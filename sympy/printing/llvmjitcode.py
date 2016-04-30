from __future__ import print_function, division

'''
Use llvmlite to create executable functions from Sympy expressions

This module requires llvmlite (https://github.com/numba/llvmlite).
'''

import ctypes

from sympy.external import import_module
from sympy.printing.printer import Printer
from sympy import S, IndexedBase
from sympy.utilities.decorator import doctest_depends_on

llvmlite = import_module('llvmlite')
if llvmlite:
    ll = import_module('llvmlite.ir').ir
    llvm = import_module('llvmlite.binding').binding
    llvm.initialize()
    llvm.initialize_native_target()
    llvm.initialize_native_asmprinter()


class LLVMJitPrinter(Printer):
    '''Convert expressions to LLVM IR'''
    def __init__(self, module, builder, fn, *args, **kwargs):
        self.func_arg_map = kwargs.pop("func_arg_map", {})
        if not llvmlite:
            raise ImportError("llvmlite is required for LLVMJITPrinter")
        super(LLVMJitPrinter, self).__init__(*args, **kwargs)
        self.fp_type = ll.DoubleType()
        self.module = module
        self.builder = builder
        self.fn = fn
        self.ext_fn = {}  # keep track of wrappers to external functions
        self.tmp_var = {}

    def _add_tmp_var(self, name, value):
        self.tmp_var[name] = value

    def _print_Number(self, n, **kwargs):
        return ll.Constant(self.fp_type, float(n))

    def _print_Integer(self, expr):
        return ll.Constant(self.fp_type, float(expr.p))

    def _print_Symbol(self, s):
        val = self.tmp_var.get(s)
        if not val:
            # look up parameter with name s
            val = self.func_arg_map.get(s)
        if not val:
            raise LookupError("Symbol not found: %s" % s)
        return val

    def _print_Pow(self, expr):
        base0 = self._print(expr.base)
        if expr.exp == S.NegativeOne:
            return self.builder.fdiv(ll.Constant(self.fp_type, 1.0), base0)
        if expr.exp == S.Half:
            fn = self.ext_fn.get("sqrt")
            if not fn:
                fn_type = ll.FunctionType(self.fp_type, [self.fp_type])
                fn = ll.Function(self.module, fn_type, "sqrt")
                self.ext_fn["sqrt"] = fn
            return self.builder.call(fn, [base0], "sqrt")
        if expr.exp == 2:
            return self.builder.fmul(base0, base0)

        exp0 = self._print(expr.exp)
        fn = self.ext_fn.get("pow")
        if not fn:
            fn_type = ll.FunctionType(self.fp_type, [self.fp_type, self.fp_type])
            fn = ll.Function(self.module, fn_type, "pow")
            self.ext_fn["pow"] = fn
        return self.builder.call(fn, [base0, exp0], "pow")

    def _print_Mul(self, expr):
        nodes = [self._print(a) for a in expr.args]
        e = nodes[0]
        for node in nodes[1:]:
            e = self.builder.fmul(e, node)
        return e

    def _print_Add(self, expr):
        nodes = [self._print(a) for a in expr.args]
        e = nodes[0]
        for node in nodes[1:]:
            e = self.builder.fadd(e, node)
        return e

    # TODO - assumes all called functions take one double precision argument.
    #        Should have a list of math library functions to validate this.
    def _print_Function(self, expr):
        name = expr.func.__name__
        e0 = self._print(expr.args[0])
        fn = self.ext_fn.get(name)
        if not fn:
            fn_type = ll.FunctionType(self.fp_type, [self.fp_type])
            fn = ll.Function(self.module, fn_type, name)
            self.ext_fn[name] = fn
        return self.builder.call(fn, [e0], name)

    def emptyPrinter(self, expr):
        raise TypeError("Unsupported type for LLVM JIT conversion: %s"
                        % type(expr))


# Used when parameters are passed by array.  Often used in callbacks to
# handle a variable number of parameters.
class LLVMJitCallbackPrinter(LLVMJitPrinter):
    def __init__(self, *args, **kwargs):
        super(LLVMJitCallbackPrinter, self).__init__(*args, **kwargs)

    def _print_Indexed(self, expr):
        array, idx = self.func_arg_map[expr.base]
        offset = int(expr.indices[0].evalf())
        array_ptr = self.builder.gep(array, [ll.Constant(ll.IntType(32), offset)])
        fp_array_ptr = self.builder.bitcast(array_ptr, ll.PointerType(self.fp_type))
        value = self.builder.load(fp_array_ptr)
        return value

    def _print_Symbol(self, s):
        val = self.tmp_var.get(s)
        if val:
            return val

        array, idx = self.func_arg_map.get(s, [None, 0])
        if not array:
            raise LookupError("Symbol not found: %s" % s)
        array_ptr = self.builder.gep(array, [ll.Constant(ll.IntType(32), idx)])
        fp_array_ptr = self.builder.bitcast(array_ptr,
                                            ll.PointerType(self.fp_type))
        value = self.builder.load(fp_array_ptr)
        return value


class _LoopCreator(object):
    def __init__(self, builder, func):
        self.builder = builder
        self.func = func

    def start(self, start_value, end_value, step_value):
        self.start_value = start_value
        self.end_value = end_value
        self.step_value = step_value
        self.pre_header_block = self.func.append_basic_block()
        self.builder.branch(self.pre_header_block)
        self.loop_block = self.func.append_basic_block('loop')
        self.exit_block = self.func.append_basic_block('afterloop')
        self.builder.position_at_end(self.pre_header_block)
        self.builder.branch(self.loop_block)

        self.builder.position_at_end(self.loop_block)
        self.index = self.builder.phi(ll.IntType(32))
        self.index.add_incoming(self.start_value, self.pre_header_block)

        return self.index

    def add_next(self):
        self.next_value = self.builder.add(self.index, self.step_value, 'next')
        self.index.add_incoming(self.next_value, self.loop_block)

    def finish(self):
        end_compare = self.builder.icmp_unsigned('<', self.next_value, self.end_value, 'loopcond')
        self.br = self.builder.cbranch(end_compare, self.loop_block, self.exit_block)
        self.builder.position_at_end(self.exit_block)


# ndim - dimension of integrate
# npts - number of points to evaluate
# x[i*ndim + j] (i in npts, j in ndim)
class LLVMJitCallbackVectorPrinter(LLVMJitPrinter):
    def __init__(self, *args, **kwargs):
        self.ndim = kwargs.pop("ndim", None)
        self.index = kwargs.pop("index", None)
        super(LLVMJitCallbackVectorPrinter, self).__init__(*args, **kwargs)

    def _print_Indexed(self, expr):
        array, idx = self.func_arg_map[expr.base]
        offset = int(expr.indices[0].evalf())
        array_ptr = self.builder.gep(array, [ll.Constant(ll.IntType(32), offset)])
        fp_array_ptr = self.builder.bitcast(array_ptr, ll.PointerType(self.fp_type))
        value = self.builder.load(fp_array_ptr)
        return value

    def _print_Symbol(self, s):
        val = self.tmp_var.get(s)
        if val:
            return val

        array, idx = self.func_arg_map.get(s, [None, 0])
        if not array:
            raise LookupError("Symbol not found: %s" % s)

        # x[i*ndim + idx] - need i and ndim
        e = self.builder.mul(self.index, self.ndim)
        offset = self.builder.add(e, ll.Constant(ll.IntType(32), idx))
        array_ptr = self.builder.gep(array, [offset], inbounds=True)
        fp_array_ptr = self.builder.bitcast(array_ptr,
                                            ll.PointerType(self.fp_type))
        value = self.builder.load(fp_array_ptr)
        return value


# ensure lifetime of the execution engine persists (else call to compiled
#   function will seg fault)
exe_engines = []

# ensure names for generated functions are unique
link_names = set()
current_link_suffix = 0


class LLVMJitCode(object):
    def __init__(self, signature):
        self.signature = signature
        self.fp_type = ll.DoubleType()
        self.module = ll.Module('mod1')
        self.fn = None
        self.llvm_arg_types = []
        self.llvm_ret_type = self.fp_type
        self.param_dict = {}  # map symbol name to LLVM function argument
        self.link_name = ''

    def _from_ctype(self, ctype):
        if ctype == ctypes.c_int:
            return ll.IntType(32)
        if ctype == ctypes.c_uint:
            return ll.IntType(32)
        if ctype == ctypes.c_size_t:
            return ll.IntType(64)
        if ctype == ctypes.c_ulong:
            return ll.IntType(64)
        if ctype == ctypes.c_double:
            return self.fp_type
        if ctype == ctypes.POINTER(ctypes.c_double):
            return ll.PointerType(self.fp_type)
        if ctype == ctypes.POINTER(ctypes.c_int):
            return ll.PointerType(ll.IntType(32))
        if ctype == ctypes.c_void_p:
            return ll.PointerType(ll.IntType(32))
        if ctype == ctypes.py_object:
            return ll.PointerType(ll.IntType(32))

        print("Unhandled ctype = %s" % str(ctype))

    def _create_args(self, func_args):
        """Create types for function arguments"""
        self.llvm_ret_type = self._from_ctype(self.signature.ret_type)
        self.llvm_arg_types = \
            [self._from_ctype(a) for a in self.signature.arg_ctypes]

    def _create_function_base(self):
        """Create function with name and type signature"""
        global link_names, current_link_suffix
        default_link_name = 'jit_func'
        current_link_suffix += 1
        self.link_name = default_link_name + str(current_link_suffix)
        link_names.add(self.link_name)

        fn_type = ll.FunctionType(self.llvm_ret_type, self.llvm_arg_types)
        self.fn = ll.Function(self.module, fn_type, name=self.link_name)

    def _create_param_dict(self, func_args):
        """Mapping of symbolic values to function arguments"""
        for i, a in enumerate(func_args):
            self.fn.args[i].name = str(a)
            self.param_dict[a] = self.fn.args[i]

    def _wrap_return(self, lj, vals):
        # Return a single double if there is one return value,
        #  else return a tuple of doubles.

        # Don't wrap return value in this case
        if self.signature.ret_type == ctypes.c_double:
            return vals[0]

        # Use this instead of a real PyObject*
        void_ptr = ll.PointerType(ll.IntType(32))

        # Create a wrapped double: PyObject* PyFloat_FromDouble(double v)
        wrap_type = ll.FunctionType(void_ptr, [self.fp_type])
        wrap_fn = ll.Function(lj.module, wrap_type, "PyFloat_FromDouble")

        wrapped_vals = [lj.builder.call(wrap_fn, [v]) for v in vals]
        if len(vals) == 1:
            final_val = wrapped_vals[0]
        else:
            # Create a tuple: PyObject* PyTuple_Pack(Py_ssize_t n, ...)

            # This should be Py_ssize_t
            tuple_arg_types = [ll.IntType(32)]

            tuple_arg_types.extend([void_ptr]*len(vals))
            tuple_type = ll.FunctionType(void_ptr, tuple_arg_types)
            tuple_fn = ll.Function(lj.module, tuple_type, "PyTuple_Pack")

            tuple_args = [ll.Constant(ll.IntType(32), len(wrapped_vals))]
            tuple_args.extend(wrapped_vals)

            final_val = lj.builder.call(tuple_fn, tuple_args)

        return final_val

    def _create_function(self, expr):
        """Create function body and return LLVM IR"""
        bb_entry = self.fn.append_basic_block('entry')
        builder = ll.IRBuilder(bb_entry)

        lj = LLVMJitPrinter(self.module, builder, self.fn,
                            func_arg_map=self.param_dict)

        ret_vals = self._convert_expr(lj, expr)

        lj.builder.ret(self._wrap_return(lj, ret_vals))

        strmod = str(self.module)
        return strmod

    def _convert_expr(self, lj, expr):
        try:
            # Match CSE return data structure.
            if len(expr) == 2:
                tmp_exprs = expr[0]
                final_exprs = expr[1]
                if len(final_exprs) != 1 and self.signature.ret_type == ctypes.c_double:
                    raise NotImplementedError("Return of multiple expressions not supported for this callback")
                for name, e in tmp_exprs:
                    val = lj._print(e)
                    lj._add_tmp_var(name, val)
        except TypeError:
            final_exprs = [expr]

        vals = [lj._print(e) for e in final_exprs]

        return vals

    def _compile_function(self, strmod):
        global exe_engines
        llmod = llvm.parse_assembly(strmod)

        pmb = llvm.create_pass_manager_builder()
        pmb.opt_level = 2
        pass_manager = llvm.create_module_pass_manager()
        pmb.populate(pass_manager)

        pass_manager.run(llmod)

        target_machine = \
            llvm.Target.from_default_triple().create_target_machine()
        exe_eng = llvm.create_mcjit_compiler(llmod, target_machine)
        exe_eng.finalize_object()
        exe_engines.append(exe_eng)

        if False:
            print("Assembly")
            print(target_machine.emit_assembly(llmod))

        fptr = exe_eng.get_pointer_to_function(
            llmod.get_function(self.link_name))

        return fptr


class LLVMJitCodeCallback(LLVMJitCode):
    def __init__(self, signature):
        super(LLVMJitCodeCallback, self).__init__(signature)

    def _create_param_dict(self, func_args):
        for i, a in enumerate(func_args):
            if isinstance(a, IndexedBase):
                self.param_dict[a] = (self.fn.args[i], i)
                self.fn.args[i].name = str(a)
            else:
                self.param_dict[a] = (self.fn.args[self.signature.input_arg],
                                      i)

    def _create_function(self, expr):
        """Create function body and return LLVM IR"""
        bb_entry = self.fn.append_basic_block('entry')
        builder = ll.IRBuilder(bb_entry)

        lj = LLVMJitCallbackPrinter(self.module, builder, self.fn,
                                    func_arg_map=self.param_dict)

        ret = self._convert_expr(lj, expr)

        if self.signature.ret_arg:
            output_fp_ptr = builder.bitcast(self.fn.args[self.signature.ret_arg],
                                            ll.PointerType(self.fp_type))
            for i, val in enumerate(ret):
                index = ll.Constant(ll.IntType(32), i)
                output_array_ptr = builder.gep(output_fp_ptr, [index])
                builder.store(val, output_array_ptr)
            builder.ret(ll.Constant(ll.IntType(32), 0))  # return success
        else:
            lj.builder.ret(self._wrap_return(lj, ret))

        strmod = str(self.module)
        return strmod


class LLVMJitCodeCallbackVector(LLVMJitCode):
    def __init__(self, signature):
        super(LLVMJitCodeCallbackVector, self).__init__(signature)

    def _create_param_dict(self, func_args):
        for i, a in enumerate(func_args):
            if isinstance(a, IndexedBase):
                self.param_dict[a] = (self.fn.args[i], i)
                self.fn.args[i].name = str(a)
            else:
                self.param_dict[a] = (self.fn.args[self.signature.input_arg],
                                      i)

    def _create_function(self, expr):
        """Create function body and return LLVM IR"""

        iargs = []
        iargs.append(self.signature.input_arg)
        iargs.append(self.signature.ret_arg)
        if self.signature.ndim_is_pointer:
            iargs.append(self.signature.ndim_arg)

        if self.signature.ncomp_is_pointer:
            iargs.append(self.signature.ncomp_arg)

        if self.signature.vector_len_is_pointer:
            iargs.append(self.signature.vector_len_arg)

        for i in iargs:
            self.fn.args[i].add_attribute("nocapture")
            self.fn.args[i].add_attribute("noalias")

        bb_entry = self.fn.append_basic_block('entry')
        builder = ll.IRBuilder(bb_entry)

        loop = _LoopCreator(builder, self.fn)

        start = ll.Constant(ll.IntType(32), 0)
        step = ll.Constant(ll.IntType(32), 1)
        vec_len = self.fn.args[self.signature.vector_len_arg]
        if self.signature.vector_len_is_pointer:
            vec_len = builder.load(vec_len)
        end = builder.trunc(vec_len, ll.IntType(32))
        index_var = loop.start(start, end, step)

        ndim = ll.Constant(ll.IntType(32), len(self.param_dict))
        lj = LLVMJitCallbackVectorPrinter(self.module, builder, self.fn,
                                          ndim=ndim, index=index_var,
                                          func_arg_map=self.param_dict)

        ret_vals = self._convert_expr(lj, expr)
        output_fp_ptr = builder.bitcast(self.fn.args[self.signature.ret_arg],
                                        ll.PointerType(self.fp_type))

        # fval[i*fdim + k]
        # i is 0..npts and k is 0..fdim
        #  index_var*ncomp + ret_vals index

        ncomp = self.fn.args[self.signature.ncomp_arg]
        if self.signature.ncomp_is_pointer:
            ncomp = builder.load(ncomp)

        # should assert that ncomp == len(ret_vals)

        idx1 = builder.mul(index_var, ll.Constant(ll.IntType(32), len(ret_vals)))
        for val_idx, ret in enumerate(ret_vals):
            idx2 = builder.add(idx1, ll.Constant(ll.IntType(32), val_idx))
            output_array_ptr = builder.gep(output_fp_ptr, [idx2], inbounds=True)
            builder.store(ret, output_array_ptr)

        loop.add_next()
        loop.finish()

        builder.ret(ll.Constant(ll.IntType(32), 0))  # return success

        strmod = str(self.module)
        return strmod


class CodeSignature(object):
    def __init__(self, ret_type):
        self.ret_type = ret_type
        self.arg_ctypes = []

        # Input argument array element index
        self.input_arg = 0

        # For the case output value is referenced through a parameter rather
        # than the return value
        self.ret_arg = None

        # Number of dimensions (arguments)
        self.ndim_arg = 0
        self.ndim_is_pointer = False

        # For vector inputs, which argument is the vector length
        self.vector_len_arg = None
        self.vector_len_is_pointer = False

        # Number of components (return values)
        self.ncomp_arg = 0
        self.ncomp_is_pointer = False


def _llvm_jit_code(args, expr, signature, callback_type):
    """Create a native code function from a Sympy expression"""
    if callback_type is None:
        jit = LLVMJitCode(signature)
    else:
        if signature.vector_len_arg is None:
            jit = LLVMJitCodeCallback(signature)
        else:
            jit = LLVMJitCodeCallbackVector(signature)

    jit._create_args(args)
    jit._create_function_base()
    jit._create_param_dict(args)
    strmod = jit._create_function(expr)
    if False:
        print("LLVM IR")
        print(strmod)
    fptr = jit._compile_function(strmod)
    return fptr


@doctest_depends_on(modules=('llvmlite', 'scipy'))
def llvm_callable(args, expr, callback_type=None):
    '''Compile function from a Sympy expression

    Expressions are evaluated using double precision arithmetic.
    Some single argument math functions (exp, sin, cos, etc.) are supported
    in expressions.

    Parameters
    ==========
    args : List of Symbol
        Arguments to the generated function.  Usually the free symbols in
        the expression.  Currently each one is assumed to convert to
        a double precision scalar.
    expr : Expr, or (Replacements, Expr) as returned from 'cse'
        Expression to compile.
    callback_type : string
        Create function with signature appropriate to use as a callback.
        Currently supported:
           'scipy.integrate'
           'scipy.integrate.test'
           'cubature'
           'cubature_v'
           'cuba'
           'cuba_v'

    Returns
    =======
    Compiled function that can evaluate the expression.

    Examples
    ========
    >>> import sympy.printing.llvmjitcode as jit
    >>> from sympy.abc import a
    >>> e = a*a + a + 1
    >>> e1 = jit.llvm_callable([a], e)
    >>> e.subs(a, 1.1)   # Evaluate via substitution
    3.31000000000000
    >>> e1(1.1)  # Evaluate using JIT-compiled code
    3.3100000000000005


    Callbacks for integration functions can be JIT compiled.
    >>> import sympy.printing.llvmjitcode as jit
    >>> from sympy.abc import a
    >>> from sympy import integrate
    >>> from scipy.integrate import quad
    >>> e = a*a
    >>> e1 = jit.llvm_callable([a], e, callback_type='scipy.integrate')
    >>> integrate(e, (a, 0.0, 2.0))
    2.66666666666667
    >>> quad(e1, 0.0, 2.0)[0]
    2.66666666666667

    The 'cubature' callback is for the Python wrapper around the
    cubature package ( https://github.com/saullocastro/cubature )
    and ( http://ab-initio.mit.edu/wiki/index.php/Cubature )

    The 'cuba' callback is for the Python wrapper around the Cuba
    library (https://github.com/JohannesBuchner/PyMultiNest , in the
    'pycuba' subdirectory.)
    The original Cuba library is here: http://www.feynarts.de/cuba/ .
    A version that builds the shared library needed for PyCuba is here:
    https://github.com/JohannesBuchner/cuba .

    There are two signatures for the SciPy integration callbacks.
    The first ('scipy.integrate') is the function to be passed to the
    integration routine, and will pass the signature checks.
    The second ('scipy.integrate.test') is only useful for directly calling
    the function using ctypes variables. It will not pass the signature checks
    for scipy.integrate.

    The return value from the cse module can also be compiled.  This
    can improve the performance of the compiled function.  If multiple
    expressions are given to cse, the compiled function returns a tuple.
    The 'cubature' and 'cuba' callbacks can handle multiple expressions (set
    `fdim` to match in the cubature integration call, or 'ncomp' to match in
    the cuba integration call.)

    For improved performance, the callback routines can operate on a vector of
    input values.  Use the '_v' variant ('cubature_v' or 'cuba_v').  The use of
    the vectorized callback also needs to be enabled in the integration call.
    For cubature, set 'vectorized=True' in the call.  For cuba, set 'nvec' to
    an integer (maximum vector length.)

    >>> import sympy.printing.llvmjitcode as jit
    >>> from sympy import cse, exp
    >>> from sympy.abc import x,y
    >>> e1 = x*x + y*y
    >>> e2 = 4*(x*x + y*y) + 8.0
    >>> after_cse = cse([e1,e2])
    >>> after_cse
    ([(x0, x**2), (x1, y**2)], [x0 + x1, 4*x0 + 4*x1 + 8.0])
    >>> j1 = jit.llvm_callable([x,y], after_cse)
    >>> j1(1.0, 2.0)
    (5.0, 28.0)
    '''

    if not llvmlite:
        raise ImportError("llvmlite is required for llvmjitcode")

    signature = CodeSignature(ctypes.py_object)

    arg_ctypes = []
    if callback_type is None:
        for arg in args:
            arg_ctype = ctypes.c_double
            arg_ctypes.append(arg_ctype)
    elif callback_type == 'scipy.integrate' or callback_type == 'scipy.integrate.test':
        signature.ret_type = ctypes.c_double
        arg_ctypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        arg_ctypes_formal = [ctypes.c_int, ctypes.c_double]
        signature.input_arg = 1
    elif callback_type == 'cubature':
        arg_ctypes = [ctypes.c_int,                     # unsigned ndim
                      ctypes.POINTER(ctypes.c_double),  # const double *x
                      ctypes.c_void_p,                  # void *fdata
                      ctypes.c_int,                     # unsigned fdim
                      ctypes.POINTER(ctypes.c_double)   # double *fval
                      ]
        signature.ret_type = ctypes.c_int
        signature.input_arg = 1
        signature.ret_arg = 4
    elif callback_type == 'cubature_v':
        arg_ctypes = [ctypes.c_int,                     # unsigned ndim
                      ctypes.c_size_t,                  # size_t npts
                      ctypes.POINTER(ctypes.c_double),  # const double *x
                      ctypes.c_void_p,                  # void *fdata
                      ctypes.c_int,                     # unsigned fdim
                      ctypes.POINTER(ctypes.c_double)   # double *fval
                      ]
        signature.ret_type = ctypes.c_int
        signature.input_arg = 2
        signature.ret_arg = 5
        signature.ncomp_arg = 4
        signature.vector_len_arg = 1
    elif callback_type == 'cuba' or callback_type == 'cuba_v':
        arg_ctypes = [ctypes.POINTER(ctypes.c_int),     # int* ndim
                      ctypes.POINTER(ctypes.c_double),  # const double *x
                      ctypes.POINTER(ctypes.c_int),     # int* ncomp
                      ctypes.POINTER(ctypes.c_double),  # double *f
                      ctypes.c_void_p,                  # void *fdata
                      ctypes.POINTER(ctypes.c_int),     # int* nvec
                      ctypes.POINTER(ctypes.c_int),     # int* core
                      ]
        signature.ret_type = ctypes.c_int
        signature.input_arg = 1

        signature.ndim_arg = 0
        signature.ndim_is_pointer = True

        signature.ncomp_arg = 2
        signature.ncomp_is_pointer = True

        signature.ret_arg = 3
        if callback_type == 'cuba_v':
            signature.vector_len_arg = 5
            signature.vector_len_is_pointer = True
    else:
        raise ValueError("Unknown callback type: %s" % callback_type)

    signature.arg_ctypes = arg_ctypes

    fptr = _llvm_jit_code(args, expr, signature, callback_type)

    if callback_type and callback_type == 'scipy.integrate':
        arg_ctypes = arg_ctypes_formal

    cfunc = ctypes.CFUNCTYPE(signature.ret_type, *arg_ctypes)(fptr)
    return cfunc
