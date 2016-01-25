from __future__ import print_function, division

'''
Use llvmlite to create executable functions from Sympy expressions

This module requires llvmlite (https://github.com/numba/llvmlite).
'''

import ctypes

from sympy.external import import_module
from sympy.printing.printer import Printer
from sympy import S
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

    def _print_Number(self, n, **kwargs):
        return ll.Constant(self.fp_type, float(n))

    def _print_Integer(self, expr):
        return ll.Constant(self.fp_type, float(expr.p))

    def _print_Symbol(self, s):
        # look up parameter with name s
        return self.func_arg_map.get(s)

    def _print_Pow(self, expr):
        base0 = self._print(expr.base)
        if expr.exp == S.NegativeOne:
            return self.builder.fdiv(ll.Constant(self.fp_type, 1.0), base0)
        if expr.exp == S.Half:
            fn_type = ll.FunctionType(self.fp_type, [self.fp_type])
            fn = ll.Function(self.module, fn_type, "sqrt")
            return self.builder.call(fn, [base0], "sqrt")
        if expr.exp == 2:
            return self.builder.fmul(base0, base0)

        exp0 = self._print(expr.exp)
        fn_type = ll.FunctionType(self.fp_type, [self.fp_type, self.fp_type])
        fn = ll.Function(self.module, fn_type, "pow")
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
        fn_type = ll.FunctionType(self.fp_type, [self.fp_type])
        fn = ll.Function(self.module, fn_type, name)
        return self.builder.call(fn, [e0], name)

    def emptyPrinter(self, expr):
        raise TypeError("Unsupported type for LLVM JIT conversion: %s"
                        % type(expr))


# ensure lifetime of the execution engine persists (else call to compiled
#   function will seg fault)
exe_engines = []

# ensure names for generated functions are unique
link_names = set()
current_link_suffix = 0


def llvm_jit_code(expr, func_args=None):
    """Create a native code function from a Sympy expression"""

    global exe_engines, current_link_suffix
    module = ll.Module('mod1')

    fp_type = ll.DoubleType()
    arg_types = []
    for arg in func_args:
        arg_type = fp_type
        arg_types.append(arg_type)

    default_link_name = 'jit_func'
    current_link_suffix += 1
    link_name = default_link_name + str(current_link_suffix)
    link_names.add(link_name)

    fn_type = ll.FunctionType(fp_type, arg_types)
    fn = ll.Function(module, fn_type, name=link_name)

    param_dict = {}
    for i, a in enumerate(func_args):
        fn.args[i].name = str(a)
        param_dict[a] = fn.args[i]

    bb_entry = fn.append_basic_block('entry')

    builder = ll.IRBuilder(bb_entry)

    lj = LLVMJitPrinter(module, builder, fn, func_arg_map=param_dict)

    ret = lj._print(expr)
    lj.builder.ret(ret)

    strmod = str(module)
    if False:
        print("LLVM IR")
        print(strmod)

    llmod = llvm.parse_assembly(strmod)

    pmb = llvm.create_pass_manager_builder()
    pmb.opt_level = 2
    pass_manager = llvm.create_module_pass_manager()
    pmb.populate(pass_manager)

    pass_manager.run(llmod)

    target_machine = llvm.Target.from_default_triple().create_target_machine()
    exe_eng = llvm.create_mcjit_compiler(llmod, target_machine)
    exe_eng.finalize_object()
    exe_engines.append(exe_eng)

    if False:
        print("Assembly")
        print(target_machine.emit_assembly(llmod))

    fptr = exe_eng.get_pointer_to_function(llmod.get_function(link_name))

    return fptr


@doctest_depends_on(modules=('llvmlite',))
def llvm_callable(expr, func_args=None):
    '''Compile function from a Sympy expression

    Expressions are evaluated using double precision arithmetic.
    Some single argument math functions (exp, sin, cos, etc.) are supported
    in expressions.

    Parameters
    ==========
    expr : Expr
        Expression to compile.
    func_args : List of Symbol
        Arguments to the generated function.  Usually the free symbols in
        the expression.  Currently each one is assumed to convert to
        a double precision scalar.

    Returns
    =======
    Compiled function that can evaluate the expression.

    Examples
    ========
    >>> import sympy.printing.llvmjitcode as jit
    >>> from sympy.abc import a
    >>> e = a*a + a + 1
    >>> e1 = jit.llvm_callable(e, func_args=[a])
    >>> e.subs('a',1.1)   # Evaluate via substitution
    3.31000000000000
    >>> e1(1.1)  # Evaluate using JIT-compiled code
    3.3100000000000005
    '''

    if not llvmlite:
        raise ImportError("llvmlite is required for llvmjitcode")

    fptr = llvm_jit_code(expr, func_args)
    arg_ctypes = []
    for arg in func_args:
        arg_ctype = ctypes.c_double
        arg_ctypes.append(arg_ctype)

    cfunc = ctypes.CFUNCTYPE(ctypes.c_double, *arg_ctypes)(fptr)
    return cfunc
