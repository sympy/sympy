from __future__ import print_function, division

'''
Use llvmlite to create executable functions from Sympy expressions

  Prerequisites
  -------------
  LLVM - http://llvm.org/
  llvmlite - https://github.com/numba/llvmlite

Examples:

Scalar expression
-----------------
  import sympy.printing.llvmlitecode as g

  a = Symbol('a')
  e = a*a + a + 1
  e1 = g.get_jit_callable(e, args=[a])
  print(e1(1.1), e.subs('a',1.1))
'''

import sympy
import ctypes

from sympy.external import import_module

ll = import_module('llvmlite.ir').ir
llvm = import_module('llvmlite.binding').binding

llvm.initialize()
llvm.initialize_native_target()
llvm.initialize_native_asmprinter()


class LLVMJitPrinter(sympy.printing.printer.Printer):
    def __init__(self, *args, **kwargs):
        kwargs.pop('args')
        super(LLVMJitPrinter, self).__init__(*args, **kwargs)
        self.fp_type = ll.DoubleType()
        self.int_type = ll.IntType(64)

    def _print_Number(self, n, **kwargs):
        return ll.Constant(self.fp_type, float(n))

    def _print_Integer(self, expr):
        return ll.Constant(self.fp_type, float(expr.p))

    def _print_Symbol(self, s):
        # look up parameter with name s
        return self.param_dict.get(str(s))

    def _print_Pow(self, expr):
        base0 = self._print(expr.base)
        if expr.exp == sympy.S.NegativeOne:
            return self.builder.fdiv(ll.Constant(self.fp_type, 1.0), base0)
        if expr.exp == sympy.S.Half:
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

# ensure lifetime of the execution engine persists (else call to compiled
#   function will seg fault)
exe_engines = []

# ensure names for generated functions are unique
link_names = set()
current_link_suffix = 0


def llvm_jit_code(expr, args=None, link_name=None):
    global exe_engines, current_link_suffix
    module = ll.Module('mod1')
    lj = LLVMJitPrinter(args=args)

    fp_type = ll.DoubleType()
    arg_types = []
    for arg in args:
        arg_type = fp_type
        arg_types.append(arg_type)

    if not link_name:
        default_link_name = 'jit_func'
        current_link_suffix += 1
        link_name = default_link_name + str(current_link_suffix)
    if link_name in link_names:
        print("Error, name already used: %s" % link_name)
        return

    link_names.add(link_name)

    fn_type = ll.FunctionType(fp_type, arg_types)
    fn = ll.Function(module, fn_type, name=link_name)
    lj.module = module
    lj.param_dict = {}
    for i, a in enumerate(args):
        name = str(a)
        fn.args[i].name = name
        lj.param_dict[name] = fn.args[i]
    bb_entry = fn.append_basic_block('entry')

    lj.builder = ll.IRBuilder(bb_entry)
    lj.fn = fn

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


def get_jit_callable(expr, args=None, link_name=None):
    '''Create an executable function from a Sympy expression'''
    fptr = llvm_jit_code(expr, args, link_name)
    arg_ctypes = []
    for arg in args:
        arg_ctype = ctypes.c_double
        arg_ctypes.append(arg_ctype)

    cfunc = ctypes.CFUNCTYPE(ctypes.c_double, *arg_ctypes)(fptr)
    return cfunc
