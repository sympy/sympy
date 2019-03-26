from matchpy import (Operation, CommutativeOperation, AssociativeOperation,
    ManyToOneReplacer, OneIdentityOperation, CustomConstraint)
from matchpy.expressions.functions import register_operation_iterator, register_operation_factory
from sympy import Pow, Add, Integral, Basic, Mul, S, Function, E
from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf,
    exp as sym_exp, gamma, acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh,
    tanh, coth, sech, csch, atan, acsc, asin, acot, acos, asec, fresnels,
    fresnelc, erfc, erfi, Ei, uppergamma, polylog, zeta, factorial, polygamma, digamma, li,
    expint, LambertW, loggamma)
from sympy.integrals.rubi.utility_function import (Gamma, exp, log, ProductLog, PolyGamma)


def operation_register_integral():
    Operation.register(Integral)
    register_operation_iterator(Integral, lambda a: (a._args[0],) + a._args[1], lambda a: len((a._args[0],) + a._args[1]))

def operation_register_pow():
    Operation.register(Pow)
    OneIdentityOperation.register(Pow)
    register_operation_iterator(Pow, lambda a: a._args, lambda a: len(a._args))

def operation_register_add():
    Operation.register(Add)
    OneIdentityOperation.register(Add)
    CommutativeOperation.register(Add)
    AssociativeOperation.register(Add)
    register_operation_iterator(Add, lambda a: a._args, lambda a: len(a._args))

def operation_register_mul():
    Operation.register(Mul)
    OneIdentityOperation.register(Mul)
    CommutativeOperation.register(Mul)
    AssociativeOperation.register(Mul)
    register_operation_iterator(Mul, lambda a: a._args, lambda a: len(a._args))

sympy_functions = [log, sin, cos, tan, cot, csc, sec, sqrt, erf,
    sym_exp, gamma, acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh,
    tanh, coth, sech, csch, atan, acsc, asin, acot, acos, asec, fresnels,
    fresnelc, erfc, erfi, Ei, uppergamma, polylog, zeta, factorial, polygamma, 
    digamma, li, expint, LambertW, loggamma, Gamma, exp, log, ProductLog, 
    PolyGamma
]

def operation_register_sympy_function(sympy_function):
    # print(sympy_function.__class__)
    Operation.register(sympy_function.__class__)
    register_operation_iterator(sympy_function, lambda a: a._args, lambda a: len(a._args))

def operation_register_sympy_functions():
    import multiprocessing as mp
    # Number of processors is equal to mp.cpu_count
    pool = mp.Pool(mp.cpu_count())
    # print(mp.cpu_count())
    pool.map(operation_register_sympy_function, sympy_functions)
    pool.close()

def sympy_op_factory(old_operation, new_operands, variable_name):
    return type(old_operation)(*new_operands)

def register_operation_factory_sympy():
    register_operation_factory(Basic, sympy_op_factory)

def operation_register():
    operation_register_integral()
    operation_register_mul()
    operation_register_add()
    operation_register_pow()

    # parallel execution
    operation_register_sympy_functions()

    register_operation_factory_sympy()
