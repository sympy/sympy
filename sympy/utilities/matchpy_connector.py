"""
The objects in this module allow the usage of the MatchPy pattern matching
library on SymPy expressions.
"""
from sympy.external import import_module
from sympy.functions import (log, sin, cos, tan, cot, csc, sec, erf, gamma, uppergamma)
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
from sympy.functions.special.error_functions import fresnelc, fresnels, erfc, erfi, Ei
from sympy import (Basic, Mul, Add, Pow, Integral, exp, Symbol)
from sympy.utilities.decorator import doctest_depends_on

matchpy = import_module("matchpy")

if matchpy:
    from matchpy import Operation, CommutativeOperation, AssociativeOperation, OneIdentityOperation
    from matchpy.expressions.functions import op_iter, create_operation_expression, op_len

    Operation.register(Integral)
    Operation.register(Pow)
    OneIdentityOperation.register(Pow)

    Operation.register(Add)
    OneIdentityOperation.register(Add)
    CommutativeOperation.register(Add)
    AssociativeOperation.register(Add)

    Operation.register(Mul)
    OneIdentityOperation.register(Mul)
    CommutativeOperation.register(Mul)
    AssociativeOperation.register(Mul)

    Operation.register(exp)
    Operation.register(log)
    Operation.register(gamma)
    Operation.register(uppergamma)
    Operation.register(fresnels)
    Operation.register(fresnelc)
    Operation.register(erf)
    Operation.register(Ei)
    Operation.register(erfc)
    Operation.register(erfi)
    Operation.register(sin)
    Operation.register(cos)
    Operation.register(tan)
    Operation.register(cot)
    Operation.register(csc)
    Operation.register(sec)
    Operation.register(sinh)
    Operation.register(cosh)
    Operation.register(tanh)
    Operation.register(coth)
    Operation.register(csch)
    Operation.register(sech)
    Operation.register(asin)
    Operation.register(acos)
    Operation.register(atan)
    Operation.register(acot)
    Operation.register(acsc)
    Operation.register(asec)
    Operation.register(asinh)
    Operation.register(acosh)
    Operation.register(atanh)
    Operation.register(acoth)
    Operation.register(acsch)
    Operation.register(asech)

    @op_iter.register(Integral)  # type: ignore
    def _(operation):
        return iter((operation._args[0],) + operation._args[1])

    @op_iter.register(Basic)  # type: ignore
    def _(operation):
        return iter(operation._args)

    @op_len.register(Integral)  # type: ignore
    def _(operation):
        return 1 + len(operation._args[1])

    @op_len.register(Basic)  # type: ignore
    def _(operation):
        return len(operation._args)

    @create_operation_expression.register(Basic)
    def sympy_op_factory(old_operation, new_operands, variable_name=True):
         return type(old_operation)(*new_operands)


if matchpy:
    from matchpy import Wildcard
else:
    class Wildcard:
        def __init__(self, min_length, fixed_size, variable_name, optional):
            pass


@doctest_depends_on(modules=('matchpy',))
class _WildAbstract(Wildcard, Symbol):
    min_length = None  # abstract field required in subclasses
    fixed_size = None  # abstract field required in subclasses

    def __init__(self, variable_name=None, optional=None, **assumptions):
        min_length = self.min_length
        fixed_size = self.fixed_size
        Wildcard.__init__(self, min_length, fixed_size, str(variable_name), optional)

    def __new__(cls, variable_name=None, optional=None, **assumptions):
        cls._sanitize(assumptions, cls)
        return _WildAbstract.__xnew__(cls, variable_name, optional, **assumptions)

    def __getnewargs__(self):
        return self.min_count, self.fixed_size, self.variable_name, self.optional

    @staticmethod
    def __xnew__(cls, variable_name=None, optional=None, **assumptions):
        obj = Symbol.__xnew__(cls, variable_name, **assumptions)
        return obj

    def _hashable_content(self):
        if self.optional:
            return super()._hashable_content() + (self.min_count, self.fixed_size, self.variable_name, self.optional)
        else:
            return super()._hashable_content() + (self.min_count, self.fixed_size, self.variable_name)

    def __copy__(self) -> '_WildAbstract':
        return type(self)(variable_name=self.variable_name, optional=self.optional)


@doctest_depends_on(modules=('matchpy',))
class WildDot(_WildAbstract):
    min_length = 1
    fixed_size = True


@doctest_depends_on(modules=('matchpy',))
class WildPlus(_WildAbstract):
    min_length = 1
    fixed_size = False


@doctest_depends_on(modules=('matchpy',))
class WildStar(_WildAbstract):
    min_length = 0
    fixed_size = False
