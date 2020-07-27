"""Common multiplication operator"""

from functools import cmp_to_key
from sympy.core import Basic, S, Tuple
from sympy.core.sympify import _sympify
from .map import AppliedMap
from .operator import (
    BinaryOperator, AppliedBinaryOperator
)

__all__ = [
    "MultiplicationOperator", "Multiplication",
    "scalar_mul", "scalar_pow",
]

class MultiplicationOperator(BinaryOperator):
    """
    Class for general scalar multiplication operator.

    Examples
    ========

    >>> from sympy import scalar_mul
    >>> from sympy.abc import x, y

    >>> scalar_mul(x, y)
    x * y

    >>> scalar_mul(x, x)
    x * x
    >>> scalar_mul(x, x, evaluate=True)
    x**2

    """
    str_name = '*'
    pretty_name = '⋅'
    latex_name = r'\times'
    is_left_divisible = is_right_divisible = True
    is_associative = True
    is_commutative = True

    def __new__(cls, domain, codomain, identity, **kwargs):
        domain, codomain, identity = _sympify(domain), _sympify(codomain), _sympify(identity)
        return super().__new__(cls, domain, codomain, identity)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    @property
    def identity(self):
        return self.args[2]

    def __call__(self, *args, evaluate=False):
        return Multiplication(self, args, evaluate=evaluate)

    def multiplication_process(self, seq):
        seq = self.flatten(seq)
        seq = self.remove_identity(seq)
        return self.assoc_comm_process(seq)

class Multiplication(AppliedBinaryOperator):
    def __new__(cls, mapping, args, evaluate=False, **kwargs):
        args = [_sympify(a) for a in args]

        if evaluate:
            args = mapping.multiplication_process(args)

            if not args:
                return mapping.identity
            elif len(args) == 1:
                return args[0]

        args.sort(key=cmp_to_key(Basic.compare))
        args = Tuple(*args)
        return super(AppliedMap, cls).__new__(cls, mapping, args)

scalar_mul = MultiplicationOperator(S.Complexes**2, S.Complexes, S.One)
scalar_pow = scalar_mul.exponent_operator()
