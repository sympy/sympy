from sympy import (
    Abs, Dummy, Eq, Gt,
    LambertW, Piecewise, Poly, Rational, S, Symbol, Matrix,
    acos, atan, atanh, cos, erf, erfinv, erfc, erfcinv,
    exp, log, pi, sin, sinh, sqrt, symbols,
    tan, tanh, atan2, arg,
    Lambda, imageset, cot, acot, I, EmptySet, Union, E, Interval, Intersection,
    oo)

from sympy.sets import FiniteSet, ConditionSet
from sympy.solvers.transformations import TR_product


def test_TR_product():
    x = Symbol('x', real=True)
    in_set = ConditionSet(Lambda(x, Eq(sin(x)*cos(x), 0)), S.Reals)
    out_set = Union(ConditionSet(Lambda(x, Eq(sin(x), 0)), S.Reals),
                    ConditionSet(Lambda(x, Eq(cos(x), 0)), S.Reals))
    assert TR_product(in_set) == out_set
