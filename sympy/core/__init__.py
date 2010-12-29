"""Core module. Provides the basic operations needed in sympy.
"""

from sympify import sympify
from cache import cacheit
from basic import Basic, Atom, C
from singleton import S
from expr import Expr, AtomicExpr
from symbol import Symbol, Wild, Dummy, symbols, var
from numbers import Number, Real, Rational, Integer, NumberSymbol,\
        RealNumber, igcd, ilcm, seterr, E, I, nan, oo, pi, zoo
from power import Pow, integer_nthroot
from mul import Mul
from add import Add
from relational import Rel, Eq, Ne, Lt, Le, Gt, Ge, \
    Equality, Inequality, Unequality, StrictInequality
from multidimensional import vectorize
from function import Lambda, WildFunction, Derivative, diff, FunctionClass, \
    Function, expand, PoleError, count_ops, \
    expand_mul, expand_log, expand_func,\
    expand_trig, expand_complex, expand_multinomial
from sets import Set, Interval, Union, EmptySet
from evalf import PrecisionExhausted, N
from containers import Tuple
from exprtools import gcd_terms

# expose singletons
Catalan = S.Catalan
EulerGamma = S.EulerGamma
GoldenRatio = S.GoldenRatio
pure = S.Pure

