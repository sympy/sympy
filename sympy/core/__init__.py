"""Core module. Provides the basic operations needed in sympy.
"""

from basic import Basic, Expr, S, C, sympify
from symbol import Symbol, Wild, symbols, var
from numbers import Number, Real, Rational, Integer, NumberSymbol,\
        RealNumber, igcd, ilcm, seterr
from power import Pow, integer_nthroot
from mul import Mul
from add import Add
from relational import Rel, Eq, Ne, Lt, Le, Gt, Ge, \
    Equality, Inequality, Unequality, StrictInequality
from multidimensional import vectorize
from function import Lambda, WildFunction, Derivative, diff, FunctionClass, \
    Function, expand, PoleError, expand_mul, expand_log, expand_func,\
    expand_trig, expand_complex
from sets import Set, Interval, Union, EmptySet
from evalf import PrecisionExhausted, N

# expose singletons like exp, log, oo, I, etc.
for _n, _cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)
