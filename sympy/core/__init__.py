"""Core module. Provides the basic operations needed in sympy.
"""

from basic import Basic, S, C, sympify
from symbol import Symbol, Wild, symbols, var
from numbers import Number, Real, Rational, Integer
from power import Pow
from mul import Mul
from add import Add
from relational import Rel, Eq, Ne, Lt, Le, Gt, Ge, \
    Equality, Inequality, Unequality, StrictInequality
from multidimensional import vectorize
from function import Lambda, WildFunction, Derivative, diff, FunctionClass, \
    Function, expand
from interval import Interval

# set repr output to pretty output:
Basic.set_repr_level(1)

# expose singletons like exp, log, oo, I, etc.
for _n, _cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)

