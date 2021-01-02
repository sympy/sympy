from itertools import combinations

from sympy import (ask, cos, exp, FiniteSet, Float, Function, I, log, oo,
    pi, Q, Rational, S, sin, sqrt, symbols)
from sympy.simplify import trigsimp
from sympy.printing import sstr, pretty, latex

x,y,z = symbols('x y z')

def test_printing():
    eq = Q.eq(x,y)
    assert sstr(eq) == "x = y"
    assert pretty(eq) == "x = y"
    assert latex(eq) == "x = y"
