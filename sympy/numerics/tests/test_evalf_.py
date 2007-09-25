import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics.evalf_ import evalf, polyfunc
from sympy.utilities.pytest import XFAIL

@XFAIL
def test_simple_evalf():
    assert evalf(2) == 2
    assert evalf(Rational(1,2)) == 0.5
    assert evalf(exp(pi*I)).ae(-1)
    assert evalf(pi+E).ae(5.8598744820488378)

def test_polyfunc():
    x = Symbol('x')
    p = polyfunc(2*x**3 - 2.5*x + 4)
    q = lambda t: 2*t**3 - 2.5*t + 4
    assert p(0).ae(q(0))
    assert p(1).ae(q(1))
    assert p(-7.5).ae(q(-7.5))
    p = x**3 - 5*x**2 + 4*x - 6
    pd = p.diff(x)
    assert polyfunc(p)(2) == p.subs(x, 2)
    assert polyfunc(p, True)(2) == (p.subs(x, 2), pd.subs(x, 2))
