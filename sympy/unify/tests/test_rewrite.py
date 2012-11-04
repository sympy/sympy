from sympy.unify.rewrite import rewriterule
from sympy import Wild, sin, cos
from sympy.abc import x, y, z

p, q = Wild('p'), Wild('q')

def test_simple():
    p1 = p**2
    p2 = p**3
    rl = rewriterule(p1, p2)

    expr = x**2
    assert list(rl(expr)) == [x**3]

def test_moderate():
    p1 = p**2 + q**3
    p2 = (p*q)**4
    rl = rewriterule(p1, p2)

    expr = x**2 + y**3
    assert list(rl(expr)) == [(x*y)**4]

def test_sincos():
    p1 = sin(p)**2 + sin(p)**2
    p2 = 1
    rl = rewriterule(p1, p2)

    assert list(rl(sin(x)**2 + sin(x)**2)) == [1]
    assert list(rl(sin(y)**2 + sin(y)**2)) == [1]

def test_Exprs_ok():
    rl = rewriterule(p+q, q+p)
    rl(x+y).next().is_commutative
    str(rl(x+y).next())
