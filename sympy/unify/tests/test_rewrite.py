from sympy.unify.rewrite import rewriterule
from sympy import Wild
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
    print list(rl(expr))
    assert list(rl(expr)) == [(x*y)**4]
