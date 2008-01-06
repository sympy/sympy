from sympy import sin, exp, cot, sqrt, I, E, pi, symbols
from sympy.simplify import cancel

x,y = symbols('xy')

def test_has():
    assert cot(x).has(x)
    assert cot(x).has(cot)
    assert not cot(x).has(sin)
    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(cot)

def test_sin_exp_rewrite():
    assert sin(x).rewrite(sin, exp) == -I/2*(exp(I*x)-exp(-I*x))

def test_cancel():
    assert cancel((x**2-2)/(x+sqrt(2))) == x - sqrt(2)
    assert cancel((x**2-y**2)/(x-y), x, y) == x + y
    assert cancel((x**2-y)/(x-y)) == 1/(x - y)*(x**2 - y)
    assert cancel((x**2-y**2)/(x-y), x) == x + y
    assert cancel((x**2-y**2)/(x-y), y) == x + y
    assert cancel((x**2-y**2)/(x-y)) == x + y
    assert cancel((E*x+2)/(x-pi)*(x-1)) == -1/(x - pi)*(1 - x)*(2 + E*x)
    assert cancel((x**3-1)/(x-1) < sin(x)) == (1 + x + x**2 < sin(x))
    assert cancel((x**2-1)/(x-1) == (x**2+1)/(x-I), x) == (1 + x == I + x)
    assert cancel((x**2-1)/(x-1) + (x**2+1)/(x-I), x) == 1 + I + 2*x
    assert cancel((x**2-1)/(x-1) + (x**2+1)/(x-I), y) == 1/(1 - x)*(1 - x**2) + 1/(x - I)*(1 + x**2)
