from sympy import symbols, DiracDelta, Heaviside, nan, oo, sqrt, pi, conjugate
x,y = symbols('x y')

def test_DiracDelta():
    assert DiracDelta(1) == 0
    assert DiracDelta(5.1) == 0
    assert DiracDelta(-pi) == 0
    assert DiracDelta(5,7) == 0
    assert DiracDelta(0) == oo
    assert DiracDelta(0,5) == oo
    assert DiracDelta(x).func == DiracDelta
    assert conjugate(DiracDelta(x)) == DiracDelta(x)
    assert conjugate(DiracDelta(x - y)) == DiracDelta(x - y)

def test_heaviside():
    assert Heaviside(0) == 0.5
    assert Heaviside(-5) == 0
    assert Heaviside(1) == 1
