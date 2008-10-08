from sympy import oo, diff, log, Symbol, Piecewise

x = Symbol('x')

def test_piecewise():
    assert Piecewise(x, (0,1,x)) == Piecewise(x, (0,1,x))

    p = Piecewise(x, (-oo, -1, -1), (-1, 0, x**2), (0, oo, log(x)))
    dp = Piecewise(x, (-oo, -1, 0), (-1, 0, 2*x), (0, oo, 1/x))
    assert diff(p,x) == dp

    p_x2 = Piecewise(x, (-oo, -1, -1), (-1, 0, x**4), (0, oo, log(x**2)))
    assert p.subs(x,x**2) == p_x2
