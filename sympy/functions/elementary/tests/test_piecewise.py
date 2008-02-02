from sympy import Piecewise, oo, log, symbols

x = symbols('x')


def test_piecewise():
    p = Piecewise(x, [[(-oo, -1), 1], [(1, 2), x*x], [(3, oo), log(x)]])


    p.fdiff()

    assert p.args[1][0] == [(-oo, -1), 0]
    assert p.args[1][1] == [(1, 2), 2*x]
    assert p.args[1][2] == [(3, oo), 1/x]
