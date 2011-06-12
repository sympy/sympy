from sympy import Symbol, diff
from sympy.physics.classical import DynamicSymbol

def test_dynamicsymbol():
    x = DynamicSymbol('x')
    xd = DynamicSymbol('xd')
    xdd = DynamicSymbol('xdd')
    y = DynamicSymbol('y')
    yd = DynamicSymbol('yd')
    ydd = DynamicSymbol('ydd')
    t = Symbol('t')
    assert diff(x, t) == xd
    assert diff(2 * x + 4, t) == 2 * xd
    assert diff(2 * x + 4 + y, t) == 2 * xd + yd
    assert diff(2 * x + 4 + y * x, t) == 2 * xd + x * yd + xd * y
    assert diff(2 * x + 4 + y * x, x) == 2 + y
    assert (diff(4 * x**2 + 3 * x + x * y, t) == 3 * xd + x * yd + xd * y +
            8 * x * xd)
    assert (diff(4 * x**2 + 3 * xd + x * y, t) ==  3 * xdd + x * yd + xd * y +
            8 * x * xd)
    assert diff(4 * x**2 + 3 * xd + x * y, xd) == 3
    assert diff(4 * x**2 + 3 * xd + x * y, xdd) == 0
