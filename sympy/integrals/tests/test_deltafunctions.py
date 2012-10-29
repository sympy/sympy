from sympy import (
    cos, DiracDelta, Heaviside, Function, Integral, integrate, oo,
    pi, S, sin, symbols,
)
from sympy.integrals.deltafunctions import deltaintegrate

f = Function("f")
x_1, x_2, x, y, z = symbols("x_1 x_2 x y z")

def test_deltaintegrate_DiracDelta():
    assert deltaintegrate(DiracDelta(x), x) == Heaviside(x)
    assert deltaintegrate(DiracDelta(-x), x) == Heaviside(x)

    assert deltaintegrate(DiracDelta(x) * f(x), x) == f(0) * Heaviside(x)
    assert deltaintegrate(DiracDelta(-x) * f(x), x) == f(0) * Heaviside(x)
    assert deltaintegrate(DiracDelta(x - 1) * f(x), x) == f(1) * Heaviside(x - 1)
    assert deltaintegrate(DiracDelta(1 - x) * f(x), x) == f(1) * Heaviside(x - 1)
    assert deltaintegrate(DiracDelta(x**2 + x - 2), x) == \
        Heaviside(x - 1)/3 + Heaviside(x + 2)/3

    p = cos(x)*(DiracDelta(x) + DiracDelta(x**2 - 1))*sin(x)*(x - pi)
    assert deltaintegrate(p, x) - (-pi*(cos(1)*Heaviside(-1 + x)*sin(1)/2 - \
        cos(1)*Heaviside(1 + x)*sin(1)/2) + \
        cos(1)*Heaviside(1 + x)*sin(1)/2 + \
        cos(1)*Heaviside(-1 + x)*sin(1)/2) == 0

    p = x_2*DiracDelta(x - x_2)*DiracDelta(x_2 - x_1)
    assert integrate(p, (x_2, -oo, oo)) == x*DiracDelta(x - x_1)

    p = x*y**2*z*DiracDelta(y - x)*DiracDelta(y - z)*DiracDelta(x - z)
    assert integrate(p, (y, -oo, oo)) == x**3*z*DiracDelta(x - z)**2
    assert deltaintegrate((x + 1)*DiracDelta(2*x), x) == S(1)/2 * Heaviside(x)
    assert deltaintegrate((x + 1)*DiracDelta(2*x/3 + 4/S(9)), x) == \
        S(1)/2 * Heaviside(x + S(2)/3)

    a, b, c = symbols('a b c', commutative=False)
    assert integrate(DiracDelta(x - y)*f(x - b)*f(x - a), (x, -oo, oo)) == \
        f(y - b)*f(y - a)

    p = f(x - a)*DiracDelta(x - y)*f(x - c)*f(x - b)
    assert integrate(p, (x, -oo, oo)) == f(y - a)*f(y - c)*f(y - b)

    assert integrate(DiracDelta(x - z)*f(x - b)*f(x - a)*DiracDelta(x - y),
                     (x, -oo, oo)) == DiracDelta(y - z)*f(y - b)*f(y - a)
