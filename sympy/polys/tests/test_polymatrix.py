from sympy.polys.polymatrix import PolyMatrix
from sympy.polys import Poly

from sympy.abc import x


def test_polymatrix():
    pm1 = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(x**3, x), Poly(-1 + x, x)]])
    v1 = PolyMatrix([[1, 0], [-1, 0]])
    assert pm1*v1 == PolyMatrix([[Poly(x**2 + x, x), Poly(0, x)], \
                                [Poly(x**3 - x + 1, x), Poly(0, x)]])
    assert v1*pm1 == PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(-x**2, x), Poly(x, x)]])

    pm2 = PolyMatrix([[Poly(x**2, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**2, x, domain='QQ'), \
                    Poly(x**3, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**3, x, domain='QQ')]])
    v2 = PolyMatrix([1, 0, 0, 0, 0, 0])
    assert pm2*v2 == PolyMatrix([[Poly(x**2, x, domain='QQ')]])
