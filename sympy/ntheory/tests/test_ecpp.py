from sympy.ntheory import hilbert, c_smith
from sympy import Poly, Symbol

x = Symbol('x')
assert hilbert(-3) == (1, Poly(x, x, domain='ZZ'), {(1, 1, 1.0)})
assert hilbert(-4) == (1, Poly(x - 1728, x, domain='ZZ'), {(1, 0, 1.0)})
assert hilbert(-15) == (2, Poly(x**2 + 191025*x - 121287375, x, domain='ZZ'), {(1, 1, 4.0), (2, 1, 2.0)})
assert hilbert(-23) == (3, Poly(x**3 + 3491750*x**2 - 5151296875*x + 12771880859375, x, domain='ZZ'), {(1, 1, 6.0), (2, -1, 3.0), (2, 1, 3.0)})

assert c_smith(1097, -92) == [(54, 4), (-54, -4)]
assert c_smith(1097, -8) == [(66, 2), (-66, -2)]
assert c_smith(1097, -20) is None
assert c_smith(10009, -100) == [(6, 20), (-6, -20)]
assert c_smith(10009, -108) == [(182, 8), (-182, -8)]
assert c_smith(10009, -115) == [(119, 15), (-119, -15)]
assert c_smith(10093, -107) == [(158, 12), (-158, -12)]


assert c_smith(565039, -139) == [(135, 127), (-135, -127)]
for i in [136, 139, 140, 143, 144, 147, 148]:
    assert c_smith(10009, -i) is None
