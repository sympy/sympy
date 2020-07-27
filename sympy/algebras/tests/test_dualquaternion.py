from sympy import symbols, simplify, exp, ln, conjugate, I
from sympy.algebras import Dual, Quaternion
from sympy.testing.pytest import SKIP

symbs = symbols('a:h')
a, b, c, d, e, f, g, h = symbs[:8]
p = Quaternion(*symbs[:4])
q = Quaternion(*symbs[4:8])
r = Quaternion(*symbs[8:12])
s = Quaternion(*symbs[12:16])

def test_dualquat_construction():
    x = Dual(q, p)
    assert x + x == Dual(2*q, 2*p)

def test_dualquat_conjugation():
    x = Dual(p, q)
    xconj = Dual(conjugate(p), conjugate(q))  # Quaternion conj
    # The norm should be purely a dual number
    assert x*xconj == Dual(Quaternion(a**2 + b**2 + c**2 + d**2),
                           Quaternion(2*a*e + 2*b*f + 2*c*g + 2*d*h))

    y = x.normalize()
    # When properly normalized a dual quaternion is unitary.
    assert simplify(y*conjugate(y)) == Dual(Quaternion(1), Quaternion(0))

def test_dual_power():
    x = Dual(p,q)
    assert x**2 == Dual(p**2, p*q + q*p)
    assert simplify(x**3 - Dual(p ** 3, p*q*p + q*p*p + p*p*q)) == Dual(Quaternion(0), Quaternion(0))
    assert x**-1 == x.inverse()

def test_dual_add():
    d1 = Dual(p,q)
    d2 = Dual(r,s)
    assert d1 + d2 == Dual(p+r, q+s)
    assert d1 + r == Dual(p+r, q)

def test_dual_mul():
    d1 = Dual(p, q)
    d2 = Dual(r, s)
    assert d1 * d2 == Dual(p*r, p*s+q*r)
    assert d1 * r == Dual(p*r, q*r)
    assert Dual._generic_mul(r, d1) == Dual(r*p, r*q)
    assert r * d1 == Dual(r * p, r*q)

def test_dualquat_inv():
    d1 = Dual(p, q)
    # Inversion
    assert simplify(d1*d1.inverse()) == Dual(Quaternion(1), Quaternion(0))
    assert d1.inverse() == Dual(p**-1, -(p**-1)*q*(p**-1))

@SKIP('Closed form exp/log is yet to be implemented for Dual Quaternions.')
def test_dual_exp():
    pass

@SKIP('Closed form exp/log is yet to be implemented for Dual Quaternions.')
def test_dual_ln():
    x = Dual(a, b)
    assert exp(ln(x)) == x
