from sympy.algebras.quaternion import Quaternion
from sympy import symbols, re, im, Add, Mul, I


x, y, z, w = symbols("x y z w")


def test_quaternion_construction():
    q = Quaternion(x, y, z, w)
    assert q + q == Quaternion(2*x, 2*y, 2*z, 2*w)


def test_quaternion_complex_real_addition():
    a = symbols("a", complex=True)
    b = symbols("b", real=True)
    # This symbol is not complex:
    c = symbols("c", commutative=False)

    q = Quaternion(x, y, z, w)
    assert a + q == Quaternion(x+re(a), y+im(a), z, w)
    assert 1 + q == Quaternion(1+x, y, z, w)
    assert I + q == Quaternion(x, 1+y, z, w)
    assert b + q == Quaternion(x + b, y, z, w)
    assert c + q == Add(c, Quaternion(x, y, z, w), evaluate=False)
    assert q * c == Mul(Quaternion(x, y, z, w), c, evaluate=False)
    assert c * q == Mul(c, Quaternion(x, y, z, w), evaluate=False)

    assert -q == Quaternion(-x, -y, -z, -w)
