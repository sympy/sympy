"""Tests for high--level polynomials manipulation functions. """

from sympy.polys.polyfuncs import (
    symmetrize, horner, interpolate,
)

from sympy.abc import a, b, c, d, e, x, y, z

def test_symmetrize():
    assert symmetrize(0, x, y, z) == (0, 0)
    assert symmetrize(1, x, y, z) == (1, 0)

    s1 = x + y + z
    s2 = x*y + x*z + y*z
    s3 = x*y*z

    assert symmetrize(x) == (x, 0)
    assert symmetrize(x + 1) == (x + 1, 0)

    assert symmetrize(x, x, y) == (x + y, -y)
    assert symmetrize(x + 1, x, y) == (x + y + 1, -y)

    assert symmetrize(x, x, y, z) == (s1, -y - z)
    assert symmetrize(x + 1, x, y, z) == (s1 + 1, -y - z)

    assert symmetrize(x**2, x, y, z) == (s1**2 - 2*s2, -y**2 - z**2)

    assert symmetrize(x**2 + y**2) == (-2*x*y + (x + y)**2, 0)
    assert symmetrize(x**2 - y**2) == (-2*x*y + (x + y)**2, -2*y**2)

    assert symmetrize(x**3 + y**2 + a*x**2 + b*y**3, x, y) == \
        (-3*x*y*(x + y) - 2*a*x*y + a*(x + y)**2 + (x + y)**3, y**2*(1 - a) - y**3*(1 - b))

def test_horner():
    assert horner(0) == 0
    assert horner(1) == 1
    assert horner(x) == x

    assert horner(x + 1) == x + 1
    assert horner(x**2 + 1) == x**2 + 1
    assert horner(x**2 + x) == (x + 1)*x
    assert horner(x**2 + x + 1) == (x + 1)*x + 1

    assert horner(9*x**4 + 8*x**3 + 7*x**2 + 6*x + 5) == (((9*x + 8)*x + 7)*x + 6)*x + 5
    assert horner(a*x**4 + b*x**3 + c*x**2 + d*x + e) == (((a*x + b)*x + c)*x + d)*x + e

    assert horner(4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y, wrt=x) == ((4*y + 2)*x*y + (2*y + 1)*y)*x
    assert horner(4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y, wrt=y) == ((4*x + 2)*y*x + (2*x + 1)*x)*y

def test_interpolate():
    assert interpolate([1,4,9,16], x) == x**2
    assert interpolate([(1, 1), (2, 4), (3, 9)], x) == x**2
    assert interpolate([(1, 2), (2, 5), (3, 10)], x) == 1 + x**2
    assert interpolate({1: 2, 2: 5, 3: 10}, x) == 1 + x**2

