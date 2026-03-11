"""Tests for polynomial GCD improvements.

Tests for trivial case fast-paths and HeuristicGCDFailed fallback
added to the sparse polynomial GCD dispatch in rings.py.

Related issue: https://github.com/sympy/sympy/issues/23131
"""
from __future__ import annotations

from sympy.polys.rings import ring
from sympy.polys.domains import ZZ, QQ, ZZ_I


# ============================================================
# Tests for _gcd_trivial: ground element (constant) cases
# ============================================================

def test_gcd_ground_element():
    """GCD where one polynomial is a constant (ground element)."""
    R, x, y = ring("x,y", ZZ)

    # gcd(2, 4*x + 6*y) = 2
    f = R(2)
    g = 4*x + 6*y

    h = f.gcd(g)
    assert h == R(2)

    # Verify cofactors
    h, cff, cfg = f.cofactors(g)
    assert h == R(2)
    assert h * cff == f
    assert h * cfg == g


def test_gcd_ground_element_coprime():
    """GCD of a constant and a polynomial that are coprime."""
    R, x, y = ring("x,y", ZZ)

    f = R(3)
    g = 4*x + 10*y + 7

    h, cff, cfg = f.cofactors(g)
    assert h == R(1)
    assert h * cff == f
    assert h * cfg == g


def test_gcd_ground_element_one():
    """GCD where the constant is 1."""
    R, x, y, z = ring("x,y,z", ZZ)

    f = R(1)
    g = x**2 + y**2 + z

    h = f.gcd(g)
    assert h == R(1)


def test_gcd_both_ground():
    """GCD where both polynomials are constants."""
    R, x, y = ring("x,y", ZZ)

    f = R(12)
    g = R(8)

    h, cff, cfg = f.cofactors(g)
    assert h == R(4)
    assert h * cff == f
    assert h * cfg == g


# ============================================================
# Tests for _gcd_trivial: no common generators
# ============================================================

def test_gcd_no_common_generators():
    """GCD of polynomials with completely disjoint variables."""
    R, x, y, z, w = ring("x,y,z,w", ZZ)

    # f uses only x,y; g uses only z,w
    f = x**2 + 2*x*y + y**2
    g = z**3 + w

    # No common variables => gcd is gcd of contents = 1
    h = f.gcd(g)
    assert h == R(1)

    # Verify cofactors
    h, cff, cfg = f.cofactors(g)
    assert h == R(1)
    assert h * cff == f
    assert h * cfg == g


def test_gcd_no_common_generators_with_content():
    """No common variables but with common content factor."""
    R, x, y, z, w = ring("x,y,z,w", ZZ)

    # f uses only x, g uses only y
    # gcd should be gcd(content(f), content(g)) = gcd(6, 4) = 2
    f = 6*x**2 + 12*x
    g = 4*y**3 + 8*y

    h, cff, cfg = f.cofactors(g)
    assert h == R(2)
    assert h * cff == f
    assert h * cfg == g


def test_gcd_no_common_generators_many_vars():
    """No common variables with many generators (the key perf case)."""
    R, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 = ring(
        "x0,x1,x2,x3,x4,x5,x6,x7,x8,x9", ZZ
    )

    # f uses x0..x4, g uses x5..x9 - totally disjoint
    f = x0 + x1 + x2 + x3 + x4
    g = x5 + x6 + x7 + x8 + x9

    h = f.gcd(g)
    assert h == R(1)


# ============================================================
# Tests for _gcd_ZZ fallback (HeuristicGCDFailed -> modular GCD)
# ============================================================

def test_gcd_zz_basic():
    """Basic ZZ GCD still works through the normal heugcd path."""
    R, x, y = ring("x,y", ZZ)

    f = x**2 + 2*x*y + y**2  # (x + y)^2
    g = x**2 + x*y            # x*(x + y)

    h, cff, cfg = f.cofactors(g)
    assert h == x + y
    assert h * cff == f
    assert h * cfg == g


def test_gcd_zz_coprime():
    """GCD of coprime polynomials over ZZ."""
    R, x = ring("x", ZZ)

    f = x**2 + 1
    g = x + 1

    h = f.gcd(g)
    assert h == R(1)


def test_gcd_zz_univariate():
    """Univariate GCD over ZZ."""
    R, x = ring("x", ZZ)

    f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
    g = x**3 + 6*x**2 + 11*x + 6

    h = f.gcd(g)
    assert h == x**2 + 3*x + 2


def test_gcd_qq_basic():
    """Basic QQ GCD works (routes through _gcd_ZZ internally)."""
    R, x, y = ring("x,y", QQ)

    f = QQ(1, 4)*x**2 + x*y + y**2
    g = QQ(1, 2)*x**2 + x*y

    h, cff, cfg = f.cofactors(g)
    assert h * cff == f
    assert h * cfg == g


# ============================================================
# Tests for GCD with multiple polynomials sharing a factor
# ============================================================

def test_gcd_shared_factor_many_vars():
    """GCD where both polynomials share a common factor in a multi-var ring."""
    R, x, y, z = ring("x,y,z", ZZ)

    common = x + y + z
    f = common * (x - y)
    g = common * (z + 1)

    h = f.gcd(g)
    assert h == common or h == -common

    # Verify via cofactors
    h, cff, cfg = f.cofactors(g)
    assert h * cff == f
    assert h * cfg == g


def test_gcd_p_and_p_squared():
    """GCD(p, p^2) = p -- the example from issue #23131."""
    R, x, y, z = ring("x,y,z", ZZ)

    p = x + y + z
    p2 = p ** 2

    h = p.gcd(p2)
    assert h == p or h == -p

    h, cff, cfg = p.cofactors(p2)
    assert h * cff == p
    assert h * cfg == p2
