from sympy.abc import x
from sympy.matrices import Matrix
from sympy.polys import Poly, cyclotomic_poly
from sympy.polys.numberfields.forms import StandardRep
from sympy.polys.numberfields.utilities import (
    AlgIntPowers, coeff_search, extract_fundamental_discriminant,
)


def test_AlgIntPowers_01():
    T = Poly(cyclotomic_poly(5))
    zeta_pow = AlgIntPowers(T)
    for e in range(10):
        if e % 5 < 4:
            assert Matrix([zeta_pow[e]]) == Matrix.eye(4).row(e % 5)
        else:
            assert zeta_pow[e] == [-1] * 4


def test_AlgIntPowers_02():
    T = Poly(x**3 + 2*x**2 + 3*x + 4)
    m = 7
    theta_pow = AlgIntPowers(T, m)
    for e in range(10):
        computed = theta_pow[e]
        s = StandardRep.from_poly(T, Poly(x**e, x))
        expected = [c % m for c in s.coeffs]
        assert computed == expected


def test_coeff_search():
    C = []
    search = coeff_search(2, 1)
    for i, c in enumerate(search):
        C.append(c)
        if i == 12:
            break
    assert C == [[1, 1], [1, 0], [1, -1], [0, 1], [2, 2], [2, 1], [2, 0], [2, -1], [2, -2], [1, 2], [1, -2], [0, 2], [3, 3]]


def test_extract_fundamental_discriminant():
    cases = (
        (0, {}, {0: 1}),
        (1, {}, {}),
        (8, {2: 3}, {}),
        (-8, {2: 3, -1: 1}, {}),
        (12, {2: 2, 3: 1}, {}),
        (36, {}, {2: 1, 3: 1}),
        (45, {5: 1}, {3: 1}),
        (1125, {5: 1}, {3: 1, 5: 1}),
    )
    for a, D_expected, F_expected in cases:
        D, F = extract_fundamental_discriminant(a)
        assert D == D_expected
        assert F == F_expected
