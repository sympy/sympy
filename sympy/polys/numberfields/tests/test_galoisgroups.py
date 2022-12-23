"""Tests for computing Galois groups. """

from sympy.combinatorics.named_groups import (
    SymmetricGroup, AlternatingGroup, DihedralGroup, CyclicGroup,
    AbelianGroup,
)
from sympy.combinatorics.permutations import Permutation
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.numberfields.galoisgroups import (
    Resolvent, tschirnhausen_transformation,
    MaxTriesException, galois_group,
    _galois_group_degree_4_simple,
    M20, G36minus,
)
from sympy.polys.numberfields.subfield import field_isomorphism
from sympy.polys.polytools import Poly
from sympy.testing.pytest import raises, slow


def test_Resolvent_roots():
    X = symbols('X0 X1 X2 X3')
    F = X[0]*X[2] + X[1]*X[3]
    s = [
        Permutation([0, 1, 2, 3]),
        Permutation([1, 0, 2, 3]),
        Permutation([3, 1, 2, 0])
    ]
    R = Resolvent(F, X, s)
    roots = [r(*X) for r in R.root_lambdas]
    assert roots == [
        S('X0*X2 + X1*X3'),
        S('X0*X3 + X1*X2'),
        S('X0*X1 + X2*X3')
    ]


def test_tschirnhausen_transformation():
    x = symbols('x')
    for T in [
        Poly(x**2 - 2),
        Poly(x**2 + x + 1),
        Poly(x**4 + 1),
        Poly(x**4 - x**3 + x**2 - x + 1),
    ]:
        try:
            _, U = tschirnhausen_transformation(T)
            assert U.degree() == T.degree()
            assert U.is_monic
            assert U.is_irreducible
            K = QQ.alg_field_from_poly(T)
            L = QQ.alg_field_from_poly(U)
            assert field_isomorphism(K.ext, L.ext) is not None
        except MaxTriesException:
            print('Max tries exceeded on Tschirnhausen transf.')


def test_G36minus():
    G = G36minus()
    assert G.order() == 36
    A6 = AlternatingGroup(6)
    assert not G.is_subgroup(A6)


def test__galois_group_degree_4_simple():
    x = symbols('x')
    for T, G, s in [
        (x**4 + x**3 + x**2 + x + 1, CyclicGroup(4), False),
        (x**4 + 1, AbelianGroup(2, 2), True),
        #(x**4 - 24*x**3 + 608*x**2 - 448*x + 1088, AbelianGroup(2, 2), True),
        (x**4 - 2, DihedralGroup(4), False),
        (x**4 + 8*x + 12, AlternatingGroup(4), True),
        (x**4 + x + 1, SymmetricGroup(4), False),
    ]:
        assert _galois_group_degree_4_simple(Poly(T)) == (G, s)


@slow
def test_galois_group():
    x = symbols('x')

    raises(ValueError, lambda: galois_group(Poly(0, x)))
    raises(ValueError, lambda: galois_group(Poly(1, x)))

    for T, G, s in [
        # Degree 1
        (x, CyclicGroup(1), True),
        # Degree 2
        (x**2 + x + 1, CyclicGroup(2), False),
        # Degree 3
        (x**3 + x**2 - 2*x - 1, AlternatingGroup(3), True),
        (x**3 + 2, SymmetricGroup(3), False),
        # Degree 4
        (x**4 + x**3 + x**2 + x + 1, CyclicGroup(4), False),
        (x**4 + 1, AbelianGroup(2, 2), True),
        (x**4 - 2, DihedralGroup(4), False),
        (x**4 + 8*x + 12, AlternatingGroup(4), True),
        (x**4 + x + 1, SymmetricGroup(4), False),
        # Degree 5
        (x**5 + x**4 - 4*x**3 - 3*x**2 + 3*x + 1, CyclicGroup(5), True),
        (x**5 - 5*x + 12, DihedralGroup(5), True),
        (x**5 + 2, M20(), False),
        (x**5 + 20*x + 16, AlternatingGroup(5), True),
        (x**5 - x + 1, SymmetricGroup(5), False),
    ]:
        assert galois_group(Poly(T)) == (G, s)
