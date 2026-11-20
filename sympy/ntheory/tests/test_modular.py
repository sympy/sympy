from __future__ import annotations
from sympy.ntheory.modular import crt, crt1, crt2, solve_congruence
from sympy.testing.pytest import raises


def test_crt():
    def mcrt(m, v, r, symmetric=False):
        assert crt(m, v, symmetric)[0] == r
        mm, e, s = crt1(m)
        assert crt2(m, v, mm, e, s, symmetric) == (r, mm)

    mcrt([2, 3, 5], [0, 0, 0], 0)
    mcrt([2, 3, 5], [1, 1, 1], 1)

    mcrt([2, 3, 5], [-1, -1, -1], -1, True)
    mcrt([2, 3, 5], [-1, -1, -1], 2*3*5 - 1, False)

    assert crt([656, 350], [811, 133], symmetric=True) == (-56917, 114800)

    raises(ValueError, lambda: crt([2, 3], [1]))
    raises(ValueError, lambda: crt([2], [1, 2]))
    raises(ValueError, lambda: crt([2, 3], [1], check=False))
    raises(ValueError, lambda: crt([2], [1, 2], check=False))

    # crt2 requires moduli and residues of equal length as well.
    mm, e, s = crt1([2, 3])
    raises(ValueError, lambda: crt2([2, 3], [1], mm, e, s))
    raises(ValueError, lambda: crt2([2], [1, 2], mm, e, s))
    raises(ValueError, lambda: crt2(iter([2, 3]), iter([1]), mm, e, s))
    raises(ValueError, lambda: crt2(iter([2]), iter([1, 2]), mm, e, s))
    assert crt2(iter([2, 3]), iter([1, 2]), mm, e, s) == (5, 6)

    # crt also accepts iterators while validating their lengths.
    raises(ValueError, lambda: crt(iter([2, 3]), iter([1])))
    assert crt(iter([2, 3]), iter([1, 2])) == (5, 6)


def test_crt_noncoprime_and_edge_cases():
    # Issue #31177: non-coprime moduli where gf_crt satisfies all congruences
    assert crt([3, 6], [0, 3]) == (3, 6)
    assert crt([3, 6], [0, 3], symmetric=True) == (3, 6)
    assert crt([2, 4, 3], [0, 0, 1]) == (4, 12)
    assert crt([2, 4, 3], [0, 0, 1], symmetric=True) == (4, 12)
    assert crt([4, 6], [0, 0]) == (0, 12)
    assert crt([4, 6], [0, 0], symmetric=True) == (0, 12)

    # Consistent and inconsistent non-coprime systems
    assert crt([3, 6], [2, 5]) == (5, 6)
    assert crt([3, 6], [2, 5], symmetric=True) == (-1, 6)
    assert crt([12, 6, 17], [3, 4, 2]) is None

    # Negative and non-canonical residues with non-coprime moduli
    assert crt([3, 6], [-3, -3]) == (3, 6)
    assert crt([3, 6], [9, 15]) == (3, 6)

    # Duplicate moduli (consistent and inconsistent)
    assert crt([5, 5], [0, 0]) == (0, 5)
    assert crt([5, 5], [2, 2]) == (2, 5)
    assert crt([5, 5], [2, 3]) is None

    # Empty, singleton, and large pairwise-coprime inputs
    assert crt([], []) == (0, 1)
    assert crt([5], [2]) == (2, 5)
    assert crt([5], [4], symmetric=True) == (-1, 5)
    m1, m2 = 10**20 + 7, 10**20 + 9
    r, mod = crt([m1, m2], [12345, 67890])
    assert mod == m1 * m2
    assert r % m1 == 12345 and r % m2 == 67890

    # Existing check=False semantics remain unchanged
    assert crt([12, 6, 17], [3, 4, 2], check=False) == (954, 1224)


def test_modular():
    assert solve_congruence(*list(zip([3, 4, 2], [12, 35, 17]))) == (1719, 7140)
    assert solve_congruence(*list(zip([3, 4, 2], [12, 6, 17]))) is None
    assert solve_congruence(*list(zip([3, 4, 2], [13, 7, 17]))) == (172, 1547)
    assert solve_congruence(*list(zip([-10, -3, -15], [13, 7, 17]))) == (172, 1547)
    assert solve_congruence(*list(zip([-10, -3, 1, -15], [13, 7, 7, 17]))) is None
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 7, 17]))) == (835, 1547)
    assert solve_congruence(
        *list(zip([-10, -5, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(
        *list(zip([-10, 2, 2, -15], [13, 7, 14, 17]))) == (2382, 3094)
    assert solve_congruence(*list(zip((1, 1, 2), (3, 2, 4)))) is None
    raises(
        ValueError, lambda: solve_congruence(*list(zip([3, 4, 2], [12.1, 35, 17]))))
