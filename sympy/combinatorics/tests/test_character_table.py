from __future__ import annotations

from sympy import Matrix, sqrt, I
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, AlternatingGroup
from sympy.combinatorics.galois import (
    S1TransitiveSubgroups,
    S2TransitiveSubgroups,
    S3TransitiveSubgroups,
    S4TransitiveSubgroups,
    S5TransitiveSubgroups,
    S6TransitiveSubgroups
)
from sympy.combinatorics.character_table import CharacterTable, dixon_character_table

import pytest

def _cmp_tbl(tbl1, tbl2):
    # compare character tables numerically up to permutations of rows
    tbl1 = Matrix(tbl1).n(15).tolist()
    tbl2 = Matrix(tbl2).n(15).tolist()
    return {tuple(row) for row in tbl1} == {tuple(row) for row in tbl2}


def test_dixon_character_table():
    # trivial group
    G = PermutationGroup(Permutation(size=2))
    tbl = dixon_character_table(G.conjugacy_classes())
    assert tbl._rep.domain.is_ZZ
    assert _cmp_tbl(tbl, [[1]])

    # S4
    G = SymmetricGroup(4)
    cc = [
        Permutation([[3]]),
        Permutation([[0,1],[3]]),
        Permutation([[0,1],[2,3]]),
        Permutation([[0,1,2],[3]]),
        Permutation([[0,1,2,3]])
    ]
    cc = [G.conjugacy_class(c) for c in cc]
    tbl = dixon_character_table(cc)
    expected = [
        [1, 1, 1, 1, 1],
        [1, -1, 1, 1, -1],
        [2, 0, 2, -1, 0],
        [3, 1, -1, 0, -1],
        [3, -1, -1, 0, 1],
    ]
    assert tbl._rep.domain.is_ZZ
    assert _cmp_tbl(tbl, expected)

    # A4
    G = AlternatingGroup(4)
    cc = [
        Permutation([[3]]),
        Permutation([[1,2,3]]),
        Permutation([[1,3,2]]),
        Permutation([[0,1],[2,3]])
    ]
    cc = [G.conjugacy_class(c) for c in cc]
    tbl = dixon_character_table(cc)
    dom = tbl._rep.domain
    assert dom.is_CyclotomicField and dom.zeta_order == 3
    w = (-1 + sqrt(3)*I)/2
    expected = [
        [1, 1, 1, 1],
        [1, w, w**2, 1],
        [1, w**2, w, 1],
        [3, 0, 0, -1],
    ]
    assert _cmp_tbl(tbl, expected)

    # C₅ ⋊ C₄
    G = PermutationGroup(
        Permutation([[0,1,3,4,2]]),
        Permutation([[1,4,2,3]])
    )
    cc = [
        Permutation([[4]]),
        Permutation([[1,2],[3,4]]),
        Permutation([[1,3,2,4]]),
        Permutation([[1,4,2,3]]),
        Permutation([[0,1,3,4,2]])
    ]
    cc = [G.conjugacy_class(c) for c in cc]
    tbl = dixon_character_table(cc)
    dom = tbl._rep.domain
    assert dom.is_CyclotomicField and dom.zeta_order == 4
    expected =  [
        [1, 1, 1, 1, 1],
        [1, 1, -1, -1, 1],
        [1, -1, -I, I, 1],
        [1, -1, I, -I, 1],
        [4, 0, 0, 0, -1],
    ]
    assert _cmp_tbl(tbl, expected)


GROUPS = [
    member
    for enum_cls in (
        S1TransitiveSubgroups,
        S2TransitiveSubgroups,
        S3TransitiveSubgroups,
        S4TransitiveSubgroups,
        S5TransitiveSubgroups,
        S6TransitiveSubgroups,
    )
    for member in enum_cls
]

@pytest.mark.parametrize("group", GROUPS, ids=lambda g: g.value)
def test_generic_character_table(group):
    G = group.get_perm_group()
    tbl = CharacterTable.from_perm_group(G)
    cc = G.conjugacy_classes()
    reps = tbl.conjugacy_class_reps()
    assert all(r in c for r, c in zip(reps, cc))

    order = G.order()
    assert tbl.shape[0] == tbl.shape[1] == len(cc)
    assert all(v == 1 for v in tbl[0, :])
    assert tbl[:, 0].dot(tbl[:, 0]) == order

    tblh = tbl.H
    x = tblh * tbl
    assert x.is_diagonal()
    assert list(x.diagonal()) == [order//len(c) for c in cc]
