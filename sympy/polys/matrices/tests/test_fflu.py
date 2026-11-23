from __future__ import annotations
from sympy.polys.matrices import DomainMatrix, DM
from sympy.polys.domains import ZZ, QQ
from sympy.external.gmpy import GROUND_TYPES
import pytest


FFLU_EXAMPLES = [
    ('zz_2x3', DM([[1, 2, 3], [4, 5, 6]], ZZ)),
    ('zz_2x2', DM([[4, 3], [6, 3]], ZZ)),
    ('zz_3x2', DM([[1, 2], [3, 4], [5, 6]], ZZ)),
    ('zz_3x3', DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]], ZZ)),
    ('zz_zero', DM([[0, 0, 0], [0, 0, 0], [0, 0, 0]], ZZ)),
    ('zz_empty', DM([], ZZ)),
    ('zz_empty_0x2', DomainMatrix([], (0, 2), ZZ)),
    ('zz_empty_2x0', DomainMatrix([[], []], (2, 0), ZZ)),
    ('zz_negative', DM([[-1, -2], [-3, -4]], ZZ)),
    ('zz_mixed_signs', DM([[1, -2], [-3, 4]], ZZ)),
    ('zz_upper_triangular', DM([[1, 2, 3], [0, 4, 5], [0, 0, 6]], ZZ)),
    ('zz_lower_triangular', DM([[1, 0, 0], [2, 3, 0], [4, 5, 6]], ZZ)),
    ('zz_diagonal', DM([[2, 0, 0], [0, 3, 0], [0, 0, 4]], ZZ)),
    ('rank_deficient_3x3', DM([[1, 2, 3], [2, 4, 6], [3, 6, 9]], ZZ)),
    ('zz_1x1', DM([[5]], ZZ)),
    ('zz_nx1_2rows', DM([[81], [54]], ZZ)),
    ('zz_nx2_3rows', DM([[2, 7], [7, 45], [25, 84]], ZZ)),
    ('zz_1x2', DM([[0, 28]], ZZ)),
    ('zz_nx3_4rows', DM([[84, 30, 9], [20, 59, 13], [53, 46, 81], [63, 48, 29]], ZZ)),
    ('fflu_row_swap', DM([[0, 1, 2], [3, 4, 5], [6, 7, 8]], ZZ)),
    ('pivots_beyond_row_count', DM([[0, 0, 1, 1], [0, 0, 1, 2]], ZZ)),
    ('wide_with_zero_row', DM([[0, 0, 0, 2], [0, 0, 0, 0], [2, -1, 2, 0]], QQ)),
    ('skipped_first_column', DM([[0, 2, 2, 2], [0, 2, 2, 0], [0, 0, 1, 0]], ZZ)),
    ('skipped_column_with_row_swap', DM([[0, 1, 0, 3], [0, 0, 0, 1],
                                         [0, 0, 1, 0], [0, 0, 0, 0]], ZZ)),
]


def _check_fflu_result(A):
    if not isinstance(A, DomainMatrix):
        A = DomainMatrix.from_rep(A)
    P, L, D, U = A.fflu()
    m, n = A.shape
    assert P.shape == (m, m)
    assert L.shape == (m, m)
    assert D.shape == (m, m)
    assert U.shape == (m, n)
    assert L.is_lower
    assert D.is_diagonal
    assert U.is_upper
    A, P, L, D, U = (matrix.to_sparse() for matrix in (A, P, L, D, U))
    di, d = D.inv_den()
    assert P.matmul(A).rmul(d) == L.matmul(di).matmul(U)


@pytest.mark.parametrize('name, A', FFLU_EXAMPLES)
def test_dm_dense_fflu(name, A):
    A = A.to_dense()
    _check_fflu_result(A)


@pytest.mark.parametrize('name, A', FFLU_EXAMPLES)
def test_dm_sparse_fflu(name, A):
    A = A.to_sparse()
    _check_fflu_result(A)


@pytest.mark.parametrize('name, A', FFLU_EXAMPLES)
def test_ddm_fflu(name, A):
    A = A.to_ddm()
    _check_fflu_result(A)


@pytest.mark.parametrize('name, A', FFLU_EXAMPLES)
def test_sdm_fflu(name, A):
    A = A.to_sdm()
    _check_fflu_result(A)


@pytest.mark.parametrize('name, A', FFLU_EXAMPLES)
@pytest.mark.skipif(GROUND_TYPES != 'flint', reason="DFM requires FLINT ground types")
def test_dfm_fflu(name, A):
    A = A.to_dfm()
    _check_fflu_result(A)


def test_fflu_empty_matrix():
    A = DomainMatrix([], (0, 0), ZZ)
    P, L, D, U = A.fflu()
    assert P.shape == (0, 0)
    assert L.shape == (0, 0)
    assert D.shape == (0, 0)
    assert U.shape == (0, 0)


def test_fflu_properties():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    assert P.shape == (2, 2)
    assert L.shape == (2, 2)
    assert D.shape == (2, 2)
    assert U.shape == (2, 2)
    assert L.is_lower
    assert U.is_upper
    assert D.is_diagonal
    di, d = D.inv_den()
    assert P.matmul(A).rmul(d) == L.matmul(di).matmul(U)


def test_fflu_rank_deficient():
    A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(2), ZZ(4)]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    assert P.shape == (2, 2)
    assert L.shape == (2, 2)
    assert D.shape == (2, 2)
    assert U.shape == (2, 2)
    assert U.getitem_sympy(1, 1) == 0
