from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import ZZ, QQ
from sympy.polys.matrices.tests.test_domainmatrix import _check_fflu
import pytest


def test_fflu_2x2_matrix():
    A = DomainMatrix([[4, 3], [6, 3]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[4, 0], [6, -6]], (2, 2), ZZ)
    assert D == DomainMatrix([[4, 0], [0, -24]], (2, 2), ZZ)
    assert U == DomainMatrix([[4, 3], [0, -6]], (2, 2), ZZ)


def test_fflu_2x3_matrix():
    A = DomainMatrix([[1, 2, 3], [4, 5, 6]], (2, 3), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [4, -3]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, -3]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 2, 3], [0, -3, -6]], (2, 3), ZZ)


def test_fflu_3x2_matrix():
    A = DomainMatrix([[1, 2], [3, 4], [5, 6]], (3, 2), ZZ)
    P, L, D, U = A.fflu()
    assert P == DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), ZZ)
    assert L == DomainMatrix([[1, 0, 0], [3, -2, 0], [5, -4, 1]], (3, 3), ZZ)
    assert D == DomainMatrix([[1, 0], [0, -2]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 2], [0, -2], [0, 0]], (3, 2), ZZ)


def test_fflu_3x3_matrix():
    A = DomainMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]], (3, 3), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), ZZ)
    assert L == DomainMatrix([[1, 0, 0], [4, -3, 0], [7, -6, 1]], (3, 3), ZZ)
    assert D == DomainMatrix([[1, 0, 0], [0, -3, 0], [0, 0, 1]], (3, 3), ZZ)
    assert U == DomainMatrix([[1, 2, 3], [0, -3, -6], [0, 0, 0]], (3, 3), ZZ)


def test_fflu_zero_matrix():
    A = DomainMatrix.zeros((3, 3), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix.eye(3, ZZ)
    assert L == DomainMatrix.eye(3, ZZ)
    assert D == DomainMatrix.eye(3, ZZ)
    assert U == DomainMatrix.zeros((3, 3), ZZ)


def test_fflu_negative_entries():
    A = DomainMatrix([[-1, -2], [-3, -4]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[-1, 0], [-3, -2]], (2, 2), ZZ)
    assert D == DomainMatrix([[-1, 0], [0, 2]], (2, 2), ZZ)
    assert U == DomainMatrix([[-1, -2], [0, -2]], (2, 2), ZZ)


def test_fflu_mixed_signs():
    A = DomainMatrix([[1, -2], [-3, 4]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [-3, -2]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, -2]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, -2], [0, -2]], (2, 2), ZZ)


def test_fflu_upper_triangular():
    A = DomainMatrix([[1, 2, 3], [0, 4, 5], [0, 0, 6]], (3, 3), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), ZZ)
    assert L == DomainMatrix([[1, 0, 0], [0, 4, 0], [0, 0, 24]], (3, 3), ZZ)
    assert D == DomainMatrix([[1, 0, 0], [0, 4, 0], [0, 0, 96]], (3, 3), ZZ)
    assert U == DomainMatrix([[1, 2, 3], [0, 4, 5], [0, 0, 24]], (3, 3), ZZ)


def test_fflu_lower_triangular():
    A = DomainMatrix([[1, 0, 0], [2, 3, 0], [4, 5, 6]], (3, 3), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), ZZ)
    assert L == DomainMatrix([[1, 0, 0], [2, 3, 0], [4, 5, 18]], (3, 3), ZZ)
    assert D == DomainMatrix([[1, 0, 0], [0, 3, 0], [0, 0, 54]], (3, 3), ZZ)
    assert U == DomainMatrix([[1, 0, 0], [0, 3, 0], [0, 0, 18]], (3, 3), ZZ)


def test_fflu_diagonal():
    A = DomainMatrix.diag([2, 3, 4], ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix.eye(3, ZZ)
    assert L == DomainMatrix({0: {0: 2}, 1: {1: 6}, 2: {2: 24}}, (3, 3), ZZ)
    assert D == DomainMatrix({0: {0: 2}, 1: {1: 12}, 2: {2: 144}}, (3, 3), ZZ)
    assert U == DomainMatrix({0: {0: 2}, 1: {1: 6}, 2: {2: 24}}, (3, 3), ZZ)


@pytest.mark.parametrize('domain', [ZZ, QQ])
def test_fflu_different_domains(domain):
    A = DomainMatrix([[1, 2], [3, 4]], (2, 2), domain)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), domain)
    assert L == DomainMatrix([[1, 0], [3, -2]], (2, 2), domain)
    assert D == DomainMatrix([[1, 0], [0, -2]], (2, 2), domain)
    assert U == DomainMatrix([[1, 2], [0, -2]], (2, 2), domain)


def test_fflu_large_matrix():
    A = DomainMatrix([
        [1, 2, 3, 4, 5],
        [2, 3, 4, 5, 6],
        [3, 4, 5, 6, 7],
        [4, 5, 6, 7, 8],
        [5, 6, 7, 8, 9]
    ], (5, 5), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]], (5, 5), ZZ)
    assert L == DomainMatrix([
        [1, 0, 0, 0, 0],
        [2, -1, 0, 0, 0],
        [3, -2, 1, 0, 0],
        [4, -3, 0, 1, 0],
        [5, -4, 0, 0, 1]
    ], (5, 5), ZZ)
    assert D == DomainMatrix([
        [1, 0, 0, 0, 0],
        [0, -1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1]
    ], (5, 5), ZZ)
    assert U == DomainMatrix([
        [1, 2, 3, 4, 5],
        [0, -1, -2, -3, -4],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0]
    ], (5, 5), ZZ)


def test_fflu_all_zero_column():
    A = DomainMatrix([[0, 1], [0, 2]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [0, 2]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, 2]], (2, 2), ZZ)
    assert U == DomainMatrix([[0, 1], [0, 2]], (2, 2), ZZ)


def test_fflu_all_zero_row():
    A = DomainMatrix([[1, 2], [0, 0]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 2], [0, 0]], (2, 2), ZZ)


def test_fflu_empty_matrix():
    A = DomainMatrix([], (0, 0), ZZ)
    P, L, D, U = A.fflu()
    assert P == DomainMatrix([], (0, 0), ZZ)
    assert L == DomainMatrix([], (0, 0), ZZ)
    assert D == DomainMatrix([], (0, 0), ZZ)
    assert U == DomainMatrix([], (0, 0), ZZ)


def test_fflu_with_permutations():
    A = DomainMatrix([[0, 1], [1, 0]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[0, 1], [1, 0]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)


def test_fflu_zero_pivot_handling():
    A = DomainMatrix([[0, 1], [1, 0]], (2, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[0, 1], [1, 0]], (2, 2), ZZ)
    assert L == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert D == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)


def test_fflu_empty_matrix_with_cols():
    A = DomainMatrix([], (0, 2), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([], (0, 0), ZZ)
    assert L == DomainMatrix([], (0, 0), ZZ)
    assert D == DomainMatrix([], (0, 0), ZZ)
    assert U == DomainMatrix([], (0, 2), ZZ)


def test_fflu_identity():
    A = DomainMatrix.eye(3, ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix.eye(3, ZZ)
    assert L == DomainMatrix.eye(3, ZZ)
    assert D == DomainMatrix.eye(3, ZZ)
    assert U == DomainMatrix.eye(3, ZZ)


def test_fflu_empty_matrix_with_rows():
    A = DomainMatrix.zeros((2, 0), ZZ)
    P, L, D, U = A.fflu()
    assert P.to_list() == [[1, 0], [0, 1]]
    assert L.to_list() == [[1, 0], [0, 1]]
    assert D.to_list() == []
    assert U.to_list() == [[], []]


def test_fflu_single_element():
    A = DomainMatrix([[5]], (1, 1), ZZ)
    P, L, D, U = A.fflu()
    _check_fflu(A, P, L, D, U)
    assert P == DomainMatrix([[1]], (1, 1), ZZ)
    assert L == DomainMatrix([[5]], (1, 1), ZZ)
    assert D == DomainMatrix([[5]], (1, 1), ZZ)
    assert U == DomainMatrix([[5]], (1, 1), ZZ)


def test_fflu_rank_deficient():
    A = DomainMatrix([[1, 2], [2, 4], [3, 6]], (3, 2), ZZ)
    P, L, D, U = A.fflu()
    assert P == DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), ZZ)
    assert L == DomainMatrix([[1, 0, 0], [2, 1, 0], [3, 0, 1]], (3, 3), ZZ)
    assert D == DomainMatrix([[1, 0], [0, 1]], (2, 2), ZZ)
    assert U == DomainMatrix([[1, 2], [0, 0], [0, 0]], (3, 2), ZZ)
