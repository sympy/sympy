from sympy import ZZ, QQ
from sympy.polys.matrices import DM, DomainMatrix
from sympy.polys.matrices.dense import ddm_irref_den, ddm_irref
from sympy.polys.matrices.sdm import sdm_irref

import pytest


RREF_EXAMPLES = [
    (
        'zz_1',
         DM([[1, 2, 3]], ZZ),
         DM([[1, 2, 3]], ZZ),
         ZZ(1),
    ),

    (
        'zz_2',
         DomainMatrix([], (0, 0), ZZ),
         DomainMatrix([], (0, 0), ZZ),
         ZZ(1),
    ),

    (
        'zz_3',
        DM([[1, 2],
            [3, 4]], ZZ),
        DM([[-2,  0],
            [ 0, -2]], ZZ),
        ZZ(-2),
    ),

    (
        'zz_4',
        DM([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]], ZZ),
        DM([[-3,  0,  3],
            [ 0, -3, -6],
            [ 0,  0,  0]], ZZ),
        ZZ(-3),
    ),

    (
        'zz_5',
        DM([[0, 0, 0],
            [0, 0, 0],
            [1, 0, 0]], ZZ),
        DM([[1, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_6',
        DM([[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        DM([[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_7',
        DM([[1, 1, 0],
            [0, 0, 2],
            [0, 0, 0]], ZZ),
        DM([[2, 2, 0],
            [0, 0, 2],
            [0, 0, 0]], ZZ),
        ZZ(2),
    ),

    (
        'zz_8',
        DM([[2, 2, 0],
            [0, 0, 2],
            [0, 0, 0]], ZZ),
        DM([[4, 4, 0],
            [0, 0, 4],
            [0, 0, 0]], ZZ),
        ZZ(4),
    ),

    (
        'zz_9',
        DM([[2, 2, 0],
            [0, 2, 2],
            [0, 0, 2]], ZZ),
        DM([[8, 0, 0],
            [0, 8, 0],
            [0, 0, 8]], ZZ),
        ZZ(8),
    ),

    (
        'zz_10',
        DM([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [10, 11, 12]], ZZ),
        DM([[-3, 0, 3],
            [0, -3, -6],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        ZZ(-3),
    ),

    (
        'zz_11',
        DM([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [10, 11, 13]], ZZ),
        DM([[-3, 0, 0],
            [0, -3, 0],
            [0, 0, -3],
            [0, 0, 0]], ZZ),
        ZZ(-3),
    ),

    (
        'zz_12',
        DM([[1, 2, 4, 3],
            [4, 5, 10, 6],
            [7, 8, 16, 9]], ZZ),
        DM([[-3, 0, 0, 3],
            [0, -3, -6, -6],
            [0, 0, 0, 0]], ZZ),
        ZZ(-3),
    ),

    (
        'zz_13',
        DM([[1, 2, 4, 3],
            [4, 5, 10, 6],
            [7, 8, 17, 9]], ZZ),
        DM([[-3, 0, 0, 3],
            [0, -3, 0, -6],
            [0, 0, -3, 0]], ZZ),
        ZZ(-3),
    ),

    (
        'qq_1',
        DM([[(1,2), 0], [0, 2]], QQ),
        DM([[1, 0], [0, 1]], QQ),
        QQ(1),
    ),
]


def _pivots(A_rref):
    return tuple(sorted(map(min, A_rref.to_sdm().values())))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_Matrix_rref(name, A, A_rref, den):
    K = A.domain
    pivots = _pivots(A_rref)
    A = A.to_Matrix()
    A_rref = A_rref.to_Matrix() / K.to_sympy(den)
    assert A.rref() == (A_rref, pivots)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_rref_den(name, A, A_rref, den):
    pivots = _pivots(A_rref)
    assert A.rref_den() == (A_rref, den, pivots)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_rref_den(name, A, A_rref, den):
    A = A.to_ddm()
    A_rref = A_rref.to_ddm()
    pivots = _pivots(A_rref)
    assert A.rref_den() == (A_rref, den, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sdm_rref_den(name, A, A_rref, den):
    A = A.to_sdm()
    A_rref = A_rref.to_sdm()
    pivots = _pivots(A_rref)
    assert A.rref_den() == (A_rref, den, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_dense_irref_den(name, A, A_rref, den):
    A = A.to_ddm().copy()
    A_rref = A_rref.to_ddm()
    pivots = _pivots(A_rref)
    K = A.domain
    (den_found, pivots_found) = ddm_irref_den(A, K)
    assert (A, den_found, pivots_found) == (A_rref, den, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_rref(name, A, A_rref, den):
    A = A.to_field()
    A_rref = A_rref.to_field() / den
    pivots = _pivots(A_rref)
    assert A.rref() == (A_rref, pivots)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_rref(name, A, A_rref, den):
    A = A.to_field().to_ddm()
    A_rref = (A_rref.to_field() / den).to_ddm()
    pivots = _pivots(A_rref)
    assert A.rref() == (A_rref, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sdm_rref(name, A, A_rref, den):
    A = A.to_field().to_sdm()
    A_rref = (A_rref.to_field() / den).to_sdm()
    pivots = _pivots(A_rref)
    assert A.rref() == (A_rref, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dense_ddm_irref(name, A, A_rref, den):
    A = A.to_field().to_ddm().copy()
    A_rref = (A_rref.to_field() / den).to_ddm()
    pivots = _pivots(A_rref)
    K = A.domain
    pivots_found = ddm_irref(A, K)
    assert (A, pivots_found) == (A_rref, list(pivots))


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sparse_sdm_rref(name, A, A_rref, den):
    A = A.to_field().to_sdm()
    A_rref = (A_rref.to_field() / den).to_sdm()
    pivots = _pivots(A_rref)
    (A_rref_found, pivots_found, _) = sdm_irref(A)
    assert (A_rref_found, pivots_found) == (A_rref, list(pivots))
