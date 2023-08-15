from sympy import ZZ
from sympy.polys.matrices import DM, DomainMatrix
from sympy.polys.matrices.dense import ddm_iinv
from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError
from sympy.matrices.common import NonInvertibleMatrixError

import pytest
from sympy.testing.pytest import raises


# Examples are given as adjugate matrix and determinant adj_det should match
# these exactly but inv_den only matches after cancel_denom.


INVERSE_EXAMPLES = [

    (
        'zz_1',
        DomainMatrix([], (0, 0), ZZ),
        DomainMatrix([], (0, 0), ZZ),
        ZZ(1),
    ),

    (
        'zz_2',
        DM([[2]], ZZ),
        DM([[1]], ZZ),
        ZZ(2),
    ),

    (
        'zz_3',
        DM([[2, 0],
            [0, 2]], ZZ),
        DM([[2, 0],
            [0, 2]], ZZ),
        ZZ(4),
    ),

    (
        'zz_4',
        DM([[1, 2],
            [3, 4]], ZZ),
        DM([[ 4, -2],
            [-3,  1]], ZZ),
        ZZ(-2),
    ),

    (
        'zz_5',
        DM([[2, 2, 0],
            [0, 2, 2],
            [0, 0, 2]], ZZ),
        DM([[4, -4, 4],
            [0, 4, -4],
            [0, 0,  4]], ZZ),
        ZZ(8),
    ),

    (
        'zz_6',
        DM([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]], ZZ),
        DM([[-3,   6, -3],
            [ 6, -12,  6],
            [-3,   6, -3]], ZZ),
        ZZ(0),
    ),
]


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_Matrix_inv(name, A, A_inv, den):

    def _check(**kwargs):
        if den != 0:
            assert A.inv(**kwargs) == A_inv
        else:
            raises(NonInvertibleMatrixError, lambda: A.inv(**kwargs))

    K = A.domain
    A = A.to_Matrix()
    A_inv = A_inv.to_Matrix() / K.to_sympy(den)
    _check()
    for method in ['GE', 'LU', 'ADJ', 'CH', 'LDL', 'QR']:
        _check(method=method)


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_dm_inv_den(name, A, A_inv, den):
    if den != 0:
        A_inv_f, den_f = A.inv_den()
        assert A_inv_f.cancel_denom(den_f) == A_inv.cancel_denom(den)
    else:
        raises(DMNonInvertibleMatrixError, lambda: A.inv_den())


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_dm_inv(name, A, A_inv, den):
    A = A.to_field()
    if den != 0:
        A_inv = A_inv.to_field() / den
        assert A.inv() == A_inv
    else:
        raises(DMNonInvertibleMatrixError, lambda: A.inv())


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_ddm_inv(name, A, A_inv, den):
    A = A.to_field().to_ddm()
    if den != 0:
        A_inv = (A_inv.to_field() / den).to_ddm()
        assert A.inv() == A_inv
    else:
        raises(DMNonInvertibleMatrixError, lambda: A.inv())


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_sdm_inv(name, A, A_inv, den):
    A = A.to_field().to_sdm()
    if den != 0:
        A_inv = (A_inv.to_field() / den).to_sdm()
        assert A.inv() == A_inv
    else:
        raises(DMNonInvertibleMatrixError, lambda: A.inv())


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_dense_ddm_iinv(name, A, A_inv, den):
    A = A.to_field().to_ddm().copy()
    K = A.domain
    A_result = A.copy()
    if den != 0:
        A_inv = (A_inv.to_field() / den).to_ddm()
        ddm_iinv(A_result, A, K)
        assert A_result == A_inv
    else:
        raises(DMNonInvertibleMatrixError, lambda: ddm_iinv(A_result, A, K))


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_Matrix_adjugate(name, A, A_inv, den):
    A = A.to_Matrix()
    A_inv = A_inv.to_Matrix()
    assert A.adjugate() == A_inv
    for method in ["bareiss", "berkowitz", "bird", "laplace", "lu"]:
        assert A.adjugate(method=method) == A_inv


@pytest.mark.parametrize('name, A, A_inv, den', INVERSE_EXAMPLES)
def test_dm_adj_det(name, A, A_inv, den):
    assert A.adj_det() == (A_inv, den)
