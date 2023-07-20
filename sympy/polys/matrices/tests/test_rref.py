from sympy import ZZ, QQ, Matrix
from sympy.polys.matrices import DM, DomainMatrix
from sympy.polys.matrices.dense import ddm_irref_den, ddm_irref
from sympy.polys.matrices.ddm import DDM
from sympy.polys.matrices.sdm import SDM, sdm_irref, sdm_rref_den

import pytest


#
# The dense and sparse implementations of rref_den are ddm_irref_den and
# sdm_irref_den. These can give results that differ by some factor and also
# give different results if the order of the rows is changed. The tests below
# show all results on lowest terms as should be returned by cancel_denom.
#


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
        DM([[1, 0],
            [0, 1]], ZZ),
        ZZ(1),
    ),

    (
        'zz_4',
        DM([[1, 0],
            [3, 4]], ZZ),
        DM([[1, 0],
            [0, 1]], ZZ),
        ZZ(1),
    ),

    (
        'zz_5',
        DM([[0, 2],
            [3, 4]], ZZ),
        DM([[1, 0],
            [0, 1]], ZZ),
        ZZ(1),
    ),

    (
        'zz_6',
        DM([[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]], ZZ),
        DM([[1, 0, -1],
            [0, 1,  2],
            [0, 0,  0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_7',
        DM([[0, 0, 0],
            [0, 0, 0],
            [1, 0, 0]], ZZ),
        DM([[1, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_8',
        DM([[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        DM([[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_9',
        DM([[1, 1, 0],
            [0, 0, 2],
            [0, 0, 0]], ZZ),
        DM([[1, 1, 0],
            [0, 0, 1],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_10',
        DM([[2, 2, 0],
            [0, 0, 2],
            [0, 0, 0]], ZZ),
        DM([[1, 1, 0],
            [0, 0, 1],
            [0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_11',
        DM([[2, 2, 0],
            [0, 2, 2],
            [0, 0, 2]], ZZ),
        DM([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]], ZZ),
        ZZ(1),
    ),

    (
        'zz_12',
        DM([[ 1,  2,  3],
            [ 4,  5,  6],
            [ 7,  8,  9],
            [10, 11, 12]], ZZ),
        DM([[1,  0, -1],
            [0,  1,  2],
            [0,  0,  0],
            [0,  0,  0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_13',
        DM([[ 1,  2,  3],
            [ 4,  5,  6],
            [ 7,  8,  9],
            [10, 11, 13]], ZZ),
        DM([[ 1,  0,  0],
            [ 0,  1,  0],
            [ 0,  0,  1],
            [ 0,  0,  0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_14',
        DM([[1, 2,  4, 3],
            [4, 5, 10, 6],
            [7, 8, 16, 9]], ZZ),
        DM([[1, 0, 0, -1],
            [0, 1, 2,  2],
            [0, 0, 0,  0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_15',
        DM([[1, 2,  4, 3],
            [4, 5, 10, 6],
            [7, 8, 17, 9]], ZZ),
        DM([[1, 0, 0, -1],
            [0, 1, 0,  2],
            [0, 0, 1,  0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_16',
        DM([[1, 2, 0, 1],
            [1, 1, 9, 0]], ZZ),
        DM([[1, 0, 18, -1],
            [0, 1, -9,  1]], ZZ),
        ZZ(1),
    ),

    (
        'zz_17',
        DM([[1, 1, 1],
            [1, 2, 2]], ZZ),
        DM([[1, 0, 0],
            [0, 1, 1]], ZZ),
        ZZ(1),
    ),

    (
        # Here the sparse implementation and dense implementation give very
        # different denominators: 4061232 and -1765176.
        'zz_18',
        DM([[94, 24,  0, 27, 0],
            [79,  0,  0,  0, 0],
            [85, 16, 71, 81, 0],
            [ 0,  0, 72, 77, 0],
            [21,  0, 34,  0, 0]], ZZ),
        DM([[ 1,  0,  0,  0, 0],
            [ 0,  1,  0,  0, 0],
            [ 0,  0,  1,  0, 0],
            [ 0,  0,  0,  1, 0],
            [ 0,  0,  0,  0, 0]], ZZ),
        ZZ(1),
    ),

    (
        # Let's have a denominator that cannot be cancelled.
        'zz_19',
        DM([[1, 2, 4],
            [4, 5, 6]], ZZ),
        DM([[3, 0, -8],
            [0, 3, 10]], ZZ),
        ZZ(3),
    ),

    (
        'zz_20',
        DM([[0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 4]], ZZ),
        DM([[0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_21',
        DM([[0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 1]], ZZ),
        DM([[1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'zz_22',
        DM([[1, 1, 1, 0, 1],
            [1, 1, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [1, 1, 0, 1, 0],
            [1, 0, 0, 0, 0]], ZZ),
        DM([[1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 1],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0]], ZZ),
        ZZ(1),
    ),

    (
        'qq_1',
        DM([[(1,2), 0], [0, 2]], QQ),
        DM([[1, 0], [0, 1]], QQ),
        QQ(1),
    ),

    (
        # Standard square case
        'qq_2',
        DM([[0, 1],
            [1, 1]], QQ),
        DM([[1, 0],
            [0, 1]], QQ),
        QQ(1),
    ),

    (
        # m < n  case
        'qq_3',
        DM([[1, 2, 1],
            [3, 4, 1]], QQ),
        DM([[1, 0, -1],
            [0, 1,  1]], QQ),
        QQ(1),
    ),

    (
        # same m < n  but reversed
        'qq_4',
        DM([[3, 4, 1],
            [1, 2, 1]], QQ),
        DM([[1, 0, -1],
            [0, 1,  1]], QQ),
        QQ(1),
    ),

    (
        # m > n case
        'qq_5',
        DM([[1, 0],
            [1, 3],
            [0, 1]], QQ),
        DM([[1, 0],
            [0, 1],
            [0, 0]], QQ),
        QQ(1),
    ),

    (
        # Example with missing pivot
        'qq_6',
        DM([[1, 0, 1],
            [3, 0, 1]], QQ),
        DM([[1, 0, 0],
            [0, 0, 1]], QQ),
        QQ(1),
    ),

    (
        # Example with missing pivot and no replacement

        # This example is just enough to show a different result from the dense
        # and sparse versions of the algorithm:
        #
        #   >>> A = Matrix([[0, 1], [0, 2], [1, 0]])
        #   >>> A.to_DM().to_sparse().rref_den()[0].to_Matrix()
        #   Matrix([
        #   [1, 0],
        #   [0, 1],
        #   [0, 0]])
        #   >>> A.to_DM().to_dense().rref_den()[0].to_Matrix()
        #   Matrix([
        #   [2, 0],
        #   [0, 2],
        #   [0, 0]])
        #
        'qq_7',
        DM([[0, 1],
            [0, 2],
            [1, 0]], QQ),
        DM([[1, 0],
            [0, 1],
            [0, 0]], QQ),
        QQ(1),
    ),

]


def _to_DM(A, ans):
    """Convert the answer to DomainMatrix."""
    if isinstance(A, DomainMatrix):
        return A.to_dense()
    elif isinstance(A, Matrix):
        return A.to_DM().to_dense()

    if not (hasattr(A, 'shape') and hasattr(A, 'domain')):
        shape, domain = ans.shape, ans.domain
    else:
        shape, domain = A.shape, A.domain

    if isinstance(A, (DDM, list)):
        return DomainMatrix(list(A), shape, domain).to_dense()
    elif isinstance(A, (SDM, dict)):
        return DomainMatrix(dict(A), shape, domain).to_dense()
    else:
        assert False # pragma: no cover


def _pivots(A_rref):
    """Return the pivots from the rref of A."""
    return tuple(sorted(map(min, A_rref.to_sdm().values())))


def _check_cancel(result, rref_ans, den_ans):
    """Check the cancelled result."""
    rref, den, pivots = result
    if isinstance(rref, (DDM, SDM, list, dict)):
        assert type(pivots) is list
        pivots = tuple(pivots)
    rref = _to_DM(rref, rref_ans)
    rref2, den2 = rref.cancel_denom(den)
    assert rref2 == rref_ans
    assert den2 == den_ans
    assert pivots == _pivots(rref)


def _check_divide(result, rref_ans, den_ans):
    """Check the divided result."""
    rref, pivots = result
    if isinstance(rref, (DDM, SDM, list, dict)):
        assert type(pivots) is list
        pivots = tuple(pivots)
    rref = _to_DM(rref, rref_ans).to_field()
    rref_ans = rref_ans.to_field() / den_ans
    assert rref == rref_ans
    assert _pivots(rref) == pivots


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_Matrix_rref(name, A, A_rref, den):
    A = A.to_Matrix()
    _check_divide(A.rref(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_dense_rref(name, A, A_rref, den):
    A = A.to_field()
    _check_divide(A.rref(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_dense_rref_den(name, A, A_rref, den):
    _check_cancel(A.rref_den(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_sparse_rref(name, A, A_rref, den):
    A = A.to_field().to_sparse()
    _check_divide(A.rref(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_dm_sparse_rref_den(name, A, A_rref, den):
    A = A.to_sparse()
    _check_cancel(A.rref_den(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_rref_den(name, A, A_rref, den):
    A = A.to_ddm()
    _check_cancel(A.rref_den(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sdm_rref_den(name, A, A_rref, den):
    A = A.to_sdm()
    _check_cancel(A.rref_den(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_rref(name, A, A_rref, den):
    A = A.to_field().to_ddm()
    _check_divide(A.rref(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sdm_rref(name, A, A_rref, den):
    A = A.to_field().to_sdm()
    _check_divide(A.rref(), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_irref(name, A, A_rref, den):
    A = A.to_field().to_ddm().copy()
    pivots_found = ddm_irref(A, A.domain)
    _check_divide((A, pivots_found), A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_ddm_irref_den(name, A, A_rref, den):
    A = A.to_ddm().copy()
    (den_found, pivots_found) = ddm_irref_den(A, A.domain)
    result = (A, den_found, pivots_found)
    _check_cancel(result, A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sparse_sdm_rref(name, A, A_rref, den):
    A = A.to_field().to_sdm()
    _check_divide(sdm_irref(A)[:2], A_rref, den)


@pytest.mark.parametrize('name, A, A_rref, den', RREF_EXAMPLES)
def test_sparse_sdm_rref_den(name, A, A_rref, den):
    A = A.to_sdm().copy()
    K = A.domain
    _check_cancel(sdm_rref_den(A, K), A_rref, den)
