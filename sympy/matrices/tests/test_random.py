# coding=utf-8
import random

from pytest import raises

from sympy import Matrix, eye, zeros, cartes, I
from sympy.matrices.random import super_elementary_matrix
from sympy.matrices.random import *
from sympy.matrices import random as _random

TEST_DIMS = dict((d, tuple(cartes(range(d), range(d)))) for d in range(2, 6))
TEST_PRECISION = 7

random.seed(1)


def _check_triangular(m):
    s, t = m.shape
    for i in range(s):
        for j in range(i):
            assert m[i, j] == 0, m


# === fundamental constructor ===

def test_super_elementary_matrix():
    m = super_elementary_matrix(3, (0, 1))
    n = Matrix([
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 1]])
    assert m == n

    m = super_elementary_matrix(3, (0, 2), value=5)
    n = Matrix([
        [1, 0, 5],
        [0, 1, 0],
        [0, 0, 1]])
    assert m == n

    m = super_elementary_matrix(3, (2, 2), value=5)
    n = Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 5]])
    assert m == n

    m = super_elementary_matrix(3, (1, 2), 1, 2, 3, 4)
    n = Matrix([
        [1, 0, 0],
        [0, 1, 2],
        [0, -3, 4]])
    assert m == n

    for d, coords in TEST_DIMS.items():
        for i, j in coords:
            m = super_elementary_matrix(d, (i, j), value=d)
            assert m[i, j] == d

    # test ValueError
    with raises(ValueError):
        super_elementary_matrix(3, (1, 2, 3))
    with raises(ValueError):
        super_elementary_matrix(3, (1., 2.))
    with raises(ValueError):
        super_elementary_matrix(3, (1, 1), 1, 3)
    with raises(ValueError):
        super_elementary_matrix(3, (1, 2), 1, 3)


# === fundamental functions ===

def test_real_complex_matrix():
    z = 4 + 2 * I
    for d in TEST_DIMS:
        c = jordan(d, (0, d - 1), scalar=z)
        r = complex_to_real(c)
        for x in r:
            assert x.is_real
        for v in r.eigenvals(multiple=True):
            assert v in (z, z.conjugate(), 1)


def test_singular_matrix():
    for d in TEST_DIMS:
        s = invertible(d)
        assert regular_to_singular(s) == s
        for r in range(1, d):
            assert regular_to_singular(s, r).rank() == r


# === base matrices ===

def test_identity():
    for d in TEST_DIMS:
        m = identity(d)
        n = eye(d)
        assert m == n


def test_diagonal():
    for d in TEST_DIMS:
        m = diagonal_normal(d, (1,) * d)
        n = eye(d)
        assert m == n

    for d in TEST_DIMS:
        m = diagonal_normal(d, (0,) * d)
        n = zeros(d)
        assert m == n


def test_jordan():
    for d in TEST_DIMS:
        m = jordan(d, (0, d), (1,))
        n = eye(d)
        assert m.diagonal() == n.diagonal()

    for d in TEST_DIMS:
        m = jordan(d, (0, d), 0)
        n = zeros(d)
        assert m.diagonal() == n.diagonal()


def test_transposition():
    for d in TEST_DIMS:
        m = transposition(d)
        n = identity(d)
        assert m * m == n

    for d, coords in TEST_DIMS.items():
        for i, j in coords:
            m = transposition(d, (i, j))
            assert m[j, i] == 1
            assert m[i, j] == 1
            if not i == j:
                assert m[i, i] == 0
                assert m[j, j] == 0


def test_permutation():
    for d in TEST_DIMS:
        perm = random.sample(range(d), d)
        m = permutation(d, perm)
        n = identity(d).permute(perm)
        assert m == n


def test_projection():
    for d in TEST_DIMS:
        for r in range(1, d):
            m = projection(d, (0, r))
            assert sum(m) == r
            assert m * m == m
            _check_triangular(m)
            _check_triangular(m.T)


def test_elementary():
    for d in TEST_DIMS:
        m = elementary(d)
        n = identity(d)
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y

        m = elementary(d, scalar=1)
        assert round(m.det(), TEST_PRECISION) == 1


def test_rotation():
    for d in TEST_DIMS:
        m = rotation(d)
        n = identity(d)
        for x, y in zip(m.T * m, n):
            assert round(x, TEST_PRECISION) == y
        assert round(m.det(), TEST_PRECISION) == 1
    for d in TEST_DIMS:
        m = rotation(d, scalar=(1, 0))
        n = identity(d)
        assert m == n, (m, n)


def test_reflection():
    for d in TEST_DIMS:
        m = reflection(d)
        n = identity(d)
        for x, y in zip(m.T * m, n):
            assert round(x, TEST_PRECISION) == y
        assert round(m.det(), TEST_PRECISION) == -1


# === normal form matrices, i.e. defined by eigenvalues ===

def test_diagonal_normal():
    scalars = 2, 3, 5, 7
    for d in TEST_DIMS:
        rank = d - 1
        s = (0,) * (d - rank) + scalars[:rank]
        m = diagonal_normal(d, s)
        _check_triangular(m)
        _check_triangular(m.T)
        for x in m.diagonal():
            if x:
                rank -= 1
                assert x in scalars
        assert rank == 0

    spec = (1, 2, 3, 4)
    _random._TEST = True
    assert _random._TEST is True
    test = dict()
    specs = dict()
    for d in TEST_DIMS:
        s = random.choices(spec, k=d + 4)
        specs[d] = s
        m = _random.diagonal_normal(d + 4, s)
        test[d] = m
        _check_triangular(m)
        _check_triangular(m.T)
        for ev in m.eigenvals(multiple=True):
            assert ev in spec

    _random._TEST = False
    assert _random._TEST is False
    for d in TEST_DIMS:
        m = _random.diagonal_normal(d + 4, specs[d])
        assert test[d] == m


def test_jordan_normal():
    scalars = 2, 3, 5, 7
    for d in TEST_DIMS:
        rank = d - 1
        s = (0, None) * (d - rank) + scalars[:rank]
        m = jordan_normal(d, s)
        _check_triangular(m)
        assert m.rank() == rank

        for v in (0, 1, 2):
            m = jordan_normal(d, ((v, d), (v, d)))
            assert m[0, 0] == v
            for i in range(1, d):
                assert m[i, i] == v
                assert m[i - 1, i] == 1

        s = ((3, 2), (3, 2))
        m = jordan_normal(d, s)
        _check_triangular(m)

        with raises(AssertionError):
            _check_triangular(m.T)

    _random._TEST = True
    assert _random._TEST is True
    test = dict()
    specs = dict()
    for d in TEST_DIMS:
        s = random.choices(scalars, k=d + 4)
        specs[d] = s
        m = _random.jordan_normal(d + 4, s)
        test[d] = m

    _random._TEST = False
    assert _random._TEST is False
    for d in TEST_DIMS:
        m = _random.jordan_normal(d + 4, specs[d])
        assert test[d] == m


def test_isometry_normal():
    z = complex(0, 1)
    for d in TEST_DIMS:
        m = isometry_normal(d, ((z,), (z,)))
        assert abs(m.det()) == 1

    m = isometry_normal(3, ((0, -1), (0, 1)))
    assert m[2, 2] == 1

    spec = (0, -1), (1,)
    _random._TEST = True
    assert _random._TEST is True
    test = dict()
    for d in TEST_DIMS:
        test[d] = _random.isometry_normal(2 * d, spec * d)

    _random._TEST = False
    assert _random._TEST is False
    for d in TEST_DIMS:
        m = _random.isometry_normal(2 * d, spec * d)

        for x in (m - test[d]).evalf():
            assert abs(x) < TEST_PRECISION


# === compound matrices ===


def test_invertible():
    for d in TEST_DIMS:
        n = identity(d)
        i = invertible(d, None, None)
        assert n == i

    for d in TEST_DIMS:
        m = invertible(d)
        n = identity(d)
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y


def test_triangular():
    for d in TEST_DIMS:
        _check_triangular(triangular(d))
        m = triangular(d, d)
        n = identity(d)
        assert m.inv() * m == n
        for r in range(1, d):
            m = triangular(d, r)
            _check_triangular(m)
            assert m.rank() == r


def test_singular():
    for d in TEST_DIMS:
        for r in range(1, d):
            assert singular(d, r).rank() == r


def test_diagonalizable():
    for d in TEST_DIMS:
        m = diagonalizable(d, spec=(2,))
        assert m.eigenvals(multiple=True) == [2] * d
        for r in range(1, d):
            s = (0,) * (d - r) + (2,) * r
            m = diagonalizable(d, spec=s)
            assert m.rank() == r
            for ev in m.eigenvals(multiple=True):
                assert ev in (0, 2), m.eigenvals(multiple=True)


def test_trigonalizable():
    for d in TEST_DIMS:
        m = trigonalizable(d, spec=(2,))
        assert m.eigenvals(multiple=True) == [2] * d
        for r in range(1, d):
            s = (0, None) * (d - r) + (2,) * r
            m = trigonalizable(d, spec=s)
            assert m.rank() == r, (m.rank(), r)
            for ev in m.eigenvals(multiple=True):
                assert ev in (0, 2), m.eigenvals(multiple=True)


# === conjugate matrices, i.e. defined by similarity to a normal from ==


def test_idempotent():
    for d in TEST_DIMS:
        for r in range(1, d):
            m = idempotent(d, r)
            assert m * m == m
            assert m.rank() == r


def test_nilpotent():
    for d in TEST_DIMS:
        for r in range(1, d):
            m = nilpotent(d, r)
            assert m.rank() == r
            assert m ** d == zeros(d)


# === matrices conjugate by isometries ==


def test_orthogonal():
    for d in TEST_DIMS:
        n = identity(d)
        i = orthogonal(d, (1,), length=0)
        assert n == i

        m = orthogonal(d).evalf()
        for x in m.T * m - n:
            assert abs(x) < TEST_PRECISION
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        m = orthogonal(d, spec=_random._rotation_scalar_set).evalf()
        for x in m.T * m - n:
            assert abs(x) < TEST_PRECISION
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION


def test_unitary():
    # z = complex(1, 1) * _sqrt(2) / 2
    for d in TEST_DIMS:
        n = identity(d)
        i = unitary(d, (1,), length=0)
        assert n == i
        m = unitary(d).evalf()
        for x in (m.adjoint() * m - n).evalf():
            assert abs(x) < TEST_PRECISION
        det = m.evalf().det()
        assert det * det.conjugate() - 1 < TEST_PRECISION

        z = complex(0, 1)
        m = unitary(d, spec=(z, z * z)).evalf()
        for x in (m.adjoint() * m - n).evalf():
            assert abs(x) < TEST_PRECISION
        det = m.evalf().det()
        assert det * det.conjugate() - 1 < TEST_PRECISION


def test_normal():
    spec = (2, 3, 4)
    z = complex(0, 1)
    for d in TEST_DIMS:
        m = normal(d, spec).evalf()
        for x in m:
            assert x.is_real
        for x in m.T * m - m * m.T:
            assert abs(x) < TEST_PRECISION

        c = normal(d, spec, scalar_set=(z,)).evalf()
        for x in c.adjoint() * c - c * c.adjoint():
            assert abs(x) < TEST_PRECISION


# === symmetric or complex adjoined matrices ===


def test_symmetric():
    for d in TEST_DIMS:
        m = symmetric(d, scalar_set=(4,))
        assert m.T == m
        assert len(m.eigenvals(multiple=True)) == d


def test_hermite():
    z = complex(1, 2)
    for d in TEST_DIMS:
        m = hermite(d, scalar_set=(z,))
        assert m.adjoint() == m
        assert len(m.eigenvals(multiple=True)) == d


def test_raise():
    # with raises(ValueError):
    #    rotation(3, scalar=4)
    with raises(ValueError):
        jordan_normal(3, (None, 1, 2))


if __name__ == '__main__':
    import sympy

    sympy.test(__file__, verbose=True, subprocess=False)
