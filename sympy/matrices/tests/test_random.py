import random

from sympy import Matrix, diag, eye, zeros, cartes, \
    I, cos, simplify, symbols, sqrt
from sympy.matrices.random import super_elementary_matrix
from sympy.matrices.random import jordan, jordan_normal, \
    complex_to_real, invertible, regular_to_singular, diagonal_normal, \
    diagonalizable, transposition, triangular, trigonalizable, \
    isometry_normal, permutation, projection, elementary, rotation, \
    reflection, normal, nilpotent, idempotent, singular, orthogonal, unitary, \
    hermite, symmetric, rand, square

from sympy.matrices import random as _random
from sympy.testing.pytest import raises

TEST_DIMS = dict((d, tuple(cartes(range(d), range(d)))) for d in range(2, 6))
TEST_PRECISION = 7
TEST_EPSILON = 1e-7

phi, psi = symbols('phi psi')

rand.seed(1)


class Domain(list):
    rnd = random.Random(1)

    def sample(self, k):
        return self.rnd.sample(self, k)


def _is_zeros(m, precision=None):
    if precision is None:
        return all(x == 0 for x in simplify(m))
    else:
        return all(abs(x) < precision for x in m.evalf())


def _is_eye(m, precision=None):
    return _is_zeros(m - eye(*m.shape), precision)


def _is_triangular(m, precision=None):
    s, t = m.shape
    if precision is None:
        return all(m[i, j] == 0 for i in range(s) for j in range(i))
    else:
        n = m.evalf()
        return all(abs(n[i, j]) < precision for i in range(s) for j in range(i))


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

def test_diagonal():
    for d in TEST_DIMS:
        m = diagonal_normal(d, (1,) * d)
        assert _is_eye(m)

    for d in TEST_DIMS:
        m = diagonal_normal(d, (0,) * d)
        assert _is_zeros(m)


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
        n = eye(d)
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
        n = eye(d).permute(perm)
        assert m == n


def test_projection():
    for d in TEST_DIMS:
        for r in range(1, d):
            m = projection(d, (0, r))
            assert sum(m) == r
            assert m * m == m
            assert _is_triangular(m)
            assert _is_triangular(m.T)


def test_elementary():
    for d in TEST_DIMS:
        m = elementary(d)
        n = eye(d)
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y

        m = elementary(d, scalar=1)
        assert round(m.det(), TEST_PRECISION) == 1


def test_rotation():
    z = 1 + 3 * I
    for d in TEST_DIMS:
        m = rotation(d)
        assert _is_eye(m.T * m)
        assert m.det() == 1

        c = sqrt(2) / 2
        m = rotation(d, scalar=c)
        assert _is_eye(m.T * m, TEST_EPSILON), repr(m.T * m)
        assert (m.det() - 1) < TEST_EPSILON, m.det()

        z = complex(3, 1)
        m = rotation(d, scalar=z / abs(z))
        assert _is_eye(m.T * m, TEST_EPSILON), repr(m.T * m)
        assert (m.det() - 1) < TEST_EPSILON, m.det()

        z = 1 + 3 * I
        m = rotation(d, scalar=z / abs(z))
        assert _is_eye(m.T * m)
        assert m.det() == 1

        m = rotation(d, scalar=(1, 0))
        assert _is_eye(m)


def test_reflection():
    for d in TEST_DIMS:
        m = reflection(d)
        assert _is_eye(m.T * m)
        assert m.det() == -1

        c = sqrt(2) / 2
        m = reflection(d, scalar=c)
        assert _is_eye(m.T * m)
        assert m.det() == -1

        z = complex(3, 1)
        m = reflection(d, scalar=z / abs(z))
        assert _is_eye(m.T * m, TEST_EPSILON), repr(m.T * m)
        assert (m.det() + 1) < TEST_EPSILON, m.det()

        z = 1 + 3 * I
        m = reflection(d, scalar=z / abs(z))
        assert _is_eye(m.T * m)
        assert m.det() == -1


# === normal form matrices, i.e. defined by eigenvalues ===

def test_diagonal_normal():
    scalars = 2, 3, 5, 7
    for d in TEST_DIMS:
        rank = d - 1
        s = (0,) * (d - rank) + scalars[:rank]
        m = diagonal_normal(d, s)
        assert _is_triangular(m)
        assert _is_triangular(m.T)
        for x in m.diagonal():
            if x:
                rank -= 1
                assert x in scalars
        assert rank == 0

    spec = (1, 2, 3, 4)
    _random._ALT = True
    assert _random._ALT is True
    test = dict()
    specs = dict()
    for d in TEST_DIMS:
        s = random.choices(spec, k=d + 4)
        specs[d] = s
        m = _random.diagonal_normal(d + 4, s)
        test[d] = m
        assert _is_triangular(m)
        assert _is_triangular(m.T)
        for ev in m.eigenvals(multiple=True):
            assert ev in spec

    _random._ALT = False
    assert _random._ALT is False
    for d in TEST_DIMS:
        m = _random.diagonal_normal(d + 4, specs[d])
        assert test[d] == m


def test_jordan_normal():
    scalars = 2, 3, 5, 7
    for d in TEST_DIMS:
        rank = d - 1
        s = (0, None) * (d - rank) + scalars[:rank]
        m = jordan_normal(d, s)
        assert _is_triangular(m)
        assert m.rank() == rank

        for v in (0, 1, 2):
            m = jordan_normal(d, ((v, d), (v, d)))
            assert m[0, 0] == v
            for i in range(1, d):
                assert m[i, i] == v
                assert m[i - 1, i] == 1

        s = ((3, 2), (3, 2))
        m = jordan_normal(d, s)
        assert _is_triangular(m)

        with raises(AssertionError):
            assert _is_triangular(m.T)

    _random._ALT = True
    assert _random._ALT is True
    test = dict()
    specs = dict()
    for d in TEST_DIMS:
        s = random.choices(scalars, k=d + 4)
        specs[d] = s
        m = _random.jordan_normal(d + 4, s)
        test[d] = m

    _random._ALT = False
    assert _random._ALT is False
    for d in TEST_DIMS:
        m = _random.jordan_normal(d + 4, specs[d])
        assert test[d] == m


def test_isometry_normal():
    spec = -1, 1, -1
    m = isometry_normal(3, spec=spec)
    assert _is_zeros(m - diag(*spec)), repr(m)
    assert _is_eye(m * m), repr(m * m)

    m = isometry_normal(4, spec=(1,))
    assert _is_eye(m), repr(m)

    m = isometry_normal(4, spec=((1, 0), (1, 0)))
    assert _is_eye(m), repr(m)

    m = isometry_normal(3, spec=((0, -1), (0, 1)))
    n = Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    assert _is_zeros(m - n), repr(m)

    for d in TEST_DIMS:
        m = isometry_normal(d, spec=(I,))
        assert abs(m.det()) == 1, m.det()
        m = isometry_normal(d, spec=(I,), real=False)
        assert m.det() == I ** d, repr(m)

        m = isometry_normal(d, spec=(sqrt(2) / 2,))
        assert abs(m.det()) == 1, m.det()

        z = 2 + 3 * I
        m = isometry_normal(d, spec=(z / abs(z),))
        assert abs(m.det()) == 1, m.det()
        m = isometry_normal(d, spec=(z / abs(z),), real=False)
        assert _is_zeros(m - diag(*(d * [z / abs(z)]))), repr(m)

        m = isometry_normal(d, spec=(complex(0, 1),))
        assert abs(m.det()) == 1, m.det()

        z = complex(3, 2)
        m = isometry_normal(d, spec=(z / abs(z),))
        assert abs(abs(m.det()) - 1) < TEST_EPSILON, m.det()

    spec = (0, -1), (1,)
    _random._ALT = True
    assert _random._ALT is True
    test = dict()
    for d in TEST_DIMS:
        test[d] = _random.isometry_normal(2 * d, spec * d)

    _random._ALT = False
    assert _random._ALT is False
    for d in TEST_DIMS:
        m = _random.isometry_normal(2 * d, spec * d)

        for x in (m - test[d]).evalf():
            assert abs(x) < TEST_PRECISION


# === compound matrices ===


def test_invertible():
    for d in TEST_DIMS:
        i = invertible(d, length=0)
        assert _is_eye(i)

        i = invertible(d, None, None)
        assert _is_eye(i)

        m = invertible(d)
        assert _is_eye(m.inv() * m, TEST_PRECISION)

        m = invertible(d, domain=(0.5, 2.), units=(0.1, 1.))
        assert _is_eye(m.inv() * m, TEST_PRECISION)

        m = invertible(d, domain=(1, phi))
        assert _is_eye(m.inv() * m)


def test_triangular():
    for d in TEST_DIMS:
        assert _is_triangular(triangular(d))
        m = triangular(d, d)
        n = eye(d)
        assert m.inv() * m == n
        for r in range(1, d):
            m = triangular(d, r)
            assert _is_triangular(m)
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
            assert _is_zeros(m ** d)


# === matrices conjugate by isometries ==


def test_orthogonal():
    for d in TEST_DIMS:
        i = orthogonal(d, (1,), length=0)
        assert _is_eye(i), repr(i)

        m = orthogonal(d).evalf()
        assert _is_eye(m.T * m, TEST_PRECISION)
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        m = orthogonal(d, spec=_random._rotation_domain).evalf()
        assert _is_eye(m.T * m, TEST_PRECISION)
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        m = orthogonal(d, spec=(cos(phi),))
        assert _is_eye(m.T * m)


def test_unitary():
    for d in list(TEST_DIMS)[1::3]:
        i = unitary(d, (1,), length=0)
        assert _is_eye(i), repr(i)

        m = unitary(d)
        assert _is_eye(m.H * m, TEST_EPSILON), repr(simplify(m.H * m))
        assert _is_eye(m.H * m), repr((m.H * m).evalf())
        assert abs(abs(m.det()) - 1) < TEST_EPSILON

        z = complex(1, 2)
        m = unitary(d, spec=(z / abs(z),), length=d)
        assert _is_eye(m.H * m, TEST_EPSILON), repr(simplify(m.H * m))
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        z = 2 + 3 * I
        m = unitary(d, spec=(z / abs(z),), length=d)
        assert _is_eye(m.H * m), repr(simplify(m.H * m))

        m = unitary(d, spec=(sqrt(2) / 2,), length=d)
        assert _is_eye(m.H * m), repr(simplify(m.H * m))


def test_normal():
    for d in list(TEST_DIMS)[1::3]:
        spec = (2, 3, 4)
        m = normal(d, spec, length=d)
        assert all(x.is_real for x in simplify(m)), repr(simplify(m))
        assert _is_zeros(m.T * m - m * m.T), \
            repr(simplify(m.T * m - m * m.T))

        z = 1 + 3 * I
        c = normal(d, spec=(1, 2 * I, 3), domain=(z/abs(z),), length=d)
        assert _is_zeros(c.H * c - c * c.H), \
            repr(simplify(m.H * m - m * m.H))

        z = complex(3, 1)
        c = normal(d, spec, domain=(z/abs(z),), length=d)
        assert _is_zeros(c.H * c - c * c.H, TEST_EPSILON), \
            repr(simplify(m.H * m - m * m.H))


# === symmetric or complex adjoined matrices ===


def test_symmetric():
    for d in TEST_DIMS:
        m = symmetric(d,)
        assert m.T == m

        m = symmetric(d, domain=(4,))
        assert m.T == m


def test_hermite():
    z = complex(1, 2)
    for d in TEST_DIMS:
        m = hermite(d, domain=(z,))
        assert m.H == m


# === other tests ===


def test_raise():
    with raises(ValueError):
        rotation(3, scalar=4)
    with raises(ValueError):
        isometry_normal(3, spec=(4,))
    with raises(ValueError):
        orthogonal(3, spec=(4,))
    with raises(ValueError):
        unitary(3, spec=(complex(1, 1),))
    with raises(ValueError):
        jordan_normal(3, (None, 1, 2))


def test_sample():
    seed = 11
    rnd = random.Random(seed)
    domain = tuple(range(10, 25))
    for d in TEST_DIMS:
        m = diagonal_normal(d, spec=Domain(domain))
        for i in range(d):
            assert m[i, i] in domain

        m = diagonal_normal(d, spec=domain, seed=rnd)
        for i in range(d):
            assert m[i, i] in domain

        m = diagonal_normal(d, spec=domain, seed=seed)
        for i in range(d):
            assert m[i, i] in domain


def test_seed():
    _all_ = projection, jordan, transposition, \
          permutation, elementary, rotation, reflection, \
          diagonal_normal, jordan_normal, isometry_normal, \
          triangular, invertible, singular, \
          idempotent, nilpotent, diagonalizable, trigonalizable, \
          orthogonal, unitary, normal, \
          symmetric, hermite, square,

    _spec_ = diagonal_normal, jordan_normal, \
          diagonalizable, trigonalizable, normal

    _isometry_spec_ = isometry_normal, orthogonal, unitary

    _elementary_domain_ = triangular, invertible, singular, \
          idempotent, nilpotent, diagonalizable, trigonalizable, \
          symmetric, hermite, square,

    _isometry_domain_ = orthogonal, unitary, normal

    seed = 101
    rnd = random.Random(seed)
    seeds = rnd.sample(range(100), 3)
    for d in TEST_DIMS:
        for matrix in _all_:
            # 1. standard rand
            first, second = list(), list()
            for s in seeds:
                rand.seed(s)
                first.append(matrix(d))
            for s in seeds:
                rand.seed(s)
                second.append(matrix(d))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)
            # for a, b in zip(first, reversed(second)):
            #     assert a != b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _all_:
            # 2. direct seed
            first, second = list(), list()
            for s in seeds:
                rand.seed(s)
                first.append(matrix(d, seed=s))
            for s in seeds:
                second.append(matrix(d, seed=s))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _all_:
            # 3. direct seed as random
            first, second = list(), list()
            rnd = random.Random()
            for s in seeds:
                rnd.seed(s)
                first.append(matrix(d, seed=rnd))
            for s in seeds:
                rnd.seed(s)
                second.append(matrix(d, seed=rnd))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _spec_:
            # 4. spec as random
            first, second = list(), list()
            spec = Domain(range(100))
            for s in seeds:
                spec.rnd.seed(s)
                first.append(matrix(d, spec=spec, seed=spec.rnd))
            for s in seeds:
                spec.rnd.seed(s)
                second.append(matrix(d, spec=spec, seed=spec.rnd))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _isometry_spec_:
            # 4. domain/units as random
            first, second = list(), list()
            spec = Domain((sqrt(2)/2, 0))
            for s in seeds:
                spec.rnd.seed(s)
                first.append(matrix(d, spec=spec, seed=spec.rnd))
            for s in seeds:
                spec.rnd.seed(s)
                second.append(matrix(d, spec=spec, seed=spec.rnd))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _elementary_domain_:
            # 4. domain/units as random
            first, second = list(), list()
            domain = Domain(range(10))
            units = Domain(range(-5, 5))
            for s in seeds:
                domain.rnd.seed(s)
                first.append(matrix(d, domain=domain, units=units, seed=domain.rnd))
            for s in seeds:
                domain.rnd.seed(s)
                second.append(matrix(d, domain=domain, units=units, seed=domain.rnd))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _isometry_domain_:
            # 4. domain/units as random
            first, second = list(), list()
            domain = Domain((sqrt(2)/2, 0))
            for s in seeds:
                domain.rnd.seed(s)
                first.append(matrix(d, domain=domain, seed=domain.rnd))
            for s in seeds:
                domain.rnd.seed(s)
                second.append(matrix(d, domain=domain, seed=domain.rnd))
            for a, b in zip(first, second):
                assert a == b, matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)
