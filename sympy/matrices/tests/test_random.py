import random as system_random

from sympy import Matrix, diag, eye, zeros, cartes, \
    I, cos, simplify, symbols, sqrt, expand
from sympy.core.random import seed, sample
from sympy.matrices.random import super_elementary_matrix, \
    complex_to_real, invertible, regular_to_singular, diagonal_normal, \
    diagonalizable, transposition, triangular, trigonalizable, \
    isometry_normal, permutation, projection, elementary, rotation, \
    reflection, normal, nilpotent, idempotent, singular, orthogonal, unitary, \
    hermite, symmetric, square, jordan, jordan_normal

from sympy.matrices import random as _random

from sympy.testing.pytest import raises

TEST_DIMS = dict((d, tuple(cartes(range(d), range(d)))) for d in range(2, 6))
TEST_PRECISION = 7
TEST_EPSILON = 1e-7

phi, psi = symbols('phi psi')

# set fixed sympy random seed for testing purposes
sympy_seed = 12
seed(sympy_seed)

# set system random number generator with fixed seed for testing purposes
system_seed = 111
RANDOM = system_random.Random(system_seed)


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
        perm = RANDOM.sample(range(d), d)
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
        s = RANDOM.choices(spec, k=d + 4)
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
        s = RANDOM.choices(scalars, k=d + 4)
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

    m = isometry_normal(3, spec=((0, 1), 1))
    n = Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    assert _is_zeros(m - n), repr(m) + repr(n)

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

        m = invertible(d, scalars=(0.5, 2.), units=(0.1, 1.))
        assert _is_eye(m.inv() * m, TEST_PRECISION)

        m = invertible(d, scalars=(1, phi))
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

        m = orthogonal(d, spec=_random._rotation_scalars).evalf()
        assert _is_eye(m.T * m, TEST_PRECISION)
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        m = orthogonal(d, spec=(cos(phi),))
        assert _is_eye(m.T * m)


def test_unitary():
    for d in list(TEST_DIMS)[1::3]:
        i = expand(unitary(d, (1,), length=0))
        assert _is_eye(i), repr(i)

        m = expand(unitary(d))
        assert _is_eye(m.H * m, TEST_EPSILON), repr(simplify(m.H * m))
        assert _is_eye(m.H * m), repr((m.H * m).evalf())
        assert abs(abs(m.det()) - 1) < TEST_EPSILON

        z = complex(1, 2)
        m = expand(unitary(d, spec=(z / abs(z),), length=d))
        assert _is_eye(m.H * m, TEST_EPSILON), repr(simplify(m.H * m))
        assert abs(m.evalf().det()) - 1 < TEST_PRECISION

        z = 2 + 3 * I
        m = expand(unitary(d, spec=(z / abs(z),), length=d))
        assert _is_eye(m.H * m), repr(simplify(m.H * m))

        m = expand(unitary(d, spec=(sqrt(2) / 2,), length=d))
        assert _is_eye(m.H * m), repr(simplify(m.H * m))


def test_normal():
    for d in list(TEST_DIMS)[1::3]:
        spec = (2, 3, 4)
        m = normal(d, spec, length=d)
        assert all(x.is_real for x in simplify(m)), repr(simplify(m))
        assert _is_zeros(m.T * m - m * m.T), \
            repr(simplify(m.T * m - m * m.T))

        z = 1 + 3 * I
        c = expand(
            normal(d, spec=(1, 2 * I, 3), scalars=(z/abs(z),), length=d))
        assert _is_zeros(c.H * c - c * c.H), \
            repr(simplify(m.H * m - m * m.H))

        z = complex(3, 1)
        c = expand(normal(d, spec, scalars=(z/abs(z),), length=d))
        assert _is_zeros(c.H * c - c * c.H, TEST_EPSILON), \
            repr(simplify(m.H * m - m * m.H))


# === symmetric or complex adjoined matrices ===


def test_symmetric():
    for d in TEST_DIMS:
        m = symmetric(d,)
        assert m.T == m

        m = symmetric(d, scalars=(4,))
        assert m.T == m


def test_hermite():
    z = complex(1, 2)
    for d in TEST_DIMS:
        m = hermite(d, scalars=(z,))
        assert m.H == m


# === seed and sample tests ===


def test_sample():
    class Scalars(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            return self.rnd.sample(self, k)

    scalars = tuple(range(10, 25))

    for d in TEST_DIMS:
        m = diagonal_normal(d, spec=scalars)
        for i in range(d):
            assert m[i, i] in scalars

        n = diagonal_normal(d, spec=Scalars(scalars))
        for i in range(d):
            assert n[i, i] in scalars

    class CyclicPowers(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            items = self.rnd.sample(self, k)
            powers = self.rnd.choices(range(2, 6), k=k)
            return [item**p for item, p in zip(items, powers)]

    scalars = [2]

    for d in TEST_DIMS:
        n = diagonal_normal(d, spec=CyclicPowers(scalars))
        for i in range(d):
            assert n[i, i] not in scalars, n

    class Rnd(object):
        rng = system_random.Random()

        def __len__(self):
            return 1

        def sample(self, k):
            return [self.rng.random() for _ in range(k)]

    for d in TEST_DIMS:
        n = diagonal_normal(d, spec=Rnd())
        for i in range(d):
            assert 0 <= n[i, i] <= 1, repr(n)


def test_seed():

    class Scalars(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            return self.rnd.sample(self, k)


    _all_ = projection, jordan, transposition, \
          permutation, elementary, rotation, reflection, \
          diagonal_normal, jordan_normal, isometry_normal, \
          triangular, invertible, singular, \
          idempotent, nilpotent, diagonalizable, trigonalizable, \
          orthogonal, normal, \
          symmetric, hermite, square, unitary,

    _spec_ = diagonal_normal, jordan_normal, \
          diagonalizable, trigonalizable, normal

    _isometry_spec_ = isometry_normal, orthogonal, unitary

    _elementary_scalars_ = triangular, invertible, singular, \
          idempotent, nilpotent, diagonalizable, trigonalizable, \
          symmetric, hermite, square,

    _isometry_scalars_ = orthogonal, unitary, normal

    seed(101)
    seeds = sample(range(100), 3)
    for d in TEST_DIMS:
        for matrix in _all_:
            first, second = list(), list()
            for s in seeds:
                seed(s-1)
                first.append((s, expand(matrix(d))))
            for s in seeds:
                seed(s-1)
                second.append((s, expand(matrix(d))))
            for (seed_a, a), (seed_b, b) in zip(first, second):
                assert a == b,  \
                    matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _spec_:
            # 4. spec as random
            first, second = list(), list()
            spec = Scalars(range(100))
            for s in seeds:
                seed(s-1)
                spec.rnd.seed(s)
                first.append(matrix(d, spec=spec))
            for s in seeds:
                seed(s-1)
                spec.rnd.seed(s)
                second.append(matrix(d, spec=spec))
            for a, b in zip(first, second):
                assert a == b, \
                    matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _isometry_spec_:
            # 4. scalars/units as random
            first, second = list(), list()
            z = sqrt(2)/2 + sqrt(2)/2 * I
            spec = Scalars((z, -1, z*z))
            for s in seeds:
                seed(s-1)
                spec.rnd.seed(s)
                first.append(matrix(d, spec=spec))
            for s in seeds:
                seed(s-1)
                spec.rnd.seed(s)
                second.append(matrix(d, spec=spec))
            for a, b in zip(first, second):
                assert a == b, \
                    matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _elementary_scalars_:
            # 4. scalars/units as random
            first, second = list(), list()
            scalars = Scalars(range(10))
            units = Scalars(range(-5, 5))
            for s in seeds:
                seed(s-1)
                scalars.rnd.seed(s)
                first.append(matrix(d, scalars=scalars, units=units))
            for s in seeds:
                seed(s-1)
                scalars.rnd.seed(s)
                second.append(matrix(d, scalars=scalars, units=units))
            for a, b in zip(first, second):
                assert a == b, \
                    matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)

        for matrix in _isometry_scalars_:
            # 4. scalars/units as random
            first, second = list(), list()
            scalars = Scalars((sqrt(2)/2, 0))
            for s in seeds:
                seed(s-1)
                scalars.rnd.seed(s)
                first.append(matrix(d, scalars=scalars))
            for s in seeds:
                seed(s-1)
                scalars.rnd.seed(s)
                second.append(matrix(d, scalars=scalars))
            for a, b in zip(first, second):
                assert a == b, \
                    matrix.__name__ + '\n' + repr(a) + '\n' + repr(b)


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
