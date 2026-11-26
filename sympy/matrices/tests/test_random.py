import random as system_random

from sympy import eye, cartes, simplify, symbols, expand
from sympy.core.random import seed, sample
from sympy.matrices.random import regular_to_singular, elementary, \
    triangular, square
from sympy.matrices.random import _sample

from sympy.testing.pytest import raises

TEST_DIMS = dict((d, tuple(cartes(range(d), range(d)))) for d in range(2, 6))
TEST_PRECISION = 7

phi, psi = symbols('phi psi')

# set fixed sympy random seed for testing purposes
sympy_seed = 12
seed(sympy_seed)


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
        return all(
            abs(n[i, j]) < precision for i in range(s) for j in range(i))


# === fundamental functions ===

def test_singular_matrix():
    for d in TEST_DIMS:
        s = square(d)
        assert regular_to_singular(s) == s
        for r in range(1, d):
            assert regular_to_singular(s, r).rank() == r


# === base matrices ===


def test_elementary():
    for d in TEST_DIMS:
        m = elementary(d)
        n = eye(d)
        for x, y in zip(m.inv() * m, n):
            assert round(x, TEST_PRECISION) == y

        m = elementary(d, scalar=1)
        assert round(m.det(), TEST_PRECISION) == 1


# === compound matrices ===


def test_invertible_square():
    for d in TEST_DIMS:
        i = square(d, length=0)
        assert _is_eye(i)

        i = square(d, scalars=None, units=None)
        assert _is_eye(i)

        m = square(d)
        assert _is_eye(m.inv() * m, TEST_PRECISION)

        m = square(d, scalars=(0.5, 2.), units=(0.1, 1.))
        assert _is_eye(m.inv() * m, TEST_PRECISION)

        m = square(d, scalars=(1, phi))
        assert _is_eye(m.inv() * m)


def test_triangular():
    for d in TEST_DIMS:
        assert _is_triangular(triangular(d))
        m = triangular(d, rank=d)
        n = eye(d)
        assert m.inv() * m == n
        for r in range(1, d):
            m = triangular(d, rank=r)
            assert _is_triangular(m)
            assert m.rank() == r


def test_singular_square():
    for d in TEST_DIMS:
        for r in range(1, d):
            assert square(d, rank=r).rank() == r


# === seed and sample tests ===


def test_sample():
    class Scalars(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            return self.rnd.sample(self, k)

    scalars = tuple(range(10, 25))

    for d in TEST_DIMS:
        for x in _sample(scalars, d):
            assert x in scalars

        for x in _sample(Scalars(scalars), d):
            assert x in scalars

    class CyclicPowers(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            items = self.rnd.sample(self, k)
            powers = self.rnd.choices(range(2, 6), k=k)
            return [item**p for item, p in zip(items, powers)]

    scalars = [2]

    for d in TEST_DIMS:

        for x in _sample(CyclicPowers(scalars * 2 * d),d):
            assert x not in scalars

    class Rnd(object):
        rng = system_random.Random()

        def __len__(self):
            return 1

        def sample(self, k):
            return [self.rng.random() for _ in range(k)]

    for d in TEST_DIMS:
        for x in _sample(Rnd(), d):
            assert 0 <= x <= 1


def test_seed():

    class Scalars(list):
        """list class with build in sample method"""
        rnd = system_random.Random()

        def sample(self, k):
            return self.rnd.sample(self, k)

    _all_ = elementary, triangular, square

    _elementary_scalars_ = triangular, square

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

# === other tests ===


def test_raise():
    with raises(ValueError):
        square(3, rank=5)
