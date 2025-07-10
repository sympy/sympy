from sympy.polys.domains import ZZ, QQ
from sympy.polys.series.powerseriesring import PowerSeriesRing
from sympy.polys.series.python_powerseriesring import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)
from sympy.external.gmpy import GROUND_TYPES
import pytest
from sympy.testing.pytest import raises

# Rings
Ring_ZZ: list[type[PowerSeriesRing]] = [PythonPowerSeriesRingZZ]
Ring_QQ: list[type[PowerSeriesRing]] = [PythonPowerSeriesRingQQ]

if GROUND_TYPES == "flint":
    from sympy.polys.series.flint_powerseriesring import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )

    Ring_ZZ.append(FlintPowerSeriesRingZZ)
    Ring_QQ.append(FlintPowerSeriesRingQQ)


def same(R, s, coeffs_l, prec):
    """
    Helper function to assert equality of two power series.
    """
    domain = R.domain
    domain_coeffs = []
    for coeff in coeffs_l:
        if isinstance(coeff, tuple):
            domain_coeffs.append(domain(*coeff))
        else:
            domain_coeffs.append(domain(coeff))

    return R.equal_repr(s, R.from_list(domain_coeffs, prec))


@pytest.fixture(params=Ring_ZZ+Ring_QQ)
def rd_int(request):
    r = request.param
    R = r()
    x = R.gen
    one = R.one

    R3 = r(3)
    x3 = R3.gen
    one3 = R3.one

    R10 = r(10)
    x10 = R10.gen
    one10 = R10.one
    return R, x, one, R3, x3, one3, R10, x10, one10, r


@pytest.fixture(params=Ring_QQ)
def rd_rational(request):
    r = request.param
    R = r()
    x = R.gen
    one = R.one

    R3 = r(3)
    x3 = R3.gen
    one3 = R3.one

    R10 = r(10)
    x10 = R10.gen
    one10 = R10.one
    return R, x, one, R3, x3, one3, R10, x10, one10, r


def test_equal(rd_int):
    R, *_ = rd_int
    assert R.equal(R.from_list([1, 2, 3]), R.from_list([1, 2, 3])) is True
    assert R.equal(R.from_list([1, 21, 3]), R.from_list([1, 2, 3])) is False
    assert R.equal(R.from_list([1, 2, 3], 3), R.from_list([1, 2, 3], 3)) is None
    assert R.equal(R.from_list([1, 2, 13], 3), R.from_list([1, 2, 3], 3)) is False
    assert R.equal(R.from_list([1, 2, 3], 3), R.from_list([1, 2, 3], 10)) is None
    assert R.equal(R.from_list([1, 21, 3], 3), R.from_list([1, 2, 3], 10)) is False
    assert R.equal(R.from_list([1, 2, 3], 3), R.from_list([1, 2, 3, 4, 5], 10)) is None
    assert R.equal(R.from_list([1, 2, 3, 4], None), R.from_list([1, 2, 3], 2)) is None
    assert R.equal(R.from_list([1, 2, 3, 4], None), R.from_list([1, 1, 3], 2)) is False
    assert R.equal(R.from_list([1, 2], None), R.from_list([1, 2, 3], 3)) is False


def test_basics(rd_int):
    R, *_, r = rd_int
    R0 = r(0)
    gen0 = R0.gen

    assert R == r(6)
    assert hash(R) == hash(r(6))
    assert R.to_list(R.gen) == [0, 1], None
    assert R0.equal_repr(R0.zero, R0.one)
    assert R0.pretty(gen0) == "O(x**0)"
    assert same(R0, gen0, [], 0)
    assert same(R0, R0.multiply(gen0, gen0), [], 0)
    assert same(R, R.add(R.from_list([2, 4, 5], 3), R.from_list([5], 2)), [7, 4], 2)


def test_positive(rd_int):
    R, *_ = rd_int
    assert same(
        R, R.positive(R.from_list([1, 2, 3, 4, 5, 6, 7], None)), [1, 2, 3, 4, 5, 6], 6
    )


def test_negative(rd_int):
    R, *_ = rd_int
    assert same(R, R.negative(R.gen), [0, -1], None)
    assert same(
        R,
        R.negative(R.from_list([1, 2, 3, 4, 5, 6, 7], None)),
        [-1, -2, -3, -4, -5, -6],
        6,
    )


def test_int_add(rd_int):
    _, _, _, R3, x3, one3, *_ = rd_int
    assert same(R3, R3.add(x3, one3), [1, 1], None)
    assert same(R3, R3.add(one3, R3.pow_int(x3, 4)), [1], 3)


def test_rational_add(rd_rational):
    R, _, one, *_ = rd_rational
    assert same(R, R.add(R.negative(one), one), [], None)


def test_int_subtract(rd_int):
    R, x, _, R3, x3, one3, *_ = rd_int
    assert same(R, R.subtract(R.multiply(x, x), x), [0, -1, 1], None)
    assert same(R3, R3.subtract(R3.pow_int(R3.add(one3, x3), 4), x3), [1, 3, 6], 3)


def test_rational_subtract(rd_rational):
    R, x, *_ = rd_rational
    assert same(R, R.subtract(R.multiply(x, x), x), [(0, 1), (-1, 1), (1, 1)], None)


def test_int_multiply(rd_int):
    R, x, _, R3, x3, one3, R10, x10, one10, _ = rd_int
    assert same(R, R.multiply(x, x), [0, 0, 1], None)
    assert same(
        R3, R3.multiply(R3.square(R3.add(x3, one3)), R3.add(x3, one3)), [1, 3, 3], 3
    )
    assert same(
        R10, R10.multiply(R10.add(x10, one10), R10.add(x10, one10)), [1, 2, 1], None
    )


def test_rational_multiply(rd_rational):
    R, x, one, R3, x3, one3, *_ = rd_rational
    assert same(R, R.multiply(x, x), [(0, 1), (0, 1), (1, 1)], None)
    assert same(
        R, R.multiply(R.add(x, one), R.add(x, one)), [(1, 1), (2, 1), (1, 1)], None
    )
    assert same(
        R,
        R.multiply(R.subtract(x, one), R.add(x, one)),
        [(-1, 1), (0, 1), (1, 1)],
        None,
    )
    assert same(
        R3,
        R3.multiply(R3.square(R3.add(x3, one3)), R3.add(x3, one3)),
        [(1, 1), (3, 1), (3, 1)],
        3,
    )


def test_int_multiply_ground(rd_int):
    R, x, one, *_ = rd_int
    assert same(R, R.multiply_ground(one, 1), [1], None)
    assert same(R, R.multiply_ground(x, -1), [0, -1], None)
    assert same(R, R.multiply_ground(x, 0), [], None)
    assert same(R, R.multiply_ground(R.square(x), 1), [0, 0, 1], None)
    assert same(R, R.multiply_ground(R.add(x, one), ZZ(3)), [3, 3], None)


def test_rational_multiply_ground(rd_rational):
    R, x, *_ = rd_rational
    assert same(R, R.multiply_ground(x, QQ(1, 2)), [(0, 1), (1, 2)], None)
    assert same(
        R,
        R.multiply(R.multiply_ground(x, QQ(1, 2)), R.multiply_ground(x, QQ(1, 3))),
        [(0, 1), (0, 1), (1, 6)],
        None,
    )


def test_int_square(rd_int):
    R, x, one, R3, x3, one3, *_ = rd_int
    assert not same(R, R.square(R.add(x, one)), [1, 1, 1], None)
    assert same(R3, R3.square(R3.multiply(x3, R3.add(x3, one3))), [0, 0, 1], 3)


def test_rational_square(rd_rational):
    R, x, one, *_, R10, x10, one10, _ = rd_rational
    assert not same(R, R.square(R.add(x, one)), [(1, 1), (1, 1), (1, 1)], None)
    assert same(
        R10, R10.square(R10.subtract(x10, one10)), [(1, 1), (-2, 1), (1, 1)], None
    )


def test_int_pow(rd_int):
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_int
    assert same(R, R.pow_int(x, 0), [1], None)
    assert same(R, R.pow_int(R.add(x, one), 6), [1, 6, 15, 20, 15, 6], 6)
    assert same(R, R.pow_int(x, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(x3, one3), 5), [1, 5, 10], 3)
    assert same(R10, R10.pow_int(x10, 7), [0, 0, 0, 0, 0, 0, 0, 1], None)
    assert same(
        R10,
        R10.pow_int(R10.add(x10, one10), 12),
        [1, 12, 66, 220, 495, 792, 924, 792, 495, 220],
        10,
    )


def test_rational_pow(rd_rational):
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_rational
    assert same(R, R.pow_int(x, 3), [(0, 1), (0, 1), (0, 1), (1, 1)], None)
    assert same(
        R,
        R.pow_int(R.add(x, one), 6),
        [(1, 1), (6, 1), (15, 1), (20, 1), (15, 1), (6, 1)],
        6,
    )
    assert same(R, R.pow_int(x, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(x3, one3), 5), [(1, 1), (5, 1), (10, 1)], 3)
    assert same(
        R10,
        R10.pow_int(R10.add(x10, one10), 12),
        [
            (1, 1),
            (12, 1),
            (66, 1),
            (220, 1),
            (495, 1),
            (792, 1),
            (924, 1),
            (792, 1),
            (495, 1),
            (220, 1),
        ],
        10,
    )


def test_truncate(rd_int):
    R, x, one, *_ = rd_int
    assert same(R, R.truncate(R.pow_int(x, 3), 4), [0, 0, 0, 1], None)
    assert same(R, R.truncate(R.pow_int(R.add(x, one), 5), 3), [1, 5, 10], 3)


def test_int_differentiate(rd_int):
    R, x, _, R3, x3, one3, R10, x10, _, r = rd_int
    assert same(R, R.differentiate(R.pow_int(x, 3)), [0, 0, 3], None)
    assert same(R, R.differentiate(R.add(R.multiply(x, x), x)), [1, 2], None)
    assert same(R3, R3.differentiate(R3.multiply(x3, x3)), [0, 2], None)
    assert same(R3, R3.differentiate(R3.add(x3, one3)), [1], None)
    assert same(R3, R3.differentiate(r(3).pow_int(R3.add(x3, one3), 4)), [4, 12], 3)
    assert same(R10, R10.differentiate(R10.pow_int(x10, 4)), [0, 0, 0, 4], None)
    assert same(
        R10, R10.differentiate(R10.add(R10.multiply(x10, x10), x10)), [1, 2], None
    )


def test_rational_differentiate(rd_rational):
    R, x, _, R3, x3, one3, R10, x10, _, _ = rd_rational
    assert same(R, R.differentiate(R.pow_int(x, 3)), [(0, 1), (0, 1), (3, 1)], None)
    assert same(R, R.differentiate(R.add(R.multiply(x, x), x)), [(1, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.multiply(x3, x3)), [(0, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.add(x3, one3)), [(1, 1)], None)
    assert same(
        R10,
        R10.differentiate(R10.pow_int(x10, 4)),
        [(0, 1), (0, 1), (0, 1), (4, 1)],
        None,
    )
    assert same(
        R10,
        R10.differentiate(R10.add(R10.multiply(x10, x10), x10)),
        [(1, 1), (2, 1)],
        None,
    )


def test_rational_integrate(rd_rational):
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_rational
    assert same(R, R.integrate(R.add(x, one)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(
        R, R.integrate(R.multiply(x, x)), [(0, 1), (0, 1), (0, 1), (1, 3)], None
    )
    assert same(R3, R3.integrate(R3.add(x3, one3)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(
        R3, R3.integrate(R3.multiply(x3, x3)), [(0, 1), (0, 1), (0, 1), (1, 3)], None
    )
    assert same(R10, R10.integrate(R10.add(x10, one10)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(
        R10, R10.integrate(R10.pow_int(x10, 2)), [(0, 1), (0, 1), (0, 1), (1, 3)], None
    )


def test_error(rd_int):
    R, *_, r = rd_int
    raises(ValueError, lambda: r(-1))
    raises(ValueError, lambda: R.pow_int(R.gen, -1))
    raises(ValueError, lambda: R.truncate(R.gen, -1))
