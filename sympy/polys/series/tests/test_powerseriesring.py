from sympy.polys.domains import ZZ, QQ
from sympy.polys.series.powerseriesring import PowerSeriesRing
from sympy.polys.series.python_powerseriesring import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)
from sympy.external.gmpy import GROUND_TYPES
import pytest
from sympy.testing.pytest import raises
from sympy.polys.polyerrors import NotReversible

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
    return R.equal_repr(s, R(coeffs_l, prec))


@pytest.fixture(params=Ring_ZZ + Ring_QQ)
def rd_int(request):
    return request.param


@pytest.fixture(params=Ring_QQ)
def rd_rational(request):
    return request.param


def test_equal(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    assert R.equal(R([1, 2, 3]), R([1, 2, 3])) is True
    assert R.equal(R([1, 21, 3]), R([1, 2, 3])) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3], 3)) is None
    assert R.equal(R([1, 2, 13], 3), R([1, 2, 3], 3)) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3], 10)) is None
    assert R.equal(R([1, 21, 3], 3), R([1, 2, 3], 10)) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3, 4, 5], 10)) is None
    assert R.equal(R([1, 2, 3, 4], None), R([1, 2, 3], 2)) is None
    assert R.equal(R([1, 2, 3, 4], None), R([1, 1, 3], 2)) is False
    assert R.equal(R([1, 2], None), R([1, 2, 3], 3)) is False


def test_basics(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R0 = SeriesRing(0)

    assert R == SeriesRing(6)
    assert hash(R) == hash(SeriesRing(6))
    assert R.to_list(R.gen) == [0, 1], None
    assert R0.equal_repr(R0.zero, R0.one)
    assert R0.pretty(R0.gen) == "O(x**0)"
    assert same(R0, R0.gen, [], 0)
    assert same(R0, R0.multiply(R0.gen, R0.gen), [], 0)
    assert same(R, R.add(R([2, 4, 5], 3), R([5], 2)), [7, 4], 2)


def test_positive(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    assert same(R, R.positive(R([1, 2, 3, 4, 5, 6, 7], None)), [1, 2, 3, 4, 5, 6], 6)


def test_negative(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    assert same(R, R.negative(R.gen), [0, -1], None)
    assert same(
        R,
        R.negative(R([1, 2, 3, 4, 5, 6, 7], None)),
        [-1, -2, -3, -4, -5, -6],
        6,
    )


def test_add(rd_int):
    SeriesRing = rd_int
    R3 = SeriesRing(3)
    assert same(R3, R3.add(R3.gen, R3.one), [1, 1], None)
    assert same(R3, R3.add(R3.one, R3.pow_int(R3.gen, 4)), [1], 3)


def test_int_subtract(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert same(R, R.subtract(R.multiply(R.gen, R.gen), R.gen), [0, -1, 1], None)
    assert same(
        R3, R3.subtract(R3.pow_int(R3.add(R3.one, R3.gen), 4), R3.gen), [1, 3, 6], 3
    )


def test_int_multiply(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.multiply(R.gen, R.gen), [0, 0, 1], None)
    assert same(
        R3,
        R3.multiply(R3.square(R3.add(R3.gen, R3.one)), R3.add(R3.gen, R3.one)),
        [1, 3, 3],
        3,
    )
    assert same(
        R10,
        R10.multiply(R10.add(R10.gen, R10.one), R10.add(R10.gen, R10.one)),
        [1, 2, 1],
        None,
    )


def test_rational_multiply(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert same(R, R.multiply(R.gen, R.gen), [(0, 1), (0, 1), (1, 1)], None)
    assert same(
        R,
        R.multiply(R.add(R.gen, R.one), R.add(R.gen, R.one)),
        [(1, 1), (2, 1), (1, 1)],
        None,
    )
    assert same(
        R,
        R.multiply(R.subtract(R.gen, R.one), R.add(R.gen, R.one)),
        [(-1, 1), (0, 1), (1, 1)],
        None,
    )
    assert same(
        R3,
        R3.multiply(R3.square(R3.add(R3.gen, R3.one)), R3.add(R3.gen, R3.one)),
        [(1, 1), (3, 1), (3, 1)],
        3,
    )


def test_int_multiply_ground(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    assert same(R, R.multiply_ground(R.one, 1), [1], None)
    assert same(R, R.multiply_ground(R.gen, -1), [0, -1], None)
    assert same(R, R.multiply_ground(R.gen, 0), [], None)
    assert same(R, R.multiply_ground(R.square(R.gen), 1), [0, 0, 1], None)
    assert same(R, R.multiply_ground(R.add(R.gen, R.one), ZZ(3)), [3, 3], None)


def test_rational_multiply_ground(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    assert same(R, R.multiply_ground(R.gen, QQ(1, 2)), [(0, 1), (1, 2)], None)
    assert same(R, R.multiply_ground(R.one, QQ(3, 4)), [(3, 4)], None)
    assert same(
        R, R.multiply_ground(R([(2, 3), (1, 5)]), QQ(1, 2)), [(1, 3), (1, 10)], None
    )
    assert same(
        R, R.multiply_ground(R.square(R.gen), QQ(2, 7)), [(0, 1), (0, 1), (2, 7)], None
    )


def test_int_divide(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.divide(R.add(R.gen, R.one), R.one), [1, 1], 6)
    assert same(
        R,
        R.divide(R.add(R.gen, R.one), R.subtract(R.one, R.gen)),
        [1, 2, 2, 2, 2, 2],
        6,
    )
    assert same(
        R,
        R.divide(R.pow_int(R.gen, 3), R.add(R.add(R.one, R.gen), R.square(R.gen))),
        [0, 0, 0, 1, -1, 0],
        6,
    )
    assert same(R3, R3.divide(R3([1, 2, 3]), R3.add(R3.one, R3.gen)), [1, 1, 2], 3)
    assert same(
        R10,
        R10.divide(R10([0, 0, 2, 4, 6, 8, 2, 14]), R10([0, 0, 1, 3, 4, 5, 1])),
        [2, -2, 4, -6, 12, -16, 26, -68],
        8,
    )

    raises(ZeroDivisionError, lambda: R.divide(R.gen, R.zero))
    raises(ValueError, lambda: R.divide(R.one, R.gen))
    raises(ValueError, lambda: R.divide(R([0, 1]), R([0, 0, 2])))


def test_rational_divide(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    assert same(
        R,
        R.divide(R.add(R.one, R.gen), R.subtract(R.one, R.gen)),
        [(1, 1), (2, 1), (2, 1), (2, 1), (2, 1), (2, 1)],
        6,
    )
    assert same(
        R,
        R.divide(R.pow_int(R.gen, 2), R.add(R.one, R.gen)),
        [(0, 1), (0, 1), (1, 1), (-1, 1), (1, 1), (-1, 1)],
        6,
    )
    assert same(
        R,
        R.divide(R([(1, 2), (3, 4), (1, 6)]), R([(1, 3), (2, 5)])),
        [(3, 2), (9, 20), (-1, 25), (6, 125), (-36, 625), (216, 3125)],
        6,
    )
    assert same(
        R,
        R.divide(R([(2, 3), (1, 4), (5, 6)]), R([(1, 2), (3, 7)])),
        [(4, 3), (-9, 14), (326, 147), (-652, 343), (3912, 2401), (-23472, 16807)],
        6,
    )


def test_int_square(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert not same(R, R.square(R.add(R.gen, R.one)), [1, 1, 1], None)
    assert same(
        R3, R3.square(R3.multiply(R3.gen, R3.add(R3.gen, R3.one))), [0, 0, 1], 3
    )


def test_rational_square(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R10 = SeriesRing(10)
    assert not same(R, R.square(R.add(R.gen, R.one)), [(1, 1), (1, 1), (1, 1)], None)
    assert same(
        R10, R10.square(R10.subtract(R10.gen, R10.one)), [(1, 1), (-2, 1), (1, 1)], None
    )


def test_int_pow(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.pow_int(R.gen, 0), [1], None)
    assert same(R, R.pow_int(R.add(R.gen, R.one), 6), [1, 6, 15, 20, 15, 6], 6)
    assert same(R, R.pow_int(R.gen, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(R3.gen, R3.one), 5), [1, 5, 10], 3)
    assert same(R10, R10.pow_int(R10.gen, 7), [0, 0, 0, 0, 0, 0, 0, 1], None)
    assert same(
        R10,
        R10.pow_int(R10.add(R10.gen, R10.one), 12),
        [1, 12, 66, 220, 495, 792, 924, 792, 495, 220],
        10,
    )


def test_rational_pow(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.pow_int(R.gen, 3), [(0, 1), (0, 1), (0, 1), (1, 1)], None)
    assert same(
        R,
        R.pow_int(R.add(R.gen, R.one), 6),
        [(1, 1), (6, 1), (15, 1), (20, 1), (15, 1), (6, 1)],
        6,
    )
    assert same(R, R.pow_int(R.gen, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(R3.gen, R3.one), 5), [(1, 1), (5, 1), (10, 1)], 3)
    assert same(
        R10,
        R10.pow_int(R10.add(R10.gen, R10.one), 12),
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
    SeriesRing = rd_int
    R = SeriesRing()
    assert same(R, R.truncate(R.pow_int(R.gen, 3), 4), [0, 0, 0, 1], None)
    assert same(R, R.truncate(R.pow_int(R.add(R.gen, R.one), 5), 3), [1, 5, 10], 3)


def test_int_differentiate(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.differentiate(R.pow_int(R.gen, 3)), [0, 0, 3], None)
    assert same(
        R, R.differentiate(R.add(R.multiply(R.gen, R.gen), R.gen)), [1, 2], None
    )
    assert same(R3, R3.differentiate(R3.multiply(R3.gen, R3.gen)), [0, 2], None)
    assert same(R3, R3.differentiate(R3.add(R3.gen, R3.one)), [1], None)
    assert same(
        R3,
        R3.differentiate(SeriesRing(3).pow_int(R3.add(R3.gen, R3.one), 4)),
        [4, 12],
        2,
    )
    assert same(R10, R10.differentiate(R10.pow_int(R10.gen, 4)), [0, 0, 0, 4], None)
    assert same(
        R10,
        R10.differentiate(R10.add(R10.multiply(R10.gen, R10.gen), R10.gen)),
        [1, 2],
        None,
    )


def test_rational_differentiate(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.differentiate(R.pow_int(R.gen, 3)), [(0, 1), (0, 1), (3, 1)], None)
    assert same(
        R,
        R.differentiate(R.add(R.multiply(R.gen, R.gen), R.gen)),
        [(1, 1), (2, 1)],
        None,
    )
    assert same(
        R3, R3.differentiate(R3.multiply(R3.gen, R3.gen)), [(0, 1), (2, 1)], None
    )
    assert same(R3, R3.differentiate(R3.add(R3.gen, R3.one)), [(1, 1)], None)
    assert same(
        R3, R3.differentiate(R.pow_int(R3.add(R3.gen, R3.one), 4)), [4, 12, 12], 3
    )
    assert same(
        R10,
        R10.differentiate(R10.pow_int(R10.gen, 4)),
        [(0, 1), (0, 1), (0, 1), (4, 1)],
        None,
    )
    assert same(
        R10,
        R10.differentiate(R10.add(R10.multiply(R10.gen, R10.gen), R10.gen)),
        [(1, 1), (2, 1)],
        None,
    )


def test_rational_integrate(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert same(R, R.integrate(R.add(R.gen, R.one)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(
        R, R.integrate(R.multiply(R.gen, R.gen)), [(0, 1), (0, 1), (0, 1), (1, 3)], None
    )
    assert same(
        R3, R3.integrate(R3.add(R3.gen, R3.one)), [(0, 1), (1, 1), (1, 2)], None
    )
    assert same(
        R3,
        R3.integrate(R3.multiply(R3.gen, R3.gen)),
        [(0, 1), (0, 1), (0, 1), (1, 3)],
        None,
    )
    assert same(
        R10, R10.integrate(R10.add(R10.gen, R10.one)), [(0, 1), (1, 1), (1, 2)], None
    )
    assert same(
        R10,
        R10.integrate(R10.pow_int(R10.gen, 2)),
        [(0, 1), (0, 1), (0, 1), (1, 3)],
        None,
    )
    assert same(
        R10,
        R10.integrate(R10([2, 3, 4], 3)),
        [(0, 1), (2, 1), (3, 2), (4, 3)],
        4,
    )


def test_error(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    raises(ValueError, lambda: SeriesRing(-1))
    raises(ValueError, lambda: R.pow_int(R.gen, -1))
    raises(ValueError, lambda: R.truncate(R.gen, -1))


def test_int_compose(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert same(R, R.compose(R([2]), R([0, 1, 1])), [2], None)
    assert same(R, R.compose(R.add(R.one, R.gen), R.gen), [1, 1], None)
    assert same(
        R, R.compose(R.multiply(R.gen, R.gen), R.add(R.one, R.gen)), [1, 2, 1], None
    )
    assert same(R, R.compose(R([1, 1, 1]), R.square(R.gen)), [1, 0, 1, 0, 1], None)
    assert same(R3, R3.compose(R3([1, 2, 3]), R3([0, 1, 1])), [1, 2, 5], 3)
    assert same(
        R10,
        R10.compose(R10([2, 4, 5, 1, 6, 2], 7), R10([0, 1, 1, 2, 3, 4], 7)),
        [2, 4, 9, 19, 46, 101, 206],
        7,
    )

    raises(ValueError, lambda: R.compose(R([1, 2], 2), R([1, 2, 3], 2)))


def test_rational_compose(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    f1 = R([(1, 2), (3, 4)])
    g1 = R([(1, 3), (2, 5)])
    assert same(R, R.compose(f1, g1), [(3, 4), (3, 10)], None)

    f3 = R([(2, 3), (5, 7)])
    g3 = R([(0, 1), (3, 4), (1, 6)])
    assert same(R, R.compose(f3, g3), [(2, 3), (15, 28), (5, 42)], None)

    f3_2 = R3([(1, 4), (1, 2), (1, 8)])
    g3_2 = R3([(0, 1), (1, 3)])
    assert same(R3, R3.compose(f3_2, g3_2), [(1, 4), (1, 6), (1, 72)], 3)

    f10 = R10([(1, 2), (3, 4), (5, 6)], 4)
    g10 = R10([(0, 1), (2, 3), (1, 5), (-3, 7)], 4)
    assert same(R10, R10.compose(f10, g10), [(1, 2), (1, 2), (281, 540), (-25, 252)], 4)


def test_int_inverse(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert same(R, R.inverse(R.add(R.one, R.gen)), [1, -1, 1, -1, 1, -1], 6)
    assert same(
        R,
        R.inverse(R.add(R.pow_int(R.multiply(R.add(R.one, R.gen), R.gen), 3), R.one)),
        [1, 0, 0, -1, -3, -3],
        6,
    )
    assert same(R3, R3.inverse(R3([1, 3, -2])), [1, -3, 11], 3)
    assert same(
        R10,
        R10.inverse(R10([1, -2, 3])),
        [1, 2, 1, -4, -11, -10, 13, 56, 73, -22],
        10,
    )

    raises(NotReversible, lambda: R.inverse(R.zero))
    raises(NotReversible, lambda: R.inverse(R([0, 1, 2])))


def test_rational_inverse(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    f_r = R([(2, 3), (-1, 4), (3, 5)])
    assert same(
        R,
        R.inverse(f_r),
        [
            (3, 2),
            (9, 16),
            (-729, 640),
            (-4779, 5120),
            (138267, 204800),
            (1791153, 1638400),
        ],
        6,
    )

    f3 = R3([(1, 2), (-3, 4), (5, 6)])
    assert same(R3, R3.inverse(f3), [(2, 1), (3, 1), (7, 6)], 3)

    f10 = R10([(3, 5), (1, 7), (-2, 9)])
    assert same(
        R10,
        R10.inverse(f10),
        [
            (5, 3),
            (-25, 63),
            (2825, 3969),
            (-26375, 83349),
            (1779875, 5250987),
            (-7274375, 36756909),
            (1199485625, 6947055801),
            (-16690759375, 145888171821),
            (838109346875, 9190954824723),
            (-12369018828125, 193010051319183),
        ],
        10,
    )


def test_int_reversion(rd_int):
    SeriesRing = rd_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert same(R, R.reversion(R.gen), [0, 1], 6)
    assert same(
        R, R.reversion(R.multiply(R.add(R.one, R.gen), R.gen)), [0, 1, -1, 2, -5, 14], 6
    )
    assert same(
        R,
        R.reversion(R([0, 1, 53, 2, 1, 3, 2])),
        [0, 1, -53, 5616, -743856, 110349083],
        6,
    )
    assert same(
        R3, R3.reversion(R3.multiply(R3.add(R3.one, R3.gen), R3.gen)), [0, 1, -1], 3
    )
    assert same(
        R10,
        R10.reversion(R10([0, 1, 2, -1, 3])),
        [0, 1, -2, 9, -53, 347, -2429, 17808, -134991, 1049422],
        10,
    )

    raises(NotReversible, lambda: R.reversion(R.zero))
    raises(NotReversible, lambda: R.reversion(R([0, 0, 2])))


def test_rational_reversion(rd_rational):
    SeriesRing = rd_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    f1 = R([(0, 1), (3, 2), (1, 4), (-2, 5)])
    assert same(
        R,
        R.reversion(f1),
        [(0, 1), (2, 3), (-2, 27), (116, 1215), (-106, 2187), (24604, 492075)],
        6,
    )

    f3 = R3([(0, 1), (5, 4), (-1, 3)])
    assert same(R3, R3.reversion(f3), [(0, 1), (4, 5), (64, 375)], 3)

    f10 = R10([(0, 1), (2, 3), (1, 5), (-3, 7), (1, 2)])
    assert same(
        R10,
        R10.reversion(f10),
        [
            (0, 1),
            (3, 2),
            (-27, 40),
            (486, 175),
            (-209709, 22400),
            (116634897, 3920000),
            (-2627149059, 22400000),
            (157012806591, 343000000),
            (-1627722349764129, 878080000000),
            (47642334773213367, 6146560000000),
        ],
        10,
    )
