from sympy.core.symbol import symbols
from sympy.series.order import O
from sympy.polys.domains import ZZ, QQ
from sympy.polys.series.base import PowerSeriesRingProto, PowerSeriesRingFieldProto
from sympy.polys.series.ring import (
    PowerSeriesElement,
    PowerSeriesRingRing,
    PowerSeriesRingField,
)
from sympy.polys.series.ringpython import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)
from sympy.external.gmpy import GROUND_TYPES
import pytest
from sympy.testing.pytest import raises
from sympy.polys.polyerrors import NotReversible
from sympy.polys.rings import ring
from sympy.polys.ring_series import (
    rs_log,
    rs_exp,
    rs_atan,
    rs_atanh,
    rs_asin,
    rs_asinh,
    rs_tan,
    rs_tanh,
    rs_sin,
    rs_sinh,
    rs_cos,
    rs_cosh,
)
from sympy.polys.densebasic import dup_random

# Rings
Ring_ZZ: list[type[PowerSeriesRingProto]] = [PythonPowerSeriesRingZZ]
Ring_QQ: list[type[PowerSeriesRingProto]] = [PythonPowerSeriesRingQQ]

ground_int: list[type[PowerSeriesRingProto]]
ground_rational: list[type[PowerSeriesRingFieldProto]]

flint: bool = False

if GROUND_TYPES == "flint":
    flint = True
    from sympy.polys.series.ringflint import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )

    Ring_ZZ.append(FlintPowerSeriesRingZZ)
    Ring_QQ.append(FlintPowerSeriesRingQQ)

    ground_int = [FlintPowerSeriesRingZZ, FlintPowerSeriesRingQQ]
    ground_rational = [FlintPowerSeriesRingQQ]

else:
    ground_int = [PythonPowerSeriesRingZZ, PythonPowerSeriesRingQQ]
    ground_rational = [PythonPowerSeriesRingQQ]


@pytest.fixture(params=Ring_ZZ + Ring_QQ)
def ring_int(request):
    return request.param


@pytest.fixture(params=Ring_QQ)
def ring_rational(request):
    return request.param


@pytest.fixture(params=ground_int)
def groundring_int(request):
    return request.param


@pytest.fixture(params=ground_rational)
def groundring_rational(request):
    return request.param


def test_equal(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    assert R.equal(R([1, 2, 3]), R([1, 2, 3])) is True
    assert R.equal(R([1, 21, 3]), R([1, 2, 3])) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3], 3)) is None
    assert R.equal(R([1, 2, 13], 3), R([1, 2, 3], 3)) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3], 10)) is None
    assert R.equal(R([1, 21, 3], 3), R([1, 2, 3], 10)) is False
    assert R.equal(R([1, 2, 3], 3), R([1, 2, 3, 4, 5], 10)) is None
    assert R.equal(R([1, 2, 3, 4], None), R([1, 2, 3], 2)) is None
    assert R.equal(R([1, 2, 3], 2), R([1, 2, 3, 4], None)) is None
    assert R.equal(R([1, 2, 3, 4], None), R([1, 1, 3], 2)) is False
    assert R.equal(R([1, 2], None), R([1, 2, 3], 3)) is False


def test_basics(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R0 = SeriesRing(0)

    assert R == SeriesRing(6)
    assert hash(R) == hash(SeriesRing(6))
    assert R.to_list(R.gen) == [0, 1], None
    assert R.to_dense(R.gen) == [1, 0], None
    assert R0.equal_repr(R0.zero, R0.one)
    assert R0.pretty(R0.gen) == "O(x**0)"
    assert R.pretty(R([1, 2, 3])) == "1 + 2*x + 3*x**2"
    assert R0.equal_repr(R0.gen, R0([], 0))
    assert not R.equal_repr(R([1, 2, 3]), R([1, 2, 3], 5))
    assert R0.equal_repr(R0.multiply(R0.gen, R0.gen), R0([], 0))
    assert R.equal_repr(R.add(R([2, 4, 5], 3), R([5], 2)), R([7, 4], 2))

    assert R.series_prec(R([1, 2, 3])) == None
    assert R.series_prec(R([1, 2, 3], None)) == None
    assert R.series_prec(R([1, 2, 3, 4, 5, 6, 7])) == 6
    assert R.series_prec(R([1, 2, 3, 4, 5, 6, 7], 10)) == 6


def test_positive(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert R.equal_repr(
        R.positive(R([1, 2, 3, 4, 5, 6, 7], None)), R([1, 2, 3, 4, 5, 6], 6)
    )
    assert R.equal_repr(
        R.positive(R([1, 2, 3, 4, 5, 6, 7, 8, 9], 8)), R([1, 2, 3, 4, 5, 6], 6)
    )
    assert R3.equal_repr(R3.positive(R([1, 2, 7, 9, 12])), R3([1, 2, 7], 3))
    assert R3.equal_repr(R3.positive(R([1, 2, 7, 9, 12], 6)), R3([1, 2, 7], 3))


def test_negative(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    assert R.equal_repr(R.negative(R.gen), R([0, -1], None))
    assert R.equal_repr(
        R.negative(R([1, 2, 3, 4, 5, 6, 7], None)),
        R([-1, -2, -3, -4, -5, -6], 6),
    )
    assert R.equal_repr(
        R.negative(R([1, 2, 3, 4, 5, 6, 7, 8, 9], 8)),
        R([-1, -2, -3, -4, -5, -6], 6),
    )


def test_add(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert R3.equal_repr(R3.add(R3.gen, R3.one), R3([1, 1], None))
    assert R3.equal_repr(R3.add(R3.one, R3.pow_int(R3.gen, 4)), R3([1], 3))
    assert R3.equal_repr(R3.add(R([1, 2, 4, 2, 1, 1]), R3.one), R3([2, 2, 4], 3))


def test_int_subtract(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert R3.equal_repr(R3.subtract(R([1, 2, 4, 2, 1, 1]), R3.one), R3([0, 2, 4], 3))
    assert R.equal_repr(
        R.subtract(R.multiply(R.gen, R.gen), R.gen), R([0, -1, 1], None)
    )
    assert R3.equal_repr(
        R3.subtract(R3.pow_int(R3.add(R3.one, R3.gen), 4), R3.gen), R3([1, 3, 6], 3)
    )


def test_int_multiply(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(R.multiply(R.gen, R.gen), R([0, 0, 1], None))
    assert R3.equal_repr(
        R3.multiply(R([1, 1, 2, 1, 1]), R([1, 1, 2, 1, 1])), R3([1, 2, 5], 3)
    )
    assert R3.equal_repr(
        R3.multiply(R3.square(R3.add(R3.gen, R3.one)), R3.add(R3.gen, R3.one)),
        R3([1, 3, 3], 3),
    )
    assert R10.equal_repr(
        R10.multiply(R10.add(R10.gen, R10.one), R10.add(R10.gen, R10.one)),
        R10([1, 2, 1], None),
    )


def test_rational_multiply(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert R.equal_repr(R.multiply(R.gen, R.gen), R([(0, 1), (0, 1), (1, 1)], None))
    assert R.equal_repr(
        R.multiply(R.add(R.gen, R.one), R.add(R.gen, R.one)),
        R([(1, 1), (2, 1), (1, 1)], None),
    )
    assert R.equal_repr(
        R.multiply(R.subtract(R.gen, R.one), R.add(R.gen, R.one)),
        R([(-1, 1), (0, 1), (1, 1)], None),
    )
    assert R3.equal_repr(
        R3.multiply(R3.square(R3.add(R3.gen, R3.one)), R3.add(R3.gen, R3.one)),
        R3([(1, 1), (3, 1), (3, 1)], 3),
    )


def test_int_multiply_ground(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert R.equal_repr(R.multiply_ground(R.one, ZZ(1)), R([1], None))
    assert R.equal_repr(R.multiply_ground(R.gen, ZZ(-1)), R([0, -1], None))
    assert R.equal_repr(R.multiply_ground(R.gen, ZZ(0)), R([], None))
    assert R.equal_repr(R.multiply_ground(R.square(R.gen), 1), R([0, 0, 1], None))
    assert R.equal_repr(R.multiply_ground(R.add(R.gen, R.one), ZZ(3)), R([3, 3], None))
    assert R.equal_repr(
        R.multiply_ground(R.inverse(R.add(R.one, R.gen)), ZZ(7)),
        R([7, -7, 7, -7, 7, -7], 6),
    )
    assert R3.equal_repr(
        R3.multiply_ground(R([1, 2, 1, 1, 1]), ZZ(2)), R3([2, 4, 2], 3)
    )


def test_rational_multiply_ground(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    assert R.equal_repr(R.multiply_ground(R.gen, QQ(1, 2)), R([(0, 1), (1, 2)], None))
    assert R.equal_repr(R.multiply_ground(R.one, QQ(3, 4)), R([(3, 4)], None))
    assert R.equal_repr(
        R.multiply_ground(R([(2, 3), (1, 5)]), QQ(1, 2)), R([(1, 3), (1, 10)], None)
    )
    assert R.equal_repr(
        R.multiply_ground(R.square(R.gen), QQ(2, 7)), R([(0, 1), (0, 1), (2, 7)], None)
    )


def test_int_divide(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(R.divide(R.add(R.gen, R.one), R.one), R([1, 1], None))
    assert R.equal_repr(
        R.divide(R.add(R.gen, R.one), R.subtract(R.one, R.gen)),
        R([1, 2, 2, 2, 2, 2], 6),
    )
    assert R.equal_repr(
        R.divide(R.pow_int(R.add(R.gen, R.one), 3), R.add(R.one, R.gen)),
        R([1, 2, 1], None),
    )
    assert R.equal_repr(
        R.divide(R.pow_int(R.gen, 3), R.add(R.add(R.one, R.gen), R.square(R.gen))),
        R([0, 0, 0, 1, -1, 0], 6),
    )
    assert R3.equal_repr(
        R3.divide(R3([1, 2, 3]), R3.add(R3.one, R3.gen)), R3([1, 1, 2], 3)
    )
    assert R10.equal_repr(
        R10.divide(R10([0, 0, 2, 4, 6, 8, 2, 14]), R10([0, 0, 1, 3, 4, 5, 1])),
        R10([2, -2, 4, -6, 12, -16, 26, -68, 168, -346], 10),
    )

    raises(ZeroDivisionError, lambda: R.divide(R.gen, R.zero))
    raises(ValueError, lambda: R.divide(R.one, R.gen))
    raises(ValueError, lambda: R.divide(R([0, 1]), R([0, 0, 2])))


def test_rational_divide(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    assert R.equal_repr(
        R.divide(R.add(R.one, R.gen), R.subtract(R.one, R.gen)),
        R([(1, 1), (2, 1), (2, 1), (2, 1), (2, 1), (2, 1)], 6),
    )
    assert R.equal_repr(
        R.divide(R.pow_int(R.gen, 2), R.add(R.one, R.gen)),
        R([(0, 1), (0, 1), (1, 1), (-1, 1), (1, 1), (-1, 1)], 6),
    )
    assert R.equal_repr(
        R.divide(R([(1, 2), (3, 4), (1, 6)]), R([(1, 3), (2, 5)])),
        R([(3, 2), (9, 20), (-1, 25), (6, 125), (-36, 625), (216, 3125)], 6),
    )
    assert R.equal_repr(
        R.divide(R([(2, 3), (1, 4), (5, 6)]), R([(1, 2), (3, 7)])),
        R(
            [(4, 3), (-9, 14), (326, 147), (-652, 343), (3912, 2401), (-23472, 16807)],
            6,
        ),
    )


def test_int_square(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    assert not R.equal_repr(R.square(R.add(R.gen, R.one)), R([1, 1, 1], None))
    assert R3.equal_repr(
        R3.square(R3.multiply(R3.gen, R3.add(R3.gen, R3.one))), R3([0, 0, 1], 3)
    )


def test_rational_square(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R10 = SeriesRing(10)
    assert not R.equal_repr(
        R.square(R.add(R.gen, R.one)), R([(1, 1), (1, 1), (1, 1)], None)
    )
    assert R10.equal_repr(
        R10.square(R10.subtract(R10.gen, R10.one)), R10([(1, 1), (-2, 1), (1, 1)], None)
    )


def test_sqrt(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R10 = SeriesRing(10)

    assert R.equal_repr(R.sqrt(R([])), R([]))
    raises(ValueError, lambda: R.sqrt(R([3, 2])))
    raises(ValueError, lambda: R.sqrt(R([0, 1, 7])))

    assert R.equal_repr(R.sqrt(R([1, 2, 1])), R([1, 1], None))
    assert R.equal_repr(R.sqrt(R([4, 4, 1])), R([2, 1], None))
    assert R.equal_repr(R.sqrt(R([9])), R([3], None))

    assert R10.equal_repr(
        R10.sqrt(R10.subtract(R10.one, R10.gen)),
        R10(
            [
                (1, 1),
                (-1, 2),
                (-1, 8),
                (-1, 16),
                (-5, 128),
                (-7, 256),
                (-21, 1024),
                (-33, 2048),
                (-429, 32768),
                (-715, 65536),
            ],
            10,
        ),
    )

    assert R10.equal_repr(
        R10.sqrt(R10([1, 0, 0, 1])),
        R10(
            [
                (1, 1),
                (0, 1),
                (0, 1),
                (1, 2),
                (0, 1),
                (0, 1),
                (-1, 8),
                (0, 1),
                (0, 1),
                (1, 16),
            ],
            10,
        ),
    )


def test_int_pow(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(R.pow_int(R.gen, 0), R([1], None))
    assert R.equal_repr(R.pow_int(R([1, 2, 1], 5), 3), R([1, 6, 15, 20, 15], 5))
    assert R.equal_repr(R.pow_int(R.add(R.gen, R.one), 6), R([1, 6, 15, 20, 15, 6], 6))
    assert R.equal_repr(R.pow_int(R.gen, 10), R([], 6))
    assert R.equal_repr(R.pow_int(R.add(R.one, R.gen), 3), R([1, 3, 3, 1], None))
    assert R3.equal_repr(R3.pow_int(R3.add(R3.gen, R3.one), 5), R3([1, 5, 10], 3))
    assert R3.equal_repr(R3.pow_int(R([1, 1, 1, 1, 1]), 2), R3([1, 2, 3], 3))
    assert R3.equal_repr(
        R3.pow_int(R3.pow_int(R3.add(R3.one, R3.gen), 2), 2), R3([1, 4, 6], 3)
    )
    assert R10.equal_repr(R10.pow_int(R10.gen, 7), R10([0, 0, 0, 0, 0, 0, 0, 1], None))
    assert R10.equal_repr(
        R10.pow_int(R10.add(R10.gen, R10.one), 12),
        R10([1, 12, 66, 220, 495, 792, 924, 792, 495, 220], 10),
    )

    raises(ValueError, lambda: R.pow_int(R.gen, -2))


def test_rational_pow(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(R.pow_int(R.gen, 3), R([(0, 1), (0, 1), (0, 1), (1, 1)], None))
    assert R.equal_repr(
        R.pow_int(R.add(R.gen, R.one), 6),
        R([(1, 1), (6, 1), (15, 1), (20, 1), (15, 1), (6, 1)], 6),
    )
    assert R.equal_repr(R.pow_int(R.gen, 10), R([], 6))
    assert R3.equal_repr(
        R3.pow_int(R3.add(R3.gen, R3.one), 5), R3([(1, 1), (5, 1), (10, 1)], 3)
    )
    assert R10.equal_repr(
        R10.pow_int(R10.add(R10.gen, R10.one), 12),
        R10(
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
        ),
    )


def test_truncate(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    assert R.equal_repr(R.truncate(R.pow_int(R.gen, 3), 4), R([0, 0, 0, 1], None))
    assert R.equal_repr(
        R.truncate(R.pow_int(R.add(R.gen, R.one), 5), 3), R([1, 5, 10], 3)
    )


def test_int_differentiate(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(R.differentiate(R.pow_int(R.gen, 3)), R([0, 0, 3], None))
    assert R.equal_repr(
        R.differentiate(R.add(R.multiply(R.gen, R.gen), R.gen)), R([1, 2], None)
    )
    assert R3.equal_repr(
        R3.differentiate(R3.multiply(R3.gen, R3.gen)), R3([0, 2], None)
    )
    assert R3.equal_repr(R3.differentiate(R3.add(R3.gen, R3.one)), R3([1], None))
    assert R3.equal_repr(
        R3.differentiate(SeriesRing(3).pow_int(R3.add(R3.gen, R3.one), 4)),
        R3([4, 12], 2),
    )
    assert R3.equal_repr(R3.differentiate(R([2, 6, 2, 1, 2, 3])), R3([6, 4, 3], 3))
    assert R10.equal_repr(
        R10.differentiate(R10.pow_int(R10.gen, 4)), R10([0, 0, 0, 4], None)
    )
    assert R10.equal_repr(
        R10.differentiate(R10.add(R10.multiply(R10.gen, R10.gen), R10.gen)),
        R10([1, 2], None),
    )


def test_rational_differentiate(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(
        R.differentiate(R.pow_int(R.gen, 3)), R([(0, 1), (0, 1), (3, 1)], None)
    )
    assert R.equal_repr(
        R.differentiate(R.add(R.multiply(R.gen, R.gen), R.gen)),
        R([(1, 1), (2, 1)], None),
    )
    assert R3.equal_repr(
        R3.differentiate(R3.multiply(R3.gen, R3.gen)), R3([(0, 1), (2, 1)], None)
    )
    assert R3.equal_repr(R3.differentiate(R3.add(R3.gen, R3.one)), R3([(1, 1)], None))
    assert R3.equal_repr(
        R3.differentiate(R.pow_int(R3.add(R3.gen, R3.one), 4)), R3([4, 12, 12], 3)
    )
    assert R10.equal_repr(
        R10.differentiate(R10.pow_int(R10.gen, 4)),
        R10([(0, 1), (0, 1), (0, 1), (4, 1)], None),
    )
    assert R10.equal_repr(
        R10.differentiate(R10.add(R10.multiply(R10.gen, R10.gen), R10.gen)),
        R10([(1, 1), (2, 1)], None),
    )


def test_rational_integrate(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)
    assert R.equal_repr(
        R.integrate(R.add(R.gen, R.one)), R([(0, 1), (1, 1), (1, 2)], None)
    )
    assert R.equal_repr(
        R.integrate(R.multiply(R.gen, R.gen)), R([(0, 1), (0, 1), (0, 1), (1, 3)], None)
    )
    assert R3.equal_repr(
        R3.integrate(R3.add(R3.gen, R3.one)), R3([(0, 1), (1, 1), (1, 2)], None)
    )
    assert R3.equal_repr(
        R3.integrate(R3.multiply(R3.gen, R3.gen)),
        R3([(0, 1), (0, 1), (0, 1), (1, 3)], None),
    )
    assert R3.equal_repr(R3.integrate(R([2, 4, 1, 1], 5)), R3([0, 2, 2], 3))
    assert R10.equal_repr(
        R10.integrate(R10.add(R10.gen, R10.one)), R10([(0, 1), (1, 1), (1, 2)], None)
    )
    assert R10.equal_repr(
        R10.integrate(R10.pow_int(R10.gen, 2)),
        R10([(0, 1), (0, 1), (0, 1), (1, 3)], None),
    )
    assert R10.equal_repr(
        R10.integrate(R10([2, 3, 4], 3)),
        R10([(0, 1), (2, 1), (3, 2), (4, 3)], 4),
    )


def test_error(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    raises(ValueError, lambda: SeriesRing(-1))
    raises(ValueError, lambda: R.pow_int(R.gen, -1))
    raises(ValueError, lambda: R.truncate(R.gen, -1))


def test_int_compose(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.compose(R([2]), R([0, 1, 1])), R([2], None))
    assert R.equal_repr(R.compose(R.add(R.one, R.gen), R.gen), R([1, 1], None))
    assert R.equal_repr(
        R.compose(R([1, 1, 1]), R.square(R.gen)), R([1, 0, 1, 0, 1], None)
    )
    assert R3.equal_repr(R3.compose(R3([1, 2, 3]), R3([0, 1, 1])), R3([1, 2, 5], 3))
    assert R3.equal_repr(R3.compose(R3([7, 5, 8], 3), R3([0, 1, 1])), R3([7, 5, 13], 3))
    assert R3.equal_repr(
        R3.compose(R3([1, 2, 3]), R3([0, 9, 1], 3)), R3([1, 18, 245], 3)
    )
    assert R3.equal_repr(
        R3.compose(R([7, 1, 1, 2, 3]), R([0, 2, 2, 1, 1])), R3([7, 2, 6], 3)
    )

    assert R10.equal_repr(
        R10.compose(R10([2, 4, 5, 1, 6, 2], 7), R10([0, 1, 1, 2, 3, 4], 7)),
        R10([2, 4, 9, 19, 46, 101, 206], 7),
    )

    raises(ValueError, lambda: R.compose(R([1, 2], 2), R([1, 2, 3], 2)))


def test_rational_compose(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    f1 = R([(1, 2), (3, 4)])
    g1 = R([(0, 1), (2, 5)])
    assert R.equal_repr(R.compose(f1, g1), R([(1, 2), (3, 10)], None))

    f3 = R([(2, 3), (5, 7)])
    g3 = R([(0, 1), (3, 4), (1, 6)])
    assert R.equal_repr(R.compose(f3, g3), R([(2, 3), (15, 28), (5, 42)], None))

    f3_2 = R3([(1, 4), (1, 2), (1, 8)])
    g3_2 = R3([(0, 1), (1, 3)])
    assert R3.equal_repr(R3.compose(f3_2, g3_2), R3([(1, 4), (1, 6), (1, 72)], None))

    f10 = R10([(1, 2), (3, 4), (5, 6)], 4)
    g10 = R10([(0, 1), (2, 3), (1, 5), (-3, 7)], 4)
    assert R10.equal_repr(
        R10.compose(f10, g10), R10([(1, 2), (1, 2), (281, 540), (-25, 252)], 4)
    )


def test_int_inverse(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.inverse(R.one), R([1], None))
    assert R.equal_repr(R.inverse(R.add(R.one, R.gen)), R([1, -1, 1, -1, 1, -1], 6))
    assert R.equal_repr(
        R.inverse(R.add(R.pow_int(R.multiply(R.add(R.one, R.gen), R.gen), 3), R.one)),
        R([1, 0, 0, -1, -3, -3], 6),
    )
    assert R3.equal_repr(R3.inverse(R3([1, 3, -2])), R3([1, -3, 11], 3))
    assert R10.equal_repr(
        R10.inverse(R10([1, -2, 3])),
        R10([1, 2, 1, -4, -11, -10, 13, 56, 73, -22], 10),
    )

    raises(NotReversible, lambda: R.inverse(R.zero))
    raises(NotReversible, lambda: R.inverse(R([0, 1, 2])))


def test_rational_inverse(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.inverse(R([(-2, 3)])), R([(-3, 2)], None))

    f_r = R([(2, 3), (-1, 4), (3, 5)])
    assert R.equal_repr(
        R.inverse(f_r),
        R(
            [
                (3, 2),
                (9, 16),
                (-729, 640),
                (-4779, 5120),
                (138267, 204800),
                (1791153, 1638400),
            ],
            6,
        ),
    )

    f3 = R3([(1, 2), (-3, 4), (5, 6)])
    assert R3.equal_repr(R3.inverse(f3), R3([(2, 1), (3, 1), (7, 6)], 3))

    f10 = R10([(3, 5), (1, 7), (-2, 9)])
    assert R10.equal_repr(
        R10.inverse(f10),
        R10(
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
        ),
    )


def test_int_reversion(ring_int):
    SeriesRing = ring_int
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.reversion(R.gen), R([0, 1], 6))
    assert R.equal_repr(
        R.reversion(R.multiply(R.add(R.one, R.gen), R.gen)),
        R([0, 1, -1, 2, -5, 14], 6),
    )
    assert R.equal_repr(
        R.reversion(R([0, 1, 53, 2, 1, 3, 2])),
        R([0, 1, -53, 5616, -743856, 110349083], 6),
    )
    assert R3.equal_repr(
        R3.reversion(R3.multiply(R3.add(R3.one, R3.gen), R3.gen)),
        R3([0, 1, -1], 3),
    )
    assert R10.equal_repr(
        R10.reversion(R10([0, 1, 2, -1, 3])),
        R10([0, 1, -2, 9, -53, 347, -2429, 17808, -134991, 1049422], 10),
    )

    raises(NotReversible, lambda: R.reversion(R.zero))
    raises(NotReversible, lambda: R.reversion(R([0, 0, 2])))


def test_rational_reversion(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    f1 = R([(0, 1), (3, 2), (1, 4), (-2, 5)])
    assert R.equal_repr(
        R.reversion(f1),
        R([(0, 1), (2, 3), (-2, 27), (116, 1215), (-106, 2187), (24604, 492075)], 6),
    )

    f3 = R3([(0, 1), (5, 4), (-1, 3)])
    assert R3.equal_repr(R3.reversion(f3), R3([(0, 1), (4, 5), (64, 375)], 3))

    f10 = R10([(0, 1), (2, 3), (1, 5), (-3, 7), (1, 2)])
    assert R10.equal_repr(
        R10.reversion(f10),
        R10(
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
        ),
    )


def test_log(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.log(R([1])), R([]))
    raises(ValueError, lambda: R.log(R([])))
    raises(ValueError, lambda: R.log(R([(0, 1)])))
    raises(ValueError, lambda: R.log(R([(2, 1)])))

    raises(ValueError, lambda: R.log1p(R([1, 1])))

    assert R3.equal_repr(
        R3.log(R3.add(R3.one, R3.gen)),
        R3([(0, 1), (1, 1), (-1, 2)], 3),
    )
    assert R10.equal_repr(
        R10.log(R10.add(R10.one, R10.gen)),
        R10(
            [
                (0, 1),
                (1, 1),
                (-1, 2),
                (1, 3),
                (-1, 4),
                (1, 5),
                (-1, 6),
                (1, 7),
                (-1, 8),
                (1, 9),
            ],
            10,
        ),
    )


def test_exp(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    raises(ValueError, lambda: R.exp(R([1, 1])))
    assert R.equal_repr(R.exp(R([])), R([1]))
    assert R.equal_repr(R.expm1(R([])), R([]))
    assert R.equal_repr(
        R.expm1(R([0, 2, 8, 1])),
        R([(0, 1), (2, 1), (10, 1), (55, 3), (152, 3), (1274, 15)], 6),
    )

    assert R3.equal_repr(
        R3.exp(R3.gen),
        R3([(1, 1), (1, 1), (1, 2)], 3),
    )
    assert R10.equal_repr(
        R10.exp(R10.gen),
        R10(
            [
                (1, 1),
                (1, 1),
                (1, 2),
                (1, 6),
                (1, 24),
                (1, 120),
                (1, 720),
                (1, 5040),
                (1, 40320),
                (1, 362880),
            ],
            10,
        ),
    )


def test_atan(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.atan(R([])), R([]))
    raises(ValueError, lambda: R.atan(R([2])))

    assert R3.equal_repr(
        R3.atan(R3.gen),
        R3([(0, 1), (1, 1)], 3),
    )

    assert R10.equal_repr(
        R10.atan(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (-1, 3),
                (0, 1),
                (1, 5),
                (0, 1),
                (-1, 7),
                (0, 1),
                (1, 9),
            ],
            10,
        ),
    )


def test_atanh(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.atanh(R([])), R([]))
    raises(ValueError, lambda: R.atanh(R([1, 1])))

    assert R3.equal_repr(
        R3.atanh(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.atanh(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (1, 3),
                (0, 1),
                (1, 5),
                (0, 1),
                (1, 7),
                (0, 1),
                (1, 9),
            ],
            10,
        ),
    )


def test_asin(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.asin(R([])), R([]))
    raises(ValueError, lambda: R.asin(R([1, 1])))

    assert R3.equal_repr(
        R3.asin(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.asin(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (1, 6),
                (0, 1),
                (3, 40),
                (0, 1),
                (5, 112),
                (0, 1),
                (35, 1152),
            ],
            10,
        ),
    )


def test_asinh(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.asinh(R([])), R([]))
    raises(ValueError, lambda: R.asinh(R([1, 1])))

    assert R3.equal_repr(
        R3.asinh(R3.gen),
        R3([(0, 1), (1, 1)], 3),
    )
    assert R10.equal_repr(
        R10.asinh(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (-1, 6),
                (0, 1),
                (3, 40),
                (0, 1),
                (-5, 112),
                (0, 1),
                (35, 1152),
            ],
            10,
        ),
    )


def test_tan(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.tan(R([])), R([]))
    raises(ValueError, lambda: R.tan(R([1, 1])))

    assert R3.equal_repr(
        R3.tan(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.tan(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (1, 3),
                (0, 1),
                (2, 15),
                (0, 1),
                (17, 315),
                (0, 1),
                (62, 2835),
            ],
            10,
        ),
    )


def test_tanh(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.tanh(R([])), R([]))
    raises(ValueError, lambda: R.tanh(R([1, 1])))

    assert R3.equal_repr(
        R3.tanh(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.tanh(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (-1, 3),
                (0, 1),
                (2, 15),
                (0, 1),
                (-17, 315),
                (0, 1),
                (62, 2835),
            ],
            10,
        ),
    )


def test_sin(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.sin(R([])), R([]))
    raises(ValueError, lambda: R.sin(R([1, 1])))

    assert R3.equal_repr(
        R3.sin(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.sin(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (-1, 6),
                (0, 1),
                (1, 120),
                (0, 1),
                (-1, 5040),
                (0, 1),
                (1, 362880),
            ],
            10,
        ),
    )


def test_sinh(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.sinh(R([])), R([]))
    raises(ValueError, lambda: R.sinh(R([1, 1])))

    assert R3.equal_repr(
        R3.sinh(R3.gen),
        R3([(0, 1), (1, 1), (0, 1)], 3),
    )
    assert R10.equal_repr(
        R10.sinh(R10.gen),
        R10(
            [
                (0, 1),
                (1, 1),
                (0, 1),
                (1, 6),
                (0, 1),
                (1, 120),
                (0, 1),
                (1, 5040),
                (0, 1),
                (1, 362880),
            ],
            10,
        ),
    )


def test_cos(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.cos(R([])), R([1]))
    raises(ValueError, lambda: R.cos(R([1, 1])))

    assert R3.equal_repr(
        R3.cos(R3.gen),
        R3([(1, 1), (0, 1), (-1, 2)], 3),
    )
    assert R10.equal_repr(
        R10.cos(R10.gen),
        R10(
            [
                (1, 1),
                (0, 1),
                (-1, 2),
                (0, 1),
                (1, 24),
                (0, 1),
                (-1, 720),
                (0, 1),
                (1, 40320),
                (0, 1),
            ],
            10,
        ),
    )


def test_cosh(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R3 = SeriesRing(3)
    R10 = SeriesRing(10)

    assert R.equal_repr(R.cosh(R([])), R([1]))
    raises(ValueError, lambda: R.cosh(R([1, 1])))

    assert R3.equal_repr(
        R3.cosh(R3.gen),
        R3([(1, 1), (0, 1), (1, 2)], 3),
    )
    assert R10.equal_repr(
        R10.cosh(R10.gen),
        R10(
            [
                (1, 1),
                (0, 1),
                (1, 2),
                (0, 1),
                (1, 24),
                (0, 1),
                (1, 720),
                (0, 1),
                (1, 40320),
                (0, 1),
            ],
            10,
        ),
    )


def test_high_deg(ring_rational):
    SeriesRing = ring_rational
    Rs = SeriesRing(30)
    Rs50 = SeriesRing(50)
    Rp, x = ring("x", QQ)

    rand = dup_random(30, 0, 10, QQ)
    rand50 = dup_random(50, 0, 10, QQ)
    rand[0] = QQ.zero
    rand50[0] = QQ.zero
    s = Rs(rand)
    s50 = Rs50(rand50)
    p = Rp.from_list(rand[::-1])
    p50 = Rp.from_list(rand50[::-1])

    assert Rs.to_dense(Rs.exp(s)) == rs_exp(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.expm1(s)) == (rs_exp(p, x, 30) - 1).to_dense()

    assert Rs.to_dense(Rs.atan(s)) == rs_atan(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.atanh(s)) == rs_atanh(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.asin(s)) == rs_asin(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.asinh(s)) == rs_asinh(p, x, 30).to_dense()

    assert Rs.to_dense(Rs.tan(s)) == rs_tan(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.tanh(s)) == rs_tanh(p, x, 30).to_dense()
    assert Rs.to_dense(Rs.sin(s)) == rs_sin(p, x, 30).to_dense()
    assert Rs.to_dense(Rs50.sinh(s50)) == rs_sinh(p50, x, 50).to_dense()
    assert Rs.to_dense(Rs.cos(s)) == rs_cos(p, x, 30).to_dense()
    assert Rs.to_dense(Rs50.cosh(s50)) == rs_cosh(p50, x, 50).to_dense()

    rand[0] = QQ.one
    s2 = Rs(rand)
    p2 = Rp.from_list(rand[::-1])
    assert Rs.to_dense(Rs.log(s2)) == rs_log(p2, x, 30).to_dense()
    assert Rs.to_dense(Rs.log1p(s)) == rs_log(p2, x, 30).to_dense()


def test_hypot(ring_rational):
    SeriesRing = ring_rational
    R = SeriesRing()
    R10 = SeriesRing(10)

    assert R.equal_repr(R.hypot(R([(3, 1)]), R([(4, 1)])), R([(5, 1)], None))

    assert R.equal_repr(
        R.hypot(R([(1, 1), (0, 1), (-1, 1)]), R([(0, 1), (2, 1)])),
        R([(1, 1), (0, 1), (1, 1)], None),
    )

    assert R10.equal_repr(
        R10.hypot(R10([(1, 1)]), R10.gen),
        R10(
            [
                (1, 1),
                (0, 1),
                (1, 2),
                (0, 1),
                (-1, 8),
                (0, 1),
                (1, 16),
                (0, 1),
                (-5, 128),
                (0, 1),
            ],
            10,
        ),
    )

    assert R.equal_repr(
        R.hypot(R([(1, 1), (1, 1)]), R([(0, 1)])), R([(1, 1), (1, 1)], None)
    )


def test_PowerSeriesRing():
    from sympy.abc import x

    RZZ = PowerSeriesRingRing(ZZ, "x")
    RQQ = PowerSeriesRingField(QQ, "x")

    assert RZZ == PowerSeriesRingRing(ZZ, "x", 6)
    assert RQQ == PowerSeriesRingField(QQ, "x", 6)
    assert RZZ != PowerSeriesRingRing(ZZ, "x", 10)
    assert RQQ != PowerSeriesRingField(QQ, "x", 10)
    assert RZZ != RQQ

    assert RZZ.domain == ZZ
    assert RQQ.domain == QQ
    assert RZZ.gen == RZZ(x)
    assert RQQ.gen == RQQ(x)
    assert isinstance(RZZ.gen, PowerSeriesElement)
    assert isinstance(RQQ.gen, PowerSeriesElement)
    assert RZZ.prec == 6
    assert RQQ.prec == 6
    assert RZZ.symbol == RQQ.symbol


def test_PowerSeriesRing_basics():
    RZZ = PowerSeriesRingRing(ZZ, "x")
    RQQ = PowerSeriesRingField(QQ, "x")
    xz = RZZ.gen
    xq = RQQ.gen
    _x = symbols("x")

    assert RZZ.is_element(xz)
    assert not RZZ.is_element(xq)

    assert RZZ.truncate(xz + xz**2 + xz**3, 3) == RZZ.from_expr(_x + _x**2 + O(_x**3))
    assert (
        RZZ.from_list([ZZ(1), ZZ(2), ZZ(3), ZZ(4), ZZ(2), ZZ(13), ZZ(2)])
        == 1
        + 2 * xz
        + 3 * xz**2
        + 4 * xz**3
        + 2 * xz**4
        + 13 * xz**5
        + RZZ.order_term()
    )
    assert (
        RQQ.from_list([QQ(1, 3), QQ(3, 2), QQ(1, 3)])
        == QQ(1, 3) + 3 * xq / 2 + xq**2 / 3
    )

    assert RZZ.to_list(xz + 2 * xz**2 + 3 * xz**3) == [0, 1, 2, 3]
    assert RZZ.to_dense(xz + 2 * xz**2 + 3 * xz**3) == [3, 2, 1, 0]


def test_PowerSeriesElement_basics():
    R = PowerSeriesRingField(QQ, "x")
    Ry = PowerSeriesRingField(QQ, "y")
    x = R.gen
    y = Ry.gen

    _x = symbols("x")

    assert (QQ(2) + x) == x + QQ(2)
    assert (QQ(2) - x) == 2 - x
    assert (x - QQ(2)) == x - 2
    assert (QQ(2) * x) == x * QQ(2)
    assert (QQ(1) / (1 + x)) == 1 / (1 + x)
    assert ((1 + x) / QQ(3)) == (1 + x) / QQ(3)

    assert (2 + x + 3 * x**2).as_expr() == 2 + _x + 3 * _x**2
    assert (0 + 3 * x**2).as_expr() == 3 * _x**2
    assert (2 + x + 3 * x**2 + R.order_term()).as_expr() == 2 + _x + 3 * _x**2 + O(
        _x**6
    )
    assert (2 + x + 3 * x**2).constant_coefficient() == 2
    assert (2 + x + R.order_term()).removeO() == 2 + x
    assert R.from_int(5).is_ground
    assert not (1 + x + x**2).is_ground

    raises(ValueError, lambda: (1 + x + y))
    raises(ValueError, lambda: (1 + x - y))
    raises(ValueError, lambda: (1 + x * y))
    raises(ValueError, lambda: (1 + x / y))

    raises(TypeError, lambda: (x + _x))
    raises(TypeError, lambda: (x - _x))
    raises(TypeError, lambda: (x * _x))
    raises(TypeError, lambda: (x / _x))
    raises(TypeError, lambda: (_x + x))
    raises(TypeError, lambda: (_x - x))
    raises(TypeError, lambda: (_x * x))
    raises(TypeError, lambda: (_x / x))


def test_PowerSeriesRing_from_expr():
    R = PowerSeriesRingField(QQ, "x", 5)
    o = R.order_term()

    x = R.gen
    _x = symbols("x")

    assert R.from_expr(_x) == x
    assert R.from_expr(1 + _x + 2 * _x**2) == 1 + x + 2 * x**2
    assert (
        R.from_expr(2 + 4 * _x + 6 * _x**2 + 6 * _x**3 + _x**5)
        == 2 + 4 * x + 6 * x**2 + 6 * x**3 + o
    )
    assert R.from_expr(_x + O(_x**7)) == x + o
    assert (
        R.from_expr(1 + _x + (3 + _x**2) ** (-2) + _x**5)
        == QQ(10, 9) + x - 2 * x**2 / 27 + x**4 / 27 + x**5 + o
    )

    raises(ValueError, lambda: R.from_expr(symbols("y") ** 4))
    raises(ValueError, lambda: R.from_expr(O(symbols("y") ** 4)))
    raises(ValueError, lambda: R.from_expr(_x ** (-2)))
    raises(ValueError, lambda: R.from_expr(_x ** (1 / 2)))


def test_PowerSeriesRing_ring_new():
    R = PowerSeriesRingRing(QQ, "x", 5)

    assert R.ring_new(QQ(7)) == R.from_ground(QQ(7))
    assert R.ring_new(3) == R.from_int(3)


def test_PowerSeriesRing_arith(groundring_int):
    SeriesRing = groundring_int
    RL = SeriesRing(5)
    if RL.domain.is_ZZ:
        RU = PowerSeriesRingRing(RL.domain, "x", 5)
    else:
        RU = PowerSeriesRingField(RL.domain, "x", 5)
    x = RU.gen

    p1 = 1 + 2 * x + 3 * x**2
    p2 = 4 + 5 * x + 6 * x**2
    assert RL.equal_repr((p1 + p2).series, RL([5, 7, 9]))
    assert RL.equal_repr((p2 - p1).series, RL([3, 3, 3]))

    assert RL.equal_repr((x * x).series, RL([0, 0, 1]))
    assert RL.equal_repr(((1 + x) * (1 + x)).series, RL([1, 2, 1]))
    assert RL.equal_repr(((1 + 2 * x) * (3 + 4 * x)).series, RL([3, 10, 8]))

    assert RL.equal_repr((x**0).series, RL([1]))
    assert RL.equal_repr((x**1).series, RL([0, 1]))
    assert RL.equal_repr((x**3).series, RL([0, 0, 0, 1]))
    assert RL.equal_repr(((1 + x) ** 2).series, RL([1, 2, 1]))
    assert RL.equal_repr(((1 + x) ** 3).series, RL([1, 3, 3, 1]))

    p3 = 1 - x + 2 * x**2 + x**6
    assert RL.equal_repr(p3.series, RL([1, -1, 2, 0, 0, 0, 1]))
    assert RL.equal_repr((-p3).series, RL([-1, 1, -2, 0, 0, 0, -1]))

    p4 = 2 + 4 * x + 6 * x**2
    assert RL.equal_repr((p4 / 2).series, RL([1, 2, 3]))
    assert RL.equal_repr((p4 / (1 + x**2 + x**5)).series, RL([2, 4, 4, -4, -4], 5))

    assert RU.square(1 + x) == (1 + x) ** 2
    assert RL.equal_repr(RU.square(1 + x).series, RL([1, 2, 1]))
    assert RL.equal_repr(RU.square(x**2 + 7 * x**3).series, RL([0, 0, 0, 0, 1], 5))


def test_PowerSeriesRing_operations_int(groundring_int):
    SeriesRing = groundring_int
    RL = SeriesRing(5)
    if RL.domain.is_ZZ:
        RU = PowerSeriesRingRing(RL.domain, "x", 5)
    else:
        RU = PowerSeriesRingField(RL.domain, "x", 5)
    x = RU.gen

    assert RL.equal_repr(
        RU.compose(1 + x**2 + 12 * x**3, x + 7 * x**3).series,
        RL([1, 0, 1, 12, 14], 5),
    )
    assert RL.equal_repr(
        RU.compose(3 * x**3, x + 2 * x**3 + 2 * x**4).series,
        RL([0, 0, 0, 3], 5),
    )

    assert RL.equal_repr(
        RU.differentiate(x**3 + 2 * x**2 + 3 * x + 4).series, RL([3, 4, 3])
    )
    assert RL.equal_repr(
        RU.differentiate(x**4 + 5 * x**3 + 2 * x).series, RL([2, 0, 15, 4])
    )

    assert RL.equal_repr(RU.inverse(1 + x).series, RL([1, -1, 1, -1, 1], 5))
    assert RL.equal_repr(RU.inverse(1 + 2 * x + x**2).series, RL([1, -2, 3, -4, 5], 5))

    assert RL.equal_repr(RU.reversion(x + x**2).series, RL([0, 1, -1, 2, -5], 5))
    assert RL.equal_repr(
        RU.reversion(x + 2 * x**2 + x**3).series, RL([0, 1, -2, 7, -30], 5)
    )


def test_PowerSeriesRing_operations_rational(groundring_rational):
    SeriesRing = groundring_rational
    R = SeriesRing(5)
    RU = PowerSeriesRingField(QQ, "x", 5)
    x = RU.gen

    assert R.equal_repr(
        RU.compose(
            QQ(1, 2) + QQ(2, 3) * x**2 + QQ(3, 4) * x**3, x + QQ(1, 5) * x**3
        ).series,
        R([QQ(1, 2), QQ(0, 1), QQ(2, 3), QQ(3, 4), QQ(4, 15)], 5),
    )
    assert R.equal_repr(
        RU.compose(
            QQ(3, 4) + QQ(1, 6) * x**2 + QQ(2, 7) * x**4,
            x + QQ(3, 8) * x**2 + QQ(1, 9) * x**3,
        ).series,
        R([QQ(3, 4), QQ(0, 1), QQ(1, 6), QQ(1, 8), QQ(8375, 24192)], 5),
    )

    assert R.equal_repr(
        RU.differentiate(
            QQ(5, 6)
            + QQ(2, 3) * x
            + QQ(4, 7) * x**2
            + QQ(1, 8) * x**3
            + QQ(3, 5) * x**4
        ).series,
        R([QQ(2, 3), QQ(8, 7), QQ(3, 8), QQ(12, 5)]),
    )
    assert R.equal_repr(
        RU.differentiate(QQ(1, 9) * x**2 + QQ(5, 11) * x**3 + QQ(2, 13) * x**4).series,
        R([QQ(0, 1), QQ(2, 9), QQ(15, 11), QQ(8, 13)]),
    )

    p5 = QQ(2, 3) + QQ(1, 4) * x + QQ(3, 5) * x**2
    assert R.equal_repr(
        (p5 * RU.inverse(p5)).series,
        R([QQ(1, 1)], 5),
    )
    p6 = QQ(3, 7) + QQ(2, 9) * x + QQ(1, 11) * x**2 + QQ(4, 13) * x**3
    assert R.equal_repr(
        (p6 * RU.inverse(p6)).series,
        R([QQ(1, 1)], 5),
    )

    p7 = x + QQ(1, 3) * x**2 + QQ(2, 5) * x**3
    assert R.equal_repr(
        RU.compose(p7, RU.reversion(p7)).series,
        R([QQ(0, 1), QQ(1, 1)], 5),
    )
    p8 = x + QQ(3, 7) * x**3 + QQ(1, 11) * x**4
    assert R.equal_repr(
        RU.compose(p8, RU.reversion(p8)).series,
        R([QQ(0, 1), QQ(1, 1)], 5),
    )

    assert R.equal_repr(
        RU.integrate(
            QQ(1, 2) + QQ(3, 4) * x + QQ(2, 3) * x**2 + QQ(5, 6) * x**3
        ).series,
        R([QQ(0, 1), QQ(1, 2), QQ(3, 8), QQ(2, 9), QQ(5, 24)]),
    )
    assert R.equal_repr(
        RU.integrate(QQ(2, 3) * x**2 + QQ(4, 5) * x**3 + QQ(1, 7) * x**4).series,
        R([QQ(0, 1), QQ(0, 1), QQ(0, 1), QQ(2, 9), QQ(1, 5), QQ(1, 35)]),
    )

    assert R.equal_repr(RU.sqrt(x**2 + 2 * x + 1).series, R([1, 1]))
    assert R.equal_repr(
        RU.sqrt(1 + 5 * x + 2 * x**2 + 12 * x**3 + 12 * x**5).series,
        R([QQ(1, 1), QQ(5, 2), QQ(-17, 8), QQ(181, 16), QQ(-3909, 128)], 5),
    )


def test_PowerSeriesRing_series():
    R = PowerSeriesRingField(QQ, "x", 10)
    x = R.gen
    o = R.order_term()

    p1 = x + 2 * x + 3 * x**4

    assert (
        R.log(1 + p1)
        == 3 * x
        - 9 * x**2 / 2
        + 9 * x**3
        - 69 * x**4 / 4
        + 198 * x**5 / 5
        - 189 * x**6 / 2
        + 1620 * x**7 / 7
        - 4653 * x**8 / 8
        + 1485 * x**9
        + o
    )
    assert R.log1p(p1) == R.log(1 + p1)
    assert (
        R.log1p(p1)
        == 3 * x
        - 9 * x**2 / 2
        + 9 * x**3
        - 69 * x**4 / 4
        + 198 * x**5 / 5
        - 189 * x**6 / 2
        + 1620 * x**7 / 7
        - 4653 * x**8 / 8
        + 1485 * x**9
        + o
    )
    assert (
        R.exp(p1)
        == 1
        + 3 * x
        + 9 * x**2 / 2
        + 9 * x**3 / 2
        + 51 * x**4 / 8
        + 441 * x**5 / 40
        + 1161 * x**6 / 80
        + 7803 * x**7 / 560
        + 66249 * x**8 / 4480
        + 87939 * x**9 / 4480
        + o
    )
    assert R.expm1(p1) == R.exp(p1) - 1
    assert (
        R.expm1(p1)
        == 3 * x
        + 9 * x**2 / 2
        + 9 * x**3 / 2
        + 51 * x**4 / 8
        + 441 * x**5 / 40
        + 1161 * x**6 / 80
        + 7803 * x**7 / 560
        + 66249 * x**8 / 4480
        + 87939 * x**9 / 4480
        + o
    )
    assert (
        R.atan(p1)
        == 3 * x
        - 9 * x**3
        + 3 * x**4
        + 243 * x**5 / 5
        - 27 * x**6
        - 2187 * x**7 / 7
        + 243 * x**8
        + 2160 * x**9
        + o
    )
    assert (
        R.atanh(p1)
        == 3 * x
        + 9 * x**3
        + 3 * x**4
        + 243 * x**5 / 5
        + 27 * x**6
        + 2187 * x**7 / 7
        + 243 * x**8
        + 2214 * x**9
        + o
    )
    assert (
        R.asin(p1)
        == 3 * x
        + 9 * x**3 / 2
        + 3 * x**4
        + 729 * x**5 / 40
        + 27 * x**6 / 2
        + 10935 * x**7 / 112
        + 729 * x**8 / 8
        + 78273 * x**9 / 128
        + o
    )
    assert (
        R.asinh(p1)
        == 3 * x
        - 9 * x**3 / 2
        + 3 * x**4
        + 729 * x**5 / 40
        - 27 * x**6 / 2
        - 10935 * x**7 / 112
        + 729 * x**8 / 8
        + 74817 * x**9 / 128
        + o
    )
    assert (
        R.tan(p1)
        == 3 * x
        + 9 * x**3
        + 3 * x**4
        + 162 * x**5 / 5
        + 27 * x**6
        + 4131 * x**7 / 35
        + 162 * x**8
        + 16011 * x**9 / 35
        + o
    )
    assert (
        R.tanh(p1)
        == 3 * x
        - 9 * x**3
        + 3 * x**4
        + 162 * x**5 / 5
        - 27 * x**6
        - 4131 * x**7 / 35
        + 162 * x**8
        + 14121 * x**9 / 35
        + o
    )
    assert (
        R.sin(p1)
        == 3 * x
        - 9 * x**3 / 2
        + 3 * x**4
        + 81 * x**5 / 40
        - 27 * x**6 / 2
        - 243 * x**7 / 560
        + 81 * x**8 / 8
        - 60237 * x**9 / 4480
        + o
    )
    assert (
        R.sinh(p1)
        == 3 * x
        + 9 * x**3 / 2
        + 3 * x**4
        + 81 * x**5 / 40
        + 27 * x**6 / 2
        + 243 * x**7 / 560
        + 81 * x**8 / 8
        + 60723 * x**9 / 4480
        + o
    )
    assert (
        R.cos(p1)
        == 1
        - 9 * x**2 / 2
        + 27 * x**4 / 8
        - 9 * x**5
        - 81 * x**6 / 80
        + 27 * x**7 / 2
        - 19431 * x**8 / 4480
        - 243 * x**9 / 40
        + o
    )
    assert (
        R.cosh(p1)
        == 1
        + 9 * x**2 / 2
        + 27 * x**4 / 8
        + 9 * x**5
        + 81 * x**6 / 80
        + 27 * x**7 / 2
        + 20889 * x**8 / 4480
        + 243 * x**9 / 40
        + o
    )
