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
    domain = R.domain
    domain_coeffs = []
    for coeff in coeffs_l:
        if isinstance(coeff, tuple):
            domain_coeffs.append(domain(*coeff))
        else:
            domain_coeffs.append(domain(coeff))

    return R.equal_repr(s, R.from_list(domain_coeffs, prec))


@pytest.fixture(params=Ring_ZZ + Ring_QQ)
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
    R, x, one, _, _, _, R10, x10, one10, _ = rd_rational
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
    assert same(R3, R3.differentiate(r(3).pow_int(R3.add(x3, one3), 4)), [4, 12], 2)
    assert same(R10, R10.differentiate(R10.pow_int(x10, 4)), [0, 0, 0, 4], None)
    assert same(
        R10, R10.differentiate(R10.add(R10.multiply(x10, x10), x10)), [1, 2], None
    )


def test_rational_differentiate(rd_rational):
    R, x, _, R3, x3, one3, R10, x10, *_ = rd_rational
    assert same(R, R.differentiate(R.pow_int(x, 3)), [(0, 1), (0, 1), (3, 1)], None)
    assert same(R, R.differentiate(R.add(R.multiply(x, x), x)), [(1, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.multiply(x3, x3)), [(0, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.add(x3, one3)), [(1, 1)], None)
    assert same(R3, R3.differentiate(R.pow_int(R3.add(x3, one3), 4)), [4, 12, 12], 3)
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
    assert same(
        R10,
        R10.integrate(R10.from_list([2, 3, 4], 3)),
        [(0, 1), (2, 1), (3, 2), (4, 3)],
        4,
    )


def test_error(rd_int):
    R, *_, r = rd_int
    raises(ValueError, lambda: r(-1))
    raises(ValueError, lambda: R.pow_int(R.gen, -1))
    raises(ValueError, lambda: R.truncate(R.gen, -1))


def test_int_compose(rd_int):
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_int

    assert same(R, R.compose(R.from_list([2]), R.from_list([0, 1, 1])), [2], None)
    assert same(R, R.compose(R.add(one, x), x), [1, 1], None)
    assert same(R, R.compose(R.multiply(x, x), R.add(one, x)), [1, 2, 1], None)
    assert same(
        R, R.compose(R.from_list([1, 1, 1]), R.square(x)), [1, 0, 1, 0, 1], None
    )
    assert same(
        R3, R3.compose(R3.from_list([1, 2, 3]), R3.from_list([0, 1, 1])), [1, 2, 5], 3
    )
    assert same(
        R10,
        R10.compose(
            R10.from_list([2, 4, 5, 1, 6, 2], 7), R10.from_list([0, 1, 1, 2, 3, 4], 7)
        ),
        [2, 4, 9, 19, 46, 101, 206],
        7,
    )

    raises(
        ValueError, lambda: R.compose(R.from_list([1, 2], 2), R.from_list([1, 2, 3], 2))
    )


def test_rational_compose(rd_rational):
    R, _, _, R3, _, _, R10, *_ = rd_rational

    f1 = R.from_list([QQ(1, 2), QQ(3, 4)])
    g1 = R.from_list([QQ(1, 3), QQ(2, 5)])
    assert same(R, R.compose(f1, g1), [(3, 4), (3, 10)], None)

    f3 = R.from_list([QQ(2, 3), QQ(5, 7)])
    g3 = R.from_list([QQ(0, 1), QQ(3, 4), QQ(1, 6)])
    assert same(R, R.compose(f3, g3), [(2, 3), (15, 28), (5, 42)], None)

    f3_2 = R3.from_list([QQ(1, 4), QQ(1, 2), QQ(1, 8)])
    g3_2 = R3.from_list([QQ(0, 1), QQ(1, 3)])
    assert same(R3, R3.compose(f3_2, g3_2), [(1, 4), (1, 6), (1, 72)], 3)

    f10 = R10.from_list([QQ(1, 2), QQ(3, 4), QQ(5, 6)], 4)
    g10 = R10.from_list([QQ(0, 1), QQ(2, 3), QQ(1, 5), QQ(-3, 7)], 4)
    assert same(R10, R10.compose(f10, g10), [(1, 2), (1, 2), (281, 540), (-25, 252)], 4)


def test_int_inversion(rd_int):
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_int

    assert same(R, R.inversion(R.add(one, x)), [1, -1, 1, -1, 1, -1], 6)
    assert same(
        R,
        R.inversion(R.add(R.pow_int(R.multiply(R.add(one, x), x), 3), one)),
        [1, 0, 0, -1, -3, -3],
        6,
    )
    assert same(R3, R3.inversion(R3.from_list([1, 3, -2])), [1, -3, 11], 3)
    assert same(
        R10,
        R10.inversion(R10.from_list([1, -2, 3])),
        [1, 2, 1, -4, -11, -10, 13, 56, 73, -22],
        10,
    )

    raises(NotReversible, lambda: R.inversion(R.zero))
    raises(NotReversible, lambda: R.inversion(R.from_list([0, 1, 2])))


def test_rational_inversion(rd_rational):
    R, _, _, R3, _, _, R10, *_ = rd_rational

    f_r = R.from_list([QQ(2, 3), QQ(-1, 4), QQ(3, 5)])
    assert same(
        R,
        R.inversion(f_r),
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

    f3 = R3.from_list([QQ(1, 2), QQ(-3, 4), QQ(5, 6)])
    assert same(R3, R3.inversion(f3), [(2, 1), (3, 1), (7, 6)], 3)

    f10 = R10.from_list([QQ(3, 5), QQ(1, 7), QQ(-2, 9)])
    assert same(
        R10,
        R10.inversion(f10),
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
    R, x, one, R3, x3, one3, R10, x10, one10, _ = rd_int

    assert same(R, R.reversion(x), [0, 1], 6)
    assert same(R, R.reversion(R.multiply(R.add(one, x), x)), [0, 1, -1, 2, -5, 14], 6)
    assert same(
        R,
        R.reversion(R.from_list([0, 1, 53, 2, 1, 3, 2])),
        [0, 1, -53, 5616, -743856, 110349083],
        6,
    )
    assert same(R3, R3.reversion(R3.multiply(R3.add(one3, x3), x3)), [0, 1, -1], 3)
    assert same(
        R10,
        R10.reversion(R10.from_list([0, 1, 2, -1, 3])),
        [0, 1, -2, 9, -53, 347, -2429, 17808, -134991, 1049422],
        10,
    )

    raises(NotReversible, lambda: R.reversion(R.zero))
    raises(NotReversible, lambda: R.reversion(R.from_list([0, 0, 2])))


def test_rational_reversion(rd_rational):
    R, _, _, R3, _, _, R10, *_ = rd_rational

    f1 = R.from_list([QQ(0, 1), QQ(3, 2), QQ(1, 4), QQ(-2, 5)])
    assert same(
        R,
        R.reversion(f1),
        [(0, 1), (2, 3), (-2, 27), (116, 1215), (-106, 2187), (24604, 492075)],
        6,
    )

    f3 = R3.from_list([QQ(0, 1), QQ(5, 4), QQ(-1, 3)])
    assert same(R3, R3.reversion(f3), [(0, 1), (4, 5), (64, 375)], 3)

    f10 = R10.from_list([QQ(0, 1), QQ(2, 3), QQ(1, 5), QQ(-3, 7), QQ(1, 2)])
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
