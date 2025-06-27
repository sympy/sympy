from sympy.polys.domains import ZZ, QQ
from sympy.polys.series.python_powerseriesring import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)
from sympy.external.gmpy import GROUND_TYPES
import pytest
from sympy.testing.pytest import raises

flint = False
if GROUND_TYPES == "flint":
    from sympy.polys.series.flint_powerseriesring import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )

    flint = True

# Rings
Ring_ZZ = [PythonPowerSeriesRingZZ]
Ring_QQ = [PythonPowerSeriesRingQQ]

if flint:
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

@pytest.mark.parametrize("r", list(Ring_ZZ + Ring_QQ))
def test_equal(r):
    R = r()

    assert R.equal(R.from_list([1, 2, 3]), R.from_list([1, 2, 3])) is True
    assert R.equal(R.from_list([1, 2, 3], 3), R.from_list([1, 2, 3], 3)) is None
    assert R.equal(R.from_list([1, 2, 13], 3), R.from_list([1, 2, 3], 3)) is False
    assert R.equal(R.from_list([1, 2, 3], 3), R.from_list([1, 2, 3], 10)) is None
    assert R.equal(R.from_list([1, 21, 3], 3), R.from_list([1, 2, 3], 10)) is False

@pytest.mark.parametrize("r", list(Ring_ZZ + Ring_QQ))
def test_basics(r):
    R = r()
    R0 = r(0)
    gen0 = R0.gen

    assert R == r(6)
    assert hash(R) == hash(r(6))
    assert R.to_list(R.gen) == [0, 1], None
    assert R0.equal_repr(R0.zero, R0.one)
    assert R0.pretty(gen0) == "O(x**0)"

    assert same(R0, gen0, [], 0)
    assert same(R0, R0.multiply(gen0, gen0), [], 0)
    assert same(R, R.from_list([1, 2, 3, 4, 5, 6, 7]), [1, 2, 3, 4, 5, 6], 6)
    assert same(R, R.add(R.from_list([2, 4, 5], 3), R.from_list([5], 2)), [7, 4], 2)


@pytest.mark.parametrize("r", list(Ring_ZZ + Ring_QQ))
def test_int(r):
    R = r()
    x = R.gen
    one = R.one
    zero = R.zero

    R3 = r(3)
    x3 = R3.gen
    one3 = R3.one

    R10 = r(10)
    x10 = R10.gen
    one10 = R10.one

    assert same(R, zero, [], None)
    assert same(R, one, [1], None)
    assert same(R, x, [0, 1], None)

    assert same(R, R.negative(R.add(x, one)), [-1, -1], None)
    assert same(R3, R3.add(x3, one3), [1, 1], None)
    assert same(R3, R3.add(one3, R3.pow_int(x3, 4)), [1], 3)

    assert same(R, R.subtract(R.multiply(x, x), x), [0, -1, 1], None)
    assert same(R3, R3.subtract(R3.pow_int(R3.add(one3, x3), 4), x3), [1, 3, 6], 3)

    assert same(R, R.multiply(x, x), [0, 0, 1], None)
    assert same(R3, R3.multiply(R3.square(R3.add(x3, one3)),
        R3.add(x3, one3)), [1, 3, 3], 3)
    assert same(R10, R10.multiply(R10.add(x10, one10),
        R10.add(x10, one10)), [1, 2, 1], None)

    assert same(R, R.multiply_ground(one, 1), [1], None)
    assert same(R, R.multiply_ground(x, -1), [0, -1], None)
    assert same(R, R.multiply_ground(x, 0), [], None)
    assert same(R, R.multiply_ground(R.square(x), 1), [0, 0, 1], None)
    assert same(R, R.multiply_ground(R.add(x, one), ZZ(3)), [3, 3], None)

    assert not same(R, R.square(R.add(x, one)), [1, 1, 1], None)
    assert same(R3, R3.square(R3.multiply(x3, R3.add(x3, one3))),
        [0, 0, 1], 3)

    assert same(R, R.pow_int(x, 0), [1], None)
    assert same(R, R.pow_int(R.add(x, one), 6), [1, 6, 15, 20, 15, 6], 6)
    assert same(R, R.pow_int(x, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(x3, one3), 5), [1, 5, 10], 3)
    assert same(R10, R10.pow_int(x10, 7), [0, 0, 0, 0, 0, 0, 0, 1], None)
    assert same(R10, R10.pow_int(R10.add(x10, one10), 12),
        [1, 12, 66, 220, 495, 792, 924, 792, 495, 220], 10)

    assert same(R, R.truncate(R.pow_int(x, 3), 4), [0, 0, 0, 1], None)
    assert same(R, R.truncate(R.pow_int(R.add(x, one), 5), 3), [1, 5, 10], 3)

    assert same(R, R.differentiate(R.pow_int(x, 3)), [0, 0, 3], None)
    assert same(R, R.differentiate(R.add(R.multiply(x, x), x)), [1, 2], None)
    assert same(R3, R3.differentiate(R3.multiply(x3, x3)), [0, 2], None)
    assert same(R3, R3.differentiate(R3.add(x3, one3)), [1], None)
    assert same(R10, R10.differentiate(R10.pow_int(x10, 4)), [0, 0, 0, 4], None)
    assert same(R10, R10.differentiate(R10.add(R10.multiply(x10, x10), x10)),
        [1, 2], None)


@pytest.mark.parametrize(
    "r", Ring_QQ
)
def test_rational(r):
    R = r()
    x = R.gen
    one = R.one

    R3 = r(3)
    x3 = R3.gen
    one3 = R3.one

    R10 = r(10)
    x10 = R10.gen
    one10 = R10.one

    assert same(R, one, [(1, 1)], None)
    assert same(R, x, [(0, 1), (1, 1)], None)

    assert same(R, R.negative(R.add(x, one)), [(-1, 1), (-1, 1)], None)
    assert same(R, R.add(R.negative(one), one), [], None)
    assert same(R, R.subtract(R.multiply(x, x), x), [(0, 1), (-1, 1), (1, 1)], None)

    assert same(R, R.multiply(x, x), [(0, 1), (0, 1), (1, 1)], None)
    assert same(R, R.multiply(R.add(x, one), R.add(x, one)),
        [(1, 1), (2, 1), (1, 1)], None)
    assert same(R, R.multiply(R.subtract(x, one), R.add(x, one)),
        [(-1, 1), (0, 1), (1, 1)], None)
    assert same(R3, R3.multiply(R3.square(R3.add(x3, one3)),
        R3.add(x3, one3)), [(1, 1), (3, 1), (3, 1)], 3)

    assert same(R, R.multiply_ground(x, QQ(1, 2)), [(0, 1), (1, 2)], None)
    assert same(R, R.multiply(R.multiply_ground(x, QQ(1, 2)),
        R.multiply_ground(x, QQ(1, 3))), [(0, 1), (0, 1), (1, 6)], None)

    assert not same(R, R.square(R.add(x, one)), [(1, 1), (1, 1), (1, 1)],
        None)
    assert same(R10, R10.square(R10.subtract(x10, one10)),
        [(1, 1), (-2, 1), (1, 1)], None)

    assert same(R, R.pow_int(x, 3), [(0, 1), (0, 1), (0, 1), (1, 1)], None)
    assert same(R, R.pow_int(R.add(x, one), 6),
        [(1, 1), (6, 1), (15, 1), (20, 1), (15, 1), (6, 1)], 6)
    assert same(R, R.pow_int(x, 10), [], 6)
    assert same(R3, R3.pow_int(R3.add(x3, one3), 5), [(1, 1), (5, 1), (10, 1)], 3)
    assert same(R10, R10.pow_int(R10.add(x10, one10), 12),
        [(1, 1), (12, 1), (66, 1), (220, 1), (495, 1), (792, 1), (924, 1),
        (792, 1), (495, 1), (220, 1)], 10)

    assert same(R, R.differentiate(R.pow_int(x, 3)), [(0, 1), (0, 1), (3, 1)], None)
    assert same(R, R.differentiate(R.add(R.multiply(x, x), x)), [(1, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.multiply(x3, x3)), [(0, 1), (2, 1)], None)
    assert same(R3, R3.differentiate(R3.add(x3, one3)), [(1, 1)], None)
    assert same(R10, R10.differentiate(R10.pow_int(x10, 4)),
        [(0, 1), (0, 1), (0, 1), (4, 1)], None)
    assert same(R10, R10.differentiate(R10.add(R10.multiply(x10, x10), x10)),
        [(1, 1), (2, 1)], None)

    assert same(R, R.integrate(R.add(x, one)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(R, R.integrate(R.multiply(x, x)), [(0, 1), (0, 1), (0, 1), (1, 3)], None)
    assert same(R3, R3.integrate(R3.add(x3, one3)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(R3, R3.integrate(R3.multiply(x3, x3)),
        [(0, 1), (0, 1), (0, 1), (1, 3)], None)
    assert same(R10, R10.integrate(R10.add(x10, one10)), [(0, 1), (1, 1), (1, 2)], None)
    assert same(R10, R10.integrate(R10.pow_int(x10, 2)),
        [(0, 1), (0, 1), (0, 1), (1, 3)], None)

@pytest.mark.parametrize("r", list(Ring_ZZ + Ring_QQ))
def test_error(r):
    R = r()

    raises(ValueError, lambda: r(-1))
    raises(ValueError, lambda: R.pow_int(R.gen, -1))
    raises(ValueError, lambda: R.truncate(R.gen, -1))
