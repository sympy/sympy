from sympy.external.gmpy import GROUND_TYPES
import pytest

flint = False
try:
    from flint import fmpq_series, fmpz_series, fmpz_poly, fmpq_poly
    if GROUND_TYPES == 'flint':
        flint = True
except ImportError:
    pass

pytestmark = pytest.mark.skipif(not flint, reason="python-flint is not available")

from sympy import QQ, ZZ
from sympy.polys.flint_powerseriesring import (FlintPowerSeriesRingZZ,
                                                FlintPowerSeriesRingQQ)
from sympy.testing.pytest import raises, XFAIL

def test_FlintPowerSeriesRingZZ():
    R = FlintPowerSeriesRingZZ(5)
    assert isinstance(R, FlintPowerSeriesRingZZ)
    assert R.domain == ZZ
    assert R.prec == 5
    assert R.equal(R.one, fmpz_poly([1]))
    assert R.equal(R.zero, fmpz_poly([0]))
    assert R.equal(R.gen, fmpz_poly([0, 1]))

def test_FlintPowerSeriesRingQQ():
    R = FlintPowerSeriesRingQQ()
    assert isinstance(R, FlintPowerSeriesRingQQ)
    assert R.domain == QQ
    assert R.prec == 6
    assert R.equal(R.one, fmpq_poly([1]))
    assert R.equal(R.zero, fmpq_poly([0]))
    assert R.equal(R.gen, fmpq_poly([0, 1]))

# Series for testing
R_ZZ = FlintPowerSeriesRingZZ(5)
R_QQ = FlintPowerSeriesRingQQ(10)

def test_negative():
    s1 = fmpz_series([1, -2, -3], prec=4)
    res = R_ZZ.negative(s1)
    assert res.coeffs() == [-1, 2, 3]
    s2 = fmpq_series([QQ(1, 2), QQ(3, 4)], prec=6)
    res = R_QQ.negative(s2)
    expected = [-c for c in s2.coeffs()]
    assert res.coeffs() == expected

def test_equal():
    s1 = fmpz_series([1, 2, 3], prec=4)
    s2 = fmpz_series([1, 2, 3], prec=4)
    assert R_ZZ.equal(s1, s2) is True
    s3 = fmpz_series([1, 2, 3], prec=4)
    p1 = fmpz_poly([1, 2])
    assert R_ZZ.equal(s3, p1) is False
    s4 = fmpq_series([QQ(1, 2), QQ(3, 4)], prec=6)
    p2 = fmpq_poly([QQ(1, 2), QQ(3, 4)])
    assert R_QQ.equal(s4, p2) is False
    p3 = fmpq_poly([QQ(1, 2), QQ(3, 4)])
    p4 = fmpq_poly([QQ(1, 2)])
    assert R_QQ.equal(p3, p4) is False

def test_add():

    QQDUP1 = [QQ(1, 2), QQ(3, 4), QQ(1, 2), QQ(7, 8), QQ(2, 1)]
    QQDUP2 = [QQ(1, 3), QQ(2, 5), QQ(3, 4)]

    # Both are series
    s1 = fmpz_series([1, 2, 3], prec=4)
    s2 = fmpz_series([4, 5], prec=3)
    exp1 = fmpz_series([5, 7, 3], prec=3)
    assert R_ZZ.equal(R_ZZ.add(s1, s2), exp1)
    s3 = fmpq_series(QQDUP1, prec=6)
    s4 = fmpq_series(QQDUP2, prec=5)
    exp2 = fmpq_series([QQ(5, 6), QQ(23, 20), QQ(5, 4), QQ(7, 8), QQ(2, 1)], prec=5)
    assert R_QQ.equal(R_QQ.add(s3, s4), exp2)

    # One is DUP and one series
    s5 = fmpz_series([7, 14, 21, 28], prec=4)
    p1 = fmpz_poly([0, 7])
    exp3 = fmpz_series([7, 21, 21, 28], prec=4)
    assert R_ZZ.equal(R_ZZ.add(s5, p1), exp3)
    s6 = fmpq_series(QQDUP1, prec=5)
    p2 = fmpq_poly(QQDUP2)
    exp4 = fmpq_series([QQ(5, 6), QQ(23, 20), QQ(5,4), QQ(7, 8), QQ(2, 1)], prec=5)
    assert R_QQ.equal(R_QQ.add(s6, p2), exp4)

    # Both DUPs
    p3 = fmpz_poly([1, 3, 1, 2, 2])
    p4 = fmpz_poly([10, 20])
    exp5 = fmpz_poly([11, 23, 1, 2, 2])
    assert R_ZZ.equal(R_ZZ.add(p3, p4), exp5)
    p5 = fmpq_poly(QQDUP1)
    p6 = fmpq_poly(QQDUP2)
    exp6 = fmpq_poly([QQ(5, 6), QQ(23, 20), QQ(5, 4), QQ(7, 8), QQ(2, 1)])
    assert R_QQ.equal(R_QQ.add(p5, p6), exp6)

def test_subtract():

    QQDUP1 = [QQ(3, 5), QQ(2, 3), QQ(4, 7), QQ(5, 9), QQ(1, 2)]
    QQDUP2 = [QQ(1, 4), QQ(1, 6), QQ(2, 5)]
    QQDUP3 = [QQ(5, 3), QQ(-2, 7), QQ(3, 11), QQ(7, 13)]

    # Both are series
    s1 = fmpz_series([5, 8, 6], prec=4)
    s2 = fmpz_series([2, 3], prec=3)
    exp1 = fmpz_series([3, 5, 6], prec=3)
    assert R_ZZ.equal(R_ZZ.subtract(s1, s2), exp1)
    s3 = fmpq_series(QQDUP1, prec=6)
    s4 = fmpq_series(QQDUP2, prec=5)
    exp2 = fmpq_series([QQ(7, 20), QQ(1, 2), QQ(6, 35), QQ(5, 9), QQ(1, 2)], prec=5)
    assert R_QQ.equal(R_QQ.subtract(s3, s4), exp2)

    # One is DUP and one series
    s5 = fmpz_series([9, 12, 15, 18], prec=4)
    p1 = fmpz_poly([3, 5])
    exp3 = fmpz_series([6, 7, 15, 18], prec=4)
    assert R_ZZ.equal(R_ZZ.subtract(s5, p1), exp3)
    s6 = fmpq_series(QQDUP1, prec=5)
    p2 = fmpq_poly(QQDUP3)
    exp4 = fmpq_series([QQ(-16, 15), QQ(20, 21), QQ(23, 77), QQ(2,117), QQ(1, 2)], prec=5)
    assert R_QQ.equal(R_QQ.subtract(s6, p2), exp4)

    # Both DUPs
    p3 = fmpz_poly([7, 4, 9, 2, 6])
    p4 = fmpz_poly([3, 2])
    exp5 = fmpz_poly([4, 2, 9, 2, 6])
    assert R_ZZ.equal(R_ZZ.subtract(p3, p4), exp5)
    p5 = fmpq_poly(QQDUP3)
    p6 = fmpq_poly(QQDUP2)
    exp6 = fmpq_poly([QQ(17, 12), QQ(-19, 42), QQ(-7, 55), QQ(7, 13)])
    assert R_QQ.equal(R_QQ.subtract(p5, p6), exp6)

def test_multiply():

    QQDUP1 = [QQ(1, 2), QQ(3, 4), QQ(1, 2), QQ(7, 8), QQ(2, 1)]
    QQDUP2 = [QQ(1, 3), QQ(2, 5), QQ(3, 4)]

    # Both are series
    s1 = fmpz_series([1, 2, 3], prec=4)
    s2 = fmpz_series([4, 5], prec=3)
    exp1 = fmpz_series([4, 13, 22], prec=3)
    assert R_ZZ.equal(R_ZZ.multiply(s1, s2), exp1)
    s3 = fmpq_series(QQDUP1, prec=6)
    s4 = fmpq_series(QQDUP2, prec=5)
    exp2 = fmpq_series([QQ(1, 6), QQ(9, 20), QQ(101, 120), QQ(253, 240), QQ(167, 120)], prec=5)
    assert R_QQ.equal(R_QQ.multiply(s3, s4), exp2)

    # One is DUP and one series
    s5 = fmpz_series([7, 14, 21, 28], prec=4)
    p1 = fmpz_poly([0, 7])
    exp3 = fmpz_series([0, 49, 98, 147], prec=4)
    assert R_ZZ.equal(R_ZZ.multiply(s5, p1), exp3)
    p2 = fmpq_poly(QQDUP1)
    s6 = fmpq_series(QQDUP2, prec=3)
    exp4 = fmpq_series([QQ(1,6), QQ(9,20), QQ(101,120)], prec=3)
    assert R_QQ.equal(R_QQ.multiply(p2, s6), exp4)

    # Both DUPs
    p3 = fmpz_poly([1, 3])
    p4 = fmpz_poly([10])
    exp5 = fmpz_poly([10,30])
    assert R_ZZ.equal(R_ZZ.multiply(p3, p4), exp5)
    p5 = fmpq_poly(QQDUP1)
    p6 = fmpq_poly(QQDUP2)
    exp6 = fmpq_poly([QQ(1, 6), QQ(9, 20), QQ(101, 120), QQ(253, 240), QQ(167, 120), QQ(233, 160), QQ(3, 2)])
    assert R_QQ.equal(R_QQ.multiply(p5, p6), exp6)

def test_multiply_ground():
    s1 = fmpz_series([1, 2, 3], prec=4)
    exp1 = fmpz_series([2, 4, 6], prec=4)
    assert R_ZZ.equal(R_ZZ.multiply_ground(s1, 2), exp1)
    s2 = fmpq_series([QQ(1, 2)], prec=6)
    exp2 = fmpq_series([QQ(3, 8)], prec=6)
    assert R_QQ.equal(R_QQ.multiply_ground(s2, QQ(3, 4)), exp2)

def test_trunc():
    s1 = fmpz_series([1, 2, 3, 4, 5], prec=6)
    exp1 = fmpz_series([1, 2, 3], prec=3)
    assert R_ZZ.equal(R_ZZ.trunc(s1, 3), exp1)
    s2 = fmpq_series([QQ(1, 2), QQ(3, 4), QQ(5, 6)], prec=6)
    exp2 = fmpq_series([QQ(1, 2), QQ(3, 4)], prec=2)
    assert R_QQ.equal(R_QQ.trunc(s2, 2), exp2)

    raises(ValueError, lambda: R_ZZ.trunc(fmpz_poly([1, 2]), -1))
    raises(ValueError, lambda: R_QQ.trunc(fmpq_poly([QQ(1, 2)]), -1))

@XFAIL
def test_special():
    R = FlintPowerSeriesRingZZ(5)

    p1 = fmpz_series([1, 2, 3], prec=4)
    p2 = fmpz_series([0, 0, 4, 5], prec=5)

    exp = fmpz_series([0, 0, 4, 13, 22], prec=5)
    assert R.equal(R.multiply(p1, p2), exp)
