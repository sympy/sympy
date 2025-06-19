from sympy import QQ, ZZ
from sympy.polys.python_powerseriesring import (PythonPowerSeriesRingZZ,
                                                PythonPowerSeriesRingQQ)
from sympy.testing.pytest import raises

def test_PowerSeriesRingZZ():
    R = PythonPowerSeriesRingZZ(5)
    assert isinstance(R, PythonPowerSeriesRingZZ)
    assert R.domain == ZZ
    assert R.prec == 5
    assert R.one == ([1], None)
    assert R.zero == ([0], None)
    assert R.gen == ([0, 1], None)

def test_PowerSeriesRingQQ():
    R = PythonPowerSeriesRingQQ()
    assert isinstance(R, PythonPowerSeriesRingQQ)
    assert R.domain == QQ
    assert R.prec == 6
    assert R.one == ([QQ(1)], None)
    assert R.zero == ([QQ(0)], None)
    assert R.gen == ([QQ(0), QQ(1)], None)

# Series for testing
R_ZZ = PythonPowerSeriesRingZZ(5)
R_QQ = PythonPowerSeriesRingQQ(10)

def test_negative():
    assert R_ZZ.negative(([1, -2, -3], 4)) == ([-1, 2, 3], 4)
    assert R_QQ.negative(([QQ(1, 2), QQ(3, 4)], 6)) == ([QQ(-1, 2), QQ(-3, 4)], 6)

def test_equal():
    assert R_ZZ.equal(([1, 2, 3], 4), ([1, 2, 3], 4)) is True
    assert R_ZZ.equal(([1, 2, 3], 4), ([1, 2], None)) is False
    assert R_QQ.equal(([QQ(1, 2), QQ(3, 4)], 6), ([QQ(1, 2), QQ(3, 4)], None)) is False
    assert R_QQ.equal(([QQ(1, 2), QQ(3, 4)], None), ([QQ(1, 2)], None)) is False

def test_add():

    QQDUP1 = [QQ(1, 2), QQ(3, 4), QQ(1, 2), QQ(7, 8), QQ(2, 1)]
    QQDUP2 = [QQ(1, 3), QQ(2, 5), QQ(3, 4)]

    # Both are series
    assert R_ZZ.add(([1, 2, 3], 4), ([4, 5], 3)) == ([5, 7, 3], 3)
    assert R_QQ.add((QQDUP1, 6), (QQDUP2, 5)) == ([QQ(5, 6), QQ(23, 20),
        QQ(5, 4), QQ(7, 8), QQ(2, 1)], 5)

    # One is DUP and one series
    assert R_ZZ.add(([7, 14, 21, 28], 4), ([0, 7], None)) == ([7, 21, 21, 28], 4)
    assert R_QQ.add((QQDUP1, 5), (QQDUP2, None)) == ([QQ(5, 6), QQ(23, 20),
        QQ(5,4), QQ(7, 8), QQ(2, 1)], 5)

    # Both DUPs
    assert R_ZZ.add(([1, 3, 1, 2, 2], None), ([10, 20], None)) == ([11, 23, 1,
                                                                    2, 2], 5)
    assert R_QQ.add((QQDUP1, None), (QQDUP2, None)) == ([QQ(5, 6), QQ(23, 20),
        QQ(5, 4), QQ(7, 8), QQ(2, 1)], None)

def test_subtract():

    QQDUP1 = [QQ(3, 5), QQ(2, 3), QQ(4, 7), QQ(5, 9), QQ(1, 2)]
    QQDUP2 = [QQ(1, 4), QQ(1, 6), QQ(2, 5)]
    QQDUP3 = [QQ(5, 3), QQ(-2, 7), QQ(3, 11), QQ(7, 13)]

    # Both are series
    assert R_ZZ.subtract(([5, 8, 6], 4), ([2, 3], 3)) == ([3, 5, 6], 3)
    assert R_QQ.subtract((QQDUP1, 6), (QQDUP2, 5)) == ([QQ(7, 20), QQ(1, 2),
        QQ(6, 35), QQ(5, 9), QQ(1, 2)], 5)

    # One is DUP and one series
    assert R_ZZ.subtract(([9, 12, 15, 18], 4), ([3, 5], None)) == ([6, 7, 15, 18], 4)
    assert R_QQ.subtract((QQDUP1, 5), (QQDUP3, None)) == ([QQ(-16, 15), QQ(20, 21),
        QQ(23, 77), QQ(2,117), QQ(1, 2)], 5)

    # Both DUPs
    assert R_ZZ.subtract(([7, 4, 9, 2, 6], None), ([3, 2], None)) == ([4, 2, 9,
                                                                    2, 6], 5)
    assert R_QQ.subtract((QQDUP3, None), (QQDUP2, None)) == ([QQ(17, 12),
        QQ(-19, 42), QQ(-7, 55), QQ(7, 13)], None)

def test_multiply():

    QQDUP1 = [QQ(1, 2), QQ(3, 4), QQ(1, 2), QQ(7, 8), QQ(2, 1)]
    QQDUP2 = [QQ(1, 3), QQ(2, 5), QQ(3, 4)]

    # Both are series
    assert R_ZZ.multiply(([1, 2, 3], 4), ([4, 5], 3)) == ([4, 13, 22], 3)
    assert R_QQ.multiply((QQDUP1, 6), (QQDUP2, 5)) == ([QQ(1, 6), QQ(9, 20),
        QQ(101, 120), QQ(253, 240), QQ(167, 120)], 5)

    # One is DUP and one series
    assert R_ZZ.multiply(([7, 14, 21, 28], 4), ([0, 7], None)) == ([0, 49, 98,
                                                                    147], 4)
    assert R_QQ.multiply((QQDUP1, None), (QQDUP2, 3)) == ([QQ(1,6), QQ(9,20),
        QQ(101,120)], 3)

    # Both DUPs
    assert R_ZZ.multiply(([1, 3], None), ([10], None)) == ([10,30], None)
    assert R_QQ.multiply((QQDUP1, None), (QQDUP2, None)) == ([QQ(1, 6), QQ(9, 20),
        QQ(101, 120), QQ(253, 240), QQ(167, 120), QQ(233, 160), QQ(3, 2)], None)

def test_multiply_ground():
    assert R_ZZ.multiply_ground(([1, 2, 3], 4), 2) == ([2, 4, 6], 4)
    assert R_QQ.multiply_ground(([QQ(1, 2)], 6), QQ(3, 4)) == ([QQ(3, 8)], 6)

def test_trunc():
    assert R_ZZ.trunc(([1, 2, 3, 4, 5], 6), 3) == ([1, 2, 3], 3)
    assert R_QQ.trunc(([QQ(1, 2), QQ(3, 4), QQ(5, 6)], 6), 2) == ([QQ(1, 2),
                                                                    QQ(3, 4)], 2)

    raises(ValueError, lambda: R_ZZ.trunc(([1, 2], None), -1))
    raises(ValueError, lambda: R_QQ.trunc(([QQ(1, 2)], None), -1))

def test_special():
    R = PythonPowerSeriesRingZZ(5)

    p1 = ([1, 2, 3], 4)
    p2 = ([0, 0, 4, 5], 5)

    assert R.multiply(p1, p2) == ([0, 0, 4, 13, 22], 5)
