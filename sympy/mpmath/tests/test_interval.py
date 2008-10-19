from sympy.mpmath import *

def test_interval_identity():
    mp.dps = 15
    assert mpi(2) == mpi(2, 2)
    assert mpi(2) != mpi(-2, 2)
    assert mpi(-1, 1) == mpi(-1, 1)
    assert str(mpi('0.1')) == '[0.099999999999999991673, 0.10000000000000000555]'
    u = mpi(-1, 3)
    assert -1 in u
    assert 2 in u
    assert 3 in u
    assert -1.1 not in u
    assert 3.1 not in u
    assert mpi(-1, 3) in u
    assert mpi(0, 1) in u
    assert mpi(-1.1, 2) not in u
    assert mpi(2.5, 3.1) not in u
    w = mpi(-inf, inf)
    assert mpi(-5, 5) in w
    assert mpi(2, inf) in w
    assert mpi(0, 2) in mpi(0, 10)
    assert not (3 in mpi(-inf, 0))

def test_interval_arithmetic():
    mp.dps = 15
    assert mpi(2) + mpi(3,4) == mpi(5,6)
    assert mpi(1, 2)**2 == mpi(1, 4)
    assert mpi(1) + mpi(0, 1e-50) == mpi(1, mpf('1.0000000000000002'))
    x = 1 / (1 / mpi(3))
    assert x.a < 3 < x.b
    x = mpi(2) ** mpi(0.5)
    mp.dps += 5
    sq = sqrt(2)
    mp.dps -= 5
    assert x.a < sq < x.b
    assert mpi(1) / mpi(1, inf)
    assert mpi(2, 3) / inf == mpi(0, 0)
    assert mpi(0) / inf == 0
    assert mpi(0) / 0 == mpi(-inf, inf)
    assert mpi(inf) / 0 == mpi(-inf, inf)
    assert mpi(0) * inf == mpi(-inf, inf)
    assert 1 / mpi(2, inf) == mpi(0, 0.5)
    assert str((mpi(50, 50) * mpi(-10, -10)) / 3) == \
        '[-166.66666666666668561, -166.66666666666665719]'
    assert mpi(0, 4) ** 3 == mpi(0, 64)

def test_interval_mul():
    assert mpi(-1, 0) * inf == mpi(-inf, 0)
    assert mpi(-1, 0) * -inf == mpi(0, inf)
    assert mpi(0, 1) * inf == mpi(0, inf)
    assert mpi(0, 1) * mpi(0, inf) == mpi(0, inf)
    assert mpi(-1, 1) * inf == mpi(-inf, inf)
    assert mpi(-1, 1) * mpi(0, inf) == mpi(-inf, inf)
    assert mpi(-1, 1) * mpi(-inf, inf) == mpi(-inf, inf)
    assert mpi(-inf, 0) * mpi(0, 1) == mpi(-inf, 0)
    assert mpi(-inf, 0) * mpi(0, 0) * mpi(-inf, 0)
    assert mpi(-inf, 0) * mpi(-inf, inf) == mpi(-inf, inf)
    # Should be undefined?
    assert mpi(inf, inf) * 0 == mpi(-inf, inf)
    assert mpi(-inf, -inf) * 0 == mpi(-inf, inf)

def test_interval_pow():
    assert mpi(3)**2 == mpi(9, 9)
    assert mpi(-3)**2 == mpi(9, 9)
    assert mpi(-3, 1)**2 == mpi(0, 9)
    assert mpi(-3, -1)**2 == mpi(1, 9)
    assert mpi(-3, -1)**3 == mpi(-27, -1)
    assert mpi(-3, 1)**3 == mpi(-27, 1)
    assert mpi(-2, 3)**2 == mpi(0, 9)
    assert mpi(-3, 2)**2 == mpi(0, 9)
    assert mpi(4) ** -1 == mpi(0.25, 0.25)
    assert mpi(-4) ** -1 == mpi(-0.25, -0.25)
    assert mpi(4) ** -2 == mpi(0.0625, 0.0625)
    assert mpi(-4) ** -2 == mpi(0.0625, 0.0625)
    assert mpi(0, 1) ** inf == mpi(0, 1)
    assert mpi(0, 1) ** -inf == mpi(1, inf)
    assert mpi(0, inf) ** inf == mpi(0, inf)
    assert mpi(0, inf) ** -inf == mpi(0, inf)
    assert mpi(1, inf) ** inf == mpi(1, inf)
    assert mpi(1, inf) ** -inf == mpi(0, 1)
    assert mpi(2, 3) ** 1 == mpi(2, 3)
    assert mpi(2, 3) ** 0 == 1

def test_interval_sqrt():
    assert mpi(4) ** 0.5 == mpi(2)

def test_interval_div():
    assert mpi(0.5, 1) / mpi(-1, 0) == mpi(-inf, -0.5)
    assert mpi(0, 1) / mpi(0, 1) == mpi(0, inf)
    assert mpi(inf, inf) / mpi(inf, inf) == mpi(0, inf)
    assert mpi(inf, inf) / mpi(2, inf) == mpi(0, inf)
    assert mpi(inf, inf) / mpi(2, 2) == mpi(inf, inf)
    assert mpi(0, inf) / mpi(2, inf) == mpi(0, inf)
    assert mpi(0, inf) / mpi(2, 2) == mpi(0, inf)
    assert mpi(2, inf) / mpi(2, 2) == mpi(1, inf)
    assert mpi(2, inf) / mpi(2, inf) == mpi(0, inf)
    assert mpi(-4, 8) / mpi(1, inf) == mpi(-4, 8)
    assert mpi(-4, 8) / mpi(0.5, inf) == mpi(-8, 16)
    assert mpi(-inf, 8) / mpi(0.5, inf) == mpi(-inf, 16)
    assert mpi(-inf, inf) / mpi(0.5, inf) == mpi(-inf, inf)
    assert mpi(8, inf) / mpi(0.5, inf) == mpi(0, inf)
    assert mpi(-8, inf) / mpi(0.5, inf) == mpi(-16, inf)
    assert mpi(-4, 8) / mpi(inf, inf) == mpi(0, 0)
    assert mpi(0, 8) / mpi(inf, inf) == mpi(0, 0)
    assert mpi(0, 0) / mpi(inf, inf) == mpi(0, 0)
    assert mpi(-inf, 0) / mpi(inf, inf) == mpi(-inf, 0)
    assert mpi(-inf, 8) / mpi(inf, inf) == mpi(-inf, 0)
    assert mpi(-inf, inf) / mpi(inf, inf) == mpi(-inf, inf)
    assert mpi(-8, inf) / mpi(inf, inf) == mpi(0, inf)
    assert mpi(0, inf) / mpi(inf, inf) == mpi(0, inf)
    assert mpi(8, inf) / mpi(inf, inf) == mpi(0, inf)
    assert mpi(inf, inf) / mpi(inf, inf) == mpi(0, inf)
    assert mpi(-1, 2) / mpi(0, 1) == mpi(-inf, +inf)
    assert mpi(0, 1) / mpi(0, 1) == mpi(0.0, +inf)
    assert mpi(-1, 0) / mpi(0, 1) == mpi(-inf, 0.0)
    assert mpi(-0.5, -0.25) / mpi(0, 1) == mpi(-inf, -0.25)
    assert mpi(0.5, 1) / mpi(0, 1) == mpi(0.5, +inf)
    assert mpi(0.5, 4) / mpi(0, 1) == mpi(0.5, +inf)
    assert mpi(-1, -0.5) / mpi(0, 1) == mpi(-inf, -0.5)
    assert mpi(-4, -0.5) / mpi(0, 1) == mpi(-inf, -0.5)
    assert mpi(-1, 2) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(0, 1) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(-1, 0) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(-0.5, -0.25) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(0.5, 1) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(0.5, 4) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(-1, -0.5) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(-4, -0.5) / mpi(-2, 0.5) == mpi(-inf, +inf)
    assert mpi(-1, 2) / mpi(-1, 0) == mpi(-inf, +inf)
    assert mpi(0, 1) / mpi(-1, 0) == mpi(-inf, 0.0)
    assert mpi(-1, 0) / mpi(-1, 0) == mpi(0.0, +inf)
    assert mpi(-0.5, -0.25) / mpi(-1, 0) == mpi(0.25, +inf)
    assert mpi(0.5, 1) / mpi(-1, 0) == mpi(-inf, -0.5)
    assert mpi(0.5, 4) / mpi(-1, 0) == mpi(-inf, -0.5)
    assert mpi(-1, -0.5) / mpi(-1, 0) == mpi(0.5, +inf)
    assert mpi(-4, -0.5) / mpi(-1, 0) == mpi(0.5, +inf)
    assert mpi(-1, 2) / mpi(0.5, 1) == mpi(-2.0, 4.0)
    assert mpi(0, 1) / mpi(0.5, 1) == mpi(0.0, 2.0)
    assert mpi(-1, 0) / mpi(0.5, 1) == mpi(-2.0, 0.0)
    assert mpi(-0.5, -0.25) / mpi(0.5, 1) == mpi(-1.0, -0.25)
    assert mpi(0.5, 1) / mpi(0.5, 1) == mpi(0.5, 2.0)
    assert mpi(0.5, 4) / mpi(0.5, 1) == mpi(0.5, 8.0)
    assert mpi(-1, -0.5) / mpi(0.5, 1) == mpi(-2.0, -0.5)
    assert mpi(-4, -0.5) / mpi(0.5, 1) == mpi(-8.0, -0.5)
    assert mpi(-1, 2) / mpi(-2, -0.5) == mpi(-4.0, 2.0)
    assert mpi(0, 1) / mpi(-2, -0.5) == mpi(-2.0, 0.0)
    assert mpi(-1, 0) / mpi(-2, -0.5) == mpi(0.0, 2.0)
    assert mpi(-0.5, -0.25) / mpi(-2, -0.5) == mpi(0.125, 1.0)
    assert mpi(0.5, 1) / mpi(-2, -0.5) == mpi(-2.0, -0.25)
    assert mpi(0.5, 4) / mpi(-2, -0.5) == mpi(-8.0, -0.25)
    assert mpi(-1, -0.5) / mpi(-2, -0.5) == mpi(0.25, 2.0)
    assert mpi(-4, -0.5) / mpi(-2, -0.5) == mpi(0.25, 8.0)
    # Should be undefined?
    assert mpi(0, 0) / mpi(0, 0) == mpi(-inf, inf)
    assert mpi(0, 0) / mpi(0, 1) == mpi(-inf, inf)

def test_interval_cos_sin():
    mp.dps = 15
    # Around 0
    assert cos(mpi(0)) == 1
    assert sin(mpi(0)) == 0
    assert cos(mpi(0,1)) == mpi(0.54030230586813965399, 1.0)
    assert sin(mpi(0,1)) == mpi(0, 0.8414709848078966159)
    assert cos(mpi(1,2)) == mpi(-0.4161468365471424069, 0.54030230586813976501)
    assert sin(mpi(1,2)) == mpi(0.84147098480789650488, 1.0)
    assert sin(mpi(1,2.5)) == mpi(0.59847214410395643824, 1.0)
    assert cos(mpi(-1, 1)) == mpi(0.54030230586813965399, 1.0)
    assert cos(mpi(-1, 0.5)) == mpi(0.54030230586813965399, 1.0)
    assert cos(mpi(-1, 1.5)) == mpi(0.070737201667702906405, 1.0)
    assert sin(mpi(-1,1)) == mpi(-0.8414709848078966159, 0.8414709848078966159)
    assert sin(mpi(-1,0.5)) == mpi(-0.8414709848078966159, 0.47942553860420300538)
    assert sin(mpi(-1,1e-100)) == mpi(-0.8414709848078966159, 1.00000000000000002e-100)
    assert sin(mpi(-2e-100,1e-100)) == mpi(-2.00000000000000004e-100, 1.00000000000000002e-100)
    # Same interval
    assert cos(mpi(2, 2.5)) == mpi(-0.80114361554693380718, -0.41614683654714235139)
    assert cos(mpi(3.5, 4)) == mpi(-0.93645668729079634129, -0.65364362086361182946)
    assert cos(mpi(5, 5.5)) == mpi(0.28366218546322624627, 0.70866977429126010168)
    assert sin(mpi(2, 2.5)) == mpi(0.59847214410395654927, 0.90929742682568170942)
    assert sin(mpi(3.5, 4)) == mpi(-0.75680249530792831347, -0.35078322768961983646)
    assert sin(mpi(5, 5.5)) == mpi(-0.95892427466313856499, -0.70554032557039181306)
    # Higher roots
    mp.dps = 55
    w = 4*10**50 + mpf(0.5)
    for p in [15, 40, 80]:
        mp.dps = p
        assert 0 in sin(4*mpi(pi))
        assert 0 in sin(4*10**50*mpi(pi))
        assert 0 in cos((4+0.5)*mpi(pi))
        assert 0 in cos(w*mpi(pi))
        assert 1 in cos(4*mpi(pi))
        assert 1 in cos(4*10**50*mpi(pi))
    mp.dps = 15
