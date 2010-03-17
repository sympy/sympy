from sympy.mpmath import *

def test_type_compare():
    assert mpf(2) == mpc(2,0)
    assert mpf(0) == mpc(0)
    assert mpf(2) != mpc(2, 0.00001)
    assert mpf(2) == 2.0
    assert mpf(2) != 3.0
    assert mpf(2) == 2
    assert mpf(2) != '2.0'
    assert mpc(2) != '2.0'

def test_add():
    assert mpf(2.5) + mpf(3) == 5.5
    assert mpf(2.5) + 3 == 5.5
    assert mpf(2.5) + 3.0 == 5.5
    assert 3 + mpf(2.5) == 5.5
    assert 3.0 + mpf(2.5) == 5.5
    assert (3+0j) + mpf(2.5) == 5.5
    assert mpc(2.5) + mpf(3) == 5.5
    assert mpc(2.5) + 3 == 5.5
    assert mpc(2.5) + 3.0 == 5.5
    assert mpc(2.5) + (3+0j) == 5.5
    assert 3 + mpc(2.5) == 5.5
    assert 3.0 + mpc(2.5) == 5.5
    assert (3+0j) + mpc(2.5) == 5.5

def test_sub():
    assert mpf(2.5) - mpf(3) == -0.5
    assert mpf(2.5) - 3 == -0.5
    assert mpf(2.5) - 3.0 == -0.5
    assert 3 - mpf(2.5) == 0.5
    assert 3.0 - mpf(2.5) == 0.5
    assert (3+0j) - mpf(2.5) == 0.5
    assert mpc(2.5) - mpf(3) == -0.5
    assert mpc(2.5) - 3 == -0.5
    assert mpc(2.5) - 3.0 == -0.5
    assert mpc(2.5) - (3+0j) == -0.5
    assert 3 - mpc(2.5) == 0.5
    assert 3.0 - mpc(2.5) == 0.5
    assert (3+0j) - mpc(2.5) == 0.5

def test_mul():
    assert mpf(2.5) * mpf(3) == 7.5
    assert mpf(2.5) * 3 == 7.5
    assert mpf(2.5) * 3.0 == 7.5
    assert 3 * mpf(2.5) == 7.5
    assert 3.0 * mpf(2.5) == 7.5
    assert (3+0j) * mpf(2.5) == 7.5
    assert mpc(2.5) * mpf(3) == 7.5
    assert mpc(2.5) * 3 == 7.5
    assert mpc(2.5) * 3.0 == 7.5
    assert mpc(2.5) * (3+0j) == 7.5
    assert 3 * mpc(2.5) == 7.5
    assert 3.0 * mpc(2.5) == 7.5
    assert (3+0j) * mpc(2.5) == 7.5

def test_div():
    assert mpf(6) / mpf(3) == 2.0
    assert mpf(6) / 3 == 2.0
    assert mpf(6) / 3.0 == 2.0
    assert 6 / mpf(3) == 2.0
    assert 6.0 / mpf(3) == 2.0
    assert (6+0j) / mpf(3.0) == 2.0
    assert mpc(6) / mpf(3) == 2.0
    assert mpc(6) / 3 == 2.0
    assert mpc(6) / 3.0 == 2.0
    assert mpc(6) / (3+0j) == 2.0
    assert 6 / mpc(3) == 2.0
    assert 6.0 / mpc(3) == 2.0
    assert (6+0j) / mpc(3) == 2.0

def test_pow():
    assert mpf(6) ** mpf(3) == 216.0
    assert mpf(6) ** 3 == 216.0
    assert mpf(6) ** 3.0 == 216.0
    assert 6 ** mpf(3) == 216.0
    assert 6.0 ** mpf(3) == 216.0
    assert (6+0j) ** mpf(3.0) == 216.0
    assert mpc(6) ** mpf(3) == 216.0
    assert mpc(6) ** 3 == 216.0
    assert mpc(6) ** 3.0 == 216.0
    assert mpc(6) ** (3+0j) == 216.0
    assert 6 ** mpc(3) == 216.0
    assert 6.0 ** mpc(3) == 216.0
    assert (6+0j) ** mpc(3) == 216.0

def test_mixed_misc():
    assert 1 + mpf(3) == mpf(3) + 1 == 4
    assert 1 - mpf(3) == -(mpf(3) - 1) == -2
    assert 3 * mpf(2) == mpf(2) * 3 == 6
    assert 6 / mpf(2) == mpf(6) / 2 == 3
    assert 1.0 + mpf(3) == mpf(3) + 1.0 == 4
    assert 1.0 - mpf(3) == -(mpf(3) - 1.0) == -2
    assert 3.0 * mpf(2) == mpf(2) * 3.0 == 6
    assert 6.0 / mpf(2) == mpf(6) / 2.0 == 3

def test_add_misc():
    mp.dps = 15
    assert mpf(4) + mpf(-70) == -66
    assert mpf(1) + mpf(1.1)/80 == 1 + 1.1/80
    assert mpf((1, 10000000000)) + mpf(3) == mpf((1, 10000000000))
    assert mpf(3) + mpf((1, 10000000000)) == mpf((1, 10000000000))
    assert mpf((1, -10000000000)) + mpf(3) == mpf(3)
    assert mpf(3) + mpf((1, -10000000000)) == mpf(3)
    assert mpf(1) + 1e-15 != 1
    assert mpf(1) + 1e-20 == 1
    assert mpf(1.07e-22) + 0 == mpf(1.07e-22)
    assert mpf(0) + mpf(1.07e-22) == mpf(1.07e-22)

def test_complex_misc():
    # many more tests needed
    assert 1 + mpc(2) == 3
    assert not mpc(2).ae(2 + 1e-13)
    assert mpc(2+1e-15j).ae(2)

def test_complex_zeros():
    for a in [0,2]:
      for b in [0,3]:
        for c in [0,4]:
          for d in [0,5]:
            assert mpc(a,b)*mpc(c,d) == complex(a,b)*complex(c,d)

def test_hash():
    for i in range(-256, 256):
        assert hash(mpf(i)) == hash(i)
    assert hash(mpf(0.5)) == hash(0.5)
    assert hash(mpc(2,3)) == hash(2+3j)
    # Check that this doesn't fail
    assert hash(inf)
    # Check that overflow doesn't assign equal hashes to large numbers
    assert hash(mpf('1e1000')) != hash('1e10000')
    assert hash(mpc(100,'1e1000')) != hash(mpc(200,'1e1000'))

def test_arithmetic_functions():
    import operator
    ops = [(operator.add, fadd), (operator.sub, fsub), (operator.mul, fmul),
        (operator.div, fdiv)]
    a = mpf(0.27)
    b = mpf(1.13)
    c = mpc(0.51+2.16j)
    d = mpc(1.08-0.99j)
    for x in [a,b,c,d]:
        for y in [a,b,c,d]:
            for op, fop in ops:
                if fop is not fdiv:
                    mp.prec = 200
                    z0 = op(x,y)
                mp.prec = 60
                z1 = op(x,y)
                mp.prec = 53
                z2 = op(x,y)
                assert fop(x, y, prec=60) == z1
                assert fop(x, y) == z2
                if fop is not fdiv:
                    assert fop(x, y, prec=inf) == z0
                    assert fop(x, y, dps=inf) == z0
                    assert fop(x, y, exact=True) == z0
                assert fneg(fneg(z1, exact=True), prec=inf) == z1
                assert fneg(z1) == -(+z1)
    mp.dps = 15
