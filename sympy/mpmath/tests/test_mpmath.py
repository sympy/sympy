from sympy.mpmath.libmp import *
from sympy.mpmath import *
import random


#----------------------------------------------------------------------------
# Low-level tests
#

# Advanced rounding test
def test_add_rounding():
    mp.dps = 15
    a = from_float(1e-50)
    assert mpf_sub(mpf_add(fone, a, 53, round_up), fone, 53, round_up) == from_float(2.2204460492503131e-16)
    assert mpf_sub(fone, a, 53, round_up) == fone
    assert mpf_sub(fone, mpf_sub(fone, a, 53, round_down), 53, round_down) == from_float(1.1102230246251565e-16)
    assert mpf_add(fone, a, 53, round_down) == fone

def test_almost_equal():
    assert mpf(1.2).ae(mpf(1.20000001), 1e-7)
    assert not mpf(1.2).ae(mpf(1.20000001), 1e-9)
    assert not mpf(-0.7818314824680298).ae(mpf(-0.774695868667929))


#----------------------------------------------------------------------------
# Test basic arithmetic
#

# Test that integer arithmetic is exact
def test_aintegers():
    # XXX: re-fix this so that all operations are tested with all rounding modes
    random.seed(0)
    for prec in [6, 10, 25, 40, 100, 250, 725]:
      for rounding in ['d', 'u', 'f', 'c', 'n']:
        mp.dps = prec
        M = 10**(prec-2)
        M2 = 10**(prec//2-2)
        for i in range(10):
            a = random.randint(-M, M)
            b = random.randint(-M, M)
            assert mpf(a, rounding=rounding) == a
            assert int(mpf(a, rounding=rounding)) == a
            assert int(mpf(str(a), rounding=rounding)) == a
            assert mpf(a) + mpf(b) == a + b
            assert mpf(a) - mpf(b) == a - b
            assert -mpf(a) == -a
            a = random.randint(-M2, M2)
            b = random.randint(-M2, M2)
            assert mpf(a) * mpf(b) == a*b
            assert mpf_mul(from_int(a), from_int(b), mp.prec, rounding) == from_int(a*b)
    mp.dps = 15

def test_odd_int_bug():
    assert to_int(from_int(3), round_nearest) == 3

def test_str_1000_digits():
    mp.dps = 1001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '9518488472'[:9]
    assert str(pi)[-10:-1] == '2164201989'[:9]
    mp.dps = 15

def test_str_10000_digits():
    mp.dps = 10001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '5873258351'[:9]
    assert str(pi)[-10:-1] == '5256375678'[:9]
    mp.dps = 15

def test_monitor():
    f = lambda x: x**2
    a = []
    b = []
    g = monitor(f, a.append, b.append)
    assert g(3) == 9
    assert g(4) == 16
    assert a[0] == ((3,), {})
    assert b[0] == 9

def test_nint_distance():
    nint_distance(mpf(-3)) == (-3, -inf)
    nint_distance(mpc(-3)) == (-3, -inf)
    nint_distance(mpf(-3.1)) == (-3, -3)
    nint_distance(mpf(-3.01)) == (-3, -6)
    nint_distance(mpf(-3.001)) == (-3, -9)
    nint_distance(mpf(-3.0001)) == (-3, -13)
    nint_distance(mpf(-2.9)) == (-3, -3)
    nint_distance(mpf(-2.99)) == (-3, -6)
    nint_distance(mpf(-2.999)) == (-3, -9)
    nint_distance(mpf(-2.9999)) == (-3, -13)
    nint_distance(mpc(-3+0.1j)) == (-3, -3)
    nint_distance(mpc(-3+0.01j)) == (-3, -6)
    nint_distance(mpc(-3.1+0.1j)) == (-3, -3)
    nint_distance(mpc(-3.01+0.01j)) == (-3, -6)
    nint_distance(mpc(-3.001+0.001j)) == (-3, -9)
    nint_distance(mpf(0)) == (0, -inf)
    nint_distance(mpf(0.01)) == (0, -6)
    nint_distance(mpf('1e-100')) == (0, -332)
