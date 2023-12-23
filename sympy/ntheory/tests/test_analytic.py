from sympy.ntheory.analytic_ntheory import _ramanujan_tau, ramanujan_tau
from sympy.testing.pytest import raises


def test__ramanujan_tau():
    A000594 = [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643,
               -115920, 534612, -370944, -577738, 401856, 1217160,
               987136, -6905934, 2727432, 10661420, -7109760,
               -4219488, -12830688, 18643272, 21288960, -25499225,
               13865712, -73279080, 24647168]
    for n in range(1, len(A000594) + 1):
        assert _ramanujan_tau(n) == A000594[n - 1]


def test_ramanujan_tau():
    raises(ValueError, lambda: ramanujan_tau(0))
    assert ramanujan_tau(1) == 1
    assert ramanujan_tau(2) == -24
    assert ramanujan_tau(2**10) == 36697722069188608
    assert ramanujan_tau(3*5*7*11*13) == 6294721545123699018240
    assert ramanujan_tau(101**2) == -4474772103901920904697
    assert ramanujan_tau(1999) == -1159913672832202000 # 1999 is a prime
