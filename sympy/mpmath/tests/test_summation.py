from sympy.mpmath import *

def test_sumem():
    f = lambda x: x**-3
    fp = lambda x, k: (-1)**k * factorial(k+2)/2 * x**(-k-3)
    for prec in [15, 50]:
        s, err = sumem(f, [1, inf], fderiv=fp, error=1)
        assert s.ae(zeta(3))
    mp.dps = 15
    assert sumem(lambda k: k**4 + 3*k + 1, [10, 100], N=5).ae(2050333103)

def test_sumsh():
    mp.dps = 15
    assert sumsh(lambda k: (-1)**(k+1) / k, [1, inf]).ae(log(2))
    assert sumsh(lambda k: (-1)**(k+1) / k**2, [1, inf]).ae(pi**2 / 12)
    assert sumsh(lambda k: (-1)**k / log(k), [2, inf]).ae(0.9242998972229388)
    assert sumsh(lambda k: 1/factorial(k), [0, inf]).ae(e)

def test_sumrich():
    mp.dps = 15
    assert sumrich(lambda k: 1/k**2, [1, inf]).ae(pi**2 / 6)
