from sympy.mpmath import *

def test_approximation():
    mp.dps = 15
    f = lambda x: cos(2-2*x)/x
    p, err = chebyfit(f, 2, 4, 8, error=True)
    assert err < 1e-5
    for i in range(10):
        x = 2 + i/5.
        assert abs(polyval(p, x) - f(x)) < err

def test_limits():
    mp.dps = 15
    assert limit(lambda x: (x-sin(x))/x**3, 0).ae(mpf(1)/6)
    assert limit(lambda n: (1+1/n)**n, inf).ae(e)

def test_polyroots():
    #this is not a real test, it only tests a specific case
    assert polyroots([1]) == []
    try:
        polyroots([0])
        assert False
    except ValueError:
        pass

