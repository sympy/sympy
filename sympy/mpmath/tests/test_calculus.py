from sympy.mpmath import *

def test_approximation():
    mp.dps = 15
    f = lambda x: cos(2-2*x)/x
    p, err = chebyfit(f, [2, 4], 8, error=True)
    assert err < 1e-5
    for i in range(10):
        x = 2 + i/5.
        assert abs(polyval(p, x) - f(x)) < err

def test_limits():
    mp.dps = 15
    assert limit(lambda x: (x-sin(x))/x**3, 0).ae(mpf(1)/6)
    assert limit(lambda n: (1+1/n)**n, inf).ae(e)

def test_polyval():
    assert polyval([], 3) == 0
    assert polyval([0], 3) == 0
    assert polyval([5], 3) == 5
    # 4x^3 - 2x + 5
    p = [4, 0, -2, 5]
    assert polyval(p,4) == 253
    assert polyval(p,4,derivative=True) == (253, 190)

def test_polyroots():
    p = polyroots([1,-4])
    assert p[0].ae(4)
    p, q = polyroots([1,2,3])
    assert p.ae(-1 - sqrt(2)*j)
    assert q.ae(-1 + sqrt(2)*j)
    #this is not a real test, it only tests a specific case
    assert polyroots([1]) == []
    try:
        polyroots([0])
        assert False
    except ValueError:
        pass

def test_pade():
    one = mpf(1)
    mp.dps = 20
    N = 10
    a = [one]
    k = 1
    for i in range(1, N+1):
        k *= i
        a.append(one/k)
    p, q = pade(a, N//2, N//2)
    for x in arange(0, 1, 0.1):
        r = polyval(p[::-1], x)/polyval(q[::-1], x)
        assert(r.ae(exp(x), 1.0e-10))
    mp.dps = 15

def test_fourier():
    mp.dps = 15
    c, s = fourier(lambda x: x+1, [-1, 2], 2)
    #plot([lambda x: x+1, lambda x: fourierval((c, s), [-1, 2], x)], [-1, 2])
    assert c[0].ae(1.5)
    assert c[1].ae(-3*sqrt(3)/(2*pi))
    assert c[2].ae(3*sqrt(3)/(4*pi))
    assert s[0] == 0
    assert s[1].ae(3/(2*pi))
    assert s[2].ae(3/(4*pi))
    assert fourierval((c, s), [-1, 2], 1).ae(1.9134966715663442)

def test_differint():
    mp.dps = 15
    assert differint(lambda t: t, 2, -0.5).ae(8*sqrt(2/pi)/3)
