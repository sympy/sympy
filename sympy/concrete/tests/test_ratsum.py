from sympy import S, Symbol, Dummy, polygamma, pi, RootSum, Lambda
from sympy import ratsum, ratsum_def


x = Symbol("x")
n = Symbol("n")


def test_ratsum_abramov():
    # Examples taken from
    # "On the Summation of Rational Functions"
    # by S. A. Abramov

    _c0 = Dummy("c0")

    s = x**3
    S = ratsum(s, x)
    Sn = ratsum_def(s, (x, 1, n))
    assert S == _c0 + x**4/4 - x**3/2 + x**2/4
    assert Sn == n**4/4 - n**3/2 + n**2/4

    s = (6*x+3)/(4*x**4+8*x**3+8*x**2+4*x+3)
    S = ratsum(s, x)
    Sn = ratsum_def(s, (x, 1, n))
    assert S ==  -3/(2*(2*x**2 + 1))
    assert Sn == S(1)/2 - 3/(2*(2*n**2 + 1))

    s = 1/(x**2+n**2-3*x+3*n-2*x*n+2)
    S = ratsum(s, x)
    Sn = ratsum_def(s, (x, 1, n))
    assert S == -1/(-n + x - 2)
    assert Sn == S(1)/2 + 1/(-n - 1)

    s = 1/x**2
    S = ratsum(s, x)
    Sn = ratsum_def(s, (x, 1, n))
    assert S == -polygamma(1, x)
    assert Sn == -polygamma(1, n) + pi**2/6


def test_ratsum_man():
    # Examples taken from
    # "On Computing Closed Forms for Indefinite Summations"
    # by Yiu-Kwong Man
    # doi : "10.1006/jsco.1993.1053"

    _c0 = Dummy("c0")

    # First examples (Table 1)
    # 1
    s = 2**x * (x**2 - 2*x - 1) / (x**2 * (x+1)**2)
    S = ratsum(s, x)
    assert S == 2**x/x**2

    # 2
    s = (3*x**2 + 3*x + 1) / (x**3*(x+1)**3)
    S = ratsum(s, x)
    assert S == -1/x**3

    # 3
    s = (6*x + 3) / (4*x**4 + 8*x**3 + 8*x**2 + 4*x + 3)
    S = ratsum(s, x)
    assert S == -3/(2*(2*x**2 + 1))

    # 4
    s = -(x**2 + 3*x + 3) / (x**4 + 2*x**3 - 3*x**2 - 4*x + 2)
    S = ratsum(s, x)
    assert S == (x + 1)/(x**2 - 2)

    # 5
    s = x**2 * 4**x / ((x+1)*(x+2))
    S = ratsum(s, x)
    assert S == 4**x/3 - 4**x/(x + 1)

    # 6
    s = (2*x - 1)**3
    S = ratsum(s, x)
    assert S == _c0 + 2*x**4 - 8*x**3 + 11*x**2 - 6*x

    # 7
    s = 3*x**2 + 3*x + 1
    S = ratsum(s, x)
    assert S == _c0 + x**3

    # 8
    s = (x**2 + 4*x + 2) * 2**x
    S = ratsum(s, x)
    assert S == 2**x*x**2

    # 9
    s = 2**x * (x**3 - 3*x**2 - 3*x - 1) / (x**3 * (x+1)**3)
    S = ratsum(s, x)
    assert S == 2**x/x**3

    # 10
    s = x * k**x
    S = ratsum(s, x)
    assert S == k**x*(-k/(k - 1)**2 + x/(k - 1))


    # More examples (Table 2)
    # 1 (same as 2 above)
    # 2 (same as 9 above)

    # 3
    s = 3**x * (2*x**3 + x**2 + 3*x + 6) / ((x**2+2) * (x**2 + 2*x + 3))
    S = ratsum(s, x)
    assert S == 3**x*x/(x**2 + 2)

    # 4
    s = 4 * (1-x) * (x**2 - 2*x - 1) / (x**2 * (x+1)**2 * (x-2)**2 * (x-3)**2)
    S = ratsum(s, x)
    assert S == 1/(x**4 - 6*x**3 + 9*x**2)

    # 5
    s = 2**x * (x**4 - 14*x**2 - 24*x -9) / (x**2 * (x+1)**2 * (x+2)**2 * (x+3)**2)
    S = ratsum(s, x)
    assert S == 2**x/(x**4 + 4*x**3 + 4*x**2)


def test_own():

    _w = Symbol("_w")
    _a = Symbol("_a")

    # An own example
    s = (1+x+x**2+x**3)/(x**7+x)
    S = ratsum(s, x)
    S == (polygamma(0, x) +
          RootSum(_w**4 - _w**2 + 1,
                  Lambda(_a, (7*_a**4/6 - _a**3/6 - 4*_a**2/3 - _a/6 + 1)*polygamma(0, -_a + x))))
