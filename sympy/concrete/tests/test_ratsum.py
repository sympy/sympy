from sympy import S, Symbol, Dummy, polygamma, pi, RootSum, Lambda, factor
from sympy import ratsum, ratsum_def


x = Symbol("x")
n = Symbol("n")


def test_ratsum_abramov():
    # Examples taken from
    # "On the Summation of Rational Functions"
    # by S. A. Abramov

    r = x**3
    R = ratsum(r, x)
    Rn = ratsum_def(r, (x, 1, n))
    assert R == x**4/4 - x**3/2 + x**2/4
    assert Rn == n**4/4 - n**3/2 + n**2/4

    r = (6*x+3)/(4*x**4+8*x**3+8*x**2+4*x+3)
    R = ratsum(r, x)
    Rn = ratsum_def(r, (x, 1, n))
    assert R == factor(-3/(2*(2*x**2 + 1)))
    assert Rn == S(1)/2 - 3/(2*(2*n**2 + 1))

    r = 1/(x**2+n**2-3*x+3*n-2*x*n+2)
    R = ratsum(r, x)
    Rn = ratsum_def(r, (x, 1, n))
    assert R == -1/(-n + x - 2)
    assert Rn == S(1)/2 + 1/(-n - 1)

    r = 1/x**2
    R = ratsum(r, x)
    Rn = ratsum_def(r, (x, 1, n))
    assert R == -polygamma(1, x)
    assert Rn == -polygamma(1, n) + pi**2/6


def test_ratsum_man():
    # Examples taken from
    # "On Computing Closed Forms for Indefinite Summations"
    # by Yiu-Kwong Man
    # doi : "10.1006/jsco.1993.1053"

    # First examples (Table 1)
    # 1
    r = 2**x * (x**2 - 2*x - 1) / (x**2 * (x+1)**2)
    R = ratsum(r, x)
    assert R == 2**x/x**2

    # 2
    r = (3*x**2 + 3*x + 1) / (x**3*(x+1)**3)
    R = ratsum(r, x)
    assert R == -1/x**3

    # 3
    r = (6*x + 3) / (4*x**4 + 8*x**3 + 8*x**2 + 4*x + 3)
    R = ratsum(r, x)
    assert R == factor(-3/(2*(2*x**2 + 1)))

    # 4
    r = -(x**2 + 3*x + 3) / (x**4 + 2*x**3 - 3*x**2 - 4*x + 2)
    R = ratsum(r, x)
    assert R == (x + 1)/(x**2 - 2)

    # 5
    r = x**2 * 4**x / ((x+1)*(x+2))
    R = ratsum(r, x)
    assert R == 4**x/3 - 4**x/(x + 1)

    # 6
    r = (2*x - 1)**3
    R = ratsum(r, x)
    assert R == 2*x**4 - 8*x**3 + 11*x**2 - 6*x

    # 7
    r = 3*x**2 + 3*x + 1
    R = ratsum(r, x)
    assert R == x**3

    # 8
    r = (x**2 + 4*x + 2) * 2**x
    R = ratsum(r, x)
    assert R == 2**x*x**2

    # 9
    r = 2**x * (x**3 - 3*x**2 - 3*x - 1) / (x**3 * (x+1)**3)
    R = ratsum(r, x)
    assert R == 2**x/x**3

    # 10
    k = Symbol("k")
    r = x * k**x
    R = ratsum(r, x)
    assert R == k**x*(-k/(k - 1)**2 + x/(k - 1))


    # More examples (Table 2)
    # 1 (rame as 2 above)
    # 2 (rame as 9 above)

    # 3
    r = 3**x * (2*x**3 + x**2 + 3*x + 6) / ((x**2+2) * (x**2 + 2*x + 3))
    R = ratsum(r, x)
    assert R == 3**x*x/(x**2 + 2)

    # 4
    r = 4 * (1-x) * (x**2 - 2*x - 1) / (x**2 * (x+1)**2 * (x-2)**2 * (x-3)**2)
    R = ratsum(r, x)
    assert R == 1/(x**4 - 6*x**3 + 9*x**2)

    # 5
    r = 2**x * (x**4 - 14*x**2 - 24*x -9) / (x**2 * (x+1)**2 * (x+2)**2 * (x+3)**2)
    R = ratsum(r, x)
    assert R == 2**x/(x**4 + 4*x**3 + 4*x**2)


def test_own():

    _w = Symbol("_w")
    _a = Symbol("_a")

    # An own example
    r = (1+x+x**2+x**3)/(x**7+x)
    R = ratsum(r, x)
    R == (polygamma(0, x) +
          RootSum(_w**4 - _w**2 + 1,
                  Lambda(_a, (7*_a**4/6 - _a**3/6 - 4*_a**2/3 - _a/6 + 1)*polygamma(0, -_a + x))))
