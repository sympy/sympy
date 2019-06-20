# A collection of failing integrals from the issues.

from sympy import (
    integrate, Integral, exp, oo, pi, sign, sqrt, sin, cos,
    tan, S, log, gamma, sinh, sec, zeta, acos, atan
)

from sympy.utilities.pytest import XFAIL, SKIP, slow, skip, ON_TRAVIS

from sympy.abc import x, k, c, y, R, b, h, a, m


@SKIP("Too slow for @slow")
@XFAIL
def test_issue_3880():
    # integrate_hyperexponential(Poly(t*2*(1 - t0**2)*t0*(x**3 + x**2), t), Poly((1 + t0**2)**2*2*(x**2 + x + 1), t), [Poly(1, x), Poly(1 + t0**2, t0), Poly(t, t)], [x, t0, t], [exp, tan])
    assert not integrate(exp(x)*cos(2*x)*sin(2*x) * (x**3 + x**2)/(2*(x**2 + x + 1)), x).has(Integral)


@XFAIL
def test_issue_4212():
    assert not integrate(sign(x), x).has(Integral)


@XFAIL
def test_issue_4326():
    assert integrate(((h*(x - R + b))/b)*sqrt(R**2 - x**2), (x, R - b, R)).has(Integral)


@XFAIL
def test_issue_4491():
    assert not integrate(x*sqrt(x**2 + 2*x + 4), x).has(Integral)


@XFAIL
def test_issue_4511():
    # This works, but gives a complicated answer.  The correct answer is x - cos(x).
    # The last one is what Maple gives.  It is also quite slow.
    assert integrate(cos(x)**2 / (1 - sin(x))) in [x - cos(x), 1 - cos(x) + x,
            -2/(tan((S(1)/2)*x)**2 + 1) + x]


@XFAIL
def test_issue_4525():
    # Warning: takes a long time
    assert not integrate((x**m * (1 - x)**n * (a + b*x + c*x**2))/(1 + x**2), (x, 0, 1)).has(Integral)



@XFAIL
@slow
def test_issue_4540():
    if ON_TRAVIS:
        skip("Too slow for travis.")
    # Note, this integral is probably nonelementary
    assert not integrate(
        (sin(1/x) - x*exp(x)) /
        ((-sin(1/x) + x*exp(x))*x + x*sin(1/x)), x).has(Integral)


@XFAIL
def test_issue_4551():
    assert integrate(1/(x*sqrt(1 - x**2)), x).has(Integral)


@XFAIL
def test_issue_4737a():
    # Implementation of Si()
    assert integrate(sin(x)/x, x).has(Integral)


@XFAIL
def test_issue_1638b():
    assert integrate(sin(x)/x, (x, -oo, oo)) == pi/2


@XFAIL
@slow
def test_issue_4891():
    # Requires the hypergeometric function.
    assert not integrate(cos(x)**y, x).has(Integral)


@XFAIL
@slow
def test_issue_1796a():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), x).has(Integral)


@XFAIL
def test_issue_4895b():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, -oo, 0)).has(Integral)


@XFAIL
def test_issue_4895c():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, -oo, oo)).has(Integral)


@XFAIL
def test_issue_4895d():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, 0, oo)).has(Integral)


@XFAIL
@slow
def test_issue_4941():
    if ON_TRAVIS:
        skip("Too slow for travis.")
    assert not integrate(sqrt(1 + sinh(x/20)**2), (x, -25, 25)).has(Integral)


@XFAIL
def test_issue_4992():
    # Nonelementary integral.  Requires hypergeometric/Meijer-G handling.
    assert not integrate(log(x) * x**(k - 1) * exp(-x) / gamma(k), (x, 0, oo)).has(Integral)

@XFAIL
def test_issue_16084():
    assert not integrate(log(sin(x)), (x, 0, pi/2)).has(Integral)
    assert integrate(log(sin(x)), (x, 0, pi/2)) == -pi*log(2)/2

@XFAIL
def test_issue_16161():
    assert not integrate(x*sec(x)**2, x).has(Integral)
    assert integrate(x*sec(x)**2, x) == x*tan(x) + log(cos(x))


@XFAIL
def test_issue_15925a():
    assert not integrate(sqrt((1+sin(x))**2+(cos(x))**2), (x, -pi/2, pi/2)).has(Integral)


@XFAIL
@slow
def test_issue_15925b():
    assert not integrate(sqrt((-12*cos(x)**2*sin(x))**2+(12*cos(x)*sin(x)**2)**2),
                         (x, 0, pi/6)).has(Integral)


@XFAIL
def test_issue_15925b_manual():
    assert not integrate(sqrt((-12*cos(x)**2*sin(x))**2+(12*cos(x)*sin(x)**2)**2),
                         (x, 0, pi/6), manual=True).has(Integral)


@XFAIL
@slow
def test_issue_15227():
    if ON_TRAVIS:
        skip("Too slow for travis.")
    assert not integrate(log(1-x)*log((1+x)**2)/x, (x, 0, 1)).has(Integral)
    assert integrate(log(1-x)*log((1+x)**2)/x, (x, 0, 1)) == -5*zeta(3)/4


@XFAIL
def test_issue_14716():
    assert not integrate(log(x + 5)*cos(pi*x),(x, S.Half, 1)).has(Integral)
    # Mathematica can not solve it either, but
    # integrate(log(x + 5)*cos(pi*x),(x, S.Half, 1)).transform(x, y - 5).doit()
    # works
    assert integrate(log(x + 5)*cos(pi*x),(x, S.Half, 1)) == \
        -log(S(11)/2)/pi - Si(11*pi/2)/pi + Si(6*pi)/pi


@XFAIL
def test_issue_14709a():
    assert not integrate(x*acos(1 - 2*x/h), (x, 0, h)).has(Integral)
    assert integrate(x*acos(1 - 2*x/h), (x, 0, h)) == 5*h**2*pi/16


@XFAIL
def test_issue_14709b():
    assert not integrate(x*acos(1 - 2*x/21323), (x, 0, 21323)).has(Integral)


@XFAIL
def test_issue_14398():
    assert not integrate(exp(x**2)*cos(x), x).has(Integral)


@XFAIL
def test_issue_14074():
    assert not integrate(log(sin(x)), (x, 0, pi/2)).has(Integral)
    assert integrate(log(sin(x)), (x, 0, pi/2)) == -pi*log(2)/2


@XFAIL
@slow
def test_issue_14078b():
    assert not integrate((atan(4*x)-atan(2*x))/x, (x, 0, oo)).has(Integral)
    assert integrate((atan(4*x)-atan(2*x))/x, (x, 0, oo)) == pi*log(2)/2
