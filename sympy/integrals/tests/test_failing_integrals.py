# A collection of failing integrals from the issues.

from __future__ import division

from sympy import (
    integrate, Integral, exp, I, oo, pi, sign, sqrt, Rational, Symbol, sin, cos,
    tan, S, log, Function, gamma, sinh,
)

from sympy.utilities.pytest import XFAIL, SKIP, slow

from sympy.abc import x, k, c, y, R, b, h, a, m, A, z, t

import signal


class TimeOutError(Exception):
    pass


def timeout(signum, frame, time):
    raise TimeOutError("Timed out after %d seconds" % time)


def run_with_timeout(test, time):
    # Set the signal handler and a 5-second alarm
    signal.signal(signal.SIGALRM, lambda s, f: timeout(s, f, time))
    signal.alarm(time)
    r = eval(test)
    signal.alarm(0)          # Disable the alarm
    return r


@SKIP("Too slow for @slow")
@XFAIL
def test_issue_781():
    # integrate_hyperexponential(Poly(t*2*(1 - t0**2)*t0*(x**3 + x**2), t), Poly((1 + t0**2)**2*2*(x**2 + x + 1), t), [Poly(1, x), Poly(1 + t0**2, t0), Poly(t, t)], [x, t0, t], [exp, tan])
    assert not integrate(exp(x)*cos(2*x)*sin(2*x) * (x**3 + x**2)/(2*(x**2 + x + 1)), x).has(Integral)


@XFAIL
def test_issue_1113():
    assert not integrate(sign(x), x).has(Integral)


@XFAIL
def test_issue_1135():
    assert not integrate(1/sqrt(1 + tan(x)**2)).has(Integral)


@XFAIL
def test_issue_1227():
    assert integrate(((h*(x - R + b))/b)*sqrt(R**2 - x**2), (x, R - b, R)).has(Integral)


@XFAIL
def test_issue_1392():
    assert not integrate(x*sqrt(x**2 + 2*x + 4), x).has(Integral)


@XFAIL
def test_issue_1393():
    assert not integrate(x**2 * sqrt(5 - x**2), x).has(Integral)


@XFAIL
@slow
def test_issue_1412():
    # This works, but gives a complicated answer.  The correct answer is x - cos(x).
    # The last one is what Maple gives.  It is also quite slow.
    assert integrate(cos(x)**2 / (1 - sin(x))) in [x - cos(x), 1 - cos(x) + x,
            -2/(tan((S(1)/2)*x)**2 + 1) + x]


@XFAIL
def test_issue_1415():
    # The correct answer is 2*sin(x)
    assert not integrate(sin(2*x)/ sin(x)).has(Integral)


@XFAIL
def test_issue_1426():
    # Warning: takes a long time
    assert not integrate((x**m * (1 - x)**n * (a + b*x + c*x**2))/(1 + x**2), (x, 0, 1)).has(Integral)


@XFAIL
@slow
def test_issue_1441():
    # Note, this integral is probably nonelementary
    assert not integrate(
        (sin(1/x) - x*exp(x)) /
        ((-sin(1/x) + x*exp(x))*x + x*sin(1/x)), x).has(Integral)


@XFAIL
def test_issue_1452():
    assert integrate(1/(x*sqrt(1 - x**2)), x).has(Integral)


@XFAIL
def test_issue_1638a():
    # Implementation of Si()
    assert integrate(sin(x)/x, x).has(Integral)


@XFAIL
def test_issue_1638b():
    assert integrate(sin(x)/x, (x, -oo, oo)) == pi/2


@XFAIL
@slow
def test_issue_1792():
    # Requires the hypergeometric function.
    assert not integrate(cos(x)**y, x).has(Integral)


@XFAIL
def test_issue_1796a():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), x).has(Integral)


@XFAIL
def test_issue_1796b():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, -oo, 0)).has(Integral)


@XFAIL
def test_issue_1796c():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, -oo, oo)).has(Integral)


@XFAIL
def test_issue_1796d():
    assert not integrate(exp(2*b*x)*exp(-a*x**2), (x, 0, oo)).has(Integral)


@XFAIL
@slow
def test_issue_1842():
    assert not integrate(sqrt(1 + sinh(x/20)**2), (x, -25, 25)).has(Integral)


@XFAIL
@slow
def test_issue_1851():
    # Problem is with exception
    assert not integrate((-60*exp(x) - 19.2*exp(4*x))*exp(4*x), x).has(Integral)


@XFAIL
@slow
def test_issue_1869():
    assert not integrate(sin(log(x**2))).has(Integral)


@XFAIL
@slow
def test_issue_1893():
    # Nonelementary integral.  Requires hypergeometric/Meijer-G handling.
    assert not integrate(log(x) * x**(k - 1) * exp(-x) / gamma(k), (x, 0, oo)).has(Integral)
