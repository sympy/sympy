# A collection of failing integrals from the issues.

from __future__ import division

from sympy import (
    integrate, Integral, exp, I, oo, pi, sign, sqrt, Rational, Symbol, sin, cos,
    tan, S, log, Function, gamma, sinh,
)

from sympy.utilities.pytest import XFAIL, skip

from sympy.abc import x, k, c, y, R, b, h, a, m, A, z, t

import signal, os

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

@XFAIL
def test_issue_459():
    # Si special function
    assert not integrate(cos(x*y), (x, -pi/2, pi/2), (y, 0, pi)).has(Integral)

@XFAIL
def test_issue_841a():
    assert not integrate(exp(-x**2)*exp(I*k*x), (x, -oo, oo)).has(Integral)

@XFAIL
def test_issue_841b():
    a = Symbol('a', positive=True)
    assert not integrate(c*exp(-x**2/a)*exp(I*k*x), (x, -oo, oo)).has(Integral)

@XFAIL
def test_issue_1001():
    # Note: issue has additional integral
    assert not integrate(sqrt(y**2-x**2), x).has(Integral)

@XFAIL
def test_issue_1020():
    assert not integrate(sqrt(y**2 - x**2)/x, x).has(Integral)

@XFAIL
def test_issue_1113():
    assert not integrate(sign(x), x).has(Integral)

@XFAIL
def test_issue_1116():
    # (singularity detection)
    assert integrate(1/(x**2), (x, -1, 1)) == oo

@XFAIL
def test_issue_1127a():
    # The integral doesn't calculate, but the real problem is an exception
    # caused the by Real (1/2 instead of S(1)/2) making it fail.
    assert not integrate(-(1 + (-x + x**2)**(1/2))/(-x + (1 + (-x + x**2)**(1/2))*x), x).has(Integral)

@XFAIL
def test_issue_1127b():
    assert not integrate(2*a/((4*a**2+x**2)*sqrt(4*a**2+x**2)),x).has(Integral)

@XFAIL
def test_issue_1127c():
    # The integral doesn't calculate, but the real problem is an exception
    # caused by the Real (3/2 instead of S(3)/2) making it fail.
    assert not integrate(2*a/((((2*a)**2+x**2))**(3/2)),x).has(Integral)

@XFAIL
def test_issue_1127d():
    # The same as test_issue_1127c with sympified exponent
    assert not integrate(2*a/((4*a**2+x**2)*sqrt(4*a**2+x**2)),x).has(Integral)

@XFAIL
def test_issue_1127e():
    # Similar to test_issue_841a above
    assert not integrate(-(8 + 2*sin(x) + 6*exp(x))*exp(2*x), x).has(Integral)

@XFAIL
def test_issue_1135():
    assert not integrate(1/sqrt(1+tan(x)**2)).has(Integral)

@XFAIL
def test_issue_1227():
    assert not (2*integrate(((h*(x-R+b))/b)*sqrt(R**2-x**2), (x, R-b, R))).has(Integral)

@XFAIL
def test_issue_1304a():
    assert not integrate(sqrt(x**2 + y**2), x).has(Integral)

@XFAIL
def test_issue_1304b():
    assert not integrate(1/(x**2 + y**2)**(Rational(3,2)),y).has(Integral)

@XFAIL
def test_issue_1323():
    assert not integrate(1/sqrt(16 + 4*x**2), x).has(Integral)

@XFAIL
def test_issue_1392():
    assert not integrate(x*sqrt(x**2+2*x+4), x).has(Integral)

@XFAIL
def test_issue_1393():
    assert not integrate(x**2 * sqrt(5-x**2), x).has(Integral)

@XFAIL
def test_issue_1394():
    assert not integrate(x*sqrt(1+2*x), x).has(Integral)

@XFAIL
def test_issue_1412():
    # This works, but gives a complicated answer.  The correct answer is x - cos(x).
    # The last one is what Maple gives.  It is also quite slow.
    t = 5 # Timeout time, requires ~220 sec.
    assert run_with_timeout("integrate(cos(x)**2 / (1-sin(x)))", t) in [x - cos(x),
        1 - cos(x) + x, -2/(tan((S(1)/2)*x)**2+1)+x]

@XFAIL
def test_issue_1415():
    # The correct answer is 2*sin(x)
    assert not integrate(sin(2*x)/ sin(x)).has(Integral)

@XFAIL
def test_issue_1418():
    # Currently give a traceback
    assert not integrate((x**Rational(1,2) - x**3)/x**Rational(1,3), x).has(Integral)

@XFAIL
def test_issue_1426():
    # Warning: takes a long time
    assert not integrate((x**m * (1 - x)**n * (a + b*x + c*x**2))/(1 + x**2), (x, 0, 1)).has(Integral)

@XFAIL
def test_issue_1441():
    # Note, this integral is probably non-elementary
    t = 5 # Timeout time
    assert not run_with_timeout("integrate((sin(1/x) - x*exp(x))/((-sin(1/x) + x*exp(x))*x + x*sin(1/x)), x)", t).has(Integral)

@XFAIL
def test_issue_1452():
    assert integrate(1/(x*sqrt(1-x**2)),x).has(Integral)

@XFAIL
def test_issue_1576():
    assert not integrate(4*pi**2*x**2*y**4*(y**2+9*x**2)/(y**2+x**2)**3, x).has(Integral)

@XFAIL
def test_issue_1604():
    g = Function('g')
    assert integrate(exp(x)*g(x), x).has(Integral)

@XFAIL
def test_issue_1638a():
    # Implementation of Si()
    assert integrate(sin(x)/x, x).has(Integral)

@XFAIL
def test_issue_1638b():
    assert integrate(sin(x)/x, (x, -oo, oo)) == pi/2

@XFAIL
def test_issue_1768():
    assert not integrate(sin(x)*tan(x), x).has(Integral)

@XFAIL
def test_issue_1791():
    # Note: Answer contains erf()
    assert not integrate(exp(-log(x)**2),x).has(Integral)

@XFAIL
def test_issue_1792():
    # Requires the hypergeometric function.
    t = 5 # Timeout
    assert not run_with_timeout("integrate(cos(x)**y, x)", t).has(Integral)

@XFAIL
def test_issue_1793():
    P1 = -A*exp(-z)
    P2 = -A/(c*t)*(sin(x)**2 + cos(y)**2)
    assert not integrate(c*(P2 - P1), t).has(Integral)


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
def test_issue_1842():
    t = 5 # Timeout time
    assert not run_with_timeout("integrate(sqrt(1+sinh(x/20)**2),(x,-25, 25))", t).has(Integral)

@XFAIL
def test_issue_1851():
    # Problem is with exception
    assert not integrate((-60*exp(x) - 19.2*exp(4*x))*exp(4*x), x).has(Integral)

@XFAIL
def test_issue_1869():
    t = 5 # timeout
    assert not run_with_timeout("integrate(sin(log(x**2)))", t).has(Integral)

@XFAIL
def test_issue_1893():
    # Non-elementary integral.  Requires hypergeometric/Meijer-G handling.
    t = 5 # Timeout
    assert not run_with_timeout("integrate(log(x) * x**(k-1) * exp(-x) / gamma(k), (x, 0, oo))", t).has(Integral)

@XFAIL
def test_issue_1888():
    f = Function('f')
    assert integrate(f(x).diff(x)**2, x).has(Integral)
