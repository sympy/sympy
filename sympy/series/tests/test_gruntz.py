from sympy import Symbol, exp, log, oo, Rational, I, sin, gamma, loggamma, S, \
    atan, acot, pi, cancel, E
from sympy.series.gruntz import compare, mrv, rewrite, mrv_leadterm, gruntz, \
    sign
from sympy.utilities.pytest import XFAIL, skip

"""
This test suite is testing the limit algorithm using the bottom up approach.
See the documentation in limits2.py. The algorithm itself is highly recursive
by nature, so "compare" is logically the lowest part of the algorithm, yet in
some sense it's the most complex part, because it needs to calculate a limit to
return the result.

Nevertheless the rest of the algorithm depends on compare that it works
correctly.
"""

x = Symbol('x', real=True)
m = Symbol('m', real=True)


def test_compare1():
    assert compare(2, x, x) == "<"
    assert compare(x, exp(x), x) == "<"
    assert compare(exp(x), exp(x**2), x) == "<"
    assert compare(exp(x**2),exp(exp(x)), x) == "<"
    assert compare(1,exp(exp(x)), x) == "<"

    assert compare(x, 2, x) == ">"
    assert compare(exp(x), x, x) == ">"
    assert compare(exp(x**2), exp(x), x) == ">"
    assert compare(exp(exp(x)), exp(x**2), x) == ">"
    assert compare(exp(exp(x)), 1, x) == ">"

    assert compare(2, 3, x) == "="
    assert compare(3, -5, x) == "="
    assert compare(2, -5, x) == "="

    assert compare(x, x**2, x) == "="
    assert compare(x**2, x**3, x) == "="
    assert compare(x**3, 1/x, x) == "="
    assert compare(1/x, x**m, x) == "="
    assert compare(x**m, -x, x) == "="

    assert compare(exp(x), exp(-x), x) == "="
    assert compare(exp(-x), exp(2*x), x) == "="
    assert compare(exp(2*x), exp(x)**2, x) == "="
    assert compare(exp(x)**2, exp(x+exp(-x)), x) == "="
    assert compare(exp(x), exp(x+exp(-x)), x) == "="

    assert compare(exp(x**2), 1/exp(x**2), x) == "="

def test_compare2():
    assert compare(exp(x),x**5,x) == ">"
    assert compare(exp(x**2),exp(x)**2,x) == ">"
    assert compare(exp(x),exp(x+exp(-x)),x) == "="
    assert compare(exp(x+exp(-x)),exp(x),x) == "="
    assert compare(exp(x+exp(-x)),exp(-x),x) == "="
    assert compare(exp(-x),x,x) ==  ">"
    assert compare(x,exp(-x),x) ==  "<"
    assert compare(exp(x+1/x),x,x) == ">"
    assert compare(exp(-exp(x)),exp(x),x) == ">"
    assert compare(exp(exp(-exp(x))+x),exp(-exp(x)),x) == "<"

def test_compare3():
    assert compare(exp(exp(x)),exp(x+exp(-exp(x))),x) == ">"

def test_sign1():
    assert sign(Rational(0), x) == 0
    assert sign(Rational(3), x) == 1
    assert sign(Rational(-5), x) == -1
    assert sign(log(x), x) == 1
    assert sign(exp(-x), x) == 1
    assert sign(exp(x), x) == 1
    assert sign(-exp(x), x) == -1
    assert sign(3-1/x, x) == 1
    assert sign(-3-1/x, x) == -1
    assert sign(sin(1/x), x) == 1

def test_sign2():
    assert sign(x, x) == 1
    assert sign(-x, x) == -1
    y = Symbol("y", positive=True)
    assert sign(y, x) == 1
    assert sign(-y, x) == -1
    assert sign(y*x, x) == 1
    assert sign(-y*x, x) == -1

def mmrv(a, b): return set(mrv(a, b)[0].keys())

def test_mrv1():
    assert mmrv(x, x) == set([x])
    assert mmrv(x+1/x, x) == set([x])
    assert mmrv(x**2, x) == set([x])
    assert mmrv(log(x), x) == set([x])
    assert mmrv(exp(x), x) == set([exp(x)])
    assert mmrv(exp(-x), x) == set([exp(-x)])
    assert mmrv(exp(x**2), x) == set([exp(x**2)])
    assert mmrv(-exp(1/x), x) == set([x])
    assert mmrv(exp(x+1/x), x) == set([exp(x+1/x)])

def test_mrv2a():
    assert mmrv(exp(x+exp(-exp(x))), x) == set([exp(-exp(x))])
    assert mmrv(exp(x+exp(-x)), x) == set([exp(x+exp(-x)), exp(-x)])
    assert mmrv(exp(1/x+exp(-x)), x) == set([exp(-x)])

#sometimes infinite recursion due to log(exp(x**2)) not simplifying
def test_mrv2b():
    assert mmrv(exp(x+exp(-x**2)), x) == set([exp(-x**2)])

#sometimes infinite recursion due to log(exp(x**2)) not simplifying
def test_mrv2c():
    assert mmrv(exp(-x+1/x**2)-exp(x+1/x), x) == set([exp(x+1/x), exp(1/x**2-x)])

#sometimes infinite recursion due to log(exp(x**2)) not simplifying
def test_mrv3():
    assert mmrv(exp(x**2)+x*exp(x)+log(x)**x/x, x) == set([exp(x**2)])
    assert mmrv(exp(x)*(exp(1/x+exp(-x))-exp(1/x)), x) == set([exp(x), exp(-x)])
    assert mmrv(log(x**2+2*exp(exp(3*x**3*log(x)))), x) == set([exp(exp(3*x**3*log(x)))])
    assert mmrv(log(x-log(x))/log(x), x) == set([x])
    assert mmrv((exp(1/x-exp(-x))-exp(1/x))*exp(x), x) == set([exp(x), exp(-x)])
    assert mmrv(1/exp(-x+exp(-x))-exp(x), x) == set([exp(x), exp(-x), exp(x-exp(-x))])
    assert mmrv(log(log(x*exp(x*exp(x))+1)), x) == set([exp(x*exp(x))])
    assert mmrv(exp(exp(log(log(x)+1/x))), x) == set([x])

def test_mrv4():
    ln = log
    assert mmrv((ln(ln(x)+ln(ln(x)))-ln(ln(x)))/ln(ln(x)+ln(ln(ln(x))))*ln(x),
            x) == set([x])
    assert mmrv(log(log(x*exp(x*exp(x))+1)) - exp(exp(log(log(x)+1/x))), x) == \
        set([exp(x*exp(x))])

def mrewrite(a, b, c): return rewrite(a[1], a[0], b, c)
def test_rewrite1():
    e = exp(x)
    assert mrewrite(mrv(e, x), x, m) == (1/m, -x)
    e = exp(x**2)
    assert mrewrite(mrv(e, x), x, m) == (1/m, -x**2)
    e = exp(x+1/x)
    assert mrewrite(mrv(e, x), x, m) == (1/m, -x-1/x)
    e = 1/exp(-x+exp(-x))-exp(x)
    assert mrewrite(mrv(e, x), x, m) == (1/(m*exp(m))-1/m, -x)

def test_rewrite2():
    e = exp(x)*log(log(exp(x)))
    assert mmrv(e, x) == set([exp(x)])
    assert mrewrite(mrv(e, x), x, m) == (1/m*log(x), -x)

#sometimes infinite recursion due to log(exp(x**2)) not simplifying
def test_rewrite3():
    e = exp(-x+1/x**2)-exp(x+1/x)
    #both of these are correct and should be equivalent:
    assert mrewrite(mrv(e, x), x, m) in [(-1/m + m*exp(1/x+1/x**2), -x-1/x), (m - 1/m*exp(1/x + x**(-2)), x**(-2) - x)]

def test_mrv_leadterm1():
    assert mrv_leadterm(-exp(1/x), x) == (-1, 0)
    assert mrv_leadterm(1/exp(-x+exp(-x))-exp(x), x) == (-1, 0)
    assert mrv_leadterm((exp(1/x-exp(-x))-exp(1/x))*exp(x), x) == (-exp(1/x), 0)

def test_mrv_leadterm2():
    #Gruntz: p51, 3.25
    assert mrv_leadterm((log(exp(x)+x)-x)/log(exp(x)+log(x))*exp(x), x) == \
            (1, 0)

def test_mrv_leadterm3():
    #Gruntz: p56, 3.27
    assert mmrv(exp(-x+exp(-x)*exp(-x*log(x))), x) == set([exp(-x-x*log(x))])
    assert mrv_leadterm(exp(-x+exp(-x)*exp(-x*log(x))), x) == (exp(-x), 0)

def test_limit1():
    assert gruntz(x, x, oo) == oo
    assert gruntz(x, x, -oo) == -oo
    assert gruntz(-x, x, oo) == -oo
    assert gruntz(x**2, x, -oo) == oo
    assert gruntz(-x**2, x, oo) == -oo
    assert gruntz(x*log(x), x, 0, dir="+") == 0
    assert gruntz(1/x,x,oo) == 0
    assert gruntz(exp(x),x,oo) == oo
    assert gruntz(-exp(x),x,oo) == -oo
    assert gruntz(exp(x)/x,x,oo) == oo
    assert gruntz(1/x-exp(-x),x,oo) == 0
    assert gruntz(x+1/x,x,oo) == oo


def test_limit2():
    assert gruntz(x**x, x, 0, dir="+") == 1
    assert gruntz((exp(x)-1)/x, x, 0) == 1
    assert gruntz(1+1/x,x,oo) == 1
    assert gruntz(-exp(1/x),x,oo) == -1
    assert gruntz(x+exp(-x),x,oo) == oo
    assert gruntz(x+exp(-x**2),x,oo) == oo
    assert gruntz(x+exp(-exp(x)),x,oo) == oo
    assert gruntz(13+1/x-exp(-x),x,oo) == 13

def test_limit3():
    a = Symbol('a')
    assert gruntz(x-log(1+exp(x)), x, oo) == 0
    assert gruntz(x-log(a+exp(x)), x, oo) == 0
    assert gruntz(exp(x)/(1+exp(x)), x, oo) == 1
    assert gruntz(exp(x)/(a+exp(x)), x, oo) == 1

@XFAIL
def test_limit4():
    skip("Takes too long")
    #issue 364
    assert gruntz((3**x+5**x)**(1/x), x, oo) == 5
    #issue 364
    assert gruntz((3**(1/x)+5**(1/x))**x, x, 0) == 5

#@XFAIL
#def test_MrvTestCase_page47_ex3_21():
#    h = exp(-x/(1+exp(-x)))
#    expr = exp(h)*exp(-x/(1+h))*exp(exp(-x+h))/h**2-exp(x)+x
#    expected = set([1/h,exp(x),exp(x-h),exp(x/(1+h))])
#    # XXX Incorrect result
#    assert mrv(expr,x).difference(expected) == set()

def test_I():
    y = Symbol("y")
    assert gruntz(I*x, x, oo) == I*oo
    assert gruntz(y*I*x, x, oo) == y*I*oo
    assert gruntz(y*3*I*x, x, oo) == y*I*oo
    assert gruntz(y*3*sin(I)*x, x, oo) == y*I*oo

def test_issue1715():
    assert gruntz((x + 1)**(1/log(x + 1)), x, oo) == E

def test_intractable():
    assert gruntz(1/gamma(x), x, oo) == 0
    assert gruntz(1/loggamma(x), x, oo) == 0
    assert gruntz(gamma(x)/loggamma(x), x, oo) == oo
    assert gruntz(exp(gamma(x))/gamma(x), x, oo) == oo
    assert gruntz(gamma(x), x, 3) == 2
    assert gruntz(gamma(S(1)/7+1/x), x, oo) == gamma(S(1)/7)
    assert gruntz(log(x**x)/log(gamma(x)), x, oo) == 1
    assert gruntz(log(gamma(gamma(x)))/exp(x), x, oo) == oo

def test_aseries_trig():
    assert cancel(gruntz(1/log(atan(x)), x, oo) \
           - 1/(log(pi) + log(S(1)/2))) == 0
    assert gruntz(1/acot(x), x, -oo) == -oo

def test_exp_log_series():
    assert gruntz(x/log(log(x*exp(x))), x, oo) == oo

def test_issue545():
    assert gruntz(((x**7+x+1)/(2**x+x**2))**(-1/x), x, oo) == 2

def test_slow():
    skip("This is very slow.")
    assert gruntz(loggamma(loggamma(x)), x, oo) == oo
