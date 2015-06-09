from sympy import (
    Abs, acos, acosh, Add, asin, asinh, atan, Ci, cos, sinh, cosh, tanh,
    Derivative, diff, DiracDelta, E, exp, erf, erfi, EulerGamma, factor, Function,
    I, Integral, integrate, Interval, Lambda, LambertW, log,
    Matrix, O, oo, pi, Piecewise, Poly, Rational, S, simplify, sin, tan, sqrt,
    sstr, Sum, Symbol, symbols, sympify, trigsimp,
    Tuple, nan, And, Eq, Ne, re, im, polar_lift, meijerg
)
from sympy.functions.elementary.complexes import periodic_argument
from sympy.integrals.risch import NonElementaryIntegral
from sympy.physics import units
from sympy.core.compatibility import range
from sympy.utilities.pytest import XFAIL, raises, slow
from sympy.utilities.randtest import verify_numerically


x, y, a, t, x_1, x_2, z, s = symbols('x y a t x_1 x_2 z s')
n = Symbol('n', integer=True)
f = Function('f')


def diff_test(i):
    """Return the set of symbols, s, which were used in testing that
    i.diff(s) agrees with i.doit().diff(s). If there is an error then
    the assertion will fail, causing the test to fail."""
    syms = i.free_symbols
    for s in syms:
        assert (i.diff(s).doit() - i.doit().diff(s)).expand() == 0
    return syms


def test_improper_integral():
    assert integrate(log(x), (x, 0, 1)) == -1
    assert integrate(x**(-2), (x, 1, oo)) == 1


def test_constructor():
    # this is shared by Sum, so testing Integral's constructor
    # is equivalent to testing Sum's
    s1 = Integral(n, n)
    assert s1.limits == (Tuple(n),)
    s2 = Integral(n, (n,))
    assert s2.limits == (Tuple(n),)
    s3 = Integral(Sum(x, (x, 1, y)))
    assert s3.limits == (Tuple(y),)
    s4 = Integral(n, Tuple(n,))
    assert s4.limits == (Tuple(n),)

    s5 = Integral(n, (n, Interval(1, 2)))
    assert s5.limits == (Tuple(n, 1, 2),)


def test_basics():

    assert Integral(0, x) != 0
    assert Integral(x, (x, 1, 1)) != 0
    assert Integral(oo, x) != oo
    assert Integral(S.NaN, x) == S.NaN

    assert diff(Integral(y, y), x) == 0
    assert diff(Integral(x, (x, 0, 1)), x) == 0
    assert diff(Integral(x, x), x) == x
    assert diff(Integral(t, (t, 0, x)), x) == x + Integral(0, (t, 0, x))

    e = (t + 1)**2
    assert diff(integrate(e, (t, 0, x)), x) == \
        diff(Integral(e, (t, 0, x)), x).doit().expand() == \
        ((1 + x)**2).expand()
    assert diff(integrate(e, (t, 0, x)), t) == \
        diff(Integral(e, (t, 0, x)), t) == 0
    assert diff(integrate(e, (t, 0, x)), a) == \
        diff(Integral(e, (t, 0, x)), a) == 0
    assert diff(integrate(e, t), a) == diff(Integral(e, t), a) == 0

    assert integrate(e, (t, a, x)).diff(x) == \
        Integral(e, (t, a, x)).diff(x).doit().expand()
    assert Integral(e, (t, a, x)).diff(x).doit() == ((1 + x)**2)
    assert integrate(e, (t, x, a)).diff(x).doit() == (-(1 + x)**2).expand()

    assert integrate(t**2, (t, x, 2*x)).diff(x) == 7*x**2

    assert Integral(x, x).atoms() == set([x])
    assert Integral(f(x), (x, 0, 1)).atoms() == set([S(0), S(1), x])

    assert diff_test(Integral(x, (x, 3*y))) == set([y])
    assert diff_test(Integral(x, (a, 3*y))) == set([x, y])

    assert integrate(x, (x, oo, oo)) == 0 #issue 8171
    assert integrate(x, (x, -oo, -oo)) == 0

    # sum integral of terms
    assert integrate(y + x + exp(x), x) == x*y + x**2/2 + exp(x)

    assert Integral(x).is_commutative
    n = Symbol('n', commutative=False)
    assert Integral(n + x, x).is_commutative is False


def test_basics_multiple():

    assert diff_test(Integral(x, (x, 3*x, 5*y), (y, x, 2*x))) == set([x])
    assert diff_test(Integral(x, (x, 5*y), (y, x, 2*x))) == set([x])
    assert diff_test(Integral(x, (x, 5*y), (y, y, 2*x))) == set([x, y])
    assert diff_test(Integral(y, y, x)) == set([x, y])
    assert diff_test(Integral(y*x, x, y)) == set([x, y])
    assert diff_test(Integral(x + y, y, (y, 1, x))) == set([x])
    assert diff_test(Integral(x + y, (x, x, y), (y, y, x))) == set([x, y])


def test_conjugate_transpose():
    A, B = symbols("A B", commutative=False)

    x = Symbol("x", complex=True)
    p = Integral(A*B, (x,))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()

    x = Symbol("x", real=True)
    p = Integral(A*B, (x,))
    assert p.adjoint().doit() == p.doit().adjoint()
    assert p.conjugate().doit() == p.doit().conjugate()
    assert p.transpose().doit() == p.doit().transpose()


def test_integration():
    assert integrate(0, (t, 0, x)) == 0
    assert integrate(3, (t, 0, x)) == 3*x
    assert integrate(t, (t, 0, x)) == x**2/2
    assert integrate(3*t, (t, 0, x)) == 3*x**2/2
    assert integrate(3*t**2, (t, 0, x)) == x**3
    assert integrate(1/t, (t, 1, x)) == log(x)
    assert integrate(-1/t**2, (t, 1, x)) == 1/x - 1
    assert integrate(t**2 + 5*t - 8, (t, 0, x)) == x**3/3 + 5*x**2/2 - 8*x
    assert integrate(x**2, x) == x**3/3
    assert integrate((3*t*x)**5, x) == (3*t)**5 * x**6 / 6

    b = Symbol("b")
    c = Symbol("c")
    assert integrate(a*t, (t, 0, x)) == a*x**2/2
    assert integrate(a*t**4, (t, 0, x)) == a*x**5/5
    assert integrate(a*t**2 + b*t + c, (t, 0, x)) == a*x**3/3 + b*x**2/2 + c*x


def test_multiple_integration():
    assert integrate((x**2)*(y**2), (x, 0, 1), (y, -1, 2)) == Rational(1)
    assert integrate((y**2)*(x**2), x, y) == Rational(1, 9)*(x**3)*(y**3)
    assert integrate(1/(x + 3)/(1 + x)**3, x) == \
        -S(1)/8*log(3 + x) + S(1)/8*log(1 + x) + x/(4 + 8*x + 4*x**2)


def test_issue_3532():
    assert integrate(exp(-x), (x, 0, oo)) == 1


def test_issue_3560():
    assert integrate(sqrt(x)**3, x) == 2*sqrt(x)**5/5
    assert integrate(sqrt(x), x) == 2*sqrt(x)**3/3
    assert integrate(1/sqrt(x)**3, x) == -2/sqrt(x)


def test_integrate_poly():
    p = Poly(x + x**2*y + y**3, x, y)

    qx = integrate(p, x)
    qy = integrate(p, y)

    assert isinstance(qx, Poly) is True
    assert isinstance(qy, Poly) is True

    assert qx.gens == (x, y)
    assert qy.gens == (x, y)

    assert qx.as_expr() == x**2/2 + x**3*y/3 + x*y**3
    assert qy.as_expr() == x*y + x**2*y**2/2 + y**4/4


def test_integrate_poly_defined():
    p = Poly(x + x**2*y + y**3, x, y)

    Qx = integrate(p, (x, 0, 1))
    Qy = integrate(p, (y, 0, pi))

    assert isinstance(Qx, Poly) is True
    assert isinstance(Qy, Poly) is True

    assert Qx.gens == (y,)
    assert Qy.gens == (x,)

    assert Qx.as_expr() == Rational(1, 2) + y/3 + y**3
    assert Qy.as_expr() == pi**4/4 + pi*x + pi**2*x**2/2


def test_integrate_omit_var():
    y = Symbol('y')

    assert integrate(x) == x**2/2

    raises(ValueError, lambda: integrate(2))
    raises(ValueError, lambda: integrate(x*y))


def test_integrate_poly_accurately():
    y = Symbol('y')
    assert integrate(x*sin(y), x) == x**2*sin(y)/2

    # when passed to risch_norman, this will be a CPU hog, so this really
    # checks, that integrated function is recognized as polynomial
    assert integrate(x**1000*sin(y), x) == x**1001*sin(y)/1001


def test_issue_3635():
    y = Symbol('y')
    assert integrate(x**2, y) == x**2*y
    assert integrate(x**2, (y, -1, 1)) == 2*x**2

# works in sympy and py.test but hangs in `setup.py test`


def test_integrate_linearterm_pow():
    # check integrate((a*x+b)^c, x)  --  issue 3499
    y = Symbol('y', positive=True)
    # TODO: Remove conds='none' below, let the assumption take care of it.
    assert integrate(x**y, x, conds='none') == x**(y + 1)/(y + 1)
    assert integrate((exp(y)*x + 1/y)**(1 + sin(y)), x, conds='none') == \
        exp(-y)*(exp(y)*x + 1/y)**(2 + sin(y)) / (2 + sin(y))


def test_issue_3618():
    assert integrate(pi*sqrt(x), x) == 2*pi*sqrt(x)**3/3
    assert integrate(pi*sqrt(x) + E*sqrt(x)**3, x) == \
        2*pi*sqrt(x)**3/3 + 2*E *sqrt(x)**5/5


def test_issue_3623():
    assert integrate(cos((n + 1)*x), x) == Piecewise(
        (x, Eq(n + 1, 0)), (sin((n + 1)*x)/(n + 1), True))
    assert integrate(cos((n - 1)*x), x) == Piecewise(
        (x, Eq(n - 1, 0)), (sin((n - 1)*x)/(n - 1), True))
    assert integrate(cos((n + 1)*x) + cos((n - 1)*x), x) == \
        Piecewise((x, Eq(n + 1, 0)), (sin((n + 1)*x)/(n + 1), True)) + \
        Piecewise((x, Eq(n - 1, 0)), (sin((n - 1)*x)/(n - 1), True))


def test_issue_3664():
    n = Symbol('n', integer=True, nonzero=True)
    assert integrate(-1./2 * x * sin(n * pi * x/2), [x, -2, 0]) == \
        2*cos(pi*n)/(pi*n)
    assert integrate(-Rational(1)/2 * x * sin(n * pi * x/2), [x, -2, 0]) == \
        2*cos(pi*n)/(pi*n)


def test_issue_3679():
    # definite integration of rational functions gives wrong answers
    assert NS(Integral(1/(x**2 - 8*x + 17), (x, 2, 4))) == '1.10714871779409'


def test_issue_3686():  # remove this when fresnel itegrals are implemented
    from sympy import expand_func, fresnels
    assert expand_func(integrate(sin(x**2), x)) == \
        sqrt(2)*sqrt(pi)*fresnels(sqrt(2)*x/sqrt(pi))/2

def test_integrate_units():
    m = units.m
    s = units.s
    assert integrate(x * m/s, (x, 1*s, 5*s)) == 12*m*s


def test_transcendental_functions():
    assert integrate(LambertW(2*x), x) == \
        -x + x*LambertW(2*x) + x/LambertW(2*x)


def test_issue_3740():
    f = 4*log(x) - 2*log(x)**2
    fid = diff(integrate(f, x), x)
    assert abs(f.subs(x, 42).evalf() - fid.subs(x, 42).evalf()) < 1e-10


def test_issue_3788():
    assert integrate(1/(1 + x**2), x) == atan(x)


def test_issue_3952():
    f = sin(x)
    assert integrate(f, x) == -cos(x)
    raises(ValueError, lambda: integrate(f, 2*x))


def test_issue_4516():
    assert integrate(2**x - 2*x, x) == 2**x/log(2) - x**2


def test_issue_7450():
    ans = integrate(exp(-(1 + I)*x), (x, 0, oo))
    assert re(ans) == S.Half and im(ans) == -S.Half


def test_matrices():
    M = Matrix(2, 2, lambda i, j: (i + j + 1)*sin((i + j + 1)*x))

    assert integrate(M, x) == Matrix([
        [-cos(x), -cos(2*x)],
        [-cos(2*x), -cos(3*x)],
    ])


def test_integrate_functions():
    # issue 4111
    assert integrate(f(x), x) == Integral(f(x), x)
    assert integrate(f(x), (x, 0, 1)) == Integral(f(x), (x, 0, 1))
    assert integrate(f(x)*diff(f(x), x), x) == f(x)**2/2
    assert integrate(diff(f(x), x) / f(x), x) == log(f(x))


def test_integrate_derivatives():
    assert integrate(Derivative(f(x), x), x) == f(x)
    assert integrate(Derivative(f(y), y), x) == x*Derivative(f(y), y)


def test_transform():
    a = Integral(x**2 + 1, (x, -1, 2))
    fx = x
    fy = 3*y + 1
    assert a.doit() == a.transform(fx, fy).doit()
    assert a.transform(fx, fy).transform(fy, fx) == a
    fx = 3*x + 1
    fy = y
    assert a.transform(fx, fy).transform(fy, fx) == a
    a = Integral(sin(1/x), (x, 0, 1))
    assert a.transform(x, 1/y) == Integral(sin(y)/y**2, (y, 1, oo))
    assert a.transform(x, 1/y).transform(y, 1/x) == a
    a = Integral(exp(-x**2), (x, -oo, oo))
    assert a.transform(x, 2*y) == Integral(2*exp(-4*y**2), (y, -oo, oo))
    # < 3 arg limit handled properly
    assert Integral(x, x).transform(x, a*y).doit() == \
        Integral(y*a**2, y).doit()
    _3 = S(3)
    assert Integral(x, (x, 0, -_3)).transform(x, 1/y).doit() == \
        Integral(-1/x**3, (x, -oo, -1/_3)).doit()
    assert Integral(x, (x, 0, _3)).transform(x, 1/y) == \
        Integral(y**(-3), (y, 1/_3, oo))
    # issue 8400
    i = Integral(x + y, (x, 1, 2), (y, 1, 2))
    assert i.transform(x, (x + 2*y, x)).doit() == \
        i.transform(x, (x + 2*z, x)).doit() == 3


def test_issue_4052():
    f = S(1)/2*asin(x) + x*sqrt(1 - x**2)/2

    assert integrate(cos(asin(x)), x) == f
    assert integrate(sin(acos(x)), x) == f


def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)


@slow
def test_evalf_integrals():
    assert NS(Integral(x, (x, 2, 5)), 15) == '10.5000000000000'
    gauss = Integral(exp(-x**2), (x, -oo, oo))
    assert NS(gauss, 15) == '1.77245385090552'
    assert NS(gauss**2 - pi + E*Rational(
        1, 10**20), 15) in ('2.71828182845904e-20', '2.71828182845905e-20')
    # A monster of an integral from http://mathworld.wolfram.com/DefiniteIntegral.html
    t = Symbol('t')
    a = 8*sqrt(3)/(1 + 3*t**2)
    b = 16*sqrt(2)*(3*t + 1)*sqrt(4*t**2 + t + 1)**3
    c = (3*t**2 + 1)*(11*t**2 + 2*t + 3)**2
    d = sqrt(2)*(249*t**2 + 54*t + 65)/(11*t**2 + 2*t + 3)**2
    f = a - b/c - d
    assert NS(Integral(f, (t, 0, 1)), 50) == \
        NS((3*sqrt(2) - 49*pi + 162*atan(sqrt(2)))/12, 50)
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(1/x))/(1 + x + x**2), (x, 0, 1)), 15) == \
        NS('pi/sqrt(3) * log(2*pi**(5/6) / gamma(1/6))', 15)
    # http://mathworld.wolfram.com/AhmedsIntegral.html
    assert NS(Integral(atan(sqrt(x**2 + 2))/(sqrt(x**2 + 2)*(x**2 + 1)), (x,
              0, 1)), 15) == NS(5*pi**2/96, 15)
    # http://mathworld.wolfram.com/AbelsIntegral.html
    assert NS(Integral(x/((exp(pi*x) - exp(
        -pi*x))*(x**2 + 1)), (x, 0, oo)), 15) == NS('log(2)/2-1/4', 15)
    # Complex part trimming
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(sin(x)/cos(x))), (x, pi/4, pi/2)), 15, chop=True) == \
        NS('pi/4*log(4*pi**3/gamma(1/4)**4)', 15)
    #
    # Endpoints causing trouble (rounding error in integration points -> complex log)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 17, chop=True) == NS(2, 17)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 20, chop=True) == NS(2, 20)
    assert NS(
        2 + Integral(log(2*cos(x/2)), (x, -pi, pi)), 22, chop=True) == NS(2, 22)
    # Needs zero handling
    assert NS(pi - 4*Integral(
        'sqrt(1-x**2)', (x, 0, 1)), 15, maxn=30, chop=True) in ('0.0', '0')
    # Oscillatory quadrature
    a = Integral(sin(x)/x**2, (x, 1, oo)).evalf(maxn=15)
    assert 0.49 < a < 0.51
    assert NS(
        Integral(sin(x)/x**2, (x, 1, oo)), quad='osc') == '0.504067061906928'
    assert NS(Integral(
        cos(pi*x + 1)/x, (x, -oo, -1)), quad='osc') == '0.276374705640365'
    # indefinite integrals aren't evaluated
    assert NS(Integral(x, x)) == 'Integral(x, x)'
    assert NS(Integral(x, (x, y))) == 'Integral(x, (x, y))'


def test_evalf_issue_939():
    # https://github.com/sympy/sympy/issues/4038

    # The output form of an integral may differ by a step function between
    # revisions, making this test a bit useless. This can't be said about
    # other two tests. For now, all values of this evaluation are used here,
    # but in future this should be reconsidered.
    assert NS(integrate(1/(x**5 + 1), x).subs(x, 4), chop=True) in \
        ['-0.000976138910649103', '0.965906660135753', '1.93278945918216']

    assert NS(Integral(1/(x**5 + 1), (x, 2, 4))) == '0.0144361088886740'
    assert NS(
        integrate(1/(x**5 + 1), (x, 2, 4)), chop=True) == '0.0144361088886740'


@XFAIL
def test_failing_integrals():
    #---
    # Double integrals not implemented
    assert NS(Integral(
        sqrt(x) + x*y, (x, 1, 2), (y, -1, 1)), 15) == '2.43790283299492'
    # double integral + zero detection
    assert NS(Integral(sin(x + x*y), (x, -1, 1), (y, -1, 1)), 15) == '0.0'


def test_integrate_DiracDelta():
    # This is here to check that deltaintegrate is being called, but also
    # to test definite integrals. More tests are in test_deltafunctions.py
    assert integrate(DiracDelta(x) * f(x), (x, -oo, oo)) == f(0)
    assert integrate(DiracDelta(x) * f(x), (x, 0, oo)) == f(0)/2
    assert integrate(DiracDelta(x)**2, (x, -oo, oo)) == DiracDelta(0)
    # issue 4522
    assert integrate(integrate((4 - 4*x + x*y - 4*y) * \
        DiracDelta(x)*DiracDelta(y - 1), (x, 0, 1)), (y, 0, 1)) == 0
    # issue 5729
    p = exp(-(x**2 + y**2))/pi
    assert integrate(p*DiracDelta(x - 10*y), (x, -oo, oo), (y, -oo, oo)) == \
        integrate(p*DiracDelta(x - 10*y), (y, -oo, oo), (x, -oo, oo)) == \
        integrate(p*DiracDelta(10*x - y), (x, -oo, oo), (y, -oo, oo)) == \
        integrate(p*DiracDelta(10*x - y), (y, -oo, oo), (x, -oo, oo)) == \
        1/sqrt(101*pi)


@XFAIL
def test_integrate_DiracDelta_fails():
    # issue 6427
    assert integrate(integrate(integrate(
        DiracDelta(x - y - z), (z, 0, oo)), (y, 0, 1)), (x, 0, 1)) == S(1)/2


def test_integrate_returns_piecewise():
    assert integrate(x**y, x) == Piecewise(
        (log(x), Eq(y, -1)), (x**(y + 1)/(y + 1), True))
    assert integrate(x**y, y) == Piecewise(
        (y, Eq(log(x), 0)), (x**y/log(x), True))
    assert integrate(exp(n*x), x) == Piecewise(
        (x, Eq(n, 0)), (exp(n*x)/n, True))
    assert integrate(x*exp(n*x), x) == Piecewise(
        (x**2/2, Eq(n**3, 0)), ((x*n**2 - n)*exp(n*x)/n**3, True))
    assert integrate(x**(n*y), x) == Piecewise(
        (log(x), Eq(n*y, -1)), (x**(n*y + 1)/(n*y + 1), True))
    assert integrate(x**(n*y), y) == Piecewise(
        (y, Eq(n*log(x), 0)), (x**(n*y)/(n*log(x)), True))
    assert integrate(cos(n*x), x) == Piecewise(
        (x, Eq(n, 0)), (sin(n*x)/n, True))
    assert integrate(cos(n*x)**2, x) == Piecewise(
        (x, Eq(n, 0)), ((n*x/2 + sin(n*x)*cos(n*x)/2)/n, True))
    assert integrate(x*cos(n*x), x) == Piecewise(
        (x**2/2, Eq(n, 0)), (x*sin(n*x)/n + cos(n*x)/n**2, True))
    assert integrate(sin(n*x), x) == Piecewise(
        (0, Eq(n, 0)), (-cos(n*x)/n, True))
    assert integrate(sin(n*x)**2, x) == Piecewise(
        (0, Eq(n, 0)), ((n*x/2 - sin(n*x)*cos(n*x)/2)/n, True))
    assert integrate(x*sin(n*x), x) == Piecewise(
        (0, Eq(n, 0)), (-x*cos(n*x)/n + sin(n*x)/n**2, True))
    assert integrate(exp(x*y),(x,0,z)) == Piecewise( \
        (z, Eq(y,0)), (exp(y*z)/y - 1/y, True))


def test_subs1():
    e = Integral(exp(x - y), x)
    assert e.subs(y, 3) == Integral(exp(x - 3), x)
    e = Integral(exp(x - y), (x, 0, 1))
    assert e.subs(y, 3) == Integral(exp(x - 3), (x, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo))


def test_subs2():
    e = Integral(exp(x - y), x, t)
    assert e.subs(y, 3) == Integral(exp(x - 3), x, t)
    e = Integral(exp(x - y), (x, 0, 1), (t, 0, 1))
    assert e.subs(y, 3) == Integral(exp(x - 3), (x, 0, 1), (t, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo), (t, 0, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs3():
    e = Integral(exp(x - y), (x, 0, y), (t, y, 1))
    assert e.subs(y, 3) == Integral(exp(x - 3), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x - y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs4():
    e = Integral(exp(x), (x, 0, y), (t, y, 1))
    assert e.subs(y, 3) == Integral(exp(x), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x: 0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))


def test_subs5():
    e = Integral(exp(-x**2), (x, -oo, oo))
    assert e.subs(x, 5) == e
    e = Integral(exp(-x**2 + y), x)
    assert e.subs(y, 5) == Integral(exp(-x**2 + 5), x)
    e = Integral(exp(-x**2 + y), (x, x))
    assert e.subs(x, 5) == Integral(exp(y - x**2), (x, 5))
    assert e.subs(y, 5) == Integral(exp(-x**2 + 5), x)
    e = Integral(exp(-x**2 + y), (y, -oo, oo), (x, -oo, oo))
    assert e.subs(x, 5) == e
    assert e.subs(y, 5) == e
    # Test evaluation of antiderivatives
    e = Integral(exp(-x**2), (x, x))
    assert e.subs(x, 5) == Integral(exp(-x**2), (x, 5))
    e = Integral(exp(x), x)
    assert (e.subs(x,1) - e.subs(x,0) - Integral(exp(x), (x, 0, 1))
        ).doit().is_zero


def test_subs6():
    a, b = symbols('a b')
    e = Integral(x*y, (x, f(x), f(y)))
    assert e.subs(x, 1) == Integral(x*y, (x, f(1), f(y)))
    assert e.subs(y, 1) == Integral(x, (x, f(x), f(1)))
    e = Integral(x*y, (x, f(x), f(y)), (y, f(x), f(y)))
    assert e.subs(x, 1) == Integral(x*y, (x, f(1), f(y)), (y, f(1), f(y)))
    assert e.subs(y, 1) == Integral(x*y, (x, f(x), f(y)), (y, f(x), f(1)))
    e = Integral(x*y, (x, f(x), f(a)), (y, f(x), f(a)))
    assert e.subs(a, 1) == Integral(x*y, (x, f(x), f(1)), (y, f(x), f(1)))


def test_subs7():
    e = Integral(x, (x, 1, y), (y, 1, 2))
    assert e.subs({x: 1, y: 2}) == e
    e = Integral(sin(x) + sin(y), (x, sin(x), sin(y)),
                                  (y, 1, 2))
    assert e.subs(sin(y), 1) == e
    assert e.subs(sin(x), 1) == Integral(sin(x) + sin(y), (x, 1, sin(y)),
                                         (y, 1, 2))

def test_expand():
    e = Integral(f(x)+f(x**2), (x, 1, y))
    assert e.expand() == Integral(f(x), (x, 1, y)) + Integral(f(x**2), (x, 1, y))

def test_integration_variable():
    raises(ValueError, lambda: Integral(exp(-x**2), 3))
    raises(ValueError, lambda: Integral(exp(-x**2), (3, -oo, oo)))


def test_expand_integral():
    assert Integral(cos(x**2)*(sin(x**2) + 1), (x, 0, 1)).expand() == \
        Integral(cos(x**2)*sin(x**2), (x, 0, 1)) + \
        Integral(cos(x**2), (x, 0, 1))
    assert Integral(cos(x**2)*(sin(x**2) + 1), x).expand() == \
        Integral(cos(x**2)*sin(x**2), x) + \
        Integral(cos(x**2), x)


def test_as_sum_midpoint1():
    e = Integral(sqrt(x**3 + 1), (x, 2, 10))
    assert e.as_sum(1, method="midpoint") == 8*sqrt(217)
    assert e.as_sum(2, method="midpoint") == 4*sqrt(65) + 12*sqrt(57)
    assert e.as_sum(3, method="midpoint") == 8*sqrt(217)/3 + \
        8*sqrt(3081)/27 + 8*sqrt(52809)/27
    assert e.as_sum(4, method="midpoint") == 2*sqrt(730) + \
        4*sqrt(7) + 4*sqrt(86) + 6*sqrt(14)
    assert abs(e.as_sum(4, method="midpoint").n() - e.n()) < 0.5

    e = Integral(sqrt(x**3 + y**3), (x, 2, 10), (y, 0, 10))
    raises(NotImplementedError, lambda: e.as_sum(4))


def test_as_sum_midpoint2():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method="midpoint").expand() == S(1)/4 + y + y**2
    assert e.as_sum(2, method="midpoint").expand() == S(5)/16 + y + y**2
    assert e.as_sum(3, method="midpoint").expand() == S(35)/108 + y + y**2
    assert e.as_sum(4, method="midpoint").expand() == S(21)/64 + y + y**2


def test_as_sum_left():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method="left").expand() == y**2
    assert e.as_sum(2, method="left").expand() == S(1)/8 + y/2 + y**2
    assert e.as_sum(3, method="left").expand() == S(5)/27 + 2*y/3 + y**2
    assert e.as_sum(4, method="left").expand() == S(7)/32 + 3*y/4 + y**2


def test_as_sum_right():
    e = Integral((x + y)**2, (x, 0, 1))
    assert e.as_sum(1, method="right").expand() == 1 + 2*y + y**2
    assert e.as_sum(2, method="right").expand() == S(5)/8 + 3*y/2 + y**2
    assert e.as_sum(3, method="right").expand() == S(14)/27 + 4*y/3 + y**2
    assert e.as_sum(4, method="right").expand() == S(15)/32 + 5*y/4 + y**2


def test_as_sum_raises():
    e = Integral((x + y)**2, (x, 0, 1))
    raises(ValueError, lambda: e.as_sum(-1))
    raises(ValueError, lambda: e.as_sum(0))
    raises(ValueError, lambda: Integral(x).as_sum(3))
    raises(NotImplementedError, lambda: e.as_sum(oo))
    raises(NotImplementedError, lambda: e.as_sum(3, method='xxxx2'))


def test_nested_doit():
    e = Integral(Integral(x, x), x)
    f = Integral(x, x, x)
    assert e.doit() == f.doit()


def test_issue_4665():
    # Allow only upper or lower limit evaluation
    e = Integral(x**2, (x, None, 1))
    f = Integral(x**2, (x, 1, None))
    assert e.doit() == Rational(1, 3)
    assert f.doit() == Rational(-1, 3)
    assert Integral(x*y, (x, None, y)).subs(y, t) == Integral(x*t, (x, None, t))
    assert Integral(x*y, (x, y, None)).subs(y, t) == Integral(x*t, (x, t, None))
    assert integrate(x**2, (x, None, 1)) == Rational(1, 3)
    assert integrate(x**2, (x, 1, None)) == Rational(-1, 3)
    assert integrate("x**2", ("x", "1", None)) == Rational(-1, 3)


def test_integral_reconstruct():
    e = Integral(x**2, (x, -1, 1))
    assert e == Integral(*e.args)


def test_doit_integrals():
    e = Integral(Integral(2*x), (x, 0, 1))
    assert e.doit() == Rational(1, 3)
    assert e.doit(deep=False) == Rational(1, 3)
    f = Function('f')
    # doesn't matter if the integral can't be performed
    assert Integral(f(x), (x, 1, 1)).doit() == 0
    # doesn't matter if the limits can't be evaluated
    assert Integral(0, (x, 1, Integral(f(x), x))).doit() == 0
    assert Integral(x, (a, 0)).doit() == 0
    limits = ((a, 1, exp(x)), (x, 0))
    assert Integral(a, *limits).doit() == S(1)/4
    assert Integral(a, *list(reversed(limits))).doit() == 0


def test_issue_4884():
    assert integrate(sqrt(x)*(1 + x)) == \
        Piecewise(
            (2*sqrt(x)*(x + 1)**2/5 - 2*sqrt(x)*(x + 1)/15 - 4*sqrt(x)/15,
            Abs(x + 1) > 1),
            (2*I*sqrt(-x)*(x + 1)**2/5 - 2*I*sqrt(-x)*(x + 1)/15 -
            4*I*sqrt(-x)/15, True))
    assert integrate(x**x*(1 + log(x))) == x**x


def test_is_number():
    from sympy.abc import x, y, z
    from sympy import cos, sin
    assert Integral(x).is_number is False
    assert Integral(1, x).is_number is False
    assert Integral(1, (x, 1)).is_number is True
    assert Integral(1, (x, 1, 2)).is_number is True
    assert Integral(1, (x, 1, y)).is_number is False
    assert Integral(1, (x, y)).is_number is False
    assert Integral(x, y).is_number is False
    assert Integral(x, (y, 1, x)).is_number is False
    assert Integral(x, (y, 1, 2)).is_number is False
    assert Integral(x, (x, 1, 2)).is_number is True
    # `foo.is_number` should always be eqivalent to `not foo.free_symbols`
    # in each of these cases, there are pseudo-free symbols
    i = Integral(x, (y, 1, 1))
    assert i.is_number is False and i.n() == 0
    i = Integral(x, (y, z, z))
    assert i.is_number is False and i.n() == 0
    i = Integral(1, (y, z, z + 2))
    assert i.is_number is False and i.n() == 2

    assert Integral(x*y, (x, 1, 2), (y, 1, 3)).is_number is True
    assert Integral(x*y, (x, 1, 2), (y, 1, z)).is_number is False
    assert Integral(x, (x, 1)).is_number is True
    assert Integral(x, (x, 1, Integral(y, (y, 1, 2)))).is_number is True
    assert Integral(Sum(z, (z, 1, 2)), (x, 1, 2)).is_number is True
    # it is possible to get a false negative if the integrand is
    # actually an unsimplified zero, but this is true of is_number in general.
    assert Integral(sin(x)**2 + cos(x)**2 - 1, x).is_number is False
    assert Integral(f(x), (x, 0, 1)).is_number is True


def test_symbols():
    from sympy.abc import x, y, z
    assert Integral(0, x).free_symbols == set([x])
    assert Integral(x).free_symbols == set([x])
    assert Integral(x, (x, None, y)).free_symbols == set([y])
    assert Integral(x, (x, y, None)).free_symbols == set([y])
    assert Integral(x, (x, 1, y)).free_symbols == set([y])
    assert Integral(x, (x, y, 1)).free_symbols == set([y])
    assert Integral(x, (x, x, y)).free_symbols == set([x, y])
    assert Integral(x, x, y).free_symbols == set([x, y])
    assert Integral(x, (x, 1, 2)).free_symbols == set()
    assert Integral(x, (y, 1, 2)).free_symbols == set([x])
    # pseudo-free in this case
    assert Integral(x, (y, z, z)).free_symbols == set([x, z])
    assert Integral(x, (y, 1, 2), (y, None, None)).free_symbols == set([x, y])
    assert Integral(x, (y, 1, 2), (x, 1, y)).free_symbols == set([y])
    assert Integral(2, (y, 1, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (y, x, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (x, 1, 2), (y, x, 2), (y, 1, 2)).free_symbols == \
        set([x])


def test_is_zero():
    from sympy.abc import x, m
    assert Integral(0, (x, 1, x)).is_zero
    assert Integral(1, (x, 1, 1)).is_zero
    assert Integral(1, (x, 1, 2), (y, 2)).is_zero is False
    assert Integral(x, (m, 0)).is_zero
    assert Integral(x + m, (m, 0)).is_zero is None
    i = Integral(m, (m, 1, exp(x)), (x, 0))
    assert i.is_zero is None
    assert Integral(m, (x, 0), (m, 1, exp(x))).is_zero is True

    assert Integral(x, (x, oo, oo)).is_zero # issue 8171
    assert Integral(x, (x, -oo, -oo)).is_zero

    # this is zero but is beyond the scope of what is_zero
    # should be doing
    assert Integral(sin(x), (x, 0, 2*pi)).is_zero is None


def test_series():
    from sympy.abc import x
    i = Integral(cos(x), (x, x))
    e = i.lseries(x)
    assert i.nseries(x, n=8).removeO() == Add(*[next(e) for j in range(4)])


def test_issue_4403():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z', positive=True)
    assert integrate(sqrt(x**2 + z**2), x) == \
        z**2*asinh(x/z)/2 + x*sqrt(x**2 + z**2)/2
    assert integrate(sqrt(x**2 - z**2), x) == \
        -z**2*acosh(x/z)/2 + x*sqrt(x**2 - z**2)/2

    x = Symbol('x', real=True)
    y = Symbol('y', nonzero=True, real=True)
    assert integrate(1/(x**2 + y**2)**S('3/2'), x) == \
        1/(y**2*sqrt(1 + y**2/x**2))


def test_issue_4403_2():
    assert integrate(sqrt(-x**2 - 4), x) == \
        -2*atan(x/sqrt(-4 - x**2)) + x*sqrt(-4 - x**2)/2


def test_issue_4100():
    R = Symbol('R', positive=True)
    assert integrate(sqrt(R**2 - x**2), (x, 0, R)) == pi*R**2/4


def test_issue_5167():
    from sympy.abc import w, x, y, z
    f = Function('f')
    assert Integral(Integral(f(x), x), x) == Integral(f(x), x, x)
    assert Integral(f(x)).args == (f(x), Tuple(x))
    assert Integral(Integral(f(x))).args == (f(x), Tuple(x), Tuple(x))
    assert Integral(Integral(f(x)), y).args == (f(x), Tuple(x), Tuple(y))
    assert Integral(Integral(f(x), z), y).args == (f(x), Tuple(z), Tuple(y))
    assert Integral(Integral(Integral(f(x), x), y), z).args == \
        (f(x), Tuple(x), Tuple(y), Tuple(z))
    assert integrate(Integral(f(x), x), x) == Integral(f(x), x, x)
    assert integrate(Integral(f(x), y), x) == y*Integral(f(x), x)
    assert integrate(Integral(f(x), x), y) in [Integral(y*f(x), x), y*Integral(f(x), x)]
    assert integrate(Integral(2, x), x) == x**2
    assert integrate(Integral(2, x), y) == 2*x*y
    # don't re-order given limits
    assert Integral(1, x, y).args != Integral(1, y, x).args
    # do as many as possibble
    assert Integral(f(x), y, x, y, x).doit() == y**2*Integral(f(x), x, x)/2
    assert Integral(f(x), (x, 1, 2), (w, 1, x), (z, 1, y)).doit() == \
        y*(x - 1)*Integral(f(x), (x, 1, 2)) - (x - 1)*Integral(f(x), (x, 1, 2))


def test_issue_4890():
    z = Symbol('z', positive=True)
    assert integrate(exp(-log(x)**2), x) == \
        sqrt(pi)*exp(S(1)/4)*erf(log(x)-S(1)/2)/2
    assert integrate(exp(log(x)**2), x) == \
        sqrt(pi)*exp(-S(1)/4)*erfi(log(x)+S(1)/2)/2
    assert integrate(exp(-z*log(x)**2), x) == \
        sqrt(pi)*exp(1/(4*z))*erf(sqrt(z)*log(x) - 1/(2*sqrt(z)))/(2*sqrt(z))


def test_issue_4376():
    n = Symbol('n', integer=True, positive=True)
    assert simplify(integrate(n*(x**(1/n) - 1), (x, 0, S.Half)) -
                (n**2 - 2**(1/n)*n**2 - n*2**(1/n))/(2**(1 + 1/n) + n*2**(1 + 1/n))) == 0


def test_issue_4517():
    assert integrate((sqrt(x) - x**3)/x**Rational(1, 3), x) == \
        6*x**Rational(7, 6)/7 - 3*x**Rational(11, 3)/11


def test_issue_4527():
    k, m = symbols('k m', integer=True)
    assert integrate(sin(k*x)*sin(m*x), (x, 0, pi)) == Piecewise(
        (0, And(Eq(k, 0), Eq(m, 0))),
        (-pi/2, Eq(k, -m)),
        (pi/2, Eq(k, m)),
        (0, True))
    assert integrate(sin(k*x)*sin(m*x), (x,)) == Piecewise(
        (0, And(Eq(k, 0), Eq(m, 0))),
        (-x*sin(m*x)**2/2 - x*cos(m*x)**2/2 + sin(m*x)*cos(m*x)/(2*m), Eq(k, -m)),
        (x*sin(m*x)**2/2 + x*cos(m*x)**2/2 - sin(m*x)*cos(m*x)/(2*m), Eq(k, m)),
        (m*sin(k*x)*cos(m*x)/(k**2 - m**2) -
         k*sin(m*x)*cos(k*x)/(k**2 - m**2), True))

def test_issue_4199():
    ypos = Symbol('y', positive=True)
    # TODO: Remove conds='none' below, let the assumption take care of it.
    assert integrate(exp(-I*2*pi*ypos*x)*x, (x, -oo, oo), conds='none') == \
        Integral(exp(-I*2*pi*ypos*x)*x, (x, -oo, oo))


@slow
def test_issue_3940():
    a, b, c, d = symbols('a:d', positive=True, finite=True)
    assert integrate(exp(-x**2 + I*c*x), x) == \
        -sqrt(pi)*exp(-c**2/4)*erf(I*c/2 - x)/2
    assert integrate(exp(a*x**2 + b*x + c), x) == \
        sqrt(pi)*exp(c)*exp(-b**2/(4*a))*erfi(sqrt(a)*x + b/(2*sqrt(a)))/(2*sqrt(a))

    from sympy import expand_mul
    from sympy.abc import k
    assert expand_mul(integrate(exp(-x**2)*exp(I*k*x), (x, -oo, oo))) == \
        sqrt(pi)*exp(-k**2/4)
    a, d = symbols('a d', positive=True)
    assert expand_mul(integrate(exp(-a*x**2 + 2*d*x), (x, -oo, oo))) == \
        sqrt(pi)*exp(d**2/a)/sqrt(a)


def test_issue_5413():
    # Note that this is not the same as testing ratint() becuase integrate()
    # pulls out the coefficient.
    assert integrate(-a/(a**2 + x**2), x) == I*log(-I*a + x)/2 - I*log(I*a + x)/2


def test_issue_4892a():
    A, z = symbols('A z')
    c = Symbol('c', nonzero=True)
    P1 = -A*exp(-z)
    P2 = -A/(c*t)*(sin(x)**2 + cos(y)**2)

    h1 = -sin(x)**2 - cos(y)**2
    h2 = -sin(x)**2 + sin(y)**2 - 1

    # there is still some non-deterministic behavior in integrate
    # or trigsimp which permits one of the following
    assert integrate(c*(P2 - P1), t) in [
        c*(-A*(-h1)*log(c*t)/c + A*t*exp(-z)),
        c*(-A*(-h2)*log(c*t)/c + A*t*exp(-z)),
        c*( A* h1 *log(c*t)/c + A*t*exp(-z)),
        c*( A* h2 *log(c*t)/c + A*t*exp(-z)),
        (A*c*t - A*(-h1)*log(t)*exp(z))*exp(-z),
        (A*c*t - A*(-h2)*log(t)*exp(z))*exp(-z),
    ]


def test_issue_4892b():
    # Issues relating to issue 4596 are making the actual result of this hard
    # to test.  The answer should be something like
    #
    # (-sin(y) + sqrt(-72 + 48*cos(y) - 8*cos(y)**2)/2)*log(x + sqrt(-72 +
    # 48*cos(y) - 8*cos(y)**2)/(2*(3 - cos(y)))) + (-sin(y) - sqrt(-72 +
    # 48*cos(y) - 8*cos(y)**2)/2)*log(x - sqrt(-72 + 48*cos(y) -
    # 8*cos(y)**2)/(2*(3 - cos(y)))) + x**2*sin(y)/2 + 2*x*cos(y)

    expr = (sin(y)*x**3 + 2*cos(y)*x**2 + 12)/(x**2 + 2)
    assert trigsimp(factor(integrate(expr, x).diff(x) - expr)) == 0


def test_issue_5178():
    assert integrate(sin(x)*f(y, z), (x, 0, pi), (y, 0, pi), (z, 0, pi)) == \
        2*Integral(f(y, z), (y, 0, pi), (z, 0, pi))


def test_integrate_series():
    f = sin(x).series(x, 0, 10)
    g = x**2/2 - x**4/24 + x**6/720 - x**8/40320 + x**10/3628800 + O(x**11)

    assert integrate(f, x) == g
    assert diff(integrate(f, x), x) == f

    assert integrate(O(x**5), x) == O(x**6)


def test_atom_bug():
    from sympy import meijerg
    from sympy.integrals.heurisch import heurisch
    assert heurisch(meijerg([], [], [1], [], x), x) is None


def test_limit_bug():
    z = Symbol('z', nonzero=True)
    assert integrate(sin(x*y*z), (x, 0, pi), (y, 0, pi)) == \
        (log(z**2) + 2*EulerGamma + 2*log(pi))/(2*z) - \
        (-log(pi*z) + log(pi**2*z**2)/2 + Ci(pi**2*z))/z + log(pi)/z


def test_issue_4703():
    g = Function('g')
    assert integrate(exp(x)*g(x), x).has(Integral)


def test_issue_1888():
    f = Function('f')
    assert integrate(f(x).diff(x)**2, x).has(Integral)

# The following tests work using meijerint.


def test_issue_3558():
    from sympy import Si
    assert integrate(cos(x*y), (x, -pi/2, pi/2), (y, 0, pi)) == 2*Si(pi**2/2)


def test_issue_4422():
    assert integrate(1/sqrt(16 + 4*x**2), x) == asinh(x/2) / 2


def test_issue_4493():
    from sympy import simplify
    assert simplify(integrate(x*sqrt(1 + 2*x), x)) == \
        sqrt(2*x + 1)*(6*x**2 + x - 1)/15


def test_issue_4737():
    assert integrate(sin(x)/x, (x, -oo, oo)) == pi
    assert integrate(sin(x)/x, (x, 0, oo)) == pi/2


def test_issue_4992():
    # Note: psi in _check_antecedents becomes NaN.
    from sympy import simplify, expand_func, polygamma, gamma
    a = Symbol('a', positive=True)
    assert simplify(expand_func(integrate(exp(-x)*log(x)*x**a, (x, 0, oo)))) == \
        (a*polygamma(0, a) + 1)*gamma(a)


def test_issue_4487():
    from sympy import lowergamma, simplify
    assert simplify(integrate(exp(-x)*x**y, x)) == lowergamma(y + 1, x)


@XFAIL
def test_issue_4215():
    x = Symbol("x")
    assert integrate(1/(x**2), (x, -1, 1)) == oo


def test_issue_4400():
    n = Symbol('n', integer=True, positive=True)
    assert integrate((x**n)*log(x), x) == \
        n*x*x**n*log(x)/(n**2 + 2*n + 1) + x*x**n*log(x)/(n**2 + 2*n + 1) - \
        x*x**n/(n**2 + 2*n + 1)


def test_issue_6253():
    # Note: this used to raise NotImplementedError
    # Note: psi in _check_antecedents becomes NaN.
    assert integrate((sqrt(1 - x) + sqrt(1 + x))**2/x, x, meijerg=True) == \
        Integral((sqrt(-x + 1) + sqrt(x + 1))**2/x, x)


def test_issue_4153():
    assert integrate(1/(1 + x + y + z), (x, 0, 1), (y, 0, 1), (z, 0, 1)) in [
        -12*log(3) - 3*log(6)/2 + 3*log(8)/2 + 5*log(2) + 7*log(4),
        6*log(2) + 8*log(4) - 27*log(3)/2, 22*log(2) - 27*log(3)/2,
        -12*log(3) - 3*log(6)/2 + 47*log(2)/2]


def test_issue_4326():
    R, b, h = symbols('R b h')
    # It doesn't matter if we can do the integral.  Just make sure the result
    # doesn't contain nan.  This is really a test against _eval_interval.
    assert not integrate(((h*(x - R + b))/b)*sqrt(R**2 - x**2), (x, R - b, R)).has(nan)


def test_powers():
    assert integrate(2**x + 3**x, x) == 2**x/log(2) + 3**x/log(3)


def test_risch_option():
    # risch=True only allowed on indefinite integrals
    raises(ValueError, lambda: integrate(1/log(x), (x, 0, oo), risch=True))
    assert integrate(exp(-x**2), x, risch=True) == NonElementaryIntegral(exp(-x**2), x)
    assert integrate(log(1/x)*y, x, y, risch=True) == y**2*(x*log(1/x)/2 + x/2)
    assert integrate(erf(x), x, risch=True) == Integral(erf(x), x)
    # TODO: How to test risch=False?

def test_issue_6828():
    f = 1/(1.08*x**2 - 4.3)
    g = integrate(f, x).diff(x)
    assert verify_numerically(f, g, tol=1e-12)

@XFAIL
def test_integrate_Piecewise_rational_over_reals():
    f = Piecewise(
        (0,                                              t - 478.515625*pi <  0),
        (13.2075145209219*pi/(0.000871222*t + 0.995)**2, t - 478.515625*pi >= 0))

    assert integrate(f, (t, 0, oo)) == 15235.9375*pi


def test_issue_4803():
    x_max = Symbol("x_max")
    assert integrate(y/pi*exp(-(x_max - x)/cos(a)), x) == \
        y*exp((x - x_max)/cos(a))*cos(a)/pi


def test_issue_4234():
    assert integrate(1/sqrt(1 + tan(x)**2)) == tan(x) / sqrt(1 + tan(x)**2)


def test_issue_4492():
    assert simplify(integrate(x**2 * sqrt(5 - x**2), x)) == Piecewise(
        (I*(2*x**5 - 15*x**3 + 25*x - 25*sqrt(x**2 - 5)*acosh(sqrt(5)*x/5)) /
            (8*sqrt(x**2 - 5)), 1 < Abs(x**2)/5),
        ((-2*x**5 + 15*x**3 - 25*x + 25*sqrt(-x**2 + 5)*asin(sqrt(5)*x/5)) /
            (8*sqrt(-x**2 + 5)), True))


def test_issue_2708():
    # This test needs to use an integration function that can
    # not be evaluated in closed form.  Update as needed.
    f = 1/(a + z + log(z))
    integral_f = NonElementaryIntegral(f, (z, 2, 3))
    assert Integral(f, (z, 2, 3)).doit() == integral_f
    assert integrate(f + exp(z), (z, 2, 3)) == integral_f - exp(2) + exp(3)


def test_issue_8368():
    assert integrate(exp(-s*x)*cosh(x), (x, 0, oo)) == \
        Piecewise(
            (   pi*Piecewise(
                    (   -s/(pi*(-s**2 + 1)),
                        Abs(s**2) < 1),
                    (   1/(pi*s*(1 - 1/s**2)),
                        Abs(s**(-2)) < 1),
                    (   meijerg(
                            ((S(1)/2,), (0, 0)),
                            ((0, S(1)/2), (0,)),
                            polar_lift(s)**2),
                        True)
                ),
                And(
                    Abs(periodic_argument(polar_lift(s)**2, oo)) < pi,
                    cos(Abs(periodic_argument(polar_lift(s)**2, oo))/2)*sqrt(Abs(s**2)) - 1 > 0,
                    Ne(s**2, 1))
            ),
            (
                Integral(exp(-s*x)*cosh(x), (x, 0, oo)),
                True))
    assert integrate(exp(-s*x)*sinh(x), (x, 0, oo)) == \
        Piecewise(
            (   -1/(s + 1)/2 - 1/(-s + 1)/2,
                And(
                    Ne(1/s, 1),
                    Abs(periodic_argument(s, oo)) < pi/2,
                    Abs(periodic_argument(s, oo)) <= pi/2,
                    cos(Abs(periodic_argument(s, oo)))*Abs(s) - 1 > 0)),
            (   Integral(exp(-s*x)*sinh(x), (x, 0, oo)),
                True))


def test_issue_8901():
    assert integrate(sinh(1.0*x)) == 1.0*cosh(1.0*x)
    assert integrate(tanh(1.0*x)) == 1.0*x - 1.0*log(tanh(1.0*x) + 1)
    assert integrate(tanh(x)) == x - log(tanh(x) + 1)


@slow
def test_issue_7130():
    i, L, a, b = symbols('i L a b')
    integrand = (cos(pi*i*x/L)**2 / (a + b*x)).rewrite(exp)
    assert x not in integrate(integrand, (x, 0, L)).free_symbols
