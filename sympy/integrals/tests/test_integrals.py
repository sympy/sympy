from sympy import (S, symbols, integrate, Integral, Derivative, exp, oo, Symbol,
        Function, Rational, log, sin, cos, pi, E, I, Poly, LambertW, diff,
        Matrix, sympify, sqrt, atan, asin, acos, atan, DiracDelta, Heaviside,
        Lambda, sstr, Add)
from sympy.utilities.pytest import XFAIL, skip, raises
from sympy.physics.units import m, s

x,y,a,t = symbols('x,y,a,t')
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

def test_basics():

    assert Integral(0, x) != 0
    assert Integral(x, (x, 1, 1)) != 0
    assert Integral(oo, x) != oo
    assert Integral(S.NaN, x) == S.NaN

    assert diff(Integral(y, y), x)       == 0
    assert diff(Integral(x, (x,0,1)), x) == 0
    assert diff(Integral(x, x), x)       == x
    assert diff(Integral(t, (t,0,x)), x) == x + Integral(0, (t, 0, x))

    e=(t+1)**2
    assert diff(integrate(e, (t,0,x)), x) == \
           diff(Integral(e, (t, 0, x)), x).doit().expand() == \
           ((1+x)**2).expand()
    assert diff(integrate(e, (t,0,x)), t) == \
           diff(Integral(e, (t,0,x)), t) == 0
    assert diff(integrate(e, (t,0,x)), a) == \
           diff(Integral(e, (t, 0, x)), a) == 0
    assert diff(integrate(e, t), a) == diff(Integral(e, t), a) == 0

    assert integrate(e, (t,a,x)).diff(x) == \
           Integral(e, (t, a, x)).diff(x).doit().expand()
    assert Integral(e, (t, a, x)).diff(x).doit() == ((1+x)**2)
    assert integrate(e, (t,x,a)).diff(x).doit() == (-(1+x)**2).expand()

    assert integrate(t**2, (t,x,2*x)).diff(x) == 7*x**2

    assert Integral(x, x).atoms() == set([x])
    assert Integral(f(x), (x, 0, 1)).atoms() == set([S(0), S(1), x])

    assert diff_test(Integral(x, (x, 3*y))) == set([y])
    assert diff_test(Integral(x, (a, 3*y))) == set([x, y])

    # sum integral of terms
    assert integrate(y + x + exp(x), x) == x*y + x**2/2 + exp(x)

    assert Integral(x).is_commutative
    n = Symbol('n', commutative=False)
    assert Integral(x, (x, n)).is_commutative is False
    assert Integral(n + x, x).is_commutative is False

def test_basics_multiple():

    assert diff_test(Integral(x, (x, 3*x, 5*y), (y, x, 2*x))) == set([x])
    assert diff_test(Integral(x, (x, 5*y), (y, x, 2*x))) == set([x])
    assert diff_test(Integral(x, (x, 5*y), (y, y, 2*x))) == set([x, y])
    assert diff_test(Integral(y, y, x)) == set([x, y])
    assert diff_test(Integral(y*x, x, y)) == set([x, y])
    assert diff_test(Integral(x + y, y, (y, 1, x))) == set([x])
    assert diff_test(Integral(x + y, (x, x, y), (y, y, x))) == set([x, y])

def test_integration():
    assert integrate(0, (t,0,x)) == 0
    assert integrate(3, (t,0,x)) == 3*x
    assert integrate(t, (t,0,x)) == x**2/2
    assert integrate(3*t, (t,0,x))== 3*x**2/2
    assert integrate(3*t**2, (t,0,x)) == x**3
    assert integrate(1/t, (t,1,x)) == log(x)
    assert integrate(-1/t**2, (t,1,x)) == 1/x-1
    assert integrate(t**2+5*t-8, (t,0,x)) == x**3/3+5*x**2/2-8*x
    assert integrate(x**2, x) == x**3/3
    assert integrate((3*t*x)**5, x) == (3*t)**5 * x**6 / 6

    b = Symbol("b")
    c = Symbol("c")
    assert integrate(a*t, (t,0,x))==a*x**2/2
    assert integrate(a*t**4, (t,0,x))==a*x**5/5
    assert integrate(a*t**2+b*t+c, (t,0,x))==a*x**3/3+b*x**2/2+c*x

def test_multiple_integration():
    assert integrate((x**2)*(y**2), (x,0,1), (y,-1,2)) == Rational(1)
    assert integrate((y**2)*(x**2), x, y) == Rational(1,9)*(x**3)*(y**3)
    assert integrate(1/(x+3)/(1+x)**3, x) == -S(1)/8*log(3 + x) + S(1)/8*log(1 + x) + x/(4 + 8*x + 4*x**2)

def test_issue433():
    assert integrate(exp(-x), (x,0,oo)) == 1

def test_issue461():
    assert integrate(x**Rational(3,2), x) == 2*x**Rational(5,2)/5
    assert integrate(x**Rational(1,2), x) == 2*x**Rational(3,2)/3
    assert integrate(x**Rational(-3,2), x) == -2*x**Rational(-1,2)

def test_integrate_poly():
    p = Poly(x + x**2*y + y**3, x, y)

    qx = integrate(p, x)
    qy = integrate(p, y)

    assert isinstance(qx, Poly) == True
    assert isinstance(qy, Poly) == True

    assert qx.gens == (x, y)
    assert qy.gens == (x, y)

    assert qx.as_expr() == x**2/2 + x**3*y/3 + x*y**3
    assert qy.as_expr() == x*y + x**2*y**2/2 + y**4/4

def test_integrate_poly_defined():
    p = Poly(x + x**2*y + y**3, x, y)

    Qx = integrate(p, (x, 0, 1))
    Qy = integrate(p, (y, 0, pi))

    assert isinstance(Qx, Poly) == True
    assert isinstance(Qy, Poly) == True

    assert Qx.gens == (y,)
    assert Qy.gens == (x,)

    assert Qx.as_expr() == Rational(1,2) + y/3 + y**3
    assert Qy.as_expr() == pi**4/4 + pi*x + pi**2*x**2/2

def test_integrate_omit_var():
    y = Symbol('y')
    raises(ValueError, "integrate(2)")
    assert integrate(x)     == x**2/2
    assert integrate(x*y)   == x**2*y**2/4

def test_integrate_poly_accurately():
    y = Symbol('y')
    assert integrate(x*sin(y), x)       == x**2*sin(y)/2

    # when passed to risch_norman, this will be a CPU hog, so this really
    # checks, that integrated function is recognized as polynomial
    assert integrate(x**1000*sin(y), x) == x**1001*sin(y)/1001

def test_issue536():
    y = Symbol('y')
    assert integrate(x**2, y) == x**2*y
    assert integrate(x**2, (y, -1, 1)) == 2*x**2

# works in sympy and py.test but hangs in `setup.py test`
def test_integrate_linearterm_pow():
    # check integrate((a*x+b)^c, x)  --  #400
    y = Symbol('y')
    assert integrate(x**y, x) == x**(y+1)/(y+1)
    assert integrate((exp(y)*x + 1/y)**(1+sin(y)), x)   == exp(-y)*(exp(y)*x + 1/y)**(2+sin(y)) / (2+sin(y))

def test_issue519():
    assert integrate(pi*x**Rational(1,2),x) == 2*pi*x**Rational(3,2)/3
    assert integrate(pi*x**Rational(1,2) + E*x**Rational(3,2),x) == \
                                               2*pi*x**Rational(3,2)/3  + \
                                               2*E *x**Rational(5,2)/5
def test_issue524():
    assert integrate(cos((n+1) * x), x)   == sin(x*(n+1)) / (n+1)
    assert integrate(cos((n-1) * x), x)   == sin(x*(n-1)) / (n-1)

    assert integrate(cos((n+1) * x) + cos((n-1) * x), x) == \
                                             sin(x*(n+1)) / (n+1)  + \
                                             sin(x*(n-1)) / (n-1)

def test_issue565():
    assert integrate(-1./2 * x * sin(n * pi * x/2), [x, -2, 0])  == 2*cos(pi*n)/(pi*n)
    assert integrate(-Rational(1)/2 * x * sin(n * pi * x/2), [x, -2, 0]) \
                                                                 == 2*cos(pi*n)/(pi*n)
def test_issue580():
    # definite integration of rational functions gives wrong answers
    assert NS(Integral(1/(x**2-8*x+17), (x, 2, 4))) == '1.10714871779409'

def test_issue587(): # remove this when fresnel itegrals are implemented
    assert integrate(sin(x**2), x) == Integral(sin(x**2), x)

def test_integrate_units():
    assert integrate(x * m/s, (x, 1*s, 5*s)) == 12*m*s

def test_transcendental_functions():
    assert integrate(LambertW(2*x), x) == -x + x*LambertW(2*x) + x/LambertW(2*x)

def test_issue641():
    f=4*log(x)-2*log(x)**2
    fid=diff(integrate(f,x),x)
    assert abs(f.subs(x,42).evalf() - fid.subs(x,42).evalf()) < 1e-10

def test_issue689():
    assert integrate(1/(1+x**2), x) == atan(x)

def test_issue853():
    f = sin(x)
    assert integrate(f, x) == -cos(x)
    raises(ValueError, "integrate(f, 2*x)")

def test_issue1417():
    assert integrate(2**x - 2*x, x) == 2**x/log(2) - x**2

def test_matrices():
    M = Matrix(2, 2, lambda i, j: (i+j+1)*sin((i+j+1)*x))

    assert integrate(M, x) == Matrix([
        [-cos(x),   -cos(2*x)],
        [-cos(2*x), -cos(3*x)],
    ])

# issue1012
def test_integrate_functions():
    assert integrate(f(x), x)       == Integral(f(x), x)
    assert integrate(f(x), (x,0,1)) == Integral(f(x), (x,0,1))
    assert integrate(f(x)*diff(f(x), x), x) == f(x)**2/2
    assert integrate(diff(f(x),x) / f(x),x) == log(f(x))

def test_integrate_derivatives():
    assert integrate(Derivative(f(x), x), x) == f(x)
    assert integrate(Derivative(f(y), y), x) == x*Derivative(f(y), y)

def test_transform():
    a = Integral(x**2+1, (x, -1, 2))
    assert a.doit() == a.transform(x, 3*x+1).doit()
    assert a.transform(x, 3*x+1).transform(x, 3*x+1, inverse=True) == a
    assert a.transform(x, 3*x+1, inverse=True).transform(x, 3*x+1) == a
    a = Integral(sin(1/x), (x, 0, 1))
    assert a.transform(x, 1/x) == Integral(sin(x)/x**2, (x, 1, oo))
    assert a.transform(x, 1/x).transform(x, 1/x) == a
    a = Integral(exp(-x**2), (x, -oo, oo))
    assert a.transform(x, 2*x) == Integral(2*exp(-4*x**2), (x, -oo, oo))
    # < 3 arg limit handled properly
    assert Integral(x, x).transform(x, a*x) == Integral(x*a**2, x)
    raises(ValueError, "a.transform(x, 1/x)")
    raises(ValueError, "a.transform(x, 1/x)")
    _3 = S(3)
    assert Integral(x, (x, 0, -_3)).transform(x, 1/x) == \
    Integral(-1/x**3, (x, -oo, -1/_3))
    assert Integral(x, (x, 0, _3)).transform(x, 1/x) == \
    Integral(x**(-3), (x, 1/_3, oo))

def test_issue953():
    f = S(1)/2*asin(x) + x*(1 - x**2)**(S(1)/2)/2

    assert integrate(cos(asin(x)), x) == f
    assert integrate(sin(acos(x)), x) == f

def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)

def test_evalf_integrals():
    assert NS(Integral(x, (x, 2, 5)), 15) == '10.5000000000000'
    gauss = Integral(exp(-x**2), (x, -oo, oo))
    assert NS(gauss, 15) == '1.77245385090552'
    assert NS(gauss**2 - pi + E*Rational(1,10**20), 15) in ('2.71828182845904e-20', '2.71828182845905e-20')
    # A monster of an integral from http://mathworld.wolfram.com/DefiniteIntegral.html
    t = Symbol('t')
    a = 8*sqrt(3)/(1+3*t**2)
    b = 16*sqrt(2)*(3*t+1)*(4*t**2+t+1)**Rational(3,2)
    c = (3*t**2+1)*(11*t**2+2*t+3)**2
    d = sqrt(2)*(249*t**2+54*t+65)/(11*t**2+2*t+3)**2
    f = a - b/c - d
    assert NS(Integral(f, (t, 0, 1)), 50) == NS((3*sqrt(2)-49*pi+162*atan(sqrt(2)))/12,50)
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(1/x))/(1+x+x**2), (x, 0, 1)), 15) == NS('pi/sqrt(3) * log(2*pi**(5/6) / gamma(1/6))', 15)
    # http://mathworld.wolfram.com/AhmedsIntegral.html
    assert NS(Integral(atan(sqrt(x**2+2))/(sqrt(x**2+2)*(x**2+1)), (x, 0, 1)), 15) == NS(5*pi**2/96, 15)
    # http://mathworld.wolfram.com/AbelsIntegral.html
    assert NS(Integral(x/((exp(pi*x)-exp(-pi*x))*(x**2+1)), (x, 0, oo)), 15) == NS('log(2)/2-1/4',15)
    # Complex part trimming
    # http://mathworld.wolfram.com/VardisIntegral.html
    assert NS(Integral(log(log(sin(x)/cos(x))), (x, pi/4, pi/2)), 15, chop=True) == \
        NS('pi/4*log(4*pi**3/gamma(1/4)**4)', 15)
    #
    # Endpoints causing trouble (rounding error in integration points -> complex log)
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 17, chop=True) == NS(2, 17)
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 20, chop=True) == NS(2, 20)
    assert NS(2+Integral(log(2*cos(x/2)), (x, -pi, pi)), 22, chop=True) == NS(2, 22)
    # Needs zero handling
    assert NS(pi - 4*Integral('sqrt(1-x**2)', (x, 0, 1)), 15, maxn=30, chop=True) in ('0.0', '0')
    # Oscillatory quadrature
    a = Integral(sin(x)/x**2, (x, 1, oo)).evalf(maxn=15)
    assert 0.49 < a < 0.51
    assert NS(Integral(sin(x)/x**2, (x, 1, oo)), quad='osc') == '0.504067061906928'
    assert NS(Integral(cos(pi*x+1)/x, (x, -oo, -1)), quad='osc') == '0.276374705640365'

@XFAIL
def test_evalf_issue_939():
    # http://code.google.com/p/sympy/issues/detail?id=939

    # The output form of an integral may differ by a step function between
    # revisions, making this test a bit useless. This can't be said about
    # other two tests. For now, all values of this evaluation are used here,
    # but in future this should be reconsidered.
    assert NS(integrate(1/(x**5+1), x).subs(x, 4), chop=True) in \
        ['-0.000976138910649103', '0.965906660135753', '1.93278945918216']

    assert NS(Integral(1/(x**5+1), (x, 2, 4))) == '0.0144361088886740'
    assert NS(integrate(1/(x**5+1), (x, 2, 4)), chop=True) == '0.0144361088886740'

def xtest_failing_integrals():
    #---
    # Double integrals not implemented
    assert NS(Integral(sqrt(x)+x*y, (x, 1, 2), (y, -1, 1)), 15) == '2.43790283299492'
    # double integral + zero detection
    assert NS(Integral(sin(x+x*y), (x, -1, 1), (y, -1, 1)), 15) == '0.0'

def test_integrate_DiracDelta():
    assert integrate(DiracDelta(x),x) == Heaviside(x)
    assert integrate(DiracDelta(x) * f(x),x) == f(0) * Heaviside(x)
    assert integrate(DiracDelta(x) * f(x),(x,-oo,oo)) == f(0)
    assert integrate(DiracDelta(x) * f(x),(x,0,oo)) == f(0)/2
    assert integrate(DiracDelta(x**2+x-2),x) - \
           (Heaviside(-1 + x)/3 + Heaviside(2 + x)/3) == 0
    assert integrate(cos(x)*(DiracDelta(x)+DiracDelta(x**2-1))*sin(x)*(x-pi),x) - \
           (-pi*(cos(1)*Heaviside(-1 + x)*sin(1)/2 - cos(1)*Heaviside(1 + x)*sin(1)/2) + \
           cos(1)*Heaviside(1 + x)*sin(1)/2 + cos(1)*Heaviside(-1 + x)*sin(1)/2) == 0


def test_subs1():
    e = Integral(exp(x-y), x)
    assert e.subs(y, 3) == Integral(exp(x-3), x)
    e = Integral(exp(x-y), (x, 0, 1))
    assert e.subs(y, 3) == Integral(exp(x-3), (x, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x-y)*f(y), (y, -oo, oo))
    assert conv.subs({x:0}) == Integral(exp(-2*y**2), (y, -oo, oo))

def test_subs2():
    e = Integral(exp(x-y), x, t)
    assert e.subs(y, 3) == Integral(exp(x-3), x, t)
    e = Integral(exp(x-y), (x, 0, 1), (t, 0, 1))
    assert e.subs(y, 3) == Integral(exp(x-3), (x, 0, 1), (t, 0, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x-y)*f(y), (y, -oo, oo), (t, 0, 1))
    assert conv.subs({x:0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))

def test_subs3():
    e = Integral(exp(x-y), (x, 0, y), (t, y, 1))
    assert e.subs(y, 3) == Integral(exp(x-3), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(x-y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x:0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))

def test_subs4():
    e = Integral(exp(x), (x, 0, y), (t, y, 1))
    assert e.subs(y, 3) == Integral(exp(x), (x, 0, 3), (t, 3, 1))
    f = Lambda(x, exp(-x**2))
    conv = Integral(f(y)*f(y), (y, -oo, oo), (t, x, 1))
    assert conv.subs({x:0}) == Integral(exp(-2*y**2), (y, -oo, oo), (t, 0, 1))

def test_subs5():
    e = Integral(exp(-x**2), x)
    assert e.subs(x, 5) == Integral(exp(-x**2), (x, 5))
    e = Integral(exp(-x**2), (x, -oo, oo))
    assert e.subs(x, 5) == e
    e = Integral(exp(-x**2+y), x)
    assert e.subs(x, 5) == Integral(exp(y - x**2), (x, 5))
    assert e.subs(y, 5) == Integral(exp(-x**2+5), x)
    e = Integral(exp(-x**2+y), (y, -oo, oo), (x, -oo, oo))
    assert e.subs(x, 5) == e
    assert e.subs(y, 5) == e

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
    assert e.subs({x:1, y:2}) == e
    e = Integral(sin(x) + sin(y), (x, sin(x), sin(y)),
                                  (y, 1, 2))
    assert e._eval_subs(sin(y), 1) == e
    assert e._eval_subs(sin(x), 1) == Integral(sin(x) + sin(y), (x, 1, sin(y)),
                                                                (y, 1, 2))

def test_integration_variable():
    raises(ValueError, "Integral(exp(-x**2), 3)")
    raises(ValueError, "Integral(exp(-x**2), (3, -oo, oo))")

def test_expand_integral():
    assert Integral(cos(x**2)*(sin(x**2)+1),(x, 0, 1)).expand() == Integral(cos(x**2)*sin(x**2) + cos(x**2), (x, 0, 1))
    assert Integral(cos(x**2)*(sin(x**2)+1),x).expand() == Integral(cos(x**2)*sin(x**2) + cos(x**2), x)

def test_as_sum_midpoint1():
    e = Integral(sqrt(x**3+1), (x, 2, 10))
    assert e.as_sum(1, method="midpoint") == 8*217**(S(1)/2)
    assert e.as_sum(2, method="midpoint") == 4*65**(S(1)/2) + 12*57**(S(1)/2)
    assert e.as_sum(3, method="midpoint") == 8*217**(S(1)/2)/3 + \
            8*3081**(S(1)/2)/27 + 8*52809**(S(1)/2)/27
    assert e.as_sum(4, method="midpoint") == 2*730**(S(1)/2) + \
            4*7**(S(1)/2) + 4*86**(S(1)/2) + 6*14**(S(1)/2)
    assert abs(e.as_sum(4, method="midpoint").n() - e.n()) < 0.5

    e = Integral(sqrt(x**3+y**3), (x, 2, 10), (y, 0, 10))
    raises(NotImplementedError, "e.as_sum(4)")

def test_as_sum_midpoint2():
    e = Integral((x+y)**2, (x, 0, 1))
    assert e.as_sum(1, method="midpoint").expand() == S(1)/4 + y + y**2
    assert e.as_sum(2, method="midpoint").expand() == S(5)/16 + y + y**2
    assert e.as_sum(3, method="midpoint").expand() == S(35)/108 + y + y**2
    assert e.as_sum(4, method="midpoint").expand() == S(21)/64 + y + y**2

def test_as_sum_left():
    e = Integral((x+y)**2, (x, 0, 1))
    assert e.as_sum(1, method="left").expand() == y**2
    assert e.as_sum(2, method="left").expand() == S(1)/8 + y/2 + y**2
    assert e.as_sum(3, method="left").expand() == S(5)/27 + 2*y/3 + y**2
    assert e.as_sum(4, method="left").expand() == S(7)/32 + 3*y/4 + y**2

def test_as_sum_right():
    e = Integral((x+y)**2, (x, 0, 1))
    assert e.as_sum(1, method="right").expand() == 1 + 2*y + y**2
    assert e.as_sum(2, method="right").expand() == S(5)/8 + 3*y/2 + y**2
    assert e.as_sum(3, method="right").expand() == S(14)/27 + 4*y/3 + y**2
    assert e.as_sum(4, method="right").expand() == S(15)/32 + 5*y/4 + y**2

def test_as_sum_raises():
    e = Integral((x+y)**2, (x, 0, 1))
    raises(ValueError, "e.as_sum(-1)")
    raises(ValueError, "e.as_sum(0)")
    raises(NotImplementedError, "e.as_sum(oo)")
    raises(NotImplementedError, "e.as_sum(3, method='xxxx2')")

def test_nested_doit():
    e = Integral(Integral(x, x), x)
    f = Integral(x, x, x)
    assert e.doit() == f.doit()

def test_issue1566():
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

def test_doit():
    e = Integral(Integral(2*x), (x, 0, 1))
    assert e.doit() == Rational(1, 3)
    f = Function('f')
    # doesn't matter if the integral can't be performed
    assert Integral(f(x), (x, 1, 1)).doit() == 0
    # doesn't matter if the limits can't be evaluated
    assert Integral(0, (x, 1, Integral(f(x), x))).doit() == 0

@XFAIL
def test_doit2():
    e = Integral(Integral(2*x), (x, 0, 1))
    # risch currently chokes on the contained integral.
    assert e.doit(deep = False) == e

def issue_1785():
    assert integrate(sqrt(x)*(1+x)) == 2*x**Rational(3, 2)/3 + 2*x**Rational(5, 2)/5
    assert integrate(x**x*(1+log(x))) is not None

@XFAIL
def issue_1785_fail():
    assert integrate(x**x*(1+log(x)).expand(mul=True)) is None

def test_is_number():
    from sympy.abc import x, y, z
    from sympy import cos, sin
    assert Integral(x).is_number == False
    assert Integral(1, x).is_number == False
    assert Integral(1, (x, 1)).is_number == True
    assert Integral(1, (x, 1, 2)).is_number == True
    assert Integral(1, (x, 1, y)).is_number == False
    assert Integral(x, y).is_number == False
    assert Integral(x, (y, 1, x)).is_number == False
    assert Integral(x, (y, 1, 2)).is_number == False
    assert Integral(x, (x, 1, 2)).is_number == True
    assert Integral(x, (y, 1, 1)).is_number == True
    assert Integral(x*y, (x, 1, 2), (y, 1, 3)).is_number == True
    assert Integral(x*y, (x, 1, 2), (y, 1, z)).is_number == False
    assert Integral(x, (x, 1)).is_number == True
    assert Integral(x, (x, 1, Integral(y, (y, 1, 2)))).is_number == True
    # it is possible to get a false negative if the integrand is
    # actually an unsimplified zero, but this is true of is_number in general.
    assert Integral(sin(x)**2 + cos(x)**2 - 1, x).is_number == False

def test_symbols():
    from sympy.abc import x, y, z
    assert Integral(0, x).free_symbols == set()
    assert Integral(x).free_symbols == set([x])
    assert Integral(x, (x, None, y)).free_symbols == set([y])
    assert Integral(x, (x, y, None)).free_symbols == set([y])
    assert Integral(x, (x, 1, y)).free_symbols == set([y])
    assert Integral(x, (x, y, 1)).free_symbols == set([y])
    assert Integral(x, (x, x, y)).free_symbols == set([x, y])
    assert Integral(x, x, y).free_symbols == set([x, y])
    assert Integral(x, (x, 1, 2)).free_symbols == set()
    assert Integral(x, (y, 1, 2)).free_symbols == set([x])
    assert Integral(x, (y, z, z)).free_symbols == set()
    assert Integral(x, (y, 1, 2), (y, None, None)).free_symbols == set([x, y])
    assert Integral(x, (y, 1, 2), (x, 1, y)).free_symbols == set([y])
    assert Integral(2, (y, 1, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (y, x, 2), (y, 1, x), (x, 1, 2)).free_symbols == set()
    assert Integral(2, (x, 1, 2), (y, x, 2), (y, 1, 2)).free_symbols == set([x])

def test_is_zero():
    from sympy.abc import x, m, n
    assert Integral(0, (x, 1, x)).is_zero
    assert Integral(1, (x, 1, 1)).is_zero
    assert Integral(1, (x, 1, 2)).is_zero is False
    assert Integral(sin(m*x)*cos(n*x), (x, 0, 2*pi)).is_zero is None

def test_series():
    from sympy.abc import x
    i = Integral(cos(x))
    e = i.lseries(x)
    assert i.nseries(x, n=8).removeO() == Add(*[e.next() for j in range(4)])
