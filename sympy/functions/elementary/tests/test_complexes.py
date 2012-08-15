from sympy import (symbols, Symbol, sqrt, oo, re, nan, im, sign, I, E, log,
        pi, arg, conjugate, expand, exp, sin, cos, Function, Abs, zoo, atan2,
        S, DiracDelta, Rational, Heaviside)
from sympy.utilities.pytest import XFAIL

from sympy.utilities.randtest import comp
def N_equals(a, b):
    """Check whether two complex numbers are numerically close"""
    return comp(a.n(), b.n(), 1.e-6)

def test_re():
    x, y = symbols('x,y')
    a, b = symbols('a,b', real=True)

    r = Symbol('r', real=True)
    i = Symbol('i', imaginary=True)

    assert re(nan) == nan

    assert re(oo) == oo
    assert re(-oo) == -oo

    assert re(0) == 0

    assert re(1) == 1
    assert re(-1) == -1

    assert re(E) == E
    assert re(-E) == -E

    assert re(x) == re(x)
    assert re(x*I) == -im(x)
    assert re(r*I) == 0
    assert re(r) == r
    assert re(i*I) == I * i
    assert re(i) == 0

    assert re(x + y) == re(x + y)
    assert re(x + r) == re(x) + r

    assert re(re(x)) == re(x)

    assert re(2 + I) == 2
    assert re(x + I) == re(x)

    assert re(x + y*I) == re(x) - im(y)
    assert re(x + r*I) == re(x)

    assert re(log(2*I)) == log(2)

    assert re((2+I)**2).expand(complex=True) == 3

    assert re(conjugate(x)) == re(x)
    assert conjugate(re(x)) == re(x)

    assert re(x).as_real_imag() == (re(x), 0)

    assert re(i*r*x).diff(r) == re(i*x)
    assert re(i*r*x).diff(i) == -I * im(r*x)

    assert re(sqrt(a + b*I)) == (a**2 + b**2)**Rational(1,4)*cos(atan2(b, a)/2)
    assert re(a * (2 + b*I)) == 2*a

    assert re((1 + sqrt(a + b*I))/2) == \
           (a**2 + b**2)**Rational(1,4)*cos(atan2(b, a)/2)/2 + Rational(1,2)

def test_im():
    x, y = symbols('x,y')
    a, b = symbols('a,b', real=True)

    r = Symbol('r', real=True)
    i = Symbol('i', imaginary=True)

    assert im(nan) == nan

    assert im(oo*I) == oo
    assert im(-oo*I) == -oo

    assert im(0) == 0

    assert im(1) == 0
    assert im(-1) == 0

    assert im(E*I) == E
    assert im(-E*I) == -E

    assert im(x) == im(x)
    assert im(x*I) == re(x)
    assert im(r*I) == r
    assert im(r) == 0
    assert im(i*I) == 0
    assert im(i) == -I * i

    assert im(x + y) == im(x + y)
    assert im(x + r) == im(x)
    assert im(x + r*I) == im(x) + r

    assert im(im(x)*I) == im(x)

    assert im(2 + I) == 1
    assert im(x + I) == im(x) + 1

    assert im(x + y*I) == im(x) + re(y)
    assert im(x + r*I) == im(x) + r

    assert im(log(2*I)) == pi/2

    assert im((2+I)**2).expand(complex=True) == 4

    assert im(conjugate(x)) == -im(x)
    assert conjugate(im(x)) == im(x)

    assert im(x).as_real_imag() == (im(x), 0)

    assert im(i*r*x).diff(r) == im(i*x)
    assert im(i*r*x).diff(i) == -I * re(r*x)

    assert im(sqrt(a + b*I)) == (a**2 + b**2)**Rational(1,4)*sin(atan2(b, a)/2)
    assert im(a * (2 + b*I)) == a*b

    assert im((1 + sqrt(a + b*I))/2) == \
           (a**2 + b**2)**Rational(1,4)*sin(atan2(b, a)/2)/2

def test_sign():
    assert sign(1.2) == 1
    assert sign(-1.2) == -1
    assert sign(3*I) == I
    assert sign(-3*I) == -I
    assert sign(0) == 0
    assert sign(nan) == nan

    x = Symbol('x')
    assert sign(x).is_zero == None
    assert sign(x).doit() == sign(x)
    assert sign(1.2*x) == sign(x)
    assert sign(2*x) == sign(x)
    assert sign(I*x) == I*sign(x)
    assert sign(-2*I*x) == -I*sign(x)
    assert sign(conjugate(x)) == conjugate(sign(x))

    p = Symbol('p', positive = True)
    n = Symbol('n', negative = True)
    m = Symbol('m', negative = True)
    assert sign(2*p*x) == sign(x)
    assert sign(n*x) == -sign(x)
    assert sign(n*m*x) == sign(x)

    x = Symbol('x', imaginary=True)
    assert sign(x).is_zero == False
    assert sign(x).diff(x) == 2*DiracDelta(-I*x)
    assert sign(x).doit() == x / Abs(x)
    assert conjugate(sign(x)) == -sign(x)

    x = Symbol('x', real=True)
    assert sign(x).is_zero == None
    assert sign(x).diff(x) == 2*DiracDelta(x)
    assert sign(x).doit() == sign(x)
    assert conjugate(sign(x)) == sign(x)

    x = Symbol('x', nonzero=True)
    assert sign(x).is_zero == False
    assert sign(x).doit() == x / Abs(x)
    assert sign(Abs(x)) == 1
    assert Abs(sign(x)) == 1

    x = Symbol('x', positive=True)
    assert sign(x).is_zero == False
    assert sign(x).doit() == x / Abs(x)
    assert sign(Abs(x)) == 1
    assert Abs(sign(x)) == 1

    x = 0
    assert sign(x).is_zero == True
    assert sign(x).doit() == 0
    assert sign(Abs(x)) == 0
    assert Abs(sign(x)) == 0

    nz = Symbol('nz', nonzero=True, integer=True)
    assert sign(nz)**2 == 1
    assert (sign(nz)**3).args == (sign(nz), 3)

def test_as_real_imag():
    n = pi**1000
    # the special code for working out the real
    # and complex parts of a power with Integer exponent
    # should not run if there is no imaginary part, hence
    # this should not hang
    assert n.as_real_imag() == (n, 0)

    # issue 3162
    x = Symbol('x')
    assert sqrt(x).as_real_imag() == \
    ((re(x)**2 + im(x)**2)**(S(1)/4)*cos(atan2(im(x), re(x))/2), \
     (re(x)**2 + im(x)**2)**(S(1)/4)*sin(atan2(im(x), re(x))/2))

    # issue 754
    a, b = symbols('a,b', real=True)
    assert ((1 + sqrt(a + b*I))/2).as_real_imag() == \
           ((a**2 + b**2)**Rational(1,4)*cos(atan2(b, a)/2)/2 + Rational(1,2), \
            (a**2 + b**2)**Rational(1,4)*sin(atan2(b, a)/2)/2)

@XFAIL
def test_sign_issue_3068():
    n = pi**1000
    i = int(n)
    assert (n - i).round() == 1 # doesn't hang
    assert sign(n - i) == 1
    # perhaps it's not possible to get the sign right when
    # only 1 digit is being requested for this situtation;
    # 2 digits works
    assert (n - x).n(1, subs={x: i}) > 0
    assert (n - x).n(2, subs={x: i}) > 0

def test_Abs():
    x, y = symbols('x,y')
    assert sign(sign(x)) == sign(x)
    assert sign(x*y).func is sign
    assert Abs(0) == 0
    assert Abs(1) == 1
    assert Abs(-1) == 1
    assert Abs(I) == 1
    assert Abs(-I) == 1
    assert Abs(nan) == nan
    assert Abs(I * pi) == pi
    assert Abs(-I * pi) == pi
    assert Abs(I * x) == Abs(x)
    assert Abs(-I * x) == Abs(x)
    assert Abs(-2*x) == 2*Abs(x)
    assert Abs(-2.0*x) == 2.0*Abs(x)
    assert Abs(2*pi*x*y) == 2*pi*Abs(x*y)
    assert Abs(conjugate(x)) == Abs(x)
    assert conjugate(Abs(x)) == Abs(x)

    a = Symbol('a', positive=True)
    assert Abs(2*pi*x*a) == 2*pi*a*Abs(x)
    assert Abs(2*pi*I*x*a) == 2*pi*a*Abs(x)

    x = Symbol('x', real=True)
    n = Symbol('n', integer=True)
    assert x**(2*n) == Abs(x)**(2*n)
    assert Abs(x).diff(x) == sign(x)
    assert abs(x) == Abs(x) # Python built-in
    assert Abs(x)**3 == x**2*Abs(x)
    assert Abs(x)**4 == x**4
    assert (Abs(x)**(3*n)).args == (Abs(x), 3*n) # leave symbolic odd unchanged
    assert (1/Abs(x)).args == (Abs(x), -1)
    assert 1/Abs(x)**3 == 1/(x**2*Abs(x))

    x = Symbol('x', imaginary=True)
    assert Abs(x).diff(x) == -sign(x)

def test_Abs_rewrite():
    x = Symbol('x', real=True)
    a = Abs(x).rewrite(Heaviside).expand()
    assert a == x*Heaviside(x) - x*Heaviside(-x)
    for i in [-2, -1, 0, 1, 2]:
        assert a.subs(x, i) == abs(i)
    y = Symbol('y')
    assert Abs(y).rewrite(Heaviside) == Abs(y)

def test_Abs_real():
    # test some properties of abs that only apply
    # to real numbers
    x = Symbol('x', complex=True)
    assert sqrt(x**2) != Abs(x)
    assert Abs(x**2) != x**2

    x = Symbol('x', real=True)
    assert sqrt(x**2) == Abs(x)
    assert Abs(x**2) == x**2

    # if the symbol is zero, the following will still apply
    nn = Symbol('nn', nonnegative=True, real=True)
    np = Symbol('np', nonpositive=True, real=True)
    assert Abs(nn) == nn
    assert Abs(np) == -np

def test_Abs_properties():
    x = Symbol('x')
    assert Abs(x).is_real == True
    assert Abs(x).is_positive == None
    assert Abs(x).is_nonnegative == True

    w = Symbol('w', complex=True, zero=False)
    assert Abs(w).is_real == True
    assert Abs(w).is_positive == True
    assert Abs(w).is_zero == False

    q = Symbol('q', positive=True)
    assert Abs(q).is_real == True
    assert Abs(q).is_positive == True
    assert Abs(q).is_zero == False

def test_abs():
    # this tests that abs calls Abs; don't rename to
    # test_Abs since that test is already above
    a = Symbol('a', positive=True)
    assert abs(I*(1 + a)**2) == (1 + a)**2

def test_arg():
    assert arg(0) == nan
    assert arg(1) == 0
    assert arg(-1) == pi
    assert arg(I) == pi/2
    assert arg(-I) == -pi/2
    assert arg(1+I) == pi/4
    assert arg(-1+I) == 3*pi/4
    assert arg(1-I) == -pi/4

    p = Symbol('p', positive=True)
    assert arg(p) == 0

    n = Symbol('n', negative=True)
    assert arg(n) == pi

    x = Symbol('x')
    assert conjugate(arg(x)) == arg(x)

def test_conjugate():
    a = Symbol('a', real=True)
    assert conjugate(a) == a
    assert conjugate(I*a) == -I*a

    x, y = symbols('x,y')
    assert conjugate(conjugate(x)) == x
    assert conjugate(x + y) == conjugate(x) + conjugate(y)
    assert conjugate(x - y) == conjugate(x) - conjugate(y)
    assert conjugate(x * y) == conjugate(x) * conjugate(y)
    assert conjugate(x / y) == conjugate(x) / conjugate(y)
    assert conjugate(-x) == -conjugate(x)

def test_issue936():
    x = Symbol('x')
    assert Abs(x).expand(trig=True)     == Abs(x)
    assert sign(x).expand(trig=True)    == sign(x)
    assert arg(x).expand(trig=True)     == arg(x)

def test_issue3206():
    x = Symbol('x')
    assert Abs(Abs(x)) == Abs(x)

def test_issue1655_derivative_conjugate():
    x = Symbol('x', real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert (f(x).conjugate()).diff(x) == (f(x).diff(x)).conjugate()
    assert (f(y).conjugate()).diff(y) == -(f(y).diff(y)).conjugate()

def test_derivatives_issue1658():
    x = Symbol('x', real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert re(f(x)).diff(x) == re(f(x).diff(x))
    assert im(f(x)).diff(x) == im(f(x).diff(x))
    assert re(f(y)).diff(y) == -I*im(f(y).diff(y))
    assert im(f(y)).diff(y) == -I*re(f(y).diff(y))
    assert Abs(f(x)).diff(x).subs(f(x), 1+I*x).doit() == x/sqrt(1 + x**2)
    assert arg(f(x)).diff(x).subs(f(x), 1+I*x**2).doit() == 2*x/(1+x**4)
    assert Abs(f(y)).diff(y).subs(f(y), 1+y).doit() == -y/sqrt(1 - y**2)
    assert arg(f(y)).diff(y).subs(f(y), I+y**2).doit() == 2*y/(1 + y**4)

def test_periodic_argument():
    from sympy import (periodic_argument, unbranched_argument, oo,
                       principal_branch, polar_lift, pi)
    x = Symbol('x')
    p = Symbol('p', positive = True)

    assert unbranched_argument(2 + I) == periodic_argument(2 + I, oo)
    assert unbranched_argument(1 + x) == periodic_argument(1 + x, oo)
    assert N_equals(unbranched_argument((1+I)**2), pi/2)
    assert N_equals(unbranched_argument((1-I)**2), -pi/2)
    assert N_equals(periodic_argument((1+I)**2, 3*pi), pi/2)
    assert N_equals(periodic_argument((1-I)**2, 3*pi), -pi/2)

    assert unbranched_argument(principal_branch(x, pi)) \
           == periodic_argument(x, pi)

    assert unbranched_argument(polar_lift(2 + I)) == unbranched_argument(2 + I)
    assert periodic_argument(polar_lift(2 + I), 2*pi) \
           == periodic_argument(2 + I, 2*pi)
    assert periodic_argument(polar_lift(2 + I), 3*pi) \
           == periodic_argument(2 + I, 3*pi)
    assert periodic_argument(polar_lift(2 + I), pi) \
           == periodic_argument(polar_lift(2 + I), pi)

    assert unbranched_argument(polar_lift(1 + I)) == pi/4
    assert periodic_argument(2*p, p) == periodic_argument(p, p)
    assert periodic_argument(pi*p, p) == periodic_argument(p, p)

@XFAIL
def test_principal_branch_fail():
    # TODO XXX why does abs(x)._eval_evalf() not fall back to global evalf?
    assert N_equals(principal_branch((1 + I)**2, pi/2), 0)

def test_principal_branch():
    from sympy import principal_branch, polar_lift, exp_polar
    p = Symbol('p', positive=True)
    x = Symbol('x')
    neg = Symbol('x', negative=True)

    assert principal_branch(polar_lift(x), p) == principal_branch(x, p)
    assert principal_branch(polar_lift(2 + I), p) == principal_branch(2 + I, p)
    assert principal_branch(2*x, p) == 2*principal_branch(x, p)
    assert principal_branch(1, pi) == exp_polar(0)
    assert principal_branch(-1, 2*pi) == exp_polar(I*pi)
    assert principal_branch(-1, pi) == exp_polar(0)
    assert principal_branch(exp_polar(3*pi*I)*x, 2*pi) == \
           principal_branch(exp_polar(I*pi)*x, 2*pi)
    assert principal_branch(neg*exp_polar(pi*I), 2*pi) == neg*exp_polar(-I*pi)

    assert N_equals(principal_branch((1 + I)**2, 2*pi), 2*I)
    assert N_equals(principal_branch((1 + I)**2, 3*pi), 2*I)
    assert N_equals(principal_branch((1 + I)**2, 1*pi), 2*I)

    # test argument sanitization
    assert principal_branch(x, I).func is principal_branch
    assert principal_branch(x, -4).func is principal_branch
    assert principal_branch(x, -oo).func is principal_branch
    assert principal_branch(x, zoo).func is principal_branch
