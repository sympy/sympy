from sympy import (symbols, expand, expand_func, erf, nan, oo, Float, conjugate,
                   sqrt, sin, cos, pi, re, im, Abs, O, factorial, exp_polar,
                   polar_lift, Symbol, I, integrate, exp, uppergamma, expint,
                   log, loggamma, limit, hyper, meijerg, gamma, S, Shi, Chi,
                   Si, Ci, E1, Ei, sin, cos, sinh, cosh, fresnels, fresnelc)

from sympy.functions.special.error_functions import _erfs

from sympy.core.function import ArgumentIndexError

from sympy.utilities.pytest import raises

x, y, z = symbols('x,y,z')
w = Symbol("w", real=True)
n = Symbol("n", integer=True)

def test_erf():
    assert erf(nan) == nan

    assert erf(oo) == 1
    assert erf(-oo) == -1

    assert erf(0) == 0

    assert erf(I*oo) == oo*I
    assert erf(-I*oo) == -oo*I

    assert erf(-2) == -erf(2)
    assert erf(-x*y) == -erf(x*y)
    assert erf(-x - y) == -erf(x + y)

    assert erf(I).is_real == False
    assert erf(0).is_real == True

    assert conjugate(erf(z)) == erf(conjugate(z))

    assert erf(x).as_leading_term(x) == 2*x/sqrt(pi)
    assert erf(1/x).as_leading_term(x) == erf(1/x)

    assert erf(z).rewrite('uppergamma') == sqrt(z**2)*erf(sqrt(z**2))/z

    assert limit(exp(x)*exp(x**2)*(erf(x+1/exp(x))-erf(x)), x, oo) == 2/sqrt(pi)
    assert limit((1-erf(z))*exp(z**2)*z, z, oo) == 1/sqrt(pi)
    assert limit((1-erf(x))*exp(x**2)*sqrt(pi)*x, x, oo) == 1
    assert limit(((1-erf(x))*exp(x**2)*sqrt(pi)*x-1)*2*x**2, x, oo) == -1

    raises(ArgumentIndexError, lambda: erf(x).fdiff(2))

def test_erf_series():
    assert erf(x).series(x, 0, 7) == 2*x/sqrt(pi) - \
        2*x**3/3/sqrt(pi) + x**5/5/sqrt(pi) + O(x**7)

def test_erf_evalf():
    assert abs( erf(Float(2.0)) - 0.995322265 )  <  1E-8  # XXX

def test__erfs():
    assert _erfs(z).diff(z) == -2/sqrt(S.Pi)+2*z*_erfs(z)

    assert _erfs(1/z).series(z) == z/sqrt(pi) - z**3/(2*sqrt(pi)) + 3*z**5/(4*sqrt(pi)) + O(z**6)

    assert expand(erf(z).rewrite('tractable').diff(z).rewrite('intractable')) == erf(z).diff(z)
    assert _erfs(z).rewrite("intractable") == (-erf(z) + 1)*exp(z**2)

# NOTE we multiply by exp_polar(I*pi) and need this to be on the principal
#      branch, hence take x in the lower half plane (d=0).
def mytn(expr1, expr2, expr3, x, d=0):
    from sympy.utilities.randtest import test_numerically, random_complex_number
    subs = {}
    for a in expr1.free_symbols:
        if a != x:
            subs[a] = random_complex_number()
    return expr2 == expr3 and test_numerically(expr1.subs(subs),
                                               expr2.subs(subs), x, d=d)

def mytd(expr1, expr2, x):
    from sympy.utilities.randtest import test_derivative_numerically, \
                                         random_complex_number
    subs = {}
    for a in expr1.free_symbols:
        if a != x:
            subs[a] = random_complex_number()
    return expr1.diff(x) == expr2 and test_derivative_numerically(expr1.subs(subs), x)

def tn_branch(func, s=None):
    from sympy import I, pi, exp_polar
    from random import uniform
    def fn(x):
        if s is None:
            return func(x)
        return func(s, x)
    c = uniform(1, 5)
    expr = fn(c*exp_polar(I*pi)) - fn(c*exp_polar(-I*pi))
    eps = 1e-15
    expr2 = fn(-c + eps*I) - fn(-c - eps*I)
    return abs(expr.n() - expr2.n()).n() < 1e-10

def test_ei():
    pos = Symbol('p', positive=True)
    neg = Symbol('n', negative=True)
    assert Ei(-pos) == Ei(polar_lift(-1)*pos) - I*pi
    assert Ei(neg) == Ei(polar_lift(neg)) - I*pi
    assert tn_branch(Ei)
    assert mytd(Ei(x), exp(x)/x, x)
    assert mytn(Ei(x), Ei(x).rewrite(uppergamma),
                -uppergamma(0, x*polar_lift(-1)) - I*pi, x)
    assert mytn(Ei(x), Ei(x).rewrite(expint),
                -expint(1, x*polar_lift(-1)) - I*pi, x)
    assert Ei(x).rewrite(expint).rewrite(Ei) == Ei(x)
    assert Ei(x*exp_polar(2*I*pi)) == Ei(x) + 2*I*pi
    assert Ei(x*exp_polar(-2*I*pi)) == Ei(x) - 2*I*pi

    assert mytn(Ei(x), Ei(x).rewrite(Shi), Chi(x) + Shi(x), x)
    assert mytn(Ei(x*polar_lift(I)), Ei(x*polar_lift(I)).rewrite(Si),
                Ci(x) + I*Si(x) + I*pi/2, x)

def test_expint():
    assert mytn(expint(x, y), expint(x, y).rewrite(uppergamma),
                y**(x-1)*uppergamma(1 - x, y), x)
    assert mytd(expint(x, y), -y**(x - 1)*meijerg([], [1, 1], [0, 0, 1-x], [], y), x)
    assert mytd(expint(x, y), -expint(x - 1, y), y)
    assert mytn(expint(1, x), expint(1, x).rewrite(Ei),
                -Ei(x*polar_lift(-1)) + I*pi, x)

    assert expint(-4, x) == exp(-x)/x + 4*exp(-x)/x**2 + 12*exp(-x)/x**3 \
                            + 24*exp(-x)/x**4 + 24*exp(-x)/x**5
    assert expint(-S(3)/2, x) == \
           exp(-x)/x + 3*exp(-x)/(2*x**2) - 3*sqrt(pi)*erf(sqrt(x))/(4*x**S('5/2')) \
           + 3*sqrt(pi)/(4*x**S('5/2'))

    assert tn_branch(expint, 1)
    assert tn_branch(expint, 2)
    assert tn_branch(expint, 3)
    assert tn_branch(expint, 1.7)
    assert tn_branch(expint, pi)

    assert expint(y, x*exp_polar(2*I*pi)) \
           == x**(y - 1)*(exp(2*I*pi*y) - 1)*gamma(-y + 1) + expint(y, x)
    assert expint(y, x*exp_polar(-2*I*pi)) \
           == x**(y - 1)*(exp(-2*I*pi*y) - 1)*gamma(-y + 1) + expint(y, x)
    assert expint(2, x*exp_polar(2*I*pi)) == 2*I*pi*x + expint(2, x)
    assert expint(2, x*exp_polar(-2*I*pi)) == -2*I*pi*x + expint(2, x)
    assert expint(1,x).rewrite(Ei).rewrite(expint) == expint(1, x)

    assert mytn(E1(x), E1(x).rewrite(Shi), Shi(x) - Chi(x), x)
    assert mytn(E1(polar_lift(I)*x), E1(polar_lift(I)*x).rewrite(Si),
                -Ci(x) + I*Si(x) - I*pi/2, x)

    assert mytn(expint(2, x), expint(2, x).rewrite(Ei).rewrite(expint),
                -x*E1(x) + exp(-x), x)
    assert mytn(expint(3, x), expint(3, x).rewrite(Ei).rewrite(expint),
                x**2*E1(x)/2 + (1-x)*exp(-x)/2, x)

def tn_arg(func):
    def test(arg, e1, e2):
        from random import uniform
        v = uniform(1, 5)
        v1 = func(arg*x).subs(x, v).n()
        v2 = func(e1*v + e2*1e-15).n()
        return abs(v1 - v2).n() < 1e-10
    return test(exp_polar(I*pi/2), I, 1) and \
           test(exp_polar(-I*pi/2), -I, 1) and \
           test(exp_polar(I*pi), -1, I) and \
           test(exp_polar(-I*pi), -1, -I)

def test_si():
    assert Si(I*x) == I*Shi(x)
    assert Shi(I*x) == I*Si(x)
    assert Si(-I*x) == -I*Shi(x)
    assert Shi(-I*x) == -I*Si(x)
    assert Si(-x) == -Si(x)
    assert Shi(-x) == -Shi(x)
    assert Si(exp_polar(2*pi*I)*x) == Si(x)
    assert Si(exp_polar(-2*pi*I)*x) == Si(x)
    assert Shi(exp_polar(2*pi*I)*x) == Shi(x)
    assert Shi(exp_polar(-2*pi*I)*x) == Shi(x)

    assert mytd(Si(x), sin(x)/x, x)
    assert mytd(Shi(x), sinh(x)/x, x)

    assert mytn(Si(x), Si(x).rewrite(Ei),
                -I*(-Ei(x*exp_polar(-I*pi/2))/2 \
                        + Ei(x*exp_polar(I*pi/2))/2 - I*pi) + pi/2, x)
    assert mytn(Si(x), Si(x).rewrite(expint),
                -I*(-expint(1, x*exp_polar(-I*pi/2))/2 + \
                    expint(1, x*exp_polar(I*pi/2))/2) + pi/2, x)
    assert mytn(Shi(x), Shi(x).rewrite(Ei),
                Ei(x)/2 - Ei(x*exp_polar(I*pi))/2 + I*pi/2, x)
    assert mytn(Shi(x), Shi(x).rewrite(expint),
                expint(1, x)/2 - expint(1, x*exp_polar(I*pi))/2 - I*pi/2, x)

    assert tn_arg(Si)
    assert tn_arg(Shi)

    from sympy import O
    assert Si(x).nseries(x, n=8) == x - x**3/18 + x**5/600 - x**7/35280 + O(x**9)
    assert Shi(x).nseries(x, n=8) == x + x**3/18 + x**5/600 + x**7/35280 + O(x**9)
    assert Si(sin(x)).nseries(x, n=5) == x - 2*x**3/9 + 17*x**5/450 + O(x**6)
    assert Si(x).nseries(x, 1, n=3) == Si(1) + x*sin(1) + x**2*(-sin(1)/2 + cos(1)/2) + O(x**3)

def test_ci():
    m1 = exp_polar(I*pi)
    m1_ = exp_polar(-I*pi)
    pI = exp_polar(I*pi/2)
    mI = exp_polar(-I*pi/2)

    assert Ci(m1*x) == Ci(x) + I*pi
    assert Ci(m1_*x) == Ci(x) - I*pi
    assert Ci(pI*x) == Chi(x) + I*pi/2
    assert Ci(mI*x) == Chi(x) - I*pi/2
    assert Chi(m1*x) == Chi(x) + I*pi
    assert Chi(m1_*x) == Chi(x) - I*pi
    assert Chi(pI*x) == Ci(x) + I*pi/2
    assert Chi(mI*x) == Ci(x) - I*pi/2
    assert Ci(exp_polar(2*I*pi)*x) == Ci(x) + 2*I*pi
    assert Chi(exp_polar(-2*I*pi)*x) == Chi(x) - 2*I*pi
    assert Chi(exp_polar(2*I*pi)*x) == Chi(x) + 2*I*pi
    assert Ci(exp_polar(-2*I*pi)*x) == Ci(x) - 2*I*pi

    assert mytd(Ci(x), cos(x)/x, x)
    assert mytd(Chi(x), cosh(x)/x, x)

    assert mytn(Ci(x), Ci(x).rewrite(Ei),
                Ei(x*exp_polar(-I*pi/2))/2 + Ei(x*exp_polar(I*pi/2))/2, x)
    assert mytn(Chi(x), Chi(x).rewrite(Ei),
                Ei(x)/2 + Ei(x*exp_polar(I*pi))/2 - I*pi/2, x)

    assert tn_arg(Ci)
    assert tn_arg(Chi)

    from sympy import O, EulerGamma, log, limit
    assert Ci(x).nseries(x, n=4) == EulerGamma + log(x) - x**2/4 + x**4/96 + O(x**5)
    assert Chi(x).nseries(x, n=4) == EulerGamma + log(x) + x**2/4 + x**4/96 + O(x**5)
    assert limit(log(x) - Ci(2*x), x, 0) == -log(2) - EulerGamma

def test_fresnel():
    assert fresnels(0) == 0
    assert fresnels(oo) == S.Half
    assert fresnels(-oo) == -S.Half

    assert fresnels(z) == fresnels(z)
    assert fresnels(-z) == -fresnels(z)
    assert fresnels(I*z) == -I*fresnels(z)
    assert fresnels(-I*z) == I*fresnels(z)

    assert conjugate(fresnels(z)) == fresnels(conjugate(z))

    assert fresnels(z).diff(z) == sin(pi*z**2/2)

    assert fresnels(z).rewrite(erf) == (S.One+I)/4 * (erf((S.One+I)/2*sqrt(pi)*z) - I*erf((S.One-I)/2*sqrt(pi)*z))

    assert fresnels(z).rewrite(hyper) == pi*z**3/6 * hyper([S(3)/4], [S(3)/2, S(7)/4], -pi**2*z**4/16)

    assert fresnels(z).series(z, n=15) == pi*z**3/6 - pi**3*z**7/336 + pi**5*z**11/42240 + O(z**15)

    assert fresnels(w).is_real is True

    assert fresnels(z).as_real_imag() == ((fresnels(re(z) - I*re(z)*Abs(im(z))/Abs(re(z)))/2 + fresnels(re(z) + I*re(z)*Abs(im(z))/Abs(re(z)))/2,
                                          I*(fresnels(re(z) - I*re(z)*Abs(im(z))/Abs(re(z))) - fresnels(re(z) + I*re(z)*Abs(im(z))/Abs(re(z))))*
                                          re(z)*Abs(im(z))/(2*im(z)*Abs(re(z)))))

    assert fresnels(2+3*I).as_real_imag() == (fresnels(2 + 3*I)/2 + fresnels(2 - 3*I)/2, I*(fresnels(2 - 3*I) - fresnels(2 + 3*I))/2)

    assert expand_func(integrate(fresnels(z), z)) == z*fresnels(z) + cos(pi*z**2/2)/pi

    assert fresnels(z).rewrite(meijerg) == sqrt(2)*pi*z**(S(9)/4)* \
           meijerg(((), (1,)), ((S(3)/4,), (S(1)/4, 0)), -pi**2*z**4/16)/(2*(-z)**(S(3)/4)*(z**2)**(S(3)/4))

    assert fresnelc(0) == 0
    assert fresnelc(oo) == S.Half
    assert fresnelc(-oo) == -S.Half

    assert fresnelc(z) == fresnelc(z)
    assert fresnelc(-z) == -fresnelc(z)
    assert fresnelc(I*z) == I*fresnelc(z)
    assert fresnelc(-I*z) == -I*fresnelc(z)

    assert conjugate(fresnelc(z)) == fresnelc(conjugate(z))

    assert fresnelc(z).diff(z) == cos(pi*z**2/2)

    assert fresnelc(z).rewrite(erf) == (S.One - I)/4 * (erf((S.One + I)/2*sqrt(pi)*z) + I*erf((S.One - I)/2*sqrt(pi)*z))

    assert fresnelc(z).rewrite(hyper) == z * hyper([S.One/4], [S.One/2, S(5)/4], -pi**2*z**4/16)

    assert fresnelc(z).series(z, n=15) == z - pi**2*z**5/40 + pi**4*z**9/3456 - pi**6*z**13/599040 + O(z**15)

    assert fresnelc(w).is_real is True

    assert fresnelc(z).as_real_imag() == ((fresnelc(re(z) - I*re(z)*Abs(im(z))/Abs(re(z)))/2 + fresnelc(re(z) + I*re(z)*Abs(im(z))/Abs(re(z)))/2,
                                           I*(fresnelc(re(z) - I*re(z)*Abs(im(z))/Abs(re(z))) - fresnelc(re(z) + I*re(z)*Abs(im(z))/Abs(re(z))))*
                                           re(z)*Abs(im(z))/(2*im(z)*Abs(re(z)))))

    assert fresnelc(2+3*I).as_real_imag() == (fresnelc(2 - 3*I)/2 + fresnelc(2 + 3*I)/2, I*(fresnelc(2 - 3*I) - fresnelc(2 + 3*I))/2)

    assert expand_func(integrate(fresnelc(z), z)) == z*fresnelc(z) - sin(pi*z**2/2)/pi

    assert fresnelc(z).rewrite(meijerg) == sqrt(2)*pi*z**(S(3)/4)* \
           meijerg(((), (1,)), ((S(1)/4,), (S(3)/4, 0)), -pi**2*z**4/16)/(2*(-z)**(S(1)/4)*(z**2)**(S(1)/4))


    from sympy.utilities.randtest import test_numerically

    test_numerically(re(fresnels(z)), fresnels(z).as_real_imag()[0], z)
    test_numerically(im(fresnels(z)), fresnels(z).as_real_imag()[1], z)
    test_numerically(fresnels(z), fresnels(z).rewrite(hyper), z)
    test_numerically(fresnels(z), fresnels(z).rewrite(meijerg), z)

    test_numerically(re(fresnelc(z)), fresnelc(z).as_real_imag()[0], z)
    test_numerically(im(fresnelc(z)), fresnelc(z).as_real_imag()[1], z)
    test_numerically(fresnelc(z), fresnelc(z).rewrite(hyper), z)
    test_numerically(fresnelc(z), fresnelc(z).rewrite(meijerg), z)
