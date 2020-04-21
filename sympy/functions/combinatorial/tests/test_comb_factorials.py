from sympy.concrete.products import Product
from sympy.core.expr import unchanged
from sympy.core.function import ArgumentIndexError, expand_func
from sympy.core.mod import Mod
from sympy.core.mul import Mul
from sympy.core.numbers import (EulerGamma, Float, I, Rational, nan,
    oo, pi, zoo)
EulerGamma = EulerGamma()  # no EulerGamma = S.EulerGamma in numbers
from sympy.core.relational import Eq
from sympy.core.singleton import S
from sympy.core.symbol import (Dummy, Symbol, symbols)
from sympy.functions.combinatorial.factorials import (
    ff, rf, binomial, factorial, factorial2,
    )
from sympy.functions.combinatorial.factorials import subfactorial
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.special.gamma_functions import (
    gamma, loggamma, polygamma, uppergamma)
from sympy.polys.polytools import Poly
from sympy.series.order import O
from sympy.sets.fancysets import Range
from sympy.simplify.simplify import simplify
from sympy.testing.pytest import XFAIL, raises, slow
from sympy.utilities.iterables import cartes

#Solves and Fixes Issue #10388 - This is the updated test for the same solved issue

def test_rf_eval_apply():
    x, y = symbols('x,y')
    n, k = symbols('n k', integer=True)
    m = Symbol('m', integer=True, nonnegative=True)

    assert rf(nan, y) is nan
    assert rf(x, nan) is nan

    assert unchanged(rf, x, y)

    assert rf(oo, 0) == 1
    assert rf(-oo, 0) == 1

    assert rf(oo, 6) is oo
    assert rf(-oo, 7) is -oo
    assert rf(-oo, 6) is oo

    assert rf(oo, -6) is oo
    assert rf(-oo, -7) is oo

    assert rf(-1, pi) == 0
    assert rf(-5, 1 + I) == 0

    assert unchanged(rf, -3, k)
    assert unchanged(rf, x, Symbol('k', integer=False))
    assert rf(-3, Symbol('k', integer=False)) == 0
    assert rf(Symbol('x', negative=True, integer=True), Symbol('k', integer=False)) == 0

    assert rf(x, 0) == 1
    assert rf(x, 1) == x
    assert rf(x, 2) == x*(x + 1)
    assert rf(x, 3) == x*(x + 1)*(x + 2)
    assert rf(x, 5) == x*(x + 1)*(x + 2)*(x + 3)*(x + 4)

    assert rf(x, -1) == 1/(x - 1)
    assert rf(x, -2) == 1/((x - 1)*(x - 2))
    assert rf(x, -3) == 1/((x - 1)*(x - 2)*(x - 3))

    assert rf(1, 100) == factorial(100)

    assert rf(x**2 + 3*x, 2) == (x**2 + 3*x)*(x**2 + 3*x + 1)
    assert isinstance(rf(x**2 + 3*x, 2), Mul)
    assert rf(x**3 + x, -2) == 1/((x**3 + x - 1)*(x**3 + x - 2))

    assert rf(Poly(x**2 + 3*x, x), 2) == Poly(x**4 + 8*x**3 + 19*x**2 + 12*x, x)
    assert isinstance(rf(Poly(x**2 + 3*x, x), 2), Poly)
    raises(ValueError, lambda: rf(Poly(x**2 + 3*x, x, y), 2))
    assert rf(Poly(x**3 + x, x), -2) == 1/(x**6 - 9*x**5 + 35*x**4 - 75*x**3 + 94*x**2 - 66*x + 20)
    raises(ValueError, lambda: rf(Poly(x**3 + x, x, y), -2))

    assert rf(x, m).is_integer is None
    assert rf(n, k).is_integer is None
    assert rf(n, m).is_integer is True
    assert rf(n, k + pi).is_integer is False
    assert rf(n, m + pi).is_integer is False
    assert rf(pi, m).is_integer is False

    assert rf(x, k).rewrite(ff) == ff(x + k - 1, k)
    assert rf(x, k).rewrite(binomial) == factorial(k)*binomial(x + k - 1, k)
    assert rf(n, k).rewrite(factorial) == \
        factorial(n + k - 1) / factorial(n - 1)
    assert rf(x, y).rewrite(factorial) == rf(x, y)
    assert rf(x, y).rewrite(binomial) == rf(x, y)
    assert rf(x, y).rewrite(gamma) == gamma(x + y)/gamma(x)

    import random
    from mpmath import rf as mpmath_rf
    for i in range(100):
        x = -500 + 500 * random.random()
        k = -500 + 500 * random.random()
        assert (abs(mpmath_rf(x, k) - rf(x, k)) < 10**(-15))


def test_ff_eval_apply():
    x, y = symbols('x,y')
    n, k = symbols('n k', integer=True)
    m = Symbol('m', integer=True, nonnegative=True)

    assert ff(nan, y) is nan
    assert ff(x, nan) is nan

    assert unchanged(ff, x, y)

    assert ff(oo, 0) == 1
    assert ff(-oo, 0) == 1

    assert ff(oo, 6) is oo
    assert ff(-oo, 7) is -oo
    assert ff(-oo, 6) is oo

    assert ff(oo, -6) is oo
    assert ff(-oo, -7) is oo

    assert ff(x, 0) == 1
    assert ff(x, 1) == x
    assert ff(x, 2) == x*(x - 1)
    assert ff(x, 3) == x*(x - 1)*(x - 2)
    assert ff(x, 5) == x*(x - 1)*(x - 2)*(x - 3)*(x - 4)

    assert ff(x, -1) == 1/(x + 1)
    assert ff(x, -2) == 1/((x + 1)*(x + 2))
    assert ff(x, -3) == 1/((x + 1)*(x + 2)*(x + 3))

    assert ff(100, 100) == factorial(100)

    assert ff(2*x**2 - 5*x, 2) == (2*x**2  - 5*x)*(2*x**2 - 5*x - 1)
    assert isinstance(ff(2*x**2 - 5*x, 2), Mul)
    assert ff(x**2 + 3*x, -2) == 1/((x**2 + 3*x + 1)*(x**2 + 3*x + 2))

    assert ff(Poly(2*x**2 - 5*x, x), 2) == Poly(4*x**4 - 28*x**3 + 59*x**2 - 35*x, x)
    assert isinstance(ff(Poly(2*x**2 - 5*x, x), 2), Poly)
    raises(ValueError, lambda: ff(Poly(2*x**2 - 5*x, x, y), 2))
    assert ff(Poly(x**2 + 3*x, x), -2) == 1/(x**4 + 12*x**3 + 49*x**2 + 78*x + 40)
    raises(ValueError, lambda: ff(Poly(x**2 + 3*x, x, y), -2))


    assert ff(x, m).is_integer is None
    assert ff(n, k).is_integer is None
    assert ff(n, m).is_integer is True
    assert ff(n, k + pi).is_integer is False
    assert ff(n, m + pi).is_integer is False
    assert ff(pi, m).is_integer is False

    assert isinstance(ff(x, x), ff)
    assert ff(n, n) == factorial(n)

    assert ff(x, k).rewrite(rf) == rf(x - k + 1, k)
    assert ff(x, k).rewrite(gamma) == (-1)**k*gamma(k - x) / gamma(-x)
    assert ff(n, k).rewrite(factorial) == factorial(n) / factorial(n - k)
    assert ff(x, k).rewrite(binomial) == factorial(k) * binomial(x, k)
    assert ff(x, y).rewrite(factorial) == ff(x, y)
    assert ff(x, y).rewrite(binomial) == ff(x, y)

    import random
    from mpmath import ff as mpmath_ff
    for i in range(100):
        x = -500 + 500 * random.random()
        k = -500 + 500 * random.random()
        a = mpmath_ff(x, k)
        b = ff(x, k)
        assert (abs(a - b) < abs(a) * 10**(-15))


def test_rf_ff_eval_hiprec():
    maple = Float('6.9109401292234329956525265438452')
    us = ff(18, Rational(2, 3)).evalf(32)
    assert abs(us - maple)/us < 1e-31

    maple = Float('6.8261540131125511557924466355367')
    us = rf(18, Rational(2, 3)).evalf(32)
    assert abs(us - maple)/us < 1e-31

    maple = Float('34.007346127440197150854651814225')
    us = rf(Float('4.4', 32), Float('2.2', 32));
    assert abs(us - maple)/us < 1e-31


def test_rf_lambdify_mpmath():
    from sympy import lambdify
    x, y = symbols('x,y')
    f = lambdify((x,y), rf(x, y), 'mpmath')
    maple = Float('34.007346127440197')
    us = f(4.4, 2.2)
    assert abs(us - maple)/us < 1e-15


def test_factorial():
    i = Symbol('i')
    x = Symbol('x')
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True, nonnegative=True)
    r = Symbol('r', integer=False)
    s = Symbol('s', integer=False, negative=True)
    t = Symbol('t', nonnegative=True)
    u = Symbol('u', noninteger=True)

    assert factorial(-2) is zoo
    assert factorial(0) == 1
    assert factorial(7) == 5040
    assert factorial(19) == 121645100408832000
    assert factorial(31) == 8222838654177922817725562880000000
    assert factorial(n).func == factorial
    assert factorial(2*n).func == factorial

    f3 = factorial(3., evaluate=False)
    assert f3.rewrite(Product) == Product(i, (i, 1, 3.0))
    assert f3.is_integer is False
    assert f3.is_positive
    assert f3.is_even is False
    assert f3.is_composite is False

    assert factorial(x).is_integer is None
    assert factorial(n).is_integer is None
    assert factorial(k).is_integer
    assert factorial(r).is_integer is None

    assert factorial(n).is_positive is None
    assert factorial(k).is_positive

    assert factorial(x).is_real is None
    assert factorial(n).is_real is None
    assert factorial(k).is_real is True
    assert factorial(r).is_real is None
    assert factorial(s).is_real is True
    assert factorial(t).is_real is True
    assert factorial(u).is_real is None  # e.g. I!

    assert factorial(x).is_composite is None
    assert factorial(n).is_composite is None
    assert factorial(k).is_composite is None
    assert factorial(k + 3).is_composite is True
    assert factorial(r).is_composite is None
    assert factorial(s).is_composite is None
    assert factorial(t).is_composite is None
    assert factorial(u).is_composite is None

    assert factorial(oo) is oo


def test_factorial_Mod():
    pr = Symbol('pr', prime=True)
    p, q = 10**9 + 9, 10**9 + 33 # prime modulo
    r, s = 10**7 + 5, 33333333 # composite modulo
    assert Mod(factorial(pr - 1), pr) == pr - 1
    assert Mod(factorial(pr - 1), -pr) == -1
    assert Mod(factorial(r - 1, evaluate=False), r) == 0
    assert Mod(factorial(s - 1, evaluate=False), s) == 0
    assert Mod(factorial(p - 1, evaluate=False), p) == p - 1
    assert Mod(factorial(q - 1, evaluate=False), q) == q - 1
    assert Mod(factorial(p - 50, evaluate=False), p) == 854928834
    assert Mod(factorial(q - 1800, evaluate=False), q) == 905504050
    assert Mod(factorial(153, evaluate=False), r) == Mod(factorial(153), r)
    assert Mod(factorial(255, evaluate=False), s) == Mod(factorial(255), s)
    assert Mod(factorial(4, evaluate=False), 3) == S.Zero
    assert Mod(factorial(5, evaluate=False), 6) == S.Zero
    assert factorial(6, evaluate=False) % 11 == 5
    assert factorial(5, evaluate=False) % 15 == 0


def test_factorial_diff():
    n = Symbol('n', integer=True)

    assert factorial(n).diff(n) == \
        gamma(1 + n)*polygamma(0, 1 + n)
    assert factorial(n**2).diff(n) == \
        2*n*gamma(1 + n**2)*polygamma(0, 1 + n**2)
    raises(ArgumentIndexError, lambda: factorial(n**2).fdiff(2))


def test_factorial_series():
    n = Symbol('n', integer=True)

    assert factorial(n).series(n, 0, 3) == \
        1 - n*EulerGamma + n**2*(EulerGamma**2/2 + pi**2/12) + O(n**3)


def test_factorial_rewrite():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True, nonnegative=True)

    assert factorial(n).rewrite(gamma) == gamma(n + 1)
    _i = Dummy('i')
    assert factorial(k).rewrite(Product).dummy_eq(Product(_i, (_i, 1, k)))
    assert factorial(n).rewrite(Product) == factorial(n)


def test_factorial2():
    n = Symbol('n', integer=True)

    assert factorial2(-1) == 1
    assert factorial2(0) == 1
    assert factorial2(7) == 105
    assert factorial2(8) == 384

    # The following is exhaustive
    tt = Symbol('tt', integer=True, nonnegative=True)
    tte = Symbol('tte', even=True, nonnegative=True)
    tpe = Symbol('tpe', even=True, positive=True)
    tto = Symbol('tto', odd=True, nonnegative=True)
    tf = Symbol('tf', integer=True, nonnegative=False)
    tfe = Symbol('tfe', even=True, nonnegative=False)
    tfo = Symbol('tfo', odd=True, nonnegative=False)
    ft = Symbol('ft', integer=False, nonnegative=True)
    ff = Symbol('ff', integer=False, nonnegative=False)
    fn = Symbol('fn', integer=False)
    nt = Symbol('nt', nonnegative=True)
    nf = Symbol('nf', nonnegative=False)
    nn = Symbol('nn')
    z = Symbol('z', zero=True)
    #Solves and Fixes Issue #10388 - This is the updated test for the same solved issue
    raises(ValueError, lambda: factorial2(oo))
    raises(ValueError, lambda: factorial2(Rational(5, 2)))
    raises(ValueError, lambda: factorial2(-4))
    assert factorial2(n).is_integer is None
    assert factorial2(tt - 1).is_integer
    assert factorial2(tte - 1).is_integer
    assert factorial2(tpe - 3).is_integer
    assert factorial2(tto - 4).is_integer
    assert factorial2(tto - 2).is_integer
    assert factorial2(tf).is_integer is None
    assert factorial2(tfe).is_integer is None
    assert factorial2(tfo).is_integer is None
    assert factorial2(ft).is_integer is None
    assert factorial2(ff).is_integer is None
    assert factorial2(fn).is_integer is None
    assert factorial2(nt).is_integer is None
    assert factorial2(nf).is_integer is None
    assert factorial2(nn).is_integer is None

    assert factorial2(n).is_positive is None
    assert factorial2(tt - 1).is_positive is True
    assert factorial2(tte - 1).is_positive is True
    assert factorial2(tpe - 3).is_positive is True
    assert factorial2(tpe - 1).is_positive is True
    assert factorial2(tto - 2).is_positive is True
    assert factorial2(tto - 1).is_positive is True
    assert factorial2(tf).is_positive is None
    assert factorial2(tfe).is_positive is None
    assert factorial2(tfo).is_positive is None
    assert factorial2(ft).is_positive is None
    assert factorial2(ff).is_positive is None
    assert factorial2(fn).is_positive is None
    assert factorial2(nt).is_positive is None
    assert factorial2(nf).is_positive is None
    assert factorial2(nn).is_positive is None

    assert factorial2(tt).is_even is None
    assert factorial2(tt).is_odd is None
    assert factorial2(tte).is_even is None
    assert factorial2(tte).is_odd is None
    assert factorial2(tte + 2).is_even is True
    assert factorial2(tpe).is_even is True
    assert factorial2(tpe).is_odd is False
    assert factorial2(tto).is_odd is True
    assert factorial2(tf).is_even is None
    assert factorial2(tf).is_odd is None
    assert factorial2(tfe).is_even is None
    assert factorial2(tfe).is_odd is None
    assert factorial2(tfo).is_even is False
    assert factorial2(tfo).is_odd is None
    assert factorial2(z).is_even is False
    assert factorial2(z).is_odd is True


def test_factorial2_rewrite():
    n = Symbol('n', integer=True)
    assert factorial2(n).rewrite(gamma) == \
        2**(n/2)*Piecewise((1, Eq(Mod(n, 2), 0)), (sqrt(2)/sqrt(pi), Eq(Mod(n, 2), 1)))*gamma(n/2 + 1)
    assert factorial2(2*n).rewrite(gamma) == 2**n*gamma(n + 1)
    assert factorial2(2*n + 1).rewrite(gamma) == \
        sqrt(2)*2**(n + S.Half)*gamma(n + Rational(3, 2))/sqrt(pi)

@XFAIL
def test_binomial_fail1():
    i, j = symbols('i, j', integer=True)
    # doesn't hold for (i, j) == (0, 0)
    assert binomial(i, j).equals(
        binomial(i - 1, j - 1) + binomial(i - 1, j)) is not True


@XFAIL
def test_binomial_fail2():
    i, j = symbols('i, j', integer=True)
    # gives None but should be True
    assert binomial(2*j, j).equals(
        (-4)**j*binomial(-S.Half, j)) is True


@slow
@XFAIL
def test_binomial_fail3():
    i, j = symbols('i j', integer=True, zero=False)
    # gives None but should be Tru
    assert binomial(i, j).equals(
        binomial(i - 1, j - 1) + binomial(i - 1, j))


def test_binomial():
    x = Symbol('x')

    n, k = symbols('n k', integer=True)
    kp = Symbol('kp', integer=True, positive=True)
    kn = Symbol('kn', integer=True, negative=True)
    kne = Symbol('kn', integer=True, negative=True, even=True)
    nt, kt = symbols('nt kt', integer=False)
    ntr = symbols('ntr', rational=False, real=True)
    a, b = symbols('a b', integer=True, nonnegative=True)
    u = Symbol('u', negative=True)
    v = Symbol('v', nonnegative=True)
    p = Symbol('p', positive=True)
    z = Symbol('z', zero=True)
    nz = Symbol('z', zero=False)

    i, j = symbols('i, j', finite=True)
    assert binomial(i, i - j).equals(binomial(i, j))
    # Consecutive Neighbours
    assert binomial(i, j).equals((i/(i - j))* binomial(i - 1, j))
    assert binomial(i, j).equals(((i - j + 1)/(i + 1))* binomial(i + 1, j))
    assert binomial(i, j).equals(((i - j + 1)/j)* binomial(i, j - 1))
    # Relations between Contiguous functions
    assert (binomial(i, j) + binomial(i, j + 1)).equals(binomial(i + 1, j + 1))

    # if not (i == j == 0) check Identity
    i, j = symbols('i', positive=True, integer=True), symbols('i', negative=True, integer=True)
    assert binomial(i, j).equals(binomial(i - 1, j - 1) + binomial(i - 1, j))
    i, j = j, i
    assert binomial(i, j).equals(binomial(i - 1, j - 1) + binomial(i - 1, j))

    i, j = symbols('i, j', positive=True, integer=True)

    # issues 14529 and 14625
    assert binomial(x, z) == binomial(x, x) == 1
    assert binomial(x, 1) == binomial(x, x - 1) == x
    assert binomial(x, 1.) == binomial(x, x - 1.) == x
    assert binomial(x + 1, x) == x + 1

    assert unchanged(binomial, x, -1)
    assert binomial(-1 + nz, -1) == 0
    assert binomial(v, v + 1) == 0
    assert binomial(nt, nt + 1) == 0
    assert unchanged(binomial, x, x + 1)  # Piecewise((0, Ne(n, -1)), (1, True))

    assert unchanged(binomial, kp, -kn)
    assert binomial(kp, 2*kp) == 0
    assert unchanged(binomial, n, u)
    assert unchanged(binomial, kp, u)
    assert unchanged(binomial, n, p)
    assert unchanged(binomial, n, k)
    assert unchanged(binomial, n, n + p)

    assert unchanged(binomial, n, 3)
    assert binomial(n, 3).expand(func=True) ==  n**3/6 - n**2/2 + n/3
    assert expand_func(binomial(n, 3)) ==  n*(n - 2)*(n - 1)/6
    assert expand_func(binomial(n, 2)) == n*(n - 1)/2
    assert expand_func(binomial(n, n - 2)) == n*(n - 1)/2
    assert expand_func(binomial(n, n - 3)) == n*(n - 2)*(n - 1)/6
    #issue #18802
    assert expand_func(binomial(x + 1, x - 1)) == x*(x + 1)/2
    assert expand_func(binomial(x**2 + 1, x**2)) == x**2 + 1

    assert binomial(n, k).is_integer
    assert binomial(nt, k).is_integer is None
    assert binomial(k, nt).is_integer is None
    assert binomial(k, ntr).is_integer is False
    assert binomial(x, nt).is_integer is None  # see next,
    assert binomial(S.Half, S.Half*3) == 0     #<--------+
    assert binomial(n, a).is_integer
    assert binomial(x, x + kp).is_integer
    assert not binomial(0, sqrt(3) + sqrt(11)*I).is_integer
    assert binomial(3, oo, evaluate=False).is_integer

    assert binomial(-oo, oo) == 0
    assert binomial(-1, oo) is S.ComplexInfinity
    assert binomial(oo, oo) is S.NaN
    assert binomial(-oo, oo, evaluate=False).is_integer
    assert binomial(-1, oo, evaluate=False).is_integer is False
    assert binomial(oo, oo, evaluate=False).is_integer is False
    assert binomial(1, 2) == 0
    assert binomial(-1, 2) == 1
    assert binomial(1, -1) == 0
    assert binomial(-1, 1) == -1
    assert binomial(-10, 1) == -10
    assert binomial(-10, 7) == -11440
    assert binomial(-S.Half, S(5)/2) == 0
    assert binomial(-S.Half, -3.) == 0

    assert binomial(gamma(25), 6) == 79232165267303928292058750056084441948572511312165380965440075720159859792344339983120618959044048198214221915637090855535036339620413440000
    assert binomial(1324, 47) == 906266255662694632984994480774946083064699457235920708992926525848438478406790323869952
    assert binomial(1735, 43) == 190910140420204130794758005450919715396159959034348676124678207874195064798202216379800
    assert binomial(2512, 53) == 213894469313832631145798303740098720367984955243020898718979538096223399813295457822575338958939834177325304000
    assert binomial(3383, 52) == 27922807788818096863529701501764372757272890613101645521813434902890007725667814813832027795881839396839287659777235
    assert binomial(4321, 51) == 124595639629264868916081001263541480185227731958274383287107643816863897851139048158022599533438936036467601690983780576

    assert binomial(a, b).is_nonnegative is True
    assert binomial(kne, kne - 1).is_zero is False
    assert binomial(kne, kne - 2).is_nonnegative is True
    assert binomial(kne, kne - 3).is_nonnegative is False
    assert binomial(-1, 2, evaluate=False).is_nonnegative is True
    assert binomial(10, 5, evaluate=False).is_nonnegative is True
    assert binomial(10, -3, evaluate=False).is_nonnegative is True
    assert binomial(-10, -3, evaluate=False).is_nonnegative is True
    assert binomial(-10, 2, evaluate=False).is_nonnegative is True
    assert binomial(-10, 1, evaluate=False).is_nonnegative is False
    assert binomial(-10, 7, evaluate=False).is_nonnegative is False

    # issue #13980 and #13981
    assert binomial(-7, -5) == 0
    assert binomial(-23, -12) == 0
    assert binomial(Rational(13, 2), -10) == 0
    assert binomial(-49, -51) == 1225

    assert binomial(19, Rational(-7, 2)) == S(-68719476736)/(911337863661225*pi)
    assert binomial(0, Rational(3, 2)) == S(-2)/(3*pi)
    assert binomial(-3, Rational(-7, 2)) is zoo
    assert binomial(kn, kt) is zoo

    assert binomial(nt, kt).func == binomial
    assert unchanged(binomial, nt, Rational(15, 6))
    assert binomial(kn - 3, kn - 1, evaluate=False).rewrite(gamma) == 0
    assert binomial(nt, Rational(15, 6)).rewrite(gamma
        ) == 8*gamma(nt + 1)/(15*sqrt(pi)*gamma(nt - Rational(3, 2)))
    assert binomial(Rational(20, 3), Rational(-10, 8)
        ).rewrite(gamma) == gamma(Rational(23, 3)
        )/(gamma(Rational(-1, 4))*gamma(Rational(107, 12)))
    assert binomial(Rational(19, 2), Rational(-7, 2)
        ).rewrite(gamma) == Rational(-1615, 8388608)
    assert binomial(Rational(-13, 5), Rational(-7, 8)
        ).rewrite(gamma) == gamma(Rational(-8, 5))/(gamma(Rational(-29, 40))*gamma(Rational(1,8)))
    assert binomial(Rational(-19, 8), Rational(-13, 5)
        ).rewrite(gamma) == gamma(Rational(-11, 8)
        )/(gamma(Rational(-8, 5))*gamma(Rational(49, 40)))

    # binomial for complexes
    assert binomial(I, Rational(-89, 8)).rewrite(gamma
        ) == gamma(1 + I)/(gamma(Rational(-81, 8)
        )*gamma(Rational(97, 8) + I))
    assert binomial(I, 2*I).rewrite(gamma
        ) == gamma(1 + I)/(gamma(1 - I)*gamma(1 + 2*I))
    assert binomial(-7, I) is zoo
    assert binomial(Rational(-7, 6), I).rewrite(gamma
        ) == gamma(Rational(-1, 6))/(gamma(Rational(-1, 6) - I)*gamma(1 + I))
    assert binomial(1 + 2*I, 1 + 3*I).rewrite(gamma
        ) == gamma(2 + 2*I)/(gamma(1 - I)*gamma(2 + 3*I))
    assert binomial(I, 5) == Rational(1, 3) - I/S(12)
    assert binomial(2*I + 3, 7) == -13*I/S(63)
    assert isinstance(binomial(I, n), binomial)
    assert expand_func(binomial(3, 2, evaluate=False)) == 3
    assert expand_func(binomial(n, 0, evaluate=False)) == 1
    assert expand_func(binomial(n, -2, evaluate=False)) == binomial(n, -2)
    assert expand_func(binomial(n, k)) == binomial(n, k)
    # expansion of numerical results is automatic
    assert binomial(3 + sqrt(11)*I, 3) == -10


def test_binomial_Mod():
    # positive k

    p, q = 10**5 + 3, 10**9 + 33 # prime modulo
    r = 10**7 + 5 # composite modulo

    # Lucas Theorem
    assert Mod(binomial(156675, 4433, evaluate=False), p
        ) == Mod(binomial(156675, 4433), p)

    # factorial Mod
    assert Mod(binomial(1234, 432, evaluate=False), q
        ) == Mod(binomial(1234, 432), q)

    # binomial factorize
    assert Mod(binomial(253, 113, evaluate=False), r
        ) == Mod(binomial(253, 113), r)

    # postive n
    lucas = 0
    fac = 0
    binfac = 0
    for k in range(1, 5):
        for n in range(k, k + 4):
            for i in Range(2, n*2 + 1):
                assert binomial(n, k, evaluate=False)%i == \
                    binomial(n, k)%i, (k, n, i)
                if i.is_prime:
                    if i <= n:
                        lucas += 1
                    else:
                        fac += 1
                else:
                    binfac += 1
    assert lucas
    assert fac
    assert binfac

    # negative n
    lucas = 0
    fac = 0
    binfac = 0
    for k in range(1, 5):
        for n in range(-k - 4, 1):
            for i in Range(2, -n*2 + 1):
                assert binomial(n, k, evaluate=False)%i == \
                    binomial(n, k)%i, (n, k, i,
                    binomial(n, k, evaluate=False)%i,
                    binomial(n, k)%i)
                if i.is_prime:
                    if i <= abs(n):
                        lucas += 1
                    else:
                        fac += 1
                else:
                    binfac += 1

    assert lucas
    assert fac
    assert binfac

    x = Symbol('x', zero=False)
    assert binomial(3, 4, evaluate=False) % x == 0
    assert binomial(3, 3, evaluate=False) % x == 1 % x
    assert binomial(-1, 3, evaluate=False) % x == -1 % x
    assert binomial(3, 1, evaluate=False) % x == 3 % x
    assert binomial(-3, 1, evaluate=False) % x == -3 % x
    assert binomial(0, -5, evaluate=False) % x == 0
    assert binomial(-3, -2, evaluate=False) % x == 0
    a = binomial(-3, -4, evaluate=False) % 3
    b = binomial(-3, 1, evaluate=False) % 3
    assert a == b

    # check int-like handling
    assert binomial(5, 2, evaluate=False) % 7 == 3
    a = binomial(5., 2, evaluate=False) % 7
    assert a == 3 and a.is_Float
    a = binomial(5, 2., evaluate=False) % 7
    assert a == 3 and a.is_Float
    a = binomial(5, 2, evaluate=False) % 7.
    assert a == 3 and a.is_Float
    assert (binomial(5, 5., evaluate=False) % 7).is_Float
    assert (binomial(5., 5, evaluate=False) % 7).is_Float
    assert (binomial(5, 5, evaluate=False) % 7.).is_Float
    assert (binomial(5, 1., evaluate=False) % 7).is_Float
    assert (binomial(5., 1, evaluate=False) % 7).is_Float
    assert (binomial(5, 1, evaluate=False) % 7.).is_Float

    raises(ValueError, lambda: binomial(pi, 2, evaluate=False) % 3)
    raises(ValueError, lambda: binomial(4, pi, evaluate=False) % 3)
    raises(ValueError, lambda: binomial(3, 2, evaluate=False) % pi)

    n = symbols('n', integer=True)
    assert unchanged(Mod, binomial(n, 2), 7)


def test_binomial_diff():
    n, k = symbols('n k', integer=True)

    assert binomial(n, k).diff(n) == \
        (-polygamma(0, 1 + n - k) + polygamma(0, 1 + n))*binomial(n, k)
    assert binomial(n**2, k**3).diff(n) == \
        2*n*(-polygamma(
            0, 1 + n**2 - k**3) + polygamma(0, 1 + n**2))*binomial(n**2, k**3)

    assert binomial(n, k).diff(k) == \
        (-polygamma(0, 1 + k) + polygamma(0, 1 + n - k))*binomial(n, k)
    assert binomial(n**2, k**3).diff(k) == \
        3*k**2*(-polygamma(
            0, 1 + k**3) + polygamma(0, 1 + n**2 - k**3))*binomial(n**2, k**3)
    raises(ArgumentIndexError, lambda: binomial(n, k).fdiff(3))


def test_simplify():
    n, k = symbols('n k', integer=True)
    assert factorial(factorial(n)/(factorial(k)*factorial(
        -k + n))).simplify() == gamma(1 + gamma(n + 1)/(
            gamma(k + 1)*gamma(-k + n + 1)))
    assert factorial(n**3 + k).as_leading_term(n) == factorial(k)


def test_binomial_rewrite():
    n, k = symbols('n k', integer=True, nonnegative=True)
    x, y = symbols('x y')
    i, i_0 = symbols('i i_0')

    assert binomial(n, k).rewrite(factorial) == factorial(n)/(
        factorial(k)*factorial(n - k))
    assert binomial(n, k).rewrite(gamma) == gamma(
        n + 1)/(gamma(k + 1)*gamma(n - k + 1))
    assert binomial(n, k).rewrite(ff) == ff(n, k) / factorial(k)
    assert binomial(n, x).rewrite(ff) == binomial(n, x)
    assert binomial(n, k).rewrite('tractable') == exp(
        -loggamma(k + 1))*exp(loggamma(n + 1))*exp(-loggamma(-k + n + 1))
    assert binomial(n, k).rewrite(factorial) == factorial(n)/(
        factorial(k)*factorial(n - k))
    assert binomial(n, k).rewrite(gamma) == gamma(
        n + 1)/(gamma(k + 1)*gamma(n - k + 1))
    assert binomial(n, k).rewrite(ff) == ff(n, k) / factorial(k)
    assert binomial(n, x).rewrite(ff) == binomial(n, x)
    assert binomial(n, k).rewrite('tractable') == exp(
        -loggamma(k + 1))*exp(loggamma(n + 1))*exp(-loggamma(-k + n + 1))

    y, i = i, i_0  # tests igen, too
    h = S.Half
    v = [-1, 0, zoo, -h, 1, 0, zoo, 1, 1, 1, zoo, 0, 0, 0, 2/pi,
        h, 1, 1, 2/pi, 1, 0, h, -2/(3*pi), 0, 1, 1, 4/pi, S(3)/2, 1,
        S(3)/2, 4/(3*pi), 1, 0, S(3)/8, -4/(15*pi), 0]
    ix = 0
    for i in range(-1, 2):
        for j in range(i - 1, i + 2):
            assert v[ix] == binomial(i, j)
            assert v[ix + 1] == binomial(i + h, j)
            assert v[ix + 2] == binomial(i, j + h)
            assert v[ix + 3] == binomial(i + h, j + h)
            ix += 4
    rw = binomial(x, y).rewrite(Piecewise)
    # make sure that special cases produce Float
    assert rw.xreplace({x: 2, y: 1.}).is_Float
    assert rw.xreplace({x: 3., y: 2}).is_Float
    for F in (Piecewise, factorial, ff, gamma, 'tractable', 'expand'):
        for i in range(-2, 2):
            for j in range(i - 1, i + 2):
                for hi, hj in cartes((0, S.Half), (0, S.Half)):
                    x, y = i + hi, j + hj
                    u = binomial(x, y, evaluate=False)
                    if F == 'expand':
                        rw = expand_func(u)
                    else:
                        rw = u.rewrite(F)
                    assert Eq(rw, binomial(x, y)), (x, y, F)



@XFAIL
def test_factorial_simplify_fail():
    # simplify(factorial(x + 1).diff(x) - ((x + 1)*factorial(x)).diff(x))) == 0
    from sympy.abc import x
    assert simplify(x*polygamma(0, x + 1) - x*polygamma(0, x + 2) +
                    polygamma(0, x + 1) - polygamma(0, x + 2) + 1) == 0


def test_subfactorial():
    assert all(subfactorial(i) == ans for i, ans in enumerate(
        [1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496]))
    assert subfactorial(oo) is oo
    assert subfactorial(nan) is nan
    assert subfactorial(23) == 9510425471055777937262
    assert unchanged(subfactorial, 2.2)

    x = Symbol('x')
    assert subfactorial(x).rewrite(uppergamma) == uppergamma(x + 1, -1)/S.Exp1

    tt = Symbol('tt', integer=True, nonnegative=True)
    tf = Symbol('tf', integer=True, nonnegative=False)
    tn = Symbol('tf', integer=True)
    ft = Symbol('ft', integer=False, nonnegative=True)
    ff = Symbol('ff', integer=False, nonnegative=False)
    fn = Symbol('ff', integer=False)
    nt = Symbol('nt', nonnegative=True)
    nf = Symbol('nf', nonnegative=False)
    nn = Symbol('nf')
    te = Symbol('te', even=True, nonnegative=True)
    to = Symbol('to', odd=True, nonnegative=True)
    assert subfactorial(tt).is_integer
    assert subfactorial(tf).is_integer is None
    assert subfactorial(tn).is_integer is None
    assert subfactorial(ft).is_integer is None
    assert subfactorial(ff).is_integer is None
    assert subfactorial(fn).is_integer is None
    assert subfactorial(nt).is_integer is None
    assert subfactorial(nf).is_integer is None
    assert subfactorial(nn).is_integer is None
    assert subfactorial(tt).is_nonnegative
    assert subfactorial(tf).is_nonnegative is None
    assert subfactorial(tn).is_nonnegative is None
    assert subfactorial(ft).is_nonnegative is None
    assert subfactorial(ff).is_nonnegative is None
    assert subfactorial(fn).is_nonnegative is None
    assert subfactorial(nt).is_nonnegative is None
    assert subfactorial(nf).is_nonnegative is None
    assert subfactorial(nn).is_nonnegative is None
    assert subfactorial(tt).is_even is None
    assert subfactorial(tt).is_odd is None
    assert subfactorial(te).is_odd is True
    assert subfactorial(to).is_even is True
