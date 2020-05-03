from sympy.core import (
    Rational, Symbol, S, Float, Integer, Mul, Number, Pow,
    Basic, I, nan, pi, symbols, oo, zoo, N)
from sympy.core.tests.test_evalf import NS
from sympy.core.function import expand_multinomial
from sympy.functions.elementary.miscellaneous import sqrt, cbrt
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.special.error_functions import erf
from sympy.functions.elementary.trigonometric import (
    sin, cos, tan, sec, csc, sinh, cosh, tanh, atan)
from sympy.series.order import O
from sympy.core.expr import unchanged


def test_rational():
    a = Rational(1, 5)

    r = sqrt(5)/5
    assert sqrt(a) == r
    assert 2*sqrt(a) == 2*r

    r = a*a**S.Half
    assert a**Rational(3, 2) == r
    assert 2*a**Rational(3, 2) == 2*r

    r = a**5*a**Rational(2, 3)
    assert a**Rational(17, 3) == r
    assert 2 * a**Rational(17, 3) == 2*r


def test_large_rational():
    e = (Rational(123712**12 - 1, 7) + Rational(1, 7))**Rational(1, 3)
    assert e == 234232585392159195136 * (Rational(1, 7)**Rational(1, 3))


def test_negative_real():
    def feq(a, b):
        return abs(a - b) < 1E-10

    assert feq(S.One / Float(-0.5), -Integer(2))


def test_expand():
    x = Symbol('x')
    assert (2**(-1 - x)).expand() == S.Half*2**(-x)


def test_issue_3449():
    #test if powers are simplified correctly
    #see also issue 3995
    x = Symbol('x')
    assert ((x**Rational(1, 3))**Rational(2)) == x**Rational(2, 3)
    assert (
        (x**Rational(3))**Rational(2, 5)) == (x**Rational(3))**Rational(2, 5)

    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    assert (a**2)**b == (abs(a)**b)**2
    assert sqrt(1/a) != 1/sqrt(a)  # e.g. for a = -1
    assert (a**3)**Rational(1, 3) != a
    assert (x**a)**b != x**(a*b)  # e.g. x = -1, a=2, b=1/2
    assert (x**.5)**b == x**(.5*b)
    assert (x**.5)**.5 == x**.25
    assert (x**2.5)**.5 != x**1.25  # e.g. for x = 5*I

    k = Symbol('k', integer=True)
    m = Symbol('m', integer=True)
    assert (x**k)**m == x**(k*m)
    assert Number(5)**Rational(2, 3) == Number(25)**Rational(1, 3)

    assert (x**.5)**2 == x**1.0
    assert (x**2)**k == (x**k)**2 == x**(2*k)

    a = Symbol('a', positive=True)
    assert (a**3)**Rational(2, 5) == a**Rational(6, 5)
    assert (a**2)**b == (a**b)**2
    assert (a**Rational(2, 3))**x == a**(x*Rational(2, 3)) != (a**x)**Rational(2, 3)


def test_issue_3866():
    assert --sqrt(sqrt(5) - 1) == sqrt(sqrt(5) - 1)


def test_negative_one():
    x = Symbol('x', complex=True)
    y = Symbol('y', complex=True)
    assert 1/x**y == x**(-y)


def test_issue_4362():
    neg = Symbol('neg', negative=True)
    nonneg = Symbol('nonneg', nonnegative=True)
    any = Symbol('any')
    num, den = sqrt(1/neg).as_numer_denom()
    assert num == sqrt(-1)
    assert den == sqrt(-neg)
    num, den = sqrt(1/nonneg).as_numer_denom()
    assert num == 1
    assert den == sqrt(nonneg)
    num, den = sqrt(1/any).as_numer_denom()
    assert num == sqrt(1/any)
    assert den == 1

    def eqn(num, den, pow):
        return (num/den)**pow
    npos = 1
    nneg = -1
    dpos = 2 - sqrt(3)
    dneg = 1 - sqrt(3)
    assert dpos > 0 and dneg < 0 and npos > 0 and nneg < 0
    # pos or neg integer
    eq = eqn(npos, dpos, 2)
    assert eq.is_Pow and eq.as_numer_denom() == (1, dpos**2)
    eq = eqn(npos, dneg, 2)
    assert eq.is_Pow and eq.as_numer_denom() == (1, dneg**2)
    eq = eqn(nneg, dpos, 2)
    assert eq.is_Pow and eq.as_numer_denom() == (1, dpos**2)
    eq = eqn(nneg, dneg, 2)
    assert eq.is_Pow and eq.as_numer_denom() == (1, dneg**2)
    eq = eqn(npos, dpos, -2)
    assert eq.is_Pow and eq.as_numer_denom() == (dpos**2, 1)
    eq = eqn(npos, dneg, -2)
    assert eq.is_Pow and eq.as_numer_denom() == (dneg**2, 1)
    eq = eqn(nneg, dpos, -2)
    assert eq.is_Pow and eq.as_numer_denom() == (dpos**2, 1)
    eq = eqn(nneg, dneg, -2)
    assert eq.is_Pow and eq.as_numer_denom() == (dneg**2, 1)
    # pos or neg rational
    pow = S.Half
    eq = eqn(npos, dpos, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq = eqn(npos, dneg, pow)
    assert eq.is_Pow is False and eq.as_numer_denom() == ((-npos)**pow, (-dneg)**pow)
    eq = eqn(nneg, dpos, pow)
    assert not eq.is_Pow or eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq = eqn(nneg, dneg, pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq = eqn(npos, dpos, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq = eqn(npos, dneg, -pow)
    assert eq.is_Pow is False and eq.as_numer_denom() == (-(-npos)**pow*(-dneg)**pow, npos)
    eq = eqn(nneg, dpos, -pow)
    assert not eq.is_Pow or eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq = eqn(nneg, dneg, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)
    # unknown exponent
    pow = 2*any
    eq = eqn(npos, dpos, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq = eqn(npos, dneg, pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-npos)**pow, (-dneg)**pow)
    eq = eqn(nneg, dpos, pow)
    assert eq.is_Pow and eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq = eqn(nneg, dneg, pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq = eqn(npos, dpos, -pow)
    assert eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq = eqn(npos, dneg, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-npos)**pow)
    eq = eqn(nneg, dpos, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq = eqn(nneg, dneg, -pow)
    assert eq.is_Pow and eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)

    x = Symbol('x')
    y = Symbol('y')
    assert ((1/(1 + x/3))**(-S.One)).as_numer_denom() == (3 + x, 3)
    notp = Symbol('notp', positive=False)  # not positive does not imply real
    b = ((1 + x/notp)**-2)
    assert (b**(-y)).as_numer_denom() == (1, b**y)
    assert (b**(-S.One)).as_numer_denom() == ((notp + x)**2, notp**2)
    nonp = Symbol('nonp', nonpositive=True)
    assert (((1 + x/nonp)**-2)**(-S.One)).as_numer_denom() == ((-nonp -
            x)**2, nonp**2)

    n = Symbol('n', negative=True)
    assert (x**n).as_numer_denom() == (1, x**-n)
    assert sqrt(1/n).as_numer_denom() == (S.ImaginaryUnit, sqrt(-n))
    n = Symbol('0 or neg', nonpositive=True)
    # if x and n are split up without negating each term and n is negative
    # then the answer might be wrong; if n is 0 it won't matter since
    # 1/oo and 1/zoo are both zero as is sqrt(0)/sqrt(-x) unless x is also
    # zero (in which case the negative sign doesn't matter):
    # 1/sqrt(1/-1) = -I but sqrt(-1)/sqrt(1) = I
    assert (1/sqrt(x/n)).as_numer_denom() == (sqrt(-n), sqrt(-x))
    c = Symbol('c', complex=True)
    e = sqrt(1/c)
    assert e.as_numer_denom() == (e, 1)
    i = Symbol('i', integer=True)
    assert ((1 + x/y)**i).as_numer_denom() == ((x + y)**i, y**i)


def test_Pow_signs():
    """Cf. issues 4595 and 5250"""
    x = Symbol('x')
    y = Symbol('y')
    n = Symbol('n', even=True)
    assert (3 - y)**2 != (y - 3)**2
    assert (3 - y)**n != (y - 3)**n
    assert (-3 + y - x)**2 != (3 - y + x)**2
    assert (y - 3)**3 != -(3 - y)**3


def test_power_with_noncommutative_mul_as_base():
    x = Symbol('x', commutative=False)
    y = Symbol('y', commutative=False)
    assert not (x*y)**3 == x**3*y**3
    assert (2*x*y)**3 == 8*(x*y)**3


def test_power_rewrite_exp():
    assert (I**I).rewrite(exp) == exp(-pi/2)

    expr = (2 + 3*I)**(4 + 5*I)
    assert expr.rewrite(exp) == exp((4 + 5*I)*(log(sqrt(13)) + I*atan(Rational(3, 2))))
    assert expr.rewrite(exp).expand() == \
        169*exp(5*I*log(13)/2)*exp(4*I*atan(Rational(3, 2)))*exp(-5*atan(Rational(3, 2)))

    assert ((6 + 7*I)**5).rewrite(exp) == 7225*sqrt(85)*exp(5*I*atan(Rational(7, 6)))

    expr = 5**(6 + 7*I)
    assert expr.rewrite(exp) == exp((6 + 7*I)*log(5))
    assert expr.rewrite(exp).expand() == 15625*exp(7*I*log(5))

    assert Pow(123, 789, evaluate=False).rewrite(exp) == 123**789
    assert (1**I).rewrite(exp) == 1**I
    assert (0**I).rewrite(exp) == 0**I

    expr = (-2)**(2 + 5*I)
    assert expr.rewrite(exp) == exp((2 + 5*I)*(log(2) + I*pi))
    assert expr.rewrite(exp).expand() == 4*exp(-5*pi)*exp(5*I*log(2))

    assert ((-2)**S(-5)).rewrite(exp) == (-2)**S(-5)

    x, y = symbols('x y')
    assert (x**y).rewrite(exp) == exp(y*log(x))
    assert (7**x).rewrite(exp) == exp(x*log(7), evaluate=False)
    assert ((2 + 3*I)**x).rewrite(exp) == exp(x*(log(sqrt(13)) + I*atan(Rational(3, 2))))
    assert (y**(5 + 6*I)).rewrite(exp) == exp(log(y)*(5 + 6*I))

    assert all((1/func(x)).rewrite(exp) == 1/(func(x).rewrite(exp)) for func in
                    (sin, cos, tan, sec, csc, sinh, cosh, tanh))


def test_zero():
    x = Symbol('x')
    y = Symbol('y')
    assert 0**x != 0
    assert 0**(2*x) == 0**x
    assert 0**(1.0*x) == 0**x
    assert 0**(2.0*x) == 0**x
    assert (0**(2 - x)).as_base_exp() == (0, 2 - x)
    assert 0**(x - 2) != S.Infinity**(2 - x)
    assert 0**(2*x*y) == 0**(x*y)
    assert 0**(-2*x*y) == S.ComplexInfinity**(x*y)


def test_pow_as_base_exp():
    x = Symbol('x')
    assert (S.Infinity**(2 - x)).as_base_exp() == (S.Infinity, 2 - x)
    assert (S.Infinity**(x - 2)).as_base_exp() == (S.Infinity, x - 2)
    p = S.Half**x
    assert p.base, p.exp == p.as_base_exp() == (S(2), -x)
    # issue 8344:
    assert Pow(1, 2, evaluate=False).as_base_exp() == (S.One, S(2))


def test_issue_6100_12942_4473():
    x = Symbol('x')
    y = Symbol('y')
    assert x**1.0 != x
    assert x != x**1.0
    assert True != x**1.0
    assert x**1.0 is not True
    assert x is not True
    assert x*y != (x*y)**1.0
    # Pow != Symbol
    assert (x**1.0)**1.0 != x
    assert (x**1.0)**2.0 != x**2
    b = Basic()
    assert Pow(b, 1.0, evaluate=False) != b
    # if the following gets distributed as a Mul (x**1.0*y**1.0 then
    # __eq__ methods could be added to Symbol and Pow to detect the
    # power-of-1.0 case.
    assert ((x*y)**1.0).func is Pow


def test_issue_6208():
    from sympy import root, Rational
    I = S.ImaginaryUnit
    assert sqrt(33**(I*Rational(9, 10))) == -33**(I*Rational(9, 20))
    assert root((6*I)**(2*I), 3).as_base_exp()[1] == Rational(1, 3)  # != 2*I/3
    assert root((6*I)**(I/3), 3).as_base_exp()[1] == I/9
    assert sqrt(exp(3*I)) == exp(I*Rational(3, 2))
    assert sqrt(-sqrt(3)*(1 + 2*I)) == sqrt(sqrt(3))*sqrt(-1 - 2*I)
    assert sqrt(exp(5*I)) == -exp(I*Rational(5, 2))
    assert root(exp(5*I), 3).exp == Rational(1, 3)


def test_issue_6990():
    x = Symbol('x')
    a = Symbol('a')
    b = Symbol('b')
    assert (sqrt(a + b*x + x**2)).series(x, 0, 3).removeO() == \
        b*x/(2*sqrt(a)) + x**2*(1/(2*sqrt(a)) - \
        b**2/(8*a**Rational(3, 2))) + sqrt(a)


def test_issue_6068():
    x = Symbol('x')
    assert sqrt(sin(x)).series(x, 0, 7) == \
        sqrt(x) - x**Rational(5, 2)/12 + x**Rational(9, 2)/1440 - \
        x**Rational(13, 2)/24192 + O(x**7)
    assert sqrt(sin(x)).series(x, 0, 9) == \
        sqrt(x) - x**Rational(5, 2)/12 + x**Rational(9, 2)/1440 - \
        x**Rational(13, 2)/24192 - 67*x**Rational(17, 2)/29030400 + O(x**9)
    assert sqrt(sin(x**3)).series(x, 0, 19) == \
        x**Rational(3, 2) - x**Rational(15, 2)/12 + x**Rational(27, 2)/1440 + O(x**19)
    assert sqrt(sin(x**3)).series(x, 0, 20) == \
        x**Rational(3, 2) - x**Rational(15, 2)/12 + x**Rational(27, 2)/1440 - \
        x**Rational(39, 2)/24192 + O(x**20)


def test_issue_6782():
    x = Symbol('x')
    assert sqrt(sin(x**3)).series(x, 0, 7) == x**Rational(3, 2) + O(x**7)
    assert sqrt(sin(x**4)).series(x, 0, 3) == x**2 + O(x**3)


def test_issue_6653():
    x = Symbol('x')
    assert (1 / sqrt(1 + sin(x**2))).series(x, 0, 3) == 1 - x**2/2 + O(x**3)


def test_issue_6429():
    x = Symbol('x')
    c = Symbol('c')
    f = (c**2 + x)**(0.5)
    assert f.series(x, x0=0, n=1) == (c**2)**0.5 + O(x)
    assert f.taylor_term(0, x) == (c**2)**0.5
    assert f.taylor_term(1, x) == 0.5*x*(c**2)**(-0.5)
    assert f.taylor_term(2, x) == -0.125*x**2*(c**2)**(-1.5)


def test_issue_7638():
    f = pi/log(sqrt(2))
    assert ((1 + I)**(I*f/2))**0.3 == (1 + I)**(0.15*I*f)
    # if 1/3 -> 1.0/3 this should fail since it cannot be shown that the
    # sign will be +/-1; for the previous "small arg" case, it didn't matter
    # that this could not be proved
    assert (1 + I)**(4*I*f) == ((1 + I)**(12*I*f))**Rational(1, 3)

    assert (((1 + I)**(I*(1 + 7*f)))**Rational(1, 3)).exp == Rational(1, 3)
    r = symbols('r', real=True)
    assert sqrt(r**2) == abs(r)
    assert cbrt(r**3) != r
    assert sqrt(Pow(2*I, 5*S.Half)) != (2*I)**Rational(5, 4)
    p = symbols('p', positive=True)
    assert cbrt(p**2) == p**Rational(2, 3)
    assert NS(((0.2 + 0.7*I)**(0.7 + 1.0*I))**(0.5 - 0.1*I), 1) == '0.4 + 0.2*I'
    assert sqrt(1/(1 + I)) == sqrt(1 - I)/sqrt(2)  # or 1/sqrt(1 + I)
    e = 1/(1 - sqrt(2))
    assert sqrt(e) == I/sqrt(-1 + sqrt(2))
    assert e**Rational(-1, 2) == -I*sqrt(-1 + sqrt(2))
    assert sqrt((cos(1)**2 + sin(1)**2 - 1)**(3 + I)).exp in [S.Half,
                                                              Rational(3, 2) + I/2]
    assert sqrt(r**Rational(4, 3)) != r**Rational(2, 3)
    assert sqrt((p + I)**Rational(4, 3)) == (p + I)**Rational(2, 3)
    assert sqrt((p - p**2*I)**2) == p - p**2*I
    assert sqrt((p + r*I)**2) != p + r*I
    e = (1 + I/5)
    assert sqrt(e**5) == e**(5*S.Half)
    assert sqrt(e**6) == e**3
    assert sqrt((1 + I*r)**6) != (1 + I*r)**3


def test_issue_8582():
    assert 1**oo is nan
    assert 1**(-oo) is nan
    assert 1**zoo is nan
    assert 1**(oo + I) is nan
    assert 1**(1 + I*oo) is nan
    assert 1**(oo + I*oo) is nan


def test_issue_8650():
    n = Symbol('n', integer=True, nonnegative=True)
    assert (n**n).is_positive is True
    x = 5*n + 5
    assert (x**(5*(n + 1))).is_positive is True


def test_issue_13914():
    b = Symbol('b')
    assert (-1)**zoo is nan
    assert 2**zoo is nan
    assert (S.Half)**(1 + zoo) is nan
    assert I**(zoo + I) is nan
    assert b**(I + zoo) is nan


def test_better_sqrt():
    n = Symbol('n', integer=True, nonnegative=True)
    assert sqrt(3 + 4*I) == 2 + I
    assert sqrt(3 - 4*I) == 2 - I
    assert sqrt(-3 - 4*I) == 1 - 2*I
    assert sqrt(-3 + 4*I) == 1 + 2*I
    assert sqrt(32 + 24*I) == 6 + 2*I
    assert sqrt(32 - 24*I) == 6 - 2*I
    assert sqrt(-32 - 24*I) == 2 - 6*I
    assert sqrt(-32 + 24*I) == 2 + 6*I

    # triple (3, 4, 5):
    # parity of 3 matches parity of 5 and
    # den, 4, is a square
    assert sqrt((3 + 4*I)/4) == 1 + I/2
    # triple (8, 15, 17)
    # parity of 8 doesn't match parity of 17 but
    # den/2, 8/2, is a square
    assert sqrt((8 + 15*I)/8) == (5 + 3*I)/4
    # handle the denominator
    assert sqrt((3 - 4*I)/25) == (2 - I)/5
    assert sqrt((3 - 4*I)/26) == (2 - I)/sqrt(26)
    # mul
    #  issue #12739
    assert sqrt((3 + 4*I)/(3 - 4*I)) == (3 + 4*I)/5
    assert sqrt(2/(3 + 4*I)) == sqrt(2)/5*(2 - I)
    assert sqrt(n/(3 + 4*I)).subs(n, 2) == sqrt(2)/5*(2 - I)
    assert sqrt(-2/(3 + 4*I)) == sqrt(2)/5*(1 + 2*I)
    assert sqrt(-n/(3 + 4*I)).subs(n, 2) == sqrt(2)/5*(1 + 2*I)
    # power
    assert sqrt(1/(3 + I*4)) == (2 - I)/5
    assert sqrt(1/(3 - I)) == sqrt(10)*sqrt(3 + I)/10
    # symbolic
    i = symbols('i', imaginary=True)
    assert sqrt(3/i) == Mul(sqrt(3), 1/sqrt(i), evaluate=False)
    # multiples of 1/2; don't make this too automatic
    assert sqrt(3 + 4*I)**3 == (2 + I)**3
    assert Pow(3 + 4*I, Rational(3, 2)) == 2 + 11*I
    assert Pow(6 + 8*I, Rational(3, 2)) == 2*sqrt(2)*(2 + 11*I)
    n, d = (3 + 4*I), (3 - 4*I)**3
    a = n/d
    assert a.args == (1/d, n)
    eq = sqrt(a)
    assert eq.args == (a, S.Half)
    assert expand_multinomial(eq) == sqrt((-117 + 44*I)*(3 + 4*I))/125
    assert eq.expand() == (7 - 24*I)/125

    # issue 12775
    # pos im part
    assert sqrt(2*I) == (1 + I)
    assert sqrt(2*9*I) == Mul(3, 1 + I, evaluate=False)
    assert Pow(2*I, 3*S.Half) == (1 + I)**3
    # neg im part
    assert sqrt(-I/2) == Mul(S.Half, 1 - I, evaluate=False)
    # fractional im part
    assert Pow(Rational(-9, 2)*I, Rational(3, 2)) == 27*(1 - I)**3/8


def test_issue_2993():
    x = Symbol('x')
    assert str((2.3*x - 4)**0.3) == '1.5157165665104*(0.575*x - 1)**0.3'
    assert str((2.3*x + 4)**0.3) == '1.5157165665104*(0.575*x + 1)**0.3'
    assert str((-2.3*x + 4)**0.3) == '1.5157165665104*(1 - 0.575*x)**0.3'
    assert str((-2.3*x - 4)**0.3) == '1.5157165665104*(-0.575*x - 1)**0.3'
    assert str((2.3*x - 2)**0.3) == '1.28386201800527*(x - 0.869565217391304)**0.3'
    assert str((-2.3*x - 2)**0.3) == '1.28386201800527*(-x - 0.869565217391304)**0.3'
    assert str((-2.3*x + 2)**0.3) == '1.28386201800527*(0.869565217391304 - x)**0.3'
    assert str((2.3*x + 2)**0.3) == '1.28386201800527*(x + 0.869565217391304)**0.3'
    assert str((2.3*x - 4)**Rational(1, 3)) == '2**(2/3)*(0.575*x - 1)**(1/3)'
    eq = (2.3*x + 4)
    assert eq**2 == 16*(0.575*x + 1)**2
    assert (1/eq).args == (eq, -1)  # don't change trivial power
    # issue 17735
    q=.5*exp(x) - .5*exp(-x) + 0.1
    assert int((q**2).subs(x, 1)) == 1
    # issue 17756
    y = Symbol('y')
    assert len(sqrt(x/(x + y)**2 + Float('0.008', 30)).subs(y, pi.n(25)).atoms(Float)) == 2
    # issue 17756
    a, b, c, d, e, f, g = symbols('a:g')
    expr = sqrt(1 + a*(c**4 + g*d - 2*g*e - f*(-g + d))**2/
        (c**3*b**2*(d - 3*e + 2*f)**2))/2
    r = [
    (a, N('0.0170992456333788667034850458615', 30)),
    (b, N('0.0966594956075474769169134801223', 30)),
    (c, N('0.390911862903463913632151616184', 30)),
    (d, N('0.152812084558656566271750185933', 30)),
    (e, N('0.137562344465103337106561623432', 30)),
    (f, N('0.174259178881496659302933610355', 30)),
    (g, N('0.220745448491223779615401870086', 30))]
    tru = expr.n(30, subs=dict(r))
    seq = expr.subs(r)
    # although `tru` is the right way to evaluate
    # expr with numerical values, `seq` will have
    # significant loss of precision if extraction of
    # the largest coefficient of a power's base's terms
    # is done improperly
    assert seq == tru

def test_issue_17450():
    assert (erf(cosh(1)**7)**I).is_real is None
    assert (erf(cosh(1)**7)**I).is_imaginary is False
    assert (Pow(exp(1+sqrt(2)), ((1-sqrt(2))*I*pi), evaluate=False)).is_real is None
    assert ((-10)**(10*I*pi/3)).is_real is False
    assert ((-5)**(4*I*pi)).is_real is False


def test_issue_18190():
    assert sqrt(1 / tan(1 + I)) == 1 / sqrt(tan(1 + I))


def test_issue_14815():
    x = Symbol('x', real=True)
    assert sqrt(x).is_extended_negative is False
    x = Symbol('x', real=False)
    assert sqrt(x).is_extended_negative is None
    x = Symbol('x', complex=True)
    assert sqrt(x).is_extended_negative is False
    x = Symbol('x', extended_real=True)
    assert sqrt(x).is_extended_negative is False
    assert sqrt(zoo, evaluate=False).is_extended_negative is None
    assert sqrt(nan, evaluate=False).is_extended_negative is None


def test_issue_18509():
    assert unchanged(Mul, oo, 1/pi**oo)
    assert (1/pi**oo).is_extended_positive == False


def test_issue_18762():
    e, p = symbols('e p')
    g0 = sqrt(1 + e**2 - 2*e*cos(p))
    assert len(g0.series(e, 1, 3).args) == 4
