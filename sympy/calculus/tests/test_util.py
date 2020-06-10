from sympy import (Symbol, S, exp, log, sqrt, oo, E, zoo, pi, tan, sin, cos,
                   cot, sec, csc, Abs, symbols, I, re, simplify,
                   expint, Rational)
from sympy.calculus.util import (function_range, continuous_domain, not_empty_in,
                                 periodicity, lcim, AccumBounds, is_convex,
                                 stationary_points, minimum, maximum)
from sympy.core import Add, Mul, Pow
from sympy.sets.sets import (Interval, FiniteSet, EmptySet, Complement,
                            Union)
from sympy.testing.pytest import raises
from sympy.abc import x

a = Symbol('a', real=True)


def test_function_range():
    x, y, a, b = symbols('x y a b')
    assert function_range(sin(x), x, Interval(-pi/2, pi/2)
        ) == Interval(-1, 1)
    assert function_range(sin(x), x, Interval(0, pi)
        ) == Interval(0, 1)
    assert function_range(tan(x), x, Interval(0, pi)
        ) == Interval(-oo, oo)
    assert function_range(tan(x), x, Interval(pi/2, pi)
        ) == Interval(-oo, 0)
    assert function_range((x + 3)/(x - 2), x, Interval(-5, 5)
        ) == Union(Interval(-oo, Rational(2, 7)), Interval(Rational(8, 3), oo))
    assert function_range(1/(x**2), x, Interval(-1, 1)
        ) == Interval(1, oo)
    assert function_range(exp(x), x, Interval(-1, 1)
        ) == Interval(exp(-1), exp(1))
    assert function_range(log(x) - x, x, S.Reals
        ) == Interval(-oo, -1)
    assert function_range(sqrt(3*x - 1), x, Interval(0, 2)
        ) == Interval(0, sqrt(5))
    assert function_range(x*(x - 1) - (x**2 - x), x, S.Reals
        ) == FiniteSet(0)
    assert function_range(x*(x - 1) - (x**2 - x) + y, x, S.Reals
        ) == FiniteSet(y)
    assert function_range(sin(x), x, Union(Interval(-5, -3), FiniteSet(4))
        ) == Union(Interval(-sin(3), 1), FiniteSet(sin(4)))
    assert function_range(cos(x), x, Interval(-oo, -4)
        ) == Interval(-1, 1)
    assert function_range(cos(x), x, S.EmptySet) == S.EmptySet
    raises(NotImplementedError, lambda : function_range(
        exp(x)*(sin(x) - cos(x))/2 - x, x, S.Reals))
    raises(NotImplementedError, lambda : function_range(
        sin(x) + x, x, S.Reals)) # issue 13273
    raises(NotImplementedError, lambda : function_range(
        log(x), x, S.Integers))
    raises(NotImplementedError, lambda : function_range(
        sin(x)/2, x, S.Naturals))


def test_continuous_domain():
    x = Symbol('x')
    assert continuous_domain(sin(x), x, Interval(0, 2*pi)) == Interval(0, 2*pi)
    assert continuous_domain(tan(x), x, Interval(0, 2*pi)) == \
        Union(Interval(0, pi/2, False, True), Interval(pi/2, pi*Rational(3, 2), True, True),
              Interval(pi*Rational(3, 2), 2*pi, True, False))
    assert continuous_domain((x - 1)/((x - 1)**2), x, S.Reals) == \
        Union(Interval(-oo, 1, True, True), Interval(1, oo, True, True))
    assert continuous_domain(log(x) + log(4*x - 1), x, S.Reals) == \
        Interval(Rational(1, 4), oo, True, True)
    assert continuous_domain(1/sqrt(x - 3), x, S.Reals) == Interval(3, oo, True, True)
    assert continuous_domain(1/x - 2, x, S.Reals) == \
        Union(Interval.open(-oo, 0), Interval.open(0, oo))
    assert continuous_domain(1/(x**2 - 4) + 2, x, S.Reals) == \
        Union(Interval.open(-oo, -2), Interval.open(-2, 2), Interval.open(2, oo))
    domain = continuous_domain(log(tan(x)**2 + 1), x, S.Reals)
    assert not domain.contains(3*pi/2)
    assert domain.contains(5)


def test_not_empty_in():
    assert not_empty_in(FiniteSet(x, 2*x).intersect(Interval(1, 2, True, False)), x) == \
        Interval(S.Half, 2, True, False)
    assert not_empty_in(FiniteSet(x, x**2).intersect(Interval(1, 2)), x) == \
        Union(Interval(-sqrt(2), -1), Interval(1, 2))
    assert not_empty_in(FiniteSet(x**2 + x, x).intersect(Interval(2, 4)), x) == \
        Union(Interval(-sqrt(17)/2 - S.Half, -2),
              Interval(1, Rational(-1, 2) + sqrt(17)/2), Interval(2, 4))
    assert not_empty_in(FiniteSet(x/(x - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(a/(a - 1)).intersect(S.Reals), a) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet((x**2 - 3*x + 2)/(x - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(1))
    assert not_empty_in(FiniteSet(3, 4, x/(x - 1)).intersect(Interval(2, 3)), x) == \
        Interval(-oo, oo)
    assert not_empty_in(FiniteSet(4, x/(x - 1)).intersect(Interval(2, 3)), x) == \
        Interval(S(3)/2, 2)
    assert not_empty_in(FiniteSet(x/(x**2 - 1)).intersect(S.Reals), x) == \
        Complement(S.Reals, FiniteSet(-1, 1))
    assert not_empty_in(FiniteSet(x, x**2).intersect(Union(Interval(1, 3, True, True),
                                                           Interval(4, 5))), x) == \
        Union(Interval(-sqrt(5), -2), Interval(-sqrt(3), -1, True, True),
              Interval(1, 3, True, True), Interval(4, 5))
    assert not_empty_in(FiniteSet(1).intersect(Interval(3, 4)), x) == S.EmptySet
    assert not_empty_in(FiniteSet(x**2/(x + 2)).intersect(Interval(1, oo)), x) == \
        Union(Interval(-2, -1, True, False), Interval(2, oo))
    raises(ValueError, lambda: not_empty_in(x))
    raises(ValueError, lambda: not_empty_in(Interval(0, 1), x))
    raises(NotImplementedError,
           lambda: not_empty_in(FiniteSet(x).intersect(S.Reals), x, a))


def test_periodicity():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z', real=True)

    assert periodicity(sin(2*x), x) == pi
    assert periodicity((-2)*tan(4*x), x) == pi/4
    assert periodicity(sin(x)**2, x) == 2*pi
    assert periodicity(3**tan(3*x), x) == pi/3
    assert periodicity(tan(x)*cos(x), x) == 2*pi
    assert periodicity(sin(x)**(tan(x)), x) == 2*pi
    assert periodicity(tan(x)*sec(x), x) == 2*pi
    assert periodicity(sin(2*x)*cos(2*x) - y, x) == pi/2
    assert periodicity(tan(x) + cot(x), x) == pi
    assert periodicity(sin(x) - cos(2*x), x) == 2*pi
    assert periodicity(sin(x) - 1, x) == 2*pi
    assert periodicity(sin(4*x) + sin(x)*cos(x), x) == pi
    assert periodicity(exp(sin(x)), x) == 2*pi
    assert periodicity(log(cot(2*x)) - sin(cos(2*x)), x) == pi
    assert periodicity(sin(2*x)*exp(tan(x) - csc(2*x)), x) == pi
    assert periodicity(cos(sec(x) - csc(2*x)), x) == 2*pi
    assert periodicity(tan(sin(2*x)), x) == pi
    assert periodicity(2*tan(x)**2, x) == pi
    assert periodicity(sin(x%4), x) == 4
    assert periodicity(sin(x)%4, x) == 2*pi
    assert periodicity(tan((3*x-2)%4), x) == Rational(4, 3)
    assert periodicity((sqrt(2)*(x+1)+x) % 3, x) == 3 / (sqrt(2)+1)
    assert periodicity((x**2+1) % x, x) is None
    assert periodicity(sin(re(x)), x) == 2*pi
    assert periodicity(sin(x)**2 + cos(x)**2, x) is S.Zero
    assert periodicity(tan(x), y) is S.Zero
    assert periodicity(sin(x) + I*cos(x), x) == 2*pi
    assert periodicity(x - sin(2*y), y) == pi

    assert periodicity(exp(x), x) is None
    assert periodicity(exp(I*x), x) == 2*pi
    assert periodicity(exp(I*z), z) == 2*pi
    assert periodicity(exp(z), z) is None
    assert periodicity(exp(log(sin(z) + I*cos(2*z)), evaluate=False), z) == 2*pi
    assert periodicity(exp(log(sin(2*z) + I*cos(z)), evaluate=False), z) == 2*pi
    assert periodicity(exp(sin(z)), z) == 2*pi
    assert periodicity(exp(2*I*z), z) == pi
    assert periodicity(exp(z + I*sin(z)), z) is None
    assert periodicity(exp(cos(z/2) + sin(z)), z) == 4*pi
    assert periodicity(log(x), x) is None
    assert periodicity(exp(x)**sin(x), x) is None
    assert periodicity(sin(x)**y, y) is None

    assert periodicity(Abs(sin(Abs(sin(x)))), x) == pi
    assert all(periodicity(Abs(f(x)), x) == pi for f in (
        cos, sin, sec, csc, tan, cot))
    assert periodicity(Abs(sin(tan(x))), x) == pi
    assert periodicity(Abs(sin(sin(x) + tan(x))), x) == 2*pi
    assert periodicity(sin(x) > S.Half, x) == 2*pi

    assert periodicity(x > 2, x) is None
    assert periodicity(x**3 - x**2 + 1, x) is None
    assert periodicity(Abs(x), x) is None
    assert periodicity(Abs(x**2 - 1), x) is None

    assert periodicity((x**2 + 4)%2, x) is None
    assert periodicity((E**x)%3, x) is None

    assert periodicity(sin(expint(1, x))/expint(1, x), x) is None


def test_periodicity_check():
    x = Symbol('x')
    y = Symbol('y')

    assert periodicity(tan(x), x, check=True) == pi
    assert periodicity(sin(x) + cos(x), x, check=True) == 2*pi
    assert periodicity(sec(x), x) == 2*pi
    assert periodicity(sin(x*y), x) == 2*pi/abs(y)
    assert periodicity(Abs(sec(sec(x))), x) == pi


def test_lcim():
    from sympy import pi

    assert lcim([S.Half, S(2), S(3)]) == 6
    assert lcim([pi/2, pi/4, pi]) == pi
    assert lcim([2*pi, pi/2]) == 2*pi
    assert lcim([S.One, 2*pi]) is None
    assert lcim([S(2) + 2*E, E/3 + Rational(1, 3), S.One + E]) == S(2) + 2*E


def test_is_convex():
    assert is_convex(1/x, x, domain=Interval(0, oo)) == True
    assert is_convex(1/x, x, domain=Interval(-oo, 0)) == False
    assert is_convex(x**2, x, domain=Interval(0, oo)) == True
    assert is_convex(log(x), x) == False
    raises(NotImplementedError, lambda: is_convex(log(x), x, a))


def test_stationary_points():
    x, y = symbols('x y')

    assert stationary_points(sin(x), x, Interval(-pi/2, pi/2)
        ) == {-pi/2, pi/2}
    assert  stationary_points(sin(x), x, Interval.Ropen(0, pi/4)
        ) == EmptySet()
    assert stationary_points(tan(x), x,
        ) == EmptySet()
    assert stationary_points(sin(x)*cos(x), x, Interval(0, pi)
        ) == {pi/4, pi*Rational(3, 4)}
    assert stationary_points(sec(x), x, Interval(0, pi)
        ) == {0, pi}
    assert stationary_points((x+3)*(x-2), x
        ) == FiniteSet(Rational(-1, 2))
    assert stationary_points((x + 3)/(x - 2), x, Interval(-5, 5)
        ) == EmptySet()
    assert stationary_points((x**2+3)/(x-2), x
        ) == {2 - sqrt(7), 2 + sqrt(7)}
    assert stationary_points((x**2+3)/(x-2), x, Interval(0, 5)
        ) == {2 + sqrt(7)}
    assert stationary_points(x**4 + x**3 - 5*x**2, x, S.Reals
        ) == FiniteSet(-2, 0, Rational(5, 4))
    assert stationary_points(exp(x), x
        ) == EmptySet()
    assert stationary_points(log(x) - x, x, S.Reals
        ) == {1}
    assert stationary_points(cos(x), x, Union(Interval(0, 5), Interval(-6, -3))
        ) == {0, -pi, pi}
    assert stationary_points(y, x, S.Reals
        ) == S.Reals
    assert stationary_points(y, x, S.EmptySet) == S.EmptySet


def test_maximum():
    x, y = symbols('x y')
    assert maximum(sin(x), x) is S.One
    assert maximum(sin(x), x, Interval(0, 1)) == sin(1)
    assert maximum(tan(x), x) is oo
    assert maximum(tan(x), x, Interval(-pi/4, pi/4)) is S.One
    assert maximum(sin(x)*cos(x), x, S.Reals) == S.Half
    assert simplify(maximum(sin(x)*cos(x), x, Interval(pi*Rational(3, 8), pi*Rational(5, 8)))
        ) == sqrt(2)/4
    assert maximum((x+3)*(x-2), x) is oo
    assert maximum((x+3)*(x-2), x, Interval(-5, 0)) == S(14)
    assert maximum((x+3)/(x-2), x, Interval(-5, 0)) == Rational(2, 7)
    assert simplify(maximum(-x**4-x**3+x**2+10, x)
        ) == 41*sqrt(41)/512 + Rational(5419, 512)
    assert maximum(exp(x), x, Interval(-oo, 2)) == exp(2)
    assert maximum(log(x) - x, x, S.Reals) is S.NegativeOne
    assert maximum(cos(x), x, Union(Interval(0, 5), Interval(-6, -3))
        ) is S.One
    assert maximum(cos(x)-sin(x), x, S.Reals) == sqrt(2)
    assert maximum(y, x, S.Reals) == y

    raises(ValueError, lambda : maximum(sin(x), x, S.EmptySet))
    raises(ValueError, lambda : maximum(log(cos(x)), x, S.EmptySet))
    raises(ValueError, lambda : maximum(1/(x**2 + y**2 + 1), x, S.EmptySet))
    raises(ValueError, lambda : maximum(sin(x), sin(x)))
    raises(ValueError, lambda : maximum(sin(x), x*y, S.EmptySet))
    raises(ValueError, lambda : maximum(sin(x), S.One))


def test_minimum():
    x, y = symbols('x y')

    assert minimum(sin(x), x) is S.NegativeOne
    assert minimum(sin(x), x, Interval(1, 4)) == sin(4)
    assert minimum(tan(x), x) is -oo
    assert minimum(tan(x), x, Interval(-pi/4, pi/4)) is S.NegativeOne
    assert minimum(sin(x)*cos(x), x, S.Reals) == Rational(-1, 2)
    assert simplify(minimum(sin(x)*cos(x), x, Interval(pi*Rational(3, 8), pi*Rational(5, 8)))
        ) == -sqrt(2)/4
    assert minimum((x+3)*(x-2), x) == Rational(-25, 4)
    assert minimum((x+3)/(x-2), x, Interval(-5, 0)) == Rational(-3, 2)
    assert minimum(x**4-x**3+x**2+10, x) == S(10)
    assert minimum(exp(x), x, Interval(-2, oo)) == exp(-2)
    assert minimum(log(x) - x, x, S.Reals) is -oo
    assert minimum(cos(x), x, Union(Interval(0, 5), Interval(-6, -3))
        ) is S.NegativeOne
    assert minimum(cos(x)-sin(x), x, S.Reals) == -sqrt(2)
    assert minimum(y, x, S.Reals) == y

    raises(ValueError, lambda : minimum(sin(x), x, S.EmptySet))
    raises(ValueError, lambda : minimum(log(cos(x)), x, S.EmptySet))
    raises(ValueError, lambda : minimum(1/(x**2 + y**2 + 1), x, S.EmptySet))
    raises(ValueError, lambda : minimum(sin(x), sin(x)))
    raises(ValueError, lambda : minimum(sin(x), x*y, S.EmptySet))
    raises(ValueError, lambda : minimum(sin(x), S.One))


def test_AccumBounds():
    assert AccumBounds(1, 2).args == (1, 2)
    assert AccumBounds(1, 2).delta is S.One
    assert AccumBounds(1, 2).mid == Rational(3, 2)
    assert AccumBounds(1, 3).is_real == True

    assert AccumBounds(1, 1) is S.One

    assert AccumBounds(1, 2) + 1 == AccumBounds(2, 3)
    assert 1 + AccumBounds(1, 2) == AccumBounds(2, 3)
    assert AccumBounds(1, 2) + AccumBounds(2, 3) == AccumBounds(3, 5)

    assert -AccumBounds(1, 2) == AccumBounds(-2, -1)

    assert AccumBounds(1, 2) - 1 == AccumBounds(0, 1)
    assert 1 - AccumBounds(1, 2) == AccumBounds(-1, 0)
    assert AccumBounds(2, 3) - AccumBounds(1, 2) == AccumBounds(0, 2)

    assert x + AccumBounds(1, 2) == Add(AccumBounds(1, 2), x)
    assert a + AccumBounds(1, 2) == AccumBounds(1 + a, 2 + a)
    assert AccumBounds(1, 2) - x == Add(AccumBounds(1, 2), -x)

    assert AccumBounds(-oo, 1) + oo == AccumBounds(-oo, oo)
    assert AccumBounds(1, oo) + oo is oo
    assert AccumBounds(1, oo) - oo == AccumBounds(-oo, oo)
    assert (-oo - AccumBounds(-1, oo)) is -oo
    assert AccumBounds(-oo, 1) - oo is -oo

    assert AccumBounds(1, oo) - oo == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, 1) - (-oo) == AccumBounds(-oo, oo)
    assert (oo - AccumBounds(1, oo)) == AccumBounds(-oo, oo)
    assert (-oo - AccumBounds(1, oo)) is -oo

    assert AccumBounds(1, 2)/2 == AccumBounds(S.Half, 1)
    assert 2/AccumBounds(2, 3) == AccumBounds(Rational(2, 3), 1)
    assert 1/AccumBounds(-1, 1) == AccumBounds(-oo, oo)

    assert abs(AccumBounds(1, 2)) == AccumBounds(1, 2)
    assert abs(AccumBounds(-2, -1)) == AccumBounds(1, 2)
    assert abs(AccumBounds(-2, 1)) == AccumBounds(0, 2)
    assert abs(AccumBounds(-1, 2)) == AccumBounds(0, 2)
    c = Symbol('c')
    raises(ValueError, lambda: AccumBounds(0, c))
    raises(ValueError, lambda: AccumBounds(1, -1))


def test_AccumBounds_mul():
    assert AccumBounds(1, 2)*2 == AccumBounds(2, 4)
    assert 2*AccumBounds(1, 2) == AccumBounds(2, 4)
    assert AccumBounds(1, 2)*AccumBounds(2, 3) == AccumBounds(2, 6)

    assert AccumBounds(1, 2)*0 == 0
    assert AccumBounds(1, oo)*0 == AccumBounds(0, oo)
    assert AccumBounds(-oo, 1)*0 == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, oo)*0 == AccumBounds(-oo, oo)

    assert AccumBounds(1, 2)*x == Mul(AccumBounds(1, 2), x, evaluate=False)

    assert AccumBounds(0, 2)*oo == AccumBounds(0, oo)
    assert AccumBounds(-2, 0)*oo == AccumBounds(-oo, 0)
    assert AccumBounds(0, 2)*(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-2, 0)*(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-1, 1)*oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, 1)*(-oo) == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, oo)*oo == AccumBounds(-oo, oo)


def test_AccumBounds_div():
    assert AccumBounds(-1, 3)/AccumBounds(3, 4) == AccumBounds(Rational(-1, 3), 1)
    assert AccumBounds(-2, 4)/AccumBounds(-3, 4) == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)/AccumBounds(-4, 0) == AccumBounds(S.Half, oo)

    # these two tests can have a better answer
    # after Union of AccumBounds is improved
    assert AccumBounds(-3, -2)/AccumBounds(-2, 1) == AccumBounds(-oo, oo)
    assert AccumBounds(2, 3)/AccumBounds(-2, 2) == AccumBounds(-oo, oo)

    assert AccumBounds(-3, -2)/AccumBounds(0, 4) == AccumBounds(-oo, Rational(-1, 2))
    assert AccumBounds(2, 4)/AccumBounds(-3, 0) == AccumBounds(-oo, Rational(-2, 3))
    assert AccumBounds(2, 4)/AccumBounds(0, 3) == AccumBounds(Rational(2, 3), oo)

    assert AccumBounds(0, 1)/AccumBounds(0, 1) == AccumBounds(0, oo)
    assert AccumBounds(-1, 0)/AccumBounds(0, 1) == AccumBounds(-oo, 0)
    assert AccumBounds(-1, 2)/AccumBounds(-2, 2) == AccumBounds(-oo, oo)

    assert 1/AccumBounds(-1, 2) == AccumBounds(-oo, oo)
    assert 1/AccumBounds(0, 2) == AccumBounds(S.Half, oo)
    assert (-1)/AccumBounds(0, 2) == AccumBounds(-oo, Rational(-1, 2))
    assert 1/AccumBounds(-oo, 0) == AccumBounds(-oo, 0)
    assert 1/AccumBounds(-1, 0) == AccumBounds(-oo, -1)
    assert (-2)/AccumBounds(-oo, 0) == AccumBounds(0, oo)
    assert 1/AccumBounds(-oo, -1) == AccumBounds(-1, 0)

    assert AccumBounds(1, 2)/a == Mul(AccumBounds(1, 2), 1/a, evaluate=False)

    assert AccumBounds(1, 2)/0 == AccumBounds(1, 2)*zoo
    assert AccumBounds(1, oo)/oo == AccumBounds(0, oo)
    assert AccumBounds(1, oo)/(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, -1)/oo == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, -1)/(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-oo, oo)/oo == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, oo)/(-oo) == AccumBounds(-oo, oo)
    assert AccumBounds(-1, oo)/oo == AccumBounds(0, oo)
    assert AccumBounds(-1, oo)/(-oo) == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, 1)/oo == AccumBounds(-oo, 0)
    assert AccumBounds(-oo, 1)/(-oo) == AccumBounds(0, oo)

def test_issue_18795():
    r = Symbol('r', real=True)
    a = AccumBounds(-1,1)
    c = AccumBounds(7, oo)
    b = AccumBounds(-oo, oo)
    assert c - tan(r) == AccumBounds(7-tan(r), oo)
    assert b + tan(r) == AccumBounds(-oo, oo)
    assert (a + r)/a == AccumBounds(-oo, oo)*AccumBounds(r - 1, r + 1)
    assert (b + a)/a == AccumBounds(-oo, oo)

def test_AccumBounds_func():
    assert (x**2 + 2*x + 1).subs(x, AccumBounds(-1, 1)) == AccumBounds(-1, 4)
    assert exp(AccumBounds(0, 1)) == AccumBounds(1, E)
    assert exp(AccumBounds(-oo, oo)) == AccumBounds(0, oo)
    assert log(AccumBounds(3, 6)) == AccumBounds(log(3), log(6))


def test_AccumBounds_pow():
    assert AccumBounds(0, 2)**2 == AccumBounds(0, 4)
    assert AccumBounds(-1, 1)**2 == AccumBounds(0, 1)
    assert AccumBounds(1, 2)**2 == AccumBounds(1, 4)
    assert AccumBounds(-1, 2)**3 == AccumBounds(-1, 8)
    assert AccumBounds(-1, 1)**0 == 1

    assert AccumBounds(1, 2)**Rational(5, 2) == AccumBounds(1, 4*sqrt(2))
    assert AccumBounds(-1, 2)**Rational(1, 3) == AccumBounds(-1, 2**Rational(1, 3))
    assert AccumBounds(0, 2)**S.Half == AccumBounds(0, sqrt(2))

    assert AccumBounds(-4, 2)**Rational(2, 3) == AccumBounds(0, 2*2**Rational(1, 3))

    assert AccumBounds(-1, 5)**S.Half == AccumBounds(0, sqrt(5))
    assert AccumBounds(-oo, 2)**S.Half == AccumBounds(0, sqrt(2))
    assert AccumBounds(-2, 3)**Rational(-1, 4) == AccumBounds(0, oo)

    assert AccumBounds(1, 5)**(-2) == AccumBounds(Rational(1, 25), 1)
    assert AccumBounds(-1, 3)**(-2) == AccumBounds(0, oo)
    assert AccumBounds(0, 2)**(-2) == AccumBounds(Rational(1, 4), oo)
    assert AccumBounds(-1, 2)**(-3) == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)**(-3) == AccumBounds(Rational(-1, 8), Rational(-1, 27))
    assert AccumBounds(-3, -2)**(-2) == AccumBounds(Rational(1, 9), Rational(1, 4))
    assert AccumBounds(0, oo)**S.Half == AccumBounds(0, oo)
    assert AccumBounds(-oo, -1)**Rational(1, 3) == AccumBounds(-oo, -1)
    assert AccumBounds(-2, 3)**(Rational(-1, 3)) == AccumBounds(-oo, oo)
    assert AccumBounds(-oo, 0)**(-2) == AccumBounds(0, oo)
    assert AccumBounds(-2, 0)**(-2) == AccumBounds(Rational(1, 4), oo)

    assert AccumBounds(Rational(1, 3), S.Half)**oo is S.Zero
    assert AccumBounds(0, S.Half)**oo is S.Zero
    assert AccumBounds(S.Half, 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(0, 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(2, 3)**oo is oo
    assert AccumBounds(1, 2)**oo == AccumBounds(0, oo)
    assert AccumBounds(S.Half, 3)**oo == AccumBounds(0, oo)
    assert AccumBounds(Rational(-1, 3), Rational(-1, 4))**oo is S.Zero
    assert AccumBounds(-1, Rational(-1, 2))**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-3, -2)**oo == FiniteSet(-oo, oo)
    assert AccumBounds(-2, -1)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-2, Rational(-1, 2))**oo == AccumBounds(-oo, oo)
    assert AccumBounds(Rational(-1, 2), S.Half)**oo is S.Zero
    assert AccumBounds(Rational(-1, 2), 1)**oo == AccumBounds(0, oo)
    assert AccumBounds(Rational(-2, 3), 2)**oo == AccumBounds(0, oo)
    assert AccumBounds(-1, 1)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, S.Half)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-1, 2)**oo == AccumBounds(-oo, oo)
    assert AccumBounds(-2, S.Half)**oo == AccumBounds(-oo, oo)

    assert AccumBounds(1, 2)**x == Pow(AccumBounds(1, 2), x)

    assert AccumBounds(2, 3)**(-oo) is S.Zero
    assert AccumBounds(0, 2)**(-oo) == AccumBounds(0, oo)
    assert AccumBounds(-1, 2)**(-oo) == AccumBounds(-oo, oo)

    assert (tan(x)**sin(2*x)).subs(x, AccumBounds(0, pi/2)) == \
        Pow(AccumBounds(-oo, oo), AccumBounds(0, 1))


def test_comparison_AccumBounds():
    assert (AccumBounds(1, 3) < 4) == S.true
    assert (AccumBounds(1, 3) < -1) == S.false
    assert (AccumBounds(1, 3) < 2).rel_op == '<'
    assert (AccumBounds(1, 3) <= 2).rel_op == '<='

    assert (AccumBounds(1, 3) > 4) == S.false
    assert (AccumBounds(1, 3) > -1) == S.true
    assert (AccumBounds(1, 3) > 2).rel_op == '>'
    assert (AccumBounds(1, 3) >= 2).rel_op == '>='

    assert (AccumBounds(1, 3) < AccumBounds(4, 6)) == S.true
    assert (AccumBounds(1, 3) < AccumBounds(2, 4)).rel_op == '<'
    assert (AccumBounds(1, 3) < AccumBounds(-2, 0)) == S.false

    assert (AccumBounds(1, 3) <= AccumBounds(4, 6)) == S.true
    assert (AccumBounds(1, 3) <= AccumBounds(-2, 0)) == S.false

    assert (AccumBounds(1, 3) > AccumBounds(4, 6)) == S.false
    assert (AccumBounds(1, 3) > AccumBounds(-2, 0)) == S.true

    assert (AccumBounds(1, 3) >= AccumBounds(4, 6)) == S.false
    assert (AccumBounds(1, 3) >= AccumBounds(-2, 0)) == S.true

    # issue 13499
    assert (cos(x) > 0).subs(x, oo) == (AccumBounds(-1, 1) > 0)

    c = Symbol('c')
    raises(TypeError, lambda: (AccumBounds(0, 1) < c))
    raises(TypeError, lambda: (AccumBounds(0, 1) <= c))
    raises(TypeError, lambda: (AccumBounds(0, 1) > c))
    raises(TypeError, lambda: (AccumBounds(0, 1) >= c))


def test_contains_AccumBounds():
    assert (1 in AccumBounds(1, 2)) == S.true
    raises(TypeError, lambda: a in AccumBounds(1, 2))
    assert 0 in AccumBounds(-1, 0)
    raises(TypeError, lambda:
        (cos(1)**2 + sin(1)**2 - 1) in AccumBounds(-1, 0))
    assert (-oo in AccumBounds(1, oo)) == S.true
    assert (oo in AccumBounds(-oo, 0)) == S.true

    # issue 13159
    assert Mul(0, AccumBounds(-1, 1)) == Mul(AccumBounds(-1, 1), 0) == 0
    import itertools
    for perm in itertools.permutations([0, AccumBounds(-1, 1), x]):
        assert Mul(*perm) == 0


def test_intersection_AccumBounds():
    assert AccumBounds(0, 3).intersection(AccumBounds(1, 2)) == AccumBounds(1, 2)
    assert AccumBounds(0, 3).intersection(AccumBounds(1, 4)) == AccumBounds(1, 3)
    assert AccumBounds(0, 3).intersection(AccumBounds(-1, 2)) == AccumBounds(0, 2)
    assert AccumBounds(0, 3).intersection(AccumBounds(-1, 4)) == AccumBounds(0, 3)
    assert AccumBounds(0, 1).intersection(AccumBounds(2, 3)) == S.EmptySet
    raises(TypeError, lambda: AccumBounds(0, 3).intersection(1))


def test_union_AccumBounds():
    assert AccumBounds(0, 3).union(AccumBounds(1, 2)) == AccumBounds(0, 3)
    assert AccumBounds(0, 3).union(AccumBounds(1, 4)) == AccumBounds(0, 4)
    assert AccumBounds(0, 3).union(AccumBounds(-1, 2)) == AccumBounds(-1, 3)
    assert AccumBounds(0, 3).union(AccumBounds(-1, 4)) == AccumBounds(-1, 4)
    raises(TypeError, lambda: AccumBounds(0, 3).union(1))


def test_issue_16469():
    x = Symbol("x", real=True)
    f = abs(x)
    assert function_range(f, x, S.Reals) == Interval(0, oo, False, True)

def test_issue_18747():
    assert periodicity(exp(pi*I*(x/4+S.Half/2)), x) == 8
