from sympy.core.numbers import (I, Rational, pi, oo)
from sympy.core.singleton import S
from sympy.core.symbol import Symbol, Dummy
from sympy.core.function import Lambda
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.trigonometric import sec, csc
from sympy.functions.elementary.hyperbolic import (coth, sech,
                                                   atanh, asech, acoth, acsch)
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.calculus.singularities import (
    singularities,
    is_increasing,
    is_strictly_increasing,
    is_decreasing,
    is_strictly_decreasing,
    is_monotonic
)
from sympy.sets import Interval, FiniteSet, Union, ImageSet
from sympy.testing.pytest import raises
from sympy.abc import x, y


def test_singularities():
    x = Symbol('x')
    assert singularities(x**2, x) == S.EmptySet
    assert singularities(x/(x**2 + 3*x + 2), x) == FiniteSet(-2, -1)
    assert singularities(1/(x**2 + 1), x) == FiniteSet(I, -I)
    assert singularities(x/(x**3 + 1), x) == \
        FiniteSet(-1, (1 - sqrt(3) * I) / 2, (1 + sqrt(3) * I) / 2)
    assert singularities(1/(y**2 + 2*I*y + 1), y) == \
        FiniteSet(-I + sqrt(2)*I, -I - sqrt(2)*I)
    _n = Dummy('n')
    assert singularities(sech(x), x).dummy_eq(Union(
        ImageSet(Lambda(_n, 2*_n*I*pi + I*pi/2), S.Integers),
        ImageSet(Lambda(_n, 2*_n*I*pi + 3*I*pi/2), S.Integers)))
    assert singularities(coth(x), x).dummy_eq(Union(
        ImageSet(Lambda(_n, 2*_n*I*pi + I*pi), S.Integers),
        ImageSet(Lambda(_n, 2*_n*I*pi), S.Integers)))
    assert singularities(atanh(x), x) == FiniteSet(-1, 1)
    assert singularities(acoth(x), x) == FiniteSet(-1, 1)
    assert singularities(asech(x), x) == FiniteSet(0)
    assert singularities(acsch(x), x) == FiniteSet(0)

    x = Symbol('x', real=True)
    assert singularities(1/(x**2 + 1), x) == S.EmptySet
    assert singularities(exp(1/x), x, S.Reals) == FiniteSet(0)
    assert singularities(exp(1/x), x, Interval(1, 2)) == S.EmptySet
    assert singularities(log((x - 2)**2), x, Interval(1, 3)) == FiniteSet(2)
    raises(NotImplementedError, lambda: singularities(x**-oo, x))
    assert singularities(sec(x), x, Interval(0, 3*pi)) == FiniteSet(
        pi/2, 3*pi/2, 5*pi/2)
    assert singularities(csc(x), x, Interval(0, 3*pi)) == FiniteSet(
        0, pi, 2*pi, 3*pi)
    assert singularities(0**x, x) == Interval.open(-oo, 0)
    assert singularities(0**(x-1), x) == Interval.open(-oo, 1)
    assert singularities(0**(x+1)/(x-2), x) == Union({2}, Interval.open(-oo, -1))


def test_is_increasing():
    """Test whether is_increasing returns correct value."""
    a = Symbol('a', negative=True)

    assert is_increasing(x**3 - 3*x**2 + 4*x, S.Reals)
    assert is_increasing(-x**2, Interval(-oo, 0))
    assert not is_increasing(-x**2, Interval(0, oo))
    assert not is_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval(-2, 3))
    assert is_increasing(x**2 + y, Interval(1, oo), x)
    assert is_increasing(-x**2*a, Interval(1, oo), x)
    assert is_increasing(1)

    assert is_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval(-2, 3)) is False


def test_is_strictly_increasing():
    """Test whether is_strictly_increasing returns correct value."""
    assert is_strictly_increasing(
        4*x**3 - 6*x**2 - 72*x + 30, Interval.Ropen(-oo, -2))
    assert is_strictly_increasing(
        4*x**3 - 6*x**2 - 72*x + 30, Interval.Lopen(3, oo))
    assert not is_strictly_increasing(
        4*x**3 - 6*x**2 - 72*x + 30, Interval.open(-2, 3))
    assert not is_strictly_increasing(-x**2, Interval(0, oo))
    assert not is_strictly_decreasing(1)

    assert is_strictly_increasing(4*x**3 - 6*x**2 - 72*x + 30, Interval.open(-2, 3)) is False


def test_is_decreasing():
    """Test whether is_decreasing returns correct value."""
    b = Symbol('b', positive=True)

    assert is_decreasing(1/(x**2 - 3*x), Interval.open(Rational(3,2), 3))
    assert is_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))
    assert is_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    assert not is_decreasing(1/(x**2 - 3*x), Interval.Ropen(-oo, Rational(3, 2)))
    assert not is_decreasing(-x**2, Interval(-oo, 0))
    assert not is_decreasing(-x**2*b, Interval(-oo, 0), x)


def test_is_strictly_decreasing():
    """Test whether is_strictly_decreasing returns correct value."""
    assert is_strictly_decreasing(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    assert not is_strictly_decreasing(
        1/(x**2 - 3*x), Interval.Ropen(-oo, Rational(3, 2)))
    assert not is_strictly_decreasing(-x**2, Interval(-oo, 0))
    assert not is_strictly_decreasing(1)
    assert is_strictly_decreasing(1/(x**2 - 3*x), Interval.open(Rational(3,2), 3))
    assert is_strictly_decreasing(1/(x**2 - 3*x), Interval.open(1.5, 3))


def test_is_monotonic():
    """Test whether is_monotonic returns correct value."""
    assert is_monotonic(1/(x**2 - 3*x), Interval.open(Rational(3,2), 3))
    assert is_monotonic(1/(x**2 - 3*x), Interval.open(1.5, 3))
    assert is_monotonic(1/(x**2 - 3*x), Interval.Lopen(3, oo))
    assert is_monotonic(x**3 - 3*x**2 + 4*x, S.Reals)
    assert not is_monotonic(-x**2, S.Reals)
    assert is_monotonic(x**2 + y + 1, Interval(1, 2), x)
    raises(NotImplementedError, lambda: is_monotonic(x**2 + y + 1))

    # Test cases for functions with derivative = 0 at isolated points.
    # x**3: derivative is 3*x**2, which is >= 0 everywhere (zero at x=0)
    assert is_monotonic(x**3, Interval(-10, 10))
    assert is_monotonic(x**3, S.Reals)

    # -x**3: derivative is -3*x**2, which is <= 0 everywhere (zero at x=0)
    assert is_monotonic(-x**3, Interval(-10, 10))
    assert is_monotonic(-x**3, S.Reals)

    # x**5: derivative is 5*x**4, which is >= 0 everywhere (zero at x=0)
    assert is_monotonic(x**5, S.Reals)

    # x**7: another odd power
    assert is_monotonic(x**7, S.Reals)
    assert is_monotonic(-x**7, S.Reals)

    # Test strictly increasing linear and polynomial functions
    assert is_monotonic(x, S.Reals)
    assert is_monotonic(2*x + 3, S.Reals)
    assert is_monotonic(x**3 + x, S.Reals)
    assert is_monotonic(x**3 + 2*x, S.Reals)

    # Test strictly decreasing functions
    assert is_monotonic(-x, S.Reals)
    assert is_monotonic(-2*x + 5, S.Reals)

    # Test non-monotonic polynomial functions
    assert not is_monotonic(x**2, S.Reals)  # Has minimum at x=0
    assert not is_monotonic(x**3 - 3*x, S.Reals)  # Has local max and min
    assert not is_monotonic(x**4 - x**2, S.Reals)
    assert not is_monotonic(x**4, S.Reals)  # Has minimum at x=0

    # Test monotonic on restricted intervals
    assert is_monotonic(x**2, Interval(0, oo))  # Increasing on [0, oo)
    assert is_monotonic(x**2, Interval.Lopen(0, oo))
    assert is_monotonic(x**2, Interval(1, 10))
    assert is_monotonic(-x**2, Interval(-oo, 0))  # Increasing on (-oo, 0]
    assert not is_monotonic(x**2, Interval(-1, 1))  # Not monotonic on [-1, 1]

    # Test constant functions (monotonic both ways)
    assert is_monotonic(5, S.Reals)
    assert is_monotonic(pi, Interval(-10, 10))
    assert is_monotonic(0, S.Reals)

    # Test rational functions on appropriate intervals
    assert is_monotonic(1/x, Interval.open(0, oo))  # Decreasing on (0, oo)
    assert is_monotonic(1/x, Interval.open(-oo, 0))  # Decreasing on (-oo, 0)
    assert is_monotonic(1/x, Interval(1, 10))
    assert is_monotonic(1/x, Interval(-10, -1))

    # Test polynomial with all stationary points
    assert is_monotonic(x**3 + 3*x**2 + 3*x + 1, S.Reals)  # (x+1)**3, always increasing
    assert is_monotonic((x + 2)**3, S.Reals)
    assert is_monotonic((x - 1)**5, S.Reals)

    # Test with specified symbol (multivariate)
    assert is_monotonic(x**2 + y, Interval(-oo, 0), x)  # Decreasing in x on (-oo, 0]
    assert is_monotonic(x**2 + y, Interval(0, oo), x)  # Increasing in x on [0, oo)
    assert is_monotonic(x + y**2, Interval(0, oo), y)  # Increasing in y on [0, oo)

    # Edge cases with single point intervals
    assert is_monotonic(x**2, Interval(5, 5))  # Single point
    assert is_monotonic(x**3 - x, Interval(0, 0))  # Single point

    # Test x**4 on restricted intervals where they're monotonic
    assert is_monotonic(x**4, Interval(0, oo))
    assert is_monotonic(x**4, Interval.Lopen(0, oo))
    assert is_monotonic(-x**4, Interval(-oo, 0))

    # Test more complex polynomials with all odd powers
    assert is_monotonic(x**7 + x**5 + x**3 + x, S.Reals)
    assert is_monotonic(-(x**7 + x**5 + x**3 + x), S.Reals)

    # Test polynomials with derivative that has zeros but maintains sign
    assert is_monotonic(x**5 + x**3, S.Reals)  # f'(x) = 5x**4 + 3x**2 >= 0
    assert is_monotonic(x**9 + x**7 + x**5 + x**3 + x, S.Reals)

def test_issue_23401():
    x = Symbol('x')
    expr = (x + 1)/(-1.0e-3*x**2 + 0.1*x + 0.1)
    assert is_increasing(expr, Interval(1,2), x)
