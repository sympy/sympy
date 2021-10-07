"""Tests on algebraic numbers. """

from sympy import S, Rational, Symbol, Poly, sqrt, I, Tuple, AlgebraicNumber, exp
from sympy.testing.pytest import raises, slow
from sympy.sets.sets import FiniteSet
from sympy import Point2D
from sympy.solvers.solveset import nonlinsolve
from sympy.geometry import Circle, intersection
from sympy.polys.partfrac import apart
from sympy.polys.numberfields.minpoly import minimal_polynomial
from sympy.polys.numberfields.numbers import (
    to_number_field,
    isolate,
    IntervalPrinter,
)
from sympy.polys.polyerrors import IsomorphismFailed, GeneratorsError, NotAlgebraic
from sympy.polys.polyclasses import DMP
from sympy.polys.domains import QQ
from sympy.abc import x, y, z

Q = Rational


def test_to_number_field():
    assert to_number_field(sqrt(2)) == AlgebraicNumber(sqrt(2))
    assert to_number_field(
        [sqrt(2), sqrt(3)]) == AlgebraicNumber(sqrt(2) + sqrt(3))

    a = AlgebraicNumber(sqrt(2) + sqrt(3), [S.Half, S.Zero, Rational(-9, 2), S.Zero])

    assert to_number_field(sqrt(2), sqrt(2) + sqrt(3)) == a
    assert to_number_field(sqrt(2), AlgebraicNumber(sqrt(2) + sqrt(3))) == a

    raises(IsomorphismFailed, lambda: to_number_field(sqrt(2), sqrt(3)))


def test_AlgebraicNumber():
    minpoly, root = x**2 - 2, sqrt(2)

    a = AlgebraicNumber(root, gen=x)

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is False

    assert a.coeffs() == [S.One, S.Zero]
    assert a.native_coeffs() == [QQ(1), QQ(0)]

    a = AlgebraicNumber(root, gen=x, alias='y')

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.alias == Symbol('y')
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is True

    a = AlgebraicNumber(root, gen=x, alias=Symbol('y'))

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.alias == Symbol('y')
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is True

    assert AlgebraicNumber(sqrt(2), []).rep == DMP([], QQ)
    assert AlgebraicNumber(sqrt(2), ()).rep == DMP([], QQ)
    assert AlgebraicNumber(sqrt(2), (0, 0)).rep == DMP([], QQ)

    assert AlgebraicNumber(sqrt(2), [8]).rep == DMP([QQ(8)], QQ)
    assert AlgebraicNumber(sqrt(2), [Rational(8, 3)]).rep == DMP([QQ(8, 3)], QQ)

    assert AlgebraicNumber(sqrt(2), [7, 3]).rep == DMP([QQ(7), QQ(3)], QQ)
    assert AlgebraicNumber(
        sqrt(2), [Rational(7, 9), Rational(3, 2)]).rep == DMP([QQ(7, 9), QQ(3, 2)], QQ)

    assert AlgebraicNumber(sqrt(2), [1, 2, 3]).rep == DMP([QQ(2), QQ(5)], QQ)

    a = AlgebraicNumber(AlgebraicNumber(root, gen=x), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is False

    assert a.coeffs() == [S.One, S(2)]
    assert a.native_coeffs() == [QQ(1), QQ(2)]

    a = AlgebraicNumber((minpoly, root), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is False

    a = AlgebraicNumber((Poly(minpoly), root), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.alias is None
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.is_aliased is False

    assert AlgebraicNumber( sqrt(3)).rep == DMP([ QQ(1), QQ(0)], QQ)
    assert AlgebraicNumber(-sqrt(3)).rep == DMP([ QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(2))

    assert a == b

    c = AlgebraicNumber(sqrt(2), gen=x)

    assert a == b
    assert a == c

    a = AlgebraicNumber(sqrt(2), [1, 2])
    b = AlgebraicNumber(sqrt(2), [1, 3])

    assert a != b and a != sqrt(2) + 3

    assert (a == x) is False and (a != x) is True

    a = AlgebraicNumber(sqrt(2), [1, 0])
    b = AlgebraicNumber(sqrt(2), [1, 0], alias=y)

    assert a.as_poly(x) == Poly(x, domain='QQ')
    assert b.as_poly() == Poly(y, domain='QQ')

    assert a.as_expr() == sqrt(2)
    assert a.as_expr(x) == x
    assert b.as_expr() == sqrt(2)
    assert b.as_expr(x) == x

    a = AlgebraicNumber(sqrt(2), [2, 3])
    b = AlgebraicNumber(sqrt(2), [2, 3], alias=y)

    p = a.as_poly()

    assert p == Poly(2*p.gen + 3)

    assert a.as_poly(x) == Poly(2*x + 3, domain='QQ')
    assert b.as_poly() == Poly(2*y + 3, domain='QQ')

    assert a.as_expr() == 2*sqrt(2) + 3
    assert a.as_expr(x) == 2*x + 3
    assert b.as_expr() == 2*sqrt(2) + 3
    assert b.as_expr(x) == 2*x + 3

    a = AlgebraicNumber(sqrt(2))
    b = to_number_field(sqrt(2))
    assert a.args == b.args == (sqrt(2), Tuple(1, 0))
    b = AlgebraicNumber(sqrt(2), alias='alpha')
    assert b.args == (sqrt(2), Tuple(1, 0), Symbol('alpha'))

    a = AlgebraicNumber(sqrt(2), [1, 2, 3])
    assert a.args == (sqrt(2), Tuple(1, 2, 3))


def test_to_algebraic_integer():
    a = AlgebraicNumber(sqrt(3), gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 3
    assert a.root == sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(2*sqrt(3), gen=x).to_algebraic_integer()
    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, [Rational(7, 19), 3], gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(7, 19), QQ(3)], QQ)


def test_IntervalPrinter():
    ip = IntervalPrinter()
    assert ip.doprint(x**Q(1, 3)) == "x**(mpi('1/3'))"
    assert ip.doprint(sqrt(x)) == "x**(mpi('1/2'))"


def test_isolate():
    assert isolate(1) == (1, 1)
    assert isolate(S.Half) == (S.Half, S.Half)

    assert isolate(sqrt(2)) == (1, 2)
    assert isolate(-sqrt(2)) == (-2, -1)

    assert isolate(sqrt(2), eps=Rational(1, 100)) == (Rational(24, 17), Rational(17, 12))
    assert isolate(-sqrt(2), eps=Rational(1, 100)) == (Rational(-17, 12), Rational(-24, 17))

    raises(NotImplementedError, lambda: isolate(I))

def test_minpoly_fraction_field():
    assert minimal_polynomial(1/x, y) == -x*y + 1
    assert minimal_polynomial(1 / (x + 1), y) == (x + 1)*y - 1

    assert minimal_polynomial(sqrt(x), y) == y**2 - x
    assert minimal_polynomial(sqrt(x + 1), y) == y**2 - x - 1
    assert minimal_polynomial(sqrt(x) / x, y) == x*y**2 - 1
    assert minimal_polynomial(sqrt(2) * sqrt(x), y) == y**2 - 2 * x
    assert minimal_polynomial(sqrt(2) + sqrt(x), y) == \
        y**4 + (-2*x - 4)*y**2 + x**2 - 4*x + 4

    assert minimal_polynomial(x**Rational(1,3), y) == y**3 - x
    assert minimal_polynomial(x**Rational(1,3) + sqrt(x), y) == \
        y**6 - 3*x*y**4 - 2*x*y**3 + 3*x**2*y**2 - 6*x**2*y - x**3 + x**2

    assert minimal_polynomial(sqrt(x) / z, y) == z**2*y**2 - x
    assert minimal_polynomial(sqrt(x) / (z + 1), y) == (z**2 + 2*z + 1)*y**2 - x

    assert minimal_polynomial(1/x, y, polys=True) == Poly(-x*y + 1, y, domain='ZZ(x)')
    assert minimal_polynomial(1 / (x + 1), y, polys=True) == \
        Poly((x + 1)*y - 1, y, domain='ZZ(x)')
    assert minimal_polynomial(sqrt(x), y, polys=True) == Poly(y**2 - x, y, domain='ZZ(x)')
    assert minimal_polynomial(sqrt(x) / z, y, polys=True) == \
        Poly(z**2*y**2 - x, y, domain='ZZ(x, z)')

    # this is (sqrt(1 + x**3)/x).integrate(x).diff(x) - sqrt(1 + x**3)/x
    a = sqrt(x)/sqrt(1 + x**(-3)) - sqrt(x**3 + 1)/x + 1/(x**Rational(5, 2)* \
        (1 + x**(-3))**Rational(3, 2)) + 1/(x**Rational(11, 2)*(1 + x**(-3))**Rational(3, 2))

    assert minimal_polynomial(a, y) == y

    raises(NotAlgebraic, lambda: minimal_polynomial(exp(x), y))
    raises(GeneratorsError, lambda: minimal_polynomial(sqrt(x), x))
    raises(GeneratorsError, lambda: minimal_polynomial(sqrt(x) - y, x))
    raises(NotImplementedError, lambda: minimal_polynomial(sqrt(x), y, compose=False))

@slow
def test_minpoly_fraction_field_slow():
    assert minimal_polynomial(minimal_polynomial(sqrt(x**Rational(1,5) - 1),
        y).subs(y, sqrt(x**Rational(1,5) - 1)), z) == z

def test_minpoly_domain():
    assert minimal_polynomial(sqrt(2), x, domain=QQ.algebraic_field(sqrt(2))) == \
        x - sqrt(2)
    assert minimal_polynomial(sqrt(8), x, domain=QQ.algebraic_field(sqrt(2))) == \
        x - 2*sqrt(2)
    assert minimal_polynomial(sqrt(Rational(3,2)), x,
        domain=QQ.algebraic_field(sqrt(2))) == 2*x**2 - 3

    raises(NotAlgebraic, lambda: minimal_polynomial(y, x, domain=QQ))


def test_issue_14831():
    a = -2*sqrt(2)*sqrt(12*sqrt(2) + 17)
    assert minimal_polynomial(a, x) == x**2 + 16*x - 8
    e = (-3*sqrt(12*sqrt(2) + 17) + 12*sqrt(2) +
         17 - 2*sqrt(2)*sqrt(12*sqrt(2) + 17))
    assert minimal_polynomial(e, x) == x


def test_issue_18248():
    assert nonlinsolve([x*y**3-sqrt(2)/3, x*y**6-4/(9*(sqrt(3)))],x,y) == \
            FiniteSet((sqrt(3)/2, sqrt(6)/3), (sqrt(3)/2, -sqrt(6)/6 - sqrt(2)*I/2),
            (sqrt(3)/2, -sqrt(6)/6 + sqrt(2)*I/2))


def test_issue_13230():
    c1 = Circle(Point2D(3, sqrt(5)), 5)
    c2 = Circle(Point2D(4, sqrt(7)), 6)
    assert intersection(c1, c2) == [Point2D(-1 + (-sqrt(7) + sqrt(5))*(-2*sqrt(7)/29
    + 9*sqrt(5)/29 + sqrt(196*sqrt(35) + 1941)/29), -2*sqrt(7)/29 + 9*sqrt(5)/29
    + sqrt(196*sqrt(35) + 1941)/29), Point2D(-1 + (-sqrt(7) + sqrt(5))*(-sqrt(196*sqrt(35)
    + 1941)/29 - 2*sqrt(7)/29 + 9*sqrt(5)/29), -sqrt(196*sqrt(35) + 1941)/29 - 2*sqrt(7)/29 + 9*sqrt(5)/29)]

def test_issue_19760():
    e = 1/(sqrt(1 + sqrt(2)) - sqrt(2)*sqrt(1 + sqrt(2))) + 1
    mp_expected = x**4 - 4*x**3 + 4*x**2 - 2

    for comp in (True, False):
        mp = Poly(minimal_polynomial(e, compose=comp))
        assert mp(x) == mp_expected, "minimal_polynomial(e, compose=%s) = %s; %s expected" % (comp, mp(x), mp_expected)


def test_issue_20163():
    assert apart(1/(x**6+1), extension=[sqrt(3), I]) == \
        (sqrt(3) + I)/(2*x + sqrt(3) + I)/6 + \
        (sqrt(3) - I)/(2*x + sqrt(3) - I)/6 - \
        (sqrt(3) - I)/(2*x - sqrt(3) + I)/6 - \
        (sqrt(3) + I)/(2*x - sqrt(3) - I)/6 + \
        I/(x + I)/6 - I/(x - I)/6
