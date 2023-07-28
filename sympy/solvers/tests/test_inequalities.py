"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy.concrete.summations import Sum
from sympy.core.function import Function
from sympy.core.numbers import I, Rational, oo, pi
from sympy.core.relational import Relational, Eq, Ge, Gt, Le, Lt, Ne
from sympy.core.singleton import S
from sympy.core.symbol import (Dummy, Symbol)
from sympy.core.sympify import sympify
from sympy.core.random import random, choice
from sympy.ntheory.generate import randprime
from sympy.assumptions.ask import Q
from sympy.matrices.dense import Matrix, eye
from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.miscellaneous import root, sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import cos, sin, tan
from sympy.integrals.integrals import Integral
from sympy.logic.boolalg import And, Or
from sympy.polys.polytools import Poly, PurePoly
from sympy.sets.sets import FiniteSet, Interval, Union
from sympy.solvers.inequalities import (reduce_inequalities,
                                        solve_poly_inequality as psolve,
                                        reduce_rational_inequalities,
                                        solve_univariate_inequality as isolve,
                                        reduce_abs_inequality,
                                        _solve_inequality,
                                        lp,
                                        UnboundedLinearProgrammingError,
                                        InfeasibleLinearProgrammingError)
from sympy.polys.rootoftools import rootof
from sympy.solvers.solvers import solve
from sympy.solvers.solveset import solveset
from sympy.abc import x, y, z
from sympy.core.symbol import symbols

from sympy.external.importtools import import_module


from sympy.core.mod import Mod

from sympy.testing.pytest import raises, XFAIL


inf = oo.evalf()


def test_solve_poly_inequality():
    assert psolve(Poly(0, x), '==') == [S.Reals]
    assert psolve(Poly(1, x), '==') == [S.EmptySet]
    assert psolve(PurePoly(x + 1, x), ">") == [Interval(-1, oo, True, False)]


def test_reduce_poly_inequalities_real_interval():
    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=False) == \
        S.Reals if x.is_real else Interval(-oo, oo)
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=False) == \
        FiniteSet(0).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=False) == \
        FiniteSet(0).complement(S.Reals)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[Le(x**2, 1)]], x, relational=False) == Interval(-1, 1)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 1)]], x, relational=False) == Interval(-1, 1, True, True)
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1)]], x, relational=False) == \
        Union(Interval(-oo, -1), Interval(1, oo))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1)]], x, relational=False) == \
        Interval(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 1)]], x, relational=False) == \
        FiniteSet(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities([[Eq(
        x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).evalf()
    assert reduce_rational_inequalities(
        [[Le(x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0)
    assert reduce_rational_inequalities([[Lt(
        x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0, True, True)
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1.0)]], x, relational=False) == \
        Union(Interval(-inf, -1.0), Interval(1.0, inf))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1.0)]], x, relational=False) == \
        Union(Interval(-inf, -1.0, right_open=True),
        Interval(1.0, inf, left_open=True))
    assert reduce_rational_inequalities([[Ne(
        x**2, 1.0)]], x, relational=False) == \
        FiniteSet(-1.0, 1.0).complement(S.Reals)

    s = sqrt(2)

    assert reduce_rational_inequalities([[Lt(
        x**2 - 1, 0), Gt(x**2 - 1, 0)]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities([[Le(x**2 - 1, 0), Ge(
        x**2 - 1, 0)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[Le(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False
        ) == Union(Interval(-s, -1, False, False), Interval(1, s, False, False))
    assert reduce_rational_inequalities(
        [[Le(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False
        ) == Union(Interval(-s, -1, False, True), Interval(1, s, True, False))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False
        ) == Union(Interval(-s, -1, True, False), Interval(1, s, False, True))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False
        ) == Union(Interval(-s, -1, True, True), Interval(1, s, True, True))
    assert reduce_rational_inequalities(
        [[Lt(x**2 - 2, 0), Ne(x**2 - 1, 0)]], x, relational=False
        ) == Union(Interval(-s, -1, True, True), Interval(-1, 1, True, True),
        Interval(1, s, True, True))

    assert reduce_rational_inequalities([[Lt(x**2, -1.)]], x) is S.false


def test_reduce_poly_inequalities_complex_relational():
    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=True) == False
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=True) == And(Lt(-oo, x), Lt(x, oo))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=True) == \
        And(Gt(x, -oo), Lt(x, oo), Ne(x, 0))
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=True) == \
        And(Gt(x, -oo), Lt(x, oo), Ne(x, 0))

    for one in (S.One, S(1.0)):
        inf = one*oo
        assert reduce_rational_inequalities(
            [[Eq(x**2, one)]], x, relational=True) == \
            Or(Eq(x, -one), Eq(x, one))
        assert reduce_rational_inequalities(
            [[Le(x**2, one)]], x, relational=True) == \
            And(And(Le(-one, x), Le(x, one)))
        assert reduce_rational_inequalities(
            [[Lt(x**2, one)]], x, relational=True) == \
            And(And(Lt(-one, x), Lt(x, one)))
        assert reduce_rational_inequalities(
            [[Ge(x**2, one)]], x, relational=True) == \
            And(Or(And(Le(one, x), Lt(x, inf)), And(Le(x, -one), Lt(-inf, x))))
        assert reduce_rational_inequalities(
            [[Gt(x**2, one)]], x, relational=True) == \
            And(Or(And(Lt(-inf, x), Lt(x, -one)), And(Lt(one, x), Lt(x, inf))))
        assert reduce_rational_inequalities(
            [[Ne(x**2, one)]], x, relational=True) == \
            Or(And(Lt(-inf, x), Lt(x, -one)),
               And(Lt(-one, x), Lt(x, one)),
               And(Lt(one, x), Lt(x, inf)))


def test_reduce_rational_inequalities_real_relational():
    assert reduce_rational_inequalities([], x) == False
    assert reduce_rational_inequalities(
        [[(x**2 + 3*x + 2)/(x**2 - 16) >= 0]], x, relational=False) == \
        Union(Interval.open(-oo, -4), Interval(-2, -1), Interval.open(4, oo))

    assert reduce_rational_inequalities(
        [[((-2*x - 10)*(3 - x))/((x**2 + 5)*(x - 2)**2) < 0]], x,
        relational=False) == \
        Union(Interval.open(-5, 2), Interval.open(2, 3))

    assert reduce_rational_inequalities([[(x + 1)/(x - 5) <= 0]], x,
        relational=False) == \
        Interval.Ropen(-1, 5)

    assert reduce_rational_inequalities([[(x**2 + 4*x + 3)/(x - 1) > 0]], x,
        relational=False) == \
        Union(Interval.open(-3, -1), Interval.open(1, oo))

    assert reduce_rational_inequalities([[(x**2 - 16)/(x - 1)**2 < 0]], x,
        relational=False) == \
        Union(Interval.open(-4, 1), Interval.open(1, 4))

    assert reduce_rational_inequalities([[(3*x + 1)/(x + 4) >= 1]], x,
        relational=False) == \
        Union(Interval.open(-oo, -4), Interval.Ropen(Rational(3, 2), oo))

    assert reduce_rational_inequalities([[(x - 8)/x <= 3 - x]], x,
        relational=False) == \
        Union(Interval.Lopen(-oo, -2), Interval.Lopen(0, 4))

    # issue sympy/sympy#10237
    assert reduce_rational_inequalities(
        [[x < oo, x >= 0, -oo < x]], x, relational=False) == Interval(0, oo)


def test_reduce_abs_inequalities():
    e = abs(x - 5) < 3
    ans = And(Lt(2, x), Lt(x, 8))
    assert reduce_inequalities(e) == ans
    assert reduce_inequalities(e, x) == ans
    assert reduce_inequalities(abs(x - 5)) == Eq(x, 5)
    assert reduce_inequalities(
        abs(2*x + 3) >= 8) == Or(And(Le(Rational(5, 2), x), Lt(x, oo)),
        And(Le(x, Rational(-11, 2)), Lt(-oo, x)))
    assert reduce_inequalities(abs(x - 4) + abs(
        3*x - 5) < 7) == And(Lt(S.Half, x), Lt(x, 4))
    assert reduce_inequalities(abs(x - 4) + abs(3*abs(x) - 5) < 7) == \
        Or(And(S(-2) < x, x < -1), And(S.Half < x, x < 4))

    nr = Symbol('nr', extended_real=False)
    raises(TypeError, lambda: reduce_inequalities(abs(nr - 5) < 3))
    assert reduce_inequalities(x < 3, symbols=[x, nr]) == And(-oo < x, x < 3)


def test_reduce_inequalities_general():
    assert reduce_inequalities(Ge(sqrt(2)*x, 1)) == And(sqrt(2)/2 <= x, x < oo)
    assert reduce_inequalities(x + 1 > 0) == And(S.NegativeOne < x, x < oo)


def test_reduce_inequalities_boolean():
    assert reduce_inequalities(
        [Eq(x**2, 0), True]) == Eq(x, 0)
    assert reduce_inequalities([Eq(x**2, 0), False]) == False
    assert reduce_inequalities(x**2 >= 0) is S.true  # issue 10196


def test_reduce_inequalities_multivariate():
    assert reduce_inequalities([Ge(x**2, 1), Ge(y**2, 1)]) == And(
        Or(And(Le(S.One, x), Lt(x, oo)), And(Le(x, -1), Lt(-oo, x))),
        Or(And(Le(S.One, y), Lt(y, oo)), And(Le(y, -1), Lt(-oo, y))))


def test_reduce_inequalities_errors():
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(sin(x) + x, 1)))
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(x**2*y + y, 1)))


def test__solve_inequalities():
    assert reduce_inequalities(x + y < 1, symbols=[x]) == (x < 1 - y)
    assert reduce_inequalities(x + y >= 1, symbols=[x]) == (x < oo) & (x >= -y + 1)
    assert reduce_inequalities(Eq(0, x - y), symbols=[x]) == Eq(x, y)
    assert reduce_inequalities(Ne(0, x - y), symbols=[x]) == Ne(x, y)


def test_issue_6343():
    eq = -3*x**2/2 - x*Rational(45, 4) + Rational(33, 2) > 0
    assert reduce_inequalities(eq) == \
        And(x < Rational(-15, 4) + sqrt(401)/4, -sqrt(401)/4 - Rational(15, 4) < x)


def test_issue_8235():
    assert reduce_inequalities(x**2 - 1 < 0) == \
        And(S.NegativeOne < x, x < 1)
    assert reduce_inequalities(x**2 - 1 <= 0) == \
        And(S.NegativeOne <= x, x <= 1)
    assert reduce_inequalities(x**2 - 1 > 0) == \
        Or(And(-oo < x, x < -1), And(x < oo, S.One < x))
    assert reduce_inequalities(x**2 - 1 >= 0) == \
        Or(And(-oo < x, x <= -1), And(S.One <= x, x < oo))

    eq = x**8 + x - 9  # we want CRootOf solns here
    sol = solve(eq >= 0)
    tru = Or(And(rootof(eq, 1) <= x, x < oo), And(-oo < x, x <= rootof(eq, 0)))
    assert sol == tru

    # recast vanilla as real
    assert solve(sqrt((-x + 1)**2) < 1) == And(S.Zero < x, x < 2)


def test_issue_5526():
    assert reduce_inequalities(0 <=
        x + Integral(y**2, (y, 1, 3)) - 1, [x]) == \
        (x >= -Integral(y**2, (y, 1, 3)) + 1)
    f = Function('f')
    e = Sum(f(x), (x, 1, 3))
    assert reduce_inequalities(0 <= x + e + y**2, [x]) == \
        (x >= -y**2 - Sum(f(x), (x, 1, 3)))


def test_solve_univariate_inequality():
    assert isolve(x**2 >= 4, x, relational=False) == Union(Interval(-oo, -2),
        Interval(2, oo))
    assert isolve(x**2 >= 4, x) == Or(And(Le(2, x), Lt(x, oo)), And(Le(x, -2),
        Lt(-oo, x)))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x, relational=False) == \
        Union(Interval(1, 2), Interval(3, oo))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x) == \
        Or(And(Le(1, x), Le(x, 2)), And(Le(3, x), Lt(x, oo)))
    assert isolve((x - 1)*(x - 2)*(x - 4) < 0, x, domain = FiniteSet(0, 3)) == \
        Or(Eq(x, 0), Eq(x, 3))
    # issue 2785:
    assert isolve(x**3 - 2*x - 1 > 0, x, relational=False) == \
        Union(Interval(-1, -sqrt(5)/2 + S.Half, True, True),
              Interval(S.Half + sqrt(5)/2, oo, True, True))
    # issue 2794:
    assert isolve(x**3 - x**2 + x - 1 > 0, x, relational=False) == \
        Interval(1, oo, True)
    #issue 13105
    assert isolve((x + I)*(x + 2*I) < 0, x) == Eq(x, 0)
    assert isolve(((x - 1)*(x - 2) + I)*((x - 1)*(x - 2) + 2*I) < 0, x) == Or(Eq(x, 1), Eq(x, 2))
    assert isolve((((x - 1)*(x - 2) + I)*((x - 1)*(x - 2) + 2*I))/(x - 2) > 0, x) == Eq(x, 1)
    raises (ValueError, lambda: isolve((x**2 - 3*x*I + 2)/x < 0, x))

    # numerical testing in valid() is needed
    assert isolve(x**7 - x - 2 > 0, x) == \
        And(rootof(x**7 - x - 2, 0) < x, x < oo)

    # handle numerator and denominator; although these would be handled as
    # rational inequalities, these test confirm that the right thing is done
    # when the domain is EX (e.g. when 2 is replaced with sqrt(2))
    assert isolve(1/(x - 2) > 0, x) == And(S(2) < x, x < oo)
    den = ((x - 1)*(x - 2)).expand()
    assert isolve((x - 1)/den <= 0, x) == \
        (x > -oo) & (x < 2) & Ne(x, 1)

    n = Dummy('n')
    raises(NotImplementedError, lambda: isolve(Abs(x) <= n, x, relational=False))
    c1 = Dummy("c1", positive=True)
    raises(NotImplementedError, lambda: isolve(n/c1 < 0, c1))
    n = Dummy('n', negative=True)
    assert isolve(n/c1 > -2, c1) == (-n/2 < c1)
    assert isolve(n/c1 < 0, c1) == True
    assert isolve(n/c1 > 0, c1) == False

    zero = cos(1)**2 + sin(1)**2 - 1
    raises(NotImplementedError, lambda: isolve(x**2 < zero, x))
    raises(NotImplementedError, lambda: isolve(
        x**2 < zero*I, x))
    raises(NotImplementedError, lambda: isolve(1/(x - y) < 2, x))
    raises(NotImplementedError, lambda: isolve(1/(x - y) < 0, x))
    raises(TypeError, lambda: isolve(x - I < 0, x))

    zero = x**2 + x - x*(x + 1)
    assert isolve(zero < 0, x, relational=False) is S.EmptySet
    assert isolve(zero <= 0, x, relational=False) is S.Reals

    # make sure iter_solutions gets a default value
    raises(NotImplementedError, lambda: isolve(
        Eq(cos(x)**2 + sin(x)**2, 1), x))


def test_trig_inequalities():
    # all the inequalities are solved in a periodic interval.
    assert isolve(sin(x) < S.Half, x, relational=False) == \
        Union(Interval(0, pi/6, False, True), Interval.open(pi*Rational(5, 6), 2*pi))
    assert isolve(sin(x) > S.Half, x, relational=False) == \
        Interval(pi/6, pi*Rational(5, 6), True, True)
    assert isolve(cos(x) < S.Zero, x, relational=False) == \
        Interval(pi/2, pi*Rational(3, 2), True, True)
    assert isolve(cos(x) >= S.Zero, x, relational=False) == \
        Union(Interval(0, pi/2), Interval.Ropen(pi*Rational(3, 2), 2*pi))

    assert isolve(tan(x) < S.One, x, relational=False) == \
        Union(Interval.Ropen(0, pi/4), Interval.open(pi/2, pi))

    assert isolve(sin(x) <= S.Zero, x, relational=False) == \
        Union(FiniteSet(S.Zero), Interval.Ropen(pi, 2*pi))

    assert isolve(sin(x) <= S.One, x, relational=False) == S.Reals
    assert isolve(cos(x) < S(-2), x, relational=False) == S.EmptySet
    assert isolve(sin(x) >= S.NegativeOne, x, relational=False) == S.Reals
    assert isolve(cos(x) > S.One, x, relational=False) == S.EmptySet


def test_issue_9954():
    assert isolve(x**2 >= 0, x, relational=False) == S.Reals
    assert isolve(x**2 >= 0, x, relational=True) == S.Reals.as_relational(x)
    assert isolve(x**2 < 0, x, relational=False) == S.EmptySet
    assert isolve(x**2 < 0, x, relational=True) == S.EmptySet.as_relational(x)


@XFAIL
def test_slow_general_univariate():
    r = rootof(x**5 - x**2 + 1, 0)
    assert solve(sqrt(x) + 1/root(x, 3) > 1) == \
        Or(And(0 < x, x < r**6), And(r**6 < x, x < oo))


def test_issue_8545():
    eq = 1 - x - abs(1 - x)
    ans = And(Lt(1, x), Lt(x, oo))
    assert reduce_abs_inequality(eq, '<', x) == ans
    eq = 1 - x - sqrt((1 - x)**2)
    assert reduce_inequalities(eq < 0) == ans


def test_issue_8974():
    assert isolve(-oo < x, x) == And(-oo < x, x < oo)
    assert isolve(oo > x, x) == And(-oo < x, x < oo)


def test_issue_10198():
    assert reduce_inequalities(
        -1 + 1/abs(1/x - 1) < 0) == (x > -oo) & (x < S(1)/2) & Ne(x, 0)

    assert reduce_inequalities(abs(1/sqrt(x)) - 1, x) == Eq(x, 1)
    assert reduce_abs_inequality(-3 + 1/abs(1 - 1/x), '<', x) == \
        Or(And(-oo < x, x < 0),
        And(S.Zero < x, x < Rational(3, 4)), And(Rational(3, 2) < x, x < oo))
    raises(ValueError,lambda: reduce_abs_inequality(-3 + 1/abs(
        1 - 1/sqrt(x)), '<', x))


def test_issue_10047():
    # issue 10047: this must remain an inequality, not True, since if x
    # is not real the inequality is invalid
    # assert solve(sin(x) < 2) == (x <= oo)

    # with PR 16956, (x <= oo) autoevaluates when x is extended_real
    # which is assumed in the current implementation of inequality solvers
    assert solve(sin(x) < 2) == True
    assert solveset(sin(x) < 2, domain=S.Reals) == S.Reals


def test_issue_10268():
    assert solve(log(x) < 1000) == And(S.Zero < x, x < exp(1000))


@XFAIL
def test_isolve_Sets():
    n = Dummy('n')
    assert isolve(Abs(x) <= n, x, relational=False) == \
        Piecewise((S.EmptySet, n < 0), (Interval(-n, n), True))


def test_integer_domain_relational_isolve():

    dom = FiniteSet(0, 3)
    x = Symbol('x',zero=False)
    assert isolve((x - 1)*(x - 2)*(x - 4) < 0, x, domain=dom) == Eq(x, 3)

    x = Symbol('x')
    assert isolve(x + 2 < 0, x, domain=S.Integers) == \
           (x <= -3) & (x > -oo) & Eq(Mod(x, 1), 0)
    assert isolve(2 * x + 3 > 0, x, domain=S.Integers) == \
           (x >= -1) & (x < oo)  & Eq(Mod(x, 1), 0)
    assert isolve((x ** 2 + 3 * x - 2) < 0, x, domain=S.Integers) == \
           (x >= -3) & (x <= 0)  & Eq(Mod(x, 1), 0)
    assert isolve((x ** 2 + 3 * x - 2) > 0, x, domain=S.Integers) == \
           ((x >= 1) & (x < oo)  & Eq(Mod(x, 1), 0)) | (
               (x <= -4) & (x > -oo)  & Eq(Mod(x, 1), 0))


def test_issue_10671_12466():
    assert solveset(sin(y), y, Interval(0, pi)) == FiniteSet(0, pi)
    i = Interval(1, 10)
    assert solveset((1/x).diff(x) < 0, x, i) == i
    assert solveset((log(x - 6)/x) <= 0, x, S.Reals) == \
        Interval.Lopen(6, 7)


def test__solve_inequality():
    for op in (Gt, Lt, Le, Ge, Eq, Ne):
        assert _solve_inequality(op(x, 1), x).lhs == x
        assert _solve_inequality(op(S.One, x), x).lhs == x
    # don't get tricked by symbol on right: solve it
    assert _solve_inequality(Eq(2*x - 1, x), x) == Eq(x, 1)
    ie = Eq(S.One, y)
    assert _solve_inequality(ie, x) == ie
    for fx in (x**2, exp(x), sin(x) + cos(x), x*(1 + x)):
        for c in (0, 1):
            e = 2*fx - c > 0
            assert _solve_inequality(e, x, linear=True) == (
                fx > c/S(2))
    assert _solve_inequality(2*x**2 + 2*x - 1 < 0, x, linear=True) == (
        x*(x + 1) < S.Half)
    assert _solve_inequality(Eq(x*y, 1), x) == Eq(x*y, 1)
    nz = Symbol('nz', nonzero=True)
    assert _solve_inequality(Eq(x*nz, 1), x) == Eq(x, 1/nz)
    assert _solve_inequality(x*nz < 1, x) == (x*nz < 1)
    a = Symbol('a', positive=True)
    assert _solve_inequality(a/x > 1, x) == (S.Zero < x) & (x < a)
    assert _solve_inequality(a/x > 1, x, linear=True) == (1/x > 1/a)
    # make sure to include conditions under which solution is valid
    e = Eq(1 - x, x*(1/x - 1))
    assert _solve_inequality(e, x) == Ne(x, 0)
    assert _solve_inequality(x < x*(1/x - 1), x) == (x < S.Half) & Ne(x, 0)


def test__pt():
    from sympy.solvers.inequalities import _pt
    assert _pt(-oo, oo) == 0
    assert _pt(S.One, S(3)) == 2
    assert _pt(S.One, oo) == _pt(oo, S.One) == 2
    assert _pt(S.One, -oo) == _pt(-oo, S.One) == S.Half
    assert _pt(S.NegativeOne, oo) == _pt(oo, S.NegativeOne) == Rational(-1, 2)
    assert _pt(S.NegativeOne, -oo) == _pt(-oo, S.NegativeOne) == -2
    assert _pt(x, oo) == _pt(oo, x) == x + 1
    assert _pt(x, -oo) == _pt(-oo, x) == x - 1
    raises(ValueError, lambda: _pt(Dummy('i', infinite=True), S.One))


def test_lp():
    np = import_module("numpy")
    scipy = import_module("scipy")

    def get_results_with_scipy(objective, constraints, variables):
        if scipy is not None and np is not None:
            from sympy.solvers.inequalities import _np
            nonpos, rep, xx = _np(constraints, [])
            assert not rep  # only testing nonneg variables
            C, _D = linear_eq_to_matrix(objective, *variables)
            A, B = linear_eq_to_matrix(nonpos, *variables)
            assert _D[0] == 0  # scipy only deals with D = 0



            A_sci = Matrix([[A], [-eye(len(variables))]])
            B_sci = Matrix([[B], [Matrix([0] * len(variables))]])
            C_sci = C
            A_sci = np.array(A_sci.tolist())
            B_sci = np.array(B_sci.tolist())
            C_sci = np.array(C_sci.tolist())
            res = scipy.optimize.linprog(C_sci, A_ub=A_sci, b_ub=B_sci)
            return res

    r1 = y+2*z <= 3
    r2 = -x - 3*z <= -2
    r3 = 2*x + y + 7*z <= 5
    constraints = [r1, r2, r3]
    objective = -x - y - 5 * z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints)
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    r1 = x - y + 2*z <= 3
    r2 = -x + 2*y - 3*z <= -2
    r3 = 2*x + y - 7*z <= -5
    constraints = [r1, r2, r3]
    objective = -x - y - 5*z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    r1 = x - y + 2*z <= -4
    r2 = -x + 2*y - 3*z <= 8
    r3 = 2*x + y - 7*z <= 10
    constraints = [r1, r2, r3]
    const = 2
    objective = -x-y-5*z+const # has constant term
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective-const, constraints, variables)
        assert optimum.evalf() == (sympify(-scipy_res.fun)+const)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # Section 4 Problem 1 from
    # http://web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
    # answer on page 55
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    r1 = x1 - x2 - 2*x3 - x4 <= 4
    r2 = 2*x1 + x3 -4*x4 <= 2
    r3 = -2*x1 + x2 + x4 <= 1
    optimum, argmax = lp(max,
        x1 - 2*x2 - 3*x3 - x4, [r1, r2, r3], [])
    assert optimum == 4
    assert list(argmax.values()) == [7, 0, 0, 3]

    # equality
    r1 = Eq(x,y)
    r2 = Eq(y,z)
    r3 = z <= 3
    constraints = [r1, r2, r3]
    objective = x
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # Binary predicate
    r1 = Q.ge(x, y)
    r2 = Q.le(x, 4)
    constraints = [r1, r2]
    objective = x + y
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert optimum == 8
    assert [x.subs(argmax), y.subs(argmax)] == [4, 4]

    # input contains Floats
    r1 = x - y + 2.0*z <= -4
    r2 = -x + 2*y - 3.0*z <= 8
    r3 = 2*x + y - 7*z <= 10
    constraints = [r1, r2, r3]
    objective = -x-y-5*z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # input contains non-float or non-Rational
    r1 = x - y + sqrt(2) * z <= -4
    r2 = -x + 2*y - 3*z <= 8
    r3 = 2*x + y - 7*z <= 10
    raises(TypeError, lambda: lp(max, -x-y-5*z, [r1, r2, r3], []))

    r1 = x >= 0
    raises(UnboundedLinearProgrammingError, lambda: lp(max, x, [r1], []))
    r2 = x <= -1
    raises(InfeasibleLinearProgrammingError, lambda: lp(max, x, [r1,r2], []))

    # strict inequalities are not allowed
    r1 = x > 0
    raises(TypeError, lambda: lp(max, x, [r1], []))

    # not equals not allowed
    r1 = Ne(x, 0)
    raises(TypeError, lambda: lp(max, x, [r1], []))

    def make_random_problem(num_variables=2, num_constraints=2, sparsity=.1):
        def rand():
            if random() < sparsity:
                return sympify(0)
            int1, int2 = [randprime(0, 200) for _ in range(2)]
            return Rational(int1, int2)*choice([-1, 1])
        variables = symbols('x1:%s' % (num_variables + 1))
        constraints = [(sum(rand()*x for x in variables) <= rand())
                       for _ in range(num_constraints)]
        objective = sum(rand() * x for x in variables)
        return objective, constraints, variables

    # testing random problems
    if scipy is not None and np is not None:
        for _ in range(50):
            objective, constraints, variables = make_random_problem()
            constraints = [c for c in constraints if isinstance(c, Relational)] # in case c auto simplifies to True or False
            if len(constraints) == 0:
                continue

            # check lp maximization
            # scipy minimizes, so negative objective for it
            scipy_res = get_results_with_scipy(-objective, constraints, variables)
            if scipy_res.status == 0:
                optimum, argmax = lp(max, objective, constraints, [])
                scipy_op = -scipy_res.fun  # negated to give actual max
                assert abs(optimum.evalf() - scipy_op) < .1**10
                assert objective.subs(argmax) == optimum
                for constr in constraints:
                    assert constr.subs(argmax) == True
            elif scipy_res.status == 2:
                # scipy: problem is infeasible
                raises(InfeasibleLinearProgrammingError,
                       lambda: lp(max, objective, constraints, []))
            elif scipy_res.status == 3:
                # scipy: problem is unbounded
                raises(UnboundedLinearProgrammingError,
                       lambda: lp(max, objective, constraints, []))
            else:
                # scipy: either iteration limit reached or numerical difficulties
                pass

            # check lp minimization
            scipy_res = get_results_with_scipy(objective, constraints, variables)
            if scipy_res.status == 0:
                optimum, argmax = lp(min, objective, constraints, [])
                scipy_op = scipy_res.fun
                assert abs(optimum.evalf() - scipy_op) < .1**10
                assert objective.subs(argmax) == optimum
                for constr in constraints:
                    assert constr.subs(argmax) == True
            elif scipy_res.status == 2:
                # scipy: problem is infeasible
                raises(InfeasibleLinearProgrammingError,
                       lambda: lp(max, objective, constraints, []))
            elif scipy_res.status == 3:
                # scipy: problem is unbounded
                raises(UnboundedLinearProgrammingError,
                       lambda: lp(max, objective, constraints, []))
            else:
                # scipy: either iteration limit reached or numerical difficulties
                pass
