from __future__ import annotations

from sympy.core.numbers import Rational
from sympy.core.relational import Eq, Ne
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, Xor
from sympy.polys.polyerrors import PolynomialError
from sympy.polys.rootoftools import CRootOf
from sympy.sets.sets import Interval, FiniteSet, Union
from sympy.polys.cad.qe import (quantifier_elimination, decide,
    sample_points, solution_set, _truth_values)
from sympy.testing.pytest import raises
from sympy.abc import a, b, c, d, e, x, y, z

qe = quantifier_elimination


def _equivalent(f, g, gens):
    """Whether two quantifier-free formulas agree on every cell of a
    decomposition sign-invariant for the polynomials of both."""
    _, _, truth = _truth_values(Equivalent(f, g), gens, [], None)
    return all(value for _, value in truth)


def test_decide():
    assert decide(x >= 22, [('forall', x)]) is False
    assert decide(x >= 22, [('exists', x)]) is True
    assert decide(x >= 22, [('exists', [x, y])]) is True
    assert decide(Eq(x, 2*y), [('forall', x), ('exists', y)]) is True
    assert decide(Eq(x, y**2), [('forall', x), ('exists', y)]) is False
    assert decide(Eq(x, y**2), [('exists', x), ('forall', y)]) is False
    assert decide(Eq(y, x**2), [('forall', x), ('exists', y)]) is True
    assert decide(Ne(x, x + 1), [('forall', x)]) is True
    assert decide(Ne(x, 2*x), [('forall', x)]) is False
    assert decide(Eq(x, 2*x), [('exists', x)]) is True
    assert decide(And(x + 1 < y, Eq(x, y)), [('exists', [x, y])]) is False
    assert decide(Or(x >= 22, x < 22), [('forall', x)]) is True
    assert decide(And(x >= 5, x < 1), [('forall', x)]) is False
    assert decide(Implies(x > 0, x >= 0), [('forall', x)]) is True
    assert decide(Implies(x > 0, x < 0), [('forall', x)]) is False
    assert decide(Equivalent(x >= 0, x**3 >= 0), [('forall', x)]) is True
    assert decide(Xor(x > 0, x <= 0), [('forall', x)]) is True
    assert decide(Not(x**2 < 0), [('forall', x)]) is True
    assert decide(Eq(y, z*x), [('forall', [x, y]), ('exists', z)]) is False
    assert decide(Eq(y, z*x), [('forall', y), ('exists', [x, z])]) is True
    assert decide(x**2 + y**2 < 0, [('exists', [x, y])]) is False
    assert decide(x**2 + y**2 < 1, [('forall', x), ('exists', y)]) is False
    assert decide(x**2 + y**2 < 1, [('exists', x), ('exists', y)]) is True
    assert decide(x**2 + y**2 - 2*x*y >= 0, [('forall', [x, y])]) is True
    assert decide(S.true, [('forall', x)]) is True
    assert decide(S.false, [('exists', x)]) is False
    # the quantifier order matters
    assert decide(x < y, [('forall', x), ('exists', y)]) is True
    assert decide(x < y, [('exists', y), ('forall', x)]) is False
    # a real root of every odd degree polynomial, none for x**2 + 1
    assert decide(Eq(x**3 + a*x + b, 0), [('forall', [a, b]), ('exists', x)]) is True
    assert decide(Eq(x**2 + 1, 0), [('exists', x)]) is False
    # positive definite quadratic form
    assert decide(x**2 + x*y + y**2 > 0, [('forall', [x, y])]) is False
    assert decide(Or(x**2 + x*y + y**2 > 0, And(Eq(x, 0), Eq(y, 0))), [('forall', [x, y])]) is True
    assert decide(x**2 - x*y + y**2 >= 0, [('forall', [x, y])]) is True
    assert decide(x**2 - 3*x*y + y**2 >= 0, [('forall', [x, y])]) is False


def test_quantifier_elimination_one_variable():
    # solvability of the monic quadratic
    assert qe(Eq(x**2 + a*x + b, 0), [('exists', x)]) == (a**2 - 4*b >= 0)
    # positivity of the monic quadratic
    assert qe(x**2 + b*x + c > 0, [('forall', x)]) == (b**2 - 4*c < 0)
    assert qe(x**2 + b*x + c >= 0, [('forall', x)]) == (b**2 - 4*c <= 0)
    # projections of an ellipse (the answers 8*x**2 - 8*x - 29 <= 0 and
    # 8*y**2 + 16*y - 85 <= 0 are due to QEPCAD)
    ellipse = Eq(3*x**2 + 2*x*y + y**2 - x + y - 7, 0)
    r0, r1 = CRootOf(8*x**2 - 8*x - 29, 0), CRootOf(8*x**2 - 8*x - 29, 1)
    assert qe(ellipse, [('exists', y)]) == And(x >= r0, x <= r1)
    assert solution_set(ellipse, x, [('exists', y)]) == Interval(r0, r1)
    assert solution_set(8*x**2 - 8*x - 29 <= 0, x) == Interval(r0, r1)
    s0, s1 = CRootOf(8*x**2 + 16*x - 85, 0), CRootOf(8*x**2 + 16*x - 85, 1)
    assert solution_set(ellipse, y, [('exists', x)]) == Interval(s0, s1)
    assert solution_set(8*y**2 + 16*y - 85 <= 0, y) == Interval(s0, s1)
    # a circle and its inside
    assert qe(Eq(x**2 + y**2, 1), [('exists', y)]) == And(x >= -1, x <= 1)
    assert qe(x**2 + y**2 < 1, [('exists', y)]) == And(x > -1, x < 1)
    assert qe(x**2 + y**2 < 1, [('forall', y)]) == S.false
    assert qe(x**2 + y**2 > -1, [('forall', y)]) == S.true
    assert qe(Or(x**2 + y**2 < 1, x**2 + y**2 > 4), [('forall', y)]) == Or(x < -2, x > 2)
    # unions of intervals and points
    assert solution_set(x**2 > 2, x) == Union(
        Interval.open(-S.Infinity, CRootOf(x**2 - 2, 0)),
        Interval.open(CRootOf(x**2 - 2, 1), S.Infinity))
    assert solution_set(x**2 >= 2, x) == Union(
        Interval(-S.Infinity, CRootOf(x**2 - 2, 0)),
        Interval(CRootOf(x**2 - 2, 1), S.Infinity))
    assert solution_set(Eq(x**2, 2), x) == FiniteSet(CRootOf(x**2 - 2, 0), CRootOf(x**2 - 2, 1))
    assert solution_set(Or(Eq(x, 1), And(x > 2, x <= 3)), x) == Union(FiniteSet(1), Interval.Lopen(2, 3))
    assert solution_set(x**2 < 0, x) == S.EmptySet
    assert solution_set(x**2 >= 0, x) == S.Reals
    assert solution_set(Ne(x, 1), x) == Union(Interval.open(-S.Infinity, 1), Interval.open(1, S.Infinity))
    assert solution_set(And(x >= 1, x <= 1), x) == FiniteSet(1)
    assert solution_set(x**3 - x > 0, x) == Union(Interval.open(-1, 0), Interval.open(1, S.Infinity))
    assert qe(x**3 - x > 0, []) == Or(And(x > -1, x < 0), x > 1)
    # for which x is there a y with y > x on the unit circle
    assert solution_set(And(Eq(x**2 + y**2, 1), y > x), x, [('exists', y)]) == \
        Interval.Ropen(-1, CRootOf(2*x**2 - 1, 1))
    # the quantified variables need not appear
    assert qe(x > 1, [('forall', y)]) == (x > 1)
    assert solution_set(x > 1, x, [('exists', [y, z])]) == Interval.open(1, S.Infinity)
    # rational coefficients
    assert solution_set(x/2 > Rational(1, 3), x) == Interval.open(Rational(2, 3), S.Infinity)


def test_quantifier_elimination_several_variables():
    # a linear polynomial has a positive value: a /= 0 \/ b > 0
    r = qe(a*x + b > 0, [('exists', x)])
    assert _equivalent(r, Or(Ne(a, 0), b > 0), [a, b])
    # solvability of the general quadratic, QEPCAD gives
    # 4 a c - b^2 <= 0 /\ [ c = 0 \/ a /= 0 \/ 4 a c - b^2 < 0 ]
    r = qe(Eq(a*x**2 + b*x + c, 0), [('exists', x)])
    ref = And(4*a*c - b**2 <= 0, Or(Eq(c, 0), Ne(a, 0), 4*a*c - b**2 < 0))
    assert _equivalent(r, ref, [a, b, c])
    assert r == Or(Eq(c, 0), 4*a*c - b**2 < 0, And(a > 0, Eq(4*a*c - b**2, 0)),
                   And(a < 0, 4*a*c - b**2 <= 0))
    # no real root, QEPCAD: 4 a c - b^2 >= 0 /\ c /= 0 /\ [ b = 0 \/ 4 a c - b^2 > 0 ]
    r = qe(Ne(a*x**2 + b*x + c, 0), [('forall', x)])
    ref = And(4*a*c - b**2 >= 0, Ne(c, 0), Or(Eq(b, 0), 4*a*c - b**2 > 0))
    assert _equivalent(r, ref, [a, b, c])
    assert r == Or(And(Eq(a, 0), Eq(b, 0), Ne(c, 0)), 4*a*c - b**2 > 0)
    # a positive value of the general quadratic; Redlog answers
    # a > 0 or (b < 0 and a = 0) or (a = 0 and (b > 0 or (c > 0 and b = 0)))
    # or (a < 0 and 4*a*c - b^2 < 0)
    r = qe(a*x**2 + b*x + c > 0, [('exists', x)])
    assert r == Or(a > 0, c > 0, 4*a*c - b**2 < 0)
    ref = Or(a > 0, And(b < 0, Eq(a, 0)), And(Eq(a, 0), Or(b > 0, And(c > 0, Eq(b, 0)))),
             And(a < 0, 4*a*c - b**2 < 0))
    assert _equivalent(r, ref, [a, b, c])
    # with a /= 0 assumed the condition is the discriminant
    r = qe(And(Ne(a, 0), Eq(a*x**2 + b*x + c, 0)), [('exists', x)])
    assert _equivalent(r, And(Ne(a, 0), 4*a*c - b**2 <= 0), [a, b, c])
    # a quantifier-free formula is described by its own polynomials
    assert qe((x**2 + y**2 < 1) & (x > y), []) == And(x - y > 0, x**2 + y**2 - 1 < 0)
    assert qe(x**2 + y**2 < 0, []) == S.false
    assert qe(x**2 + y**2 >= 0, []) == S.true
    assert qe(Or(x > y, x <= y), [], free=[x, y]) == S.true
    # existence of a point of the circle over a
    r = qe(And(Eq(x**2 + y**2, 1), Eq(a, x)), [('exists', [x, y])])
    assert r == And(a >= -1, a <= 1)
    # x**2 + a*x + 1 = 0 has a real root iff |a| >= 2
    assert qe(Eq(x**2 + a*x + 1, 0), [('exists', x)]) == Or(a >= 2, a <= -2)
    # x**4 + d*x + e = 0 has a real root iff 256 e^3 <= 27 d^4
    r = qe(Eq(x**4 + d*x + e, 0), [('exists', x)])
    assert _equivalent(r, 256*e**3 - 27*d**4 <= 0, [d, e])
    # the depressed cubic has a real root for every p, q
    p, q = symbols('p q')
    assert qe(Eq(x**3 + p*x + q, 0), [('exists', x)], free=[p, q]) == S.true
    # the projection factors of two free variables may not suffice: for
    # x >= 0 the condition is y < sqrt(x), and no factor vanishes there
    raises(NotImplementedError, lambda: qe(Eq(z**2, x) & (z > y), [('exists', z)], free=[x, y]))


def test_sample_points():
    assert sample_points(x**2 + y**2 < 0, [x, y]) == []
    assert sample_points((x**2 + y**2 < 1) & (x > y), [x, y]) == [
        {x: 0, y: -Rational(1, 2)}, {x: CRootOf(2*x**2 - 1, 1), y: 0}, {x: Rational(3, 4), y: 0}]
    pts = sample_points(Eq(x**2 + y**2, 1), [x, y])
    assert len(pts) == 4 and all((p[x]**2 + p[y]**2 - 1).equals(0) for p in pts)
    assert sample_points(x > 2, [x]) == [{x: 3}]
    assert sample_points(S.true, [x, y]) == [{x: 0, y: 0}]
    assert sample_points((x > 0) & (y > 0) & (x + y < 1), [x, y]) == [{x: Rational(1, 2), y: Rational(1, 3)}]
    f = And(x**2 + y**2 - 7 > 0, x + y - 2 > 0, x**2 + y - 10 > 0)
    pts = sample_points(f, [x, y])
    assert len(pts) == 13 and all(f.subs(pt) for pt in pts)
    assert sample_points(And(x**2 + y**2 - 7 > 0, x + y - 2 > 0, x**2 + y - 10 > 0,
                             x**2 + y**2 - 7 < 0), [x, y]) == []


def test_errors():
    raises(ValueError, lambda: qe(x > 0, [('some', x)]))
    raises(ValueError, lambda: qe(x > 0, [('exists', x)], free=[x]))
    raises(ValueError, lambda: qe(x + y > 0, [('exists', x)], free=[]))
    raises(ValueError, lambda: qe(x + y > 0, [('exists', x)], free=[z]))
    raises(PolynomialError, lambda: qe(sqrt(x) > 0, [('exists', x)]))
    raises(ValueError, lambda: qe(x, [('exists', x)]))
    from sympy.logic.boolalg import ITE
    raises(ValueError, lambda: qe(ITE(x > 0, y > 0, y < 0), [('exists', [x, y])]))


def test_sign_formula_minimization():
    # the disc, described by one factor
    assert qe(x**2 + y**2 < 1, [], free=[x, y]) == (x**2 + y**2 - 1 < 0)
    # a half plane with a redundant polynomial in the formula
    assert qe(Or(x - y > 0, And(x - y > 0, x > 0)), [], free=[x, y]) == (x - y > 0)
    # complement of a line
    assert qe(Not(Eq(x, y)), [], free=[x, y]) == Ne(x - y, 0)
    # a formula that needs the sign of two factors
    assert qe(x*y > 0, [], free=[x, y]) == Or(And(x < 0, y < 0), And(x > 0, y > 0))
    r = qe(x*y >= 0, [], free=[x, y])
    assert r == Or(Eq(y, 0), And(x >= 0, y > 0), And(x <= 0, y <= 0))
    assert _equivalent(r, x*y >= 0, [x, y])
    # merged signs
    r = qe(Ne(x*y, 0), [], free=[x, y])
    assert _equivalent(r, And(Ne(x, 0), Ne(y, 0)), [x, y])
    assert r.count(Ne) == 2 or r == And(Ne(x, 0), Ne(y, 0))
