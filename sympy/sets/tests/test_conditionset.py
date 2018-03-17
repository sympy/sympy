from sympy.sets import (ConditionSet, Intersection, FiniteSet,
    EmptySet, Union)
from sympy import (Symbol, Eq, Lt, S, Abs, sin, pi, Lambda, Interval,
    And, Mod)
from sympy.utilities.pytest import raises


w = Symbol('w')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
L = Symbol('lambda')


def test_CondSet():
    sin_sols_principal = ConditionSet(x, Eq(sin(x), 0),
                                      Interval(0, 2*pi, False, True))
    assert pi in sin_sols_principal
    assert pi/2 not in sin_sols_principal
    assert 3*pi not in sin_sols_principal
    assert 5 in ConditionSet(x, x**2 > 4, S.Reals)
    assert 1 not in ConditionSet(x, x**2 > 4, S.Reals)

    assert isinstance(ConditionSet(x, x < 1, {x, y}).base_set, FiniteSet)
    raises(TypeError, lambda: ConditionSet(x, x + 1, {x, y}))
    raises(TypeError, lambda: ConditionSet(x, x, 1))

    I = S.Integers
    C = ConditionSet
    assert C(x, x < 1, C(x, x < 2, I)
        ) == C(x, (x < 1) & (x < 2), I)
    assert C(y, y < 1, C(x, y < 2, I)
        ) == C(x, (x < 1) & (y < 2), I)
    assert C(y, y < 1, C(x, x < 2, I)
        ) == C(y, (y < 1) & (y < 2), I)
    assert C(y, y < 1, C(x, y < x, I)
        ) == C(x, (x < 1) & (y < x), I)
    assert C(y, x < 1, C(x, y < x, I)
        ) == C(L, (x < 1) & (y < L), I)
    c = C(y, x < 1, C(x, L < y, I))
    assert c == C(c.sym, (L < y) & (x < 1), I)
    assert c.sym not in (x, y, L)
    c = C(y, x < 1, C(x, y < x, FiniteSet(L)))
    assert c == C(
        c.sym, c.condition.xreplace({L: c.sym}), FiniteSet(L))
    assert c.sym not in (x, y, L)


def test_CondSet_intersect():
    input_conditionset = ConditionSet(x, x**2 > 4, Interval(1, 4, False, False))
    other_domain = Interval(0, 3, False, False)
    output_conditionset = ConditionSet(x, x**2 > 4, Interval(1, 3, False, False))
    assert Intersection(input_conditionset, other_domain) == output_conditionset


def test_issue_9849():
    assert ConditionSet(x, Eq(x, x), S.Naturals) == S.Naturals
    assert ConditionSet(x, Eq(Abs(sin(x)), -1), S.Naturals) == S.EmptySet


def test_simplified_FiniteSet_in_CondSet():
    assert ConditionSet(x, And(x < 1, x > -3), FiniteSet(0, 1, 2)) == FiniteSet(0)
    assert ConditionSet(x, x < 0, FiniteSet(0, 1, 2)) == EmptySet()
    assert ConditionSet(x, And(x < -3), EmptySet()) == EmptySet()
    y = Symbol('y')
    assert (ConditionSet(x, And(x > 0), FiniteSet(-1, 0, 1, y)) ==
        Union(FiniteSet(1), ConditionSet(x, And(x > 0), FiniteSet(y))))
    assert (ConditionSet(x, Eq(Mod(x, 3), 1), FiniteSet(1, 4, 2, y)) ==
        Union(FiniteSet(1, 4), ConditionSet(x, Eq(Mod(x, 3), 1), FiniteSet(y))))


def test_free_symbols():
    assert ConditionSet(x, Eq(y, 0), FiniteSet(z)
        ).free_symbols == {y, z}
    assert ConditionSet(x, Eq(x, 0), FiniteSet(z)
        ).free_symbols == {z}
    assert ConditionSet(x, Eq(x, 0), FiniteSet(x, z)
        ).free_symbols == {x, z}


def test_subs_CondSet():
    s = FiniteSet(z, y)
    c = ConditionSet(x, x < 2, s)
    # you can only replace sym with a symbol that is not in
    # the free symbols
    assert c.subs(x, 1) == c
    assert c.subs(x, y) == c
    assert c.subs(x, w) == ConditionSet(w, w < 2, s)
    assert ConditionSet(x, x < y, s
        ).subs(y, w) == ConditionSet(x, x < w, s.subs(y, w))

    # to eventually be removed
    c = ConditionSet((x, y), {x + 1, x + y}, S.Reals)
    assert c.subs(x, z) == c


def test_dummy_eq():
    C = ConditionSet
    I = S.Integers
    c = C(x, x < 1, I)
    assert c.dummy_eq(C(y, y < 1, I))
    assert c.dummy_eq(1) == False
    assert c.dummy_eq(C(x, x < 1, S.Reals)) == False
    raises(ValueError, lambda: c.dummy_eq(C(x, x < 1, S.Reals), z))

    # to eventually be removed
    c1 = ConditionSet((x, y), {x + 1, x + y}, S.Reals)
    c2 = ConditionSet((x, y), {x + 1, x + y}, S.Reals)
    c3 = ConditionSet((x, y), {x + 1, x + y}, S.Complexes)
    assert c1.dummy_eq(c2)
    assert c1.dummy_eq(c3) is False
    assert c.dummy_eq(c1) is False
    assert c1.dummy_eq(c) is False
