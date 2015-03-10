from sympy.utilities.pytest import XFAIL, raises
from sympy import (S, Symbol, symbols, nan, oo, I, pi, Float, And, Or, Not,
                   Implies, Xor, zoo, sqrt, Rational, simplify, Function)
from sympy.core.compatibility import range
from sympy.core.relational import (Relational, Equality, Unequality,
                                   GreaterThan, LessThan, StrictGreaterThan,
                                   StrictLessThan, Rel, Eq, Lt, Le,
                                   Gt, Ge, Ne)
from sympy.sets.sets import Interval, FiniteSet

x, y, z, t = symbols('x,y,z,t')


def test_rel_ne():
    assert Relational(x, y, '!=') == Ne(x, y)


def test_rel_subs():
    e = Relational(x, y, '==')
    e = e.subs(x, z)

    assert isinstance(e, Equality)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '>=')
    e = e.subs(x, z)

    assert isinstance(e, GreaterThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '<=')
    e = e.subs(x, z)

    assert isinstance(e, LessThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '>')
    e = e.subs(x, z)

    assert isinstance(e, StrictGreaterThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '<')
    e = e.subs(x, z)

    assert isinstance(e, StrictLessThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Eq(x, 0)
    assert e.subs(x, 0) is S.true
    assert e.subs(x, 1) is S.false


def test_wrappers():
    e = x + x**2

    res = Relational(y, e, '==')
    assert Rel(y, x + x**2, '==') == res
    assert Eq(y, x + x**2) == res

    res = Relational(y, e, '<')
    assert Lt(y, x + x**2) == res

    res = Relational(y, e, '<=')
    assert Le(y, x + x**2) == res

    res = Relational(y, e, '>')
    assert Gt(y, x + x**2) == res

    res = Relational(y, e, '>=')
    assert Ge(y, x + x**2) == res

    res = Relational(y, e, '!=')
    assert Ne(y, x + x**2) == res


def test_Eq():
    assert Eq(x**2) == Eq(x**2, 0)
    assert Eq(x**2) != Eq(x**2, 1)

    assert Eq(x, x)  # issue 5719

    # issue 6116
    p = Symbol('p', positive=True)
    assert Eq(p, 0) is S.false


def test_rel_Infinity():
    # NOTE: All of these are actually handled by sympy.core.Number, and do
    # not create Relational objects.
    assert (oo > oo) is S.false
    assert (oo > -oo) is S.true
    assert (oo > 1) is S.true
    assert (oo < oo) is S.false
    assert (oo < -oo) is S.false
    assert (oo < 1) is S.false
    assert (oo >= oo) is S.true
    assert (oo >= -oo) is S.true
    assert (oo >= 1) is S.true
    assert (oo <= oo) is S.true
    assert (oo <= -oo) is S.false
    assert (oo <= 1) is S.false
    assert (-oo > oo) is S.false
    assert (-oo > -oo) is S.false
    assert (-oo > 1) is S.false
    assert (-oo < oo) is S.true
    assert (-oo < -oo) is S.false
    assert (-oo < 1) is S.true
    assert (-oo >= oo) is S.false
    assert (-oo >= -oo) is S.true
    assert (-oo >= 1) is S.false
    assert (-oo <= oo) is S.true
    assert (-oo <= -oo) is S.true
    assert (-oo <= 1) is S.true


def test_bool():
    assert Eq(0, 0) is S.true
    assert Eq(1, 0) is S.false
    assert Ne(0, 0) is S.false
    assert Ne(1, 0) is S.true
    assert Lt(0, 1) is S.true
    assert Lt(1, 0) is S.false
    assert Le(0, 1) is S.true
    assert Le(1, 0) is S.false
    assert Le(0, 0) is S.true
    assert Gt(1, 0) is S.true
    assert Gt(0, 1) is S.false
    assert Ge(1, 0) is S.true
    assert Ge(0, 1) is S.false
    assert Ge(1, 1) is S.true
    assert Eq(I, 2) is S.false
    assert Ne(I, 2) is S.true
    raises(TypeError, lambda: Gt(I, 2))
    raises(TypeError, lambda: Ge(I, 2))
    raises(TypeError, lambda: Lt(I, 2))
    raises(TypeError, lambda: Le(I, 2))
    a = Float('.000000000000000000001', '')
    b = Float('.0000000000000000000001', '')
    assert Eq(pi + a, pi + b) is S.false


def test_rich_cmp():
    assert (x < y) == Lt(x, y)
    assert (x <= y) == Le(x, y)
    assert (x > y) == Gt(x, y)
    assert (x >= y) == Ge(x, y)


def test_doit():
    from sympy import Symbol
    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    np = Symbol('np', nonpositive=True)
    nn = Symbol('nn', nonnegative=True)

    assert Gt(p, 0).doit() is S.true
    assert Gt(p, 1).doit() == Gt(p, 1)
    assert Ge(p, 0).doit() is S.true
    assert Le(p, 0).doit() is S.false
    assert Lt(n, 0).doit() is S.true
    assert Le(np, 0).doit() is S.true
    assert Gt(nn, 0).doit() == Gt(nn, 0)
    assert Lt(nn, 0).doit() is S.false

    assert Eq(x, 0).doit() == Eq(x, 0)


def test_new_relational():
    x = Symbol('x')

    assert Eq(x) == Relational(x, 0)       # None ==> Equality
    assert Eq(x) == Relational(x, 0, '==')
    assert Eq(x) == Relational(x, 0, 'eq')
    assert Eq(x) == Equality(x, 0)
    assert Eq(x, -1) == Relational(x, -1)       # None ==> Equality
    assert Eq(x, -1) == Relational(x, -1, '==')
    assert Eq(x, -1) == Relational(x, -1, 'eq')
    assert Eq(x, -1) == Equality(x, -1)
    assert Eq(x) != Relational(x, 1)       # None ==> Equality
    assert Eq(x) != Relational(x, 1, '==')
    assert Eq(x) != Relational(x, 1, 'eq')
    assert Eq(x) != Equality(x, 1)
    assert Eq(x, -1) != Relational(x, 1)       # None ==> Equality
    assert Eq(x, -1) != Relational(x, 1, '==')
    assert Eq(x, -1) != Relational(x, 1, 'eq')
    assert Eq(x, -1) != Equality(x, 1)

    assert Ne(x, 0) == Relational(x, 0, '!=')
    assert Ne(x, 0) == Relational(x, 0, '<>')
    assert Ne(x, 0) == Relational(x, 0, 'ne')
    assert Ne(x, 0) == Unequality(x, 0)
    assert Ne(x, 0) != Relational(x, 1, '!=')
    assert Ne(x, 0) != Relational(x, 1, '<>')
    assert Ne(x, 0) != Relational(x, 1, 'ne')
    assert Ne(x, 0) != Unequality(x, 1)

    assert Ge(x, 0) == Relational(x, 0, '>=')
    assert Ge(x, 0) == Relational(x, 0, 'ge')
    assert Ge(x, 0) == GreaterThan(x, 0)
    assert Ge(x, 1) != Relational(x, 0, '>=')
    assert Ge(x, 1) != Relational(x, 0, 'ge')
    assert Ge(x, 1) != GreaterThan(x, 0)
    assert (x >= 1) == Relational(x, 1, '>=')
    assert (x >= 1) == Relational(x, 1, 'ge')
    assert (x >= 1) == GreaterThan(x, 1)
    assert (x >= 0) != Relational(x, 1, '>=')
    assert (x >= 0) != Relational(x, 1, 'ge')
    assert (x >= 0) != GreaterThan(x, 1)

    assert Le(x, 0) == Relational(x, 0, '<=')
    assert Le(x, 0) == Relational(x, 0, 'le')
    assert Le(x, 0) == LessThan(x, 0)
    assert Le(x, 1) != Relational(x, 0, '<=')
    assert Le(x, 1) != Relational(x, 0, 'le')
    assert Le(x, 1) != LessThan(x, 0)
    assert (x <= 1) == Relational(x, 1, '<=')
    assert (x <= 1) == Relational(x, 1, 'le')
    assert (x <= 1) == LessThan(x, 1)
    assert (x <= 0) != Relational(x, 1, '<=')
    assert (x <= 0) != Relational(x, 1, 'le')
    assert (x <= 0) != LessThan(x, 1)

    assert Gt(x, 0) == Relational(x, 0, '>')
    assert Gt(x, 0) == Relational(x, 0, 'gt')
    assert Gt(x, 0) == StrictGreaterThan(x, 0)
    assert Gt(x, 1) != Relational(x, 0, '>')
    assert Gt(x, 1) != Relational(x, 0, 'gt')
    assert Gt(x, 1) != StrictGreaterThan(x, 0)
    assert (x > 1) == Relational(x, 1, '>')
    assert (x > 1) == Relational(x, 1, 'gt')
    assert (x > 1) == StrictGreaterThan(x, 1)
    assert (x > 0) != Relational(x, 1, '>')
    assert (x > 0) != Relational(x, 1, 'gt')
    assert (x > 0) != StrictGreaterThan(x, 1)

    assert Lt(x, 0) == Relational(x, 0, '<')
    assert Lt(x, 0) == Relational(x, 0, 'lt')
    assert Lt(x, 0) == StrictLessThan(x, 0)
    assert Lt(x, 1) != Relational(x, 0, '<')
    assert Lt(x, 1) != Relational(x, 0, 'lt')
    assert Lt(x, 1) != StrictLessThan(x, 0)
    assert (x < 1) == Relational(x, 1, '<')
    assert (x < 1) == Relational(x, 1, 'lt')
    assert (x < 1) == StrictLessThan(x, 1)
    assert (x < 0) != Relational(x, 1, '<')
    assert (x < 0) != Relational(x, 1, 'lt')
    assert (x < 0) != StrictLessThan(x, 1)

    # finally, some fuzz testing
    from random import randint
    from sympy.core.compatibility import unichr
    for i in range(100):
        while 1:
            strtype, length = (unichr, 65535) if randint(0, 1) else (chr, 255)
            relation_type = strtype(randint(0, length))
            if randint(0, 1):
                relation_type += strtype(randint(0, length))
            if relation_type not in ('==', 'eq', '!=', '<>', 'ne', '>=', 'ge',
                                     '<=', 'le', '>', 'gt', '<', 'lt'):
                break

        raises(ValueError, lambda: Relational(x, 1, relation_type))


def test_relational_bool_output():
    # https://github.com/sympy/sympy/issues/5931
    raises(TypeError, lambda: bool(x > 3))
    raises(TypeError, lambda: bool(x >= 3))
    raises(TypeError, lambda: bool(x < 3))
    raises(TypeError, lambda: bool(x <= 3))
    raises(TypeError, lambda: bool(Eq(x, 3)))
    raises(TypeError, lambda: bool(Ne(x, 3)))


def test_relational_logic_symbols():
    # See issue 6204
    assert (x < y) & (z < t) == And(x < y, z < t)
    assert (x < y) | (z < t) == Or(x < y, z < t)
    assert ~(x < y) == Not(x < y)
    assert (x < y) >> (z < t) == Implies(x < y, z < t)
    assert (x < y) << (z < t) == Implies(z < t, x < y)
    assert (x < y) ^ (z < t) == Xor(x < y, z < t)

    assert isinstance((x < y) & (z < t), And)
    assert isinstance((x < y) | (z < t), Or)
    assert isinstance(~(x < y), GreaterThan)
    assert isinstance((x < y) >> (z < t), Implies)
    assert isinstance((x < y) << (z < t), Implies)
    assert isinstance((x < y) ^ (z < t), (Or, Xor))


def test_univariate_relational_as_set():
    assert (x > 0).as_set() == Interval(0, oo, True, True)
    assert (x >= 0).as_set() == Interval(0, oo)
    assert (x < 0).as_set() == Interval(-oo, 0, True, True)
    assert (x <= 0).as_set() == Interval(-oo, 0)
    assert Eq(x, 0).as_set() == FiniteSet(0)
    assert Ne(x, 0).as_set() == Interval(-oo, 0, True, True) + \
        Interval(0, oo, True, True)

    assert (x**2 >= 4).as_set() == Interval(-oo, -2) + Interval(2, oo)


@XFAIL
def test_multivariate_relational_as_set():
    assert (x*y >= 0).as_set() == Interval(0, oo)*Interval(0, oo) + \
        Interval(-oo, 0)*Interval(-oo, 0)


def test_Not():
    assert Not(Equality(x, y)) == Unequality(x, y)
    assert Not(Unequality(x, y)) == Equality(x, y)
    assert Not(StrictGreaterThan(x, y)) == LessThan(x, y)
    assert Not(StrictLessThan(x, y)) == GreaterThan(x, y)
    assert Not(GreaterThan(x, y)) == StrictLessThan(x, y)
    assert Not(LessThan(x, y)) == StrictGreaterThan(x, y)


def test_evaluate():
    assert str(Eq(x, x, evaluate=False)) == 'Eq(x, x)'
    assert Eq(x, x, evaluate=False).doit() == S.true
    assert str(Ne(x, x, evaluate=False)) == 'Ne(x, x)'
    assert Ne(x, x, evaluate=False).doit() == S.false

    assert str(Ge(x, x, evaluate=False)) == 'x >= x'
    assert str(Le(x, x, evaluate=False)) == 'x <= x'
    assert str(Gt(x, x, evaluate=False)) == 'x > x'
    assert str(Lt(x, x, evaluate=False)) == 'x < x'


def assert_all_ineq_raise_TypeError(a, b):
    raises(TypeError, lambda: a > b)
    raises(TypeError, lambda: a >= b)
    raises(TypeError, lambda: a < b)
    raises(TypeError, lambda: a <= b)
    raises(TypeError, lambda: b > a)
    raises(TypeError, lambda: b >= a)
    raises(TypeError, lambda: b < a)
    raises(TypeError, lambda: b <= a)


def assert_all_ineq_give_class_Inequality(a, b):
    """All inequality operations on `a` and `b` result in class Inequality."""
    from sympy.core.relational import _Inequality as Inequality
    assert isinstance(a > b,  Inequality)
    assert isinstance(a >= b, Inequality)
    assert isinstance(a < b,  Inequality)
    assert isinstance(a <= b, Inequality)
    assert isinstance(b > a,  Inequality)
    assert isinstance(b >= a, Inequality)
    assert isinstance(b < a,  Inequality)
    assert isinstance(b <= a, Inequality)


def test_imaginary_compare_raises_TypeError():
    # See issue #5724
    assert_all_ineq_raise_TypeError(I, x)


def test_complex_compare_not_real():
    # two cases which are not real
    y = Symbol('y', imaginary=True)
    z = Symbol('z', complex=True, real=False)
    for w in (y, z):
        assert_all_ineq_raise_TypeError(2, w)
    # some cases which should remain un-evaluated
    t = Symbol('t')
    x = Symbol('x', real=True)
    z = Symbol('z', complex=True)
    for w in (x, z, t):
        assert_all_ineq_give_class_Inequality(2, w)


def test_imaginary_and_inf_compare_raises_TypeError():
    # See pull request #7835
    y = Symbol('y', imaginary=True)
    assert_all_ineq_raise_TypeError(oo, y)
    assert_all_ineq_raise_TypeError(-oo, y)


def test_complex_pure_imag_not_ordered():
    raises(TypeError, lambda: 2*I < 3*I)

    # more generally
    x = Symbol('x', real=True, nonzero=True)
    y = Symbol('y', imaginary=True)
    z = Symbol('z', complex=True)
    assert_all_ineq_raise_TypeError(I, y)

    t = I*x   # an imaginary number, should raise errors
    assert_all_ineq_raise_TypeError(2, t)

    t = -I*y   # a real number, so no errors
    assert_all_ineq_give_class_Inequality(2, t)

    t = I*z   # unknown, should be unevaluated
    assert_all_ineq_give_class_Inequality(2, t)


def test_x_minus_y_not_same_as_x_lt_y():
    """
    A consequence of pull request #7792 is that `x - y < 0` and `x < y`
    are not synonymous.
    """
    x = I + 2
    y = I + 3
    raises(TypeError, lambda: x < y)
    assert x - y < 0

    ineq = Lt(x, y, evaluate=False)
    raises(TypeError, lambda: ineq.doit())
    assert ineq.lhs - ineq.rhs < 0

    t = Symbol('t', imaginary=True)
    x = 2 + t
    y = 3 + t
    ineq = Lt(x, y, evaluate=False)
    raises(TypeError, lambda: ineq.doit())
    assert ineq.lhs - ineq.rhs < 0

    # this one should give error either way
    x = I + 2
    y = 2*I + 3
    raises(TypeError, lambda: x < y)
    raises(TypeError, lambda: x - y < 0)


def test_nan_equality_exceptions():
    # See issue #7774
    import random
    assert Equality(nan, nan) is S.false
    assert Unequality(nan, nan) is S.true

    # See issue #7773
    A = (x, S(0), S(1)/3, pi, oo, -oo)
    assert Equality(nan, random.choice(A)) is S.false
    assert Equality(random.choice(A), nan) is S.false
    assert Unequality(nan, random.choice(A)) is S.true
    assert Unequality(random.choice(A), nan) is S.true


def test_nan_inequality_raise_errors():
    # See discussion in pull request #7776.  We test inequalities with
    # a set including examples of various classes.
    for q in (x, S(0), S(10), S(1)/3, pi, S(1.3), oo, -oo, nan):
        assert_all_ineq_raise_TypeError(q, nan)


def test_nan_complex_inequalities():
    # Comparisons of NaN with non-real raise errors, we're not too
    # fussy whether its the NaN error or complex error.
    for r in (I, zoo, Symbol('z', imaginary=True)):
        assert_all_ineq_raise_TypeError(r, nan)


def test_complex_infinity_inequalities():
    raises(TypeError, lambda: zoo > 0)
    raises(TypeError, lambda: zoo >= 0)
    raises(TypeError, lambda: zoo < 0)
    raises(TypeError, lambda: zoo <= 0)


def test_inequalities_symbol_name_same():
    """Using the operator and functional forms should give same results."""
    # We test all combinations from a set
    # FIXME: could replace with random selection after test passes
    A = (x, y, S(0), S(1)/3, pi, oo, -oo)
    for a in A:
        for b in A:
            assert Gt(a, b) == (a > b)
            assert Lt(a, b) == (a < b)
            assert Ge(a, b) == (a >= b)
            assert Le(a, b) == (a <= b)

    for b in (y, S(0), S(1)/3, pi, oo, -oo):
        assert Gt(x, b, evaluate=False) == (x > b)
        assert Lt(x, b, evaluate=False) == (x < b)
        assert Ge(x, b, evaluate=False) == (x >= b)
        assert Le(x, b, evaluate=False) == (x <= b)

    for b in (y, S(0), S(1)/3, pi, oo, -oo):
        assert Gt(b, x, evaluate=False) == (b > x)
        assert Lt(b, x, evaluate=False) == (b < x)
        assert Ge(b, x, evaluate=False) == (b >= x)
        assert Le(b, x, evaluate=False) == (b <= x)


def test_inequalities_symbol_name_same_complex():
    """Using the operator and functional forms should give same results.
    With complex non-real numbers, both should raise errors.
    """
    # FIXME: could replace with random selection after test passes
    for a in (x, S(0), S(1)/3, pi, oo):
        raises(TypeError, lambda: Gt(a, I))
        raises(TypeError, lambda: a > I)
        raises(TypeError, lambda: Lt(a, I))
        raises(TypeError, lambda: a < I)
        raises(TypeError, lambda: Ge(a, I))
        raises(TypeError, lambda: a >= I)
        raises(TypeError, lambda: Le(a, I))
        raises(TypeError, lambda: a <= I)


def test_inequalities_cant_sympify_other():
    # see issue 7833
    from operator import gt, lt, ge, le

    bar = "foo"

    for a in (x, S(0), S(1)/3, pi, I, zoo, oo, -oo, nan):
        for op in (lt, gt, le, ge):
            raises(TypeError, lambda: op(a, bar))


def test_ineq_avoid_wild_symbol_flip():
    # see issue #7951, we try to avoid this internally, e.g., by using
    # __lt__ instead of "<".
    from sympy.core.symbol import Wild
    p = symbols('p', cls=Wild)
    # x > p might flip, but Gt should not:
    assert Gt(x, p) == Gt(x, p, evaluate=False)
    # Previously failed as 'p > x':
    e = Lt(x, y).subs({y: p})
    assert e == Lt(x, p, evaluate=False)
    # Previously failed as 'p <= x':
    e = Ge(x, p).doit()
    assert e == Ge(x, p, evaluate=False)


def test_issue_8245():
    a = S("6506833320952669167898688709329/5070602400912917605986812821504")
    q = a.n(10)
    assert (a == q) is True
    assert (a != q) is False
    assert (a > q) == False
    assert (a < q) == False
    assert (a >= q) == True
    assert (a <= q) == True

    a = sqrt(2)
    r = Rational(str(a.n(30)))
    assert (r == a) is False
    assert (r != a) is True
    assert (r > a) == True
    assert (r < a) == False
    assert (r >= a) == True
    assert (r <= a) == False
    a = sqrt(2)
    r = Rational(str(a.n(29)))
    assert (r == a) is False
    assert (r != a) is True
    assert (r > a) == False
    assert (r < a) == True
    assert (r >= a) == False
    assert (r <= a) == True


def test_issue_8449():
    p = Symbol('p', nonnegative=True)
    assert Lt(-oo, p)
    assert Ge(-oo, p) is S.false
    assert Gt(oo, -p)
    assert Le(oo, -p) is S.false


def test_simplify():
    assert simplify(x*(y + 1) - x*y - x + 1 < x) == (x > 1)
    assert simplify(S(1) < -x) == (x < -1)


def test_equals():
    w, x, y, z = symbols('w:z')
    f = Function('f')
    assert Eq(x, 1).equals(Eq(x*(y + 1) - x*y - x + 1, x))
    assert Eq(x, y).equals(x < y, True) == False
    assert Eq(x, f(1)).equals(Eq(x, f(2)), True) == f(1) - f(2)
    assert Eq(f(1), y).equals(Eq(f(2), y), True) == f(1) - f(2)
    assert Eq(x, f(1)).equals(Eq(f(2), x), True) == f(1) - f(2)
    assert Eq(f(1), x).equals(Eq(x, f(2)), True) == f(1) - f(2)
    assert Eq(w, x).equals(Eq(y, z), True) == False
    assert Eq(f(1), f(2)).equals(Eq(f(3), f(4)), True) == f(1) - f(3)
    assert (x < y).equals(y > x, True) == True
    assert (x < y).equals(y >= x, True) == False
    assert (x < y).equals(z < y, True) == False
    assert (x < y).equals(x < z, True) == False
    assert (x < f(1)).equals(x < f(2), True) == f(1) - f(2)
    assert (f(1) < x).equals(f(2) < x, True) == f(1) - f(2)


def test_reversed():
    assert (x < y).reversed == (y > x)
    assert (x <= y).reversed == (y >= x)
    assert Eq(x, y, evaluate=False).reversed == Eq(y, x, evaluate=False)
    assert Ne(x, y, evaluate=False).reversed == Ne(y, x, evaluate=False)
    assert (x >= y).reversed == (y <= x)
    assert (x > y).reversed == (y < x)


def test_canonical():
    one = S(1)

    def unchanged(v):
        c = v.canonical
        return v.is_Relational and c.is_Relational and v == c

    def isreversed(v):
        return v.canonical == v.reversed

    assert unchanged(x < one)
    assert unchanged(x <= one)
    assert isreversed(Eq(one, x, evaluate=False))
    assert unchanged(Eq(x, one, evaluate=False))
    assert isreversed(Ne(one, x, evaluate=False))
    assert unchanged(Ne(x, one, evaluate=False))
    assert unchanged(x >= one)
    assert unchanged(x > one)

    assert unchanged(x < y)
    assert unchanged(x <= y)
    assert isreversed(Eq(y, x, evaluate=False))
    assert unchanged(Eq(x, y, evaluate=False))
    assert isreversed(Ne(y, x, evaluate=False))
    assert unchanged(Ne(x, y, evaluate=False))
    assert isreversed(x >= y)
    assert isreversed(x > y)
    assert (-x < 1).canonical == (x > -1)
    assert isreversed(-x > y)


@XFAIL
def test_issue_8444():
    x = symbols('x', real=True)
    assert (x <= oo) == (x >= -oo) == True

    x = symbols('x')
    assert x >= floor(x)
    assert (x < floor(x)) == False
    assert Gt(x, floor(x)) == Gt(x, floor(x), evaluate=False)
    assert Ge(x, floor(x)) == Ge(x, floor(x), evaluate=False)
    assert x <= ceiling(x)
    assert (x > ceiling(x)) == False
    assert Lt(x, ceiling(x)) == Lt(x, ceiling(x), evaluate=False)
    assert Le(x, ceiling(x)) == Le(x, ceiling(x), evaluate=False)
    i = symbols('i', integer=True)
    assert (i > floor(i)) == False
    assert (i < ceiling(i)) == False
