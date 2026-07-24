"""Tests for ``EML`` and the ``.rewrite(EML)`` conversion."""
from __future__ import annotations

from sympy.core.function import ArgumentIndexError
from sympy.core.numbers import E, Rational, nan, oo, zoo
from sympy.core.symbol import Symbol, symbols
from sympy.functions.elementary.exponential import (
    EML, exp, log, to_eml, from_eml)
from sympy.functions.elementary.hyperbolic import (
    sinh, cosh, tanh, coth, sech, csch,
    asinh, acosh, atanh, acoth)
from sympy.functions.elementary.trigonometric import (
    sin, cos, tan, cot, sec, csc,
    asin, acos, atan)
from sympy.sets.sets import FiniteSet
from sympy.simplify.simplify import simplify
from sympy.testing.pytest import raises
import pickle
from sympy import S, srepr, sympify, I


x, y, z = symbols('x y z')


# ---------------------------------------------------------------------------
# Construction and argument validation
# ---------------------------------------------------------------------------

def test_eml_nargs():
    assert EML.nargs == FiniteSet(2)


def test_eml_arg_count():
    # SymPy derives nargs from the eval() signature and raises TypeError on
    # wrong arg counts.
    raises(TypeError, lambda: EML(x))
    raises(TypeError, lambda: EML(x, y, z))
    raises(TypeError, lambda: EML())


def test_eml_basic_construction():
    # Symbolic stays unevaluated.
    e = EML(x, y)
    assert isinstance(e, EML)
    assert e.args == (x, y)


# ---------------------------------------------------------------------------
# Eval simplifications
# ---------------------------------------------------------------------------

def test_eml_eval_simplifications():
    # EML does not auto-evaluate into exp/log form: that conversion is manual
    # (via rewrite / expand_func / from_eml), so these all stay as EML nodes.
    assert EML(x, 1) == EML(x, 1)
    assert isinstance(EML(x, 1), EML)
    assert isinstance(EML(0, 1), EML)
    assert isinstance(EML(1, 1), EML)
    assert isinstance(EML(2, 1), EML)
    assert isinstance(EML(-oo, y), EML)
    assert isinstance(EML(-oo, 2), EML)
    assert isinstance(EML(x, 0), EML)
    # Even two concrete Numbers stay unevaluated.
    assert isinstance(EML(2, 3), EML)
    assert isinstance(EML(0, 2), EML)

    # The manual conversion still gives the expected exp/log form.
    assert EML(x, 1).rewrite(exp) == exp(x)
    assert EML(-oo, y).rewrite(exp) == -log(y)
    assert EML(2, 3).rewrite(exp) == exp(2) - log(3)

    # NaN must still propagate from either argument...
    assert EML(nan, y) is nan
    assert EML(x, nan) is nan
    # ... and from genuine NaN results: exp(zoo) is NaN.
    assert EML(zoo, 1) is nan
    assert EML(zoo, x) is nan


def test_eml_no_overzealous_eval():
    # EML(x, y) with symbolic args stays as EML.
    assert isinstance(EML(x, y), EML)
    # E is a NumberSymbol (not is_Number); leave EML(2, E) symbolic.
    assert isinstance(EML(2, E), EML) or EML(2, E) == exp(2) - log(E)


# ---------------------------------------------------------------------------
# Differentiation
# ---------------------------------------------------------------------------

def test_eml_fdiff():
    # d/dx EML(x, y) = exp(x)
    assert EML(x, y).diff(x) == exp(x)
    # d/dy EML(x, y) = -1/y
    assert EML(x, y).diff(y) == -1/y
    # Mixed
    assert EML(x, y).diff(x, y) == 0
    # Higher order in x
    assert EML(x, y).diff(x, 2) == exp(x)
    # Differentiate when one arg is a function of the symbol
    assert EML(x**2, y).diff(x) == 2*x*exp(x**2)
    raises(ArgumentIndexError, lambda: EML(x, y).fdiff(3))
    raises(ArgumentIndexError, lambda: EML(x, y).fdiff(0))


def test_eml_diff_chain():
    # Chain rule across the second argument.
    assert EML(x, x**2).diff(x) == exp(x) - 2/x


# ---------------------------------------------------------------------------
# Rewrites of EML itself
# ---------------------------------------------------------------------------

def test_eml_rewrite_as_exp_log():
    assert EML(x, y).rewrite(exp) == exp(x) - log(y)
    assert EML(x, y).rewrite(log) == exp(x) - log(y)
    # Round-trip via simplify.
    assert simplify(EML(x, y).rewrite(exp) - (exp(x) - log(y))) == 0


# ---------------------------------------------------------------------------
# Numerical evaluation
# ---------------------------------------------------------------------------

def test_eml_evalf():
    # EML(1, 2) = e - log(2)
    val = EML(1, 2).evalf()
    expected = (exp(1) - log(2)).evalf()
    assert abs(val - expected) < 1e-12

    val2 = EML(0, 1).evalf()
    assert abs(val2 - 1) < 1e-12


# ---------------------------------------------------------------------------
# Assumption hooks
# ---------------------------------------------------------------------------

def test_eml_is_extended_real():
    a = Symbol('a', real=True)
    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    assert EML(a, p).is_extended_real is True
    # y not known positive: cannot conclude real
    assert EML(a, n).is_extended_real is not True


def test_eml_is_finite():
    a = Symbol('a', finite=True)
    p = Symbol('p', finite=True, positive=True)
    assert EML(a, p).is_finite is True


def test_eml_conjugate():
    from sympy.functions.elementary.complexes import conjugate
    assert conjugate(EML(x, y)) == EML(conjugate(x), conjugate(y))


# ---------------------------------------------------------------------------
# Conversion method: expr.rewrite(EML)
# ---------------------------------------------------------------------------

def _equivalent(a, b, sample=None):
    """Test mathematical equivalence of two expressions.

    Tries a symbolic round-trip via simplify(); if that fails (e.g. because
    SymPy's simplifier cannot fold ``-log(1/y) - log(y)`` without knowing
    ``y > 0``), falls back to a numerical check at a sample point.
    """
    delta = simplify((a - b).rewrite(exp))
    if delta == 0:
        return True
    # Try with expand_log(force=True) which assumes branches conveniently.
    from sympy import expand_log
    if expand_log(delta, force=True) == 0:
        return True
    # Numerical fallback at a positive real sample point.
    if sample is None:
        sample = {s: Rational(7, 10) + 1 for s in delta.free_symbols}
    val = delta.subs(sample).evalf()
    return abs(complex(val)) < 1e-10


def test_rewrite_exp():
    # exp(x).rewrite(EML) must contain EML, not collapse back.
    r = exp(x).rewrite(EML)
    assert r.has(EML)
    assert r == EML(x, 1, evaluate=False)
    # Round-trip.
    assert r.rewrite(exp) == exp(x)


def test_rewrite_log():
    r = log(x).rewrite(EML)
    assert r.has(EML)
    # Pure-EML grammar (Odrzywolek 2026):
    #     log(x) = EML(1, EML(EML(1, x), 1))
    # Constructed with evaluate=False to preserve the tree structure.
    expected_inner = EML(1, x, evaluate=False)
    expected_middle = EML(expected_inner, 1, evaluate=False)
    expected = EML(1, expected_middle, evaluate=False)
    assert r == expected
    # Round-trip via exp for a positive symbol is mathematically equal to log.
    xp = Symbol('x', positive=True)
    assert simplify(log(xp).rewrite(EML).rewrite(exp) - log(xp)) == 0


def test_rewrite_log_pure_eml_tree_shape():
    # The result of log(x).rewrite(EML) is a pure binary tree:
    # only EML nodes and the constants 1 / x as leaves. No Add, Mul, Pow.
    r = log(x).rewrite(EML)
    from sympy.core.add import Add
    from sympy.core.mul import Mul
    from sympy.core.power import Pow
    assert not r.has(Add)
    assert not r.has(Mul)
    assert not r.has(Pow)
    # Three nested EML nodes.
    assert isinstance(r, EML)
    assert isinstance(r.args[1], EML)
    assert isinstance(r.args[1].args[0], EML)


def test_rewrite_log_recursive_nesting_depth():
    # log(log(x)) should produce a doubly-nested pure-EML tree.
    e = log(log(x))
    r = e.rewrite(EML)
    assert r.has(EML)
    # Round-trip for a sufficiently large positive value where log(log(x)) > 0.
    xp = Symbol('x', positive=True)
    r_pos = log(log(xp)).rewrite(EML)
    diff_val = (r_pos - log(log(xp))).subs(xp, 10).rewrite(exp).evalf()
    assert abs(complex(diff_val)) < 1e-10


def test_rewrite_sums_and_products():
    e = exp(x) + log(y)
    r = e.rewrite(EML)
    assert r.has(EML)
    assert _equivalent(r, e)

    e2 = exp(x) * log(y)
    r2 = e2.rewrite(EML)
    assert r2.has(EML)
    assert _equivalent(r2, e2)

    e3 = exp(x) - 3*log(y) + 7
    assert _equivalent(e3.rewrite(EML), e3)


def test_rewrite_trig_direct():
    # Direct trig: numerical round-trip at a sample point.
    for fn in (sin, cos, tan, cot, sec, csc):
        e = fn(x)
        r = e.rewrite(EML)
        assert r.has(EML), f'{fn.__name__}.rewrite(EML) lost EML'
        diff_val = (r - e).subs(x, Rational(7, 10)).rewrite(exp).evalf()
        assert abs(complex(diff_val)) < 1e-10, f'{fn.__name__} numerical mismatch'


def test_rewrite_hyperbolic_direct():
    for fn in (sinh, cosh, tanh, coth, sech, csch):
        e = fn(x)
        r = e.rewrite(EML)
        assert r.has(EML), f'{fn.__name__}.rewrite(EML) lost EML'
        diff_val = (r - e).subs(x, Rational(7, 10)).rewrite(exp).evalf()
        assert abs(complex(diff_val)) < 1e-10, f'{fn.__name__} numerical mismatch'


def test_rewrite_inverse_trig():
    for fn in (asin, acos, atan):
        e = fn(x)
        r = e.rewrite(EML)
        assert r.has(EML), f'{fn.__name__}.rewrite(EML) lost EML'
        diff_val = (r - e).subs(x, Rational(1, 2)).rewrite(exp).evalf()
        assert abs(complex(diff_val)) < 1e-10, f'{fn.__name__} numerical mismatch'


def test_rewrite_inverse_hyperbolic():
    # acosh has a different domain; test at a point > 1 for it separately.
    for fn in (asinh, atanh, acoth):
        e = fn(x)
        r = e.rewrite(EML)
        assert r.has(EML), f'{fn.__name__}.rewrite(EML) lost EML'
        sample = Rational(1, 3) if fn is atanh else Rational(1, 2)
        if fn is acoth:
            sample = Rational(3, 2)  # |x| > 1
        diff_val = (r - e).subs(x, sample).rewrite(exp).evalf()
        assert abs(complex(diff_val)) < 1e-10, f'{fn.__name__} numerical mismatch'

    e = acosh(x)
    r = e.rewrite(EML)
    assert r.has(EML)
    diff_val = (r - e).subs(x, Rational(3, 2)).rewrite(exp).evalf()
    assert abs(complex(diff_val)) < 1e-10


def test_rewrite_compositions():
    # Pythagorean identity expressed via EML still equals 1.
    e = sin(x)**2 + cos(x)**2
    r = e.rewrite(EML)
    assert r.has(EML)
    # simplify of the rewritten form should reach 1.
    assert simplify(r.rewrite(exp)) == 1

    # Nested composition. exp(log(x)) auto-simplifies to x, so use a form
    # that does not collapse: log of a composite. SymPy's simplify cannot
    # always fold ``-log(1/u) - log(u)`` without ``u > 0``, so verify
    # numerically instead.
    e2 = log(exp(x) + 1)
    r2 = e2.rewrite(EML)
    assert r2.has(EML)
    sample_val = (r2 - e2).subs(x, Rational(1, 2)).rewrite(exp).evalf()
    assert abs(complex(sample_val)) < 1e-10


def test_rewrite_constants():
    # EML(0, 1) stays an EML node; its exp/log form is 1.
    assert isinstance(EML(0, 1), EML)
    assert EML(0, 1).rewrite(exp) == 1
    # exp(1) auto-evaluates to E (NumberSymbol singleton) at construction, so
    # `.rewrite(EML)` on the bare literal `exp(1)` cannot reach an EML form.
    # Within a non-collapsing expression, exp(x) rewrites to EML(x, 1).
    r = (exp(x) + 2).rewrite(EML)
    assert r.has(EML)
    # log of a compound argument is rewritten via the pure-form tree.
    rl = log(2 * y).rewrite(EML)
    assert rl.has(EML)


def test_paper_log_identity():
    """The paper's explicit pure-EML form for ln(x)."""
    r = log(x).rewrite(EML)
    # Step through the paper's formula: ln x = eml(1, eml(eml(1, x), 1)).
    t1 = EML(1, x, evaluate=False)            # exp(1) - log(x) = e - log(x)
    t2 = EML(t1, 1, evaluate=False)           # exp(t1) - 0    = exp(e - log(x))
    t3 = EML(1, t2, evaluate=False)           # e - log(t2)    = log(x)
    assert r == t3
    # Numeric verification at a positive sample.
    xp = Symbol('x', positive=True)
    assert simplify(log(xp).rewrite(EML).rewrite(exp) - log(xp)) == 0


def test_rewrite_passthrough_for_unsupported():
    # A symbol has no rewrite; passes through unchanged.
    assert x.rewrite(EML) == x

    # A function with no exp/log rewrite (LambertW) should pass through.
    from sympy.functions.elementary.exponential import LambertW
    expr = LambertW(x) + exp(x)
    r = expr.rewrite(EML)
    # exp(x) part rewrites; LambertW stays.
    assert r.has(EML)
    assert r.has(LambertW)


def test_rewrite_polynomial_in_symbol():
    # With the Pow -> EML rewrite, x**2 = exp(2*log(x)) is converted to EML
    # form while the linear term and the constant are left untouched.
    r = (x**2 + 3*x + 2).rewrite(EML)
    assert r.has(EML)
    assert _equivalent(r, x**2 + 3*x + 2)


# ---------------------------------------------------------------------------
# Convenience accessors
# ---------------------------------------------------------------------------

def test_eml_arg_properties():
    e = EML(x + 1, y * 2)
    assert e.exp_arg == x + 1
    assert e.log_arg == y * 2
    # args still primary tuple.
    assert e.args == (x + 1, 2*y)


# ---------------------------------------------------------------------------
# Extra rewrites
# ---------------------------------------------------------------------------

def test_eml_rewrite_as_sinh_cosh():
    r = EML(x, y).rewrite(sinh)
    # exp(x) = sinh(x) + cosh(x); so result should contain those.
    assert r.has(sinh) and r.has(cosh)
    assert simplify(r - (exp(x) - log(y))) == 0

    r2 = EML(x, y).rewrite(cosh)
    assert r2.has(sinh) and r2.has(cosh)
    assert simplify(r2 - (exp(x) - log(y))) == 0


def test_eml_rewrite_as_pow():
    # E**x and exp(x) are the same in SymPy's canonical form.
    r = EML(x, y).rewrite('Pow')
    assert simplify(r - (exp(x) - log(y))) == 0


def test_eml_rewrite_as_tractable():
    # tractable form is used internally by limit/gruntz machinery; it should
    # equal the exp/log form.
    r = EML(x, y).rewrite('tractable')
    assert r == exp(x) - log(y)


# ---------------------------------------------------------------------------
# Expansion
# ---------------------------------------------------------------------------

def test_eml_expand_func():
    assert EML(x, y).expand(func=True) == exp(x) - log(y)
    # expand also works through compound expressions.
    e = EML(x, y) + EML(0, x)
    assert e.expand(func=True) == exp(x) - log(y) + 1 - log(x)


# ---------------------------------------------------------------------------
# nth derivative
# ---------------------------------------------------------------------------

def test_eml_nth_derivative():
    # Derivatives wrt the first argument: always exp(x).
    assert EML(x, y).diff(x, 4) == exp(x)
    assert EML(x, y).diff(x, 10) == exp(x)
    # Derivatives wrt the second argument: (-1)^n (n-1)! / y^n
    assert EML(x, y).diff(y, 1) == -1/y
    assert EML(x, y).diff(y, 2) == 1/y**2
    assert EML(x, y).diff(y, 3) == -2/y**3
    assert EML(x, y).diff(y, 4) == 6/y**4


# ---------------------------------------------------------------------------
# as_real_imag
# ---------------------------------------------------------------------------

def test_eml_as_real_imag():
    from sympy.functions.elementary.complexes import Abs, arg as arg_f
    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    re_part, im_part = EML(a, b).as_real_imag()
    # For real a, b: re = exp(a) - log(|b|), im = -arg(b).
    expected = (exp(a) - log(Abs(b)), -arg_f(b))
    assert (re_part, im_part) == expected


# ---------------------------------------------------------------------------
# Series and leading term
# ---------------------------------------------------------------------------

def test_eml_leading_term():
    # As x -> 0: exp(x) -> 1, log(y) stays. Leading term = 1 - log(y).
    assert EML(x, y).as_leading_term(x) == 1 - log(y)


def test_eml_nseries():
    from sympy.series.order import Order
    s = EML(x, y).series(x, 0, 4)
    # Should expand exp(x) to Taylor series; log(y) is constant.
    assert s.removeO() == 1 - log(y) + x + x**2/2 + x**3/6
    assert s.getO() == Order(x**4, x)


# ---------------------------------------------------------------------------
# Refined assumption hooks
# ---------------------------------------------------------------------------

def test_eml_is_extended_positive():
    # 0 < y < 1 => -log(y) > 0 => EML(0, y) > 0
    assert EML(0, Rational(1, 2)).is_extended_positive is True


def test_eml_is_zero():
    # Trivial zero.
    assert EML(0, 1).is_zero is False  # EML(0, 1) is 1, not zero
    # Symbolic root: y = exp(exp(x)).
    assert EML(1, exp(E)).is_zero is True


def test_eml_is_rational_algebraic():
    # EML(0, 1) = 1 => rational & algebraic.
    assert EML(0, 1).is_rational is True
    assert EML(0, 1).is_algebraic is True
    # exp of nonzero algebraic is transcendental => EML(a, 1) is not algebraic.
    a = Symbol('a', algebraic=True, nonzero=True)
    assert EML(a, 1).is_algebraic is False
    # log of algebraic != 1 is transcendental => EML(0, b) is not algebraic
    # when b is algebraic and b != 1.
    b = Symbol('b', algebraic=True, positive=True)
    # Cannot conclude without ruling out b == 1; with rational b != 1 we can.
    b2 = Symbol('b2', rational=True, positive=True)
    # If b2 != 1 we'd know; without that information, the answer is None.
    # Just check that the framework doesn't crash and gives a fuzzy result.
    res = EML(0, b).is_algebraic
    assert res in (False, None)
    assert EML(0, b2).is_algebraic in (False, None)


# ---------------------------------------------------------------------------
# Extra edge-case tests added for deep EML trees, branch cuts, and serialisation
# ---------------------------------------------------------------------------


def test_eml_complex_branch_extra():
    x = Symbol('x')
    z = -1 + 0.1 * I
    r = log(x).rewrite(EML)
    v1 = r.rewrite('exp').subs(x, z).evalf()
    v2 = log(z).evalf()
    assert abs(complex(v1 - v2)) < 1e-10


def test_eml_log_zero_behavior_extra():
    x = Symbol('x')
    assert EML(x, 0).rewrite('exp') == exp(x) + zoo
    expr = EML(S.NegativeInfinity, x)
    val = expr.subs(x, 0.1).rewrite('exp').evalf()
    assert abs(complex(val - (-log(0.1).evalf()))) < 1e-10


def test_eml_pickle_srepr_roundtrip_extra():
    a = Symbol('a')
    e = EML(1, EML(1, a, evaluate=False), evaluate=False)
    p = pickle.loads(pickle.dumps(e))
    assert p == e
    s = srepr(e)
    assert sympify(s) == e


def test_eml_deep_tree_smoke_extra():
    def make_deep(depth: int):
        node = S.One
        for _ in range(depth):
            node = EML(node, S.One, evaluate=False)
        return node

    for d in (10, 50, 100):
        e = make_deep(d)
        assert e is not None
        if d > 0:
            assert e.has(EML)


# ---------------------------------------------------------------------------
# Pow -> EML conversion (a**b = exp(b*log(a)))
# ---------------------------------------------------------------------------

def test_pow_rewrite_as_eml():
    a = Symbol('a', positive=True)
    b = Symbol('b')
    r = (a**b).rewrite(EML)
    assert r.has(EML)
    # a**b = exp(b*log(a)) = EML(b*log(a)-in-EML, 1)
    assert r == EML(b*log(a).rewrite(EML), 1, evaluate=False)


def test_pow_rewrite_as_eml_roundtrip():
    a = Symbol('a', positive=True)
    b = Symbol('b')
    # Round-trips back to a**b (checked numerically; SymPy cannot fold the
    # EML encoding of log symbolically without forcing branch assumptions).
    assert _equivalent((a**b).rewrite(EML), a**b)


def test_sqrt_rewrite_as_eml():
    a = Symbol('a', positive=True)
    from sympy import sqrt
    r = sqrt(a).rewrite(EML)
    assert r.has(EML)
    assert _equivalent(r, sqrt(a))


def test_pow_numeric_base_rewrite_as_eml():
    # 2**x = exp(x*log(2)) -> EML form.
    r = (2**x).rewrite(EML)
    assert r.has(EML)
    assert _equivalent(r, 2**x)


# ---------------------------------------------------------------------------
# Public to_eml() transform
# ---------------------------------------------------------------------------

def test_to_eml_basic():
    assert to_eml(exp(x)) == EML(x, 1, evaluate=False)
    assert to_eml(log(x)) == log(x).rewrite(EML)
    assert to_eml(exp(x)) == exp(x).rewrite(EML)


def test_to_eml_keeps_arithmetic_structure():
    from sympy.core.add import Add
    r = to_eml(exp(x) + log(y))
    assert isinstance(r, Add)
    assert r == exp(x).rewrite(EML) + log(y).rewrite(EML)


def test_to_eml_power():
    a = Symbol('a', positive=True)
    b = Symbol('b')
    assert to_eml(a**b) == (a**b).rewrite(EML)
    assert to_eml(a**b).has(EML)


def test_to_eml_iterable():
    out = to_eml([exp(x), log(x)])
    assert out == [EML(x, 1, evaluate=False), log(x).rewrite(EML)]
    assert isinstance(out, list)
    out_t = to_eml((exp(x), log(x)))
    assert isinstance(out_t, tuple)


def test_to_eml_deep_flag():
    # With deep=False only the outermost function is rewritten; the inner
    # sin is left intact. deep=True rewrites everything.
    e = exp(sin(x))
    shallow = to_eml(e, deep=False)
    deep = to_eml(e, deep=True)
    assert shallow.has(EML) and shallow.has(sin)
    assert not deep.has(sin)
    assert deep.count(EML) > shallow.count(EML)


def test_to_eml_sympifies_input():
    assert to_eml(1) == S.One
    assert to_eml('exp(x)') == exp(x).rewrite(EML)


def test_to_eml_numbers_flag():
    # e = EML(1, 1) is the one numeric constant with a finite EML encoding.
    assert to_eml(E, numbers=True) == EML(S.One, S.One, evaluate=False)
    assert to_eml(E, numbers=True).expand(func=True) == E
    # Off by default: e is left as-is.
    assert to_eml(E) == E
    # General integers have no finite EML encoding and are left untouched.
    assert to_eml(2, numbers=True) == S(2)
    # e nested inside an expression is encoded too.
    r = to_eml(exp(x) + E, numbers=True)
    assert r.has(EML) and not r.has(E)


# ---------------------------------------------------------------------------
# from_eml() inverse transform
# ---------------------------------------------------------------------------

def test_from_eml_basic():
    assert from_eml(EML(x, y)) == exp(x) - log(y)
    assert from_eml(EML(x, 1, evaluate=False)) == exp(x)


def test_from_eml_no_eml_remains():
    expr = log(x) + exp(y)
    back = from_eml(to_eml(expr))
    assert not back.has(EML)


def test_from_eml_roundtrip_exp():
    # exp round-trips exactly; log round-trips up to simplify for x > 0.
    xp = Symbol('xp', positive=True)
    assert from_eml(to_eml(exp(x))) == exp(x)
    assert simplify(from_eml(to_eml(log(xp))) - log(xp)) == 0


def test_from_eml_iterable():
    out = from_eml([EML(x, 1, evaluate=False), EML(1, y, evaluate=False)])
    assert out == [exp(x), E - log(y)]
    assert isinstance(out, list)
    assert isinstance(from_eml((EML(x, y),)), tuple)


def test_from_eml_equals_rewrite_exp():
    expr = sin(x).rewrite(EML)
    assert from_eml(expr) == expr.rewrite(exp)


def test_to_from_eml_numbers_roundtrip():
    # to_eml(..., numbers=True) then from_eml recovers e.
    assert from_eml(to_eml(E, numbers=True)) == E


def test_from_eml_only_rewrites_eml():
    # from_eml must convert *only* EML nodes; non-EML functions (here sin)
    # are left untouched rather than expanded into exponential form.
    assert from_eml(sin(y) + EML(x, y)) == sin(y) + exp(x) - log(y)
    assert from_eml(cos(x) * EML(x, 1, evaluate=False)) == cos(x) * exp(x)


def test_eml_rewrite_then_subs():
    # With manual (non-eager) evaluation, substituting into the pure-EML form
    # of log keeps the tree intact instead of mangling it into exp(E)/2.
    r = log(x).rewrite(EML).subs(x, 2)
    assert r == EML(1, EML(EML(1, 2, evaluate=False), 1, evaluate=False),
                    evaluate=False)
    # ... and it expands back to log(2).
    assert from_eml(r) == log(2)
    assert from_eml(-EML(-oo, 2, evaluate=False)) == log(2)


def test_eml_nested_recursive_reduction():
    # A nested EML expression must keep its structure (no eager collapse into
    # exp(EML(...))) and round-trip cleanly back to log(x).
    nested = EML(1, EML(EML(1, x), 1))
    assert nested == log(x).rewrite(EML)
    assert from_eml(nested) == log(x)
    # Full round trip: to_eml(ln(x)) -> from_eml -> log(x).
    assert to_eml(log(x)) == EML(1, EML(EML(1, x, evaluate=False),
                                        1, evaluate=False), evaluate=False)
    assert from_eml(to_eml(log(x))) == log(x)


# ---------------------------------------------------------------------------
# Numeric evaluation and code generation (lambdify / code printers)
# ---------------------------------------------------------------------------

def test_eml_lambdify():
    from sympy import lambdify
    f = lambdify((x, y), EML(x, y))
    assert abs(f(0.0, 1.0) - 1.0) < 1e-12
    # e - log(2)
    expected = float((exp(1) - log(2)).evalf())
    assert abs(f(1.0, 2.0) - expected) < 1e-9


# NB: the code-printer assertions for EML now live alongside each printer's
# own tests -- see test_ccode_EML (printing/tests/test_c.py),
# test_fcode_EML (test_fortran.py), test_octave_EML (test_octave.py) and
# test_numpy_EML (test_numpy.py).
