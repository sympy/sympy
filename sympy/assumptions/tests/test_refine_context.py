from __future__ import annotations

from sympy import Q, Abs, Max, Min, sign, sqrt
from sympy.abc import x, y
from sympy.assumptions.refine import RefineContext, refine, refine_abs, refine_sign
from sympy.functions.elementary.integers import ceiling, floor
from sympy.testing.pytest import XFAIL


def test_refine_context_caches_ask_results():
    ctx = RefineContext(Q.positive(x))

    assert ctx.ask(Q.positive(x)) is True
    assert Q.positive(x) in ctx.cache
    assert ctx.ask(Q.positive(x)) is True


def test_refine_handlers_still_accept_plain_assumptions():
    assert refine_abs(Abs(x), Q.positive(x)) == x
    assert refine_sign(sign(x), Q.positive(x)) == 1


def test_refine_preserves_existing_recursive_behavior():
    assert refine(1 + Abs(x), Q.positive(x)) == 1 + x
    assert refine(sqrt(x**2), Q.real(x)) == Abs(x)
    assert refine(sqrt(x**2), Q.positive(x)) == x


def test_refine_preserves_existing_integer_floor_ceiling_behavior():
    assert refine(floor(x), Q.integer(x)) == x
    assert refine(ceiling(x), Q.integer(x)) == x
    assert refine(floor(x + y), Q.integer(x)) == x + floor(y)
    assert refine(ceiling(x + y), Q.integer(x)) == x + ceiling(y)


def test_refine_should_derive_positive_sum_from_positive_terms():
    assert refine(sqrt((x + y)**2), Q.positive(x) & Q.positive(y)) == x + y


@XFAIL
def test_refine_should_use_relative_assumptions_for_abs():
    assert refine(Abs(x), Q.gt(x, 1)) == x


@XFAIL
def test_refine_should_use_relative_assumptions_for_shifted_abs():
    assert refine(Abs(x - 2), Q.gt(x, 3)) == x - 2


@XFAIL
def test_refine_should_use_relative_assumptions_for_sign():
    assert refine(sign(x), Q.gt(x, 0)) == 1


@XFAIL
def test_refine_should_use_ordering_assumptions_for_max_min():
    assert refine(Max(x, y), Q.gt(x, y)) == x
    assert refine(Min(x, y), Q.gt(x, y)) == y
