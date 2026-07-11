"""
Helpers for cross-checking the various backends that can answer an
assumptions query against each other, including the old assumptions
system.

``ask()`` is really a chain of increasingly expensive backends -- the
fact/recursive-handler based :func:`~sympy.assumptions.ask._ask_recursive`,
then :func:`~sympy.assumptions.satask.satask`, then
:func:`~sympy.assumptions.lra_satask.lra_satask` -- each of which is tried
in turn until one of them returns something other than ``None``. Because
of this, a bug in one backend can be hidden by another backend picking up
the slack for a particular query, which lets bugs slip through the
old, backend-agnostic ``assert ask(...) is ...`` style tests.

This module exposes each backend individually (see :data:`BACKENDS`),
together with :func:`old_ask`, a converter that answers (a restricted
subset of) queries using the old assumptions system, and :func:`check_ask`,
a test helper that runs a query through every backend and checks that:

- backends that are expected to know the answer give exactly the
  expected (``True``/``False``) result,
- every other backend either agrees or admits it doesn't know
  (returns ``None``), and
- the top-level :func:`~.ask` itself is consistent with the above.

See https://github.com/sympy/sympy/issues/28191 for the motivating
discussion.
"""
from sympy.assumptions.ask import Q, ask, _ask_recursive
from sympy.assumptions.assume import AppliedPredicate
from sympy.assumptions.cnf import CNF
from sympy.assumptions.satask import satask
from sympy.assumptions.lra_satask import lra_satask
from sympy.core.symbol import Symbol
from sympy.logic.algorithms.lra_theory import UnhandledInput


class BackendUnsupported(Exception):
    """
    Raised by a backend helper (see :data:`BACKENDS`) to signal that the
    given proposition/assumptions can't be evaluated by that backend at
    all, as opposed to the backend evaluating them and simply not being
    able to determine a result (which is reported as ``None``, not by
    raising).

    :func:`check_ask` silently skips backends that raise this.
    """


class OldAssumptionsUnsupported(BackendUnsupported):
    """Raised by :func:`old_ask` when a query cannot be represented using
    the old assumptions system."""


# Predicates that have a directly corresponding keyword argument in the
# old assumptions system (see ``sympy.core.assumptions._assume_defined``).
# Predicates that are absent here (matrix predicates, relations like
# ``Q.eq``, ``Q.hermitian``, ``Q.is_true``, ...) have no old-assumptions
# analogue, so ``old_ask`` refuses to evaluate them.
_OLD_ASSUMPTIONS_NAMES = {getattr(Q, name): name for name in (
    'commutative', 'complex', 'algebraic', 'transcendental',
    'extended_real', 'real', 'imaginary',
    'integer', 'noninteger', 'rational', 'irrational',
    'finite', 'infinite',
    'positive', 'negative', 'zero', 'nonzero', 'nonpositive', 'nonnegative',
    'extended_positive', 'extended_negative', 'extended_nonzero',
    'extended_nonpositive', 'extended_nonnegative',
    'even', 'odd', 'prime', 'composite',
)}


def _old_kwargs_from_assumptions(assumptions):
    """
    Convert *assumptions* -- a conjunction of unary predicates applied to
    plain ``Symbol``s -- into a ``{Symbol: {kwarg: bool}}`` mapping
    suitable for building old-assumptions ``Symbol``s.

    Raises :class:`OldAssumptionsUnsupported` if *assumptions* can't be
    represented this way, e.g. because it contains a disjunction (the old
    assumptions system can't express "x is real or x is a matrix"), a
    predicate with no old-assumptions analogue, or a predicate applied to
    something other than a bare ``Symbol``.
    """
    cnf = CNF.from_prop(assumptions)
    kwargs = {}
    for clause in cnf.clauses:
        if len(clause) != 1:
            raise OldAssumptionsUnsupported(f"disjunction in {assumptions}")
        (literal,) = clause
        lit = literal.lit
        if lit is True and not literal.is_Not:
            # trivially true literal, e.g. from CNF.from_prop(True)
            continue
        if not isinstance(lit, AppliedPredicate) or len(lit.arguments) != 1:
            raise OldAssumptionsUnsupported(f"non-unary literal {lit}")
        sym = lit.arguments[0]
        if not isinstance(sym, Symbol):
            raise OldAssumptionsUnsupported(f"assumption on non-symbol {sym}")
        name = _OLD_ASSUMPTIONS_NAMES.get(lit.function)
        if name is None:
            raise OldAssumptionsUnsupported(
                f"no old-assumptions analogue for {lit.function}")
        kwargs.setdefault(sym, {})[name] = not literal.is_Not
    return kwargs


def old_ask(proposition, assumptions=True):
    """
    Evaluate *proposition* given *assumptions* using the old assumptions
    system, as a (much more limited) analogue of :func:`~.ask`.

    This works by building fresh ``Symbol``s that carry the old-assumptions
    keyword arguments implied by *assumptions*, substituting them into
    *proposition*'s argument, and reading off the corresponding
    ``is_<name>`` property.

    Raises :class:`OldAssumptionsUnsupported` if *proposition* or
    *assumptions* can't be expressed in the old assumptions system, e.g.
    the proposition isn't a unary predicate with a corresponding
    ``is_<name>`` attribute, or a free symbol of the proposition already
    carries old-style assumptions of its own that aren't accounted for by
    *assumptions*.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.abc import x
    >>> from sympy.assumptions.tests.common import old_ask
    >>> old_ask(Q.positive(x), Q.positive(x))
    True
    >>> old_ask(Q.even(x**2), Q.even(x))
    True
    """
    if not isinstance(proposition, AppliedPredicate) or len(proposition.arguments) != 1:
        raise OldAssumptionsUnsupported(f"non-unary proposition {proposition}")

    name = _OLD_ASSUMPTIONS_NAMES.get(proposition.function)
    if name is None:
        raise OldAssumptionsUnsupported(
            f"no old-assumptions analogue for {proposition.function}")

    kwargs = _old_kwargs_from_assumptions(assumptions)

    expr = proposition.arguments[0]
    for sym in expr.free_symbols:
        if not isinstance(sym, Symbol):
            raise OldAssumptionsUnsupported(f"non-symbol atom {sym}")
        if sym not in kwargs and sym.assumptions0.keys() - {'commutative'}:
            # sym carries old-style assumptions of its own that aren't
            # accounted for by *assumptions*; too risky to guess intent.
            raise OldAssumptionsUnsupported(
                f"{sym} already carries old assumptions {sym.assumptions0}")

    replacements = {sym: Symbol(sym.name, **sym_kwargs)
                    for sym, sym_kwargs in kwargs.items()}
    new_expr = expr.xreplace(replacements)
    return getattr(new_expr, 'is_' + name)


def _lra_satask_backend(proposition, assumptions):
    try:
        return lra_satask(proposition, assumptions=assumptions)
    except UnhandledInput as e:
        raise BackendUnsupported(str(e)) from e


# The individual backends that ``ask()`` chains together, plus ``old_ask``
# for the old assumptions system. Each maps a name to a callable with
# signature ``(proposition, assumptions) -> True | False | None`` that may
# raise ``BackendUnsupported`` if it can't process the given input at all.
BACKENDS = {
    'recursive': _ask_recursive,
    'satask': satask,
    'lra_satask': _lra_satask_backend,
    'old': old_ask,
}

# The subset of BACKENDS that ask() itself is built out of, in the order
# ask() tries them. 'old' is a separate system entirely and is not part of
# the ask() chain.
_ASK_BACKENDS = ('recursive', 'satask', 'lra_satask')


def check_ask(proposition, assumptions=True, ideal_result=None,
              backends_with_ideal_result=()):
    """
    Check *proposition* given *assumptions* against every backend in
    :data:`BACKENDS`, as proposed in
    https://github.com/sympy/sympy/issues/28191.

    Parameters
    ==========

    proposition, assumptions : Boolean
        As for :func:`~.ask`.

    ideal_result : True, False or None
        The mathematically correct answer for *proposition* given
        *assumptions*.

    backends_with_ideal_result : iterable of str
        Names of the backends (a subset of ``BACKENDS``) that are
        expected to be strong enough to return *ideal_result* for this
        query. Every other backend that supports this query is expected
        to return ``None``, i.e. to correctly admit that it doesn't know
        rather than getting the wrong answer. Backends that raise
        :class:`BackendUnsupported` for the given input (e.g. ``old_ask``
        on a query it can't represent) are silently skipped.

        Since :func:`~.ask` itself is built out of the ``'recursive'``,
        ``'satask'`` and ``'lra_satask'`` backends (in that order), this
        also asserts that ``ask(proposition, assumptions)`` returns
        *ideal_result* whenever any of those three backends does, and
        ``None`` otherwise.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.abc import x
    >>> from sympy.assumptions.tests.common import check_ask

    A query that every backend capable of handling it should get right:

    >>> check_ask(Q.zero(x**2), Q.zero(x), True,
    ...     ['recursive', 'satask', 'old'])

    A query where only the recursive handlers and the old assumptions
    system currently know the answer; ``satask`` is (for now) expected
    to correctly return ``None`` instead of a wrong answer, rather than
    silently being bailed out by the recursive backend that ``ask()``
    tries first:

    >>> check_ask(Q.even(x**2), Q.even(x), True, ['recursive', 'old'])
    """
    backends_with_ideal_result = set(backends_with_ideal_result)
    unknown = backends_with_ideal_result - set(BACKENDS)
    if unknown:
        raise ValueError(f"unknown backend(s): {unknown}")

    for name, backend in BACKENDS.items():
        try:
            result = backend(proposition, assumptions)
        except BackendUnsupported:
            continue
        expected = ideal_result if name in backends_with_ideal_result else None
        assert result is expected, (
            f"backend {name!r} gave {result!r} for "
            f"ask({proposition}, {assumptions}), expected {expected!r}")

    expected_ask = (ideal_result
                     if backends_with_ideal_result & set(_ASK_BACKENDS)
                     else None)
    result = ask(proposition, assumptions)
    assert result is expected_ask, (
        f"ask({proposition}, {assumptions}) gave {result!r}, "
        f"expected {expected_ask!r}")
