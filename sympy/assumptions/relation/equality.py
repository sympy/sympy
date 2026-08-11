"""
Module for mathematical equality [1] and inequalities [2].

The purpose of this module is to provide the instances which represent the
binary predicates in order to combine the relationals into logical inference
system. Objects such as ``Q.eq``, ``Q.lt`` should remain internal to
assumptions module, and user must use the classes such as :obj:`~.Eq()`,
:obj:`~.Lt()` instead to construct the relational expressions.

References
==========

.. [1] https://en.wikipedia.org/wiki/Equality_(mathematics)
.. [2] https://en.wikipedia.org/wiki/Inequality_(mathematics)
"""
from __future__ import annotations
from sympy.assumptions import Q

from .binrel import BinaryRelation

__all__ = ['EqualityPredicate', 'UnequalityPredicate', 'StrictGreaterThanPredicate',
    'GreaterThanPredicate', 'StrictLessThanPredicate', 'LessThanPredicate']


def _query(predicate_ctor, expr, assumptions):
    """Resolve a *unary* predicate (Q.zero, Q.infinite, Q.extended_real, ...)
    on `expr`.

    First tries the cheap, already-known fuzzy property on the expression
    itself (e.g. `expr.is_zero`) - this alone answers the large majority of
    cases (structural facts, and facts implied by assumptions declared on the
    symbols that make up `expr`, e.g. `symbols('m', integer=True)`).

    If that's inconclusive, falls back to ``_ask_recursive`` - the fast
    fact/dispatch-only lookup introduced by GH-28129/#29921 - to also take
    local *assumptions* (the second argument to ask()) into account. Unlike
    ask(), ``_ask_recursive`` never invokes the expensive satask solver, so
    this cannot trigger another costly top-level ask() cycle - which is the
    exact complaint in GH-29863 (citing GH-28129) about calling is_eq/is_ge
    from here.
    """
    prop = getattr(expr, 'is_' + predicate_ctor.name, None)
    if prop is not None:
        return prop
    from sympy.assumptions.ask import _ask_recursive
    return _ask_recursive(predicate_ctor(expr), assumptions)


def _fuzzy_eq(lhs, rhs, assumptions=True):
    """Fuzzy equality check that avoids calling is_eq(), and therefore avoids
    the ask() <-> is_eq() mutual recursion described in GH-29863 (is_eq routes
    through AssumptionsWrapper straight back into ask(Q.eq(...))).
    """
    from sympy.core.numbers import NaN
    from sympy.core.logic import fuzzy_xor
    from sympy.logic.boolalg import BooleanAtom, Boolean

    # NaN is never equal to anything, not even itself - must be checked
    # before the structural `lhs == rhs` shortcut below, since sympy defines
    # NaN == NaN as structurally True.
    if isinstance(lhs, NaN) or isinstance(rhs, NaN):
        return False

    if lhs == rhs:
        return True

    if all(isinstance(i, BooleanAtom) for i in (lhs, rhs)):
        return False  # True != False (and they weren't equal above)

    if not (lhs.is_Symbol or rhs.is_Symbol) and (
            isinstance(lhs, Boolean) != isinstance(rhs, Boolean)):
        return False  # only Booleans can equal Booleans; avoids lhs - rhs below

    lhs_inf = _query(Q.infinite, lhs, assumptions)
    rhs_inf = _query(Q.infinite, rhs, assumptions)
    if lhs_inf or rhs_inf:
        if fuzzy_xor([lhs_inf, rhs_inf]):
            return False
        lhs_real = _query(Q.extended_real, lhs, assumptions)
        rhs_real = _query(Q.extended_real, rhs, assumptions)
        if fuzzy_xor([lhs_real, rhs_real]):
            return False
        # both infinite & both real (or both known non-real): fall through
        # to the diff-based check below, which handles +-oo vs +-oo.

    diff = lhs - rhs
    if isinstance(diff, NaN):
        # e.g. zoo - zoo, or an indeterminate difference between two
        # complex infinities that aren't structurally identical
        return None

    return _query(Q.zero, diff, assumptions)


def _fuzzy_ge(lhs, rhs, assumptions=True):
    """Fuzzy >= check that avoids calling is_ge(), and therefore avoids the
    ask() <-> is_ge() mutual recursion described in GH-29863.
    """
    from sympy.core.numbers import NaN

    if isinstance(lhs, NaN) or isinstance(rhs, NaN):
        return None

    if lhs == rhs:
        return True

    lhs_real = _query(Q.extended_real, lhs, assumptions)
    rhs_real = _query(Q.extended_real, rhs, assumptions)
    if not (lhs_real and rhs_real):
        # is_ge also bails out unless both sides are confirmed extended real
        return None

    if (_query(Q.infinite, lhs, assumptions) and
            _query(Q.extended_positive, lhs, assumptions)) or \
       (_query(Q.infinite, rhs, assumptions) and
            _query(Q.extended_negative, rhs, assumptions)):
        return True

    diff = lhs - rhs
    if isinstance(diff, NaN):
        return None

    return _query(Q.extended_nonnegative, diff, assumptions)


class EqualityPredicate(BinaryRelation):
    """
    Binary predicate for $=$.

    The purpose of this class is to provide the instance which represent
    the equality predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Eq()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_eq`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.eq(0, 0)
    Q.eq(0, 0)
    >>> ask(_)
    True

    See Also
    ========

    sympy.core.relational.Eq

    """
    is_reflexive = True
    is_symmetric = True

    name = 'eq'
    handler = None  # Do not allow dispatching by this predicate

    @property
    def negated(self):
        return Q.ne

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_eq() cycle (GH-29863)
        return _fuzzy_eq(*args, assumptions)


class UnequalityPredicate(BinaryRelation):
    r"""
    Binary predicate for $\neq$.

    The purpose of this class is to provide the instance which represent
    the inequation predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Ne()` instead to construct the inequation expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_neq`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.ne(0, 0)
    Q.ne(0, 0)
    >>> ask(_)
    False

    See Also
    ========

    sympy.core.relational.Ne

    """
    is_reflexive = False
    is_symmetric = True

    name = 'ne'
    handler = None

    @property
    def negated(self):
        return Q.eq

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_neq() cycle (GH-29863)
        res = _fuzzy_eq(*args, assumptions)
        return None if res is None else not res


class StrictGreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for $>$.

    The purpose of this class is to provide the instance which represent
    the ">" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Gt()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_gt`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.gt(0, 0)
    Q.gt(0, 0)
    >>> ask(_)
    False

    See Also
    ========

    sympy.core.relational.Gt

    """
    is_reflexive = False
    is_symmetric = False

    name = 'gt'
    handler = None

    @property
    def reversed(self):
        return Q.lt

    @property
    def negated(self):
        return Q.le

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_gt() cycle (GH-29863)
        # gt(a,b) := not ge(b,a)
        res = _fuzzy_ge(args[1], args[0], assumptions)
        return None if res is None else not res


class GreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for $>=$.

    The purpose of this class is to provide the instance which represent
    the ">=" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Ge()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_ge`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.ge(0, 0)
    Q.ge(0, 0)
    >>> ask(_)
    True

    See Also
    ========

    sympy.core.relational.Ge

    """
    is_reflexive = True
    is_symmetric = False

    name = 'ge'
    handler = None

    @property
    def reversed(self):
        return Q.le

    @property
    def negated(self):
        return Q.lt

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_ge() cycle (GH-29863)
        return _fuzzy_ge(*args, assumptions)


class StrictLessThanPredicate(BinaryRelation):
    """
    Binary predicate for $<$.

    The purpose of this class is to provide the instance which represent
    the "<" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Lt()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_lt`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.lt(0, 0)
    Q.lt(0, 0)
    >>> ask(_)
    False

    See Also
    ========

    sympy.core.relational.Lt

    """
    is_reflexive = False
    is_symmetric = False

    name = 'lt'
    handler = None

    @property
    def reversed(self):
        return Q.gt

    @property
    def negated(self):
        return Q.ge

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_lt() cycle (GH-29863)
        # lt(a,b) := not ge(a,b)
        res = _fuzzy_ge(*args, assumptions)
        return None if res is None else not res


class LessThanPredicate(BinaryRelation):
    """
    Binary predicate for $<=$.

    The purpose of this class is to provide the instance which represent
    the "<=" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Le()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_le`

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.le(0, 0)
    Q.le(0, 0)
    >>> ask(_)
    True

    See Also
    ========

    sympy.core.relational.Le

    """
    is_reflexive = True
    is_symmetric = False

    name = 'le'
    handler = None

    @property
    def reversed(self):
        return Q.ge

    @property
    def negated(self):
        return Q.gt

    def eval(self, args, assumptions=True):
        # Use non-recursive helper; avoids ask()<->is_le() cycle (GH-29863)
        # le(a,b) := ge(b,a)
        return _fuzzy_ge(args[1], args[0], assumptions)
