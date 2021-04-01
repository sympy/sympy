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
from sympy.assumptions import ask, Q
from sympy.core import Add, S, Expr
from sympy.core.logic import fuzzy_and, fuzzy_bool, fuzzy_not, fuzzy_xor
from sympy.core.relational import (
    is_eq, is_neq, is_gt, is_ge, is_lt, is_le, _n2
)
from sympy.core.sympify import _sympify
from sympy.functions.elementary.complexes import arg
from sympy.simplify.simplify import clear_coefficients
from sympy.utilities.iterables import sift

from .binrel import BinaryRelation

__all__ = ['EqualityPredicate', 'UnequalityPredicate', 'StrictGreaterThanPredicate',
    'GreaterThanPredicate', 'StrictLessThanPredicate', 'LessThanPredicate']


class EqualityPredicate(BinaryRelation):
    """
    Binary predicate for $=$.

    The purpose of this class is to provide the instance which represent
    the equality predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Eq()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_eq()`

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
        if assumptions == True:
            # if no assumption is made, use is_eq for fast computation
            return is_eq(*args)

        # compute numerical equality with assumptions
        # note that this is parallel with is_eq
        lhs, rhs = args

        lhs_inf = ask(Q.infinite(lhs), assumptions)
        rhs_inf = ask(Q.infinite(rhs), assumptions)
        if lhs_inf or rhs_inf:
            if fuzzy_xor([rhs_inf, rhs_inf]):
                return False
            lhs_extreal = ask(Q.extended_real(lhs), assumptions)
            rhs_extreal = ask(Q.extended_real(rhs), assumptions)
            if fuzzy_xor([lhs_extreal, rhs_extreal]):
                return False
            if fuzzy_and([rhs_extreal, rhs_extreal]):
                lhs_extpos = ask(Q.extended_positive(lhs), assumptions)
                rhs_extpos = ask(Q.extended_positive(rhs), assumptions)
                return fuzzy_xor([lhs_extpos, fuzzy_not(rhs_extpos)])

            # Try to split real/imaginary parts and equate them
            I = S.ImaginaryUnit

            def split_real_imag(expr):
                real_imag = lambda t: (
                    'real' if ask(Q.extended_real(t), assumptions) else
                    'imag' if ask(Q.extended_real(I*t), assumptions) else None)
                return sift(Add.make_args(expr), real_imag)

            lhs_ri = split_real_imag(lhs)
            if not lhs_ri[None]:
                rhs_ri = split_real_imag(rhs)
                if not rhs_ri[None]:
                    eq_real = ask(Q.eq(Add(*lhs_ri['real']), Add(*rhs_ri['real'])), assumptions)
                    eq_imag = ask(Q.eq(I * Add(*lhs_ri['imag']), I * Add(*rhs_ri['imag'])), assumptions)
                    return fuzzy_and(map(fuzzy_bool, [eq_real, eq_imag]))

            # Compare e.g. zoo with 1+I*oo by comparing args
            arglhs = arg(lhs)
            argrhs = arg(rhs)
            # Guard against Eq(nan, nan) -> Falsesymp
            if not (arglhs == S.NaN and argrhs == S.NaN):
                return fuzzy_bool(ask(Q.eq(arglhs, argrhs), assumptions))

        if all(isinstance(i, Expr) for i in (lhs, rhs)):
            # see if the difference evaluates
            dif = lhs - rhs
            z = ask(Q.zero(dif), assumptions)
            if z is not None:
                if z is False and ask(Q.commutative(dif), assumptions):  # issue 10728
                    return False
                if z:
                    return True

            n2 = _n2(lhs, rhs)
            if n2 is not None:
                return _sympify(n2 == 0)

            # see if the ratio evaluates
            n, d = dif.as_numer_denom()
            rv = None

            n_zero = ask(Q.zero(n), assumptions)
            n_finite = ask(Q.finite(n), assumptions)
            if n_zero:
                rv = ask(Q.nonzero(d), assumptions)
            elif n_finite:
                d_infinite = ask(Q.infinite(d), assumptions)
                if d_infinite:
                    rv = True
                elif n_zero is False:
                    rv = d_infinite
                    if rv is None:
                        # if the condition that makes the denominator
                        # infinite does not make the original expression
                        # True then False can be returned
                        l, r = clear_coefficients(d, S.Infinity)
                        args = [_.subs(l, r) for _ in (lhs, rhs)]
                        if args != [lhs, rhs]:
                            rv = fuzzy_bool(ask(Q.eq(*args), assumptions))
                            if rv is True:
                                rv = None
            elif any(ask(Q.infinite(a), assumptions) for a in Add.make_args(n)):
                # (inf or nan)/x != 0
                rv = False
            if rv is not None:
                return rv


class UnequalityPredicate(BinaryRelation):
    r"""
    Binary predicate for $\neq$.

    The purpose of this class is to provide the instance which represent
    the inequation predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Ne()` instead to construct the inequation expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_neq()`

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
        if assumptions == True:
            return is_neq(*args)
        return fuzzy_not(ask(Q.eq(*args), assumptions))


class StrictGreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for $>$.

    The purpose of this class is to provide the instance which represent
    the ">" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Gt()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_gt()`

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
        return is_gt(*args)


class GreaterThanPredicate(BinaryRelation):
    """
    Binary predicate for $>=$.

    The purpose of this class is to provide the instance which represent
    the ">=" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Ge()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_ge()`

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
        return is_ge(*args)


class StrictLessThanPredicate(BinaryRelation):
    """
    Binary predicate for $<$.

    The purpose of this class is to provide the instance which represent
    the "<" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Lt()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_lt()`

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
        return is_lt(*args)


class LessThanPredicate(BinaryRelation):
    """
    Binary predicate for $<=$.

    The purpose of this class is to provide the instance which represent
    the "<=" predicate in order to allow the logical inference.
    This class must remain internal to assumptions module and user must
    use :obj:`~.Le()` instead to construct the equality expression.

    Evaluating this predicate to ``True`` or ``False`` is done by
    :func:`~.core.relational.is_le()`

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
        return is_le(*args)
