"""
Module for mathematical equality [1] and inequation [2].

References
==========

.. [1] https://en.wikipedia.org/wiki/Equality_(mathematics)
.. [2] https://en.wikipedia.org/wiki/Inequation
"""
from sympy import S
from sympy.assumptions import Q
from sympy.core.relational import is_eq, is_neq, _eval_is_eq

from .binrel import BinaryRelation

__all__ = ['EqualityPredicate', 'UnequalityPredicate']


class EqualityPredicate(BinaryRelation):
    """
    Binary predicate for ``==``

    """
    # TODO: Add examples

    is_reflexive = True
    is_symmetric = True

    name = 'eq'
    str_name = latex_name = "="
    handler = _eval_is_eq   # this allows registering via Q.eq

    @property
    def negated(self):
        return Q.ne

    def eval(self, args, assumptions=True):
        return is_eq(*args)

    def _simplify_applied(self, lhs, rhs, **kwargs):
        return eqsimp(self(lhs, rhs), **kwargs)

    def _eval_binary_symbols(self, lhs, rhs):
        args = (lhs, rhs)
        if S.true in args or S.false in args:
            if lhs.is_Symbol:
                return {lhs}
            elif rhs.is_Symbol:
                return {rhs}
        return set()


class UnequalityPredicate(BinaryRelation):
    """
    Binary predicate for ``!=``.

    """
    # TODO: Add examples

    is_reflexive = False
    is_symmetric = True

    name = 'ne'
    str_name = "!="
    latex_name = r"\neq"
    handler = None

    @property
    def negated(self):
        return Q.eq

    def eval(self, args, assumptions=True):
        return is_neq(*args)

    def _simplify_applied(self, lhs, rhs, **kwargs):
        eq = Q.eq(lhs, rhs).simplify(**kwargs)
        return self(*eq.arguments)

    def _eval_binary_symbols(self, lhs, rhs):
        args = (lhs, rhs)
        if S.true in args or S.false in args:
            if lhs.is_Symbol:
                return {lhs}
            elif rhs.is_Symbol:
                return {rhs}
        return set()

from .simplify import eqsimp
