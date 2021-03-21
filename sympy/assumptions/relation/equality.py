"""
Module for mathematical equality [1] and inequation [2].

References
==========

.. [1] https://en.wikipedia.org/wiki/Equality_(mathematics)
.. [2] https://en.wikipedia.org/wiki/Inequation
"""
from sympy.assumptions import Q
from sympy.core.relational import is_eq, is_neq, _eval_is_eq

from .binrel import BinaryRelation

__all__ = ['EqualityPredicate', 'UnequalityPredicate']


class EqualityPredicate(BinaryRelation):
    """
    Binary predicate for $=$.

    This uses :func:`sympy.core.relational.is_eq()` for evaluation. Dispatching
    the handler via ``Q.eq`` is supported.

    Examples
    ========

    >>> from sympy import ask, Q
    >>> Q.eq(0, 0)
    Q.eq(0, 0)
    >>> ask(_)
    True

    New types can be supported by registration.

    >>> from sympy import Basic
    >>> class MyBasic(Basic):
    ...     def __new__(cls, arg):
    ...         return super().__new__(cls, arg)
    >>> @Q.eq.register(MyBasic, MyBasic)
    ... def _(lhs, rhs):
    ...     return ask(Q.eq(lhs.args[0], rhs.args[0]))

    >>> ask(Q.eq(MyBasic(1), MyBasic(1)))
    True
    >>> ask(Q.eq(MyBasic(2), MyBasic(1)))
    False

    By dispatching to ``Q.eq``, ``MyBasic`` is supported by ``Q.ne`` as well.

    >>> ask(Q.ne(MyBasic(1), MyBasic(1)))
    False
    >>> ask(Q.ne(MyBasic(2), MyBasic(1)))
    True

    """

    is_reflexive = True
    is_symmetric = True

    name = 'eq'
    handler = _eval_is_eq   # this allows registering via Q.eq

    @property
    def negated(self):
        return Q.ne

    def eval(self, args, assumptions=True):
        return is_eq(*args)


class UnequalityPredicate(BinaryRelation):
    r"""
    Binary predicate for $\neq$.

    This predicate delegates evaluation logic to ``Q.eq`` and does not support
    multipledispatch handler. To support new types, dispatch them to ``Q.eq``.
    See the docstring of :obj:`sympy.assumptions.relation.equality.EqualityPredicate`.

    """
    # TODO: Add examples

    is_reflexive = False
    is_symmetric = True

    name = 'ne'
    handler = None

    @property
    def negated(self):
        return Q.eq

    def eval(self, args, assumptions=True):
        return is_neq(*args)
