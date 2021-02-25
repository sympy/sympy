"""
Module for mathematical equality [1] and inequation [2].

References
==========

.. [1] https://en.wikipedia.org/wiki/Equality_(mathematics)
.. [2] https://en.wikipedia.org/wiki/Inequation
"""
from sympy import S
from sympy.assumptions import Q
from sympy.core.relational import Eq, Ne, is_eq, is_neq, _eval_is_eq

from .binrel import BinaryRelation, _DeprecatedRelational

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
    rel_op = "=="
    latex_name = "="

    handler = _eval_is_eq   # this allows registering via Q.eq

    @property
    def negated(self):
        return Q.ne

    @property
    def as_old_relational(self):
        return Eq

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
    rel_op = "!="
    latex_name = r"\neq"
    fortran_name = "/="

    handler = None

    @property
    def negated(self):
        return Q.eq

    @property
    def as_old_relational(self):
        return Ne

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


# compatibility for old relational


class _DeprecatedEq(_DeprecatedRelational):

    is_Equality = True

    def __new__(cls, *args, **options):
        if len(args) == 2:
            # allow rel.__class__(rel.lhs, rel.rhs)
            args = (Q.eq, args[0], args[1])
        return super().__new__(cls, *args, **options)

    def _eval_rewrite_as_Add(self, *args, **kwargs):
        """
        return Eq(L, R) as L - R. To control the evaluation of
        the result set pass `evaluate=True` to give L - R;
        if `evaluate=None` then terms in L and R will not cancel
        but they will be listed in canonical order; otherwise
        non-canonical args will be returned.

        Examples
        ========

        >>> from sympy import Eq, Add
        >>> from sympy.abc import b, x
        >>> eq = Eq(x + b, x - b)
        >>> eq.rewrite(Add)
        2*b
        >>> eq.rewrite(Add, evaluate=None).args
        (b, b, x, -x)
        >>> eq.rewrite(Add, evaluate=False).args
        (b, x, b, -x)
        """
        if self.function is not Q.eq:
            return self

        from sympy.core.add import _unevaluated_Add, Add
        L, R = self.arguments
        evaluate = kwargs.get('evaluate', True)
        if evaluate:
            # allow cancellation of args
            return L - R
        args = Add.make_args(L) + Add.make_args(-R)
        if evaluate is None:
            # no cancellation, but canonical
            return _unevaluated_Add(*args)
        # no cancellation, not canonical
        return Add._from_args(args)

    def integrate(self, *args, **kwargs):
        """See the integrate function in sympy.integrals"""
        from sympy.integrals import integrate
        return integrate(self, *args, **kwargs)

    def as_poly(self, *gens, **kwargs):
        '''Returns lhs-rhs as a Poly

        Examples
        ========

        >>> from sympy import Eq
        >>> from sympy.abc import x
        >>> Eq(x**2, 1).as_poly(x)
        Poly(x**2 - 1, x, domain='ZZ')
        '''
        if self.function is not Q.eq:
            raise AttributeError("%s does not support attribute 'as_poly'" % self.function)
        return (self.lhs - self.rhs).as_poly(*gens, **kwargs)


class _DeprecatedNe(_DeprecatedRelational):
    def __new__(cls, *args, **options):
        if len(args) == 2:
            # allow rel.__class__(rel.lhs, rel.rhs)
            args = (Q.ne, args[0], args[1])
        return super().__new__(cls, *args, **options)



from .simplify import eqsimp
