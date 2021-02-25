"""
Supports old relational classes.
"""

from sympy.assumptions import Q
from .binrel import _DeprecatedRelational


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
