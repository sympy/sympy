"""
Module for mathematical equality [1] and inequation [2].

References
==========

.. [1] https://en.wikipedia.org/wiki/Equality_(mathematics)
.. [2] https://en.wikipedia.org/wiki/Inequation
"""
from sympy import S
from sympy.assumptions import ask, Q
from sympy.core import Add, Expr
from sympy.core.logic import fuzzy_bool, fuzzy_and, fuzzy_xor, fuzzy_not
from sympy.core.sympify import _sympify
from sympy.functions import arg
from sympy.multipledispatch import Dispatcher
from sympy.simplify.simplify import clear_coefficients
from sympy.utilities.iterables import sift

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
    handler = Dispatcher("EqualityHandler", doc="Handler for key 'eq'.")

    @property
    def negated(self):
        return Q.ne

    def eval(self, args, assumptions=True):

        # NaN never equals to anything
        if S.NaN in args:
            return False

        # here, _eval_Eq is only called for backwards compatibility
        # new code should use is_eq with multiple dispatch as
        # outlined in the docstring
        lhs, rhs = args
        for side1, side2 in (lhs, rhs), (rhs, lhs):
            eval_func = getattr(side1, '_eval_Eq', None)
            if eval_func is not None:
                retval = eval_func(side2)
                if retval is not None:
                    return retval

        # evaulate by multipledispatch. if lhs or rhs is nonnumeric, this step
        # must return evaluated result.
        ret = super().eval(args, assumptions)
        if ret is None:
            # evaluation of numeric arguments
            ret = is_eq(*args, assumptions)
        return ret

    def _simplify_applied(self, lhs, rhs, **kwargs):
        return eqsimp(self(lhs, rhs), **kwargs)


def is_eq(lhs, rhs, assumptions):
    """
    Compare numeric *lhs* and *rhs*. Used by ``EqualityPredicate.eval()``.
    """
    # This is not dispatched as (Expr, Expr) to handler because too many
    # unambiguity is raised.

    if lhs.is_infinite or rhs.is_infinite:
        if fuzzy_xor([lhs.is_infinite, rhs.is_infinite]):
            return False
        if fuzzy_xor([lhs.is_extended_real, rhs.is_extended_real]):
            return False
        if fuzzy_and([lhs.is_extended_real, rhs.is_extended_real]):
            return fuzzy_xor([lhs.is_extended_positive, fuzzy_not(rhs.is_extended_positive)])

        # Try to split real/imaginary parts and equate them
        I = S.ImaginaryUnit

        def split_real_imag(expr):
            real_imag = lambda t: (
                'real' if t.is_extended_real else
                'imag' if (I * t).is_extended_real else None)
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
        # Guard against Q.eq(nan, nan) -> Falsesymp
        if not (arglhs == S.NaN and argrhs == S.NaN):
            return fuzzy_bool(ask(Q.eq(arglhs, argrhs), assumptions))

    if all(isinstance(i, Expr) for i in (lhs, rhs)):
        # see if the difference evaluates
        dif = lhs - rhs
        z = dif.is_zero
        if z is not None:
            if z is False and dif.is_commutative:  # issue 10728
                return False
            if z:
                return True

        n2 = _n2(lhs, rhs)
        if n2 is not None:
            return _sympify(n2 == 0)

        # see if the ratio evaluates
        n, d = dif.as_numer_denom()
        rv = None
        if n.is_zero:
            rv = d.is_nonzero
        elif n.is_finite:
            if d.is_infinite:
                rv = True
            elif n.is_zero is False:
                rv = d.is_infinite
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
        elif any(a.is_infinite for a in Add.make_args(n)):
            # (inf or nan)/x != 0
            rv = False
        if rv is not None:
            return rv

    return None

def _n2(a, b):
    """Return (a - b).evalf(2) if a and b are comparable, else None.
    This should only be used when a and b are already sympified.
    """
    # /!\ it is very important (see issue 8245) not to
    # use a re-evaluated number in the calculation of dif
    if a.is_comparable and b.is_comparable:
        dif = (a - b).evalf(2)
        if dif.is_comparable:
            return dif


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
    handler = Dispatcher("UnequalityHandler", doc="Handler for key 'ne'.")

    @property
    def negated(self):
        return Q.eq

    def _simplify_applied(self, lhs, rhs, **kwargs):
        eq = Q.eq(lhs, rhs).simplify(**kwargs)
        return self(*eq.arguments)

from .simplify import eqsimp
