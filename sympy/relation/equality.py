"""
Module for mathematical equality.
"""
from sympy.assumptions import Q
from sympy.core import Add, Expr
from sympy.core.add import _unevaluated_Add
from sympy.solvers.solveset import linear_coeffs
from .binrel import BinaryRelation


class Equal(BinaryRelation):
    """
    Binary equality.
    """

    is_reflexive = True
    is_symmetric = True

    name = 'eq'
    str_name = latex_name = "="

    @property
    def handler(self):
        from .handlers.equality import EqualHandler
        return EqualHandler

    @property
    def reversed(self):
        return Q.eq

    def _eval_relation(self, lhs, rhs):
        # logic for simple numbers
        return lhs == rhs

    def _eval_simplify(self, **kwargs):
        lhs, rhs = self.arguments

        # standard simplify
        e = super()._eval_simplify(**kwargs)

        if not isinstance(e.function, Equal):
            return e
        if not isinstance(e.lhs, Expr) or not isinstance(e.rhs, Expr):
            return e
        free = self.free_symbols
        if len(free) == 1:
            try:
                x = free.pop()
                m, b = linear_coeffs(
                    _convert_to_Add(e, evaluate=False), x)
                if m.is_zero is False:
                    enew = e.function(x, -b / m)
                else:
                    enew = e.function(m * x, -b)
                measure = kwargs['measure']
                if measure(enew) <= kwargs['ratio'] * measure(e):
                    e = enew
            except ValueError:
                pass
        return e.canonical

def _convert_to_Add(rel, evaluate=True):
    # Mimic Relational._eval_rewrite_as_Add.
    args = rel.lhs, rel.rhs
    L, R = args
    if evaluate:
        # allow cancellation of args
        return L - R
    args = Add.make_args(L) + Add.make_args(-R)
    if evaluate is None:
        # no cancellation, but canonical
        return _unevaluated_Add(*args)
    # no cancellation, not canonical
    return Add._from_args(args)
