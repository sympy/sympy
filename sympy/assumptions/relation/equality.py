"""
Module for mathematical equality.
"""
from sympy.assumptions import Q
from sympy.multipledispatch import Dispatcher

from .binrel import BinaryRelation


class Equal(BinaryRelation):
    """
    Binary equality.
    """

    is_reflexive = True
    is_symmetric = True

    name = 'eq'
    str_name = latex_name = "="

    handler = Dispatcher(
        "EqualHandler",
        doc="Handler for Q.eq. Test that two expressions are equal."
    )

    @property
    def reversed(self):
        return Q.eq

    def _eval_relation(self, lhs, rhs):
        # logic for simple numbers
        return lhs == rhs

    def _simplify_applied(self, lhs, rhs, **kwargs):
        return eqsimp(self(lhs, rhs), **kwargs)


from .simplify import eqsimp
