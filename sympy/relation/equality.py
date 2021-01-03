"""
Module for mathematical equality.
"""
from sympy.assumptions import Q
from sympy.core import Equality
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

    @property
    def as_Relational(self):
        return Equality

    def _eval_relation(self, lhs, rhs):
        # logic for simple real numbers
        return lhs == rhs
