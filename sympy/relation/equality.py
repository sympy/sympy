"""
Module for mathematical equality.
"""
from sympy.assumptions import ask, Q
from sympy.core import Add, Equality, Expr, S
from sympy.core.logic import fuzzy_and, fuzzy_bool, fuzzy_xor, fuzzy_not
from sympy.core.relational import _n2
from sympy.functions import arg
from sympy.logic.boolalg import Boolean, BooleanAtom
from sympy.simplify.simplify import clear_coefficients
from sympy.utilities.iterables import sift
from .binrel import BinaryRelation, AppliedBinaryRelation
from .relop import relop_add, relop_mul, relop_pow


class Equal(BinaryRelation):
    """
    Binary equality.
    """

    is_reflexive = True

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


@relop_add.register(Equal, Equal)
def eq_add(rel1, rel2, assumptions=True):
    lhs = rel1.lhs + rel2.lhs
    rhs = rel1.rhs + rel2.rhs
    return Q.eq(lhs, rhs)

@relop_add.register(Equal, Expr)
def eq_expr_add(arg1, arg2, assumptions=True):
    if isinstance(arg1, AppliedBinaryRelation):
        lhs = arg1.lhs + arg2
        rhs = arg1.rhs + arg2
    return Q.eq(lhs, rhs)

@relop_add.register(Expr, Equal)
def expr_eq_add(arg1, arg2, assumptions=True):
    lhs = arg1 + arg2.lhs
    rhs = arg1 + arg2.rhs
    return Q.eq(lhs, rhs)

@relop_mul.register(Equal, Equal)
def eq_mul(rel1, rel2, assumptions=True):
    lhs = rel1.lhs*rel2.lhs
    rhs = rel1.rhs*rel2.rhs
    return Q.eq(lhs, rhs)

@relop_mul.register(Equal, Expr)
def eq_expr_mul(arg1, arg2, assumptions=True):
    lhs = arg1.lhs*arg2
    rhs = arg1.rhs*arg2
    return Q.eq(lhs, rhs)

@relop_mul.register(Expr, Equal)
def expr_eq_mul(arg1, arg2, assumptions=True):
    lhs = arg1*arg2.lhs
    rhs = arg1*arg2.rhs
    return Q.eq(lhs, rhs)

@relop_pow.register(Equal, Equal)
def eq_pow(rel1, rel2, assumptions=True):
    lhs = rel1.lhs**rel2.lhs
    rhs = rel1.rhs**rel2.rhs
    return Q.eq(lhs, rhs)

@relop_pow.register(Equal, Expr)
def eq_expr_pow(arg1, arg2, assumptions=True):
    lhs = arg1.lhs**arg2
    rhs = arg1.rhs**arg2
    return Q.eq(lhs, rhs)

@relop_pow.register(Expr, Equal)
def eq_expr_pow(arg1, arg2, assumptions=True):
    lhs = arg1**arg2.lhs
    rhs = arg1**arg2.rhs
    return Q.eq(lhs, rhs)
