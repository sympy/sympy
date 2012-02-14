""" Optimizations of the expression tree representation for better CSE
opportunities.
"""
from sympy.core import Add, Mul, Expr, S
from sympy.core.exprtools import factor_terms
from sympy.utilities.iterables import preorder_traversal

class Neg(Expr):
    """ Stub to hold negated expression.
    """
    __slots__ = []

def sub_pre(e):
    """ Replace Add(x, Mul(NegativeOne(-1), y)) with Sub(x, y).
    """
    # make canonical, first
    adds = {}
    for a in e.atoms(Add):
        adds[a] = a.could_extract_minus_sign()
    e = e.subs([(a, Mul(-1, -a, evaluate=False)
                    if adds[a] else a) for a in adds])
    # now replace any persisting Adds, a, that can have -1 extracted with Neg(-a)
    reps = dict([(a, Neg(-a)) for a in e.atoms(Add)
           if adds.get(a, a.could_extract_minus_sign())])
    e = e.xreplace(reps)
    return e

def sub_post(e):
    """ Replace Neg(x) with -x.
    """
    replacements = []
    for node in preorder_traversal(e):
        if isinstance(node, Neg):
            replacements.append((node, -node.args[0]))
    for node, replacement in replacements:
        e = e.xreplace({node: replacement})

    return e

default_optimizations = [
    (sub_pre, sub_post),
    (factor_terms, None),
]
