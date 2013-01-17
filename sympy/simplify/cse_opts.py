""" Optimizations of the expression tree representation for better CSE
opportunities.
"""
from sympy.core import Add, Basic, Expr, Mul
from sympy.core.basic import preorder_traversal
from sympy.core.exprtools import factor_terms
from sympy.utilities.iterables import default_sort_key


class Neg(Expr):
    """ Stub to hold negated expression.
    """
    __slots__ = []


def sub_pre(e):
    """ Replace y - x with Neg(x - y) if -1 can be extracted from y - x.
    """
    reps = [a for a in e.atoms(Add) if a.could_extract_minus_sign()]

    # make it canonical
    reps.sort(key=default_sort_key)

    e = e.subs([(a, Mul(-1, -a, evaluate=False)) for a in reps])
    # now replace any persisting Adds, a, that can have -1 extracted with Neg(-a)
    if isinstance(e, Basic):
        negs = {}
        for a in sorted(e.atoms(Add), key=default_sort_key):
            if a in reps or a.could_extract_minus_sign():
                negs[a] = Neg(-a)
        e = e.xreplace(negs)
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
