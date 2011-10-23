""" Optimizations of the expression tree representation for better CSE
opportunities.
"""
from sympy.core import Add, Mul, Expr
from sympy.utilities.iterables import preorder_traversal

class Sub(Expr):
    """ Stub of a Sub operator to replace Add(x, Mul(NegativeOne(-1), y)).
    """
    __slots__ = []

def sub_pre(e):
    """ Replace Add(x, Mul(NegativeOne(-1), y)) with Sub(x, y).
    """
    replacements = []
    for node in preorder_traversal(e):
        if node.is_Add:
            positives = []
            negatives = []
            for arg in node.args:
                if arg.is_Mul:
                    a, b = arg.as_two_terms()
                    if (a.is_number and a.is_negative):
                        negatives.append(Mul(-a, b))
                        continue
                positives.append(arg)
            if len(negatives) > 0:
                replacement = Sub(Add(*positives), Add(*negatives))
                replacements.append((node, replacement))
    for node, replacement in replacements:
        e = e.subs(node, replacement)

    return e

def sub_post(e):
    """ Replace Sub(x,y) with the canonical form Add(x, Mul(NegativeOne(-1), y)).
    """
    replacements = []
    for node in preorder_traversal(e):
        if isinstance(node, Sub):
            replacements.append((node, Add(node.args[0], Mul(-1, node.args[1]))))
    for node, replacement in replacements:
        e = e.subs(node, replacement)

    return e

default_optimizations = [
    (sub_pre, sub_post),
]
