""" Optimizations of the expression tree representation for better CSE
opportunities.
"""

from sympy.core.basic import Basic
from sympy.core.operations import AssocOp
from sympy.utilities.iterables import preorder_traversal

from sympy import Add, Mul


def assumed(e, name):
    """ Return True if the given assumption is true about the sympy expression.

    Examples
    --------
    >>> from sympy import symbols
    >>> from sympy.simplify.cse_opts import assumed
    >>> from sympy.abc import x, y
    >>> assumed(x+y, 'is_Add')
    True
    >>> assumed(x+y, 'is_Mul')
    False

    """
    return getattr(e, name, False)

class Sub(AssocOp):
    """ Stub of a Sub operator to replace Add(x, Mul(NegativeOne(-1), y)).
    """
    __slots__ = []
    is_Add = False
    is_Sub = True

    def _eval_subs(self, old, new):
        if self==old:
            return new
        else:
            return self.__class__(*[s._eval_subs(old, new) for s in self.args ])

def sub_pre(e):
    """ Replace Add(x, Mul(NegativeOne(-1), y)) with Sub(x, y).
    """
    replacements = []
    for node in preorder_traversal(e):
        if assumed(node, 'is_Add'):
            positives = []
            negatives = []
            for arg in node.args:
                if (assumed(arg, 'is_Mul') and
                    assumed(arg.args[0], 'is_number') and
                    assumed(arg.args[0], 'is_negative')):
                    negatives.append(Mul(-arg.args[0], *arg.args[1:]))
                else:
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
        if assumed(node, 'is_Sub'):
            replacements.append((node, Add(node.args[0], Mul(-1, node.args[1]))))
    for node, replacement in replacements:
        e = e.subs(node, replacement)

    return e

default_optimizations = [
    (sub_pre, sub_post),
]
