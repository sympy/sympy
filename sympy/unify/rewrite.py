""" Functions to support rewriting of SymPy expressions """

from sympy.unify.usympy import unify
from sympy.unify.usympy import rebuild
from sympy.rules.tools import subs
from sympy import Expr

def rewriterule(source, target, variables=(), condition=None):
    """ Rewrite rule

    Transform expressions that match source into expressions that match target
    treating all `variables` as wilds.

    >>> from sympy.abc import x, y, z
    >>> from sympy.unify.rewrite import rewriterule
    >>> rl = rewriterule(x + y, x**y, [x, y])
    >>> list(rl(z + 3))
    [3**z, z**3]

    Use ``condition`` to specify additional requirements.  Inputs are taken in
    the same order as is found in variables.

    >>> rl = rewriterule(x + y, x**y, [x, y], lambda x, y: x.is_integer)
    >>> list(rl(z + 3))
    [3**z]
    """

    def rewrite_rl(expr):
        for match in unify(source, expr, {}, variables=variables):
            if (condition is None or
                condition(*[match.get(var, var) for var in variables]) == True):
                expr2 = subs(match)(target)
                if isinstance(expr2, Expr):
                    expr2 = rebuild(expr2)
                yield expr2
    return rewrite_rl
