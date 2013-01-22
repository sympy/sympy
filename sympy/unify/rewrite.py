""" Functions to support rewriting of SymPy expressions """

from sympy.unify.usympy import unify
from sympy.unify.usympy import rebuild
from sympy.rules.tools import subs
from sympy import Expr

def rewriterule(p1, p2, variables=()):
    def rewrite_rl(expr):
        for match in unify(p1, expr, {}, variables=variables):
            expr2 = subs(match)(p2)
            if isinstance(expr2, Expr):
                expr2 = rebuild(expr2)
            yield expr2
    return rewrite_rl
