from sympy.unify.unify_sympy import unify
from sympy.unify.unify_sympy import rebuild
from sympy.rules.tools import subs
from sympy import Expr

def rewriterule(p1, p2):
    def rewrite_rl(expr):
        match = unify(p1, expr, {})
        for m in match:
            expr2 = subs(m)(p2)
            if isinstance(expr2, Expr):
                expr2 = rebuild(expr2)
            yield expr2
    return rewrite_rl
