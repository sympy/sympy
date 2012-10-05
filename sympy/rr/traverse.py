# Strategies to traverse a SymPy Tree
from sympy import Basic

def top_down(rule):
    def top_down_rl(expr):
        newexpr = rule(expr)
        if newexpr.is_Atom:
            return newexpr
        return Basic.__new__(newexpr.__class__, *map(top_down_rl, newexpr.args))
    return top_down_rl

def bottom_up(rule):
    def bottom_up_rl(expr):
        if expr.is_Atom:
            return rule(expr)
        else:
            return rule(Basic.__new__(expr.__class__,
                                      *map(bottom_up_rl, expr.args)))
    return bottom_up_rl
