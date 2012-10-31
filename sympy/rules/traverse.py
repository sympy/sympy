# Strategies to traverse a Tree
from util import new, is_leaf, children

def top_down(rule):
    """ Apply a rule down an AST running it on the top nodes first """
    def top_down_rl(expr):
        newexpr = rule(expr)
        if is_leaf(newexpr):
            return newexpr
        return new(type(newexpr), *map(top_down_rl, children(newexpr)))
    return top_down_rl

def bottom_up(rule):
    """ Apply a rule down an AST running it on the bottom nodes first """
    def bottom_up_rl(expr):
        if is_leaf(expr):
            return rule(expr)
        else:
            return rule(new(type(expr), *map(bottom_up_rl, children(expr))))
    return bottom_up_rl
