""" Branching Strategies to Traverse a Tree """

from sympy.rules.util import new, is_leaf, children
from strat_pure import notempty
from itertools import product

def top_down(brule):
    """ Apply a rule down a tree running it on the top nodes first """
    def top_down_rl(expr):
        brl = notempty(brule)
        for newexpr in brl(expr):
            if is_leaf(newexpr):
                yield newexpr
            else:
                for args in product(*map(top_down_rl, children(newexpr))):
                    yield new(type(newexpr), *args)
    return top_down_rl
