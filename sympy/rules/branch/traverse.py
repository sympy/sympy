""" Branching Strategies to Traverse a Tree """

from sympy.rules.util import new, is_leaf
from strat_pure import chain, identity, do_one
from sympy.core.compatibility import product

def top_down(brule):
    """ Apply a rule down a tree running it on the top nodes first """
    return chain(do_one(brule, identity), lambda expr: sall(top_down(brule))(expr))

def sall(brule):
    """ Strategic all - apply rule to args """
    def all_rl(expr):
        if is_leaf(expr):
            yield expr
        else:
            op = type(expr)
            argss = product(*map(brule, expr.args))
            for args in argss:
                yield new(op, *args)
    return all_rl
