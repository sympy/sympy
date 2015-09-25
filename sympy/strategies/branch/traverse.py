""" Branching Strategies to Traverse a Tree """

from __future__ import print_function, division

from itertools import product

from sympy.strategies.util import expr_fns
from .core import chain, identity, do_one

def top_down(brule, fns=expr_fns):
    """ Apply a rule down a tree running it on the top nodes first """
    return chain(do_one(brule, identity),
                 lambda expr: sall(top_down(brule, fns), fns)(expr))

def sall(brule, fns=expr_fns):
    """ Strategic all - apply rule to args """
    op, new, children, leaf = map(fns.get, ('op', 'new', 'children', 'leaf'))
    def all_rl(expr):
        if leaf(expr):
            yield expr
        else:
            myop = op(expr)
            argss = product(*map(brule, children(expr)))
            for args in argss:
                yield new(myop, *args)
    return all_rl
