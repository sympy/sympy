""" Strategies to Traverse a Tree """
from util import new, is_leaf
from sympy.strategies.core import chain, do_one

def top_down(rule):
    """ Apply a rule down a tree running it on the top nodes first """
    return chain(rule, lambda expr: sall(top_down(rule))(expr))

def bottom_up(rule):
    """ Apply a rule down a tree running it on the bottom nodes first """
    return chain(lambda expr: sall(bottom_up(rule))(expr), rule)

def top_down_once(rule):
    """ Apply a rule down a tree - stop on success """
    return do_one(rule, lambda expr: sall(top_down(rule))(expr))

def bottom_up_once(rule):
    """ Apply a rule up a tree - stop on success """
    return do_one(lambda expr: sall(bottom_up(rule))(expr), rule)

def sall(rule):
    """ Strategic all - apply rule to args """
    def all_rl(expr):
        if is_leaf(expr):
            return expr
        else:
            args = map(rule, expr.args)
            return new(type(expr), *args)
    return all_rl
