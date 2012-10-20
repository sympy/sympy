from strat_pure import (exhaust, memoize, condition, chain, debug, try_safe,
        null_safe, do_one)
from traverse import bottom_up, top_down

def canon(*rules):
    """ Strategy for canonicalization

    Apply each rule in a bottom_up fashion through the tree.
    Do each one in turn.
    Keep doing this until there is no change.
    """
    return exhaust(chain(*map(bottom_up, rules)))

def typed(ruletypes):
    """ Apply rules based on the expression type

    inputs:
        ruletypes -- a dict mapping {Type: rule}

    >>> from sympy.rules import rm_id, typed
    >>> from sympy import Add, Mul
    >>> rm_zeros = rm_id(lambda x: x==0)
    >>> rm_ones  = rm_id(lambda x: x==1)
    >>> remove_idents = typed({Add: rm_zeros, Mul: rm_ones})
    """
    def typed_rl(expr):
        rl = ruletypes.get(type(expr), lambda x:x)
        return rl(expr)
    return typed_rl
