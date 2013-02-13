from strat_pure import (exhaust, memoize, condition, chain, debug, tryit,
        null_safe, do_one, switch, minimize)
from traverse import bottom_up, top_down

def canon(*rules):
    """ Strategy for canonicalization

    Apply each rule in a bottom_up fashion through the tree.
    Do each one in turn.
    Keep doing this until there is no change.
    """
    return exhaust(top_down(exhaust(do_one(*rules))))

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
    return switch(type, ruletypes)
