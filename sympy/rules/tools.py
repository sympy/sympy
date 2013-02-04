import rl
from strat import do_one
from traverse import top_down

def subs(d):
    """ Full simultaneous exact substitution

    Example
    =======

    >>> from sympy.rules.tools import subs
    >>> from sympy import Basic
    >>> mapping = {1: 4, 4: 1, Basic(5): Basic(6, 7)}
    >>> expr = Basic(1, Basic(2, 3), Basic(4, Basic(5)))
    >>> subs(mapping)(expr)
    Basic(4, Basic(2, 3), Basic(1, Basic(6, 7)))
    """
    if d:
        return top_down(do_one(*map(rl.subs, *zip(*d.items()))))
    else:
        return lambda x: x
