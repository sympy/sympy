"""
Module for equation solving with respect to certain variable.
"""

from sympy import FiniteSet, Or, Q, S, solveset


# solveeqn

def solveeqn(eqn, symbol=None, domain=S.Complexes):
    """
    Rewrite the equation with respect to a variable. This function
    just runs solveset and convert the result to equation form.
    """
    if symbol is None:
        free_symbols = eqn.free_symbols
        if len(free_symbols) == 1:
            symbol = free_symbols.pop()
        else:
            raise ValueError("symbol must be specified.")
    rel = eqn.as_Relational()
    solve_result = solveset(rel, symbol, domain)
    if isinstance(solve_result, FiniteSet):
        rels = [Q.eq(symbol, elem) for elem in solve_result]
        return Or(*rels)
    raise NotImplementedError
