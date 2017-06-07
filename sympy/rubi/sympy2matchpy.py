'''
Contains Sympy to MatchPy converter function.
'''

import matchpy
from .operation import *
from sympy.core import sympify
from .symbol import VariableSymbol, ConstantSymbol


get_matchpy_node = {
    "Add": Add,
    "Mul": Mul,
    "Pow": Pow,
    "log": Log,
    'Int': Int,
    'ZeroQ': ZeroQ,
    'NonzeroQ': NonzeroQ,
    'FreeQ': FreeQ,
    'And': And,
    'Or': Or,
    'RemoveContent': RemoveContent,
    'Log': Log,
    'List': List,
}


def sympy2matchpy(expr, parse_rules=False):
    '''
    Converts a SymPy expression into a MatchPy expression
    `parse_rules` is used internally to parse rules.
    '''
    if expr.is_Atom:
        if expr.is_Number:
            if parse_rules:
                return ConstantSymbol_parse(expr)
            else:
                return ConstantSymbol(expr)
        return matchpy.Symbol(str(expr))
    elif type(expr).__name__ in get_matchpy_node.keys():
        matchpy_node = get_matchpy_node[type(expr).__name__]
        args = [sympy2matchpy(i, parse_rules) for i in expr.args]
        return matchpy_node(*args)
    else:
        raise ValueError(('Could not parse: {}').format(expr))
