'''
Contains Sympy to MatchPy converter function.
'''

import matchpy
from operation import *
from sympy import sympify

get_matchpy_node = {
    "Add": Add,
    "Mul": Mul,
    "Pow": Pow,
    "log": Log,
    'Int': Int,
}

get_matchpy_function = {
    'NonzeroQ': NonzeroQ,
    'FreeQ': FreeQ,
    'And': And,
    'Or': Or,
    'RemoveContent': RemoveContent,
    'Log': Log_parse,
    'List': List_parse
}

class ConstantSymbol(matchpy.Symbol):
    def __init__(self, value):
        super(self.__class__, self).__init__(str(value))
        self.value = value


def sympy2matchpy(expr):
    '''
    Converts a SymPy expression into a MatchPy expression
    '''
    print(expr)
    if expr.is_Atom:
        if expr.is_Number:
            return ConstantSymbol_parse(str(expr))
        return matchpy.Symbol(str(expr))
    elif type(expr).__name__ in get_matchpy_node.keys():
        matchpy_node = get_matchpy_node[type(expr).__name__]
        args = [sympy2matchpy(i) for i in expr.args]
        return matchpy_node(*args)
    elif type(expr).__name__ in get_matchpy_function.keys():
        matchpy_node = get_matchpy_function[type(expr).__name__]
        args = [sympy2matchpy(i) for i in expr.args]
        return matchpy_node(*args)
    else:
        raise ValueError(('Could not parse: {}').format(expr))
