'''
Contains Sympy to MatchPy converter function.
'''
from sympy.integrals.rubi.operation import *
from sympy.core import sympify
from sympy.integrals.rubi.symbol import VariableSymbol, Integer

get_matchpy_node = {
    "Add": Add,
    "Mul": Mul,
    "Pow": Pow,
    "log": Log,
    'Int': Int,
    'Log': Log,
    'sqrt': Sqrt
}

def sympy2matchpy(expr):
    '''
    Converts a SymPy expression into a MatchPy expression
    '''

    if expr == None:
        return None
    if isinstance(expr, int):
        return Integer(expr)
    elif expr.is_Atom:
        if expr.is_Number:
            return Integer(expr)
        return VariableSymbol(str(expr))
    elif type(expr).__name__ in get_matchpy_node.keys():
        matchpy_node = get_matchpy_node[type(expr).__name__]
        args = [sympy2matchpy(i) for i in expr.args]
        return matchpy_node(*args)
    else:
        raise ValueError(('Could not parse: {}').format(expr))
