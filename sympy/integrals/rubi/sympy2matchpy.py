'''
Contains Sympy to MatchPy converter function.
'''
from sympy.integrals.rubi.operation import (Int, Mul, Add, Pow, And, Or, ZeroQ, NonzeroQ, List, Log, RemoveContent, PositiveIntegerQ, NegativeIntegerQ, PositiveQ, IntegerQ, IntegersQ, PosQ, NegQ, FracPart, IntPart, RationalQ, Subst, LinearQ, Sqrt, NegativeQ, ArcCosh, Rational, Less, Not, Simplify, Denominator, Coefficient, SumSimplerQ, Equal, Unequal, SimplerQ, LessEqual, IntLinearcQ, Greater, GreaterEqual, FractionQ, ExpandIntegrand, With, Set, Hypergeometric2F1, TogetherSimplify, Inequality)
from sympy.core import sympify
from sympy.integrals.rubi.symbol import VariableSymbol, matchpyInteger

get_matchpy_node = {
    "Add": Add,
    "Mul": Mul,
    "Pow": Pow,
    "log": Log,
    'Int': Int,
    'Log': Log,
    'sqrt': Sqrt,
}

def sympy2matchpy(expr):
    '''
    Converts a SymPy expression into a MatchPy expression
    '''

    if expr == None:
        return None
    if isinstance(expr, int):
        return matchpyInteger(expr)
    elif expr.is_Atom:
        if expr.is_Number:
            return matchpyInteger(expr)
        return VariableSymbol(str(expr))
    elif type(expr).__name__ in get_matchpy_node.keys():
        matchpy_node = get_matchpy_node[type(expr).__name__]
        args = [sympy2matchpy(i) for i in expr.args]
        return matchpy_node(*args)
    else:
        raise ValueError(('Could not parse: {}').format(expr))
