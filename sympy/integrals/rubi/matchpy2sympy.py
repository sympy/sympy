from sympy.integrals.rubi.symbol import VariableSymbol, Integer
from sympy.integrals.rubi.utility_function import (Int, ZeroQ, NonzeroQ, List, Log, RemoveContent,
    PositiveIntegerQ, NegativeIntegerQ, PositiveQ, IntegerQ, PosQ, NegQ, FracPart,
    IntPart, RationalQ, Subst, LinearQ, Sqrt, NegativeQ, ArcCosh, RationalQ, Less,
    Not, Simplify, Denominator, Coefficient, SumSimplerQ, Equal, Unequal, SimplerQ,
    LessEqual, IntLinearcQ, Greater, GreaterEqual, FractionQ, ExpandIntegrand,
    With, Set, Hypergeometric2F1, TogetherSimplify)
from sympy import Add, Mul, Pow, Symbol, S, I, And, Or

m2s_function = {
    'Mul': Mul,
    'Pow': Pow,
    'Add': Add,
    'And': And,
    'Or': Or,
    'Int':  Int,
    'ZeroQ':  ZeroQ,
    'NonzeroQ':  NonzeroQ,
    'List':  List,
    'Log':  Log,
    'RemoveContent':  RemoveContent,
    'PositiveIntegerQ':  PositiveIntegerQ,
    'NegativeIntegerQ':  NegativeIntegerQ,
    'PositiveQ':  PositiveQ,
    'IntegerQ':  IntegerQ,
    'PosQ':  PosQ,
    'NegQ':  NegQ,
    'FracPart':  FracPart,
    'IntPart':  IntPart,
    'RationalQ':  RationalQ,
    'Subst':  Subst,
    'LinearQ':  LinearQ,
    'Sqrt':  Sqrt,
    'NegativeQ':  NegativeQ,
    'ArcCosh':  ArcCosh,
    'RationalQ':  RationalQ,
    'Less':  Less,
    'Not':  Not,
    'Simplify':  Simplify,
    'Denominator':  Denominator,
    'Coefficient':  Coefficient,
    'SumSimplerQ':  SumSimplerQ,
    'Equal':  Equal,
    'Unequal':  Unequal,
    'SimplerQ':  SimplerQ,
    'LessEqual':  LessEqual,
    'IntLinearcQ':  IntLinearcQ,
    'Greater':  Greater,
    'GreaterEqual':  GreaterEqual,
    'FractionQ':  FractionQ,
    'ExpandIntegrand':  ExpandIntegrand,
    'With':  With,
    'Set':  Set,
    'Hypergeometric2F1':  Hypergeometric2F1,
    'TogetherSimplify':  TogetherSimplify,
}

def matchpy2sympy(expr):
    if isinstance(expr, VariableSymbol):
        name = expr.name
        if name == 'I':
            return I
        return Symbol(name, real=True, imaginary=False)
    elif isinstance(expr, Integer):
        return S(expr.name)
    args = [matchpy2sympy(i) for i in expr.operands]
    return m2s_function[expr.name](*args)
