'''
Contains all the `ReplacementRule` for Rubi integration
'''

from matchpy import Wildcard, Pattern, ReplacementRule, is_match, replace_all, ManyToOneReplacer
from .operation import *
from .symbol import VariableSymbol, ConstantSymbol
from .constraint import cons

a, b, c, d, e, f, g, h, x = map(VariableSymbol, 'abcdefghx')
n, m = map(VariableSymbol, 'nm')


a_, b_, c_, d_, e_, f_, g_, h_ = map(Wildcard.dot, 'abcdefgh')
n_, m_ = map(Wildcard.dot, 'nm')
x_ = Wildcard.symbol('x')
u_ = Wildcard.symbol('u')

one = Wildcard.dot('x1')

def rubi_object():
    '''
    Function which compiles Rubi rules and returns `rubi` object
    '''

    rubi = ManyToOneReplacer()

    pattern2 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(ConstantSymbol(1), m_))), (m, x)))
    rule2 = ReplacementRule(pattern2, lambda m, x : Mul(Pow(Add(ConstantSymbol(1), m), ConstantSymbol(-1)), Pow(x, Add(ConstantSymbol(1), m))))

    rubi.add(rule2)

    return rubi
