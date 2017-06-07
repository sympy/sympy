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
x_, u_ = map(Wildcard.symbol, 'xu')


def rubi_object():
    '''
    Function which compiles Rubi rules and returns `rubi` object
    '''

    rubi = ManyToOneReplacer()

    pattern0 = Pattern(Int(Pow(Add(Mul(b_, x_), a_), ConstantSymbol(-1)), x_), cons(FreeQ(List(a_, b_), x_), (b, a, x)))
    rule0 = ReplacementRule(pattern0, lambda b, a, x : Mul(log(RemoveContent(Add(Mul(b, x), a), x)), Pow(b, ConstantSymbol(-1))))

    pattern1 = Pattern(Int(Pow(Add(Mul(b_, x_), a_), m_), x_), cons(And(NonzeroQ(Add(ConstantSymbol(1), m_)), FreeQ(List(a_, b_, m_), x_)), (a, b, x, m)))
    rule1 = ReplacementRule(pattern1, lambda a, b, x, m : Mul(Pow(Add(ConstantSymbol(1), m), ConstantSymbol(-1)), Pow(Add(Mul(b, x), a), Add(ConstantSymbol(1), m)), Pow(b, ConstantSymbol(-1))))

    pattern2 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(ConstantSymbol(1), m_))), (m, x)))
    rule2 = ReplacementRule(pattern2, lambda m, x : Mul(Pow(Add(ConstantSymbol(1), m), ConstantSymbol(-1)), Pow(x, Add(ConstantSymbol(1), m))))

    rubi.add(rule0)
    rubi.add(rule1)
    rubi.add(rule2)

    return rubi
