'''
Contains all the `ReplacementRule` for Rubi integration
'''

from matchpy import Wildcard, Pattern, ReplacementRule, is_match, replace_all, ManyToOneReplacer
from operation import Int, Mul, Add, Pow, Log, FreeQ, NonzeroQ, List
from symbol import VariableSymbol, ConstantSymbol
from constraint import cons

a, b, c, d, e, f, g, h, x = map(VariableSymbol, 'abcdefghx')
n, m = map(VariableSymbol, 'nm')
a_, b_, c_, d_, e_, f_, g_, h_ = map(Wildcard.dot, 'abcdefgh')
n_, m_ = map(Wildcard.dot, 'nm')

x_ = Wildcard.symbol('x')
u_ = Wildcard.symbol('u')

one = ConstantSymbol(1)
m_one = ConstantSymbol(-1)

def rubi_object():
    '''
    Function which compiles Rubi rules and returns `rubi` object
    '''
    pattern1 = Pattern(Int(Mul(one, Pow(Add(a_, Mul(b_, x_)), m_one)), x_), cons(FreeQ(List(a_, b_), x_), (b, a, x)))
    rule1 = ReplacementRule(pattern1, lambda a, x: Mul(a, Log(x)))

    #pattern2 = Pattern(Int(Mul(1, Pow(Add(a_, Mul(b_, x_)), -1)), x_), cons(FreeQ(List(a, b), x), (b, a, x)))
    #pattern2 = Pattern(Int(Pow(x_, m_), x), FreeQ((m,), x), NonzeroQ(Add(m_, one)))
    #rule2 = ReplacementRule(pattern2, lambda m, x: Mul(Pow(x, Add(m, one)), Pow(Add(m, one), m_one)))

    rubi = ManyToOneReplacer(rule1)
    #rubi.add(rule2)

    return rubi
