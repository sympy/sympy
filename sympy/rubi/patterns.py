'''
Contains all the `ReplacementRule` for Rubi integration
'''

from matchpy import Wildcard, Pattern, ReplacementRule, ManyToOneReplacer
from .operation import *
from .symbol import VariableSymbol, Integer
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

    pattern1 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), x_), cons(FreeQ(List(a_, b_), x_), (b, a, x)))
    rule1 = ReplacementRule(pattern1, lambda b, a, x : Mul(Pow(b, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))))
    rubi.add(rule1)

    pattern2 = Pattern(Int(Pow(Add(a_, x_), Integer(-1)), x_), cons(FreeQ(List(a_, Integer(1)), x_), (a, x)))
    rule2 = ReplacementRule(pattern2, lambda a, x : Log(RemoveContent(Add(a, x), x)))
    rubi.add(rule2)

    pattern3 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, b_, m_), x_)), (a, b, x, m)))
    rule3 = ReplacementRule(pattern3, lambda a, b, x, m : Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule3)

    pattern4 = Pattern(Int(Pow(Mul(b_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), b_, m_), x_)), (b, x, m)))
    rule4 = ReplacementRule(pattern4, lambda b, x, m : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule4)

    pattern5 = Pattern(Int(Pow(Add(a_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, Integer(1), m_), x_)), (a, x, m)))
    rule5 = ReplacementRule(pattern5, lambda a, x, m : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule5)

    pattern6 = Pattern(Int(Pow(x_, m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), Integer(1), m_), x_)), (x, m)))
    rule6 = ReplacementRule(pattern6, lambda x, m : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule6)

    pattern7 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(m_, Integer(1)))), (m, x)))
    rule7 = ReplacementRule(pattern7, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule7)

    pattern8 = Pattern(Int(x_, x_), cons(And(NonzeroQ(Integer(2)), FreeQ(Integer(1), x_)), (x,)))
    rule8 = ReplacementRule(pattern8, lambda x : Mul(Rational(1, 2), Pow(x, Integer(2))))
    rubi.add(rule8)

    return rubi
