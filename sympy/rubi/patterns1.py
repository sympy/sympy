
from matchpy import Wildcard, Pattern, ReplacementRule, ManyToOneReplacer, is_match
from .operation import *
from .symbol import VariableSymbol, Integer
from .constraint import cons
from .sympy2matchpy import sympy2matchpy
from sympy import sympify

a, b, c, d, e, f, g, h, x, u = map(VariableSymbol, 'abcdefghxu')
n, m = map(VariableSymbol, 'nm')


a_, b_, c_, d_, e_, f_, g_, h_ = map(Wildcard.dot, 'abcdefgh')
n_, m_ = map(Wildcard.dot, 'nm')
x_, u_ = map(Wildcard.symbol, 'xu')


def rubi_object():
    rubi = ManyToOneReplacer()

    #pattern533 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, Mul(b_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, a, n, b)))
    pattern533 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(x_, n_)), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_)), (x, a, n, b)))
    pattern533 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(x_, n_)), x_))
    rule533 = ReplacementRule(pattern533, lambda x, a, n, b : a)
    rubi.add(rule533)

    #sub = Int(Mul(Add(Mul(b, x), a), Pow(x, Integer(3))), x)
    #print(is_match(sub, pattern533))
    #print()
    #print()

    return rubi
