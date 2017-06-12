'''
Contains all the `ReplacementRule` for Rubi integration
'''

from matchpy import Wildcard, Pattern, ReplacementRule, ManyToOneReplacer
from .operation import *
from .symbol import VariableSymbol, Integer
from .constraint import cons

a, b, c, d, e, f, g, h, x, u = map(VariableSymbol, 'abcdefghxu')
n, m = map(VariableSymbol, 'nm')


a_, b_, c_, d_, e_, f_, g_, h_ = map(Wildcard.dot, 'abcdefgh')
n_, m_ = map(Wildcard.dot, 'nm')
x_, u_ = map(Wildcard.symbol, 'xu')


def rubi_object():
    '''
    Function which compiles Rubi rules and returns `rubi` object
    '''

    rubi = ManyToOneReplacer()

    pattern1 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (d, x, c, a, b)))
    rule1 = ReplacementRule(pattern1, lambda d, x, c, a, b : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule1)

    pattern2 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (d, x, c, a)))
    rule2 = ReplacementRule(pattern2, lambda d, x, c, a : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule2)

    pattern3 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule3 = ReplacementRule(pattern3, lambda x, c, a, b : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule3)

    pattern4 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule4 = ReplacementRule(pattern4, lambda x, c, a : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), Integer(-1)), x))
    rubi.add(rule4)

    pattern5 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-2), Integer(3)))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (d, x, a, c, b)))
    rule5 = ReplacementRule(pattern5, lambda d, x, a, c, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule5)

    pattern6 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (d, x, c, b)))
    rule6 = ReplacementRule(pattern6, lambda d, x, c, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule6)

    pattern7 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-2), Integer(3)))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Add(Mul(Integer(-1), a_, d_), c_))), (d, x, a, c)))
    rule7 = ReplacementRule(pattern7, lambda d, x, a, c : With(List(Set(q, Rt(Add(Mul(Integer(-1), a, d), c), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule7)

    pattern8 = Pattern(Int(Mul(Pow(Mul(d_, x_), Rational(Integer(-2), Integer(3))), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (d, x, a, b)))
    rule8 = ReplacementRule(pattern8, lambda d, x, a, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule8)

    pattern9 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, a, c, b)))
    rule9 = ReplacementRule(pattern9, lambda x, a, c, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule9)

    pattern10 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (d, x, c)))
    rule10 = ReplacementRule(pattern10, lambda d, x, c : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule10)

    pattern11 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (d, x, b)))
    rule11 = ReplacementRule(pattern11, lambda d, x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule11)

    pattern12 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule12 = ReplacementRule(pattern12, lambda x, c, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule12)

    pattern13 = Pattern(Int(Mul(Pow(Mul(d_, x_), Rational(Integer(-2), Integer(3))), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (d, x, a)))
    rule13 = ReplacementRule(pattern13, lambda d, x, a : With(List(Set(q, Rt(Mul(Integer(-1), a, d), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule13)

    pattern14 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, a, c)))
    rule14 = ReplacementRule(pattern14, lambda x, a, c : With(List(Set(q, Rt(Add(Mul(Integer(-1), a), c), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule14)

    pattern15 = Pattern(Int(Mul(Pow(x_, Rational(Integer(-2), Integer(3))), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule15 = ReplacementRule(pattern15, lambda x, a, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule15)

    pattern16 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (d, x)))
    rule16 = ReplacementRule(pattern16, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule16)

    pattern17 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-2), Integer(3)))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule17 = ReplacementRule(pattern17, lambda x, c : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule17)

    pattern18 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Rational(Integer(-5), Integer(3)))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule18 = ReplacementRule(pattern18, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule18)

    pattern19 = Pattern(Int(Mul(Pow(x_, Rational(Integer(-2), Integer(3))), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule19 = ReplacementRule(pattern19, lambda x, a : With(List(Set(q, Rt(Mul(Integer(-1), a), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule19)

    pattern20 = Pattern(Int(Pow(x_, Rational(Integer(-5), Integer(3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule20 = ReplacementRule(pattern20, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule20)

    pattern21 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(3)))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (d, x, a, c, b)))
    rule21 = ReplacementRule(pattern21, lambda d, x, a, c, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule21)

    pattern22 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (d, x, c, b)))
    rule22 = ReplacementRule(pattern22, lambda d, x, c, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule22)

    pattern23 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(3)))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NegQ(Add(Mul(Integer(-1), a_, d_), c_))), (d, x, a, c)))
    rule23 = ReplacementRule(pattern23, lambda d, x, a, c : With(List(Set(q, Rt(Add(Mul(a, d), Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule23)

    pattern24 = Pattern(Int(Mul(Pow(Mul(d_, x_), Rational(Integer(-1), Integer(3))), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (d, x, a, b)))
    rule24 = ReplacementRule(pattern24, lambda d, x, a, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule24)

    pattern25 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, a, c, b)))
    rule25 = ReplacementRule(pattern25, lambda x, a, c, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule25)

    pattern26 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (d, x, c)))
    rule26 = ReplacementRule(pattern26, lambda d, x, c : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule26)

    pattern27 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (d, x, b)))
    rule27 = ReplacementRule(pattern27, lambda d, x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule27)

    pattern28 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule28 = ReplacementRule(pattern28, lambda x, c, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule28)

    pattern29 = Pattern(Int(Mul(Pow(Mul(d_, x_), Rational(Integer(-1), Integer(3))), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (d, x, a)))
    rule29 = ReplacementRule(pattern29, lambda d, x, a : With(List(Set(q, Rt(Mul(a, d), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule29)

    pattern30 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, a, c)))
    rule30 = ReplacementRule(pattern30, lambda x, a, c : With(List(Set(q, Rt(Add(a, Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule30)

    pattern31 = Pattern(Int(Mul(Pow(x_, Rational(Integer(-1), Integer(3))), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule31 = ReplacementRule(pattern31, lambda x, a, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule31)

    pattern32 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (d, x)))
    rule32 = ReplacementRule(pattern32, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule32)

    pattern33 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Rational(Integer(-1), Integer(3)))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule33 = ReplacementRule(pattern33, lambda x, c : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Rational(Integer(1), Integer(3))))))))
    rubi.add(rule33)

    pattern34 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Rational(Integer(-4), Integer(3)))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule34 = ReplacementRule(pattern34, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule34)

    pattern35 = Pattern(Int(Mul(Pow(x_, Rational(Integer(-1), Integer(3))), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule35 = ReplacementRule(pattern35, lambda x, a : With(List(Set(q, Rt(a, Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule35)

    pattern36 = Pattern(Int(Pow(x_, Rational(Integer(-4), Integer(3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule36 = ReplacementRule(pattern36, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Rational(Integer(3), Integer(2)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))), Mul(Rational(Integer(1), Integer(2)), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Rational(Integer(3), Integer(2)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Rational(Integer(1), Integer(3))))))))
    rubi.add(rule36)

    pattern37 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (d, x, a, c, b)))
    rule37 = ReplacementRule(pattern37, lambda d, x, a, c, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule37)

    pattern38 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (d, x, c, b)))
    rule38 = ReplacementRule(pattern38, lambda d, x, c, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule38)

    pattern39 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (d, x, a, c)))
    rule39 = ReplacementRule(pattern39, lambda d, x, a, c : Add(Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule39)

    pattern40 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (d, x, a, b)))
    rule40 = ReplacementRule(pattern40, lambda d, x, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule40)

    pattern41 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, a, c, b)))
    rule41 = ReplacementRule(pattern41, lambda x, a, c, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule41)

    pattern42 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (d, x, c)))
    rule42 = ReplacementRule(pattern42, lambda d, x, c : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule42)

    pattern43 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (d, x, b)))
    rule43 = ReplacementRule(pattern43, lambda d, x, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule43)

    pattern44 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule44 = ReplacementRule(pattern44, lambda x, c, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule44)

    pattern45 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (d, x, a)))
    rule45 = ReplacementRule(pattern45, lambda d, x, a : Add(Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule45)

    pattern46 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, a, c)))
    rule46 = ReplacementRule(pattern46, lambda x, a, c : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule46)

    pattern47 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule47 = ReplacementRule(pattern47, lambda x, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule47)

    pattern48 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (d, x)))
    rule48 = ReplacementRule(pattern48, lambda d, x : Add(Mul(zoo, d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule48)

    pattern49 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule49 = ReplacementRule(pattern49, lambda x, c : Add(Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule49)

    pattern50 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule50 = ReplacementRule(pattern50, lambda x, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule50)

    pattern51 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule51 = ReplacementRule(pattern51, lambda x, a : Add(Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule51)

    pattern52 = Pattern(Int(Pow(x_, Integer(-2)), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule52 = ReplacementRule(pattern52, lambda x : Mul(Integer(2), zoo, Int(Pow(x, Integer(-1)), x)))
    rubi.add(rule52)

    pattern53 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), x_), cons(FreeQ(List(a_, b_), x_), (a, b, x)))
    rule53 = ReplacementRule(pattern53, lambda a, b, x : Mul(Pow(b, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))))
    rubi.add(rule53)

    pattern54 = Pattern(Int(Pow(Add(a_, x_), Integer(-1)), x_), cons(FreeQ(List(a_, Integer(1)), x_), (a, x)))
    rule54 = ReplacementRule(pattern54, lambda a, x : Log(RemoveContent(Add(a, x), x)))
    rubi.add(rule54)

    pattern55 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(True, (x,)))
    rule55 = ReplacementRule(pattern55, lambda x : Log(x))
    rubi.add(rule55)

    pattern56 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (d, x, n, a, c, b, m)))
    rule56 = ReplacementRule(pattern56, lambda d, x, n, a, c, b, m : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule56)

    pattern57 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, c_, d_, m_, n_), x_)), (d, x, n, c, b, m)))
    rule57 = ReplacementRule(pattern57, lambda d, x, n, c, b, m : Add(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule57)

    pattern58 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_)), (d, x, n, a, c, m)))
    rule58 = ReplacementRule(pattern58, lambda d, x, n, a, c, m : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))
    rubi.add(rule58)

    pattern59 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, Integer(0), d_, m_, n_), x_)), (d, x, n, a, b, m)))
    rule59 = ReplacementRule(pattern59, lambda d, x, n, a, b, m : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule59)

    pattern60 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_)), (x, n, a, c, b, m)))
    rule60 = ReplacementRule(pattern60, lambda x, n, a, c, b, m : Add(Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule60)

    pattern61 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, n_), x_)), (d, x, n, c, m)))
    rule61 = ReplacementRule(pattern61, lambda d, x, n, c, m : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule61)

    pattern62 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, Integer(0), d_, m_, n_), x_)), (d, x, n, b, m)))
    rule62 = ReplacementRule(pattern62, lambda d, x, n, b, m : Add(Mul(zoo, Pow(b, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule62)

    pattern63 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, n_), x_)), (x, n, c, b, m)))
    rule63 = ReplacementRule(pattern63, lambda x, n, c, b, m : Add(Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule63)

    pattern64 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), Integer(0), d_, m_, n_), x_)), (d, x, n, a, m)))
    rule64 = ReplacementRule(pattern64, lambda d, x, n, a, m : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule64)

    pattern65 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_)), (x, n, a, c, m)))
    rule65 = ReplacementRule(pattern65, lambda x, n, a, c, m : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule65)

    pattern66 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, Integer(0), Integer(1), m_, n_), x_)), (x, n, a, b, m)))
    rule66 = ReplacementRule(pattern66, lambda x, n, a, b, m : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule66)

    pattern67 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, m_, n_), x_)), (d, x, n, m)))
    rule67 = ReplacementRule(pattern67, lambda d, x, n, m : Add(Mul(zoo, d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule67)

    pattern68 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, n_), x_)), (x, n, c, m)))
    rule68 = ReplacementRule(pattern68, lambda x, n, c, m : Add(Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule68)

    pattern69 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), m_, n_), x_)), (x, n, b, m)))
    rule69 = ReplacementRule(pattern69, lambda x, n, b, m : Add(Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(b, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Mul(b, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule69)

    pattern70 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, n, a, m)))
    rule70 = ReplacementRule(pattern70, lambda x, n, a, m : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule70)

    pattern71 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, n, m)))
    rule71 = ReplacementRule(pattern71, lambda x, n, m : Add(Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(x, Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule71)

    pattern72 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, c_, d_), x_), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (d, x, n, a, c, b, m)))
    rule72 = ReplacementRule(pattern72, lambda d, x, n, a, c, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule72)

    pattern73 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, c_, d_), x_), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_)), (d, x, n, c, b, m)))
    rule73 = ReplacementRule(pattern73, lambda d, x, n, c, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule73)

    pattern74 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_)), (d, x, n, a, c, m)))
    rule74 = ReplacementRule(pattern74, lambda d, x, n, a, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), c, Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule74)

    pattern75 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_)), (d, x, n, a, b, m)))
    rule75 = ReplacementRule(pattern75, lambda d, x, n, a, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule75)

    pattern76 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, c_, Integer(1)), x_), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, a, c, b, m)))
    rule76 = ReplacementRule(pattern76, lambda x, n, a, c, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule76)

    pattern77 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_)), (d, x, n, c, m)))
    rule77 = ReplacementRule(pattern77, lambda d, x, n, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(d, Pow(x, p))), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule77)

    pattern78 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_)), (d, x, n, b, m)))
    rule78 = ReplacementRule(pattern78, lambda d, x, n, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), d, Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule78)

    pattern79 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_)), (x, n, c, b, m)))
    rule79 = ReplacementRule(pattern79, lambda x, n, c, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule79)

    pattern80 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_)), (d, x, n, a, m)))
    rule80 = ReplacementRule(pattern80, lambda d, x, n, a, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule80)

    pattern81 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_)), (x, n, a, c, m)))
    rule81 = ReplacementRule(pattern81, lambda x, n, a, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), c, Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule81)

    pattern82 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_)), (x, n, a, b, m)))
    rule82 = ReplacementRule(pattern82, lambda x, n, a, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule82)

    pattern83 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_)), (d, x, n, m)))
    rule83 = ReplacementRule(pattern83, lambda d, x, n, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(d, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule83)

    pattern84 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_)), (x, n, c, m)))
    rule84 = ReplacementRule(pattern84, lambda x, n, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule84)

    pattern85 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_)), (x, n, b, m)))
    rule85 = ReplacementRule(pattern85, lambda x, n, b, m : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule85)

    pattern86 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_)), (x, n, a, m)))
    rule86 = ReplacementRule(pattern86, lambda x, n, a, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule86)

    pattern87 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_)), (x, n, m)))
    rule87 = ReplacementRule(pattern87, lambda x, n, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Pow(x, p), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule87)

    pattern88 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (d, x, n, a, c, b, m)))
    rule88 = ReplacementRule(pattern88, lambda d, x, n, a, c, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule88)

    pattern89 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (d, x, n, c, b, m)))
    rule89 = ReplacementRule(pattern89, lambda d, x, n, c, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule89)

    pattern90 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (d, x, n, a, c, m)))
    rule90 = ReplacementRule(pattern90, lambda d, x, n, a, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule90)

    pattern91 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (d, x, n, a, b, m)))
    rule91 = ReplacementRule(pattern91, lambda d, x, n, a, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule91)

    pattern92 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, a, c, b, m)))
    rule92 = ReplacementRule(pattern92, lambda x, n, a, c, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule92)

    pattern93 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (d, x, n, c, m)))
    rule93 = ReplacementRule(pattern93, lambda d, x, n, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule93)

    pattern94 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (d, x, n, b, m)))
    rule94 = ReplacementRule(pattern94, lambda d, x, n, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule94)

    pattern95 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, n, c, b, m)))
    rule95 = ReplacementRule(pattern95, lambda x, n, c, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule95)

    pattern96 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (d, x, n, a, m)))
    rule96 = ReplacementRule(pattern96, lambda d, x, n, a, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule96)

    pattern97 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, n, a, c, m)))
    rule97 = ReplacementRule(pattern97, lambda x, n, a, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule97)

    pattern98 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, n, a, b, m)))
    rule98 = ReplacementRule(pattern98, lambda x, n, a, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule98)

    pattern99 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (d, x, n, m)))
    rule99 = ReplacementRule(pattern99, lambda d, x, n, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule99)

    pattern100 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, n, c, m)))
    rule100 = ReplacementRule(pattern100, lambda x, n, c, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule100)

    pattern101 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, n, b, m)))
    rule101 = ReplacementRule(pattern101, lambda x, n, b, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Mul(b, x), Pow(p, Integer(-1))))))))
    rubi.add(rule101)

    pattern102 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, n, a, m)))
    rule102 = ReplacementRule(pattern102, lambda x, n, a, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule102)

    pattern103 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x, n, m)))
    rule103 = ReplacementRule(pattern103, lambda x, n, m : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Integer(1)))))
    rubi.add(rule103)

    pattern104 = Pattern(Int(Pow(Add(a_, Mul(b_, u_)), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (a, b, u, m)))
    rule104 = ReplacementRule(pattern104, lambda a, b, u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, Mul(b, x)), m), x), x, u)))
    rubi.add(rule104)

    pattern105 = Pattern(Int(Pow(Mul(b_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (b, u, m)))
    rule105 = ReplacementRule(pattern105, lambda b, u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Mul(b, x), m), x), x, u)))
    rubi.add(rule105)

    pattern106 = Pattern(Int(Pow(Add(a_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (a, u, m)))
    rule106 = ReplacementRule(pattern106, lambda a, u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, x), m), x), x, u)))
    rubi.add(rule106)

    pattern107 = Pattern(Int(Pow(u_, m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (u, m)))
    rule107 = ReplacementRule(pattern107, lambda u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, m), x), x, u)))
    rubi.add(rule107)

    pattern108 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, b_, m_), x_)), (a, b, x, m)))
    rule108 = ReplacementRule(pattern108, lambda a, b, x, m : Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule108)

    pattern109 = Pattern(Int(Pow(Mul(b_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), b_, m_), x_)), (b, x, m)))
    rule109 = ReplacementRule(pattern109, lambda b, x, m : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule109)

    pattern110 = Pattern(Int(Pow(Add(a_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, Integer(1), m_), x_)), (a, x, m)))
    rule110 = ReplacementRule(pattern110, lambda a, x, m : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule110)

    pattern111 = Pattern(Int(Pow(x_, m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), Integer(1), m_), x_)), (x, m)))
    rule111 = ReplacementRule(pattern111, lambda x, m : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule111)

    pattern112 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(m_, Integer(1)))), (x, m)))
    rule112 = ReplacementRule(pattern112, lambda x, m : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule112)

    pattern113 = Pattern(Int(x_, x_), cons(And(NonzeroQ(Integer(2)), FreeQ(Integer(1), x_)), (x,)))
    rule113 = ReplacementRule(pattern113, lambda x : Mul(Rational(Integer(1), Integer(2)), Pow(x, Integer(2))))
    rubi.add(rule113)


    return rubi
