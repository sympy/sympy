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

    pattern1 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-9), Integer(4))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(4)))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (b, x, c, a, d)))
    rule1 = ReplacementRule(pattern1, lambda b, x, c, a, d : Add(Mul(Integer(-1), Rational(Integer(1), Integer(5)), Pow(b, Integer(-1)), d, Int(Mul(Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Rational(Integer(4), Integer(5)), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule1)

    pattern2 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-9), Integer(4))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (x, c, a, d)))
    rule2 = ReplacementRule(pattern2, lambda x, c, a, d : Add(Mul(Integer(-1), Rational(Integer(1), Integer(5)), d, Int(Mul(Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Rational(Integer(4), Integer(5)), Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule2)

    pattern3 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-9), Integer(4))), Pow(Add(c_, x_), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (b, x, c, a)))
    rule3 = ReplacementRule(pattern3, lambda b, x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(5)), Pow(b, Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Rational(Integer(4), Integer(5)), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule3)

    pattern4 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-9), Integer(4))), Pow(Add(c_, x_), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (x, c, a)))
    rule4 = ReplacementRule(pattern4, lambda x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(5)), Int(Mul(Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Rational(Integer(4), Integer(5)), Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule4)

    pattern5 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-5), Integer(4))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(4)))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (b, x, c, a, d)))
    rule5 = ReplacementRule(pattern5, lambda b, x, c, a, d : Add(Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Rational(Integer(-1), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule5)

    pattern6 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-5), Integer(4))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (x, c, a, d)))
    rule6 = ReplacementRule(pattern6, lambda x, c, a, d : Add(Mul(Rational(Integer(1), Integer(2)), Add(Mul(Integer(-1), a, d), c), Int(Mul(Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Rational(Integer(-1), Integer(4))), Pow(Add(c, Mul(d, x)), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule6)

    pattern7 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-5), Integer(4))), Pow(Add(c_, x_), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (b, x, c, a)))
    rule7 = ReplacementRule(pattern7, lambda b, x, c, a : Add(Mul(Rational(Integer(1), Integer(2)), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Rational(Integer(-1), Integer(4))), Pow(Add(c, x), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule7)

    pattern8 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-5), Integer(4))), Pow(Add(c_, x_), Rational(Integer(-1), Integer(4)))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (x, c, a)))
    rule8 = ReplacementRule(pattern8, lambda x, c, a : Add(Mul(Rational(Integer(1), Integer(2)), Add(Mul(Integer(-1), a), c), Int(Mul(Pow(Add(a, x), Rational(Integer(-5), Integer(4))), Pow(Add(c, x), Rational(Integer(-5), Integer(4)))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Rational(Integer(-1), Integer(4))), Pow(Add(c, x), Rational(Integer(-1), Integer(4))))))
    rubi.add(rule8)

    pattern9 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-3), Integer(2))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-3), Integer(2)))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (b, x, c, a, d)))
    rule9 = ReplacementRule(pattern9, lambda b, x, c, a, d : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule9)

    pattern10 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-3), Integer(2))), Pow(Add(c_, Mul(d_, x_)), Rational(Integer(-3), Integer(2)))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule10 = ReplacementRule(pattern10, lambda x, c, a, d : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule10)

    pattern11 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Rational(Integer(-3), Integer(2))), Pow(Add(c_, x_), Rational(Integer(-3), Integer(2)))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (b, x, c, a)))
    rule11 = ReplacementRule(pattern11, lambda b, x, c, a : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule11)

    pattern12 = Pattern(Int(Mul(Pow(Add(a_, x_), Rational(Integer(-3), Integer(2))), Pow(Add(c_, x_), Rational(Integer(-3), Integer(2)))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule12 = ReplacementRule(pattern12, lambda x, c, a : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule12)

    pattern13 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, Mul(b_, x_))), Integer(-1)), Pow(Sqrt(Add(c_, Mul(d_, x_))), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (b, x, c, a, d)))
    rule13 = ReplacementRule(pattern13, lambda b, x, c, a, d : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule13)

    pattern14 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, x_)), Integer(-1)), Pow(Sqrt(Add(c_, Mul(d_, x_))), Integer(-1))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule14 = ReplacementRule(pattern14, lambda x, c, a, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule14)

    pattern15 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, Mul(b_, x_))), Integer(-1)), Pow(Sqrt(Add(c_, x_)), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (b, x, c, a)))
    rule15 = ReplacementRule(pattern15, lambda b, x, c, a : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule15)

    pattern16 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, x_)), Integer(-1)), Pow(Sqrt(Add(c_, x_)), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule16 = ReplacementRule(pattern16, lambda x, c, a : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule16)

    pattern17 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, Mul(b_, x_))), Integer(-1)), Pow(Sqrt(Add(c_, Mul(d_, x_))), Integer(-1))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (b, x, c, a, d)))
    rule17 = ReplacementRule(pattern17, lambda b, x, c, a, d : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule17)

    pattern18 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, x_)), Integer(-1)), Pow(Sqrt(Add(c_, Mul(d_, x_))), Integer(-1))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule18 = ReplacementRule(pattern18, lambda x, c, a, d : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule18)

    pattern19 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, Mul(b_, x_))), Integer(-1)), Pow(Sqrt(Add(c_, x_)), Integer(-1))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (b, x, c, a)))
    rule19 = ReplacementRule(pattern19, lambda b, x, c, a : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule19)

    pattern20 = Pattern(Int(Mul(Pow(Sqrt(Add(a_, x_)), Integer(-1)), Pow(Sqrt(Add(c_, x_)), Integer(-1))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule20 = ReplacementRule(pattern20, lambda x, c, a : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule20)

    pattern21 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (b, x, c, a, d)))
    rule21 = ReplacementRule(pattern21, lambda b, x, c, a, d : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule21)

    pattern22 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule22 = ReplacementRule(pattern22, lambda x, c, a, d : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule22)

    pattern23 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (b, x, c, a)))
    rule23 = ReplacementRule(pattern23, lambda b, x, c, a : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule23)

    pattern24 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule24 = ReplacementRule(pattern24, lambda x, c, a : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), Integer(-1)), x))
    rubi.add(rule24)

    pattern25 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, x, c, a, d)))
    rule25 = ReplacementRule(pattern25, lambda b, x, c, a, d : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule25)

    pattern26 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (b, x, c, d)))
    rule26 = ReplacementRule(pattern26, lambda b, x, c, d : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule26)

    pattern27 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule27 = ReplacementRule(pattern27, lambda x, c, a, d : Add(Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule27)

    pattern28 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (b, x, a, d)))
    rule28 = ReplacementRule(pattern28, lambda b, x, a, d : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule28)

    pattern29 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, x, c, a)))
    rule29 = ReplacementRule(pattern29, lambda b, x, c, a : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule29)

    pattern30 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule30 = ReplacementRule(pattern30, lambda x, c, d : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule30)

    pattern31 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (b, x, d)))
    rule31 = ReplacementRule(pattern31, lambda b, x, d : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule31)

    pattern32 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (b, x, c)))
    rule32 = ReplacementRule(pattern32, lambda b, x, c : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule32)

    pattern33 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule33 = ReplacementRule(pattern33, lambda x, a, d : Add(Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule33)

    pattern34 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule34 = ReplacementRule(pattern34, lambda x, c, a : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule34)

    pattern35 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (b, x, a)))
    rule35 = ReplacementRule(pattern35, lambda b, x, a : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule35)

    pattern36 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule36 = ReplacementRule(pattern36, lambda x, d : Add(Mul(zoo, d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule36)

    pattern37 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule37 = ReplacementRule(pattern37, lambda x, c : Add(Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule37)

    pattern38 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (b, x)))
    rule38 = ReplacementRule(pattern38, lambda b, x : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule38)

    pattern39 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule39 = ReplacementRule(pattern39, lambda x, a : Add(Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule39)

    pattern40 = Pattern(Int(Pow(x_, Integer(-2)), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule40 = ReplacementRule(pattern40, lambda x : Mul(Integer(2), zoo, Int(Pow(x, Integer(-1)), x)))
    rubi.add(rule40)

    pattern41 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), x_), cons(FreeQ(List(a_, b_), x_), (a, b, x)))
    rule41 = ReplacementRule(pattern41, lambda a, b, x : Mul(Pow(b, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))))
    rubi.add(rule41)

    pattern42 = Pattern(Int(Pow(Add(a_, x_), Integer(-1)), x_), cons(FreeQ(List(a_, Integer(1)), x_), (a, x)))
    rule42 = ReplacementRule(pattern42, lambda a, x : Log(RemoveContent(Add(a, x), x)))
    rubi.add(rule42)

    pattern43 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(True, (x,)))
    rule43 = ReplacementRule(pattern43, lambda x : Log(x))
    rubi.add(rule43)

    pattern44 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, b, x, c, a, d)))
    rule44 = ReplacementRule(pattern44, lambda m, b, x, c, a, d : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule44)

    pattern45 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (m, x, c, a, d)))
    rule45 = ReplacementRule(pattern45, lambda m, x, c, a, d : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule45)

    pattern46 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (m, b, x, c, a)))
    rule46 = ReplacementRule(pattern46, lambda m, b, x, c, a : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule46)

    pattern47 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (m, x, c, a)))
    rule47 = ReplacementRule(pattern47, lambda m, x, c, a : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule47)

    pattern48 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(m_, n_, Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, b, n, x, c, a, d)))
    rule48 = ReplacementRule(pattern48, lambda m, b, n, x, c, a, d : Add(Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule48)

    pattern49 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(m_, n_, Integer(0)), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (m, n, x, c, a, d)))
    rule49 = ReplacementRule(pattern49, lambda m, n, x, c, a, d : Add(Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule49)

    pattern50 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (m, b, n, x, c, a)))
    rule50 = ReplacementRule(pattern50, lambda m, b, n, x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule50)

    pattern51 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (m, n, x, c, a)))
    rule51 = ReplacementRule(pattern51, lambda m, n, x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule51)

    pattern52 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(Integer(0), m_, n_), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, b, n, x, c, a, d)))
    rule52 = ReplacementRule(pattern52, lambda m, b, n, x, c, a, d : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule52)

    pattern53 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(Integer(0), m_, n_), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (m, n, x, c, a, d)))
    rule53 = ReplacementRule(pattern53, lambda m, n, x, c, a, d : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule53)

    pattern54 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (m, b, n, x, c, a)))
    rule54 = ReplacementRule(pattern54, lambda m, b, n, x, c, a : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule54)

    pattern55 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Rational(Integer(1), Integer(2)))), IntegerQ(Add(n_, Rational(Integer(1), Integer(2)))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (m, n, x, c, a)))
    rule55 = ReplacementRule(pattern55, lambda m, n, x, c, a : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule55)

    pattern56 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), FreeQ(List(a_, b_, c_, d_, m_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, b, x, c, a, d)))
    rule56 = ReplacementRule(pattern56, lambda m, b, x, c, a, d : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule56)

    pattern57 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_), x_)), (m, x, c, a, d)))
    rule57 = ReplacementRule(pattern57, lambda m, x, c, a, d : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule57)

    pattern58 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_), x_)), (m, b, x, c, a)))
    rule58 = ReplacementRule(pattern58, lambda m, b, x, c, a : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule58)

    pattern59 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, c_)), Not(IntegerQ(Mul(Integer(2), m_))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_), x_)), (m, x, c, a)))
    rule59 = ReplacementRule(pattern59, lambda m, x, c, a : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Pow(x, Integer(2))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x)))
    rubi.add(rule59)

    pattern60 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Rational(Integer(3), Integer(2)))), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, b, x, c, a, d)))
    rule60 = ReplacementRule(pattern60, lambda m, b, x, c, a, d : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule60)

    pattern61 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Rational(Integer(3), Integer(2)))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (m, x, c, a, d)))
    rule61 = ReplacementRule(pattern61, lambda m, x, c, a, d : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule61)

    pattern62 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Rational(Integer(3), Integer(2)))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (m, b, x, c, a)))
    rule62 = ReplacementRule(pattern62, lambda m, b, x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule62)

    pattern63 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Rational(Integer(3), Integer(2)))), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (m, x, c, a)))
    rule63 = ReplacementRule(pattern63, lambda m, x, c, a : Add(Mul(Integer(-1), Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Rational(Integer(1), Integer(2)), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule63)

    pattern64 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(FreeQ(List(a_, b_, c_, d_, m_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, b, x, c, a, d)))
    rule64 = ReplacementRule(pattern64, lambda m, b, x, c, a, d : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule64)

    pattern65 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), cons(And(FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (b, x, c, a, d)))
    rule65 = ReplacementRule(pattern65, lambda b, x, c, a, d : Int(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), x))
    rubi.add(rule65)

    pattern66 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, x, c, a, d)))
    rule66 = ReplacementRule(pattern66, lambda m, x, c, a, d : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule66)

    pattern67 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, b, x, c, a)))
    rule67 = ReplacementRule(pattern67, lambda m, b, x, c, a : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x))
    rubi.add(rule67)

    pattern68 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a, d)))
    rule68 = ReplacementRule(pattern68, lambda x, c, a, d : Int(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), x))
    rubi.add(rule68)

    pattern69 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (b, x, c, a)))
    rule69 = ReplacementRule(pattern69, lambda b, x, c, a : Int(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), x))
    rubi.add(rule69)

    pattern70 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, x, c, a)))
    rule70 = ReplacementRule(pattern70, lambda m, x, c, a : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x))
    rubi.add(rule70)

    pattern71 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a)))
    rule71 = ReplacementRule(pattern71, lambda x, c, a : Int(Add(Mul(a, c), Pow(x, Integer(2))), x))
    rubi.add(rule71)

    pattern72 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, n, x, c, a, d)))
    rule72 = ReplacementRule(pattern72, lambda m, b, n, x, c, a, d : Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule72)

    pattern73 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(a_, b_, c_, d_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, n, x, c, a, d)))
    rule73 = ReplacementRule(pattern73, lambda b, n, x, c, a, d : Mul(Rational(Integer(1), Integer(2)), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule73)

    pattern74 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), b_, c_, d_, m_, n_), x_)), (m, b, n, x, c, d)))
    rule74 = ReplacementRule(pattern74, lambda m, b, n, x, c, d : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule74)

    pattern75 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_)), (m, n, x, c, a, d)))
    rule75 = ReplacementRule(pattern75, lambda m, n, x, c, a, d : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule75)

    pattern76 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_)), (m, b, n, x, c, a)))
    rule76 = ReplacementRule(pattern76, lambda m, b, n, x, c, a : Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule76)

    pattern77 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), b_, c_, d_, Integer(1), n_), x_)), (b, n, x, c, d)))
    rule77 = ReplacementRule(pattern77, lambda b, n, x, c, d : Mul(Rational(Integer(1), Integer(2)), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule77)

    pattern78 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1), n_), x_)), (n, x, c, a, d)))
    rule78 = ReplacementRule(pattern78, lambda n, x, c, a, d : Mul(Rational(Integer(1), Integer(2)), Pow(Add(a, x), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule78)

    pattern79 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1), n_), x_)), (b, n, x, c, a)))
    rule79 = ReplacementRule(pattern79, lambda b, n, x, c, a : Mul(Rational(Integer(1), Integer(2)), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule79)

    pattern80 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, n_), x_)), (m, n, x, c, d)))
    rule80 = ReplacementRule(pattern80, lambda m, n, x, c, d : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule80)

    pattern81 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, n_), x_)), (m, b, n, x, c)))
    rule81 = ReplacementRule(pattern81, lambda m, b, n, x, c : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule81)

    pattern82 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_)), (m, n, x, c, a)))
    rule82 = ReplacementRule(pattern82, lambda m, n, x, c, a : Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule82)

    pattern83 = Pattern(Int(Mul(x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1), n_), x_)), (n, x, c, d)))
    rule83 = ReplacementRule(pattern83, lambda n, x, c, d : Mul(Rational(Integer(1), Integer(2)), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule83)

    pattern84 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1), n_), x_)), (b, n, x, c)))
    rule84 = ReplacementRule(pattern84, lambda b, n, x, c : Mul(Rational(Integer(1), Integer(2)), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule84)

    pattern85 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1), n_), x_)), (n, x, c, a)))
    rule85 = ReplacementRule(pattern85, lambda n, x, c, a : Mul(Rational(Integer(1), Integer(2)), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule85)

    pattern86 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, n_), x_)), (m, n, x, c)))
    rule86 = ReplacementRule(pattern86, lambda m, n, x, c : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule86)

    pattern87 = Pattern(Int(Mul(x_, Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1), n_), x_)), (n, x, c)))
    rule87 = ReplacementRule(pattern87, lambda n, x, c : Mul(Rational(Integer(1), Integer(2)), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule87)

    pattern88 = Pattern(Int(Pow(Add(a_, Mul(b_, u_)), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (u, a, b, m)))
    rule88 = ReplacementRule(pattern88, lambda u, a, b, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, Mul(b, x)), m), x), x, u)))
    rubi.add(rule88)

    pattern89 = Pattern(Int(Pow(Mul(b_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (u, b, m)))
    rule89 = ReplacementRule(pattern89, lambda u, b, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Mul(b, x), m), x), x, u)))
    rubi.add(rule89)

    pattern90 = Pattern(Int(Pow(Add(a_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (u, a, m)))
    rule90 = ReplacementRule(pattern90, lambda u, a, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, x), m), x), x, u)))
    rubi.add(rule90)

    pattern91 = Pattern(Int(Pow(u_, m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (u, m)))
    rule91 = ReplacementRule(pattern91, lambda u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, m), x), x, u)))
    rubi.add(rule91)

    pattern92 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, b_, m_), x_)), (a, b, m, x)))
    rule92 = ReplacementRule(pattern92, lambda a, b, m, x : Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule92)

    pattern93 = Pattern(Int(Pow(Mul(b_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), b_, m_), x_)), (b, m, x)))
    rule93 = ReplacementRule(pattern93, lambda b, m, x : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule93)

    pattern94 = Pattern(Int(Pow(Add(a_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, Integer(1), m_), x_)), (a, m, x)))
    rule94 = ReplacementRule(pattern94, lambda a, m, x : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule94)

    pattern95 = Pattern(Int(Pow(x_, m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), Integer(1), m_), x_)), (m, x)))
    rule95 = ReplacementRule(pattern95, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule95)

    pattern96 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(m_, Integer(1)))), (m, x)))
    rule96 = ReplacementRule(pattern96, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule96)

    pattern97 = Pattern(Int(x_, x_), cons(And(NonzeroQ(Integer(2)), FreeQ(Integer(1), x_)), (x,)))
    rule97 = ReplacementRule(pattern97, lambda x : Mul(Rational(Integer(1), Integer(2)), Pow(x, Integer(2))))
    rubi.add(rule97)

    return rubi
