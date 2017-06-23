
from matchpy import Wildcard, Pattern, ReplacementRule, ManyToOneReplacer
from .operation import *
from .symbol import VariableSymbol, Integer
from .constraint import cons

a, b, c, d, e, f, g, h, x, u, p = map(VariableSymbol, 'abcdefghxup')
n, m = map(VariableSymbol, 'nm')
zoo = VariableSymbol('zoo')

a_, b_, c_, d_, e_, f_, g_, h_, _p = map(Wildcard.dot, 'abcdefghp')
n_, m_ = map(Wildcard.dot, 'nm')
x_, u_ = map(Wildcard.symbol, 'xu')


def rubi_object():
    rubi = ManyToOneReplacer()
    pattern1 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(True, (x,)))
    rule1 = ReplacementRule(pattern1, lambda x : Log(x))
    rubi.add(rule1)

    pattern2 = Pattern(Int(Pow(x_, m_), x_), cons(And(FreeQ(m_, x_), NonzeroQ(Add(m_, Integer(1)))), (m, x)))
    rule2 = ReplacementRule(pattern2, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule2)

    pattern3 = Pattern(Int(x_, x_), cons(And(NonzeroQ(Integer(2)), FreeQ(Integer(1), x_)), (x,)))
    rule3 = ReplacementRule(pattern3, lambda x : Mul(Integer(1/2), Pow(x, Integer(2))))
    rubi.add(rule3)

    pattern4 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), x_), cons(FreeQ(List(a_, b_), x_), (b, a, x)))
    rule4 = ReplacementRule(pattern4, lambda b, a, x : Mul(Pow(b, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))))
    rubi.add(rule4)

    pattern5 = Pattern(Int(Pow(Add(a_, x_), Integer(-1)), x_), cons(FreeQ(List(a_, Integer(1)), x_), (a, x)))
    rule5 = ReplacementRule(pattern5, lambda a, x : Log(RemoveContent(Add(a, x), x)))
    rubi.add(rule5)

    pattern6 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, b_, m_), x_)), (m, b, a, x)))
    rule6 = ReplacementRule(pattern6, lambda m, b, a, x : Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule6)

    pattern7 = Pattern(Int(Pow(Mul(b_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), b_, m_), x_)), (m, b, x)))
    rule7 = ReplacementRule(pattern7, lambda m, b, x : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule7)

    pattern8 = Pattern(Int(Pow(Add(a_, x_), m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(a_, Integer(1), m_), x_)), (m, a, x)))
    rule8 = ReplacementRule(pattern8, lambda m, a, x : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule8)

    pattern9 = Pattern(Int(Pow(x_, m_), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), FreeQ(List(Integer(0), Integer(1), m_), x_)), (m, x)))
    rule9 = ReplacementRule(pattern9, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule9)

    pattern10 = Pattern(Int(Pow(Add(a_, Mul(b_, u_)), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (x, u, m, a, b)))
    rule10 = ReplacementRule(pattern10, lambda x, u, m, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, Mul(b, x)), m), x), x, u)))
    rubi.add(rule10)

    pattern11 = Pattern(Int(Pow(Mul(b_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), b_, m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (x, u, m, b)))
    rule11 = ReplacementRule(pattern11, lambda x, u, m, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Mul(b, x), m), x), x, u)))
    rubi.add(rule11)

    pattern12 = Pattern(Int(Pow(Add(a_, u_), m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(a_, Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (x, u, m, a)))
    rule12 = ReplacementRule(pattern12, lambda x, u, m, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, x), m), x), x, u)))
    rubi.add(rule12)

    pattern13 = Pattern(Int(Pow(u_, m_), x_), cons(And(LinearQ(u_, x_), FreeQ(List(Integer(0), Integer(1), m_), x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (x, u, m)))
    rule13 = ReplacementRule(pattern13, lambda x, u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, m), x), x, u)))
    rubi.add(rule13)

    pattern14 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule14 = ReplacementRule(pattern14, lambda x, c, a, d, b : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule14)

    pattern15 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule15 = ReplacementRule(pattern15, lambda x, c, a, d : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule15)

    pattern16 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule16 = ReplacementRule(pattern16, lambda x, c, a, b : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule16)

    pattern17 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule17 = ReplacementRule(pattern17, lambda x, c, a : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), Integer(-1)), x))
    rubi.add(rule17)

    pattern18 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule18 = ReplacementRule(pattern18, lambda x, c, a, d, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule18)

    pattern19 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule19 = ReplacementRule(pattern19, lambda x, c, d, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule19)

    pattern20 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule20 = ReplacementRule(pattern20, lambda x, c, a, d : Add(Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule20)

    pattern21 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule21 = ReplacementRule(pattern21, lambda x, a, d, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule21)

    pattern22 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule22 = ReplacementRule(pattern22, lambda x, c, a, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule22)

    pattern23 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule23 = ReplacementRule(pattern23, lambda x, c, d : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule23)

    pattern24 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule24 = ReplacementRule(pattern24, lambda x, d, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule24)

    pattern25 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule25 = ReplacementRule(pattern25, lambda x, c, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule25)

    pattern26 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule26 = ReplacementRule(pattern26, lambda x, a, d : Add(Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule26)

    pattern27 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule27 = ReplacementRule(pattern27, lambda x, c, a : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule27)

    pattern28 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule28 = ReplacementRule(pattern28, lambda x, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule28)

    pattern29 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule29 = ReplacementRule(pattern29, lambda x, d : Add(Mul(zoo, d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule29)

    pattern30 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule30 = ReplacementRule(pattern30, lambda x, c : Add(Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule30)

    pattern31 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule31 = ReplacementRule(pattern31, lambda x, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule31)

    pattern32 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule32 = ReplacementRule(pattern32, lambda x, a : Add(Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule32)

    pattern33 = Pattern(Int(Pow(x_, Integer(-2)), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule33 = ReplacementRule(pattern33, lambda x : Mul(Integer(2), zoo, Int(Pow(x, Integer(-1)), x)))
    rubi.add(rule33)

    pattern34 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule34 = ReplacementRule(pattern34, lambda x, m, n, c, a, d, b : Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule34)

    pattern35 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(a_, b_, c_, d_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, n, c, a, d, b)))
    rule35 = ReplacementRule(pattern35, lambda x, n, c, a, d, b : Mul(Integer(1/2), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule35)

    pattern36 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), b_, c_, d_, m_, n_), x_)), (x, m, n, c, d, b)))
    rule36 = ReplacementRule(pattern36, lambda x, m, n, c, d, b : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule36)

    pattern37 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_)), (x, m, n, c, a, d)))
    rule37 = ReplacementRule(pattern37, lambda x, m, n, c, a, d : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule37)

    pattern38 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_)), (x, m, n, c, a, b)))
    rule38 = ReplacementRule(pattern38, lambda x, m, n, c, a, b : Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule38)

    pattern39 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), b_, c_, d_, Integer(1), n_), x_)), (x, n, c, d, b)))
    rule39 = ReplacementRule(pattern39, lambda x, n, c, d, b : Mul(Integer(1/2), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule39)

    pattern40 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1), n_), x_)), (x, n, c, a, d)))
    rule40 = ReplacementRule(pattern40, lambda x, n, c, a, d : Mul(Integer(1/2), Pow(Add(a, x), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule40)

    pattern41 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1), n_), x_)), (x, n, c, a, b)))
    rule41 = ReplacementRule(pattern41, lambda x, n, c, a, b : Mul(Integer(1/2), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule41)

    pattern42 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, n_), x_)), (x, m, n, c, d)))
    rule42 = ReplacementRule(pattern42, lambda x, m, n, c, d : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule42)

    pattern43 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, n_), x_)), (x, m, n, c, b)))
    rule43 = ReplacementRule(pattern43, lambda x, m, n, c, b : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule43)

    pattern44 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_)), (x, m, n, c, a)))
    rule44 = ReplacementRule(pattern44, lambda x, m, n, c, a : Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule44)

    pattern45 = Pattern(Int(Mul(x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1), n_), x_)), (x, n, c, d)))
    rule45 = ReplacementRule(pattern45, lambda x, n, c, d : Mul(Integer(1/2), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule45)

    pattern46 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1), n_), x_)), (x, n, c, b)))
    rule46 = ReplacementRule(pattern46, lambda x, n, c, b : Mul(Integer(1/2), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule46)

    pattern47 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1), n_), x_)), (x, n, c, a)))
    rule47 = ReplacementRule(pattern47, lambda x, n, c, a : Mul(Integer(1/2), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule47)

    pattern48 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, n_), x_)), (x, m, n, c)))
    rule48 = ReplacementRule(pattern48, lambda x, m, n, c : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule48)

    pattern49 = Pattern(Int(Mul(x_, Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1), n_), x_)), (x, n, c)))
    rule49 = ReplacementRule(pattern49, lambda x, n, c : Mul(Integer(1/2), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule49)

    pattern50 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, m, c, a, d, b)))
    rule50 = ReplacementRule(pattern50, lambda x, m, c, a, d, b : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule50)

    pattern51 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, m, c, a, d)))
    rule51 = ReplacementRule(pattern51, lambda x, m, c, a, d : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule51)

    pattern52 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, m, c, a, b)))
    rule52 = ReplacementRule(pattern52, lambda x, m, c, a, b : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule52)

    pattern53 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, c, a)))
    rule53 = ReplacementRule(pattern53, lambda x, m, c, a : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule53)

    pattern54 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-3/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-3/2))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule54 = ReplacementRule(pattern54, lambda x, c, a, d, b : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule54)

    pattern55 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-3/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-3/2))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule55 = ReplacementRule(pattern55, lambda x, c, a, d : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule55)

    pattern56 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-3/2)), Pow(Add(c_, x_), Integer(-3/2))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule56 = ReplacementRule(pattern56, lambda x, c, a, b : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule56)

    pattern57 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-3/2)), Pow(Add(c_, x_), Integer(-3/2))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule57 = ReplacementRule(pattern57, lambda x, c, a : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule57)

    pattern58 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, m, c, a, d, b)))
    rule58 = ReplacementRule(pattern58, lambda x, m, c, a, d, b : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule58)

    pattern59 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, m, c, a, d)))
    rule59 = ReplacementRule(pattern59, lambda x, m, c, a, d : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule59)

    pattern60 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, m, c, a, b)))
    rule60 = ReplacementRule(pattern60, lambda x, m, c, a, b : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule60)

    pattern61 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, c, a)))
    rule61 = ReplacementRule(pattern61, lambda x, m, c, a : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule61)

    pattern62 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(FreeQ(List(a_, b_, c_, d_, m_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (x, m, c, a, d, b)))
    rule62 = ReplacementRule(pattern62, lambda x, m, c, a, d, b : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule62)

    pattern63 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), cons(And(FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a, d, b)))
    rule63 = ReplacementRule(pattern63, lambda x, c, a, d, b : Int(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), x))
    rubi.add(rule63)

    pattern64 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (x, m, c, a, d)))
    rule64 = ReplacementRule(pattern64, lambda x, m, c, a, d : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule64)

    pattern65 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (x, m, c, a, b)))
    rule65 = ReplacementRule(pattern65, lambda x, m, c, a, b : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x))
    rubi.add(rule65)

    pattern66 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a, d)))
    rule66 = ReplacementRule(pattern66, lambda x, c, a, d : Int(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), x))
    rubi.add(rule66)

    pattern67 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a, b)))
    rule67 = ReplacementRule(pattern67, lambda x, c, a, b : Int(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), x))
    rubi.add(rule67)

    pattern68 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_), x_), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (x, m, c, a)))
    rule68 = ReplacementRule(pattern68, lambda x, m, c, a : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x))
    rubi.add(rule68)

    pattern69 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (x, c, a)))
    rule69 = ReplacementRule(pattern69, lambda x, c, a : Int(Add(Mul(a, c), Pow(x, Integer(2))), x))
    rubi.add(rule69)

    pattern70 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule70 = ReplacementRule(pattern70, lambda x, c, a, d, b : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule70)

    pattern71 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule71 = ReplacementRule(pattern71, lambda x, c, a, d : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule71)

    pattern72 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule72 = ReplacementRule(pattern72, lambda x, c, a, b : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule72)

    pattern73 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule73 = ReplacementRule(pattern73, lambda x, c, a : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule73)

    pattern74 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule74 = ReplacementRule(pattern74, lambda x, c, a, d, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule74)

    pattern75 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule75 = ReplacementRule(pattern75, lambda x, c, a, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule75)

    pattern76 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule76 = ReplacementRule(pattern76, lambda x, c, a, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule76)

    pattern77 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule77 = ReplacementRule(pattern77, lambda x, c, a : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule77)

    pattern78 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), FreeQ(List(a_, b_, c_, d_, m_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, m, c, a, d, b)))
    rule78 = ReplacementRule(pattern78, lambda x, m, c, a, d, b : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule78)

    pattern79 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_), x_)), (x, m, c, a, d)))
    rule79 = ReplacementRule(pattern79, lambda x, m, c, a, d : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule79)

    pattern80 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_), x_)), (x, m, c, a, b)))
    rule80 = ReplacementRule(pattern80, lambda x, m, c, a, b : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule80)

    pattern81 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(ZeroQ(Add(a_, c_)), Not(IntegerQ(Mul(Integer(2), m_))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_), x_)), (x, m, c, a)))
    rule81 = ReplacementRule(pattern81, lambda x, m, c, a : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Pow(x, Integer(2))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x)))
    rubi.add(rule81)

    pattern82 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-5/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (x, c, a, d, b)))
    rule82 = ReplacementRule(pattern82, lambda x, c, a, d, b : Add(Mul(Integer(1/2), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-1/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule82)

    pattern83 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-5/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (x, c, a, d)))
    rule83 = ReplacementRule(pattern83, lambda x, c, a, d : Add(Mul(Integer(1/2), Add(Mul(Integer(-1), a, d), c), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Integer(-1/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule83)

    pattern84 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-5/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (x, c, a, b)))
    rule84 = ReplacementRule(pattern84, lambda x, c, a, b : Add(Mul(Integer(1/2), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-1/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule84)

    pattern85 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-5/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (x, c, a)))
    rule85 = ReplacementRule(pattern85, lambda x, c, a : Add(Mul(Integer(1/2), Add(Mul(Integer(-1), a), c), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Integer(-1/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule85)

    pattern86 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-9/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (x, c, a, d, b)))
    rule86 = ReplacementRule(pattern86, lambda x, c, a, d, b : Add(Mul(Integer(-1), Integer(1/5), Pow(b, Integer(-1)), d, Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule86)

    pattern87 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-9/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (x, c, a, d)))
    rule87 = ReplacementRule(pattern87, lambda x, c, a, d : Add(Mul(Integer(-1), Integer(1/5), d, Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule87)

    pattern88 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-9/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (x, c, a, b)))
    rule88 = ReplacementRule(pattern88, lambda x, c, a, b : Add(Mul(Integer(-1), Integer(1/5), Pow(b, Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule88)

    pattern89 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-9/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), cons(And(ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (x, c, a)))
    rule89 = ReplacementRule(pattern89, lambda x, c, a : Add(Mul(Integer(-1), Integer(1/5), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule89)

    pattern90 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule90 = ReplacementRule(pattern90, lambda x, m, n, c, a, d, b : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule90)

    pattern91 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, m, n, c, a, d)))
    rule91 = ReplacementRule(pattern91, lambda x, m, n, c, a, d : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule91)

    pattern92 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, m, n, c, a, b)))
    rule92 = ReplacementRule(pattern92, lambda x, m, n, c, a, b : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule92)

    pattern93 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, n, c, a)))
    rule93 = ReplacementRule(pattern93, lambda x, m, n, c, a : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule93)

    pattern94 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule94 = ReplacementRule(pattern94, lambda x, m, n, c, a, d, b : Add(Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule94)

    pattern95 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(Mul(a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, m, n, c, a, d)))
    rule95 = ReplacementRule(pattern95, lambda x, m, n, c, a, d : Add(Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule95)

    pattern96 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, m, n, c, a, b)))
    rule96 = ReplacementRule(pattern96, lambda x, m, n, c, a, b : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule96)

    pattern97 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, n, c, a)))
    rule97 = ReplacementRule(pattern97, lambda x, m, n, c, a : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule97)

    pattern98 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(m_), FreeQ(List(a_, b_, c_, d_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, a, d, b)))
    rule98 = ReplacementRule(pattern98, lambda x, m, n, c, a, d, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule98)

    pattern99 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, d_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, a, d, b)))
    rule99 = ReplacementRule(pattern99, lambda x, n, c, a, d, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule99)

    pattern100 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, d, b)))
    rule100 = ReplacementRule(pattern100, lambda x, m, n, c, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule100)

    pattern101 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, a, d)))
    rule101 = ReplacementRule(pattern101, lambda x, m, n, c, a, d : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule101)

    pattern102 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(m_), FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, a, d, b)))
    rule102 = ReplacementRule(pattern102, lambda x, m, c, a, d, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule102)

    pattern103 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, a, d, b)))
    rule103 = ReplacementRule(pattern103, lambda x, m, n, a, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule103)

    pattern104 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(m_), FreeQ(List(a_, b_, c_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, a, b)))
    rule104 = ReplacementRule(pattern104, lambda x, m, n, c, a, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule104)

    pattern105 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, d, b)))
    rule105 = ReplacementRule(pattern105, lambda x, n, c, d, b : Int(ExpandIntegrand(Mul(b, x, Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule105)

    pattern106 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, a, d)))
    rule106 = ReplacementRule(pattern106, lambda x, n, c, a, d : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule106)

    pattern107 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, a, d, b)))
    rule107 = ReplacementRule(pattern107, lambda x, c, a, d, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x))
    rubi.add(rule107)

    pattern108 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, Mul(b_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, a, d, b)))
    rule108 = ReplacementRule(pattern108, lambda x, n, a, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule108)

    pattern109 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, a, b)))
    rule109 = ReplacementRule(pattern109, lambda x, n, c, a, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x))
    rubi.add(rule109)

    pattern110 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, d)))
    rule110 = ReplacementRule(pattern110, lambda x, m, n, c, d : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule110)

    pattern111 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, d, b)))
    rule111 = ReplacementRule(pattern111, lambda x, m, c, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule111)

    pattern112 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), b_, Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, d, b)))
    rule112 = ReplacementRule(pattern112, lambda x, m, n, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Mul(d, x), n)), x), x))
    rubi.add(rule112)

    pattern113 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, b)))
    rule113 = ReplacementRule(pattern113, lambda x, m, n, c, b : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule113)

    pattern114 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, a, d)))
    rule114 = ReplacementRule(pattern114, lambda x, m, c, a, d : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule114)

    pattern115 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, a, d)))
    rule115 = ReplacementRule(pattern115, lambda x, m, n, a, d : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule115)

    pattern116 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c, a)))
    rule116 = ReplacementRule(pattern116, lambda x, m, n, c, a : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule116)

    pattern117 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, a, d, b)))
    rule117 = ReplacementRule(pattern117, lambda x, m, a, d, b : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule117)

    pattern118 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, x_)), x_), cons(And(PositiveIntegerQ(m_), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, a, b)))
    rule118 = ReplacementRule(pattern118, lambda x, m, c, a, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x))
    rubi.add(rule118)

    pattern119 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, a, b)))
    rule119 = ReplacementRule(pattern119, lambda x, m, n, a, b : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule119)

    pattern120 = Pattern(Int(Mul(x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), c_, d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, d)))
    rule120 = ReplacementRule(pattern120, lambda x, n, c, d : Int(ExpandIntegrand(Mul(x, Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule120)

    pattern121 = Pattern(Int(Mul(b_, x_, Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, d, b)))
    rule121 = ReplacementRule(pattern121, lambda x, c, d, b : Int(ExpandIntegrand(Mul(b, x, Add(c, Mul(d, x))), x), x))
    rubi.add(rule121)

    pattern122 = Pattern(Int(Mul(b_, x_, Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), b_, Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, d, b)))
    rule122 = ReplacementRule(pattern122, lambda x, n, d, b : Int(ExpandIntegrand(Mul(b, x, Pow(Mul(d, x), n)), x), x))
    rubi.add(rule122)

    pattern123 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, b)))
    rule123 = ReplacementRule(pattern123, lambda x, n, c, b : Int(ExpandIntegrand(Mul(b, x, Pow(Add(c, x), n)), x), x))
    rubi.add(rule123)

    pattern124 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, a, d)))
    rule124 = ReplacementRule(pattern124, lambda x, c, a, d : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, Mul(d, x))), x), x))
    rubi.add(rule124)

    pattern125 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, a, d)))
    rule125 = ReplacementRule(pattern125, lambda x, n, a, d : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x))
    rubi.add(rule125)

    pattern126 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c, a)))
    rule126 = ReplacementRule(pattern126, lambda x, n, c, a : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, x), n)), x), x))
    rubi.add(rule126)

    pattern127 = Pattern(Int(Mul(d_, x_, Add(a_, Mul(b_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, a, d, b)))
    rule127 = ReplacementRule(pattern127, lambda x, a, d, b : Int(ExpandIntegrand(Mul(d, x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule127)

    pattern128 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, a, b)))
    rule128 = ReplacementRule(pattern128, lambda x, c, a, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x))
    rubi.add(rule128)

    pattern129 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, Mul(b_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, a, b)))
    rule129 = ReplacementRule(pattern129, lambda x, n, a, b : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Int(Mul(Pow(x_, m_), Add(c_, Mul(d_, x_))), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, d)))
    rule130 = ReplacementRule(pattern130, lambda x, m, c, d : Int(ExpandIntegrand(Mul(Pow(x, m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, d)))
    rule131 = ReplacementRule(pattern131, lambda x, m, n, d : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Mul(d, x), n)), x), x))
    rubi.add(rule131)

    pattern132 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, c)))
    rule132 = ReplacementRule(pattern132, lambda x, m, n, c : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule132)

    pattern133 = Pattern(Int(Mul(d_, x_, Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), b_, Integer(0), d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, d, b)))
    rule133 = ReplacementRule(pattern133, lambda x, m, d, b : Int(ExpandIntegrand(Mul(d, x, Pow(Mul(b, x), m)), x), x))
    rubi.add(rule133)

    pattern134 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Add(c_, x_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, b)))
    rule134 = ReplacementRule(pattern134, lambda x, m, c, b : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Add(c, x)), x), x))
    rubi.add(rule134)

    pattern135 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, b)))
    rule135 = ReplacementRule(pattern135, lambda x, m, n, b : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Mul(b, x), m)), x), x))
    rubi.add(rule135)

    pattern136 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, x_), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, a, d)))
    rule136 = ReplacementRule(pattern136, lambda x, m, a, d : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule136)

    pattern137 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, x_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c, a)))
    rule137 = ReplacementRule(pattern137, lambda x, m, c, a : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, x)), x), x))
    rubi.add(rule137)

    pattern138 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n, a)))
    rule138 = ReplacementRule(pattern138, lambda x, m, n, a : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule138)

    pattern139 = Pattern(Int(Mul(x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, a, b)))
    rule139 = ReplacementRule(pattern139, lambda x, m, a, b : Int(ExpandIntegrand(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule139)

    pattern140 = Pattern(Int(Mul(x_, Add(c_, Mul(d_, x_))), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, d)))
    rule140 = ReplacementRule(pattern140, lambda x, c, d : Int(ExpandIntegrand(Mul(x, Add(c, Mul(d, x))), x), x))
    rubi.add(rule140)

    pattern141 = Pattern(Int(Mul(x_, Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, d)))
    rule141 = ReplacementRule(pattern141, lambda x, n, d : Int(ExpandIntegrand(Mul(x, Pow(Mul(d, x), n)), x), x))
    rubi.add(rule141)

    pattern142 = Pattern(Int(Mul(x_, Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, c)))
    rule142 = ReplacementRule(pattern142, lambda x, n, c : Int(ExpandIntegrand(Mul(x, Pow(Add(c, x), n)), x), x))
    rubi.add(rule142)

    pattern143 = Pattern(Int(Mul(b_, d_, Pow(x_, Integer(2))), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), b_, Integer(0), d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, d, b)))
    rule143 = ReplacementRule(pattern143, lambda x, d, b : Int(ExpandIntegrand(Mul(b, d, Pow(x, Integer(2))), x), x))
    rubi.add(rule143)

    pattern144 = Pattern(Int(Mul(b_, x_, Add(c_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, b)))
    rule144 = ReplacementRule(pattern144, lambda x, c, b : Int(ExpandIntegrand(Mul(b, x, Add(c, x)), x), x))
    rubi.add(rule144)

    pattern145 = Pattern(Int(Mul(b_, x_, Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, b)))
    rule145 = ReplacementRule(pattern145, lambda x, n, b : Int(ExpandIntegrand(Mul(b, x, Pow(x, n)), x), x))
    rubi.add(rule145)

    pattern146 = Pattern(Int(Mul(d_, x_, Add(a_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, a, d)))
    rule146 = ReplacementRule(pattern146, lambda x, a, d : Int(ExpandIntegrand(Mul(d, x, Add(a, x)), x), x))
    rubi.add(rule146)

    pattern147 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c, a)))
    rule147 = ReplacementRule(pattern147, lambda x, c, a : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, x)), x), x))
    rubi.add(rule147)

    pattern148 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n, a)))
    rule148 = ReplacementRule(pattern148, lambda x, n, a : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, x)), x), x))
    rubi.add(rule148)

    pattern149 = Pattern(Int(Mul(x_, Add(a_, Mul(b_, x_))), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, a, b)))
    rule149 = ReplacementRule(pattern149, lambda x, a, b : Int(ExpandIntegrand(Mul(x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule149)

    pattern150 = Pattern(Int(Mul(d_, x_, Pow(x_, m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, d)))
    rule150 = ReplacementRule(pattern150, lambda x, m, d : Int(ExpandIntegrand(Mul(d, x, Pow(x, m)), x), x))
    rubi.add(rule150)

    pattern151 = Pattern(Int(Mul(Pow(x_, m_), Add(c_, x_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, c)))
    rule151 = ReplacementRule(pattern151, lambda x, m, c : Int(ExpandIntegrand(Mul(Pow(x, m), Add(c, x)), x), x))
    rubi.add(rule151)

    pattern152 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (x, m, n)))
    rule152 = ReplacementRule(pattern152, lambda x, m, n : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(x, n)), x), x))
    rubi.add(rule152)

    pattern153 = Pattern(Int(Mul(x_, Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, b)))
    rule153 = ReplacementRule(pattern153, lambda x, m, b : Int(ExpandIntegrand(Mul(x, Pow(Mul(b, x), m)), x), x))
    rubi.add(rule153)

    pattern154 = Pattern(Int(Mul(x_, Pow(Add(a_, x_), m_)), x_), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m, a)))
    rule154 = ReplacementRule(pattern154, lambda x, m, a : Int(ExpandIntegrand(Mul(x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule154)

    pattern155 = Pattern(Int(Mul(d_, Pow(x_, Integer(2))), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, d)))
    rule155 = ReplacementRule(pattern155, lambda x, d : Int(ExpandIntegrand(Mul(d, Pow(x, Integer(2))), x), x))
    rubi.add(rule155)

    pattern156 = Pattern(Int(Mul(x_, Add(c_, x_)), x_), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (x, c)))
    rule156 = ReplacementRule(pattern156, lambda x, c : Int(ExpandIntegrand(Mul(x, Add(c, x)), x), x))
    rubi.add(rule156)

    pattern157 = Pattern(Int(Mul(x_, Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), n_), x_), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n)))
    rule157 = ReplacementRule(pattern157, lambda x, n : Int(ExpandIntegrand(Mul(x, Pow(x, n)), x), x))
    rubi.add(rule157)

    pattern158 = Pattern(Int(Mul(b_, Pow(x_, Integer(2))), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, b)))
    rule158 = ReplacementRule(pattern158, lambda x, b : Int(ExpandIntegrand(Mul(b, Pow(x, Integer(2))), x), x))
    rubi.add(rule158)

    pattern159 = Pattern(Int(Mul(x_, Add(a_, x_)), x_), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x, a)))
    rule159 = ReplacementRule(pattern159, lambda x, a : Int(ExpandIntegrand(Mul(x, Add(a, x)), x), x))
    rubi.add(rule159)

    pattern160 = Pattern(Int(Mul(x_, Pow(x_, m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), Integer(1)), x_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (x, m)))
    rule160 = ReplacementRule(pattern160, lambda x, m : Int(ExpandIntegrand(Mul(x, Pow(x, m)), x), x))
    rubi.add(rule160)

    pattern161 = Pattern(Int(Pow(x_, Integer(2)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), Integer(1)), x_), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x,)))
    rule161 = ReplacementRule(pattern161, lambda x : Int(ExpandIntegrand(Pow(x, Integer(2)), x), x))
    rubi.add(rule161)

    pattern162 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), FreeQ(List(a_, b_, c_, d_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, c, a, d, b)))
    rule162 = ReplacementRule(pattern162, lambda x, n, m, c, a, d, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule162)

    pattern163 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, d_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, c, a, d, b)))
    rule163 = ReplacementRule(pattern163, lambda x, n, c, a, d, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule163)

    pattern164 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, c, a, d)))
    rule164 = ReplacementRule(pattern164, lambda x, n, m, c, a, d : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule164)

    pattern165 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, Mul(d_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, c, a, d, b)))
    rule165 = ReplacementRule(pattern165, lambda x, m, c, a, d, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule165)

    pattern166 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, a, d, b)))
    rule166 = ReplacementRule(pattern166, lambda x, n, m, a, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule166)

    pattern167 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), FreeQ(List(a_, b_, c_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, c, a, b)))
    rule167 = ReplacementRule(pattern167, lambda x, n, m, c, a, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule167)

    pattern168 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, c, a, d)))
    rule168 = ReplacementRule(pattern168, lambda x, n, c, a, d : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule168)

    pattern169 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, b_, c_, d_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule169 = ReplacementRule(pattern169, lambda x, c, a, d, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x))
    rubi.add(rule169)

    pattern170 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, Mul(b_, x_))), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, a, d, b)))
    rule170 = ReplacementRule(pattern170, lambda x, n, a, d, b : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule170)

    pattern171 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), FreeQ(List(a_, b_, c_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, c, a, b)))
    rule171 = ReplacementRule(pattern171, lambda x, n, c, a, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x))
    rubi.add(rule171)

    pattern172 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, Mul(d_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, c, a, d)))
    rule172 = ReplacementRule(pattern172, lambda x, m, c, a, d : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule172)

    pattern173 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, a, d)))
    rule173 = ReplacementRule(pattern173, lambda x, n, m, a, d : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule173)

    pattern174 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, c, a)))
    rule174 = ReplacementRule(pattern174, lambda x, n, m, c, a : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule174)

    pattern175 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, a, d, b)))
    rule175 = ReplacementRule(pattern175, lambda x, m, a, d, b : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule175)

    pattern176 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, c, a, b)))
    rule176 = ReplacementRule(pattern176, lambda x, m, c, a, b : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x))
    rubi.add(rule176)

    pattern177 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, a, b)))
    rule177 = ReplacementRule(pattern177, lambda x, n, m, a, b : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule177)

    pattern178 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, Integer(1), c_, d_, Integer(1)), x_)), (x, c, a, d)))
    rule178 = ReplacementRule(pattern178, lambda x, c, a, d : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, Mul(d, x))), x), x))
    rubi.add(rule178)

    pattern179 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, x_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, a, d)))
    rule179 = ReplacementRule(pattern179, lambda x, n, a, d : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x))
    rubi.add(rule179)

    pattern180 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, c, a)))
    rule180 = ReplacementRule(pattern180, lambda x, n, c, a : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, x), n)), x), x))
    rubi.add(rule180)

    pattern181 = Pattern(Int(Mul(d_, x_, Add(a_, Mul(b_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, b_, Integer(0), d_, Integer(1)), x_)), (x, a, d, b)))
    rule181 = ReplacementRule(pattern181, lambda x, a, d, b : Int(ExpandIntegrand(Mul(d, x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule181)

    pattern182 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule182 = ReplacementRule(pattern182, lambda x, c, a, b : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x))
    rubi.add(rule182)

    pattern183 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, Mul(b_, x_))), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, a, b)))
    rule183 = ReplacementRule(pattern183, lambda x, n, a, b : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule183)

    pattern184 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, x_), m_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, a, d)))
    rule184 = ReplacementRule(pattern184, lambda x, m, a, d : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule184)

    pattern185 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, c, a)))
    rule185 = ReplacementRule(pattern185, lambda x, m, c, a : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, x)), x), x))
    rubi.add(rule185)

    pattern186 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (x, n, m, a)))
    rule186 = ReplacementRule(pattern186, lambda x, n, m, a : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule186)

    pattern187 = Pattern(Int(Mul(x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, a, b)))
    rule187 = ReplacementRule(pattern187, lambda x, m, a, b : Int(ExpandIntegrand(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule187)

    pattern188 = Pattern(Int(Mul(d_, x_, Add(a_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1)), x_)), (x, a, d)))
    rule188 = ReplacementRule(pattern188, lambda x, a, d : Int(ExpandIntegrand(Mul(d, x, Add(a, x)), x), x))
    rubi.add(rule188)

    pattern189 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1)), x_)), (x, c, a)))
    rule189 = ReplacementRule(pattern189, lambda x, c, a : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, x)), x), x))
    rubi.add(rule189)

    pattern190 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, x_)), x_), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), n_), x_), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (x, n, a)))
    rule190 = ReplacementRule(pattern190, lambda x, n, a : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, x)), x), x))
    rubi.add(rule190)

    pattern191 = Pattern(Int(Mul(x_, Add(a_, Mul(b_, x_))), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1)), x_)), (x, a, b)))
    rule191 = ReplacementRule(pattern191, lambda x, a, b : Int(ExpandIntegrand(Mul(x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule191)

    pattern192 = Pattern(Int(Mul(x_, Pow(Add(a_, x_), m_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1)), x_), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (x, m, a)))
    rule192 = ReplacementRule(pattern192, lambda x, m, a : Int(ExpandIntegrand(Mul(x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule192)

    pattern193 = Pattern(Int(Mul(x_, Add(a_, x_)), x_), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1)), x_)), (x, a)))
    rule193 = ReplacementRule(pattern193, lambda x, a : Int(ExpandIntegrand(Mul(x, Add(a, x)), x), x))
    rubi.add(rule193)

    pattern194 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, n, c, a, d, b)))
    rule194 = ReplacementRule(pattern194, lambda x, n, c, a, d, b : Add(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule194)

    pattern195 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, n, c, d, b)))
    rule195 = ReplacementRule(pattern195, lambda x, n, c, d, b : Add(Mul(c, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule195)

    pattern196 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, n, c, a, d)))
    rule196 = ReplacementRule(pattern196, lambda x, n, c, a, d : Add(Mul(Add(Mul(Integer(-1), a, d), c), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule196)

    pattern197 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, n, a, d, b)))
    rule197 = ReplacementRule(pattern197, lambda x, n, a, d, b : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d, Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Mul(d, x), n))))
    rubi.add(rule197)

    pattern198 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, c, a, b)))
    rule198 = ReplacementRule(pattern198, lambda x, n, c, a, b : Add(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule198)

    pattern199 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, n, c, d)))
    rule199 = ReplacementRule(pattern199, lambda x, n, c, d : Add(Mul(c, Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule199)

    pattern200 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, n, d, b)))
    rule200 = ReplacementRule(pattern200, lambda x, n, d, b : Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Mul(d, x), n)))
    rubi.add(rule200)

    pattern201 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, n, c, b)))
    rule201 = ReplacementRule(pattern201, lambda x, n, c, b : Add(Mul(c, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule201)

    pattern202 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, n, a, d)))
    rule202 = ReplacementRule(pattern202, lambda x, n, a, d : Add(Mul(Integer(-1), a, d, Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), Integer(-1))), x)), Mul(Pow(n, Integer(-1)), Pow(Mul(d, x), n))))
    rubi.add(rule202)

    pattern203 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, n, c, a)))
    rule203 = ReplacementRule(pattern203, lambda x, n, c, a : Add(Mul(Add(Mul(Integer(-1), a), c), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule203)

    pattern204 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, n, a, b)))
    rule204 = ReplacementRule(pattern204, lambda x, n, a, b : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(x, n))))
    rubi.add(rule204)

    pattern205 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, n, d)))
    rule205 = ReplacementRule(pattern205, lambda x, n, d : Mul(Pow(n, Integer(-1)), Pow(Mul(d, x), n)))
    rubi.add(rule205)

    pattern206 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, n, c)))
    rule206 = ReplacementRule(pattern206, lambda x, n, c : Add(Mul(c, Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule206)

    pattern207 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, n, b)))
    rule207 = ReplacementRule(pattern207, lambda x, n, b : Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(x, n)))
    rubi.add(rule207)

    pattern208 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, n, a)))
    rule208 = ReplacementRule(pattern208, lambda x, n, a : Add(Mul(Integer(-1), a, Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), Integer(-1))), x)), Mul(Pow(n, Integer(-1)), Pow(x, n))))
    rubi.add(rule208)

    pattern209 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x, n)))
    rule209 = ReplacementRule(pattern209, lambda x, n : Mul(Pow(n, Integer(-1)), Pow(x, n)))
    rubi.add(rule209)

    pattern210 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, n, c, a, d, b)))
    rule210 = ReplacementRule(pattern210, lambda x, n, c, a, d, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule210)

    pattern211 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, n, c, d, b)))
    rule211 = ReplacementRule(pattern211, lambda x, n, c, d, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule211)

    pattern212 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, n, c, a, d)))
    rule212 = ReplacementRule(pattern212, lambda x, n, c, a, d : Add(Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x))))
    rubi.add(rule212)

    pattern213 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, n, a, d, b)))
    rule213 = ReplacementRule(pattern213, lambda x, n, a, d, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule213)

    pattern214 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, c, a, b)))
    rule214 = ReplacementRule(pattern214, lambda x, n, c, a, b : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule214)

    pattern215 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, n, c, d)))
    rule215 = ReplacementRule(pattern215, lambda x, n, c, d : Add(Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(c, Integer(-1)), Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x))))
    rubi.add(rule215)

    pattern216 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, n, d, b)))
    rule216 = ReplacementRule(pattern216, lambda x, n, d, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1)))), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule216)

    pattern217 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, n, c, b)))
    rule217 = ReplacementRule(pattern217, lambda x, n, c, b : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule217)

    pattern218 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, n, a, d)))
    rule218 = ReplacementRule(pattern218, lambda x, n, a, d : Add(Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Integer(-1))), x))))
    rubi.add(rule218)

    pattern219 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, n, c, a)))
    rule219 = ReplacementRule(pattern219, lambda x, n, c, a : Add(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x))))
    rubi.add(rule219)

    pattern220 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, n, a, b)))
    rule220 = ReplacementRule(pattern220, lambda x, n, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Mul(Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule220)

    pattern221 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, n, d)))
    rule221 = ReplacementRule(pattern221, lambda x, n, d : Add(Mul(zoo, Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(zoo, Int(Mul(Pow(x, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1)))), x))))
    rubi.add(rule221)

    pattern222 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, n, c)))
    rule222 = ReplacementRule(pattern222, lambda x, n, c : Add(Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(c, Integer(-1)), Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x))))
    rubi.add(rule222)

    pattern223 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, n, b)))
    rule223 = ReplacementRule(pattern223, lambda x, n, b : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(x, Add(n, Integer(1)))), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule223)

    pattern224 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, n, a)))
    rule224 = ReplacementRule(pattern224, lambda x, n, a : Add(Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Integer(-1))), x))))
    rubi.add(rule224)

    pattern225 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x, n)))
    rule225 = ReplacementRule(pattern225, lambda x, n : Add(Mul(zoo, Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(zoo, Int(Mul(Pow(x, Integer(-1)), Pow(x, Add(n, Integer(1)))), x))))
    rubi.add(rule225)

    pattern226 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (x, c, a, d, b)))
    rule226 = ReplacementRule(pattern226, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule226)

    pattern227 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule227 = ReplacementRule(pattern227, lambda x, c, d, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule227)

    pattern228 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule228 = ReplacementRule(pattern228, lambda x, c, a, d : With(List(Set(q, Rt(Add(Mul(Integer(-1), a, d), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule228)

    pattern229 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (x, a, d, b)))
    rule229 = ReplacementRule(pattern229, lambda x, a, d, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule229)

    pattern230 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, c, a, b)))
    rule230 = ReplacementRule(pattern230, lambda x, c, a, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule230)

    pattern231 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule231 = ReplacementRule(pattern231, lambda x, c, d : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule231)

    pattern232 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule232 = ReplacementRule(pattern232, lambda x, d, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule232)

    pattern233 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule233 = ReplacementRule(pattern233, lambda x, c, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule233)

    pattern234 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule234 = ReplacementRule(pattern234, lambda x, a, d : With(List(Set(q, Rt(Mul(Integer(-1), a, d), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule234)

    pattern235 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(PosQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule235 = ReplacementRule(pattern235, lambda x, c, a : With(List(Set(q, Rt(Add(Mul(Integer(-1), a), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule235)

    pattern236 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule236 = ReplacementRule(pattern236, lambda x, a, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule236)

    pattern237 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule237 = ReplacementRule(pattern237, lambda x, d : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule237)

    pattern238 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule238 = ReplacementRule(pattern238, lambda x, c : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule238)

    pattern239 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-4/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule239 = ReplacementRule(pattern239, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule239)

    pattern240 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule240 = ReplacementRule(pattern240, lambda x, a : With(List(Set(q, Rt(Mul(Integer(-1), a), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule240)

    pattern241 = Pattern(Int(Pow(x_, Integer(-4/3)), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule241 = ReplacementRule(pattern241, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule241)

    pattern242 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (x, c, a, d, b)))
    rule242 = ReplacementRule(pattern242, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(a, d), Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule242)

    pattern243 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule243 = ReplacementRule(pattern243, lambda x, c, d, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule243)

    pattern244 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NegQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule244 = ReplacementRule(pattern244, lambda x, c, a, d : With(List(Set(q, Rt(Add(Mul(a, d), Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule244)

    pattern245 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (x, a, d, b)))
    rule245 = ReplacementRule(pattern245, lambda x, a, d, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule245)

    pattern246 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, c, a, b)))
    rule246 = ReplacementRule(pattern246, lambda x, c, a, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(a, Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule246)

    pattern247 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule247 = ReplacementRule(pattern247, lambda x, c, d : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule247)

    pattern248 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule248 = ReplacementRule(pattern248, lambda x, d, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule248)

    pattern249 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule249 = ReplacementRule(pattern249, lambda x, c, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule249)

    pattern250 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule250 = ReplacementRule(pattern250, lambda x, a, d : With(List(Set(q, Rt(Mul(a, d), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule250)

    pattern251 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(NegQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule251 = ReplacementRule(pattern251, lambda x, c, a : With(List(Set(q, Rt(Add(a, Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule251)

    pattern252 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule252 = ReplacementRule(pattern252, lambda x, a, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule252)

    pattern253 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule253 = ReplacementRule(pattern253, lambda x, d : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule253)

    pattern254 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule254 = ReplacementRule(pattern254, lambda x, c : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule254)

    pattern255 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-4/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule255 = ReplacementRule(pattern255, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule255)

    pattern256 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule256 = ReplacementRule(pattern256, lambda x, a : With(List(Set(q, Rt(a, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule256)

    pattern257 = Pattern(Int(Pow(x_, Integer(-4/3)), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule257 = ReplacementRule(pattern257, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule257)

    pattern258 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (x, c, a, d, b)))
    rule258 = ReplacementRule(pattern258, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule258)

    pattern259 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule259 = ReplacementRule(pattern259, lambda x, c, d, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule259)

    pattern260 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), PosQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule260 = ReplacementRule(pattern260, lambda x, c, a, d : With(List(Set(q, Rt(Add(Mul(Integer(-1), a, d), c), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule260)

    pattern261 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (x, a, d, b)))
    rule261 = ReplacementRule(pattern261, lambda x, a, d, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule261)

    pattern262 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, c, a, b)))
    rule262 = ReplacementRule(pattern262, lambda x, c, a, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule262)

    pattern263 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule263 = ReplacementRule(pattern263, lambda x, c, d : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule263)

    pattern264 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule264 = ReplacementRule(pattern264, lambda x, d, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule264)

    pattern265 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule265 = ReplacementRule(pattern265, lambda x, c, b : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule265)

    pattern266 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule266 = ReplacementRule(pattern266, lambda x, a, d : With(List(Set(q, Rt(Mul(Integer(-1), a, d), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule266)

    pattern267 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(PosQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule267 = ReplacementRule(pattern267, lambda x, c, a : With(List(Set(q, Rt(Add(Mul(Integer(-1), a), c), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule267)

    pattern268 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule268 = ReplacementRule(pattern268, lambda x, a, b : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule268)

    pattern269 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule269 = ReplacementRule(pattern269, lambda x, d : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule269)

    pattern270 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(PosQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule270 = ReplacementRule(pattern270, lambda x, c : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule270)

    pattern271 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-5/3))), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule271 = ReplacementRule(pattern271, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule271)

    pattern272 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(PosQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule272 = ReplacementRule(pattern272, lambda x, a : With(List(Set(q, Rt(Mul(Integer(-1), a), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule272)

    pattern273 = Pattern(Int(Pow(x_, Integer(-5/3)), x_), cons(And(PosQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule273 = ReplacementRule(pattern273, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule273)

    pattern274 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))))), (x, c, a, d, b)))
    rule274 = ReplacementRule(pattern274, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(a, d), Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule274)

    pattern275 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule275 = ReplacementRule(pattern275, lambda x, c, d, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule275)

    pattern276 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NegQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule276 = ReplacementRule(pattern276, lambda x, c, a, d : With(List(Set(q, Rt(Add(Mul(a, d), Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule276)

    pattern277 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), d_), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_))), (x, a, d, b)))
    rule277 = ReplacementRule(pattern277, lambda x, a, d, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule277)

    pattern278 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_))))), (x, c, a, b)))
    rule278 = ReplacementRule(pattern278, lambda x, c, a, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(a, Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule278)

    pattern279 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule279 = ReplacementRule(pattern279, lambda x, c, d : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule279)

    pattern280 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule280 = ReplacementRule(pattern280, lambda x, d, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule280)

    pattern281 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule281 = ReplacementRule(pattern281, lambda x, c, b : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule281)

    pattern282 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule282 = ReplacementRule(pattern282, lambda x, a, d : With(List(Set(q, Rt(Mul(a, d), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule282)

    pattern283 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule283 = ReplacementRule(pattern283, lambda x, c, a : With(List(Set(q, Rt(Add(a, Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule283)

    pattern284 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1))))), (x, a, b)))
    rule284 = ReplacementRule(pattern284, lambda x, a, b : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule284)

    pattern285 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule285 = ReplacementRule(pattern285, lambda x, d : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule285)

    pattern286 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule286 = ReplacementRule(pattern286, lambda x, c : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule286)

    pattern287 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-5/3))), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule287 = ReplacementRule(pattern287, lambda x, b : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule287)

    pattern288 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(NegQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule288 = ReplacementRule(pattern288, lambda x, a : With(List(Set(q, Rt(a, Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule288)

    pattern289 = Pattern(Int(Pow(x_, Integer(-5/3)), x_), cons(And(NegQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule289 = ReplacementRule(pattern289, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule289)

    pattern290 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, n, c, a, d, b)))
    rule290 = ReplacementRule(pattern290, lambda x, n, c, a, d, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule290)

    pattern291 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, n, c, d, b)))
    rule291 = ReplacementRule(pattern291, lambda x, n, c, d, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule291)

    pattern292 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, n, c, a, d)))
    rule292 = ReplacementRule(pattern292, lambda x, n, c, a, d : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule292)

    pattern293 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, n, a, d, b)))
    rule293 = ReplacementRule(pattern293, lambda x, n, a, d, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule293)

    pattern294 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, c, a, b)))
    rule294 = ReplacementRule(pattern294, lambda x, n, c, a, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule294)

    pattern295 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, n, c, d)))
    rule295 = ReplacementRule(pattern295, lambda x, n, c, d : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule295)

    pattern296 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, n, d, b)))
    rule296 = ReplacementRule(pattern296, lambda x, n, d, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(b, Integer(-1)), Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule296)

    pattern297 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, n, c, b)))
    rule297 = ReplacementRule(pattern297, lambda x, n, c, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule297)

    pattern298 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, n, a, d)))
    rule298 = ReplacementRule(pattern298, lambda x, n, a, d : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Pow(x, p)), Integer(-1))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule298)

    pattern299 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, n, c, a)))
    rule299 = ReplacementRule(pattern299, lambda x, n, c, a : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule299)

    pattern300 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, n, a, b)))
    rule300 = ReplacementRule(pattern300, lambda x, n, a, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule300)

    pattern301 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, n, d)))
    rule301 = ReplacementRule(pattern301, lambda x, n, d : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule301)

    pattern302 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, n, c)))
    rule302 = ReplacementRule(pattern302, lambda x, n, c : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule302)

    pattern303 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, n, b)))
    rule303 = ReplacementRule(pattern303, lambda x, n, b : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(b, Integer(-1)), Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule303)

    pattern304 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, n, a)))
    rule304 = ReplacementRule(pattern304, lambda x, n, a : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Pow(x, p)), Integer(-1))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule304)

    pattern305 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x, n)))
    rule305 = ReplacementRule(pattern305, lambda x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule305)

    pattern306 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(c_, d_, n_), x_)), (d, c, n, x)))
    rule306 = ReplacementRule(pattern306, lambda d, c, n, x : Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule306)

    pattern307 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(c_, Integer(1), n_), x_)), (c, n, x)))
    rule307 = ReplacementRule(pattern307, lambda c, n, x : Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule307)

    pattern308 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(a_, b_, c_, d_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, n, c, a, d, b)))
    rule308 = ReplacementRule(pattern308, lambda x, n, c, a, d, b : Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(b, Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))))
    rubi.add(rule308)

    pattern309 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, n_), x_)), (x, n, c, a, d)))
    rule309 = ReplacementRule(pattern309, lambda x, n, c, a, d : Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))))
    rubi.add(rule309)

    pattern310 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_, n_), x_)), (x, n, a, d, b)))
    rule310 = ReplacementRule(pattern310, lambda x, n, a, d, b : Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), b, x)))))
    rubi.add(rule310)

    pattern311 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(a_, b_, c_, Integer(1), n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, n, c, a, b)))
    rule311 = ReplacementRule(pattern311, lambda x, n, c, a, b : Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(c, x))))))
    rubi.add(rule311)

    pattern312 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_, n_), x_)), (x, n, a, d)))
    rule312 = ReplacementRule(pattern312, lambda x, n, a, d : Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), x)))))
    rubi.add(rule312)

    pattern313 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), n_), x_)), (x, n, c, a)))
    rule313 = ReplacementRule(pattern313, lambda x, n, c, a : Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(c, x))))))
    rubi.add(rule313)

    pattern314 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1), n_), x_)), (x, n, a, b)))
    rule314 = ReplacementRule(pattern314, lambda x, n, a, b : Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), b, x)))))
    rubi.add(rule314)

    pattern315 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), n_), x_)), (x, n, a)))
    rule315 = ReplacementRule(pattern315, lambda x, n, a : Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), x)))))
    rubi.add(rule315)

    pattern316 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, b_, c_, d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, a, d, b)))
    rule316 = ReplacementRule(pattern316, lambda x, m, n, c, a, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule316)

    pattern317 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, d, b)))
    rule317 = ReplacementRule(pattern317, lambda x, m, n, c, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule317)

    pattern318 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, a, d)))
    rule318 = ReplacementRule(pattern318, lambda x, m, n, c, a, d : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule318)

    pattern319 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, a, d, b)))
    rule319 = ReplacementRule(pattern319, lambda x, m, n, a, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule319)

    pattern320 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, b_, c_, Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, a, b)))
    rule320 = ReplacementRule(pattern320, lambda x, m, n, c, a, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule320)

    pattern321 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, d)))
    rule321 = ReplacementRule(pattern321, lambda x, m, n, c, d : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule321)

    pattern322 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), b_, Integer(0), d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, d, b)))
    rule322 = ReplacementRule(pattern322, lambda x, m, n, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule322)

    pattern323 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, b)))
    rule323 = ReplacementRule(pattern323, lambda x, m, n, c, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule323)

    pattern324 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, a, d)))
    rule324 = ReplacementRule(pattern324, lambda x, m, n, a, d : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule324)

    pattern325 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c, a)))
    rule325 = ReplacementRule(pattern325, lambda x, m, n, c, a : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule325)

    pattern326 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, a, b)))
    rule326 = ReplacementRule(pattern326, lambda x, m, n, a, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule326)

    pattern327 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, d)))
    rule327 = ReplacementRule(pattern327, lambda x, m, n, d : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule327)

    pattern328 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, c)))
    rule328 = ReplacementRule(pattern328, lambda x, m, n, c : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule328)

    pattern329 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, b)))
    rule329 = ReplacementRule(pattern329, lambda x, m, n, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Mul(b, x), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule329)

    pattern330 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n, a)))
    rule330 = ReplacementRule(pattern330, lambda x, m, n, a : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule330)

    pattern331 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (x, m, n)))
    rule331 = ReplacementRule(pattern331, lambda x, m, n : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(-1)))), x)), Mul(Pow(x, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule331)

    pattern332 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, b_, c_, d_), x_), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, a, d, b)))
    rule332 = ReplacementRule(pattern332, lambda x, m, n, c, a, d, b : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule332)

    pattern333 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, d, b)))
    rule333 = ReplacementRule(pattern333, lambda x, m, n, c, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule333)

    pattern334 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, a, d)))
    rule334 = ReplacementRule(pattern334, lambda x, m, n, c, a, d : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))
    rubi.add(rule334)

    pattern335 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, a, d, b)))
    rule335 = ReplacementRule(pattern335, lambda x, m, n, a, d, b : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule335)

    pattern336 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(a_, b_, c_, Integer(1)), x_), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, a, b)))
    rule336 = ReplacementRule(pattern336, lambda x, m, n, c, a, b : Add(Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule336)

    pattern337 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, d)))
    rule337 = ReplacementRule(pattern337, lambda x, m, n, c, d : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule337)

    pattern338 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), b_, Integer(0), d_), x_), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, d, b)))
    rule338 = ReplacementRule(pattern338, lambda x, m, n, d, b : Add(Mul(zoo, Pow(b, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n)), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule338)

    pattern339 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, b)))
    rule339 = ReplacementRule(pattern339, lambda x, m, n, c, b : Add(Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule339)

    pattern340 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, a, d)))
    rule340 = ReplacementRule(pattern340, lambda x, m, n, a, d : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule340)

    pattern341 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c, a)))
    rule341 = ReplacementRule(pattern341, lambda x, m, n, c, a : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule341)

    pattern342 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, a, b)))
    rule342 = ReplacementRule(pattern342, lambda x, m, n, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule342)

    pattern343 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, d)))
    rule343 = ReplacementRule(pattern343, lambda x, m, n, d : Add(Mul(zoo, d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n)), x)), Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule343)

    pattern344 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (x, m, n, c)))
    rule344 = ReplacementRule(pattern344, lambda x, m, n, c : Add(Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule344)

    pattern345 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, b)))
    rule345 = ReplacementRule(pattern345, lambda x, m, n, b : Add(Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(b, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1)))), x))))
    rubi.add(rule345)

    pattern346 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n, a)))
    rule346 = ReplacementRule(pattern346, lambda x, m, n, a : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1)))), x))))
    rubi.add(rule346)

    pattern347 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (x, m, n)))
    rule347 = ReplacementRule(pattern347, lambda x, m, n : Add(Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(x, Add(m, Integer(1)))), x))))
    rubi.add(rule347)

    pattern348 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, a, d, b)))
    rule348 = ReplacementRule(pattern348, lambda x, m, n, c, a, d, b : Add(Mul(Pow(b, Integer(-1)), n, Add(Mul(Integer(-1), a, d), Mul(b, c)), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule348)

    pattern349 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, d_), x_), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, d, b)))
    rule349 = ReplacementRule(pattern349, lambda x, m, n, c, d, b : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule349)

    pattern350 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, a, d)))
    rule350 = ReplacementRule(pattern350, lambda x, m, n, c, a, d : Add(Mul(n, Add(Mul(Integer(-1), a, d), c), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule350)

    pattern351 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), d_), x_), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, a, d, b)))
    rule351 = ReplacementRule(pattern351, lambda x, m, n, a, d, b : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), m)), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule351)

    pattern352 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, Integer(1)), x_), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, a, b)))
    rule352 = ReplacementRule(pattern352, lambda x, m, n, c, a, b : Add(Mul(Pow(b, Integer(-1)), n, Add(Mul(Integer(-1), a), Mul(b, c)), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule352)

    pattern353 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, d)))
    rule353 = ReplacementRule(pattern353, lambda x, m, n, c, d : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule353)

    pattern354 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, d, b)))
    rule354 = ReplacementRule(pattern354, lambda x, m, n, d, b : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule354)

    pattern355 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, b)))
    rule355 = ReplacementRule(pattern355, lambda x, m, n, c, b : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule355)

    pattern356 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, a, d)))
    rule356 = ReplacementRule(pattern356, lambda x, m, n, a, d : Add(Mul(Integer(-1), a, d, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), m)), x)), Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule356)

    pattern357 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c, a)))
    rule357 = ReplacementRule(pattern357, lambda x, m, n, c, a : Add(Mul(n, Add(Mul(Integer(-1), a), c), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule357)

    pattern358 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, a, b)))
    rule358 = ReplacementRule(pattern358, lambda x, m, n, a, b : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), m)), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule358)

    pattern359 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, d)))
    rule359 = ReplacementRule(pattern359, lambda x, m, n, d : Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule359)

    pattern360 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, c)))
    rule360 = ReplacementRule(pattern360, lambda x, m, n, c : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule360)

    pattern361 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, b)))
    rule361 = ReplacementRule(pattern361, lambda x, m, n, b : Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule361)

    pattern362 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n, a)))
    rule362 = ReplacementRule(pattern362, lambda x, m, n, a : Add(Mul(Integer(-1), a, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), m)), x)), Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule362)

    pattern363 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (x, m, n)))
    rule363 = ReplacementRule(pattern363, lambda x, m, n : Mul(Pow(x, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule363)

    pattern364 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(b_, d_)), FreeQ(List(a_, b_, c_, d_), x_)), (x, c, a, d, b)))
    rule364 = ReplacementRule(pattern364, lambda x, c, a, d, b : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))), Mul(Integer(-1), b, x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule364)

    pattern365 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(d_, Integer(1))), FreeQ(List(a_, Integer(1), c_, d_), x_)), (x, c, a, d)))
    rule365 = ReplacementRule(pattern365, lambda x, c, a, d : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(x, Integer(2))), Mul(Integer(-1), x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule365)

    pattern366 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(b_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule366 = ReplacementRule(pattern366, lambda x, a, d, b : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, b, x), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule366)

    pattern367 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(b_, Integer(1))), FreeQ(List(a_, b_, c_, Integer(1)), x_)), (x, c, a, b)))
    rule367 = ReplacementRule(pattern367, lambda x, c, a, b : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))), Mul(Integer(-1), b, x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule367)

    pattern368 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(d_, Integer(1))), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule368 = ReplacementRule(pattern368, lambda x, a, d : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, x), Mul(Integer(-1), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule368)

    pattern369 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(ZeroQ(Integer(2)), PositiveQ(Add(a_, c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule369 = ReplacementRule(pattern369, lambda x, c, a : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(x, Integer(2))), Mul(Integer(-1), x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule369)

    pattern370 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Add(b_, Integer(1))), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule370 = ReplacementRule(pattern370, lambda x, a, b : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, b, x), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule370)

    pattern371 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(a_), ZeroQ(Integer(2)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule371 = ReplacementRule(pattern371, lambda x, a : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, x), Mul(Integer(-1), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule371)

    pattern372 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(b_), FreeQ(List(a_, b_, c_, d_), x_), PositiveQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule372 = ReplacementRule(pattern372, lambda x, c, a, d, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(b, c), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule372)

    pattern373 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(b_), PositiveQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule373 = ReplacementRule(pattern373, lambda x, c, d, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(b, c), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule373)

    pattern374 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), FreeQ(List(a_, Integer(1), c_, d_), x_), PositiveQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule374 = ReplacementRule(pattern374, lambda x, c, a, d : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), c, Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule374)

    pattern375 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(b_), PositiveQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule375 = ReplacementRule(pattern375, lambda x, a, d, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule375)

    pattern376 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(b_), FreeQ(List(a_, b_, c_, Integer(1)), x_), PositiveQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule376 = ReplacementRule(pattern376, lambda x, c, a, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Mul(b, c), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule376)

    pattern377 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), PositiveQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule377 = ReplacementRule(pattern377, lambda x, c, d : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule377)

    pattern378 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(0)), PositiveQ(b_), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule378 = ReplacementRule(pattern378, lambda x, d, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Mul(d, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule378)

    pattern379 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(b_), PositiveQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule379 = ReplacementRule(pattern379, lambda x, c, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(b, c), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule379)

    pattern380 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), PositiveQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule380 = ReplacementRule(pattern380, lambda x, a, d : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule380)

    pattern381 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), PositiveQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule381 = ReplacementRule(pattern381, lambda x, c, a : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule381)

    pattern382 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(PositiveQ(b_), PositiveQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule382 = ReplacementRule(pattern382, lambda x, a, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule382)

    pattern383 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(0)), PositiveQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule383 = ReplacementRule(pattern383, lambda x, d : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Mul(d, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule383)

    pattern384 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), PositiveQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule384 = ReplacementRule(pattern384, lambda x, c : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule384)

    pattern385 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(b_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(0)), PositiveQ(b_), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule385 = ReplacementRule(pattern385, lambda x, b : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Pow(x, Integer(2))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule385)

    pattern386 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(PositiveQ(Integer(1)), PositiveQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule386 = ReplacementRule(pattern386, lambda x, a : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule386)

    pattern387 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(PositiveQ(Integer(0)), PositiveQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule387 = ReplacementRule(pattern387, lambda x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Pow(x, Integer(2))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule387)

    pattern388 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(ZeroQ(Add(b_, Mul(Integer(-1), d_))), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule388 = ReplacementRule(pattern388, lambda x, c, a, d, b : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule388)

    pattern389 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(b_, c_)), ZeroQ(Add(b_, Mul(Integer(-1), d_))), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule389 = ReplacementRule(pattern389, lambda x, c, d, b : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule389)

    pattern390 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(ZeroQ(Add(Mul(Integer(-1), d_), Integer(1))), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule390 = ReplacementRule(pattern390, lambda x, c, a, d : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule390)

    pattern391 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(ZeroQ(Add(b_, Integer(-1))), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule391 = ReplacementRule(pattern391, lambda x, c, a, b : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule391)

    pattern392 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(c_), ZeroQ(Add(Mul(Integer(-1), d_), Integer(1))), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule392 = ReplacementRule(pattern392, lambda x, c, d : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule392)

    pattern393 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(b_, c_)), ZeroQ(Add(b_, Integer(-1))), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule393 = ReplacementRule(pattern393, lambda x, c, b : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule393)

    pattern394 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(ZeroQ(Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule394 = ReplacementRule(pattern394, lambda x, c, a : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule394)

    pattern395 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(c_), ZeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule395 = ReplacementRule(pattern395, lambda x, c : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule395)

    pattern396 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule396 = ReplacementRule(pattern396, lambda x, c, a, d, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule396)

    pattern397 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule397 = ReplacementRule(pattern397, lambda x, c, d, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule397)

    pattern398 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule398 = ReplacementRule(pattern398, lambda x, c, a, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule398)

    pattern399 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule399 = ReplacementRule(pattern399, lambda x, a, d, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(Mul(d, x)), Integer(-1)), Sqrt(Add(a, Mul(b, x)))))))
    rubi.add(rule399)

    pattern400 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule400 = ReplacementRule(pattern400, lambda x, c, a, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule400)

    pattern401 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule401 = ReplacementRule(pattern401, lambda x, c, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule401)

    pattern402 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule402 = ReplacementRule(pattern402, lambda x, d, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Mul(d, x)), Integer(-1))))))
    rubi.add(rule402)

    pattern403 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule403 = ReplacementRule(pattern403, lambda x, c, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule403)

    pattern404 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule404 = ReplacementRule(pattern404, lambda x, a, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Pow(Sqrt(Mul(d, x)), Integer(-1)), Sqrt(Add(a, x))))))
    rubi.add(rule404)

    pattern405 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule405 = ReplacementRule(pattern405, lambda x, c, a : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule405)

    pattern406 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule406 = ReplacementRule(pattern406, lambda x, a, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Add(a, Mul(b, x)))))))
    rubi.add(rule406)

    pattern407 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule407 = ReplacementRule(pattern407, lambda x, d : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Mul(d, x)), Integer(-1))))))
    rubi.add(rule407)

    pattern408 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule408 = ReplacementRule(pattern408, lambda x, c : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule408)

    pattern409 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(b_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule409 = ReplacementRule(pattern409, lambda x, b : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Mul(b, x))))))
    rubi.add(rule409)

    pattern410 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule410 = ReplacementRule(pattern410, lambda x, a : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Add(a, x))))))
    rubi.add(rule410)

    pattern411 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule411 = ReplacementRule(pattern411, lambda x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Integer(1))))
    rubi.add(rule411)

    pattern412 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, c, a, d, b)))
    rule412 = ReplacementRule(pattern412, lambda x, m, c, a, d, b : Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), Mul(b, c)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), Mul(b, c)))), m), x)))
    rubi.add(rule412)

    pattern413 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, m, c, d, b)))
    rule413 = ReplacementRule(pattern413, lambda x, m, c, d, b : Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(b, c, x), Mul(b, d, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(b, c, x), Mul(b, d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule413)

    pattern414 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, m, c, a, d)))
    rule414 = ReplacementRule(pattern414, lambda x, m, c, a, d : Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), c))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), c))), m), x)))
    rubi.add(rule414)

    pattern415 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, m, c, a, b)))
    rule415 = ReplacementRule(pattern415, lambda x, m, c, a, b : Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), m), Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2))), Mul(x, Add(a, Mul(b, c)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2))), Mul(x, Add(a, Mul(b, c)))), m), x)))
    rubi.add(rule415)

    pattern416 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, m, c, d)))
    rule416 = ReplacementRule(pattern416, lambda x, m, c, d : Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(c, x), Mul(d, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(c, x), Mul(d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule416)

    pattern417 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, m, c, b)))
    rule417 = ReplacementRule(pattern417, lambda x, m, c, b : Mul(Pow(Mul(b, x), m), Pow(Add(c, x), m), Pow(Add(Mul(b, c, x), Mul(b, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(b, c, x), Mul(b, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule417)

    pattern418 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, c, a)))
    rule418 = ReplacementRule(pattern418, lambda x, m, c, a : Mul(Pow(Add(a, x), m), Pow(Add(c, x), m), Pow(Add(Mul(a, c), Pow(x, Integer(2)), Mul(x, Add(a, c))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Pow(x, Integer(2)), Mul(x, Add(a, c))), m), x)))
    rubi.add(rule418)

    pattern419 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), m_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, m, c)))
    rule419 = ReplacementRule(pattern419, lambda x, m, c : Mul(Pow(x, m), Pow(Add(c, x), m), Pow(Add(Mul(c, x), Pow(x, Integer(2))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(c, x), Pow(x, Integer(2))), m), x)))
    rubi.add(rule419)

    pattern420 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(PosQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule420 = ReplacementRule(pattern420, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule420)

    pattern421 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NonzeroQ(Mul(b_, c_)), PosQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule421 = ReplacementRule(pattern421, lambda x, c, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule421)

    pattern422 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(PosQ(d_), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule422 = ReplacementRule(pattern422, lambda x, c, a, d : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule422)

    pattern423 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), PosQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule423 = ReplacementRule(pattern423, lambda x, a, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule423)

    pattern424 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(PosQ(Pow(b_, Integer(-1))), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule424 = ReplacementRule(pattern424, lambda x, c, a, b : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule424)

    pattern425 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NonzeroQ(c_), PosQ(d_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule425 = ReplacementRule(pattern425, lambda x, c, d : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule425)

    pattern426 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NonzeroQ(Integer(0)), PosQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule426 = ReplacementRule(pattern426, lambda x, d, b : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule426)

    pattern427 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NonzeroQ(Mul(b_, c_)), PosQ(Pow(b_, Integer(-1))), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule427 = ReplacementRule(pattern427, lambda x, c, b : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule427)

    pattern428 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), cons(And(PosQ(d_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule428 = ReplacementRule(pattern428, lambda x, a, d : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule428)

    pattern429 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(PosQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule429 = ReplacementRule(pattern429, lambda x, c, a : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule429)

    pattern430 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), PosQ(Pow(b_, Integer(-1))), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule430 = ReplacementRule(pattern430, lambda x, a, b : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule430)

    pattern431 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NonzeroQ(Integer(0)), PosQ(d_), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule431 = ReplacementRule(pattern431, lambda x, d : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule431)

    pattern432 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NonzeroQ(c_), PosQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule432 = ReplacementRule(pattern432, lambda x, c : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule432)

    pattern433 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Mul(b_, x_), Integer(-1/3))), x_), cons(And(NonzeroQ(Integer(0)), PosQ(Pow(b_, Integer(-1))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule433 = ReplacementRule(pattern433, lambda x, b : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule433)

    pattern434 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), cons(And(PosQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule434 = ReplacementRule(pattern434, lambda x, a : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule434)

    pattern435 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(NonzeroQ(Integer(0)), PosQ(Integer(1)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule435 = ReplacementRule(pattern435, lambda x : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(q, Integer(-1)))))))
    rubi.add(rule435)

    pattern436 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NegQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, c, a, d, b)))
    rule436 = ReplacementRule(pattern436, lambda x, c, a, d, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule436)

    pattern437 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NonzeroQ(Mul(b_, c_)), NegQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, c, d, b)))
    rule437 = ReplacementRule(pattern437, lambda x, c, d, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule437)

    pattern438 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NegQ(d_), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, c, a, d)))
    rule438 = ReplacementRule(pattern438, lambda x, c, a, d : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule438)

    pattern439 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), NegQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, a, d, b)))
    rule439 = ReplacementRule(pattern439, lambda x, a, d, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(1)))))))
    rubi.add(rule439)

    pattern440 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(Pow(b_, Integer(-1))), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, c, a, b)))
    rule440 = ReplacementRule(pattern440, lambda x, c, a, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule440)

    pattern441 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), cons(And(NegQ(d_), NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, c, d)))
    rule441 = ReplacementRule(pattern441, lambda x, c, d : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule441)

    pattern442 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NonzeroQ(Integer(0)), NegQ(Mul(Pow(b_, Integer(-1)), d_)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, d, b)))
    rule442 = ReplacementRule(pattern442, lambda x, d, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule442)

    pattern443 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(Pow(b_, Integer(-1))), NonzeroQ(Mul(b_, c_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, c, b)))
    rule443 = ReplacementRule(pattern443, lambda x, c, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule443)

    pattern444 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), cons(And(NegQ(d_), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, a, d)))
    rule444 = ReplacementRule(pattern444, lambda x, a, d : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule444)

    pattern445 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, c, a)))
    rule445 = ReplacementRule(pattern445, lambda x, c, a : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule445)

    pattern446 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), cons(And(NegQ(Pow(b_, Integer(-1))), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, a, b)))
    rule446 = ReplacementRule(pattern446, lambda x, a, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(1)))))))
    rubi.add(rule446)

    pattern447 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), cons(And(NegQ(d_), NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, d)))
    rule447 = ReplacementRule(pattern447, lambda x, d : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule447)

    pattern448 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), cons(And(NegQ(Integer(1)), NonzeroQ(c_), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, c)))
    rule448 = ReplacementRule(pattern448, lambda x, c : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule448)

    pattern449 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Mul(b_, x_), Integer(-1/3))), x_), cons(And(NonzeroQ(Integer(0)), NegQ(Pow(b_, Integer(-1))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, b)))
    rule449 = ReplacementRule(pattern449, lambda x, b : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule449)

    pattern450 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), cons(And(NegQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, a)))
    rule450 = ReplacementRule(pattern450, lambda x, a : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule450)

    pattern451 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(NegQ(Integer(1)), NonzeroQ(Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x,)))
    rule451 = ReplacementRule(pattern451, lambda x : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(q, Integer(1)))))))
    rubi.add(rule451)

    pattern452 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule452 = ReplacementRule(pattern452, lambda x, m, n, c, a, d, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule452)

    pattern453 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, d_), x_)), (x, m, n, c, d, b)))
    rule453 = ReplacementRule(pattern453, lambda x, m, n, c, d, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule453)

    pattern454 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (x, m, n, c, a, d)))
    rule454 = ReplacementRule(pattern454, lambda x, m, n, c, a, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule454)

    pattern455 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), d_), x_)), (x, m, n, a, d, b)))
    rule455 = ReplacementRule(pattern455, lambda x, m, n, a, d, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule455)

    pattern456 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, c_, Integer(1)), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, m, n, c, a, b)))
    rule456 = ReplacementRule(pattern456, lambda x, m, n, c, a, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule456)

    pattern457 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_)), (x, m, n, c, d)))
    rule457 = ReplacementRule(pattern457, lambda x, m, n, c, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule457)

    pattern458 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_)), (x, m, n, d, b)))
    rule458 = ReplacementRule(pattern458, lambda x, m, n, d, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule458)

    pattern459 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_)), (x, m, n, c, b)))
    rule459 = ReplacementRule(pattern459, lambda x, m, n, c, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule459)

    pattern460 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_)), (x, m, n, a, d)))
    rule460 = ReplacementRule(pattern460, lambda x, m, n, a, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule460)

    pattern461 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_)), (x, m, n, c, a)))
    rule461 = ReplacementRule(pattern461, lambda x, m, n, c, a : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule461)

    pattern462 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_)), (x, m, n, a, b)))
    rule462 = ReplacementRule(pattern462, lambda x, m, n, a, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule462)

    pattern463 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_)), (x, m, n, d)))
    rule463 = ReplacementRule(pattern463, lambda x, m, n, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule463)

    pattern464 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_)), (x, m, n, c)))
    rule464 = ReplacementRule(pattern464, lambda x, m, n, c : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule464)

    pattern465 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_)), (x, m, n, b)))
    rule465 = ReplacementRule(pattern465, lambda x, m, n, b : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Mul(b, x), Pow(p, Integer(-1))))))))
    rubi.add(rule465)

    pattern466 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_)), (x, m, n, a)))
    rule466 = ReplacementRule(pattern466, lambda x, m, n, a : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule466)

    pattern467 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_)), (x, m, n)))
    rule467 = ReplacementRule(pattern467, lambda x, m, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Integer(1)))))
    rubi.add(rule467)

    pattern468 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, c_, d_), x_), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule468 = ReplacementRule(pattern468, lambda x, m, n, c, a, d, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule468)

    pattern469 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, c_, d_), x_), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_)), (x, m, n, c, d, b)))
    rule469 = ReplacementRule(pattern469, lambda x, m, n, c, d, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule469)

    pattern470 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, Integer(1), c_, d_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_)), (x, m, n, c, a, d)))
    rule470 = ReplacementRule(pattern470, lambda x, m, n, c, a, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), c, Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule470)

    pattern471 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, b_, Integer(0), d_), x_), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_)), (x, m, n, a, d, b)))
    rule471 = ReplacementRule(pattern471, lambda x, m, n, a, d, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule471)

    pattern472 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, c_, Integer(1)), x_), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, m, n, c, a, b)))
    rule472 = ReplacementRule(pattern472, lambda x, m, n, c, a, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule472)

    pattern473 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), c_, d_), x_), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_)), (x, m, n, c, d)))
    rule473 = ReplacementRule(pattern473, lambda x, m, n, c, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(d, Pow(x, p))), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule473)

    pattern474 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, Integer(0), d_), x_), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_)), (x, m, n, d, b)))
    rule474 = ReplacementRule(pattern474, lambda x, m, n, d, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), d, Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule474)

    pattern475 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, c_, Integer(1)), x_), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_)), (x, m, n, c, b)))
    rule475 = ReplacementRule(pattern475, lambda x, m, n, c, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule475)

    pattern476 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), FreeQ(List(a_, Integer(1), Integer(0), d_), x_), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_)), (x, m, n, a, d)))
    rule476 = ReplacementRule(pattern476, lambda x, m, n, a, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule476)

    pattern477 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1)), x_), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_)), (x, m, n, c, a)))
    rule477 = ReplacementRule(pattern477, lambda x, m, n, c, a : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), c, Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule477)

    pattern478 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, b_, Integer(0), Integer(1)), x_), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_)), (x, m, n, a, b)))
    rule478 = ReplacementRule(pattern478, lambda x, m, n, a, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule478)

    pattern479 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), Integer(0), d_), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_)), (x, m, n, d)))
    rule479 = ReplacementRule(pattern479, lambda x, m, n, d : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(d, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule479)

    pattern480 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), c_, Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_)), (x, m, n, c)))
    rule480 = ReplacementRule(pattern480, lambda x, m, n, c : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule480)

    pattern481 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), b_, Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_)), (x, m, n, b)))
    rule481 = ReplacementRule(pattern481, lambda x, m, n, b : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule481)

    pattern482 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(a_, Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_)), (x, m, n, a)))
    rule482 = ReplacementRule(pattern482, lambda x, m, n, a : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule482)

    pattern483 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1)), x_), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_)), (x, m, n)))
    rule483 = ReplacementRule(pattern483, lambda x, m, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Pow(x, p), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule483)

    pattern484 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule484 = ReplacementRule(pattern484, lambda x, m, n, c, a, d, b : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule484)

    pattern485 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, c_, d_, m_, n_), x_)), (x, m, n, c, d, b)))
    rule485 = ReplacementRule(pattern485, lambda x, m, n, c, d, b : Add(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule485)

    pattern486 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_)), (x, m, n, c, a, d)))
    rule486 = ReplacementRule(pattern486, lambda x, m, n, c, a, d : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))
    rubi.add(rule486)

    pattern487 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, Integer(0), d_, m_, n_), x_)), (x, m, n, a, d, b)))
    rule487 = ReplacementRule(pattern487, lambda x, m, n, a, d, b : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule487)

    pattern488 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_)), (x, m, n, c, a, b)))
    rule488 = ReplacementRule(pattern488, lambda x, m, n, c, a, b : Add(Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule488)

    pattern489 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, n_), x_)), (x, m, n, c, d)))
    rule489 = ReplacementRule(pattern489, lambda x, m, n, c, d : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule489)

    pattern490 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, Integer(0), d_, m_, n_), x_)), (x, m, n, d, b)))
    rule490 = ReplacementRule(pattern490, lambda x, m, n, d, b : Add(Mul(zoo, Pow(b, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule490)

    pattern491 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, n_), x_)), (x, m, n, c, b)))
    rule491 = ReplacementRule(pattern491, lambda x, m, n, c, b : Add(Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule491)

    pattern492 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), Integer(0), d_, m_, n_), x_)), (x, m, n, a, d)))
    rule492 = ReplacementRule(pattern492, lambda x, m, n, a, d : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule492)

    pattern493 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_)), (x, m, n, c, a)))
    rule493 = ReplacementRule(pattern493, lambda x, m, n, c, a : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule493)

    pattern494 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, b_, Integer(0), Integer(1), m_, n_), x_)), (x, m, n, a, b)))
    rule494 = ReplacementRule(pattern494, lambda x, m, n, a, b : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule494)

    pattern495 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, m_, n_), x_)), (x, m, n, d)))
    rule495 = ReplacementRule(pattern495, lambda x, m, n, d : Add(Mul(zoo, d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule495)

    pattern496 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, n_), x_)), (x, m, n, c)))
    rule496 = ReplacementRule(pattern496, lambda x, m, n, c : Add(Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule496)

    pattern497 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), m_, n_), x_)), (x, m, n, b)))
    rule497 = ReplacementRule(pattern497, lambda x, m, n, b : Add(Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(b, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Mul(b, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule497)

    pattern498 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, m, n, a)))
    rule498 = ReplacementRule(pattern498, lambda x, m, n, a : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule498)

    pattern499 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, m, n)))
    rule499 = ReplacementRule(pattern499, lambda x, m, n : Add(Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(x, Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule499)

    pattern500 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), FreeQ(List(b_, c_, d_, m_, n_), x_), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (x, m, n, c, d, b)))
    rule500 = ReplacementRule(pattern500, lambda x, m, n, c, d, b : Mul(Pow(b, Integer(-1)), Pow(c, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), d, x))))
    rubi.add(rule500)

    pattern501 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), FreeQ(List(Integer(1), c_, d_, m_, n_), x_), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (x, m, n, c, d)))
    rule501 = ReplacementRule(pattern501, lambda x, m, n, c, d : Mul(Pow(c, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), d, x))))
    rubi.add(rule501)

    pattern502 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), FreeQ(List(b_, c_, Integer(1), m_, n_), x_), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1))), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1))))))))), (x, m, n, c, b)))
    rule502 = ReplacementRule(pattern502, lambda x, m, n, c, b : Mul(Pow(b, Integer(-1)), Pow(c, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), x))))
    rubi.add(rule502)

    pattern503 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), FreeQ(List(Integer(1), c_, Integer(1), m_, n_), x_), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (x, m, n, c)))
    rule503 = ReplacementRule(pattern503, lambda x, m, n, c : Mul(Pow(c, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), x))))
    rubi.add(rule503)

    pattern504 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(b_, c_, d_, m_, n_), x_), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)))), (x, m, n, c, d, b)))
    rule504 = ReplacementRule(pattern504, lambda x, m, n, c, d, b : Mul(Pow(d, Integer(-1)), Pow(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d), Mul(Integer(-1), m)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule504)

    pattern505 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(Integer(1), c_, d_, m_, n_), x_), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)))), (x, m, n, c, d)))
    rule505 = ReplacementRule(pattern505, lambda x, m, n, c, d : Mul(Pow(d, Integer(-1)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d), Mul(Integer(-1), m)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule505)

    pattern506 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(b_, c_, Integer(1), m_, n_), x_), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)))))), (x, m, n, c, b)))
    rule506 = ReplacementRule(pattern506, lambda x, m, n, c, b : Mul(Pow(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1))), Mul(Integer(-1), m)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule506)

    pattern507 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(n_)), FreeQ(List(Integer(1), c_, Integer(1), m_, n_), x_), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)))))), (x, m, n, c)))
    rule507 = ReplacementRule(pattern507, lambda x, m, n, c : Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1))), Mul(Integer(-1), m)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule507)

    pattern508 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(b_, c_, d_, m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (x, m, n, c, d, b)))
    rule508 = ReplacementRule(pattern508, lambda x, m, n, c, d, b : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), n)), x)))
    rubi.add(rule508)

    pattern509 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(Integer(1), c_, d_, m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (x, m, n, c, d)))
    rule509 = ReplacementRule(pattern509, lambda x, m, n, c, d : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(x, m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), n)), x)))
    rubi.add(rule509)

    pattern510 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(b_, c_, Integer(1), m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1))))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (x, m, n, c, b)))
    rule510 = ReplacementRule(pattern510, lambda x, m, n, c, b : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), n)), x)))
    rubi.add(rule510)

    pattern511 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1))))), FreeQ(List(Integer(1), c_, Integer(1), m_, n_), x_), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (x, m, n, c)))
    rule511 = ReplacementRule(pattern511, lambda x, m, n, c : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(x, m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), n)), x)))
    rubi.add(rule511)

    pattern512 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(b_, c_, d_, m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)))), (x, m, n, c, d, b)))
    rule512 = ReplacementRule(pattern512, lambda x, m, n, c, d, b : Mul(Pow(Mul(b, x), FracPart(m)), Pow(Mul(Integer(-1), b, c, Pow(d, Integer(-1))), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), m), Pow(Add(c, Mul(d, x)), n)), x)))
    rubi.add(rule512)

    pattern513 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(Integer(1), c_, d_, m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)))), (x, m, n, c, d)))
    rule513 = ReplacementRule(pattern513, lambda x, m, n, c, d : Mul(Pow(x, FracPart(m)), Pow(Mul(Integer(-1), c, Pow(d, Integer(-1))), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), m), Pow(Add(c, Mul(d, x)), n)), x)))
    rubi.add(rule513)

    pattern514 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), FreeQ(List(b_, c_, Integer(1), m_, n_), x_), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)))))), (x, m, n, c, b)))
    rule514 = ReplacementRule(pattern514, lambda x, m, n, c, b : Mul(Pow(Mul(b, x), FracPart(m)), Pow(Mul(Integer(-1), b, c), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), m), Pow(Add(c, x), n)), x)))
    rubi.add(rule514)

    pattern515 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1))))), FreeQ(List(Integer(1), c_, Integer(1), m_, n_), x_)), (x, m, n, c)))
    rule515 = ReplacementRule(pattern515, lambda x, m, n, c : Mul(Pow(x, FracPart(m)), Pow(Mul(Integer(-1), c), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), m), Pow(Add(c, x), n)), x)))
    rubi.add(rule515)

    pattern516 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), FreeQ(List(a_, b_, c_, d_, m_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (x, m, n, c, a, d, b)))
    rule516 = ReplacementRule(pattern516, lambda x, m, n, c, a, d, b : Mul(Pow(b, Add(Mul(Integer(-1), n), Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), n), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, Mul(b, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))))
    rubi.add(rule516)

    pattern517 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_), x_)), (x, m, n, c, a, d)))
    rule517 = ReplacementRule(pattern517, lambda x, m, n, c, a, d : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), n), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, x), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))))
    rubi.add(rule517)

    pattern518 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), FreeQ(List(a_, b_, c_, Integer(1), m_), x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (x, m, n, c, a, b)))
    rule518 = ReplacementRule(pattern518, lambda x, m, n, c, a, b : Mul(Pow(b, Add(Mul(Integer(-1), n), Integer(-1))), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(a, Mul(b, x))))))
    rubi.add(rule518)

    pattern519 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_), x_)), (x, m, n, c, a)))
    rule519 = ReplacementRule(pattern519, lambda x, m, n, c, a : Mul(Pow(Add(Mul(Integer(-1), a), c), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(a, x)))))
    rubi.add(rule519)

    pattern520 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), PositiveQ(Mul(b_, Pow(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)), Integer(-1)))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), d_, Pow(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)), Integer(-1)))))))), (x, m, n, c, a, d, b)))
    rule520 = ReplacementRule(pattern520, lambda x, m, n, c, a, d, b : Mul(Pow(b, Integer(-1)), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), n)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, Mul(b, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))))
    rubi.add(rule520)

    pattern521 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_), PositiveQ(Pow(Add(Mul(Integer(-1), a_, d_), c_), Integer(-1))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), d_, Pow(Add(Mul(Integer(-1), a_, d_), c_), Integer(-1)))))))), (x, m, n, c, a, d)))
    rule521 = ReplacementRule(pattern521, lambda x, m, n, c, a, d : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Mul(Integer(-1), n)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, x), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))))
    rubi.add(rule521)

    pattern522 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_), PositiveQ(Mul(b_, Pow(Add(Mul(Integer(-1), a_), Mul(b_, c_)), Integer(-1)))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a_), Mul(b_, c_)), Integer(-1)))))))), (x, m, n, c, a, b)))
    rule522 = ReplacementRule(pattern522, lambda x, m, n, c, a, b : Mul(Pow(b, Integer(-1)), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(Integer(-1), n)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(a, Mul(b, x))))))
    rubi.add(rule522)

    pattern523 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), PositiveQ(Pow(Add(Mul(Integer(-1), a_), c_), Integer(-1))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a_), c_), Integer(-1)))))))), (x, m, n, c, a)))
    rule523 = ReplacementRule(pattern523, lambda x, m, n, c, a : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Mul(Integer(-1), n)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(a, x)))))
    rubi.add(rule523)

    pattern524 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), FreeQ(List(a_, b_, c_, d_, m_, n_), x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (x, m, n, c, a, d, b)))
    rule524 = ReplacementRule(pattern524, lambda x, m, n, c, a, d, b : Mul(Pow(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), IntPart(n))), Pow(Mul(b, Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(Mul(b, c, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(b, d, x, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))), n)), x)))
    rubi.add(rule524)

    pattern525 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (x, m, n, c, a, d)))
    rule525 = ReplacementRule(pattern525, lambda x, m, n, c, a, d : Mul(Pow(Mul(Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Pow(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Mul(Integer(-1), IntPart(n))), Int(Mul(Pow(Add(a, x), m), Pow(Add(Mul(c, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(d, x, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))), n)), x)))
    rubi.add(rule525)

    pattern526 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (x, m, n, c, a, b)))
    rule526 = ReplacementRule(pattern526, lambda x, m, n, c, a, b : Mul(Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(Integer(-1), IntPart(n))), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(c, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(Mul(b, c, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(b, x, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)))), n)), x)))
    rubi.add(rule526)

    pattern527 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (x, m, n, c, a)))
    rule527 = ReplacementRule(pattern527, lambda x, m, n, c, a : Mul(Pow(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(c, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Pow(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Mul(Integer(-1), IntPart(n))), Int(Mul(Pow(Add(a, x), m), Pow(Add(Mul(c, Pow(Add(Mul(Integer(-1), a), c), Integer(-1))), Mul(x, Pow(Add(Mul(Integer(-1), a), c), Integer(-1)))), n)), x)))
    rubi.add(rule527)

    pattern528 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, d_, m_, n_), x_)), (x, u, m, n, c, a, d, b)))
    rule528 = ReplacementRule(pattern528, lambda x, u, m, n, c, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule528)

    pattern529 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, d_, Integer(1), n_), x_)), (x, u, n, c, a, d, b)))
    rule529 = ReplacementRule(pattern529, lambda x, u, n, c, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule529)

    pattern530 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, d_, m_, n_), x_)), (x, u, m, n, c, d, b)))
    rule530 = ReplacementRule(pattern530, lambda x, u, m, n, c, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule530)

    pattern531 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, d_, m_, n_), x_)), (x, u, m, n, c, a, d)))
    rule531 = ReplacementRule(pattern531, lambda x, u, m, n, c, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule531)

    pattern532 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, d_, m_, Integer(1)), x_)), (x, u, m, c, a, d, b)))
    rule532 = ReplacementRule(pattern532, lambda x, u, m, c, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule532)

    pattern533 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Pow(Add(a_, Mul(b_, u_)), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), d_, m_, n_), x_)), (x, u, m, n, a, d, b)))
    rule533 = ReplacementRule(pattern533, lambda x, u, m, n, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule533)

    pattern534 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, Integer(1), m_, n_), x_)), (x, u, m, n, c, a, b)))
    rule534 = ReplacementRule(pattern534, lambda x, u, m, n, c, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule534)

    pattern535 = Pattern(Int(Mul(b_, u_, Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, d_, Integer(1), n_), x_)), (x, u, n, c, d, b)))
    rule535 = ReplacementRule(pattern535, lambda x, u, n, c, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule535)

    pattern536 = Pattern(Int(Mul(Add(a_, u_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, d_, Integer(1), n_), x_)), (x, u, n, c, a, d)))
    rule536 = ReplacementRule(pattern536, lambda x, u, n, c, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule536)

    pattern537 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, d_, Integer(1), Integer(1)), x_)), (x, u, c, a, d, b)))
    rule537 = ReplacementRule(pattern537, lambda x, u, c, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule537)

    pattern538 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Add(a_, Mul(b_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), d_, Integer(1), n_), x_)), (x, u, n, a, d, b)))
    rule538 = ReplacementRule(pattern538, lambda x, u, n, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule538)

    pattern539 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1), n_), x_)), (x, u, n, c, a, b)))
    rule539 = ReplacementRule(pattern539, lambda x, u, n, c, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule539)

    pattern540 = Pattern(Int(Mul(Pow(u_, m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, n_), x_)), (x, u, m, n, c, d)))
    rule540 = ReplacementRule(pattern540, lambda x, u, m, n, c, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule540)

    pattern541 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, d_, m_, Integer(1)), x_)), (x, u, m, c, d, b)))
    rule541 = ReplacementRule(pattern541, lambda x, u, m, c, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule541)

    pattern542 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Mul(d_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), d_, m_, n_), x_)), (x, u, m, n, d, b)))
    rule542 = ReplacementRule(pattern542, lambda x, u, m, n, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule542)

    pattern543 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, n_), x_)), (x, u, m, n, c, b)))
    rule543 = ReplacementRule(pattern543, lambda x, u, m, n, c, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule543)

    pattern544 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, d_, m_, Integer(1)), x_)), (x, u, m, c, a, d)))
    rule544 = ReplacementRule(pattern544, lambda x, u, m, c, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule544)

    pattern545 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Pow(Add(a_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), d_, m_, n_), x_)), (x, u, m, n, a, d)))
    rule545 = ReplacementRule(pattern545, lambda x, u, m, n, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule545)

    pattern546 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, n_), x_)), (x, u, m, n, c, a)))
    rule546 = ReplacementRule(pattern546, lambda x, u, m, n, c, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule546)

    pattern547 = Pattern(Int(Mul(d_, u_, Pow(Add(a_, Mul(b_, u_)), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), d_, m_, Integer(1)), x_)), (x, u, m, a, d, b)))
    rule547 = ReplacementRule(pattern547, lambda x, u, m, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule547)

    pattern548 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, Integer(1), m_, Integer(1)), x_)), (x, u, m, c, a, b)))
    rule548 = ReplacementRule(pattern548, lambda x, u, m, c, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x, u)))
    rubi.add(rule548)

    pattern549 = Pattern(Int(Mul(Pow(u_, n_), Pow(Add(a_, Mul(b_, u_)), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), Integer(1), m_, n_), x_)), (x, u, m, n, a, b)))
    rule549 = ReplacementRule(pattern549, lambda x, u, m, n, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule549)

    pattern550 = Pattern(Int(Mul(u_, Pow(Add(c_, Mul(d_, u_)), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1), n_), x_)), (x, u, n, c, d)))
    rule550 = ReplacementRule(pattern550, lambda x, u, n, c, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule550)

    pattern551 = Pattern(Int(Mul(b_, u_, Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, d_, Integer(1), Integer(1)), x_)), (x, u, c, d, b)))
    rule551 = ReplacementRule(pattern551, lambda x, u, c, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule551)

    pattern552 = Pattern(Int(Mul(b_, u_, Pow(Mul(d_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), d_, Integer(1), n_), x_)), (x, u, n, d, b)))
    rule552 = ReplacementRule(pattern552, lambda x, u, n, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule552)

    pattern553 = Pattern(Int(Mul(b_, u_, Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1), n_), x_)), (x, u, n, c, b)))
    rule553 = ReplacementRule(pattern553, lambda x, u, n, c, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule553)

    pattern554 = Pattern(Int(Mul(Add(a_, u_), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, d_, Integer(1), Integer(1)), x_)), (x, u, c, a, d)))
    rule554 = ReplacementRule(pattern554, lambda x, u, c, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule554)

    pattern555 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Add(a_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1), n_), x_)), (x, u, n, a, d)))
    rule555 = ReplacementRule(pattern555, lambda x, u, n, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x, u)))
    rubi.add(rule555)

    pattern556 = Pattern(Int(Mul(Add(a_, u_), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1), n_), x_)), (x, u, n, c, a)))
    rule556 = ReplacementRule(pattern556, lambda x, u, n, c, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule556)

    pattern557 = Pattern(Int(Mul(d_, u_, Add(a_, Mul(b_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), d_, Integer(1), Integer(1)), x_)), (x, u, a, d, b)))
    rule557 = ReplacementRule(pattern557, lambda x, u, a, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule557)

    pattern558 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, c_, Integer(1), Integer(1), Integer(1)), x_)), (x, u, c, a, b)))
    rule558 = ReplacementRule(pattern558, lambda x, u, c, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x, u)))
    rubi.add(rule558)

    pattern559 = Pattern(Int(Mul(Pow(u_, n_), Add(a_, Mul(b_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1), n_), x_)), (x, u, n, a, b)))
    rule559 = ReplacementRule(pattern559, lambda x, u, n, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule559)

    pattern560 = Pattern(Int(Mul(Pow(u_, m_), Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, d_, m_, Integer(1)), x_)), (x, u, m, c, d)))
    rule560 = ReplacementRule(pattern560, lambda x, u, m, c, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule560)

    pattern561 = Pattern(Int(Mul(Pow(u_, m_), Pow(Mul(d_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, m_, n_), x_)), (x, u, m, n, d)))
    rule561 = ReplacementRule(pattern561, lambda x, u, m, n, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule561)

    pattern562 = Pattern(Int(Mul(Pow(u_, m_), Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, n_), x_)), (x, u, m, n, c)))
    rule562 = ReplacementRule(pattern562, lambda x, u, m, n, c : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule562)

    pattern563 = Pattern(Int(Mul(d_, u_, Pow(Mul(b_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), d_, m_, Integer(1)), x_)), (x, u, m, d, b)))
    rule563 = ReplacementRule(pattern563, lambda x, u, m, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule563)

    pattern564 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, Integer(1), m_, Integer(1)), x_)), (x, u, m, c, b)))
    rule564 = ReplacementRule(pattern564, lambda x, u, m, c, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Add(c, x)), x), x, u)))
    rubi.add(rule564)

    pattern565 = Pattern(Int(Mul(Pow(u_, n_), Pow(Mul(b_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), m_, n_), x_)), (x, u, m, n, b)))
    rule565 = ReplacementRule(pattern565, lambda x, u, m, n, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule565)

    pattern566 = Pattern(Int(Mul(d_, u_, Pow(Add(a_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), d_, m_, Integer(1)), x_)), (x, u, m, a, d)))
    rule566 = ReplacementRule(pattern566, lambda x, u, m, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule566)

    pattern567 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, Integer(1), m_, Integer(1)), x_)), (x, u, m, c, a)))
    rule567 = ReplacementRule(pattern567, lambda x, u, m, c, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Add(c, x)), x), x, u)))
    rubi.add(rule567)

    pattern568 = Pattern(Int(Mul(Pow(u_, n_), Pow(Add(a_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, u, m, n, a)))
    rule568 = ReplacementRule(pattern568, lambda x, u, m, n, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule568)

    pattern569 = Pattern(Int(Mul(u_, Pow(Add(a_, Mul(b_, u_)), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), Integer(1), m_, Integer(1)), x_)), (x, u, m, a, b)))
    rule569 = ReplacementRule(pattern569, lambda x, u, m, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule569)

    pattern570 = Pattern(Int(Mul(u_, Add(c_, Mul(d_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, d_, Integer(1), Integer(1)), x_)), (x, u, c, d)))
    rule570 = ReplacementRule(pattern570, lambda x, u, c, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule570)

    pattern571 = Pattern(Int(Mul(u_, Pow(Mul(d_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, Integer(1), n_), x_)), (x, u, n, d)))
    rule571 = ReplacementRule(pattern571, lambda x, u, n, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule571)

    pattern572 = Pattern(Int(Mul(u_, Pow(Add(c_, u_), n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1), n_), x_)), (x, u, n, c)))
    rule572 = ReplacementRule(pattern572, lambda x, u, n, c : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule572)

    pattern573 = Pattern(Int(Mul(b_, d_, Pow(u_, Integer(2))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), d_, Integer(1), Integer(1)), x_)), (x, u, d, b)))
    rule573 = ReplacementRule(pattern573, lambda x, u, d, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, d, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule573)

    pattern574 = Pattern(Int(Mul(b_, u_, Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, c_, Integer(1), Integer(1), Integer(1)), x_)), (x, u, c, b)))
    rule574 = ReplacementRule(pattern574, lambda x, u, c, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Add(c, x)), x), x, u)))
    rubi.add(rule574)

    pattern575 = Pattern(Int(Mul(b_, u_, Pow(u_, n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), Integer(1), n_), x_)), (x, u, n, b)))
    rule575 = ReplacementRule(pattern575, lambda x, u, n, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(x, n)), x), x, u)))
    rubi.add(rule575)

    pattern576 = Pattern(Int(Mul(d_, u_, Add(a_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), d_, Integer(1), Integer(1)), x_)), (x, u, a, d)))
    rule576 = ReplacementRule(pattern576, lambda x, u, a, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Add(a, x)), x), x, u)))
    rubi.add(rule576)

    pattern577 = Pattern(Int(Mul(Add(a_, u_), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), c_, Integer(1), Integer(1), Integer(1)), x_)), (x, u, c, a)))
    rule577 = ReplacementRule(pattern577, lambda x, u, c, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Add(c, x)), x), x, u)))
    rubi.add(rule577)

    pattern578 = Pattern(Int(Mul(Pow(u_, n_), Add(a_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1), n_), x_)), (x, u, n, a)))
    rule578 = ReplacementRule(pattern578, lambda x, u, n, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Add(a, x)), x), x, u)))
    rubi.add(rule578)

    pattern579 = Pattern(Int(Mul(u_, Add(a_, Mul(b_, u_))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, b_, Integer(0), Integer(1), Integer(1), Integer(1)), x_)), (x, u, a, b)))
    rule579 = ReplacementRule(pattern579, lambda x, u, a, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule579)

    pattern580 = Pattern(Int(Mul(d_, u_, Pow(u_, m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, m_, Integer(1)), x_)), (x, u, m, d)))
    rule580 = ReplacementRule(pattern580, lambda x, u, m, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(x, m)), x), x, u)))
    rubi.add(rule580)

    pattern581 = Pattern(Int(Mul(Pow(u_, m_), Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), m_, Integer(1)), x_)), (x, u, m, c)))
    rule581 = ReplacementRule(pattern581, lambda x, u, m, c : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Add(c, x)), x), x, u)))
    rubi.add(rule581)

    pattern582 = Pattern(Int(Mul(Pow(u_, m_), Pow(u_, n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_), x_)), (x, u, m, n)))
    rule582 = ReplacementRule(pattern582, lambda x, u, m, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(x, n)), x), x, u)))
    rubi.add(rule582)

    pattern583 = Pattern(Int(Mul(u_, Pow(Mul(b_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), m_, Integer(1)), x_)), (x, u, m, b)))
    rule583 = ReplacementRule(pattern583, lambda x, u, m, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule583)

    pattern584 = Pattern(Int(Mul(u_, Pow(Add(a_, u_), m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), m_, Integer(1)), x_)), (x, u, m, a)))
    rule584 = ReplacementRule(pattern584, lambda x, u, m, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule584)

    pattern585 = Pattern(Int(Mul(d_, Pow(u_, Integer(2))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), d_, Integer(1), Integer(1)), x_)), (x, u, d)))
    rule585 = ReplacementRule(pattern585, lambda x, u, d : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule585)

    pattern586 = Pattern(Int(Mul(u_, Add(c_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), c_, Integer(1), Integer(1), Integer(1)), x_)), (x, u, c)))
    rule586 = ReplacementRule(pattern586, lambda x, u, c : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(c, x)), x), x, u)))
    rubi.add(rule586)

    pattern587 = Pattern(Int(Mul(u_, Pow(u_, n_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), Integer(1), n_), x_)), (x, u, n)))
    rule587 = ReplacementRule(pattern587, lambda x, u, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(x, n)), x), x, u)))
    rubi.add(rule587)

    pattern588 = Pattern(Int(Mul(b_, Pow(u_, Integer(2))), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), b_, Integer(0), Integer(1), Integer(1), Integer(1)), x_)), (x, u, b)))
    rule588 = ReplacementRule(pattern588, lambda x, u, b : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule588)

    pattern589 = Pattern(Int(Mul(u_, Add(a_, u_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(a_, Integer(1), Integer(0), Integer(1), Integer(1), Integer(1)), x_)), (x, u, a)))
    rule589 = ReplacementRule(pattern589, lambda x, u, a : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(a, x)), x), x, u)))
    rubi.add(rule589)

    pattern590 = Pattern(Int(Mul(u_, Pow(u_, m_)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), m_, Integer(1)), x_)), (x, u, m)))
    rule590 = ReplacementRule(pattern590, lambda x, u, m : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(x, m)), x), x, u)))
    rubi.add(rule590)

    pattern591 = Pattern(Int(Pow(u_, Integer(2)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0))), FreeQ(List(Integer(0), Integer(1), Integer(0), Integer(1), Integer(1), Integer(1)), x_)), (x, u)))
    rule591 = ReplacementRule(pattern591, lambda x, u : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, Integer(2)), x), x, u)))
    rubi.add(rule591)

    return rubi
