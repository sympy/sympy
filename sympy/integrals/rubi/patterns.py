from sympy.external import import_module
matchpy = import_module("matchpy")

if matchpy:
    Wildcard, Pattern, ReplacementRule, ManyToOneReplacer = matchpy.Wildcard, matchpy.Pattern, matchpy.ReplacementRule, matchpy.ManyToOneReplacer
else:
    raise ImportError('MatchPy could not be imported')

from sympy.integrals.rubi.operation import *
from sympy.integrals.rubi.symbol import VariableSymbol, Integer
from sympy.integrals.rubi.constraint import cons, FreeQ

a, b, c, d, e, f, g, h, x, u, p = map(VariableSymbol, 'abcdefghxup')
n, m = map(VariableSymbol, 'nm')
zoo = VariableSymbol('zoo')

a_, b_, c_, d_, e_, f_, g_, h_, p_ = map(Wildcard.dot, 'abcdefghp')
n_, m_ = map(Wildcard.dot, 'nm')
x_, u_ = map(Wildcard.symbol, 'xu')


def rubi_object():
    rubi = ManyToOneReplacer()
    pattern1 = Pattern(Int(Mul(Pow(Mul(c_, x_), m_), Pow(Add(a_, Mul(b_, Pow(x_, n_))), p_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (m, a, b, c, p, x, n)))
    rule1 = ReplacementRule(pattern1, lambda m, a, b, c, p, x, n : Mul(Pow(a, p), Pow(c, Integer(-1)), Pow(Mul(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Pow(n, Integer(-1)), Add(m, Integer(1))), Add(Integer(1), Mul(Pow(n, Integer(-1)), Add(m, Integer(1)))), Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(x, n)))))
    rubi.add(rule1)

    pattern2 = Pattern(Int(Mul(c_, x_, Pow(Add(a_, Mul(b_, Pow(x_, n_))), p_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (a, b, c, p, x, n)))
    rule2 = ReplacementRule(pattern2, lambda a, b, c, p, x, n : Mul(Integer(1/2), Pow(a, p), c, Pow(x, Integer(2)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Integer(2), Pow(n, Integer(-1))), Add(Integer(1), Mul(Integer(2), Pow(n, Integer(-1)))), Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(x, n)))))
    rubi.add(rule2)

    pattern3 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(a_, Mul(b_, Pow(x_, n_))), p_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (m, a, b, p, x, n)))
    rule3 = ReplacementRule(pattern3, lambda m, a, b, p, x, n : Mul(Pow(a, p), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Pow(n, Integer(-1)), Add(m, Integer(1))), Add(Integer(1), Mul(Pow(n, Integer(-1)), Add(m, Integer(1)))), Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(x, n)))))
    rubi.add(rule3)

    pattern4 = Pattern(Int(Mul(Pow(Mul(c_, x_), m_), Pow(Add(a_, Pow(x_, n_)), p_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (m, a, c, p, x, n)))
    rule4 = ReplacementRule(pattern4, lambda m, a, c, p, x, n : Mul(Pow(a, p), Pow(c, Integer(-1)), Pow(Mul(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Pow(n, Integer(-1)), Add(m, Integer(1))), Add(Integer(1), Mul(Pow(n, Integer(-1)), Add(m, Integer(1)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, n)))))
    rubi.add(rule4)

    pattern5 = Pattern(Int(Mul(x_, Pow(Add(a_, Mul(b_, Pow(x_, n_))), p_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (a, b, p, x, n)))
    rule5 = ReplacementRule(pattern5, lambda a, b, p, x, n : Mul(Integer(1/2), Pow(a, p), Pow(x, Integer(2)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Integer(2), Pow(n, Integer(-1))), Add(Integer(1), Mul(Integer(2), Pow(n, Integer(-1)))), Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(x, n)))))
    rubi.add(rule5)

    pattern6 = Pattern(Int(Mul(c_, x_, Pow(Add(a_, Pow(x_, n_)), p_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (a, c, p, x, n)))
    rule6 = ReplacementRule(pattern6, lambda a, c, p, x, n : Mul(Integer(1/2), Pow(a, p), c, Pow(x, Integer(2)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Integer(2), Pow(n, Integer(-1))), Add(Integer(1), Mul(Integer(2), Pow(n, Integer(-1)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, n)))))
    rubi.add(rule6)

    pattern7 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(a_, Pow(x_, n_)), p_)), x_), FreeQ(a, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (m, a, p, x, n)))
    rule7 = ReplacementRule(pattern7, lambda m, a, p, x, n : Mul(Pow(a, p), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Pow(n, Integer(-1)), Add(m, Integer(1))), Add(Integer(1), Mul(Pow(n, Integer(-1)), Add(m, Integer(1)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, n)))))
    rubi.add(rule7)

    pattern8 = Pattern(Int(Mul(x_, Pow(Add(a_, Pow(x_, n_)), p_)), x_), FreeQ(a, x), FreeQ(n, x), FreeQ(p, x), cons(And(Not(PositiveIntegerQ(p_)), Or(NegativeIntegerQ(p_), PositiveQ(a_))), (a, p, x, n)))
    rule8 = ReplacementRule(pattern8, lambda a, p, x, n : Mul(Integer(1/2), Pow(a, p), Pow(x, Integer(2)), Hypergeometric2F1(Mul(Integer(-1), p), Mul(Integer(2), Pow(n, Integer(-1))), Add(Integer(1), Mul(Integer(2), Pow(n, Integer(-1)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, n)))))
    rubi.add(rule8)

    pattern9 = Pattern(Int(a, x_), FreeQ(a, x), cons(True, (x,)))
    rule9 = ReplacementRule(pattern9, lambda x : Mul(a, x))
    rubi.add(rule9)

    pattern10 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(True, (x,)))
    rule10 = ReplacementRule(pattern10, lambda x : Log(x))
    rubi.add(rule10)

    pattern11 = Pattern(Int(Pow(x_, m_), x_), FreeQ(m, x), cons(NonzeroQ(Add(m_, Integer(1))), (m, x)))
    rule11 = ReplacementRule(pattern11, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule11)

    pattern12 = Pattern(Int(x_, x_), cons(NonzeroQ(Integer(2)), (x,)))
    rule12 = ReplacementRule(pattern12, lambda x : Mul(Integer(1/2), Pow(x, Integer(2))))
    rubi.add(rule12)

    pattern13 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), x_), FreeQ(a, x), FreeQ(b, x), cons(True, (a, x, b)))
    rule13 = ReplacementRule(pattern13, lambda a, x, b : Mul(Pow(b, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))))
    rubi.add(rule13)

    pattern14 = Pattern(Int(Pow(Add(a_, x_), Integer(-1)), x_), FreeQ(a, x), cons(True, (a, x)))
    rule14 = ReplacementRule(pattern14, lambda a, x : Log(RemoveContent(Add(a, x), x)))
    rubi.add(rule14)

    pattern15 = Pattern(Int(Pow(Add(a_, Mul(b_, x_)), m_), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), cons(NonzeroQ(Add(m_, Integer(1))), (m, b, a, x)))
    rule15 = ReplacementRule(pattern15, lambda m, b, a, x : Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule15)

    pattern16 = Pattern(Int(Pow(Mul(b_, x_), m_), x_), FreeQ(b, x), FreeQ(m, x), cons(NonzeroQ(Add(m_, Integer(1))), (m, b, x)))
    rule16 = ReplacementRule(pattern16, lambda m, b, x : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule16)

    pattern17 = Pattern(Int(Pow(Add(a_, x_), m_), x_), FreeQ(a, x), FreeQ(m, x), cons(NonzeroQ(Add(m_, Integer(1))), (m, a, x)))
    rule17 = ReplacementRule(pattern17, lambda m, a, x : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule17)

    pattern18 = Pattern(Int(Pow(x_, m_), x_), FreeQ(m, x), cons(NonzeroQ(Add(m_, Integer(1))), (m, x)))
    rule18 = ReplacementRule(pattern18, lambda m, x : Mul(Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule18)

    pattern19 = Pattern(Int(Pow(Add(a_, Mul(b_, u_)), m_), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (m, b, a, u, x)))
    rule19 = ReplacementRule(pattern19, lambda m, b, a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, Mul(b, x)), m), x), x, u)))
    rubi.add(rule19)

    pattern20 = Pattern(Int(Pow(Mul(b_, u_), m_), x_), FreeQ(b, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (m, b, u, x)))
    rule20 = ReplacementRule(pattern20, lambda m, b, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Mul(b, x), m), x), x, u)))
    rubi.add(rule20)

    pattern21 = Pattern(Int(Pow(Add(a_, u_), m_), x_), FreeQ(a, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (m, a, u, x)))
    rule21 = ReplacementRule(pattern21, lambda m, a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(Add(a, x), m), x), x, u)))
    rubi.add(rule21)

    pattern22 = Pattern(Int(Pow(u_, m_), x_), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Add(u_, Mul(Integer(-1), x_)))), (m, u, x)))
    rule22 = ReplacementRule(pattern22, lambda m, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, m), x), x, u)))
    rubi.add(rule22)

    pattern23 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), (a, b, c, d, x)))
    rule23 = ReplacementRule(pattern23, lambda a, b, c, d, x : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule23)

    pattern24 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), c_)), (a, c, d, x)))
    rule24 = ReplacementRule(pattern24, lambda a, c, d, x : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule24)

    pattern25 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(ZeroQ(Add(a_, Mul(b_, c_))), (a, b, c, x)))
    rule25 = ReplacementRule(pattern25, lambda a, b, c, x : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Integer(-1)), x))
    rubi.add(rule25)

    pattern26 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(c, x), cons(ZeroQ(Add(a_, c_)), (a, c, x)))
    rule26 = ReplacementRule(pattern26, lambda a, c, x : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), Integer(-1)), x))
    rubi.add(rule26)

    pattern27 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), (b, a, c, d, x)))
    rule27 = ReplacementRule(pattern27, lambda b, a, c, d, x : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule27)

    pattern28 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Mul(b_, c_)), (b, c, d, x)))
    rule28 = ReplacementRule(pattern28, lambda b, c, d, x : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x))))
    rubi.add(rule28)

    pattern29 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule29 = ReplacementRule(pattern29, lambda a, c, d, x : Add(Mul(Integer(-1), d, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule29)

    pattern30 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(NonzeroQ(Mul(Integer(-1), a_, d_)), (b, a, d, x)))
    rule30 = ReplacementRule(pattern30, lambda b, a, d, x : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule30)

    pattern31 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), (b, a, c, x)))
    rule31 = ReplacementRule(pattern31, lambda b, a, c, x : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule31)

    pattern32 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1))), x_), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(c_), (c, d, x)))
    rule32 = ReplacementRule(pattern32, lambda c, d, x : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Int(Pow(Add(c, Mul(d, x)), Integer(-1)), x)), Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule32)

    pattern33 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), FreeQ(b, x), FreeQ(d, x), cons(NonzeroQ(Integer(0)), (b, d, x)))
    rule33 = ReplacementRule(pattern33, lambda b, d, x : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x))))
    rubi.add(rule33)

    pattern34 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(b, x), FreeQ(c, x), cons(NonzeroQ(Mul(b_, c_)), (b, c, x)))
    rule34 = ReplacementRule(pattern34, lambda b, c, x : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule34)

    pattern35 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(NonzeroQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule35 = ReplacementRule(pattern35, lambda a, d, x : Add(Mul(Pow(a, Integer(-1)), Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule35)

    pattern36 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(c, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule36 = ReplacementRule(pattern36, lambda a, c, x : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule36)

    pattern37 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(NonzeroQ(Mul(Integer(-1), a_)), (b, a, x)))
    rule37 = ReplacementRule(pattern37, lambda b, a, x : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Pow(Add(a, Mul(b, x)), Integer(-1)), x)), Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule37)

    pattern38 = Pattern(Int(Mul(Pow(d_, Integer(-1)), Pow(x_, Integer(-2))), x_), FreeQ(d, x), cons(NonzeroQ(Integer(0)), (d, x)))
    rule38 = ReplacementRule(pattern38, lambda d, x : Add(Mul(zoo, d, Int(Mul(Pow(d, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule38)

    pattern39 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1))), x_), FreeQ(c, x), cons(NonzeroQ(c_), (c, x)))
    rule39 = ReplacementRule(pattern39, lambda c, x : Add(Mul(Pow(c, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(c, Integer(-1)), Int(Pow(Add(c, x), Integer(-1)), x))))
    rubi.add(rule39)

    pattern40 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-2))), x_), FreeQ(b, x), cons(NonzeroQ(Integer(0)), (b, x)))
    rule40 = ReplacementRule(pattern40, lambda b, x : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1))), x)), Mul(zoo, Pow(b, Integer(-1)), Int(Pow(x, Integer(-1)), x))))
    rubi.add(rule40)

    pattern41 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(NonzeroQ(Mul(Integer(-1), a_)), (a, x)))
    rule41 = ReplacementRule(pattern41, lambda a, x : Add(Mul(Pow(a, Integer(-1)), Int(Pow(x, Integer(-1)), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Pow(Add(a, x), Integer(-1)), x))))
    rubi.add(rule41)

    pattern42 = Pattern(Int(Pow(x_, Integer(-2)), x_), cons(NonzeroQ(Integer(0)), (x,)))
    rule42 = ReplacementRule(pattern42, lambda x : Mul(Integer(2), zoo, Int(Pow(x, Integer(-1)), x)))
    rubi.add(rule42)

    pattern43 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, a, c, d, x, n)))
    rule43 = ReplacementRule(pattern43, lambda m, b, a, c, d, x, n : Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule43)

    pattern44 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x, n)))
    rule44 = ReplacementRule(pattern44, lambda b, a, c, d, x, n : Mul(Integer(1/2), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))
    rubi.add(rule44)

    pattern45 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2)))), (m, b, c, d, x, n)))
    rule45 = ReplacementRule(pattern45, lambda m, b, c, d, x, n : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule45)

    pattern46 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (m, a, c, d, x, n)))
    rule46 = ReplacementRule(pattern46, lambda m, a, c, d, x, n : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule46)

    pattern47 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (m, b, a, c, x, n)))
    rule47 = ReplacementRule(pattern47, lambda m, b, a, c, x, n : Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule47)

    pattern48 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3)))), (b, c, d, x, n)))
    rule48 = ReplacementRule(pattern48, lambda b, c, d, x, n : Mul(Integer(1/2), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule48)

    pattern49 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x, n)))
    rule49 = ReplacementRule(pattern49, lambda a, c, d, x, n : Mul(Integer(1/2), Pow(Add(a, x), Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))
    rubi.add(rule49)

    pattern50 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x, n)))
    rule50 = ReplacementRule(pattern50, lambda b, a, c, x, n : Mul(Integer(1/2), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule50)

    pattern51 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2)))), (m, c, d, x, n)))
    rule51 = ReplacementRule(pattern51, lambda m, c, d, x, n : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule51)

    pattern52 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2)))), (m, b, c, x, n)))
    rule52 = ReplacementRule(pattern52, lambda m, b, c, x, n : Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule52)

    pattern53 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2))), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (m, a, c, x, n)))
    rule53 = ReplacementRule(pattern53, lambda m, a, c, x, n : Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule53)

    pattern54 = Pattern(Int(Mul(x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3)))), (c, d, x, n)))
    rule54 = ReplacementRule(pattern54, lambda c, d, x, n : Mul(Integer(1/2), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))))
    rubi.add(rule54)

    pattern55 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), NonzeroQ(Mul(b_, c_)), ZeroQ(Add(n_, Integer(3)))), (b, c, x, n)))
    rule55 = ReplacementRule(pattern55, lambda b, c, x, n : Mul(Integer(1/2), b, Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule55)

    pattern56 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), ZeroQ(Add(n_, Integer(3))), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x, n)))
    rule56 = ReplacementRule(pattern56, lambda a, c, x, n : Mul(Integer(1/2), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule56)

    pattern57 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), ZeroQ(Add(m_, n_, Integer(2)))), (m, c, x, n)))
    rule57 = ReplacementRule(pattern57, lambda m, c, x, n : Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))))
    rubi.add(rule57)

    pattern58 = Pattern(Int(Mul(x_, Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(2)), NonzeroQ(c_), ZeroQ(Add(n_, Integer(3)))), (c, x, n)))
    rule58 = ReplacementRule(pattern58, lambda c, x, n : Mul(Integer(1/2), Pow(c, Integer(-1)), Pow(x, Integer(2)), Pow(Add(c, x), Add(n, Integer(1)))))
    rubi.add(rule58)

    pattern59 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x)))
    rule59 = ReplacementRule(pattern59, lambda m, a, b, c, d, x : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule59)

    pattern60 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(Mul(a_, d_), c_))), (m, a, c, d, x)))
    rule60 = ReplacementRule(pattern60, lambda m, a, c, d, x : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule60)

    pattern61 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(a_, Mul(b_, c_)))), (m, a, b, c, x)))
    rule61 = ReplacementRule(pattern61, lambda m, a, b, c, x : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule61)

    pattern62 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PositiveIntegerQ(Add(m_, Integer(1/2))), ZeroQ(Add(a_, c_))), (m, a, c, x)))
    rule62 = ReplacementRule(pattern62, lambda m, a, c, x : Add(Mul(Integer(2), a, c, m, Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(-1))), Pow(Add(c, x), Add(m, Integer(-1)))), x)), Mul(x, Pow(Add(a, x), m), Pow(Add(c, x), m), Pow(Add(Mul(Integer(2), m), Integer(1)), Integer(-1)))))
    rubi.add(rule62)

    pattern63 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-3/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-3/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), (a, b, c, d, x)))
    rule63 = ReplacementRule(pattern63, lambda a, b, c, d, x : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule63)

    pattern64 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-3/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-3/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), c_)), (a, c, d, x)))
    rule64 = ReplacementRule(pattern64, lambda a, c, d, x : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))
    rubi.add(rule64)

    pattern65 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-3/2)), Pow(Add(c_, x_), Integer(-3/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(ZeroQ(Add(a_, Mul(b_, c_))), (a, b, c, x)))
    rule65 = ReplacementRule(pattern65, lambda a, b, c, x : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, Mul(b, x))), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule65)

    pattern66 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-3/2)), Pow(Add(c_, x_), Integer(-3/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(ZeroQ(Add(a_, c_)), (a, c, x)))
    rule66 = ReplacementRule(pattern66, lambda a, c, x : Mul(Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Sqrt(Add(a, x)), Integer(-1)), Pow(Sqrt(Add(c, x)), Integer(-1))))
    rubi.add(rule66)

    pattern67 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x)))
    rule67 = ReplacementRule(pattern67, lambda m, a, b, c, d, x : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule67)

    pattern68 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(Mul(a_, d_), c_))), (m, a, c, d, x)))
    rule68 = ReplacementRule(pattern68, lambda m, a, c, d, x : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule68)

    pattern69 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(a_, Mul(b_, c_)))), (m, a, b, c, x)))
    rule69 = ReplacementRule(pattern69, lambda m, a, b, c, x : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule69)

    pattern70 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(NegativeIntegerQ(Add(m_, Integer(3/2))), ZeroQ(Add(a_, c_))), (m, a, c, x)))
    rule70 = ReplacementRule(pattern70, lambda m, a, c, x : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), x, Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(Mul(Integer(2), m), Integer(3)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(m, Integer(1)))), x))))
    rubi.add(rule70)

    pattern71 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, a, b, c, d, x)))
    rule71 = ReplacementRule(pattern71, lambda m, a, b, c, d, x : Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule71)

    pattern72 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (a, b, c, d, x)))
    rule72 = ReplacementRule(pattern72, lambda a, b, c, d, x : Int(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), x))
    rubi.add(rule72)

    pattern73 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, a, c, d, x)))
    rule73 = ReplacementRule(pattern73, lambda m, a, c, d, x : Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x))
    rubi.add(rule73)

    pattern74 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, a, b, c, x)))
    rule74 = ReplacementRule(pattern74, lambda m, a, b, c, x : Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x))
    rubi.add(rule74)

    pattern75 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (a, c, d, x)))
    rule75 = ReplacementRule(pattern75, lambda a, c, d, x : Int(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), x))
    rubi.add(rule75)

    pattern76 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (a, b, c, x)))
    rule76 = ReplacementRule(pattern76, lambda a, b, c, x : Int(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), x))
    rubi.add(rule76)

    pattern77 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), cons(And(ZeroQ(Add(a_, c_)), Or(IntegerQ(m_), And(PositiveQ(a_), PositiveQ(c_)))), (m, a, c, x)))
    rule77 = ReplacementRule(pattern77, lambda m, a, c, x : Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x))
    rubi.add(rule77)

    pattern78 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, c_)), Or(IntegerQ(Integer(1)), And(PositiveQ(a_), PositiveQ(c_)))), (a, c, x)))
    rule78 = ReplacementRule(pattern78, lambda a, c, x : Int(Add(Mul(a, c), Pow(x, Integer(2))), x))
    rubi.add(rule78)

    pattern79 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (a, b, c, d, x)))
    rule79 = ReplacementRule(pattern79, lambda a, b, c, d, x : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule79)

    pattern80 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(Mul(a_, d_), c_))), (a, c, d, x)))
    rule80 = ReplacementRule(pattern80, lambda a, c, d, x : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule80)

    pattern81 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_)), ZeroQ(Add(a_, Mul(b_, c_)))), (a, b, c, x)))
    rule81 = ReplacementRule(pattern81, lambda a, b, c, x : Mul(Pow(b, Integer(-1)), ArcCosh(Mul(Pow(a, Integer(-1)), b, x))))
    rubi.add(rule81)

    pattern82 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PositiveQ(a_), ZeroQ(Add(a_, c_))), (a, c, x)))
    rule82 = ReplacementRule(pattern82, lambda a, c, x : ArcCosh(Mul(Pow(a, Integer(-1)), x)))
    rubi.add(rule82)

    pattern83 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), (a, b, c, d, x)))
    rule83 = ReplacementRule(pattern83, lambda a, b, c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule83)

    pattern84 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(ZeroQ(Add(Mul(a_, d_), c_)), (a, c, d, x)))
    rule84 = ReplacementRule(pattern84, lambda a, c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule84)

    pattern85 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(ZeroQ(Add(a_, Mul(b_, c_))), (a, b, c, x)))
    rule85 = ReplacementRule(pattern85, lambda a, b, c, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule85)

    pattern86 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(ZeroQ(Add(a_, c_)), (a, c, x)))
    rule86 = ReplacementRule(pattern86, lambda a, c, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule86)

    pattern87 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x)))
    rule87 = ReplacementRule(pattern87, lambda m, a, b, c, d, x : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule87)

    pattern88 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(Mul(a_, d_), c_))), (m, a, c, d, x)))
    rule88 = ReplacementRule(pattern88, lambda m, a, c, d, x : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, Mul(d, x)), FracPart(m)), Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule88)

    pattern89 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), cons(And(Not(IntegerQ(Mul(Integer(2), m_))), ZeroQ(Add(a_, Mul(b_, c_)))), (m, a, b, c, x)))
    rule89 = ReplacementRule(pattern89, lambda m, a, b, c, x : Mul(Pow(Add(a, Mul(b, x)), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule89)

    pattern90 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), cons(And(ZeroQ(Add(a_, c_)), Not(IntegerQ(Mul(Integer(2), m_)))), (m, a, c, x)))
    rule90 = ReplacementRule(pattern90, lambda m, a, c, x : Mul(Pow(Add(a, x), FracPart(m)), Pow(Add(c, x), FracPart(m)), Pow(Add(Mul(a, c), Pow(x, Integer(2))), Mul(Integer(-1), FracPart(m))), Int(Pow(Add(Mul(a, c), Pow(x, Integer(2))), m), x)))
    rubi.add(rule90)

    pattern91 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-5/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (a, b, c, d, x)))
    rule91 = ReplacementRule(pattern91, lambda a, b, c, d, x : Add(Mul(Integer(1/2), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-1/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule91)

    pattern92 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-5/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (a, c, d, x)))
    rule92 = ReplacementRule(pattern92, lambda a, c, d, x : Add(Mul(Integer(1/2), Add(Mul(Integer(-1), a, d), c), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Integer(-1/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule92)

    pattern93 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-5/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (a, b, c, x)))
    rule93 = ReplacementRule(pattern93, lambda a, b, c, x : Add(Mul(Integer(1/2), Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-1/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule93)

    pattern94 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-5/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, c_)), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (a, c, x)))
    rule94 = ReplacementRule(pattern94, lambda a, c, x : Add(Mul(Integer(1/2), Add(Mul(Integer(-1), a), c), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(2), Pow(Add(a, x), Integer(-1/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule94)

    pattern95 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-9/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1)), d_))), (a, b, c, d, x)))
    rule95 = ReplacementRule(pattern95, lambda a, b, c, d, x : Add(Mul(Integer(-1), Integer(1/5), Pow(b, Integer(-1)), d, Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule95)

    pattern96 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-9/4)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(a_, d_), c_)), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1)), d_))), (a, c, d, x)))
    rule96 = ReplacementRule(pattern96, lambda a, c, d, x : Add(Mul(Integer(-1), Integer(1/5), d, Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, Mul(d, x)), Integer(-1/4)))))
    rubi.add(rule96)

    pattern97 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-9/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, Mul(b_, c_))), PosQ(Mul(Pow(a_, Integer(-1)), b_, Pow(c_, Integer(-1))))), (a, b, c, x)))
    rule97 = ReplacementRule(pattern97, lambda a, b, c, x : Add(Mul(Integer(-1), Integer(1/5), Pow(b, Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Integer(-5/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule97)

    pattern98 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-9/4)), Pow(Add(c_, x_), Integer(-1/4))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(ZeroQ(Add(a_, c_)), PosQ(Mul(Pow(a_, Integer(-1)), Pow(c_, Integer(-1))))), (a, c, x)))
    rule98 = ReplacementRule(pattern98, lambda a, c, x : Add(Mul(Integer(-1), Integer(1/5), Int(Mul(Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-5/4))), x)), Mul(Integer(-1), Integer(4/5), Pow(Add(a, x), Integer(-5/4)), Pow(Add(c, x), Integer(-1/4)))))
    rubi.add(rule98)

    pattern99 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x, n)))
    rule99 = ReplacementRule(pattern99, lambda m, a, b, c, d, x, n : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule99)

    pattern100 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(Mul(a_, d_), c_))), (m, a, c, d, x, n)))
    rule100 = ReplacementRule(pattern100, lambda m, a, c, d, x, n : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule100)

    pattern101 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, Mul(b_, c_)))), (m, a, b, c, x, n)))
    rule101 = ReplacementRule(pattern101, lambda m, a, b, c, x, n : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule101)

    pattern102 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(Integer(0), m_, n_), ZeroQ(Add(a_, c_))), (m, a, c, x, n)))
    rule102 = ReplacementRule(pattern102, lambda m, a, c, x, n : Add(Mul(Integer(2), c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule102)

    pattern103 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(Mul(a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x, n)))
    rule103 = ReplacementRule(pattern103, lambda m, a, b, c, d, x, n : Add(Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule103)

    pattern104 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(Mul(a_, d_), c_))), (m, a, c, d, x, n)))
    rule104 = ReplacementRule(pattern104, lambda m, a, c, d, x, n : Add(Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule104)

    pattern105 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, Mul(b_, c_)))), (m, a, b, c, x, n)))
    rule105 = ReplacementRule(pattern105, lambda m, a, b, c, x, n : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule105)

    pattern106 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(IntegerQ(Add(m_, Integer(1/2))), IntegerQ(Add(n_, Integer(1/2))), Less(m_, n_, Integer(0)), ZeroQ(Add(a_, c_))), (m, a, c, x, n)))
    rule106 = ReplacementRule(pattern106, lambda m, a, c, x, n : Add(Mul(Integer(-1), Integer(1/2), Pow(a, Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(1/2), Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule106)

    pattern107 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, a, c, d, x, n)))
    rule107 = ReplacementRule(pattern107, lambda m, b, a, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule107)

    pattern108 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, a, c, d, x, n)))
    rule108 = ReplacementRule(pattern108, lambda b, a, c, d, x, n : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule108)

    pattern109 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, c, d, x, n)))
    rule109 = ReplacementRule(pattern109, lambda m, b, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule109)

    pattern110 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, a, c, d, x, n)))
    rule110 = ReplacementRule(pattern110, lambda m, a, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule110)

    pattern111 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, a, c, d, x)))
    rule111 = ReplacementRule(pattern111, lambda m, b, a, c, d, x : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule111)

    pattern112 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, a, d, x, n)))
    rule112 = ReplacementRule(pattern112, lambda m, b, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule112)

    pattern113 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, a, c, x, n)))
    rule113 = ReplacementRule(pattern113, lambda m, b, a, c, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule113)

    pattern114 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, c, d, x, n)))
    rule114 = ReplacementRule(pattern114, lambda b, c, d, x, n : Int(ExpandIntegrand(Mul(b, x, Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule114)

    pattern115 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (a, c, d, x, n)))
    rule115 = ReplacementRule(pattern115, lambda a, c, d, x, n : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule115)

    pattern116 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (b, a, c, d, x)))
    rule116 = ReplacementRule(pattern116, lambda b, a, c, d, x : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x))
    rubi.add(rule116)

    pattern117 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, a, d, x, n)))
    rule117 = ReplacementRule(pattern117, lambda b, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule117)

    pattern118 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, a, c, x, n)))
    rule118 = ReplacementRule(pattern118, lambda b, a, c, x, n : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x))
    rubi.add(rule118)

    pattern119 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, c, d, x, n)))
    rule119 = ReplacementRule(pattern119, lambda m, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule119)

    pattern120 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Add(c_, Mul(d_, x_))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, c, d, x)))
    rule120 = ReplacementRule(pattern120, lambda m, b, c, d, x : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule120)

    pattern121 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, d, x, n)))
    rule121 = ReplacementRule(pattern121, lambda m, b, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Mul(d, x), n)), x), x))
    rubi.add(rule121)

    pattern122 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, c, x, n)))
    rule122 = ReplacementRule(pattern122, lambda m, b, c, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule122)

    pattern123 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, a, c, d, x)))
    rule123 = ReplacementRule(pattern123, lambda m, a, c, d, x : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule123)

    pattern124 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, a, d, x, n)))
    rule124 = ReplacementRule(pattern124, lambda m, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule124)

    pattern125 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, a, c, x, n)))
    rule125 = ReplacementRule(pattern125, lambda m, a, c, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule125)

    pattern126 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, a, d, x)))
    rule126 = ReplacementRule(pattern126, lambda m, b, a, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule126)

    pattern127 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, a, c, x)))
    rule127 = ReplacementRule(pattern127, lambda m, b, a, c, x : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x))
    rubi.add(rule127)

    pattern128 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, a, x, n)))
    rule128 = ReplacementRule(pattern128, lambda m, b, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule128)

    pattern129 = Pattern(Int(Mul(x_, Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (c, d, x, n)))
    rule129 = ReplacementRule(pattern129, lambda c, d, x, n : Int(ExpandIntegrand(Mul(x, Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule129)

    pattern130 = Pattern(Int(Mul(b_, x_, Add(c_, Mul(d_, x_))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (b, c, d, x)))
    rule130 = ReplacementRule(pattern130, lambda b, c, d, x : Int(ExpandIntegrand(Mul(b, x, Add(c, Mul(d, x))), x), x))
    rubi.add(rule130)

    pattern131 = Pattern(Int(Mul(b_, x_, Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, d, x, n)))
    rule131 = ReplacementRule(pattern131, lambda b, d, x, n : Int(ExpandIntegrand(Mul(b, x, Pow(Mul(d, x), n)), x), x))
    rubi.add(rule131)

    pattern132 = Pattern(Int(Mul(b_, x_, Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, c, x, n)))
    rule132 = ReplacementRule(pattern132, lambda b, c, x, n : Int(ExpandIntegrand(Mul(b, x, Pow(Add(c, x), n)), x), x))
    rubi.add(rule132)

    pattern133 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (a, c, d, x)))
    rule133 = ReplacementRule(pattern133, lambda a, c, d, x : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, Mul(d, x))), x), x))
    rubi.add(rule133)

    pattern134 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, x_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (a, d, x, n)))
    rule134 = ReplacementRule(pattern134, lambda a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x))
    rubi.add(rule134)

    pattern135 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (a, c, x, n)))
    rule135 = ReplacementRule(pattern135, lambda a, c, x, n : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, x), n)), x), x))
    rubi.add(rule135)

    pattern136 = Pattern(Int(Mul(d_, x_, Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (b, a, d, x)))
    rule136 = ReplacementRule(pattern136, lambda b, a, d, x : Int(ExpandIntegrand(Mul(d, x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule136)

    pattern137 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (b, a, c, x)))
    rule137 = ReplacementRule(pattern137, lambda b, a, c, x : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x))
    rubi.add(rule137)

    pattern138 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, a, x, n)))
    rule138 = ReplacementRule(pattern138, lambda b, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule138)

    pattern139 = Pattern(Int(Mul(Pow(x_, m_), Add(c_, Mul(d_, x_))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, c, d, x)))
    rule139 = ReplacementRule(pattern139, lambda m, c, d, x : Int(ExpandIntegrand(Mul(Pow(x, m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule139)

    pattern140 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, d, x, n)))
    rule140 = ReplacementRule(pattern140, lambda m, d, x, n : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Mul(d, x), n)), x), x))
    rubi.add(rule140)

    pattern141 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, c, x, n)))
    rule141 = ReplacementRule(pattern141, lambda m, c, x, n : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule141)

    pattern142 = Pattern(Int(Mul(d_, x_, Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, d, x)))
    rule142 = ReplacementRule(pattern142, lambda m, b, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(Mul(b, x), m)), x), x))
    rubi.add(rule142)

    pattern143 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Add(c_, x_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(b_, c_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, c, x)))
    rule143 = ReplacementRule(pattern143, lambda m, b, c, x : Int(ExpandIntegrand(Mul(Pow(Mul(b, x), m), Add(c, x)), x), x))
    rubi.add(rule143)

    pattern144 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, b, x, n)))
    rule144 = ReplacementRule(pattern144, lambda m, b, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Mul(b, x), m)), x), x))
    rubi.add(rule144)

    pattern145 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, a, d, x)))
    rule145 = ReplacementRule(pattern145, lambda m, a, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule145)

    pattern146 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, a, c, x)))
    rule146 = ReplacementRule(pattern146, lambda m, a, c, x : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, x)), x), x))
    rubi.add(rule146)

    pattern147 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(n, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, a, x, n)))
    rule147 = ReplacementRule(pattern147, lambda m, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule147)

    pattern148 = Pattern(Int(Mul(x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, a, x)))
    rule148 = ReplacementRule(pattern148, lambda m, b, a, x : Int(ExpandIntegrand(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule148)

    pattern149 = Pattern(Int(Mul(x_, Add(c_, Mul(d_, x_))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (c, d, x)))
    rule149 = ReplacementRule(pattern149, lambda c, d, x : Int(ExpandIntegrand(Mul(x, Add(c, Mul(d, x))), x), x))
    rubi.add(rule149)

    pattern150 = Pattern(Int(Mul(x_, Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (d, x, n)))
    rule150 = ReplacementRule(pattern150, lambda d, x, n : Int(ExpandIntegrand(Mul(x, Pow(Mul(d, x), n)), x), x))
    rubi.add(rule150)

    pattern151 = Pattern(Int(Mul(x_, Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(n, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (c, x, n)))
    rule151 = ReplacementRule(pattern151, lambda c, x, n : Int(ExpandIntegrand(Mul(x, Pow(Add(c, x), n)), x), x))
    rubi.add(rule151)

    pattern152 = Pattern(Int(Mul(b_, d_, Pow(x_, Integer(2))), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (b, d, x)))
    rule152 = ReplacementRule(pattern152, lambda b, d, x : Int(ExpandIntegrand(Mul(b, d, Pow(x, Integer(2))), x), x))
    rubi.add(rule152)

    pattern153 = Pattern(Int(Mul(b_, x_, Add(c_, x_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(b_, c_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (b, c, x)))
    rule153 = ReplacementRule(pattern153, lambda b, c, x : Int(ExpandIntegrand(Mul(b, x, Add(c, x)), x), x))
    rubi.add(rule153)

    pattern154 = Pattern(Int(Mul(b_, x_, Pow(x_, n_)), x_), FreeQ(b, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (b, x, n)))
    rule154 = ReplacementRule(pattern154, lambda b, x, n : Int(ExpandIntegrand(Mul(b, x, Pow(x, n)), x), x))
    rubi.add(rule154)

    pattern155 = Pattern(Int(Mul(d_, x_, Add(a_, x_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (a, d, x)))
    rule155 = ReplacementRule(pattern155, lambda a, d, x : Int(ExpandIntegrand(Mul(d, x, Add(a, x)), x), x))
    rubi.add(rule155)

    pattern156 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (a, c, x)))
    rule156 = ReplacementRule(pattern156, lambda a, c, x : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, x)), x), x))
    rubi.add(rule156)

    pattern157 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, x_)), x_), FreeQ(a, x), FreeQ(n, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (a, x, n)))
    rule157 = ReplacementRule(pattern157, lambda a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, x)), x), x))
    rubi.add(rule157)

    pattern158 = Pattern(Int(Mul(x_, Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (b, a, x)))
    rule158 = ReplacementRule(pattern158, lambda b, a, x : Int(ExpandIntegrand(Mul(x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule158)

    pattern159 = Pattern(Int(Mul(d_, x_, Pow(x_, m_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, d, x)))
    rule159 = ReplacementRule(pattern159, lambda m, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(x, m)), x), x))
    rubi.add(rule159)

    pattern160 = Pattern(Int(Mul(Pow(x_, m_), Add(c_, x_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(c_), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, c, x)))
    rule160 = ReplacementRule(pattern160, lambda m, c, x : Int(ExpandIntegrand(Mul(Pow(x, m), Add(c, x)), x), x))
    rubi.add(rule160)

    pattern161 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Greater(Add(m_, n_, Integer(2)), Integer(0)), Less(Add(Mul(Integer(9), m_), Mul(Integer(5), n_), Integer(5)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Mul(Integer(4), n_)), Integer(0))))), (m, x, n)))
    rule161 = ReplacementRule(pattern161, lambda m, x, n : Int(ExpandIntegrand(Mul(Pow(x, m), Pow(x, n)), x), x))
    rubi.add(rule161)

    pattern162 = Pattern(Int(Mul(x_, Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, b, x)))
    rule162 = ReplacementRule(pattern162, lambda m, b, x : Int(ExpandIntegrand(Mul(x, Pow(Mul(b, x), m)), x), x))
    rubi.add(rule162)

    pattern163 = Pattern(Int(Mul(x_, Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(PositiveIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, a, x)))
    rule163 = ReplacementRule(pattern163, lambda m, a, x : Int(ExpandIntegrand(Mul(x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule163)

    pattern164 = Pattern(Int(Mul(d_, Pow(x_, Integer(2))), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (d, x)))
    rule164 = ReplacementRule(pattern164, lambda d, x : Int(ExpandIntegrand(Mul(d, Pow(x, Integer(2))), x), x))
    rubi.add(rule164)

    pattern165 = Pattern(Int(Mul(x_, Add(c_, x_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(c_), LessEqual(Integer(11), Integer(0))))), (c, x)))
    rule165 = ReplacementRule(pattern165, lambda c, x : Int(ExpandIntegrand(Mul(x, Add(c, x)), x), x))
    rubi.add(rule165)

    pattern166 = Pattern(Int(Mul(x_, Pow(x_, n_)), x_), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Not(IntegerQ(n_)), Greater(Add(n_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(5), n_), Integer(14)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(4), n_), Integer(7)), Integer(0))))), (x, n)))
    rule166 = ReplacementRule(pattern166, lambda x, n : Int(ExpandIntegrand(Mul(x, Pow(x, n)), x), x))
    rubi.add(rule166)

    pattern167 = Pattern(Int(Mul(b_, Pow(x_, Integer(2))), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (b, x)))
    rule167 = ReplacementRule(pattern167, lambda b, x : Int(ExpandIntegrand(Mul(b, Pow(x, Integer(2))), x), x))
    rubi.add(rule167)

    pattern168 = Pattern(Int(Mul(x_, Add(a_, x_)), x_), FreeQ(a, x), cons(And(PositiveIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (a, x)))
    rule168 = ReplacementRule(pattern168, lambda a, x : Int(ExpandIntegrand(Mul(x, Add(a, x)), x), x))
    rubi.add(rule168)

    pattern169 = Pattern(Int(Mul(x_, Pow(x_, m_)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(m_), Or(Not(IntegerQ(Integer(1))), Greater(Add(m_, Integer(3)), Integer(0)), Less(Add(Mul(Integer(9), m_), Integer(10)), Integer(0)), And(ZeroQ(Integer(0)), LessEqual(Add(Mul(Integer(7), m_), Integer(4)), Integer(0))))), (m, x)))
    rule169 = ReplacementRule(pattern169, lambda m, x : Int(ExpandIntegrand(Mul(x, Pow(x, m)), x), x))
    rubi.add(rule169)

    pattern170 = Pattern(Int(Pow(x_, Integer(2)), x_), cons(And(NonzeroQ(Integer(0)), PositiveIntegerQ(Integer(1)), Or(Greater(Integer(4), Integer(0)), Less(Integer(19), Integer(0)), Not(IntegerQ(Integer(1))), And(ZeroQ(Integer(0)), LessEqual(Integer(11), Integer(0))))), (x,)))
    rule170 = ReplacementRule(pattern170, lambda x : Int(ExpandIntegrand(Pow(x, Integer(2)), x), x))
    rubi.add(rule170)

    pattern171 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, b, a, c, d, x, n)))
    rule171 = ReplacementRule(pattern171, lambda m, b, a, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule171)

    pattern172 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (b, a, c, d, x, n)))
    rule172 = ReplacementRule(pattern172, lambda b, a, c, d, x, n : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule172)

    pattern173 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, a, c, d, x, n)))
    rule173 = ReplacementRule(pattern173, lambda m, a, c, d, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule173)

    pattern174 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, b, a, c, d, x)))
    rule174 = ReplacementRule(pattern174, lambda m, b, a, c, d, x : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule174)

    pattern175 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, b, a, d, x, n)))
    rule175 = ReplacementRule(pattern175, lambda m, b, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule175)

    pattern176 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, b, a, c, x, n)))
    rule176 = ReplacementRule(pattern176, lambda m, b, a, c, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule176)

    pattern177 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (a, c, d, x, n)))
    rule177 = ReplacementRule(pattern177, lambda a, c, d, x, n : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x))
    rubi.add(rule177)

    pattern178 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule178 = ReplacementRule(pattern178, lambda b, a, c, d, x : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x))
    rubi.add(rule178)

    pattern179 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (b, a, d, x, n)))
    rule179 = ReplacementRule(pattern179, lambda b, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule179)

    pattern180 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (b, a, c, x, n)))
    rule180 = ReplacementRule(pattern180, lambda b, a, c, x, n : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x))
    rubi.add(rule180)

    pattern181 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, a, c, d, x)))
    rule181 = ReplacementRule(pattern181, lambda m, a, c, d, x : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x))
    rubi.add(rule181)

    pattern182 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, a, d, x, n)))
    rule182 = ReplacementRule(pattern182, lambda m, a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule182)

    pattern183 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, a, c, x, n)))
    rule183 = ReplacementRule(pattern183, lambda m, a, c, x, n : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x))
    rubi.add(rule183)

    pattern184 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, b, a, d, x)))
    rule184 = ReplacementRule(pattern184, lambda m, b, a, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule184)

    pattern185 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, b, a, c, x)))
    rule185 = ReplacementRule(pattern185, lambda m, b, a, c, x : Int(ExpandIntegrand(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x))
    rubi.add(rule185)

    pattern186 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, b, a, x, n)))
    rule186 = ReplacementRule(pattern186, lambda m, b, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule186)

    pattern187 = Pattern(Int(Mul(Add(a_, x_), Add(c_, Mul(d_, x_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (a, c, d, x)))
    rule187 = ReplacementRule(pattern187, lambda a, c, d, x : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, Mul(d, x))), x), x))
    rubi.add(rule187)

    pattern188 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Add(a_, x_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (a, d, x, n)))
    rule188 = ReplacementRule(pattern188, lambda a, d, x, n : Int(ExpandIntegrand(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x))
    rubi.add(rule188)

    pattern189 = Pattern(Int(Mul(Add(a_, x_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (a, c, x, n)))
    rule189 = ReplacementRule(pattern189, lambda a, c, x, n : Int(ExpandIntegrand(Mul(Add(a, x), Pow(Add(c, x), n)), x), x))
    rubi.add(rule189)

    pattern190 = Pattern(Int(Mul(d_, x_, Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (b, a, d, x)))
    rule190 = ReplacementRule(pattern190, lambda b, a, d, x : Int(ExpandIntegrand(Mul(d, x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule190)

    pattern191 = Pattern(Int(Mul(Add(a_, Mul(b_, x_)), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0)))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule191 = ReplacementRule(pattern191, lambda b, a, c, x : Int(ExpandIntegrand(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x))
    rubi.add(rule191)

    pattern192 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (b, a, x, n)))
    rule192 = ReplacementRule(pattern192, lambda b, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x))
    rubi.add(rule192)

    pattern193 = Pattern(Int(Mul(d_, x_, Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, a, d, x)))
    rule193 = ReplacementRule(pattern193, lambda m, a, d, x : Int(ExpandIntegrand(Mul(d, x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule193)

    pattern194 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, a, c, x)))
    rule194 = ReplacementRule(pattern194, lambda m, a, c, x : Int(ExpandIntegrand(Mul(Pow(Add(a, x), m), Add(c, x)), x), x))
    rubi.add(rule194)

    pattern195 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(n_), Less(Add(m_, n_, Integer(2)), Integer(0))))), (m, a, x, n)))
    rule195 = ReplacementRule(pattern195, lambda m, a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x))
    rubi.add(rule195)

    pattern196 = Pattern(Int(Mul(x_, Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, b, a, x)))
    rule196 = ReplacementRule(pattern196, lambda m, b, a, x : Int(ExpandIntegrand(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x))
    rubi.add(rule196)

    pattern197 = Pattern(Int(Mul(d_, x_, Add(a_, x_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (a, d, x)))
    rule197 = ReplacementRule(pattern197, lambda a, d, x : Int(ExpandIntegrand(Mul(d, x, Add(a, x)), x), x))
    rubi.add(rule197)

    pattern198 = Pattern(Int(Mul(Add(a_, x_), Add(c_, x_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (a, c, x)))
    rule198 = ReplacementRule(pattern198, lambda a, c, x : Int(ExpandIntegrand(Mul(Add(a, x), Add(c, x)), x), x))
    rubi.add(rule198)

    pattern199 = Pattern(Int(Mul(Pow(x_, n_), Add(a_, x_)), x_), FreeQ(a, x), FreeQ(n, x), cons(And(IntegerQ(n_), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(n_), Less(Add(n_, Integer(3)), Integer(0))))), (a, x, n)))
    rule199 = ReplacementRule(pattern199, lambda a, x, n : Int(ExpandIntegrand(Mul(Pow(x, n), Add(a, x)), x), x))
    rubi.add(rule199)

    pattern200 = Pattern(Int(Mul(x_, Add(a_, Mul(b_, x_))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (b, a, x)))
    rule200 = ReplacementRule(pattern200, lambda b, a, x : Int(ExpandIntegrand(Mul(x, Add(a, Mul(b, x))), x), x))
    rubi.add(rule200)

    pattern201 = Pattern(Int(Mul(x_, Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(m_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Add(m_, Integer(3)), Integer(0))))), (m, a, x)))
    rule201 = ReplacementRule(pattern201, lambda m, a, x : Int(ExpandIntegrand(Mul(x, Pow(Add(a, x), m)), x), x))
    rubi.add(rule201)

    pattern202 = Pattern(Int(Mul(x_, Add(a_, x_)), x_), FreeQ(a, x), cons(And(IntegerQ(Integer(1)), NegativeIntegerQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_)), Not(And(PositiveIntegerQ(Integer(1)), Less(Integer(4), Integer(0))))), (a, x)))
    rule202 = ReplacementRule(pattern202, lambda a, x : Int(ExpandIntegrand(Mul(x, Add(a, x)), x), x))
    rubi.add(rule202)

    pattern203 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x, n)))
    rule203 = ReplacementRule(pattern203, lambda b, a, c, d, x, n : Add(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule203)

    pattern204 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(b_, c_))), (b, c, d, x, n)))
    rule204 = ReplacementRule(pattern204, lambda b, c, d, x, n : Add(Mul(c, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule204)

    pattern205 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x, n)))
    rule205 = ReplacementRule(pattern205, lambda a, c, d, x, n : Add(Mul(Add(Mul(Integer(-1), a, d), c), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule205)

    pattern206 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_))), (b, a, d, x, n)))
    rule206 = ReplacementRule(pattern206, lambda b, a, d, x, n : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d, Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Mul(d, x), n))))
    rubi.add(rule206)

    pattern207 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x, n)))
    rule207 = ReplacementRule(pattern207, lambda b, a, c, x, n : Add(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule207)

    pattern208 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(n_), Greater(n_, Integer(0))), (c, d, x, n)))
    rule208 = ReplacementRule(pattern208, lambda c, d, x, n : Add(Mul(c, Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, Mul(d, x)), n))))
    rubi.add(rule208)

    pattern209 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0))), (b, d, x, n)))
    rule209 = ReplacementRule(pattern209, lambda b, d, x, n : Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Mul(d, x), n)))
    rubi.add(rule209)

    pattern210 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(b_, c_))), (b, c, x, n)))
    rule210 = ReplacementRule(pattern210, lambda b, c, x, n : Add(Mul(c, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule210)

    pattern211 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x, n)))
    rule211 = ReplacementRule(pattern211, lambda a, d, x, n : Add(Mul(Integer(-1), a, d, Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), Integer(-1))), x)), Mul(Pow(n, Integer(-1)), Pow(Mul(d, x), n))))
    rubi.add(rule211)

    pattern212 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x, n)))
    rule212 = ReplacementRule(pattern212, lambda a, c, x, n : Add(Mul(Add(Mul(Integer(-1), a), c), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule212)

    pattern213 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_))), (b, a, x, n)))
    rule213 = ReplacementRule(pattern213, lambda b, a, x, n : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(x, n))))
    rubi.add(rule213)

    pattern214 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0))), (d, x, n)))
    rule214 = ReplacementRule(pattern214, lambda d, x, n : Mul(Pow(n, Integer(-1)), Pow(Mul(d, x), n)))
    rubi.add(rule214)

    pattern215 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(n_), Greater(n_, Integer(0))), (c, x, n)))
    rule215 = ReplacementRule(pattern215, lambda c, x, n : Add(Mul(c, Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(n, Integer(-1)), Pow(Add(c, x), n))))
    rubi.add(rule215)

    pattern216 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0))), (b, x, n)))
    rule216 = ReplacementRule(pattern216, lambda b, x, n : Mul(Pow(b, Integer(-1)), Pow(n, Integer(-1)), Pow(x, n)))
    rubi.add(rule216)

    pattern217 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(And(RationalQ(n_), Greater(n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_))), (a, x, n)))
    rule217 = ReplacementRule(pattern217, lambda a, x, n : Add(Mul(Integer(-1), a, Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), Integer(-1))), x)), Mul(Pow(n, Integer(-1)), Pow(x, n))))
    rubi.add(rule217)

    pattern218 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Greater(n_, Integer(0))), (x, n)))
    rule218 = ReplacementRule(pattern218, lambda x, n : Mul(Pow(n, Integer(-1)), Pow(x, n)))
    rubi.add(rule218)

    pattern219 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x, n)))
    rule219 = ReplacementRule(pattern219, lambda b, a, c, d, x, n : Add(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule219)

    pattern220 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(b_, c_))), (b, c, d, x, n)))
    rule220 = ReplacementRule(pattern220, lambda b, c, d, x, n : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule220)

    pattern221 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x, n)))
    rule221 = ReplacementRule(pattern221, lambda a, c, d, x, n : Add(Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x))))
    rubi.add(rule221)

    pattern222 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_, d_))), (b, a, d, x, n)))
    rule222 = ReplacementRule(pattern222, lambda b, a, d, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Pow(d, Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule222)

    pattern223 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x, n)))
    rule223 = ReplacementRule(pattern223, lambda b, a, c, x, n : Add(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule223)

    pattern224 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(n_), Less(n_, Integer(-1))), (c, d, x, n)))
    rule224 = ReplacementRule(pattern224, lambda c, d, x, n : Add(Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(c, Integer(-1)), Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1)))), x))))
    rubi.add(rule224)

    pattern225 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1))), (b, d, x, n)))
    rule225 = ReplacementRule(pattern225, lambda b, d, x, n : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1)))), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule225)

    pattern226 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(b_, c_))), (b, c, x, n)))
    rule226 = ReplacementRule(pattern226, lambda b, c, x, n : Add(Mul(Pow(c, Integer(-1)), Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x)), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule226)

    pattern227 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x, n)))
    rule227 = ReplacementRule(pattern227, lambda a, d, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Integer(-1))), x))))
    rubi.add(rule227)

    pattern228 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x, n)))
    rule228 = ReplacementRule(pattern228, lambda a, c, x, n : Add(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x))))
    rubi.add(rule228)

    pattern229 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_))), (b, a, x, n)))
    rule229 = ReplacementRule(pattern229, lambda b, a, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), b, Int(Mul(Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Integer(-1))), x)), Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule229)

    pattern230 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1))), (d, x, n)))
    rule230 = ReplacementRule(pattern230, lambda d, x, n : Add(Mul(zoo, Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(zoo, Int(Mul(Pow(x, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1)))), x))))
    rubi.add(rule230)

    pattern231 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(n_), Less(n_, Integer(-1))), (c, x, n)))
    rule231 = ReplacementRule(pattern231, lambda c, x, n : Add(Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Pow(c, Integer(-1)), Int(Mul(Pow(x, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1)))), x))))
    rubi.add(rule231)

    pattern232 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1))), (b, x, n)))
    rule232 = ReplacementRule(pattern232, lambda b, x, n : Add(Mul(zoo, Int(Mul(Pow(b, Integer(-1)), Pow(x, Integer(-1)), Pow(x, Add(n, Integer(1)))), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)))))
    rubi.add(rule232)

    pattern233 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(And(RationalQ(n_), Less(n_, Integer(-1)), NonzeroQ(Mul(Integer(-1), a_))), (a, x, n)))
    rule233 = ReplacementRule(pattern233, lambda a, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(a, Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Integer(-1))), x))))
    rubi.add(rule233)

    pattern234 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(n_, Integer(-1))), (x, n)))
    rule234 = ReplacementRule(pattern234, lambda x, n : Add(Mul(zoo, Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1))), Mul(zoo, Int(Mul(Pow(x, Integer(-1)), Pow(x, Add(n, Integer(1)))), x))))
    rubi.add(rule234)

    pattern235 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule235 = ReplacementRule(pattern235, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule235)

    pattern236 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(c_), (b, c, d, x)))
    rule236 = ReplacementRule(pattern236, lambda b, c, d, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule236)

    pattern237 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule237 = ReplacementRule(pattern237, lambda a, c, d, x : With(List(Set(q, Rt(Add(Mul(Integer(-1), a, d), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule237)

    pattern238 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_)), (b, a, d, x)))
    rule238 = ReplacementRule(pattern238, lambda b, a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule238)

    pattern239 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule239 = ReplacementRule(pattern239, lambda b, a, c, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule239)

    pattern240 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(PosQ(c_), (c, d, x)))
    rule240 = ReplacementRule(pattern240, lambda c, d, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule240)

    pattern241 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(PosQ(Integer(0)), (b, d, x)))
    rule241 = ReplacementRule(pattern241, lambda b, d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule241)

    pattern242 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(PosQ(c_), (b, c, x)))
    rule242 = ReplacementRule(pattern242, lambda b, c, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule242)

    pattern243 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(PosQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule243 = ReplacementRule(pattern243, lambda a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), a, d), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule243)

    pattern244 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(PosQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule244 = ReplacementRule(pattern244, lambda a, c, x : With(List(Set(q, Rt(Add(Mul(Integer(-1), a), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule244)

    pattern245 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)))), (b, a, x)))
    rule245 = ReplacementRule(pattern245, lambda b, a, x : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule245)

    pattern246 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), FreeQ(d, x), cons(PosQ(Integer(0)), (d, x)))
    rule246 = ReplacementRule(pattern246, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule246)

    pattern247 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(c, x), cons(PosQ(c_), (c, x)))
    rule247 = ReplacementRule(pattern247, lambda c, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule247)

    pattern248 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-4/3))), x_), FreeQ(b, x), cons(PosQ(Integer(0)), (b, x)))
    rule248 = ReplacementRule(pattern248, lambda b, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule248)

    pattern249 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(PosQ(Mul(Integer(-1), a_)), (a, x)))
    rule249 = ReplacementRule(pattern249, lambda a, x : With(List(Set(q, Rt(Mul(Integer(-1), a), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule249)

    pattern250 = Pattern(Int(Pow(x_, Integer(-4/3)), x_), cons(PosQ(Integer(0)), (x,)))
    rule250 = ReplacementRule(pattern250, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule250)

    pattern251 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule251 = ReplacementRule(pattern251, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(a, d), Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule251)

    pattern252 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(c_), (b, c, d, x)))
    rule252 = ReplacementRule(pattern252, lambda b, c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule252)

    pattern253 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule253 = ReplacementRule(pattern253, lambda a, c, d, x : With(List(Set(q, Rt(Add(Mul(a, d), Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule253)

    pattern254 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_)), (b, a, d, x)))
    rule254 = ReplacementRule(pattern254, lambda b, a, d, x : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule254)

    pattern255 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule255 = ReplacementRule(pattern255, lambda b, a, c, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(a, Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule255)

    pattern256 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(NegQ(c_), (c, d, x)))
    rule256 = ReplacementRule(pattern256, lambda c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule256)

    pattern257 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(NegQ(Integer(0)), (b, d, x)))
    rule257 = ReplacementRule(pattern257, lambda b, d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule257)

    pattern258 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(NegQ(c_), (b, c, x)))
    rule258 = ReplacementRule(pattern258, lambda b, c, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule258)

    pattern259 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(NegQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule259 = ReplacementRule(pattern259, lambda a, d, x : With(List(Set(q, Rt(Mul(a, d), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule259)

    pattern260 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(NegQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule260 = ReplacementRule(pattern260, lambda a, c, x : With(List(Set(q, Rt(Add(a, Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule260)

    pattern261 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)))), (b, a, x)))
    rule261 = ReplacementRule(pattern261, lambda b, a, x : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule261)

    pattern262 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-1/3))), x_), FreeQ(d, x), cons(NegQ(Integer(0)), (d, x)))
    rule262 = ReplacementRule(pattern262, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule262)

    pattern263 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-1/3))), x_), FreeQ(c, x), cons(NegQ(c_), (c, x)))
    rule263 = ReplacementRule(pattern263, lambda c, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule263)

    pattern264 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-4/3))), x_), FreeQ(b, x), cons(NegQ(Integer(0)), (b, x)))
    rule264 = ReplacementRule(pattern264, lambda b, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule264)

    pattern265 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(NegQ(Mul(Integer(-1), a_)), (a, x)))
    rule265 = ReplacementRule(pattern265, lambda a, x : With(List(Set(q, Rt(a, Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule265)

    pattern266 = Pattern(Int(Pow(x_, Integer(-4/3)), x_), cons(NegQ(Integer(0)), (x,)))
    rule266 = ReplacementRule(pattern266, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(1/2), Pow(q, Integer(-1)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule266)

    pattern267 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule267 = ReplacementRule(pattern267, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a, d), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule267)

    pattern268 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(c_), (b, c, d, x)))
    rule268 = ReplacementRule(pattern268, lambda b, c, d, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule268)

    pattern269 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(PosQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule269 = ReplacementRule(pattern269, lambda a, c, d, x : With(List(Set(q, Rt(Add(Mul(Integer(-1), a, d), c), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule269)

    pattern270 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_)), (b, a, d, x)))
    rule270 = ReplacementRule(pattern270, lambda b, a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule270)

    pattern271 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(PosQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule271 = ReplacementRule(pattern271, lambda b, a, c, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(Integer(-1), a), Mul(b, c))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule271)

    pattern272 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(PosQ(c_), (c, d, x)))
    rule272 = ReplacementRule(pattern272, lambda c, d, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule272)

    pattern273 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(PosQ(Integer(0)), (b, d, x)))
    rule273 = ReplacementRule(pattern273, lambda b, d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule273)

    pattern274 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(PosQ(c_), (b, c, x)))
    rule274 = ReplacementRule(pattern274, lambda b, c, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule274)

    pattern275 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(PosQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule275 = ReplacementRule(pattern275, lambda a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), a, d), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule275)

    pattern276 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(PosQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule276 = ReplacementRule(pattern276, lambda a, c, x : With(List(Set(q, Rt(Add(Mul(Integer(-1), a), c), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule276)

    pattern277 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(PosQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)))), (b, a, x)))
    rule277 = ReplacementRule(pattern277, lambda b, a, x : With(List(Set(q, Rt(Mul(Integer(-1), a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule277)

    pattern278 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(d, x), cons(PosQ(Integer(0)), (d, x)))
    rule278 = ReplacementRule(pattern278, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule278)

    pattern279 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(c, x), cons(PosQ(c_), (c, x)))
    rule279 = ReplacementRule(pattern279, lambda c, x : With(List(Set(q, Rt(c, Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule279)

    pattern280 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-5/3))), x_), FreeQ(b, x), cons(PosQ(Integer(0)), (b, x)))
    rule280 = ReplacementRule(pattern280, lambda b, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(-1), Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule280)

    pattern281 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(PosQ(Mul(Integer(-1), a_)), (a, x)))
    rule281 = ReplacementRule(pattern281, lambda a, x : With(List(Set(q, Rt(Mul(Integer(-1), a), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule281)

    pattern282 = Pattern(Int(Pow(x_, Integer(-5/3)), x_), cons(PosQ(Integer(0)), (x,)))
    rule282 = ReplacementRule(pattern282, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(-1), Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, Mul(Integer(-1), x)), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule282)

    pattern283 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule283 = ReplacementRule(pattern283, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(Mul(a, d), Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule283)

    pattern284 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(c_), (b, c, d, x)))
    rule284 = ReplacementRule(pattern284, lambda b, c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule284)

    pattern285 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(NegQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule285 = ReplacementRule(pattern285, lambda a, c, d, x : With(List(Set(q, Rt(Add(Mul(a, d), Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule285)

    pattern286 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)), d_)), (b, a, d, x)))
    rule286 = ReplacementRule(pattern286, lambda b, a, d, x : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule286)

    pattern287 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(NegQ(Mul(Pow(b_, Integer(-1)), Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule287 = ReplacementRule(pattern287, lambda b, a, c, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), Add(a, Mul(Integer(-1), b, c))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule287)

    pattern288 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(NegQ(c_), (c, d, x)))
    rule288 = ReplacementRule(pattern288, lambda c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, Mul(d, x)), Integer(1/3)))))))
    rubi.add(rule288)

    pattern289 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(NegQ(Integer(0)), (b, d, x)))
    rule289 = ReplacementRule(pattern289, lambda b, d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule289)

    pattern290 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(NegQ(c_), (b, c, x)))
    rule290 = ReplacementRule(pattern290, lambda b, c, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule290)

    pattern291 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(NegQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule291 = ReplacementRule(pattern291, lambda a, d, x : With(List(Set(q, Rt(Mul(a, d), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule291)

    pattern292 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(NegQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule292 = ReplacementRule(pattern292, lambda a, c, x : With(List(Set(q, Rt(Add(a, Mul(Integer(-1), c)), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule292)

    pattern293 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(NegQ(Mul(Integer(-1), a_, Pow(b_, Integer(-1)))), (b, a, x)))
    rule293 = ReplacementRule(pattern293, lambda b, a, x : With(List(Set(q, Rt(Mul(a, Pow(b, Integer(-1))), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, Mul(b, x)), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule293)

    pattern294 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(d, x), cons(NegQ(Integer(0)), (d, x)))
    rule294 = ReplacementRule(pattern294, lambda d, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Mul(d, x), Integer(1/3)))))))
    rubi.add(rule294)

    pattern295 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(c, x), cons(NegQ(c_), (c, x)))
    rule295 = ReplacementRule(pattern295, lambda c, x : With(List(Set(q, Rt(Mul(Integer(-1), c), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(Add(c, x), Integer(1/3)))))))
    rubi.add(rule295)

    pattern296 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-5/3))), x_), FreeQ(b, x), cons(NegQ(Integer(0)), (b, x)))
    rule296 = ReplacementRule(pattern296, lambda b, x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Log(RemoveContent(Mul(b, x), x))), Mul(Integer(3/2), Pow(b, Integer(-1)), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule296)

    pattern297 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(NegQ(Mul(Integer(-1), a_)), (a, x)))
    rule297 = ReplacementRule(pattern297, lambda a, x : With(List(Set(q, Rt(a, Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(Add(a, x), x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule297)

    pattern298 = Pattern(Int(Pow(x_, Integer(-5/3)), x_), cons(NegQ(Integer(0)), (x,)))
    rule298 = ReplacementRule(pattern298, lambda x : With(List(Set(q, Rt(Integer(0), Integer(3)))), Add(Mul(Integer(3/2), Pow(q, Integer(-1)), Subst(Int(Pow(Add(Pow(q, Integer(2)), Mul(Integer(-1), q, x), Pow(x, Integer(2))), Integer(-1)), x), x, Pow(x, Integer(1/3)))), Mul(Integer(-1), Integer(1/2), Pow(q, Integer(-2)), Log(RemoveContent(x, x))), Mul(Integer(3/2), Pow(q, Integer(-2)), Subst(Int(Pow(Add(q, x), Integer(-1)), x), x, Pow(x, Integer(1/3)))))))
    rubi.add(rule298)

    pattern299 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x, n)))
    rule299 = ReplacementRule(pattern299, lambda b, a, c, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule299)

    pattern300 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_))), (b, c, d, x, n)))
    rule300 = ReplacementRule(pattern300, lambda b, c, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule300)

    pattern301 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x, n)))
    rule301 = ReplacementRule(pattern301, lambda a, c, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule301)

    pattern302 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_))), (b, a, d, x, n)))
    rule302 = ReplacementRule(pattern302, lambda b, a, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule302)

    pattern303 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x, n)))
    rule303 = ReplacementRule(pattern303, lambda b, a, c, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule303)

    pattern304 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (c, d, x, n)))
    rule304 = ReplacementRule(pattern304, lambda c, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, Mul(d, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule304)

    pattern305 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (b, d, x, n)))
    rule305 = ReplacementRule(pattern305, lambda b, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(b, Integer(-1)), Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule305)

    pattern306 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_))), (b, c, x, n)))
    rule306 = ReplacementRule(pattern306, lambda b, c, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), b, c), Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule306)

    pattern307 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x, n)))
    rule307 = ReplacementRule(pattern307, lambda a, d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(a, d), Pow(x, p)), Integer(-1))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule307)

    pattern308 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x, n)))
    rule308 = ReplacementRule(pattern308, lambda a, c, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule308)

    pattern309 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_))), (b, a, x, n)))
    rule309 = ReplacementRule(pattern309, lambda b, a, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Mul(b, Pow(x, p))), Integer(-1))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule309)

    pattern310 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (d, x, n)))
    rule310 = ReplacementRule(pattern310, lambda d, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(Mul(d, x), Pow(p, Integer(-1)))))))
    rubi.add(rule310)

    pattern311 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (c, x, n)))
    rule311 = ReplacementRule(pattern311, lambda c, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), c), Pow(x, p)), Integer(-1))), x), x, Pow(Add(c, x), Pow(p, Integer(-1)))))))
    rubi.add(rule311)

    pattern312 = Pattern(Int(Mul(Pow(b_, Integer(-1)), Pow(x_, Integer(-1)), Pow(x_, n_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (b, x, n)))
    rule312 = ReplacementRule(pattern312, lambda b, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(b, Integer(-1)), Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule312)

    pattern313 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), cons(And(RationalQ(n_), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_))), (a, x, n)))
    rule313 = ReplacementRule(pattern313, lambda a, x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1))), Pow(Add(a, Pow(x, p)), Integer(-1))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule313)

    pattern314 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(n_), Less(Integer(-1), n_, Integer(0))), (x, n)))
    rule314 = ReplacementRule(pattern314, lambda x, n : With(List(Set(p, Denominator(n))), Mul(p, Subst(Int(Mul(Pow(x, Mul(Integer(-1), p)), Pow(x, Add(Mul(p, Add(n, Integer(1))), Integer(-1)))), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule314)

    pattern315 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(Not(IntegerQ(n_)), (d, n, c, x)))
    rule315 = ReplacementRule(pattern315, lambda d, n, c, x : Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule315)

    pattern316 = Pattern(Int(Mul(Pow(x_, Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(n, x), cons(Not(IntegerQ(n_)), (n, c, x)))
    rule316 = ReplacementRule(pattern316, lambda n, c, x : Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule316)

    pattern317 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (a, b, c, d, x, n)))
    rule317 = ReplacementRule(pattern317, lambda a, b, c, d, x, n : Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(b, Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))))
    rubi.add(rule317)

    pattern318 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x, n)))
    rule318 = ReplacementRule(pattern318, lambda a, c, d, x, n : Mul(Integer(-1), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))))
    rubi.add(rule318)

    pattern319 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, b, d, x, n)))
    rule319 = ReplacementRule(pattern319, lambda a, b, d, x, n : Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), b, x)))))
    rubi.add(rule319)

    pattern320 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (a, b, c, x, n)))
    rule320 = ReplacementRule(pattern320, lambda a, b, c, x, n : Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(c, x))))))
    rubi.add(rule320)

    pattern321 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x, n)))
    rule321 = ReplacementRule(pattern321, lambda a, d, x, n : Mul(Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), x)))))
    rubi.add(rule321)

    pattern322 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1)), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x, n)))
    rule322 = ReplacementRule(pattern322, lambda a, c, x, n : Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(c, x))))))
    rubi.add(rule322)

    pattern323 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), Integer(-1))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_))), (a, b, x, n)))
    rule323 = ReplacementRule(pattern323, lambda a, b, x, n : Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), b, x)))))
    rubi.add(rule323)

    pattern324 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), Integer(-1))), x_), FreeQ(a, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), NonzeroQ(Mul(Integer(-1), a_))), (a, x, n)))
    rule324 = ReplacementRule(pattern324, lambda a, x, n : Mul(Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Integer(1), Add(n, Integer(1)), Add(n, Integer(2)), TogetherSimplify(Mul(Integer(-1), Pow(a, Integer(-1)), x)))))
    rubi.add(rule324)

    pattern325 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, a, c, d, x, n)))
    rule325 = ReplacementRule(pattern325, lambda m, b, a, c, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule325)

    pattern326 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, c, d, x, n)))
    rule326 = ReplacementRule(pattern326, lambda m, b, c, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule326)

    pattern327 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, a, c, d, x, n)))
    rule327 = ReplacementRule(pattern327, lambda m, a, c, d, x, n : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule327)

    pattern328 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, a, d, x, n)))
    rule328 = ReplacementRule(pattern328, lambda m, b, a, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule328)

    pattern329 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, a, c, x, n)))
    rule329 = ReplacementRule(pattern329, lambda m, b, a, c, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule329)

    pattern330 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, c, d, x, n)))
    rule330 = ReplacementRule(pattern330, lambda m, c, d, x, n : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule330)

    pattern331 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, d, x, n)))
    rule331 = ReplacementRule(pattern331, lambda m, b, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule331)

    pattern332 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, c, x, n)))
    rule332 = ReplacementRule(pattern332, lambda m, b, c, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule332)

    pattern333 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, a, d, x, n)))
    rule333 = ReplacementRule(pattern333, lambda m, a, d, x, n : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule333)

    pattern334 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, a, c, x, n)))
    rule334 = ReplacementRule(pattern334, lambda m, a, c, x, n : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule334)

    pattern335 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, a, x, n)))
    rule335 = ReplacementRule(pattern335, lambda m, b, a, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule335)

    pattern336 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, d, x, n)))
    rule336 = ReplacementRule(pattern336, lambda m, d, x, n : Add(Mul(Integer(-1), d, n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule336)

    pattern337 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, c, x, n)))
    rule337 = ReplacementRule(pattern337, lambda m, c, x, n : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule337)

    pattern338 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, b, x, n)))
    rule338 = ReplacementRule(pattern338, lambda m, b, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Mul(b, x), Add(m, Integer(1)))), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule338)

    pattern339 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, a, x, n)))
    rule339 = ReplacementRule(pattern339, lambda m, a, x, n : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule339)

    pattern340 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), Not(And(IntegerQ(n_), Not(IntegerQ(m_)))), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), LessEqual(Add(m_, n_, Integer(2)), Integer(0)), Or(FractionQ(m_), GreaterEqual(Add(m_, Mul(Integer(2), n_), Integer(1)), Integer(0)))))), (m, x, n)))
    rule340 = ReplacementRule(pattern340, lambda m, x, n : Add(Mul(Integer(-1), n, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(-1)))), x)), Mul(Pow(x, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule340)

    pattern341 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, b, a, c, d, x, n)))
    rule341 = ReplacementRule(pattern341, lambda m, b, a, c, d, x, n : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule341)

    pattern342 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, b, c, d, x, n)))
    rule342 = ReplacementRule(pattern342, lambda m, b, c, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule342)

    pattern343 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, a, c, d, x, n)))
    rule343 = ReplacementRule(pattern343, lambda m, a, c, d, x, n : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))
    rubi.add(rule343)

    pattern344 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, b, a, d, x, n)))
    rule344 = ReplacementRule(pattern344, lambda m, b, a, d, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule344)

    pattern345 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, b, a, c, x, n)))
    rule345 = ReplacementRule(pattern345, lambda m, b, a, c, x, n : Add(Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule345)

    pattern346 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, c, d, x, n)))
    rule346 = ReplacementRule(pattern346, lambda m, c, d, x, n : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n)), x)), Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule346)

    pattern347 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, b, d, x, n)))
    rule347 = ReplacementRule(pattern347, lambda m, b, d, x, n : Add(Mul(zoo, Pow(b, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n)), x)), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule347)

    pattern348 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, b, c, x, n)))
    rule348 = ReplacementRule(pattern348, lambda m, b, c, x, n : Add(Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule348)

    pattern349 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, a, d, x, n)))
    rule349 = ReplacementRule(pattern349, lambda m, a, d, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1)))), x)), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule349)

    pattern350 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, a, c, x, n)))
    rule350 = ReplacementRule(pattern350, lambda m, a, c, x, n : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule350)

    pattern351 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, b, a, x, n)))
    rule351 = ReplacementRule(pattern351, lambda m, b, a, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1)))), x))))
    rubi.add(rule351)

    pattern352 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, d, x, n)))
    rule352 = ReplacementRule(pattern352, lambda m, d, x, n : Add(Mul(zoo, d, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n)), x)), Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule352)

    pattern353 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(c_), Less(m_, n_)))))), (m, c, x, n)))
    rule353 = ReplacementRule(pattern353, lambda m, c, x, n : Add(Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n)), x))))
    rubi.add(rule353)

    pattern354 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, b, x, n)))
    rule354 = ReplacementRule(pattern354, lambda m, b, x, n : Add(Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(b, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1)))), x))))
    rubi.add(rule354)

    pattern355 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(Less(m_, Integer(-1)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(a_), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, a, x, n)))
    rule355 = ReplacementRule(pattern355, lambda m, a, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1)))), x))))
    rubi.add(rule355)

    pattern356 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Less(m_, Integer(-1)), RationalQ(m_, n_), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(Less(n_, Integer(-1)), Or(ZeroQ(Integer(0)), And(IntegerQ(n_), NonzeroQ(Integer(0)), Less(m_, n_)))))), (m, x, n)))
    rule356 = ReplacementRule(pattern356, lambda m, x, n : Add(Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(Add(m, Integer(1)), Integer(-1)), Add(m, n, Integer(2)), Int(Mul(Pow(x, n), Pow(x, Add(m, Integer(1)))), x))))
    rubi.add(rule356)

    pattern357 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, a, c, d, x, n)))
    rule357 = ReplacementRule(pattern357, lambda m, b, a, c, d, x, n : Add(Mul(Pow(b, Integer(-1)), n, Add(Mul(Integer(-1), a, d), Mul(b, c)), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule357)

    pattern358 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, c, d, x, n)))
    rule358 = ReplacementRule(pattern358, lambda m, b, c, d, x, n : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule358)

    pattern359 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, a, c, d, x, n)))
    rule359 = ReplacementRule(pattern359, lambda m, a, c, d, x, n : Add(Mul(n, Add(Mul(Integer(-1), a, d), c), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule359)

    pattern360 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, a, d, x, n)))
    rule360 = ReplacementRule(pattern360, lambda m, b, a, d, x, n : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), m)), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule360)

    pattern361 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, a, c, x, n)))
    rule361 = ReplacementRule(pattern361, lambda m, b, a, c, x, n : Add(Mul(Pow(b, Integer(-1)), n, Add(Mul(Integer(-1), a), Mul(b, c)), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule361)

    pattern362 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, c, d, x, n)))
    rule362 = ReplacementRule(pattern362, lambda m, c, d, x, n : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule362)

    pattern363 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, d, x, n)))
    rule363 = ReplacementRule(pattern363, lambda m, b, d, x, n : Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule363)

    pattern364 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(b_, c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, c, x, n)))
    rule364 = ReplacementRule(pattern364, lambda m, b, c, x, n : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule364)

    pattern365 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_, d_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, a, d, x, n)))
    rule365 = ReplacementRule(pattern365, lambda m, a, d, x, n : Add(Mul(Integer(-1), a, d, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), Add(n, Integer(-1))), Pow(Add(a, x), m)), x)), Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule365)

    pattern366 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, a, c, x, n)))
    rule366 = ReplacementRule(pattern366, lambda m, a, c, x, n : Add(Mul(n, Add(Mul(Integer(-1), a), c), Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule366)

    pattern367 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, a, x, n)))
    rule367 = ReplacementRule(pattern367, lambda m, b, a, x, n : Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, Mul(b, x)), m)), x)), Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule367)

    pattern368 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, d, x, n)))
    rule368 = ReplacementRule(pattern368, lambda m, d, x, n : Mul(Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), n), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule368)

    pattern369 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, c, x, n)))
    rule369 = ReplacementRule(pattern369, lambda m, c, x, n : Add(Mul(c, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, m), Pow(Add(c, x), Add(n, Integer(-1)))), x)), Mul(Pow(x, Add(m, Integer(1))), Pow(Add(c, x), n), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule369)

    pattern370 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, b, x, n)))
    rule370 = ReplacementRule(pattern370, lambda m, b, x, n : Mul(Pow(b, Integer(-1)), Pow(x, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule370)

    pattern371 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(Greater(n_, Integer(0)), RationalQ(m_, n_), NonzeroQ(Mul(Integer(-1), a_)), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, a, x, n)))
    rule371 = ReplacementRule(pattern371, lambda m, a, x, n : Add(Mul(Integer(-1), a, n, Pow(Add(m, n, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Add(n, Integer(-1))), Pow(Add(a, x), m)), x)), Mul(Pow(x, n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1)))))
    rubi.add(rule371)

    pattern372 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), Greater(n_, Integer(0)), RationalQ(m_, n_), Unequal(Add(m_, n_, Integer(1)), Integer(0)), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_), Not(And(IntegerQ(Add(m_, n_)), Less(Add(m_, n_, Integer(2)), Integer(0)))), Not(And(PositiveIntegerQ(m_), Or(Not(IntegerQ(n_)), Less(Integer(0), m_, n_))))), (m, x, n)))
    rule372 = ReplacementRule(pattern372, lambda m, x, n : Mul(Pow(x, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, n, Integer(1)), Integer(-1))))
    rubi.add(rule372)

    pattern373 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(b_, d_))), (a, b, c, d, x)))
    rule373 = ReplacementRule(pattern373, lambda a, b, c, d, x : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))), Mul(Integer(-1), b, x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule373)

    pattern374 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(d_, Integer(1)))), (a, c, d, x)))
    rule374 = ReplacementRule(pattern374, lambda a, c, d, x : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(x, Integer(2))), Mul(Integer(-1), x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule374)

    pattern375 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(PositiveQ(a_), ZeroQ(Add(b_, d_))), (a, b, d, x)))
    rule375 = ReplacementRule(pattern375, lambda a, b, d, x : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, b, x), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule375)

    pattern376 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveQ(Add(a_, c_)), ZeroQ(Add(b_, Integer(1)))), (a, b, c, x)))
    rule376 = ReplacementRule(pattern376, lambda a, b, c, x : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))), Mul(Integer(-1), b, x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule376)

    pattern377 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(PositiveQ(a_), ZeroQ(Add(d_, Integer(1)))), (a, d, x)))
    rule377 = ReplacementRule(pattern377, lambda a, d, x : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, x), Mul(Integer(-1), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule377)

    pattern378 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(ZeroQ(Integer(2)), PositiveQ(Add(a_, c_))), (a, c, x)))
    rule378 = ReplacementRule(pattern378, lambda a, c, x : Int(Pow(Sqrt(Add(Mul(a, c), Mul(Integer(-1), Pow(x, Integer(2))), Mul(Integer(-1), x, Add(a, Mul(Integer(-1), c))))), Integer(-1)), x))
    rubi.add(rule378)

    pattern379 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(PositiveQ(a_), ZeroQ(Add(b_, Integer(1)))), (a, b, x)))
    rule379 = ReplacementRule(pattern379, lambda a, b, x : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, b, x), Mul(Integer(-1), Pow(b, Integer(2)), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule379)

    pattern380 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), cons(And(PositiveQ(a_), ZeroQ(Integer(2))), (a, x)))
    rule380 = ReplacementRule(pattern380, lambda a, x : Int(Pow(Sqrt(Add(Mul(Integer(-1), a, x), Mul(Integer(-1), Pow(x, Integer(2))))), Integer(-1)), x))
    rubi.add(rule380)

    pattern381 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(b_), PositiveQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule381 = ReplacementRule(pattern381, lambda b, a, c, d, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(b, c), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule381)

    pattern382 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(b_), PositiveQ(Mul(b_, c_))), (b, c, d, x)))
    rule382 = ReplacementRule(pattern382, lambda b, c, d, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(b, c), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule382)

    pattern383 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(Integer(1)), PositiveQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x)))
    rule383 = ReplacementRule(pattern383, lambda a, c, d, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), c, Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule383)

    pattern384 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(PositiveQ(b_), PositiveQ(Mul(Integer(-1), a_, d_))), (b, a, d, x)))
    rule384 = ReplacementRule(pattern384, lambda b, a, d, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule384)

    pattern385 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveQ(b_), PositiveQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule385 = ReplacementRule(pattern385, lambda b, a, c, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Mul(b, c), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule385)

    pattern386 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(PositiveQ(Integer(1)), PositiveQ(c_)), (c, d, x)))
    rule386 = ReplacementRule(pattern386, lambda c, d, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule386)

    pattern387 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(d, x), cons(And(PositiveQ(Integer(0)), PositiveQ(b_)), (b, d, x)))
    rule387 = ReplacementRule(pattern387, lambda b, d, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Mul(d, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule387)

    pattern388 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), cons(And(PositiveQ(b_), PositiveQ(Mul(b_, c_))), (b, c, x)))
    rule388 = ReplacementRule(pattern388, lambda b, c, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(b, c), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule388)

    pattern389 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(PositiveQ(Integer(1)), PositiveQ(Mul(Integer(-1), a_, d_))), (a, d, x)))
    rule389 = ReplacementRule(pattern389, lambda a, d, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, Integer(2))))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule389)

    pattern390 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PositiveQ(Integer(1)), PositiveQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x)))
    rule390 = ReplacementRule(pattern390, lambda a, c, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule390)

    pattern391 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(PositiveQ(b_), PositiveQ(Mul(Integer(-1), a_))), (b, a, x)))
    rule391 = ReplacementRule(pattern391, lambda b, a, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule391)

    pattern392 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), FreeQ(d, x), cons(And(PositiveQ(Integer(0)), PositiveQ(Integer(1))), (d, x)))
    rule392 = ReplacementRule(pattern392, lambda d, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Mul(d, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule392)

    pattern393 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(c, x), cons(And(PositiveQ(Integer(1)), PositiveQ(c_)), (c, x)))
    rule393 = ReplacementRule(pattern393, lambda c, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule393)

    pattern394 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(b_, x_), Integer(-1/2))), x_), FreeQ(b, x), cons(And(PositiveQ(Integer(0)), PositiveQ(b_)), (b, x)))
    rule394 = ReplacementRule(pattern394, lambda b, x : Mul(Integer(2), Pow(Sqrt(b), Integer(-1)), Subst(Int(Pow(Sqrt(Pow(x, Integer(2))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule394)

    pattern395 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), cons(And(PositiveQ(Integer(1)), PositiveQ(Mul(Integer(-1), a_))), (a, x)))
    rule395 = ReplacementRule(pattern395, lambda a, x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule395)

    pattern396 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(PositiveQ(Integer(0)), PositiveQ(Integer(1))), (x,)))
    rule396 = ReplacementRule(pattern396, lambda x : Mul(Integer(2), Pow(Sqrt(Integer(1)), Integer(-1)), Subst(Int(Pow(Sqrt(Pow(x, Integer(2))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule396)

    pattern397 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(b_, Mul(Integer(-1), d_))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule397 = ReplacementRule(pattern397, lambda b, a, c, d, x : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule397)

    pattern398 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(Mul(b_, c_)), ZeroQ(Add(b_, Mul(Integer(-1), d_)))), (b, c, d, x)))
    rule398 = ReplacementRule(pattern398, lambda b, c, d, x : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule398)

    pattern399 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(ZeroQ(Add(Mul(Integer(-1), d_), Integer(1))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x)))
    rule399 = ReplacementRule(pattern399, lambda a, c, d, x : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule399)

    pattern400 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(ZeroQ(Add(b_, Integer(-1))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule400 = ReplacementRule(pattern400, lambda b, a, c, x : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, Mul(b, x))))))
    rubi.add(rule400)

    pattern401 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), ZeroQ(Add(Mul(Integer(-1), d_), Integer(1)))), (c, d, x)))
    rule401 = ReplacementRule(pattern401, lambda c, d, x : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule401)

    pattern402 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), cons(And(NonzeroQ(Mul(b_, c_)), ZeroQ(Add(b_, Integer(-1)))), (b, c, x)))
    rule402 = ReplacementRule(pattern402, lambda b, c, x : Mul(Integer(2), Pow(b, Integer(-1)), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Mul(b, x)))))
    rubi.add(rule402)

    pattern403 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(ZeroQ(Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x)))
    rule403 = ReplacementRule(pattern403, lambda a, c, x : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(Mul(Integer(-1), a), c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(Add(a, x)))))
    rubi.add(rule403)

    pattern404 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), ZeroQ(Integer(0))), (c, x)))
    rule404 = ReplacementRule(pattern404, lambda c, x : Mul(Integer(2), Subst(Int(Pow(Sqrt(Add(c, Pow(x, Integer(2)))), Integer(-1)), x), x, Sqrt(x))))
    rubi.add(rule404)

    pattern405 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), (b, a, c, d, x)))
    rule405 = ReplacementRule(pattern405, lambda b, a, c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule405)

    pattern406 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Mul(b_, c_)), (b, c, d, x)))
    rule406 = ReplacementRule(pattern406, lambda b, c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule406)

    pattern407 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), (a, c, d, x)))
    rule407 = ReplacementRule(pattern407, lambda a, c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule407)

    pattern408 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(NonzeroQ(Mul(Integer(-1), a_, d_)), (b, a, d, x)))
    rule408 = ReplacementRule(pattern408, lambda b, a, d, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(Mul(d, x)), Integer(-1)), Sqrt(Add(a, Mul(b, x)))))))
    rubi.add(rule408)

    pattern409 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), (b, a, c, x)))
    rule409 = ReplacementRule(pattern409, lambda b, a, c, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Add(a, Mul(b, x))), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule409)

    pattern410 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, Mul(d_, x_)), Integer(-1/2))), x_), FreeQ(c, x), FreeQ(d, x), cons(NonzeroQ(c_), (c, d, x)))
    rule410 = ReplacementRule(pattern410, lambda c, d, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Add(c, Mul(d, x))), Integer(-1))))))
    rubi.add(rule410)

    pattern411 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(d, x), cons(NonzeroQ(Integer(0)), (b, d, x)))
    rule411 = ReplacementRule(pattern411, lambda b, d, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), d, Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Mul(d, x)), Integer(-1))))))
    rubi.add(rule411)

    pattern412 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(b, x), FreeQ(c, x), cons(NonzeroQ(Mul(b_, c_)), (b, c, x)))
    rule412 = ReplacementRule(pattern412, lambda b, c, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Sqrt(Mul(b, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule412)

    pattern413 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(d, x), cons(NonzeroQ(Mul(Integer(-1), a_, d_)), (a, d, x)))
    rule413 = ReplacementRule(pattern413, lambda a, d, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Pow(Sqrt(Mul(d, x)), Integer(-1)), Sqrt(Add(a, x))))))
    rubi.add(rule413)

    pattern414 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(c, x), cons(NonzeroQ(Add(Mul(Integer(-1), a_), c_)), (a, c, x)))
    rule414 = ReplacementRule(pattern414, lambda a, c, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(Add(a, x)), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule414)

    pattern415 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/2))), x_), FreeQ(a, x), FreeQ(b, x), cons(NonzeroQ(Mul(Integer(-1), a_)), (b, a, x)))
    rule415 = ReplacementRule(pattern415, lambda b, a, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Add(a, Mul(b, x)))))))
    rubi.add(rule415)

    pattern416 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(d_, x_), Integer(-1/2))), x_), FreeQ(d, x), cons(NonzeroQ(Integer(0)), (d, x)))
    rule416 = ReplacementRule(pattern416, lambda d, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), d, Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Mul(d, x)), Integer(-1))))))
    rubi.add(rule416)

    pattern417 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(c_, x_), Integer(-1/2))), x_), FreeQ(c, x), cons(NonzeroQ(c_), (c, x)))
    rule417 = ReplacementRule(pattern417, lambda c, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Sqrt(x), Pow(Sqrt(Add(c, x)), Integer(-1))))))
    rubi.add(rule417)

    pattern418 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Mul(b_, x_), Integer(-1/2))), x_), FreeQ(b, x), cons(NonzeroQ(Integer(0)), (b, x)))
    rule418 = ReplacementRule(pattern418, lambda b, x : Mul(Integer(2), Subst(Int(Pow(Add(b, Mul(Integer(-1), Pow(x, Integer(2)))), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Mul(b, x))))))
    rubi.add(rule418)

    pattern419 = Pattern(Int(Mul(Pow(x_, Integer(-1/2)), Pow(Add(a_, x_), Integer(-1/2))), x_), FreeQ(a, x), cons(NonzeroQ(Mul(Integer(-1), a_)), (a, x)))
    rule419 = ReplacementRule(pattern419, lambda a, x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Mul(Pow(Sqrt(x), Integer(-1)), Sqrt(Add(a, x))))))
    rubi.add(rule419)

    pattern420 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(NonzeroQ(Integer(0)), (x,)))
    rule420 = ReplacementRule(pattern420, lambda x : Mul(Integer(2), Subst(Int(Pow(Add(Mul(Integer(-1), Pow(x, Integer(2))), Integer(1)), Integer(-1)), x), x, Integer(1))))
    rubi.add(rule420)

    pattern421 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, a, c, d, x)))
    rule421 = ReplacementRule(pattern421, lambda m, b, a, c, d, x : Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), Mul(b, c)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(b, d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), Mul(b, c)))), m), x)))
    rubi.add(rule421)

    pattern422 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Integer(3), Denominator(m_), Integer(4))), (m, b, c, d, x)))
    rule422 = ReplacementRule(pattern422, lambda m, b, c, d, x : Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(b, c, x), Mul(b, d, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(b, c, x), Mul(b, d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule422)

    pattern423 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (m, a, c, d, x)))
    rule423 = ReplacementRule(pattern423, lambda m, a, c, d, x : Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), c))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(d, Pow(x, Integer(2))), Mul(x, Add(Mul(a, d), c))), m), x)))
    rubi.add(rule423)

    pattern424 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (m, b, a, c, x)))
    rule424 = ReplacementRule(pattern424, lambda m, b, a, c, x : Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), m), Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2))), Mul(x, Add(a, Mul(b, c)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Mul(b, Pow(x, Integer(2))), Mul(x, Add(a, Mul(b, c)))), m), x)))
    rubi.add(rule424)

    pattern425 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), m_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4))), (m, c, d, x)))
    rule425 = ReplacementRule(pattern425, lambda m, c, d, x : Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), m), Pow(Add(Mul(c, x), Mul(d, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(c, x), Mul(d, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule425)

    pattern426 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Integer(3), Denominator(m_), Integer(4))), (m, b, c, x)))
    rule426 = ReplacementRule(pattern426, lambda m, b, c, x : Mul(Pow(Mul(b, x), m), Pow(Add(c, x), m), Pow(Add(Mul(b, c, x), Mul(b, Pow(x, Integer(2)))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(b, c, x), Mul(b, Pow(x, Integer(2)))), m), x)))
    rubi.add(rule426)

    pattern427 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), m_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (m, a, c, x)))
    rule427 = ReplacementRule(pattern427, lambda m, a, c, x : Mul(Pow(Add(a, x), m), Pow(Add(c, x), m), Pow(Add(Mul(a, c), Pow(x, Integer(2)), Mul(x, Add(a, c))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(a, c), Pow(x, Integer(2)), Mul(x, Add(a, c))), m), x)))
    rubi.add(rule427)

    pattern428 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), m_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(m_), Less(Integer(-1), m_, Integer(0)), LessEqual(Integer(3), Denominator(m_), Integer(4))), (m, c, x)))
    rule428 = ReplacementRule(pattern428, lambda m, c, x : Mul(Pow(x, m), Pow(Add(c, x), m), Pow(Add(Mul(c, x), Pow(x, Integer(2))), Mul(Integer(-1), m)), Int(Pow(Add(Mul(c, x), Pow(x, Integer(2))), m), x)))
    rubi.add(rule428)

    pattern429 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(PosQ(Mul(Pow(b_, Integer(-1)), d_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule429 = ReplacementRule(pattern429, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule429)

    pattern430 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(Mul(b_, c_)), PosQ(Mul(Pow(b_, Integer(-1)), d_))), (b, c, d, x)))
    rule430 = ReplacementRule(pattern430, lambda b, c, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule430)

    pattern431 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(PosQ(d_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x)))
    rule431 = ReplacementRule(pattern431, lambda a, c, d, x : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule431)

    pattern432 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), PosQ(Mul(Pow(b_, Integer(-1)), d_))), (b, a, d, x)))
    rule432 = ReplacementRule(pattern432, lambda b, a, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule432)

    pattern433 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(PosQ(Pow(b_, Integer(-1))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule433 = ReplacementRule(pattern433, lambda b, a, c, x : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule433)

    pattern434 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), PosQ(d_)), (c, d, x)))
    rule434 = ReplacementRule(pattern434, lambda c, d, x : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule434)

    pattern435 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PosQ(Mul(Pow(b_, Integer(-1)), d_))), (b, d, x)))
    rule435 = ReplacementRule(pattern435, lambda b, d, x : With(List(Set(q, Rt(Mul(Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule435)

    pattern436 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(And(NonzeroQ(Mul(b_, c_)), PosQ(Pow(b_, Integer(-1)))), (b, c, x)))
    rule436 = ReplacementRule(pattern436, lambda b, c, x : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule436)

    pattern437 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(PosQ(d_), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x)))
    rule437 = ReplacementRule(pattern437, lambda a, d, x : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule437)

    pattern438 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(PosQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x)))
    rule438 = ReplacementRule(pattern438, lambda a, c, x : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule438)

    pattern439 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(NonzeroQ(Mul(Integer(-1), a_)), PosQ(Pow(b_, Integer(-1)))), (b, a, x)))
    rule439 = ReplacementRule(pattern439, lambda b, a, x : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule439)

    pattern440 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), PosQ(d_)), (d, x)))
    rule440 = ReplacementRule(pattern440, lambda d, x : With(List(Set(q, Rt(d, Integer(3)))), Add(Mul(Integer(-1), Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(-1), Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule440)

    pattern441 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), PosQ(Integer(1))), (c, x)))
    rule441 = ReplacementRule(pattern441, lambda c, x : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(Add(c, x))), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(-1)))))))
    rubi.add(rule441)

    pattern442 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Mul(b_, x_), Integer(-1/3))), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), PosQ(Pow(b_, Integer(-1)))), (b, x)))
    rule442 = ReplacementRule(pattern442, lambda b, x : With(List(Set(q, Rt(Pow(b, Integer(-1)), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule442)

    pattern443 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), FreeQ(a, x), cons(And(PosQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_))), (a, x)))
    rule443 = ReplacementRule(pattern443, lambda a, x : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(-1)))))))
    rubi.add(rule443)

    pattern444 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(NonzeroQ(Integer(0)), PosQ(Integer(1))), (x,)))
    rule444 = ReplacementRule(pattern444, lambda x : With(List(Set(q, Rt(Integer(1), Integer(3)))), Add(Mul(Integer(-1), q, ArcTan(Add(Mul(Integer(2), q, Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(-1), Integer(1/2), q, Log(x)), Mul(Integer(-1), Integer(3/2), q, Log(Add(q, Integer(-1)))))))
    rubi.add(rule444)

    pattern445 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(NegQ(Mul(Pow(b_, Integer(-1)), d_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (b, a, c, d, x)))
    rule445 = ReplacementRule(pattern445, lambda b, a, c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule445)

    pattern446 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(Mul(b_, c_)), NegQ(Mul(Pow(b_, Integer(-1)), d_))), (b, c, d, x)))
    rule446 = ReplacementRule(pattern446, lambda b, c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule446)

    pattern447 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(NegQ(d_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (a, c, d, x)))
    rule447 = ReplacementRule(pattern447, lambda a, c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule447)

    pattern448 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Mul(Integer(-1), a_, d_)), NegQ(Mul(Pow(b_, Integer(-1)), d_))), (b, a, d, x)))
    rule448 = ReplacementRule(pattern448, lambda b, a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(1)))))))
    rubi.add(rule448)

    pattern449 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(NegQ(Pow(b_, Integer(-1))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (b, a, c, x)))
    rule449 = ReplacementRule(pattern449, lambda b, a, c, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule449)

    pattern450 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, Mul(d_, x_)), Integer(-2/3))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NegQ(d_), NonzeroQ(c_)), (c, d, x)))
    rule450 = ReplacementRule(pattern450, lambda c, d, x : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Add(c, Mul(d, x)))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, Mul(d, x)), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule450)

    pattern451 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), NegQ(Mul(Pow(b_, Integer(-1)), d_))), (b, d, x)))
    rule451 = ReplacementRule(pattern451, lambda b, d, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1)), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule451)

    pattern452 = Pattern(Int(Mul(Pow(Mul(b_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(b, x), FreeQ(c, x), cons(And(NegQ(Pow(b_, Integer(-1))), NonzeroQ(Mul(b_, c_))), (b, c, x)))
    rule452 = ReplacementRule(pattern452, lambda b, c, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Mul(b, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule452)

    pattern453 = Pattern(Int(Mul(Pow(Mul(d_, x_), Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(d, x), cons(And(NegQ(d_), NonzeroQ(Mul(Integer(-1), a_, d_))), (a, d, x)))
    rule453 = ReplacementRule(pattern453, lambda a, d, x : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(Mul(d, x), Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule453)

    pattern454 = Pattern(Int(Mul(Pow(Add(a_, x_), Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(a, x), FreeQ(c, x), cons(And(NegQ(Integer(1)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (a, c, x)))
    rule454 = ReplacementRule(pattern454, lambda a, c, x : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(Add(a, x), Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule454)

    pattern455 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, Mul(b_, x_)), Integer(-1/3))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(NegQ(Pow(b_, Integer(-1))), NonzeroQ(Mul(Integer(-1), a_))), (b, a, x)))
    rule455 = ReplacementRule(pattern455, lambda b, a, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, Mul(b, x)), Integer(1/3))), Integer(1)))))))
    rubi.add(rule455)

    pattern456 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Mul(d_, x_), Integer(-2/3))), x_), FreeQ(d, x), cons(And(NegQ(d_), NonzeroQ(Integer(0))), (d, x)))
    rule456 = ReplacementRule(pattern456, lambda d, x : With(List(Set(q, Rt(Mul(Integer(-1), d), Integer(3)))), Add(Mul(Pow(d, Integer(-1)), q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), Pow(d, Integer(-1)), q, Log(Mul(d, x))), Mul(Integer(3/2), Pow(d, Integer(-1)), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Mul(d, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule456)

    pattern457 = Pattern(Int(Mul(Pow(x_, Integer(-1/3)), Pow(Add(c_, x_), Integer(-2/3))), x_), FreeQ(c, x), cons(And(NegQ(Integer(1)), NonzeroQ(c_)), (c, x)))
    rule457 = ReplacementRule(pattern457, lambda c, x : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(Add(c, x))), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(1/3)), Pow(Add(c, x), Integer(-1/3))), Integer(1)))))))
    rubi.add(rule457)

    pattern458 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Mul(b_, x_), Integer(-1/3))), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), NegQ(Pow(b_, Integer(-1)))), (b, x)))
    rule458 = ReplacementRule(pattern458, lambda b, x : With(List(Set(q, Rt(Mul(Integer(-1), Pow(b, Integer(-1))), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Mul(b, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule458)

    pattern459 = Pattern(Int(Mul(Pow(x_, Integer(-2/3)), Pow(Add(a_, x_), Integer(-1/3))), x_), FreeQ(a, x), cons(And(NegQ(Integer(1)), NonzeroQ(Mul(Integer(-1), a_))), (a, x)))
    rule459 = ReplacementRule(pattern459, lambda a, x : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3)), Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(Mul(q, Pow(x, Integer(-1/3)), Pow(Add(a, x), Integer(1/3))), Integer(1)))))))
    rubi.add(rule459)

    pattern460 = Pattern(Int(Pow(x_, Integer(-1)), x_), cons(And(NegQ(Integer(1)), NonzeroQ(Integer(0))), (x,)))
    rule460 = ReplacementRule(pattern460, lambda x : With(List(Set(q, Rt(Integer(-1), Integer(3)))), Add(Mul(q, ArcTan(Add(Mul(Integer(-1), Integer(2), q, Pow(Sqrt(Integer(3)), Integer(-1))), Pow(Sqrt(Integer(3)), Integer(-1)))), Sqrt(Integer(3))), Mul(Integer(1/2), q, Log(x)), Mul(Integer(3/2), q, Log(Add(q, Integer(1)))))))
    rubi.add(rule460)

    pattern461 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, a, c, d, x, n)))
    rule461 = ReplacementRule(pattern461, lambda m, b, a, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule461)

    pattern462 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, c, d, x, n)))
    rule462 = ReplacementRule(pattern462, lambda m, b, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule462)

    pattern463 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (m, a, c, d, x, n)))
    rule463 = ReplacementRule(pattern463, lambda m, a, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule463)

    pattern464 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, a, d, x, n)))
    rule464 = ReplacementRule(pattern464, lambda m, b, a, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule464)

    pattern465 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (m, b, a, c, x, n)))
    rule465 = ReplacementRule(pattern465, lambda m, b, a, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule465)

    pattern466 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, c, d, x, n)))
    rule466 = ReplacementRule(pattern466, lambda m, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, Mul(d, x)), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule466)

    pattern467 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, d, x, n)))
    rule467 = ReplacementRule(pattern467, lambda m, b, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), d, Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule467)

    pattern468 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(b_, c_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, c, x, n)))
    rule468 = ReplacementRule(pattern468, lambda m, b, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(Mul(b, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule468)

    pattern469 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_, d_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, a, d, x, n)))
    rule469 = ReplacementRule(pattern469, lambda m, a, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule469)

    pattern470 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (m, a, c, x, n)))
    rule470 = ReplacementRule(pattern470, lambda m, a, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(Add(a, x), Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule470)

    pattern471 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, a, x, n)))
    rule471 = ReplacementRule(pattern471, lambda m, b, a, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1))))))))
    rubi.add(rule471)

    pattern472 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, d, x, n)))
    rule472 = ReplacementRule(pattern472, lambda m, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), d, Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Mul(d, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule472)

    pattern473 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, c, x, n)))
    rule473 = ReplacementRule(pattern473, lambda m, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Pow(p, Integer(-1))), Pow(Add(c, x), Mul(Integer(-1), Pow(p, Integer(-1)))))))))
    rubi.add(rule473)

    pattern474 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, b, x, n)))
    rule474 = ReplacementRule(pattern474, lambda m, b, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(b, Mul(Integer(-1), Pow(x, p))), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Mul(b, x), Pow(p, Integer(-1))))))))
    rubi.add(rule474)

    pattern475 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, a, x, n)))
    rule475 = ReplacementRule(pattern475, lambda m, a, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Mul(Pow(x, Mul(Integer(-1), Pow(p, Integer(-1)))), Pow(Add(a, x), Pow(p, Integer(-1))))))))
    rubi.add(rule475)

    pattern476 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Equal(Add(m_, n_, Integer(1)), Integer(0))), (m, x, n)))
    rule476 = ReplacementRule(pattern476, lambda m, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), Pow(x, p)), Integer(1)), Integer(-1))), x), x, Integer(1)))))
    rubi.add(rule476)

    pattern477 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(a_, b_, c_, d_, m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, a, c, d, x, n)))
    rule477 = ReplacementRule(pattern477, lambda m, b, a, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule477)

    pattern478 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), b_, c_, d_, m_, n_, x_)), (m, b, c, d, x, n)))
    rule478 = ReplacementRule(pattern478, lambda m, b, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule478)

    pattern479 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), IntLinearcQ(a_, Integer(1), c_, d_, m_, n_, x_)), (m, a, c, d, x, n)))
    rule479 = ReplacementRule(pattern479, lambda m, a, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), c, Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule479)

    pattern480 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), IntLinearcQ(a_, b_, Integer(0), d_, m_, n_, x_)), (m, b, a, d, x, n)))
    rule480 = ReplacementRule(pattern480, lambda m, b, a, d, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1)), d), Mul(Pow(b, Integer(-1)), d, Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule480)

    pattern481 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(a_, b_, c_, Integer(1), m_, n_, x_), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (m, b, a, c, x, n)))
    rule481 = ReplacementRule(pattern481, lambda m, b, a, c, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule481)

    pattern482 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), Integer(1), c_, d_, m_, n_, x_)), (m, c, d, x, n)))
    rule482 = ReplacementRule(pattern482, lambda m, c, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(d, Pow(x, p))), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule482)

    pattern483 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), b_, Integer(0), d_, m_, n_, x_)), (m, b, d, x, n)))
    rule483 = ReplacementRule(pattern483, lambda m, b, d, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), d, Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule483)

    pattern484 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(b_, c_)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), b_, c_, Integer(1), m_, n_, x_)), (m, b, c, x, n)))
    rule484 = ReplacementRule(pattern484, lambda m, b, c, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule484)

    pattern485 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Mul(Integer(-1), a_, d_)), IntLinearcQ(a_, Integer(1), Integer(0), d_, m_, n_, x_)), (m, a, d, x, n)))
    rule485 = ReplacementRule(pattern485, lambda m, a, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, d), Mul(d, Pow(x, p))), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule485)

    pattern486 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), IntLinearcQ(a_, Integer(1), c_, Integer(1), m_, n_, x_)), (m, a, c, x, n)))
    rule486 = ReplacementRule(pattern486, lambda m, a, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), c, Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule486)

    pattern487 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(a_, b_, Integer(0), Integer(1), m_, n_, x_)), (m, b, a, x, n)))
    rule487 = ReplacementRule(pattern487, lambda m, b, a, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a, Pow(b, Integer(-1))), Mul(Pow(b, Integer(-1)), Pow(x, p))), n)), x), x, Pow(Add(a, Mul(b, x)), Pow(p, Integer(-1)))))))
    rubi.add(rule487)

    pattern488 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), Integer(1), Integer(0), d_, m_, n_, x_)), (m, d, x, n)))
    rule488 = ReplacementRule(pattern488, lambda m, d, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(d, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule488)

    pattern489 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), cons(And(NonzeroQ(c_), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), Integer(1), c_, Integer(1), m_, n_, x_)), (m, c, x, n)))
    rule489 = ReplacementRule(pattern489, lambda m, c, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(c, Pow(x, p)), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule489)

    pattern490 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), b_, Integer(0), Integer(1), m_, n_, x_)), (m, b, x, n)))
    rule490 = ReplacementRule(pattern490, lambda m, b, x, n : With(List(Set(p, Denominator(m))), Mul(Pow(b, Integer(-1)), p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Mul(Pow(b, Integer(-1)), Pow(x, p)), n)), x), x, Pow(Mul(b, x), Pow(p, Integer(-1)))))))
    rubi.add(rule490)

    pattern491 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), cons(And(RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), NonzeroQ(Mul(Integer(-1), a_)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(a_, Integer(1), Integer(0), Integer(1), m_, n_, x_)), (m, a, x, n)))
    rule491 = ReplacementRule(pattern491, lambda m, a, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Add(Mul(Integer(-1), a), Pow(x, p)), n)), x), x, Pow(Add(a, x), Pow(p, Integer(-1)))))))
    rubi.add(rule491)

    pattern492 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), cons(And(NonzeroQ(Integer(0)), RationalQ(m_, n_), Less(Integer(-1), m_, Integer(0)), Less(Integer(-1), n_, Integer(0)), LessEqual(Denominator(n_), Denominator(m_)), IntLinearcQ(Integer(0), Integer(1), Integer(0), Integer(1), m_, n_, x_)), (m, x, n)))
    rule492 = ReplacementRule(pattern492, lambda m, x, n : With(List(Set(p, Denominator(m))), Mul(p, Subst(Int(Mul(Pow(x, Add(Mul(p, Add(m, Integer(1))), Integer(-1))), Pow(Pow(x, p), n)), x), x, Pow(x, Pow(p, Integer(-1)))))))
    rubi.add(rule492)

    pattern493 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1)))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, b, a, c, d, x, n)))
    rule493 = ReplacementRule(pattern493, lambda m, b, a, c, d, x, n : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))))
    rubi.add(rule493)

    pattern494 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, c, d, x, n)))
    rule494 = ReplacementRule(pattern494, lambda m, b, c, d, x, n : Add(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule494)

    pattern495 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, a, c, d, x, n)))
    rule495 = ReplacementRule(pattern495, lambda m, a, c, d, x, n : Add(Mul(Integer(-1), d, Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))))
    rubi.add(rule495)

    pattern496 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, a, d, x, n)))
    rule496 = ReplacementRule(pattern496, lambda m, b, a, d, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule496)

    pattern497 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, a, c, x, n)))
    rule497 = ReplacementRule(pattern497, lambda m, b, a, c, x, n : Add(Mul(Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule497)

    pattern498 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, c, d, x, n)))
    rule498 = ReplacementRule(pattern498, lambda m, c, d, x, n : Add(Mul(Integer(-1), Pow(c, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, Mul(d, x)), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule498)

    pattern499 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, d, x, n)))
    rule499 = ReplacementRule(pattern499, lambda m, b, d, x, n : Add(Mul(zoo, Pow(b, Integer(-1)), d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(b, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule499)

    pattern500 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(b_, c_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, c, x, n)))
    rule500 = ReplacementRule(pattern500, lambda m, b, c, x, n : Add(Mul(Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(b, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule500)

    pattern501 = Pattern(Int(Mul(Pow(Mul(d_, x_), n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NonzeroQ(Mul(Integer(-1), a_, d_)), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, a, d, x, n)))
    rule501 = ReplacementRule(pattern501, lambda m, a, d, x, n : Add(Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2)))), Mul(Integer(-1), Pow(a, Integer(-1)), Pow(d, Integer(-1)), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule501)

    pattern502 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, a, c, x, n)))
    rule502 = ReplacementRule(pattern502, lambda m, a, c, x, n : Add(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(Add(a, x), Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule502)

    pattern503 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, Mul(b_, x_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, a, x, n)))
    rule503 = ReplacementRule(pattern503, lambda m, b, a, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule503)

    pattern504 = Pattern(Int(Mul(Pow(x_, m_), Pow(Mul(d_, x_), n_)), x_), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, d, x, n)))
    rule504 = ReplacementRule(pattern504, lambda m, d, x, n : Add(Mul(zoo, d, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Mul(d, x), n)), x), Simplify(Add(m, n, Integer(2)))), Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(Mul(d, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)))))
    rubi.add(rule504)

    pattern505 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(c_), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, c, x, n)))
    rule505 = ReplacementRule(pattern505, lambda m, c, x, n : Add(Mul(Pow(c, Integer(-1)), Pow(x, Add(m, Integer(1))), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Integer(-1), Pow(c, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, Simplify(Add(m, Integer(1)))), Pow(Add(c, x), n)), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule505)

    pattern506 = Pattern(Int(Mul(Pow(x_, n_), Pow(Mul(b_, x_), m_)), x_), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, b, x, n)))
    rule506 = ReplacementRule(pattern506, lambda m, b, x, n : Add(Mul(zoo, Pow(b, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(b, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Mul(b, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule506)

    pattern507 = Pattern(Int(Mul(Pow(x_, n_), Pow(Add(a_, x_), m_)), x_), FreeQ(a, x), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Mul(Integer(-1), a_)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, a, x, n)))
    rule507 = ReplacementRule(pattern507, lambda m, a, x, n : Add(Mul(Integer(-1), Pow(a, Integer(-1)), Pow(x, Add(n, Integer(1))), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(Pow(a, Integer(-1)), Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(Add(a, x), Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule507)

    pattern508 = Pattern(Int(Mul(Pow(x_, m_), Pow(x_, n_)), x_), FreeQ(m, x), FreeQ(n, x), cons(And(NonzeroQ(Integer(0)), NonzeroQ(Add(m_, Integer(1))), NegativeIntegerQ(Simplify(Add(m_, n_, Integer(2)))), Or(SumSimplerQ(m_, Integer(1)), Not(SumSimplerQ(n_, Integer(1))))), (m, x, n)))
    rule508 = ReplacementRule(pattern508, lambda m, x, n : Add(Mul(zoo, Pow(x, Add(m, Integer(1))), Pow(x, Add(n, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1))), Mul(zoo, Pow(Add(m, Integer(1)), Integer(-1)), Int(Mul(Pow(x, n), Pow(x, Simplify(Add(m, Integer(1))))), x), Simplify(Add(m, n, Integer(2))))))
    rubi.add(rule508)

    pattern509 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (m, b, c, d, x, n)))
    rule509 = ReplacementRule(pattern509, lambda m, b, c, d, x, n : Mul(Pow(b, Integer(-1)), Pow(c, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), d, x))))
    rubi.add(rule509)

    pattern510 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (m, c, d, x, n)))
    rule510 = ReplacementRule(pattern510, lambda m, c, d, x, n : Mul(Pow(c, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), d, x))))
    rubi.add(rule510)

    pattern511 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1))), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1))))))))), (m, b, c, x, n)))
    rule511 = ReplacementRule(pattern511, lambda m, b, c, x, n : Mul(Pow(b, Integer(-1)), Pow(c, n), Pow(Mul(b, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), x))))
    rubi.add(rule511)

    pattern512 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Or(IntegerQ(n_), And(PositiveQ(c_), Not(And(ZeroQ(Add(n_, Integer(1/2))), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (m, c, x, n)))
    rule512 = ReplacementRule(pattern512, lambda m, c, x, n : Mul(Pow(c, n), Pow(x, Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(c, Integer(-1)), x))))
    rubi.add(rule512)

    pattern513 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)))), (m, b, c, d, x, n)))
    rule513 = ReplacementRule(pattern513, lambda m, b, c, d, x, n : Mul(Pow(d, Integer(-1)), Pow(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1)), d), Mul(Integer(-1), m)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule513)

    pattern514 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)))), (m, c, d, x, n)))
    rule514 = ReplacementRule(pattern514, lambda m, c, d, x, n : Mul(Pow(d, Integer(-1)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d), Mul(Integer(-1), m)), Pow(Add(c, Mul(d, x)), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)))))
    rubi.add(rule514)

    pattern515 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)))))), (m, b, c, x, n)))
    rule515 = ReplacementRule(pattern515, lambda m, b, c, x, n : Mul(Pow(Mul(Integer(-1), Pow(b, Integer(-1)), Pow(c, Integer(-1))), Mul(Integer(-1), m)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule515)

    pattern516 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(n_)), Or(IntegerQ(m_), PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)))))), (m, c, x, n)))
    rule516 = ReplacementRule(pattern516, lambda m, c, x, n : Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1))), Mul(Integer(-1), m)), Pow(Add(c, x), Add(n, Integer(1))), Pow(Add(n, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), m), Add(n, Integer(1)), Add(n, Integer(2)), Add(Integer(1), Mul(Pow(c, Integer(-1)), x)))))
    rubi.add(rule516)

    pattern517 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (m, b, c, d, x, n)))
    rule517 = ReplacementRule(pattern517, lambda m, b, c, d, x, n : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), n)), x)))
    rubi.add(rule517)

    pattern518 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Mul(Integer(-1), Pow(d_, Integer(2)))))))))), (m, c, d, x, n)))
    rule518 = ReplacementRule(pattern518, lambda m, c, d, x, n : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(x, m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), d, x)), n)), x)))
    rubi.add(rule518)

    pattern519 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1))))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (m, b, c, x, n)))
    rule519 = ReplacementRule(pattern519, lambda m, b, c, x, n : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(Mul(b, x), m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), n)), x)))
    rubi.add(rule519)

    pattern520 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1))))), Or(Not(RationalQ(n_)), And(RationalQ(m_), Not(And(ZeroQ(Add(n_, Integer(1/2))), ZeroQ(Add(Pow(c_, Integer(2)), Integer(-1)))))))), (m, c, x, n)))
    rule520 = ReplacementRule(pattern520, lambda m, c, x, n : Mul(Pow(c, IntPart(n)), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(x, m), Pow(Add(Integer(1), Mul(Pow(c, Integer(-1)), x)), n)), x)))
    rubi.add(rule520)

    pattern521 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)), d_)))), (m, b, c, d, x, n)))
    rule521 = ReplacementRule(pattern521, lambda m, b, c, d, x, n : Mul(Pow(Mul(b, x), FracPart(m)), Pow(Mul(Integer(-1), b, c, Pow(d, Integer(-1))), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), m), Pow(Add(c, Mul(d, x)), n)), x)))
    rubi.add(rule521)

    pattern522 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)), d_)))), (m, c, d, x, n)))
    rule522 = ReplacementRule(pattern522, lambda m, c, d, x, n : Mul(Pow(x, FracPart(m)), Pow(Mul(Integer(-1), c, Pow(d, Integer(-1))), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), d, x), m), Pow(Add(c, Mul(d, x)), n)), x)))
    rubi.add(rule522)

    pattern523 = Pattern(Int(Mul(Pow(Mul(b_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(b_, Integer(-1)), Pow(c_, Integer(-1)))))), (m, b, c, x, n)))
    rule523 = ReplacementRule(pattern523, lambda m, b, c, x, n : Mul(Pow(Mul(b, x), FracPart(m)), Pow(Mul(Integer(-1), b, c), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), m), Pow(Add(c, x), n)), x)))
    rubi.add(rule523)

    pattern524 = Pattern(Int(Mul(Pow(x_, m_), Pow(Add(c_, x_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), Not(PositiveQ(c_)), Not(PositiveQ(Mul(Integer(-1), Pow(c_, Integer(-1)))))), (m, c, x, n)))
    rule524 = ReplacementRule(pattern524, lambda m, c, x, n : Mul(Pow(x, FracPart(m)), Pow(Mul(Integer(-1), c), IntPart(m)), Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), Mul(Integer(-1), FracPart(m))), Int(Mul(Pow(Mul(Integer(-1), Pow(c, Integer(-1)), x), m), Pow(Add(c, x), n)), x)))
    rubi.add(rule524)

    pattern525 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)))), (m, a, b, c, d, x, n)))
    rule525 = ReplacementRule(pattern525, lambda m, a, b, c, d, x, n : Mul(Pow(b, Add(Mul(Integer(-1), n), Integer(-1))), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), n), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, Mul(b, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))))
    rubi.add(rule525)

    pattern526 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_))), (m, a, c, d, x, n)))
    rule526 = ReplacementRule(pattern526, lambda m, a, c, d, x, n : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Add(Mul(Integer(-1), a, d), c), n), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, x), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))))
    rubi.add(rule526)

    pattern527 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_)))), (m, a, b, c, x, n)))
    rule527 = ReplacementRule(pattern527, lambda m, a, b, c, x, n : Mul(Pow(b, Add(Mul(Integer(-1), n), Integer(-1))), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), n), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(a, Mul(b, x))))))
    rubi.add(rule527)

    pattern528 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), cons(And(IntegerQ(n_), Not(IntegerQ(m_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_))), (m, a, c, x, n)))
    rule528 = ReplacementRule(pattern528, lambda m, a, c, x, n : Mul(Pow(Add(Mul(Integer(-1), a), c), n), Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(a, x)))))
    rubi.add(rule528)

    pattern529 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), PositiveQ(Mul(b_, Pow(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)), Integer(-1)))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), d_, Pow(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_)), Integer(-1)))))))), (m, a, b, c, d, x, n)))
    rule529 = ReplacementRule(pattern529, lambda m, a, b, c, d, x, n : Mul(Pow(b, Integer(-1)), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), n)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, Mul(b, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))))))
    rubi.add(rule529)

    pattern530 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), PositiveQ(Pow(Add(Mul(Integer(-1), a_, d_), c_), Integer(-1))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), d_, Pow(Add(Mul(Integer(-1), a_, d_), c_), Integer(-1)))))))), (m, a, c, d, x, n)))
    rule530 = ReplacementRule(pattern530, lambda m, a, c, d, x, n : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Mul(Integer(-1), n)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), d, Add(a, x), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))))))
    rubi.add(rule530)

    pattern531 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), PositiveQ(Mul(b_, Pow(Add(Mul(Integer(-1), a_), Mul(b_, c_)), Integer(-1)))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a_), Mul(b_, c_)), Integer(-1)))))))), (m, a, b, c, x, n)))
    rule531 = ReplacementRule(pattern531, lambda m, a, b, c, x, n : Mul(Pow(b, Integer(-1)), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(Integer(-1), n)), Pow(Add(a, Mul(b, x)), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(a, Mul(b, x))))))
    rubi.add(rule531)

    pattern532 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), PositiveQ(Pow(Add(Mul(Integer(-1), a_), c_), Integer(-1))), Or(RationalQ(m_), Not(And(RationalQ(n_), PositiveQ(Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a_), c_), Integer(-1)))))))), (m, a, c, x, n)))
    rule532 = ReplacementRule(pattern532, lambda m, a, c, x, n : Mul(Pow(Add(a, x), Add(m, Integer(1))), Pow(Add(m, Integer(1)), Integer(-1)), Pow(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Mul(Integer(-1), n)), Hypergeometric2F1(Mul(Integer(-1), n), Add(m, Integer(1)), Add(m, Integer(2)), Mul(Integer(-1), Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(a, x)))))
    rubi.add(rule532)

    pattern533 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), Mul(b_, c_))), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (m, a, b, c, d, x, n)))
    rule533 = ReplacementRule(pattern533, lambda m, a, b, c, d, x, n : Mul(Pow(Mul(b, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), IntPart(n))), Pow(Mul(b, Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(Mul(b, c, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1))), Mul(b, d, x, Pow(Add(Mul(Integer(-1), a, d), Mul(b, c)), Integer(-1)))), n)), x)))
    rubi.add(rule533)

    pattern534 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, Mul(d_, x_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_, d_), c_)), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (m, a, c, d, x, n)))
    rule534 = ReplacementRule(pattern534, lambda m, a, c, d, x, n : Mul(Pow(Mul(Add(c, Mul(d, x)), Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(Integer(-1), FracPart(n))), Pow(Add(c, Mul(d, x)), FracPart(n)), Pow(Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)), Mul(Integer(-1), IntPart(n))), Int(Mul(Pow(Add(a, x), m), Pow(Add(Mul(c, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1))), Mul(d, x, Pow(Add(Mul(Integer(-1), a, d), c), Integer(-1)))), n)), x)))
    rubi.add(rule534)

    pattern535 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, x_)), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), Mul(b_, c_))), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (m, a, b, c, x, n)))
    rule535 = ReplacementRule(pattern535, lambda m, a, b, c, x, n : Mul(Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(Integer(-1), IntPart(n))), Pow(Mul(b, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)), Add(c, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(Mul(b, c, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1))), Mul(b, x, Pow(Add(Mul(Integer(-1), a), Mul(b, c)), Integer(-1)))), n)), x)))
    rubi.add(rule535)

    pattern536 = Pattern(Int(Mul(Pow(Add(a_, x_), m_), Pow(Add(c_, x_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(Not(IntegerQ(m_)), Not(IntegerQ(n_)), NonzeroQ(Add(Mul(Integer(-1), a_), c_)), Or(RationalQ(m_), Not(SimplerQ(Add(n_, Integer(1)), Add(m_, Integer(1)))))), (m, a, c, x, n)))
    rule536 = ReplacementRule(pattern536, lambda m, a, c, x, n : Mul(Pow(Mul(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Add(c, x)), Mul(Integer(-1), FracPart(n))), Pow(Add(c, x), FracPart(n)), Pow(Pow(Add(Mul(Integer(-1), a), c), Integer(-1)), Mul(Integer(-1), IntPart(n))), Int(Mul(Pow(Add(a, x), m), Pow(Add(Mul(c, Pow(Add(Mul(Integer(-1), a), c), Integer(-1))), Mul(x, Pow(Add(Mul(Integer(-1), a), c), Integer(-1)))), n)), x)))
    rubi.add(rule536)

    pattern537 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, c, u, d, x, n)))
    rule537 = ReplacementRule(pattern537, lambda m, b, a, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule537)

    pattern538 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, c, u, d, x, n)))
    rule538 = ReplacementRule(pattern538, lambda b, a, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule538)

    pattern539 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, c, u, d, x, n)))
    rule539 = ReplacementRule(pattern539, lambda m, b, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule539)

    pattern540 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, c, u, d, x, n)))
    rule540 = ReplacementRule(pattern540, lambda m, a, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule540)

    pattern541 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Add(c_, Mul(d_, u_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, c, u, d, x)))
    rule541 = ReplacementRule(pattern541, lambda m, b, a, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule541)

    pattern542 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Pow(Add(a_, Mul(b_, u_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, u, d, x, n)))
    rule542 = ReplacementRule(pattern542, lambda m, b, a, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule542)

    pattern543 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Pow(Add(c_, u_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, c, u, x, n)))
    rule543 = ReplacementRule(pattern543, lambda m, b, a, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule543)

    pattern544 = Pattern(Int(Mul(b_, u_, Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, c, u, d, x, n)))
    rule544 = ReplacementRule(pattern544, lambda b, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule544)

    pattern545 = Pattern(Int(Mul(Add(a_, u_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, c, u, d, x, n)))
    rule545 = ReplacementRule(pattern545, lambda a, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule545)

    pattern546 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Add(c_, Mul(d_, u_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, c, u, d, x)))
    rule546 = ReplacementRule(pattern546, lambda b, a, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule546)

    pattern547 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Add(a_, Mul(b_, u_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, u, d, x, n)))
    rule547 = ReplacementRule(pattern547, lambda b, a, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule547)

    pattern548 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Pow(Add(c_, u_), n_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, c, u, x, n)))
    rule548 = ReplacementRule(pattern548, lambda b, a, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule548)

    pattern549 = Pattern(Int(Mul(Pow(u_, m_), Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, c, u, d, x, n)))
    rule549 = ReplacementRule(pattern549, lambda m, c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule549)

    pattern550 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Add(c_, Mul(d_, u_))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, c, u, d, x)))
    rule550 = ReplacementRule(pattern550, lambda m, b, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule550)

    pattern551 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Mul(d_, u_), n_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, u, d, x, n)))
    rule551 = ReplacementRule(pattern551, lambda m, b, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule551)

    pattern552 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Pow(Add(c_, u_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, c, u, x, n)))
    rule552 = ReplacementRule(pattern552, lambda m, b, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule552)

    pattern553 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Add(c_, Mul(d_, u_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, c, u, d, x)))
    rule553 = ReplacementRule(pattern553, lambda m, a, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule553)

    pattern554 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Pow(Add(a_, u_), m_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, u, d, x, n)))
    rule554 = ReplacementRule(pattern554, lambda m, a, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule554)

    pattern555 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Pow(Add(c_, u_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, c, u, x, n)))
    rule555 = ReplacementRule(pattern555, lambda m, a, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule555)

    pattern556 = Pattern(Int(Mul(d_, u_, Pow(Add(a_, Mul(b_, u_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, u, d, x)))
    rule556 = ReplacementRule(pattern556, lambda m, b, a, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule556)

    pattern557 = Pattern(Int(Mul(Pow(Add(a_, Mul(b_, u_)), m_), Add(c_, u_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, c, u, x)))
    rule557 = ReplacementRule(pattern557, lambda m, b, a, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, Mul(b, x)), m), Add(c, x)), x), x, u)))
    rubi.add(rule557)

    pattern558 = Pattern(Int(Mul(Pow(u_, n_), Pow(Add(a_, Mul(b_, u_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, u, x, n)))
    rule558 = ReplacementRule(pattern558, lambda m, b, a, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule558)

    pattern559 = Pattern(Int(Mul(u_, Pow(Add(c_, Mul(d_, u_)), n_)), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (c, u, d, x, n)))
    rule559 = ReplacementRule(pattern559, lambda c, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(c, Mul(d, x)), n)), x), x, u)))
    rubi.add(rule559)

    pattern560 = Pattern(Int(Mul(b_, u_, Add(c_, Mul(d_, u_))), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, c, u, d, x)))
    rule560 = ReplacementRule(pattern560, lambda b, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule560)

    pattern561 = Pattern(Int(Mul(b_, u_, Pow(Mul(d_, u_), n_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, u, d, x, n)))
    rule561 = ReplacementRule(pattern561, lambda b, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule561)

    pattern562 = Pattern(Int(Mul(b_, u_, Pow(Add(c_, u_), n_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, c, u, x, n)))
    rule562 = ReplacementRule(pattern562, lambda b, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule562)

    pattern563 = Pattern(Int(Mul(Add(a_, u_), Add(c_, Mul(d_, u_))), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, c, u, d, x)))
    rule563 = ReplacementRule(pattern563, lambda a, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule563)

    pattern564 = Pattern(Int(Mul(Pow(Mul(d_, u_), n_), Add(a_, u_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, u, d, x, n)))
    rule564 = ReplacementRule(pattern564, lambda a, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(d, x), n), Add(a, x)), x), x, u)))
    rubi.add(rule564)

    pattern565 = Pattern(Int(Mul(Add(a_, u_), Pow(Add(c_, u_), n_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, c, u, x, n)))
    rule565 = ReplacementRule(pattern565, lambda a, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule565)

    pattern566 = Pattern(Int(Mul(d_, u_, Add(a_, Mul(b_, u_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, u, d, x)))
    rule566 = ReplacementRule(pattern566, lambda b, a, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule566)

    pattern567 = Pattern(Int(Mul(Add(a_, Mul(b_, u_)), Add(c_, u_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, c, u, x)))
    rule567 = ReplacementRule(pattern567, lambda b, a, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, Mul(b, x)), Add(c, x)), x), x, u)))
    rubi.add(rule567)

    pattern568 = Pattern(Int(Mul(Pow(u_, n_), Add(a_, Mul(b_, u_))), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, u, x, n)))
    rule568 = ReplacementRule(pattern568, lambda b, a, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule568)

    pattern569 = Pattern(Int(Mul(Pow(u_, m_), Add(c_, Mul(d_, u_))), x_), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, c, u, d, x)))
    rule569 = ReplacementRule(pattern569, lambda m, c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule569)

    pattern570 = Pattern(Int(Mul(Pow(u_, m_), Pow(Mul(d_, u_), n_)), x_), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, u, d, x, n)))
    rule570 = ReplacementRule(pattern570, lambda m, u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule570)

    pattern571 = Pattern(Int(Mul(Pow(u_, m_), Pow(Add(c_, u_), n_)), x_), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, c, u, x, n)))
    rule571 = ReplacementRule(pattern571, lambda m, c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule571)

    pattern572 = Pattern(Int(Mul(d_, u_, Pow(Mul(b_, u_), m_)), x_), FreeQ(b, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, u, d, x)))
    rule572 = ReplacementRule(pattern572, lambda m, b, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule572)

    pattern573 = Pattern(Int(Mul(Pow(Mul(b_, u_), m_), Add(c_, u_)), x_), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, c, u, x)))
    rule573 = ReplacementRule(pattern573, lambda m, b, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Mul(b, x), m), Add(c, x)), x), x, u)))
    rubi.add(rule573)

    pattern574 = Pattern(Int(Mul(Pow(u_, n_), Pow(Mul(b_, u_), m_)), x_), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, u, x, n)))
    rule574 = ReplacementRule(pattern574, lambda m, b, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule574)

    pattern575 = Pattern(Int(Mul(d_, u_, Pow(Add(a_, u_), m_)), x_), FreeQ(a, x), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, u, d, x)))
    rule575 = ReplacementRule(pattern575, lambda m, a, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule575)

    pattern576 = Pattern(Int(Mul(Pow(Add(a_, u_), m_), Add(c_, u_)), x_), FreeQ(a, x), FreeQ(c, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, c, u, x)))
    rule576 = ReplacementRule(pattern576, lambda m, a, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(Add(a, x), m), Add(c, x)), x), x, u)))
    rubi.add(rule576)

    pattern577 = Pattern(Int(Mul(Pow(u_, n_), Pow(Add(a_, u_), m_)), x_), FreeQ(a, x), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, u, x, n)))
    rule577 = ReplacementRule(pattern577, lambda m, a, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule577)

    pattern578 = Pattern(Int(Mul(u_, Pow(Add(a_, Mul(b_, u_)), m_)), x_), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, a, u, x)))
    rule578 = ReplacementRule(pattern578, lambda m, b, a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(a, Mul(b, x)), m)), x), x, u)))
    rubi.add(rule578)

    pattern579 = Pattern(Int(Mul(u_, Add(c_, Mul(d_, u_))), x_), FreeQ(c, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (c, u, d, x)))
    rule579 = ReplacementRule(pattern579, lambda c, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(c, Mul(d, x))), x), x, u)))
    rubi.add(rule579)

    pattern580 = Pattern(Int(Mul(u_, Pow(Mul(d_, u_), n_)), x_), FreeQ(d, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (u, d, x, n)))
    rule580 = ReplacementRule(pattern580, lambda u, d, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Mul(d, x), n)), x), x, u)))
    rubi.add(rule580)

    pattern581 = Pattern(Int(Mul(u_, Pow(Add(c_, u_), n_)), x_), FreeQ(c, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (c, u, x, n)))
    rule581 = ReplacementRule(pattern581, lambda c, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(c, x), n)), x), x, u)))
    rubi.add(rule581)

    pattern582 = Pattern(Int(Mul(b_, d_, Pow(u_, Integer(2))), x_), FreeQ(b, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, u, d, x)))
    rule582 = ReplacementRule(pattern582, lambda b, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, d, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule582)

    pattern583 = Pattern(Int(Mul(b_, u_, Add(c_, u_)), x_), FreeQ(b, x), FreeQ(c, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, c, u, x)))
    rule583 = ReplacementRule(pattern583, lambda b, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Add(c, x)), x), x, u)))
    rubi.add(rule583)

    pattern584 = Pattern(Int(Mul(b_, u_, Pow(u_, n_)), x_), FreeQ(b, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, u, x, n)))
    rule584 = ReplacementRule(pattern584, lambda b, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, x, Pow(x, n)), x), x, u)))
    rubi.add(rule584)

    pattern585 = Pattern(Int(Mul(d_, u_, Add(a_, u_)), x_), FreeQ(a, x), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, u, d, x)))
    rule585 = ReplacementRule(pattern585, lambda a, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Add(a, x)), x), x, u)))
    rubi.add(rule585)

    pattern586 = Pattern(Int(Mul(Add(a_, u_), Add(c_, u_)), x_), FreeQ(a, x), FreeQ(c, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, c, u, x)))
    rule586 = ReplacementRule(pattern586, lambda a, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Add(a, x), Add(c, x)), x), x, u)))
    rubi.add(rule586)

    pattern587 = Pattern(Int(Mul(Pow(u_, n_), Add(a_, u_)), x_), FreeQ(a, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, u, x, n)))
    rule587 = ReplacementRule(pattern587, lambda a, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, n), Add(a, x)), x), x, u)))
    rubi.add(rule587)

    pattern588 = Pattern(Int(Mul(u_, Add(a_, Mul(b_, u_))), x_), FreeQ(a, x), FreeQ(b, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, a, u, x)))
    rule588 = ReplacementRule(pattern588, lambda b, a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(a, Mul(b, x))), x), x, u)))
    rubi.add(rule588)

    pattern589 = Pattern(Int(Mul(d_, u_, Pow(u_, m_)), x_), FreeQ(d, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, u, d, x)))
    rule589 = ReplacementRule(pattern589, lambda m, u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, x, Pow(x, m)), x), x, u)))
    rubi.add(rule589)

    pattern590 = Pattern(Int(Mul(Pow(u_, m_), Add(c_, u_)), x_), FreeQ(c, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, c, u, x)))
    rule590 = ReplacementRule(pattern590, lambda m, c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Add(c, x)), x), x, u)))
    rubi.add(rule590)

    pattern591 = Pattern(Int(Mul(Pow(u_, m_), Pow(u_, n_)), x_), FreeQ(m, x), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, u, x, n)))
    rule591 = ReplacementRule(pattern591, lambda m, u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(Pow(x, m), Pow(x, n)), x), x, u)))
    rubi.add(rule591)

    pattern592 = Pattern(Int(Mul(u_, Pow(Mul(b_, u_), m_)), x_), FreeQ(b, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, b, u, x)))
    rule592 = ReplacementRule(pattern592, lambda m, b, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Mul(b, x), m)), x), x, u)))
    rubi.add(rule592)

    pattern593 = Pattern(Int(Mul(u_, Pow(Add(a_, u_), m_)), x_), FreeQ(a, x), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, a, u, x)))
    rule593 = ReplacementRule(pattern593, lambda m, a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(Add(a, x), m)), x), x, u)))
    rubi.add(rule593)

    pattern594 = Pattern(Int(Mul(d_, Pow(u_, Integer(2))), x_), FreeQ(d, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (u, d, x)))
    rule594 = ReplacementRule(pattern594, lambda u, d, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(d, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule594)

    pattern595 = Pattern(Int(Mul(u_, Add(c_, u_)), x_), FreeQ(c, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (c, u, x)))
    rule595 = ReplacementRule(pattern595, lambda c, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(c, x)), x), x, u)))
    rubi.add(rule595)

    pattern596 = Pattern(Int(Mul(u_, Pow(u_, n_)), x_), FreeQ(n, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (u, x, n)))
    rule596 = ReplacementRule(pattern596, lambda u, x, n : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(x, n)), x), x, u)))
    rubi.add(rule596)

    pattern597 = Pattern(Int(Mul(b_, Pow(u_, Integer(2))), x_), FreeQ(b, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (b, u, x)))
    rule597 = ReplacementRule(pattern597, lambda b, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(b, Pow(x, Integer(2))), x), x, u)))
    rubi.add(rule597)

    pattern598 = Pattern(Int(Mul(u_, Add(a_, u_)), x_), FreeQ(a, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (a, u, x)))
    rule598 = ReplacementRule(pattern598, lambda a, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Add(a, x)), x), x, u)))
    rubi.add(rule598)

    pattern599 = Pattern(Int(Mul(u_, Pow(u_, m_)), x_), FreeQ(m, x), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (m, u, x)))
    rule599 = ReplacementRule(pattern599, lambda m, u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Mul(x, Pow(x, m)), x), x, u)))
    rubi.add(rule599)

    pattern600 = Pattern(Int(Pow(u_, Integer(2)), x_), cons(And(LinearQ(u_, x_), NonzeroQ(Coefficient(u_, x_, Integer(0)))), (u, x)))
    rule600 = ReplacementRule(pattern600, lambda u, x : Mul(Pow(Coefficient(u, x, Integer(1)), Integer(-1)), Subst(Int(Pow(x, Integer(2)), x), x, u)))
    rubi.add(rule600)

    return rubi
