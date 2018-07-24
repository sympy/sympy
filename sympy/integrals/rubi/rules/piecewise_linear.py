from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Module, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec, ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma, FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

def piecewise_linear(rubi):

    pattern1 = Pattern(Integral(u_**WC('m', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: PiecewiseLinearQ(u, x)), )
    def With1(x, m, u):
        c = D(u, x)
        return Subst(Int(x**m, x), x, u)/c
    rule1 = ReplacementRule(pattern1, lambda x, m, u : With1(x, m, u))
    rubi.add(rule1)

    pattern2 = Pattern(Integral(v_/u_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, x, u, a: NonzeroQ(-a*v + b*u)))
    def With2(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return b*x/a - (-a*v + b*u)*Int(1/u, x)/a
    rule2 = ReplacementRule(pattern2, lambda x, v, u : With2(x, v, u))
    rubi.add(rule2)

    pattern3 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda n: Unequal(n, S(1))), CustomConstraint(lambda v, b, n, x, u, a: NonzeroQ(-a*v + b*u)))
    def With3(x, n, v, u):
        a = D(u, x)
        b = D(v, x)
        return -(-a*v + b*u)*Int(v**(n + S(-1))/u, x)/a + v**n/(a*n)
    rule3 = ReplacementRule(pattern3, lambda x, n, v, u : With3(x, n, v, u))
    rubi.add(rule3)

    pattern4 = Pattern(Integral(S(1)/(u_*v_), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, x, u, a: NonzeroQ(-a*v + b*u)))
    def With4(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return -a*Int(1/u, x)/(-a*v + b*u) + b*Int(1/v, x)/(-a*v + b*u)
    rule4 = ReplacementRule(pattern4, lambda x, v, u : With4(x, v, u))
    rubi.add(rule4)

    pattern5 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, u, a: NonzeroQ(-a*v + b*u) & PosQ((-a*v + b*u)/a)))
    def With5(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return S(2)*ArcTan(sqrt(v)/Rt((-a*v + b*u)/a, S(2)))/(a*Rt((-a*v + b*u)/a, S(2)))
    rule5 = ReplacementRule(pattern5, lambda x, v, u : With5(x, v, u))
    rubi.add(rule5)

    pattern6 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, u, a: NonzeroQ(-a*v + b*u) & NegQ((-a*v + b*u)/a)))
    def With6(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return -S(2)*atanh(sqrt(v)/Rt((a*v - b*u)/a, S(2)))/(a*Rt((a*v - b*u)/a, S(2)))
    rule6 = ReplacementRule(pattern6, lambda x, v, u : With6(x, v, u))
    rubi.add(rule6)

    pattern7 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Less(n, S(-1))), CustomConstraint(lambda v, b, n, x, u, a: NonzeroQ(-a*v + b*u)))
    def With7(x, n, v, u):
        a = D(u, x)
        b = D(v, x)
        return -a*Int(v**(n + S(1))/u, x)/(-a*v + b*u) + v**(n + S(1))/((n + S(1))*(-a*v + b*u))
    rule7 = ReplacementRule(pattern7, lambda x, n, v, u : With7(x, n, v, u))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(v_**n_/u_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda v, b, n, u, a: NonzeroQ(-a*v + b*u)))
    def With8(x, n, v, u):
        a = D(u, x)
        b = D(v, x)
        return v**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), -a*v/(-a*v + b*u))/((n + S(1))*(-a*v + b*u))
    rule8 = ReplacementRule(pattern8, lambda x, n, v, u : With8(x, n, v, u))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, u, a: PosQ(a*b) & NonzeroQ(-a*v + b*u)))
    def With9(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return S(2)*atanh(sqrt(u)*Rt(a*b, S(2))/(a*sqrt(v)))/Rt(a*b, S(2))
    rule9 = ReplacementRule(pattern9, lambda x, v, u : With9(x, v, u))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda v, b, u, a: NegQ(a*b) & NonzeroQ(-a*v + b*u)))
    def With10(x, v, u):
        a = D(u, x)
        b = D(v, x)
        return S(2)*ArcTan(sqrt(u)*Rt(-a*b, S(2))/(a*sqrt(v)))/Rt(-a*b, S(2))
    rule10 = ReplacementRule(pattern10, lambda x, v, u : With10(x, v, u))
    rubi.add(rule10)

    pattern11 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n, m: ZeroQ(m + n + S(2))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda v, b, n, m, u, a: NonzeroQ(-a*v + b*u)))
    def With11(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return -u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u))
    rule11 = ReplacementRule(pattern11, lambda x, n, v, u, m : With11(x, n, v, u, m))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda n, x: FreeQ(n, x)), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: NonzeroQ(m + S(1))), CustomConstraint(lambda n, m: (NegativeIntegerQ(m) & Not(IntegerQ(n))) | (PositiveIntegerQ(n) & Not(IntegerQ(m))) | (LessEqual(n, m) & PositiveIntegerQ(n, m)) | (Greater(n, S(0)) & Less(m, S(-1)) & RationalQ(m, n) & Not(IntegerQ(m + n) & Less(m + n + S(2), S(0)) & (FractionQ(m) | GreaterEqual(m + S(2)*n + S(1), S(0)))))), CustomConstraint(lambda v, b, n, x, m, u, a: NonzeroQ(-a*v + b*u)))
    def With12(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return -b*n*Int(u**(m + S(1))*v**(n + S(-1)), x)/(a*(m + S(1))) + u**(m + S(1))*v**n/(a*(m + S(1)))
    rule12 = ReplacementRule(pattern12, lambda x, n, v, u, m : With12(x, n, v, u, m))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(2))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(1))), CustomConstraint(lambda n, m: Not(PositiveIntegerQ(m) & (Not(IntegerQ(n)) | Less(S(0), m, n)))), CustomConstraint(lambda n, m: Not(IntegerQ(m + n) & Less(m + n + S(2), S(0)))), CustomConstraint(lambda v, b, n, x, m, u, a: NonzeroQ(-a*v + b*u)))
    def With13(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return -n*(-a*v + b*u)*Int(u**m*v**(n + S(-1)), x)/(a*(m + n + S(1))) + u**(m + S(1))*v**n/(a*(m + n + S(1)))
    rule13 = ReplacementRule(pattern13, lambda x, n, v, u, m : With13(x, n, v, u, m))
    rubi.add(rule13)

    pattern14 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(1))), CustomConstraint(lambda n: Not(RationalQ(n))), CustomConstraint(lambda n: SumSimplerQ(n, S(-1))), CustomConstraint(lambda v, b, n, x, m, u, a: NonzeroQ(-a*v + b*u)))
    def With14(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return -n*(-a*v + b*u)*Int(u**m*v**(n + S(-1)), x)/(a*(m + n + S(1))) + u**(m + S(1))*v**n/(a*(m + n + S(1)))
    rule14 = ReplacementRule(pattern14, lambda x, n, v, u, m : With14(x, n, v, u, m))
    rubi.add(rule14)

    pattern15 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda n, m: NonzeroQ(m + n + S(2))), CustomConstraint(lambda m: RationalQ(m)), CustomConstraint(lambda m: Less(m, S(-1))), CustomConstraint(lambda v, b, n, x, m, u, a: NonzeroQ(-a*v + b*u)))
    def With15(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return b*(m + n + S(2))*Int(u**(m + S(1))*v**n, x)/((m + S(1))*(-a*v + b*u)) - u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u))
    rule15 = ReplacementRule(pattern15, lambda x, n, v, u, m : With15(x, n, v, u, m))
    rubi.add(rule15)

    pattern16 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: Not(RationalQ(m))), CustomConstraint(lambda m: SumSimplerQ(m, S(1))), CustomConstraint(lambda v, b, n, x, m, u, a: NonzeroQ(-a*v + b*u)))
    def With16(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return b*(m + n + S(2))*Int(u**(m + S(1))*v**n, x)/((m + S(1))*(-a*v + b*u)) - u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u))
    rule16 = ReplacementRule(pattern16, lambda x, n, v, u, m : With16(x, n, v, u, m))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(u_**m_*v_**n_, x_), CustomConstraint(lambda x, v, u: PiecewiseLinearQ(u, v, x)), CustomConstraint(lambda m: Not(IntegerQ(m))), CustomConstraint(lambda n: Not(IntegerQ(n))), CustomConstraint(lambda v, b, n, m, u, a: NonzeroQ(-a*v + b*u)))
    def With17(x, n, v, u, m):
        a = D(u, x)
        b = D(v, x)
        return u**m*v**(n + S(1))*(b*u/(-a*v + b*u))**(-m)*Hypergeometric2F1(-m, n + S(1), n + S(2), -a*v/(-a*v + b*u))/(b*(n + S(1)))
    rule17 = ReplacementRule(pattern17, lambda x, n, v, u, m : With17(x, n, v, u, m))
    rubi.add(rule17)

    pattern18 = Pattern(Integral(u_**WC('n', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda x, u: PiecewiseLinearQ(u, x)), CustomConstraint(lambda x, u: Not(LinearQ(u, x))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), )
    def With18(x, n, b, u, a):
        c = D(u, x)
        return -Int(u**n, x) - c*n*Int(u**(n + S(-1))*(a + b*x)*log(a + b*x), x)/b + u**n*(a + b*x)*log(a + b*x)/b
    rule18 = ReplacementRule(pattern18, lambda x, n, b, u, a : With18(x, n, b, u, a))
    rubi.add(rule18)

    pattern19 = Pattern(Integral(u_**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), CustomConstraint(lambda a, x: FreeQ(a, x)), CustomConstraint(lambda b, x: FreeQ(b, x)), CustomConstraint(lambda m, x: FreeQ(m, x)), CustomConstraint(lambda x, u: PiecewiseLinearQ(u, x)), CustomConstraint(lambda x, u: Not(LinearQ(u, x))), CustomConstraint(lambda n: RationalQ(n)), CustomConstraint(lambda n: Greater(n, S(0))), CustomConstraint(lambda m: NonzeroQ(m + S(1))), )
    def With19(x, n, b, u, a, m):
        c = D(u, x)
        return -Int(u**n*(a + b*x)**m, x)/(m + S(1)) - c*n*Int(u**(n + S(-1))*(a + b*x)**(m + S(1))*log(a + b*x), x)/(b*(m + S(1))) + u**n*(a + b*x)**(m + S(1))*log(a + b*x)/(b*(m + S(1)))
    rule19 = ReplacementRule(pattern19, lambda x, n, b, u, a, m : With19(x, n, b, u, a, m))
    rubi.add(rule19)

    return rubi
