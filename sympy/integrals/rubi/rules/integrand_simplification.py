'''
This code is automatically generated. Never edit it manually.
For details of generating the code see `rubi_parsing_guide.md` in `parsetools`.
'''

from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint, is_match
    from sympy.integrals.rubi.utility_function import (
        Int, Sum, Set, With, Module, Scan, MapAnd, FalseQ,
        ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ,
        PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ,
        ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ,
        NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart,
        FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest,
        SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient,
        Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart,
        IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan,
        ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec,
        ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less,
        Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ,
        PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ,
        ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ,
        Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ,
        SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator,
        NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ,
        InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ,
        EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree,
        PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts,
        TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ,
        NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll,
        ExpandLinearProduct, GCD, ContentFactor, NumericFactor,
        NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst,
        ExpandExpression, Apart, SmartApart, MatchQ,
        PolynomialQuotientRemainder, FreeFactors, NonfreeFactors,
        RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms,
        ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup,
        AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor,
        RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon,
        MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ,
        GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList,
        PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ,
        RationalFunctionFactors, NonrationalFunctionFactors, Reverse,
        RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand,
        SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree,
        CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree,
        GeneralizedBinomialParts, GeneralizedTrinomialDegree,
        GeneralizedTrinomialParts, MonomialQ, MonomialSumQ,
        MinimumMonomialExponent, MonomialExponent, LinearMatchQ,
        PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ,
        TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ,
        QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms,
        NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial,
        PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD,
        AlgebraicFunctionFactors, NonalgebraicFunctionFactors,
        QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ,
        Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors,
        NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop,
        CombineExponents, FactorInteger, FactorAbsurdNumber,
        SubstForInverseFunction, SubstForFractionalPower,
        SubstForFractionalPowerOfQuotientOfLinears,
        FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ,
        SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ,
        FractionalPowerSubexpressionQ, Apply, FactorNumericGcd,
        MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ,
        TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest,
        OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors,
        PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn,
        PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree,
        FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify,
        FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand,
        NormalizeIntegrandAux, NormalizeIntegrandFactor,
        NormalizeIntegrandFactorBase, NormalizeTogether,
        NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors,
        SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm,
        TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum,
        UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear,
        PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ,
        IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor,
        FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ,
        FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator,
        SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand,
        SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM,
        SubstForFractionalPowerOfLinear, FractionalPowerOfLinear,
        InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig,
        FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ,
        PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ,
        KnownTangentIntegrandQ, KnownCotangentIntegrandQ,
        KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst,
        AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand,
        ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp,
        ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ,
        FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ,
        PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ,
        FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ,
        FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ,
        FunctionOfLog, PowerVariableExpn, PowerVariableDegree,
        PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic,
        SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ,
        Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ,
        SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2,
        ConstantFactor, SameQ, ReplacePart, CommonFactors,
        MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential,
        FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux,
        FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev,
        rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent,
        RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct,
        SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma,
        FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ,
        _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify,
        _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum,
        _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux,
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist, Sum_doit, PolynomialQuotient, Floor,
        PolynomialRemainder, Factor, PolyLog, CosIntegral, SinIntegral, LogIntegral, SinhIntegral,
        CoshIntegral, Rule, Erf, PolyGamma, ExpIntegralEi, ExpIntegralE, LogGamma , UtilityOperator, Factorial,
        Zeta, ProductLog, DerivativeDivides, HypergeometricPFQ, IntHide, OneQ, Null, rubi_exp as exp, rubi_log as log, Discriminant,
        Negative, Quotient
    )
    from sympy import (Integral, S, sqrt, And, Or, Integer, Float, Mod, I, Abs, simplify, Mul,
    Add, Pow, sign, EulerGamma)
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec, atan2)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii , Pqq, Q, R, r, C, k, u = symbols('i ii Pqq Q R r C k u')
    _UseGamma = False
    ShowSteps = False
    StepCounter = None

def integrand_simplification(rubi, matcher):
    from sympy.integrals.rubi.constraints import cons1, cons2, cons3, cons4, cons5, cons6, cons7, cons8, cons9, cons10, cons11, cons12, cons13, cons14, cons15, cons16, cons17, cons18, cons19, cons20, cons21, cons22, cons23, cons24, cons25, cons26, cons27, cons28, cons29, cons30, cons31, cons32, cons33, cons34, cons35, cons36, cons37, cons38, cons39, cons40, cons41, cons42, cons43, cons44, cons45, cons46, cons47, cons48, cons49, cons50, cons51, cons52, cons53, cons54, cons55, cons56, cons57, cons58, cons59, cons60, cons61, cons_with_33, cons_with_34, cons62, cons63, cons64, cons65, cons_with_35, cons_with_36

    pattern1 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1)
    def replacement1(u, a, x, p, n, b):
        rubi.append(1)
        return Int(u*(b*x**n)**p, x)
    rule1 = ReplacementRule(pattern1, replacement1)

    matcher.add(pattern1, 1)
    pattern2 = Pattern(Integral((a_ + x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons1)
    def replacement2(u, c, a, x, p, j, n, b):
        rubi.append(2)
        return Int(u*(b*x**n + c*x**(S(2)*n))**p, x)
    rule2 = ReplacementRule(pattern2, replacement2)

    matcher.add(pattern2, 2)
    pattern3 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons8)
    def replacement3(u, c, a, x, p, j, n, b):
        rubi.append(3)
        return Int(u*(a + c*x**(S(2)*n))**p, x)
    rule3 = ReplacementRule(pattern3, replacement3)

    matcher.add(pattern3, 3)
    pattern4 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons9)
    def replacement4(u, c, a, x, p, j, n, b):
        rubi.append(4)
        return Int(u*(a + b*x**n)**p, x)
    rule4 = ReplacementRule(pattern4, replacement4)

    matcher.add(pattern4, 4)
    pattern5 = Pattern(Integral((v_*WC('a', S(1)) + v_*WC('b', S(1)) + WC('w', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons10)
    def replacement5(u, a, v, x, p, w, b):
        rubi.append(5)
        return Int(u*(v*(a + b) + w)**p, x)
    rule5 = ReplacementRule(pattern5, replacement5)

    matcher.add(pattern5, 5)
    pattern6 = Pattern(Integral(Pm_**p_*WC('u', S(1)), x_), cons11, cons12, cons13)
    def replacement6(Pm, u, x, p):
        rubi.append(6)
        return Int(Pm**p*u, x)
    rule6 = ReplacementRule(pattern6, replacement6)

    matcher.add(pattern6, 6)
    pattern7 = Pattern(Integral(a_, x_), cons2, cons2)
    def replacement7(x, a):
        rubi.append(7)
        return Simp(a*x, x)
    rule7 = ReplacementRule(pattern7, replacement7)

    matcher.add(pattern7, 7)
    pattern8 = Pattern(Integral(a_*(b_ + x_*WC('c', S(1))), x_), cons2, cons3, cons7, cons14)
    def replacement8(x, c, a, b):
        rubi.append(8)
        return Simp(a*(b + c*x)**S(2)/(S(2)*c), x)
    rule8 = ReplacementRule(pattern8, replacement8)

    matcher.add(pattern8, 8)
    pattern9 = Pattern(Integral(-u_, x_))
    def replacement9(x, u):
        rubi.append(9)
        return Dist(S(-1), Int(u, x), x)
    rule9 = ReplacementRule(pattern9, replacement9)

    matcher.add(pattern9, 9)
    pattern10 = Pattern(Integral(u_*Complex(S(0), a_), x_), cons2, cons15)
    def replacement10(x, u, a):
        rubi.append(10)
        return Dist(Complex(S(0), a), Int(u, x), x)
    rule10 = ReplacementRule(pattern10, replacement10)

    matcher.add(pattern10, 10)
    pattern11 = Pattern(Integral(a_*u_, x_), cons2, cons2)
    def replacement11(x, u, a):
        rubi.append(11)
        return Dist(a, Int(u, x), x)
    rule11 = ReplacementRule(pattern11, replacement11)

    matcher.add(pattern11, 11)
    pattern12 = Pattern(Integral(u_, x_), cons16)
    def replacement12(x, u):
        rubi.append(12)
        return Simp(IntSum(u, x), x)
    rule12 = ReplacementRule(pattern12, replacement12)

    matcher.add(pattern12, 12)
    pattern13 = Pattern(Integral(v_**WC('m', S(1))*(b_*v_)**n_*WC('u', S(1)), x_), cons3, cons4, cons17)
    def replacement13(u, m, v, x, n, b):
        rubi.append(13)
        return Dist(b**(-m), Int(u*(b*v)**(m + n), x), x)
    rule13 = ReplacementRule(pattern13, replacement13)

    matcher.add(pattern13, 13)
    pattern14 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons18, cons19, cons20)
    def replacement14(u, m, a, v, x, n, b):
        rubi.append(14)
        return Dist(a**(m + S(1)/2)*b**(n + S(-1)/2)*sqrt(b*v)/sqrt(a*v), Int(u*v**(m + n), x), x)
    rule14 = ReplacementRule(pattern14, replacement14)

    matcher.add(pattern14, 14)
    pattern15 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons18, cons22, cons20)
    def replacement15(u, m, a, v, x, n, b):
        rubi.append(15)
        return Dist(a**(m + S(-1)/2)*b**(n + S(1)/2)*sqrt(a*v)/sqrt(b*v), Int(u*v**(m + n), x), x)
    rule15 = ReplacementRule(pattern15, replacement15)

    matcher.add(pattern15, 15)
    pattern16 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons18, cons23, cons20)
    def replacement16(u, m, a, v, x, n, b):
        rubi.append(16)
        return Dist(a**(m + n)*(a*v)**(-n)*(b*v)**n, Int(u*v**(m + n), x), x)
    rule16 = ReplacementRule(pattern16, replacement16)

    matcher.add(pattern16, 16)
    pattern17 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons18, cons23, cons24)
    def replacement17(u, m, a, v, x, n, b):
        rubi.append(17)
        return Dist(a**(-IntPart(n))*b**IntPart(n)*(a*v)**(-FracPart(n))*(b*v)**FracPart(n), Int(u*(a*v)**(m + n), x), x)
    rule17 = ReplacementRule(pattern17, replacement17)

    matcher.add(pattern17, 17)
    pattern18 = Pattern(Integral((a_ + v_*WC('b', S(1)))**WC('m', S(1))*(c_ + v_*WC('d', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons25, cons17, cons26)
    def replacement18(u, c, m, a, v, x, b, n, d):
        rubi.append(18)
        return Dist((b/d)**m, Int(u*(c + d*v)**(m + n), x), x)
    rule18 = ReplacementRule(pattern18, replacement18)

    matcher.add(pattern18, 18)
    pattern19 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons25, cons28, cons29)
    def replacement19(u, c, a, m, v, x, b, n, d):
        rubi.append(19)
        return Dist((b/d)**m, Int(u*(c + d*v)**(m + n), x), x)
    rule19 = ReplacementRule(pattern19, replacement19)

    matcher.add(pattern19, 19)
    pattern20 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons25, cons30)
    def replacement20(u, c, a, m, v, x, b, n, d):
        rubi.append(20)
        return Dist((a + b*v)**m*(c + d*v)**(-m), Int(u*(c + d*v)**(m + n), x), x)
    rule20 = ReplacementRule(pattern20, replacement20)

    matcher.add(pattern20, 20)
    pattern21 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_**S(2)*WC('c', S(1)) + v_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons7, cons31, cons32)
    def replacement21(u, c, a, m, v, x, b):
        rubi.append(21)
        return Dist(S(1)/a, Int(u*(a*v)**(m + S(1))*(b + c*v), x), x)
    rule21 = ReplacementRule(pattern21, replacement21)

    matcher.add(pattern21, 21)
    pattern22 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(v_**S(2)*WC('C', S(1)) + v_*WC('B', S(1)) + WC('A', S(0)))*WC('u', S(1)), x_), cons2, cons3, cons34, cons35, cons36, cons33, cons31, cons32)
    def replacement22(u, C, B, a, m, v, x, A, b):
        rubi.append(22)
        return Dist(b**(S(-2)), Int(u*(a + b*v)**(m + S(1))*Simp(B*b - C*a + C*b*v, x), x), x)
    rule22 = ReplacementRule(pattern22, replacement22)

    matcher.add(pattern22, 22)
    pattern23 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**WC('q', S(1))*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons37, cons38, cons39, cons40)
    def replacement23(u, q, c, m, a, x, b, p, n, d):
        rubi.append(23)
        return Dist((d/a)**p, Int(u*x**(-n*p)*(a + b*x**n)**(m + p), x), x)
    rule23 = ReplacementRule(pattern23, replacement23)

    matcher.add(pattern23, 23)
    pattern24 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**j_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons6, cons41, cons42, cons43, cons44)
    def replacement24(u, c, m, a, x, b, p, j, n, d):
        rubi.append(24)
        return Dist((-b**S(2)/d)**m, Int(u*(a - b*x**n)**(-m), x), x)
    rule24 = ReplacementRule(pattern24, replacement24)

    matcher.add(pattern24, 24)
    pattern25 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons45, cons38)
    def replacement25(u, c, a, x, p, b):
        rubi.append(25)
        return Int(S(2)**(-S(2)*p)*c**(-p)*u*(b + S(2)*c*x)**(S(2)*p), x)
    rule25 = ReplacementRule(pattern25, replacement25)

    matcher.add(pattern25, 25)
    pattern26 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38)
    def replacement26(u, c, a, x, p, n, n2, b):
        rubi.append(26)
        return Dist(c**(-p), Int(u*(b/S(2) + c*x**n)**(S(2)*p), x), x)
    rule26 = ReplacementRule(pattern26, replacement26)

    matcher.add(pattern26, 26)
    pattern27 = Pattern(Integral((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons47)
    def replacement27(e, c, a, d, x, p, b):
        rubi.append(27)
        return Dist(d/b, Subst(Int(x**p, x), x, a + b*x + c*x**S(2)), x)
    rule27 = ReplacementRule(pattern27, replacement27)

    matcher.add(pattern27, 27)
    pattern28 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons5, cons50, cons17, cons49)
    def replacement28(u, q, a, m, x, p, b):
        rubi.append(28)
        return Int(u*x**(m*p)*(a + b*x**(-p + q))**m, x)
    rule28 = ReplacementRule(pattern28, replacement28)

    matcher.add(pattern28, 28)
    pattern29 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons5, cons50, cons52, cons17, cons49, cons51)
    def replacement29(u, q, c, a, m, x, p, r, b):
        rubi.append(29)
        return Int(u*x**(m*p)*(a + b*x**(-p + q) + c*x**(-p + r))**m, x)
    rule29 = ReplacementRule(pattern29, replacement29)

    matcher.add(pattern29, 29)
    pattern30 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons53)
    def replacement30(a, m, x, n, b):
        rubi.append(30)
        return Simp(rubi_log(RemoveContent(a + b*x**n, x))/(b*n), x)
    rule30 = ReplacementRule(pattern30, replacement30)

    matcher.add(pattern30, 30)
    pattern31 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons53, cons54)
    def replacement31(a, m, x, p, n, b):
        rubi.append(31)
        return Simp((a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
    rule31 = ReplacementRule(pattern31, replacement31)

    matcher.add(pattern31, 31)
    pattern32 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons5, cons55, cons56, cons54)
    def replacement32(b1, m, a1, x, p, n, a2, b2):
        rubi.append(32)
        return Simp((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))), x)
    rule32 = ReplacementRule(pattern32, replacement32)

    matcher.add(pattern32, 32)

    pattern33 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons11, cons61, cons_with_33)
    def replacement33(a, x, p, n, Pm, Qm, b):

        m = Expon(Pm, x)
        rubi.append(33)
        return Dist(Coeff(Qm, x, m + S(-1))/(m*Coeff(Pm, x, m)), Subst(Int((a + b*x**n)**p, x), x, Pm), x)
    rule33 = ReplacementRule(pattern33, replacement33)

    matcher.add(pattern33, 33)

    pattern34 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + Pm_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons11, cons61, cons_with_34)
    def replacement34(c, a, x, p, n, Pm, n2, Qm, b):

        m = Expon(Pm, x)
        rubi.append(34)
        return Dist(Coeff(Qm, x, m + S(-1))/(m*Coeff(Pm, x, m)), Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, Pm), x)
    rule34 = ReplacementRule(pattern34, replacement34)

    matcher.add(pattern34, 34)

    pattern35 = Pattern(Integral(Pq_**m_*Qr_**p_*WC('u', S(1)), x_), cons62, cons63, cons64, cons65, cons_with_35)
    def replacement35(u, m, Qr, x, p, Pq):

        gcd = PolyGCD(Pq, Qr, x)
        rubi.append(35)
        return Int(gcd**(m + p)*u*PolynomialQuotient(Pq, gcd, x)**m*PolynomialQuotient(Qr, gcd, x)**p, x)
    rule35 = ReplacementRule(pattern35, replacement35)

    matcher.add(pattern35, 35)

    pattern36 = Pattern(Integral(Pq_*Qr_**p_*WC('u', S(1)), x_), cons63, cons64, cons65, cons_with_36)
    def replacement36(u, Qr, x, p, Pq):

        gcd = PolyGCD(Pq, Qr, x)
        rubi.append(36)
        return Int(gcd**(p + S(1))*u*PolynomialQuotient(Pq, gcd, x)*PolynomialQuotient(Qr, gcd, x)**p, x)
    rule36 = ReplacementRule(pattern36, replacement36)

    matcher.add(pattern36, 36)
    return matcher, [rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9, rule10, rule11, rule12, rule13, rule14, rule15, rule16, rule17, rule18, rule19, rule20, rule21, rule22, rule23, rule24, rule25, rule26, rule27, rule28, rule29, rule30, rule31, rule32, rule33, rule34, rule35, rule36, ]
