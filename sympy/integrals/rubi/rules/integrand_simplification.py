from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (
        sympy_op_factory, Int, Sum, Set, With, Module, Scan, MapAnd, FalseQ,
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
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist, PolynomialQuotient
    )
    from sympy import Integral, S, sqrt, And, Or
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

def integrand_simplification(rubi):
    from sympy.integrals.rubi.constraints import cons1, cons2, cons3, cons4, cons5, cons6, cons7, cons8, cons9, cons10, cons11, cons12, cons13, cons14, cons15, cons16, cons17, cons18, cons19, cons20, cons21, cons22, cons23, cons24, cons25, cons26, cons27, cons28, cons29, cons30, cons31, cons32, cons33, cons34, cons35, cons36, cons37, cons38, cons39, cons40, cons41, cons42, cons43, cons44, cons45, cons46, cons47, cons48, cons49, cons50, cons51, cons52, cons53, cons54, cons55, cons56, cons57, cons58, cons59, cons60, cons61, cons62, cons63, cons64, cons65

    pattern1 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1)
    def replacement1(a, x, b, p, u, n):
        rubi.append(1)
        return Int(u*(b*x**n)**p, x)
    rule1 = ReplacementRule(pattern1, replacement1)
    pattern2 = Pattern(Integral((a_ + x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons1)
    def replacement2(a, x, c, b, p, u, j, n):
        rubi.append(2)
        return Int(u*(b*x**n + c*x**(S(2)*n))**p, x)
    rule2 = ReplacementRule(pattern2, replacement2)
    pattern3 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons8)
    def replacement3(a, x, c, b, p, u, j, n):
        rubi.append(3)
        return Int(u*(a + c*x**(S(2)*n))**p, x)
    rule3 = ReplacementRule(pattern3, replacement3)
    pattern4 = Pattern(Integral((x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons6, cons9)
    def replacement4(a, x, c, b, p, u, j, n):
        rubi.append(4)
        return Int(u*(a + b*x**n)**p, x)
    rule4 = ReplacementRule(pattern4, replacement4)
    pattern5 = Pattern(Integral((v_*WC('a', S(1)) + v_*WC('b', S(1)) + WC('w', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons10)
    def replacement5(a, x, v, b, p, u, w):
        rubi.append(5)
        return Int(u*(v*(a + b) + w)**p, x)
    rule5 = ReplacementRule(pattern5, replacement5)
    pattern6 = Pattern(Integral(Pm_**p_*WC('u', S(1)), x_), cons5, cons11, cons12, cons13)
    def replacement6(x, Pm, u, p):
        rubi.append(6)
        return Int(Pm**p*u, x)
    rule6 = ReplacementRule(pattern6, replacement6)
    pattern7 = Pattern(Integral(a_, x_), cons2, cons2)
    def replacement7(a, x):
        rubi.append(7)
        return a*x
    rule7 = ReplacementRule(pattern7, replacement7)
    pattern8 = Pattern(Integral(a_*(b_ + x_*WC('c', S(1))), x_), cons2, cons3, cons7, cons14)
    def replacement8(a, b, x, c):
        rubi.append(8)
        return a*(b + c*x)**S(2)/(S(2)*c)
    rule8 = ReplacementRule(pattern8, replacement8)
    pattern9 = Pattern(Integral(-u_, x_))
    def replacement9(x, u):
        rubi.append(9)
        return -Int(u, x)
    rule9 = ReplacementRule(pattern9, replacement9)
    pattern10 = Pattern(Integral(u_*Complex(S(0), a_), x_), cons2, cons15)
    def replacement10(a, x, u):
        rubi.append(10)
        return Complex(S(0), a)*Int(u, x)
    rule10 = ReplacementRule(pattern10, replacement10)
    pattern11 = Pattern(Integral(a_*u_, x_), cons2, cons2)
    def replacement11(a, x, u):
        rubi.append(11)
        return a*Int(u, x)
    rule11 = ReplacementRule(pattern11, replacement11)
    pattern12 = Pattern(Integral(u_, x_), cons16)
    def replacement12(x, u):
        rubi.append(12)
        return IntSum(u, x)
    rule12 = ReplacementRule(pattern12, replacement12)
    pattern13 = Pattern(Integral(v_**WC('m', S(1))*(b_*v_)**n_*WC('u', S(1)), x_), cons3, cons4, cons17)
    def replacement13(x, v, b, u, n, m):
        rubi.append(13)
        return b**(-m)*Int(u*(b*v)**(m + n), x)
    rule13 = ReplacementRule(pattern13, replacement13)
    pattern14 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons18, cons19, cons20)
    def replacement14(a, x, v, b, u, n, m):
        rubi.append(14)
        return a**(m + S(1)/2)*b**(n + S(-1)/2)*sqrt(b*v)*Int(u*v**(m + n), x)/sqrt(a*v)
    rule14 = ReplacementRule(pattern14, replacement14)
    pattern15 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons18, cons22, cons20)
    def replacement15(a, x, v, b, u, n, m):
        rubi.append(15)
        return a**(m + S(-1)/2)*b**(n + S(1)/2)*sqrt(a*v)*Int(u*v**(m + n), x)/sqrt(b*v)
    rule15 = ReplacementRule(pattern15, replacement15)
    pattern16 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons18, cons23, cons20)
    def replacement16(a, x, v, b, u, n, m):
        rubi.append(16)
        return a**(m + n)*(a*v)**(-n)*(b*v)**n*Int(u*v**(m + n), x)
    rule16 = ReplacementRule(pattern16, replacement16)
    pattern17 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons18, cons23, cons24)
    def replacement17(a, x, v, b, u, n, m):
        rubi.append(17)
        return a**(-IntPart(n))*b**IntPart(n)*(a*v)**(-FracPart(n))*(b*v)**FracPart(n)*Int(u*(a*v)**(m + n), x)
    rule17 = ReplacementRule(pattern17, replacement17)
    pattern18 = Pattern(Integral((a_ + v_*WC('b', S(1)))**WC('m', S(1))*(c_ + v_*WC('d', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons25, cons17, cons26)
    def replacement18(a, x, c, v, b, u, d, n, m):
        rubi.append(18)
        return (b/d)**m*Int(u*(c + d*v)**(m + n), x)
    rule18 = ReplacementRule(pattern18, replacement18)
    pattern19 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons25, cons28, cons29)
    def replacement19(a, x, c, v, b, u, d, n, m):
        rubi.append(19)
        return (b/d)**m*Int(u*(c + d*v)**(m + n), x)
    rule19 = ReplacementRule(pattern19, replacement19)
    pattern20 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(c_ + v_*WC('d', S(1)))**n_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons25, cons30)
    def replacement20(a, x, c, v, b, u, d, n, m):
        rubi.append(20)
        return (a + b*v)**m*(c + d*v)**(-m)*Int(u*(c + d*v)**(m + n), x)
    rule20 = ReplacementRule(pattern20, replacement20)
    pattern21 = Pattern(Integral((v_*WC('a', S(1)))**m_*(v_**S(2)*WC('c', S(1)) + v_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons7, cons31, cons32)
    def replacement21(a, x, c, v, b, u, m):
        rubi.append(21)
        return Int(u*(a*v)**(m + S(1))*(b + c*v), x)/a
    rule21 = ReplacementRule(pattern21, replacement21)
    pattern22 = Pattern(Integral((a_ + v_*WC('b', S(1)))**m_*(v_**S(2)*WC('C', S(1)) + v_*WC('B', S(1)) + WC('A', S(0)))*WC('u', S(1)), x_), cons2, cons3, cons34, cons35, cons36, cons33, cons31, cons32)
    def replacement22(a, x, A, v, b, C, u, B, m):
        rubi.append(22)
        return Int(u*(a + b*v)**(m + S(1))*Simp(B*b - C*a + C*b*v, x), x)/b**S(2)
    rule22 = ReplacementRule(pattern22, replacement22)
    pattern23 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**WC('q', S(1))*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons37, cons38, cons39, cons40)
    def replacement23(a, x, c, b, p, q, u, d, n, m):
        rubi.append(23)
        return (d/a)**p*Int(u*x**(-n*p)*(a + b*x**n)**(m + p), x)
    rule23 = ReplacementRule(pattern23, replacement23)
    pattern24 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('m', S(1))*(c_ + x_**j_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons6, cons41, cons42, cons43, cons44)
    def replacement24(a, x, c, b, p, u, j, d, n, m):
        rubi.append(24)
        return (-b**S(2)/d)**m*Int(u*(a - b*x**n)**(-m), x)
    rule24 = ReplacementRule(pattern24, replacement24)
    pattern25 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons45, cons38)
    def replacement25(a, x, c, b, p, u):
        rubi.append(25)
        return Int(S(2)**(-S(2)*p)*c**(-p)*u*(b + S(2)*c*x)**(S(2)*p), x)
    rule25 = ReplacementRule(pattern25, replacement25)
    pattern26 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38)
    def replacement26(a, x, c, b, p, u, n2, n):
        rubi.append(26)
        return c**(-p)*Int(u*(b/S(2) + c*x**n)**(S(2)*p), x)
    rule26 = ReplacementRule(pattern26, replacement26)
    pattern27 = Pattern(Integral((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons47)
    def replacement27(a, x, c, e, b, p, d):
        rubi.append(27)
        return d*Subst(Int(x**p, x), x, a + b*x + c*x**S(2))/b
    rule27 = ReplacementRule(pattern27, replacement27)
    pattern28 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons5, cons50, cons17, cons49)
    def replacement28(a, x, b, p, q, u, m):
        rubi.append(28)
        return Int(u*x**(m*p)*(a + b*x**(-p + q))**m, x)
    rule28 = ReplacementRule(pattern28, replacement28)
    pattern29 = Pattern(Integral((x_**WC('p', S(1))*WC('a', S(1)) + x_**WC('q', S(1))*WC('b', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons5, cons50, cons52, cons17, cons49, cons51)
    def replacement29(a, x, c, b, p, q, u, r, m):
        rubi.append(29)
        return Int(u*x**(m*p)*(a + b*x**(-p + q) + c*x**(-p + r))**m, x)
    rule29 = ReplacementRule(pattern29, replacement29)
    pattern30 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons53)
    def replacement30(a, x, b, n, m):
        rubi.append(30)
        return log(RemoveContent(a + b*x**n, x))/(b*n)
    rule30 = ReplacementRule(pattern30, replacement30)
    pattern31 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons53, cons54)
    def replacement31(a, x, b, p, n, m):
        rubi.append(31)
        return (a + b*x**n)**(p + S(1))/(b*n*(p + S(1)))
    rule31 = ReplacementRule(pattern31, replacement31)
    pattern32 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), cons57, cons58, cons59, cons60, cons21, cons4, cons5, cons55, cons56, cons54)
    def replacement32(x, p, a1, b1, a2, n, b2, m):
        rubi.append(32)
        return (a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1)))
    rule32 = ReplacementRule(pattern32, replacement32)
    def With33(a, x, b, Pm, Qm, p, n):
        m = Expon(Pm, x)
        if And(Equal(Expon(Qm, x), m + S(-1)), ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x))):
            return True
        return False
    pattern33 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons11, cons61, CustomConstraint(With33))
    def replacement33(a, x, b, Pm, Qm, p, n):
        
        m = Expon(Pm, x)
        rubi.append(33)
        return Coeff(Qm, x, m + S(-1))*Subst(Int((a + b*x**n)**p, x), x, Pm)/(m*Coeff(Pm, x, m))
    rule33 = ReplacementRule(pattern33, replacement33)
    def With34(a, x, c, b, Pm, Qm, p, n2, n):
        m = Expon(Pm, x)
        if And(Equal(Expon(Qm, x), m + S(-1)), ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x))):
            return True
        return False
    pattern34 = Pattern(Integral(Qm_*(Pm_**WC('n', S(1))*WC('b', S(1)) + Pm_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons11, cons61, CustomConstraint(With34))
    def replacement34(a, x, c, b, Pm, Qm, p, n2, n):
        
        m = Expon(Pm, x)
        rubi.append(34)
        return Coeff(Qm, x, m + S(-1))*Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, Pm)/(m*Coeff(Pm, x, m))
    rule34 = ReplacementRule(pattern34, replacement34)
    def With35(x, Pq, p, Qr, u, m):
        gcd = PolyGCD(Pq, Qr, x)
        if NonzeroQ(gcd + S(-1)):
            return True
        return False
    pattern35 = Pattern(Integral(Pq_**m_*Qr_**p_*WC('u', S(1)), x_), cons62, cons63, cons64, cons65, CustomConstraint(With35))
    def replacement35(x, Pq, p, Qr, u, m):
        
        gcd = PolyGCD(Pq, Qr, x)
        rubi.append(35)
        return Int(gcd**(m + p)*u*PolynomialQuotient(Pq, gcd, x)**m*PolynomialQuotient(Qr, gcd, x)**p, x)
    rule35 = ReplacementRule(pattern35, replacement35)
    def With36(x, Pq, p, Qr, u):
        gcd = PolyGCD(Pq, Qr, x)
        if NonzeroQ(gcd + S(-1)):
            return True
        return False
    pattern36 = Pattern(Integral(Pq_*Qr_**p_*WC('u', S(1)), x_), cons63, cons64, cons65, CustomConstraint(With36))
    def replacement36(x, Pq, p, Qr, u):
        
        gcd = PolyGCD(Pq, Qr, x)
        rubi.append(36)
        return Int(gcd**(p + S(1))*u*PolynomialQuotient(Pq, gcd, x)*PolynomialQuotient(Qr, gcd, x)**p, x)
    rule36 = ReplacementRule(pattern36, replacement36)
    return[rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9, rule10, rule11, rule12, rule13, rule14, rule15, rule16, rule17, rule18, rule19, rule20, rule21, rule22, rule23, rule24, rule25, rule26, rule27, rule28, rule29, rule30, rule31, rule32, rule33, rule34, rule35, rule36, ]
