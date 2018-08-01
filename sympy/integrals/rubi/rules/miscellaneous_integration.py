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
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist, Sum_doit, PolynomialQuotient, Floor,
        PolynomialRemainder, Factor, PolyLog, CosIntegral, SinIntegral, LogIntegral, SinhIntegral,
        CoshIntegral, Rule, Erf, PolyGamma, ExpIntegralEi, ExpIntegralE, LogGamma , UtilityOperator, Factorial,
        Zeta, ProductLog, DerivativeDivides, HypergeometricPFQ, IntHide, OneQ, Null, exp, log, Discriminant,
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

def miscellaneous_integration(rubi, matcher, load_rule):
    from sympy.integrals.rubi.constraints import cons147, cons2002, cons2, cons3, cons7, cons4, cons5, cons386, cons27, cons50, cons2003, cons2004, cons2005, cons2006, cons48, cons125, cons208, cons34, cons35, cons36, cons1099, cons2007, cons_with_6938, cons_with_6939, cons66, cons21, cons_with_6940, cons_with_6941, cons_with_6942, cons84, cons1037, cons1036, cons38, cons2008, cons10, cons2009, cons_with_6946, cons2010, cons_with_6947, cons2011, cons209, cons_with_6948, cons1831, cons_with_6949, cons1244, cons_with_6950, cons2012, cons_with_6951, cons46, cons_with_6952, cons_with_6953, cons_with_6954, cons_with_6955, cons_with_6956, cons_with_6957, cons_with_6958, cons_with_6959, cons_with_6960, cons2013, cons_with_6961, cons_with_6962, cons2014, cons2015, cons2016, cons52, cons_with_6963, cons2017, cons800, cons_with_6964, cons2018, cons17, cons_with_6965, cons2019, cons586, cons_with_6966, cons2020, cons2021, cons2022, cons_with_6967, cons2023, cons_with_6968, cons2024, cons_with_6970, cons_with_6971, cons2025, cons2026, cons2027, cons2028, cons667, cons196, cons2029, cons840, cons2030, cons18, cons2031, cons_with_6978, cons148, cons_with_6979, cons45, cons2032, cons_with_6982, cons1854, cons_with_6984, cons1247, cons_with_6985, cons261, cons2033, cons_with_6986, cons367, cons2034, cons_with_6988, cons67, cons1479, cons744, cons1482, cons165, cons2035, cons_with_6993, cons_with_6994, cons2036, cons1676, cons1255, cons2037, cons347, cons_with_6997

    pattern6931 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons147, cons2002)
    def replacement6931(p, u, x, b, c, a, n):
        rubi.append(6931)
        return Dist(c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p)), Int(u*(a + b*x)**(n*p), x), x)
    rule6931 = ReplacementRule(pattern6931, replacement6931)

    if load_rule:
        matcher.add(pattern6931, 6931)
    pattern6932 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons147, cons386)
    def replacement6932(d, q, p, u, x, b, c, a):
        rubi.append(6932)
        return Dist((c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q), Int(u*(a + b*x)**(p*q), x), x)
    rule6932 = ReplacementRule(pattern6932, replacement6932)

    if load_rule:
        matcher.add(pattern6932, 6932)
    pattern6933 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons147, cons386)
    def replacement6933(d, q, p, u, x, b, c, a, n):
        rubi.append(6933)
        return Dist((c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q), Int(u*(a + b*x)**(n*p*q), x), x)
    rule6933 = ReplacementRule(pattern6933, replacement6933)

    if load_rule:
        matcher.add(pattern6933, 6933)
    pattern6934 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1099, cons2003, cons2004, cons2005, cons2006)
    def replacement6934(B, d, F, c, f, e, A, C, x, n, b, a, g):
        rubi.append(6934)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule6934 = ReplacementRule(pattern6934, replacement6934)

    if load_rule:
        matcher.add(pattern6934, 6934)
    pattern6935 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1099, cons2003, cons2007)
    def replacement6935(F, c, e, A, x, b, n, C, a, g):
        rubi.append(6935)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule6935 = ReplacementRule(pattern6935, replacement6935)

    if load_rule:
        matcher.add(pattern6935, 6935)
    pattern6936 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1099, cons2003, cons2004, cons2005, cons2006)
    def replacement6936(d, B, F, c, f, e, A, C, x, n, b, a, g):
        rubi.append(6936)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule6936 = ReplacementRule(pattern6936, replacement6936)

    if load_rule:
        matcher.add(pattern6936, 6936)
    pattern6937 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1099, cons2003, cons2007)
    def replacement6937(F, c, e, A, C, x, n, b, a, g):
        rubi.append(6937)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule6937 = ReplacementRule(pattern6937, replacement6937)

    if load_rule:
        matcher.add(pattern6937, 6937)

    pattern6938 = Pattern(Integral(u_/y_, x_), cons_with_6938)
    def replacement6938(y, x, u):

        q = DerivativeDivides(y, u, x)
        rubi.append(6938)
        return Simp(q*log(RemoveContent(y, x)), x)
    rule6938 = ReplacementRule(pattern6938, replacement6938)

    if load_rule:
        matcher.add(pattern6938, 6938)

    pattern6939 = Pattern(Integral(u_/(w_*y_), x_), cons_with_6939)
    def replacement6939(y, x, w, u):

        q = DerivativeDivides(w*y, u, x)
        rubi.append(6939)
        return Simp(q*log(RemoveContent(w*y, x)), x)
    rule6939 = ReplacementRule(pattern6939, replacement6939)

    if load_rule:
        matcher.add(pattern6939, 6939)

    pattern6940 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons21, cons66, cons_with_6940)
    def replacement6940(y, m, u, x):

        q = DerivativeDivides(y, u, x)
        rubi.append(6940)
        return Simp(q*y**(m + S(1))/(m + S(1)), x)
    rule6940 = ReplacementRule(pattern6940, replacement6940)

    if load_rule:
        matcher.add(pattern6940, 6940)

    pattern6941 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons21, cons4, cons66, cons_with_6941)
    def replacement6941(z, u, y, x, m, n):

        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        rubi.append(6941)
        return Simp(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), x)
    rule6941 = ReplacementRule(pattern6941, replacement6941)

    if load_rule:
        matcher.add(pattern6941, 6941)

    pattern6942 = Pattern(Integral(u_, x_), cons_with_6942)
    def replacement6942(x, u):

        v = SimplifyIntegrand(u, x)
        rubi.append(6942)
        return Int(v, x)
    rule6942 = ReplacementRule(pattern6942, replacement6942)

    if load_rule:
        matcher.add(pattern6942, 6942)
    pattern6943 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1037)
    def replacement6943(d, c, f, e, u, m, x, b, a, n):
        rubi.append(6943)
        return Dist((a*e**S(2) - c*f**S(2))**m, Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule6943 = ReplacementRule(pattern6943, replacement6943)

    if load_rule:
        matcher.add(pattern6943, 6943)
    pattern6944 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1036)
    def replacement6944(d, c, f, e, u, m, x, b, a, n):
        rubi.append(6944)
        return Dist((b*e**S(2) - d*f**S(2))**m, Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule6944 = ReplacementRule(pattern6944, replacement6944)

    if load_rule:
        matcher.add(pattern6944, 6944)
    pattern6945 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), cons2, cons21, cons4, cons38, cons2008, cons10)
    def replacement6945(v, p, u, m, x, w, a, n):
        rubi.append(6945)
        return Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x)
    rule6945 = ReplacementRule(pattern6945, replacement6945)

    if load_rule:
        matcher.add(pattern6945, 6945)

    pattern6946 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons2009, cons_with_6946)
    def replacement6946(d, v, u, m, y, x, b, c, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6946)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), x)
    rule6946 = ReplacementRule(pattern6946, replacement6946)

    if load_rule:
        matcher.add(pattern6946, 6946)

    pattern6947 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons2009, cons2010, cons_with_6947)
    def replacement6947(d, v, p, f, e, u, m, y, x, w, b, c, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6947)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), x)
    rule6947 = ReplacementRule(pattern6947, replacement6947)

    if load_rule:
        matcher.add(pattern6947, 6947)

    pattern6948 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons2009, cons2010, cons2011, cons_with_6948)
    def replacement6948(d, q, p, v, z, f, e, h, g, u, m, y, x, w, b, c, a, n):

        r = DerivativeDivides(y, u, x)
        rubi.append(6948)
        return Dist(r, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), x)
    rule6948 = ReplacementRule(pattern6948, replacement6948)

    if load_rule:
        matcher.add(pattern6948, 6948)

    pattern6949 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1831, cons_with_6949)
    def replacement6949(u, y, x, b, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6949)
        return Dist(a, Int(u, x), x) + Dist(b*q, Subst(Int(x**n, x), x, y), x)
    rule6949 = ReplacementRule(pattern6949, replacement6949)

    if load_rule:
        matcher.add(pattern6949, 6949)

    pattern6950 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1244, cons_with_6950)
    def replacement6950(p, u, y, x, b, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6950)
        return Dist(q, Subst(Int((a + b*x**n)**p, x), x, y), x)
    rule6950 = ReplacementRule(pattern6950, replacement6950)

    if load_rule:
        matcher.add(pattern6950, 6950)

    pattern6951 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons2012, cons_with_6951)
    def replacement6951(v, p, u, m, y, x, b, a, n):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(6951)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n)**p, x), x, y), x)
    rule6951 = ReplacementRule(pattern6951, replacement6951)

    if load_rule:
        matcher.add(pattern6951, 6951)

    pattern6952 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons2009, cons_with_6952)
    def replacement6952(v, p, u, n2, y, x, b, c, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6952)
        return Dist(q, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule6952 = ReplacementRule(pattern6952, replacement6952)

    if load_rule:
        matcher.add(pattern6952, 6952)

    pattern6953 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons46, cons2009, cons2010, cons_with_6953)
    def replacement6953(B, v, p, c, u, n2, A, y, w, x, b, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6953)
        return Dist(q, Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule6953 = ReplacementRule(pattern6953, replacement6953)

    if load_rule:
        matcher.add(pattern6953, 6953)

    pattern6954 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons4, cons5, cons46, cons2010, cons_with_6954)
    def replacement6954(B, p, u, n2, A, y, w, x, c, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6954)
        return Dist(q, Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule6954 = ReplacementRule(pattern6954, replacement6954)

    if load_rule:
        matcher.add(pattern6954, 6954)

    pattern6955 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons2010, cons_with_6955)
    def replacement6955(v, p, u, m, n2, y, w, x, b, c, a, n):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(6955)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule6955 = ReplacementRule(pattern6955, replacement6955)

    if load_rule:
        matcher.add(pattern6955, 6955)

    pattern6956 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons2009, cons2010, cons_with_6956)
    def replacement6956(B, v, p, c, z, u, n2, m, A, w, y, x, b, a, n):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(6956)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule6956 = ReplacementRule(pattern6956, replacement6956)

    if load_rule:
        matcher.add(pattern6956, 6956)

    pattern6957 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons2010, cons_with_6957)
    def replacement6957(B, z, p, u, n2, m, A, w, y, x, c, a, n):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(6957)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule6957 = ReplacementRule(pattern6957, replacement6957)

    if load_rule:
        matcher.add(pattern6957, 6957)

    pattern6958 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons2009, cons_with_6958)
    def replacement6958(d, v, p, u, m, y, x, b, c, a, n):

        q = DerivativeDivides(y, u, x)
        rubi.append(6958)
        return Dist(q, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), x)
    rule6958 = ReplacementRule(pattern6958, replacement6958)

    if load_rule:
        matcher.add(pattern6958, 6958)

    pattern6959 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons2009, cons2010, cons_with_6959)
    def replacement6959(d, q, p, v, f, e, u, m, y, x, w, b, c, a, n):

        r = DerivativeDivides(y, u, x)
        rubi.append(6959)
        return Dist(r, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), x)
    rule6959 = ReplacementRule(pattern6959, replacement6959)

    if load_rule:
        matcher.add(pattern6959, 6959)

    pattern6960 = Pattern(Integral(F_**v_*u_, x_), cons1099, cons1099, cons_with_6960)
    def replacement6960(x, v, F, u):

        q = DerivativeDivides(v, u, x)
        rubi.append(6960)
        return Simp(F**v*q/log(F), x)
    rule6960 = ReplacementRule(pattern6960, replacement6960)

    if load_rule:
        matcher.add(pattern6960, 6960)

    pattern6961 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), cons1099, cons21, cons2013, cons_with_6961)
    def replacement6961(v, F, u, x, w, m):

        q = DerivativeDivides(v, u, x)
        rubi.append(6961)
        return Dist(q, Subst(Int(F**x*x**m, x), x, v), x)
    rule6961 = ReplacementRule(pattern6961, replacement6961)

    if load_rule:
        matcher.add(pattern6961, 6961)

    pattern6962 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons38, cons_with_6962)
    def replacement6962(v, p, u, m, x, w, b, a):

        c = u/(v*D(w, x) + w*D(v, x))
        rubi.append(6962)
        return Dist(c, Subst(Int((a + b*x**p)**m, x), x, v*w), x)
    rule6962 = ReplacementRule(pattern6962, replacement6962)

    if load_rule:
        matcher.add(pattern6962, 6962)

    pattern6963 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons2014, cons2015, cons2016, cons_with_6963)
    def replacement6963(q, p, v, u, x, w, r, b, a, m):

        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(6963)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w), x)
    rule6963 = ReplacementRule(pattern6963, replacement6963)

    if load_rule:
        matcher.add(pattern6963, 6963)

    pattern6964 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons2017, cons2015, cons2016, cons_with_6964)
    def replacement6964(q, p, v, u, x, w, r, b, a, m, s):

        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(6964)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1))), x)
    rule6964 = ReplacementRule(pattern6964, replacement6964)

    if load_rule:
        matcher.add(pattern6964, 6964)

    pattern6965 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons2018, cons38, cons17, cons_with_6965)
    def replacement6965(q, p, v, u, m, x, w, b, a):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(6965)
        return Dist(c*p, Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))), x)
    rule6965 = ReplacementRule(pattern6965, replacement6965)

    if load_rule:
        matcher.add(pattern6965, 6965)

    pattern6966 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons2019, cons586, cons17, cons_with_6966)
    def replacement6966(q, p, v, u, m, x, w, r, b, a):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(6966)
        return -Dist(c*q, Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w), x)
    rule6966 = ReplacementRule(pattern6966, replacement6966)

    if load_rule:
        matcher.add(pattern6966, 6966)

    pattern6967 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons800, cons2020, cons2021, cons2022, cons17, cons_with_6967)
    def replacement6967(q, p, v, u, m, x, w, b, a, s):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(6967)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1))), x)
    rule6967 = ReplacementRule(pattern6967, replacement6967)

    if load_rule:
        matcher.add(pattern6967, 6967)

    pattern6968 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons2023, cons2021, cons2022, cons17, cons_with_6968)
    def replacement6968(q, p, v, u, m, x, w, r, b, a, s):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(6968)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1))), x)
    rule6968 = ReplacementRule(pattern6968, replacement6968)

    if load_rule:
        matcher.add(pattern6968, 6968)
    pattern6969 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons21, cons66, cons2024)
    def replacement6969(m, u, x):
        rubi.append(6969)
        return Dist(S(1)/(m + S(1)), Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1))), x)
    rule6969 = ReplacementRule(pattern6969, replacement6969)

    if load_rule:
        matcher.add(pattern6969, 6969)

    pattern6970 = Pattern(Integral(u_, x_), cons_with_6970)
    def replacement6970(x, u):

        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(6970)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule6970 = ReplacementRule(pattern6970, replacement6970)

    if load_rule:
        matcher.add(pattern6970, 6970)

    pattern6971 = Pattern(Integral(u_, x_), cons_with_6971)
    def replacement6971(x, u):

        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        rubi.append(6971)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule6971 = ReplacementRule(pattern6971, replacement6971)

    if load_rule:
        matcher.add(pattern6971, 6971)
    pattern6972 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons50, cons147, cons10, cons2025, cons2026)
    def replacement6972(q, v, z, p, u, x, w, a, m, n):
        rubi.append(6972)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p), Int(u*v**(m*p)*w**(n*p)*z**(p*q), x), x)
    rule6972 = ReplacementRule(pattern6972, replacement6972)

    if load_rule:
        matcher.add(pattern6972, 6972)
    pattern6973 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons147, cons10, cons2025)
    def replacement6973(v, p, u, x, w, a, m, n):
        rubi.append(6973)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p), Int(u*v**(m*p)*w**(n*p), x), x)
    rule6973 = ReplacementRule(pattern6973, replacement6973)

    if load_rule:
        matcher.add(pattern6973, 6973)
    pattern6974 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons5, cons147, cons10, cons2027, cons2028)
    def replacement6974(v, p, u, x, a, m):
        rubi.append(6974)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p), Int(u*v**(m*p), x), x)
    rule6974 = ReplacementRule(pattern6974, replacement6974)

    if load_rule:
        matcher.add(pattern6974, 6974)
    pattern6975 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons667, cons196, cons2029)
    def replacement6975(p, u, x, b, a, n):
        rubi.append(6975)
        return Dist(FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b)), Int(u*x**(n*p)*(a*x**(-n) + b)**p, x), x)
    rule6975 = ReplacementRule(pattern6975, replacement6975)

    if load_rule:
        matcher.add(pattern6975, 6975)
    pattern6976 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons147, cons196, cons840, cons2030)
    def replacement6976(v, p, u, x, b, a, n):
        rubi.append(6976)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b)**p, x), x)
    rule6976 = ReplacementRule(pattern6976, replacement6976)

    if load_rule:
        matcher.add(pattern6976, 6976)
    pattern6977 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons21, cons5, cons147, cons196, cons840)
    def replacement6977(v, p, u, m, x, b, a, n):
        rubi.append(6977)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x), x)
    rule6977 = ReplacementRule(pattern6977, replacement6977)

    if load_rule:
        matcher.add(pattern6977, 6977)

    pattern6978 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons21, cons52, cons800, cons18, cons2031, cons_with_6978)
    def replacement6978(u, m, x, r, b, a, s):

        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        rubi.append(6978)
        return Dist(v, Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), x)
    rule6978 = ReplacementRule(pattern6978, replacement6978)

    if load_rule:
        matcher.add(pattern6978, 6978)

    pattern6979 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons148, cons_with_6979)
    def replacement6979(u, x, b, a, n):

        v = RationalFunctionExpand(u/(a + b*x**n), x)
        rubi.append(6979)
        return Int(v, x)
    rule6979 = ReplacementRule(pattern6979, replacement6979)

    if load_rule:
        matcher.add(pattern6979, 6979)
    pattern6980 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38, cons2032)
    def replacement6980(p, c, u, n2, x, b, a, n):
        rubi.append(6980)
        return Dist(S(4)**(-p)*c**(-p), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule6980 = ReplacementRule(pattern6980, replacement6980)

    if load_rule:
        matcher.add(pattern6980, 6980)
    pattern6981 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons2032)
    def replacement6981(p, c, u, n2, x, b, a, n):
        rubi.append(6981)
        return Dist((b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p, Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule6981 = ReplacementRule(pattern6981, replacement6981)

    if load_rule:
        matcher.add(pattern6981, 6981)

    pattern6982 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons46, cons148, cons_with_6982)
    def replacement6982(c, u, n2, x, b, a, n):

        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        rubi.append(6982)
        return Int(v, x)
    rule6982 = ReplacementRule(pattern6982, replacement6982)

    if load_rule:
        matcher.add(pattern6982, 6982)
    pattern6983 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons1854)
    def replacement6983(c, u, x, b, a, m, n):
        rubi.append(6983)
        return Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x)
    rule6983 = ReplacementRule(pattern6983, replacement6983)

    if load_rule:
        matcher.add(pattern6983, 6983)

    pattern6984 = Pattern(Integral(u_, x_), cons_with_6984)
    def replacement6984(x, u):

        lst = FunctionOfLinear(u, x)
        rubi.append(6984)
        return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)
    rule6984 = ReplacementRule(pattern6984, replacement6984)

    if load_rule:
        matcher.add(pattern6984, 6984)

    pattern6985 = Pattern(Integral(u_/x_, x_), cons1247, cons2029, cons_with_6985)
    def replacement6985(x, u):

        lst = PowerVariableExpn(u, S(0), x)
        rubi.append(6985)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule6985 = ReplacementRule(pattern6985, replacement6985)

    if load_rule:
        matcher.add(pattern6985, 6985)

    pattern6986 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons17, cons261, cons1247, cons2033, cons_with_6986)
    def replacement6986(m, u, x):

        lst = PowerVariableExpn(u, m + S(1), x)
        rubi.append(6986)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule6986 = ReplacementRule(pattern6986, replacement6986)

    if load_rule:
        matcher.add(pattern6986, 6986)
    def With6987(m, x, u):
        k = Denominator(m)
        rubi.append(6987)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(S(1)/k)), x)
    pattern6987 = Pattern(Integral(u_*x_**m_, x_), cons367)
    rule6987 = ReplacementRule(pattern6987, With6987)

    if load_rule:
        matcher.add(pattern6987, 6987)

    pattern6988 = Pattern(Integral(u_, x_), cons2034, cons_with_6988)
    def replacement6988(x, u):

        lst = FunctionOfSquareRootOfQuadratic(u, x)
        rubi.append(6988)
        return Dist(S(2), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), x)
    rule6988 = ReplacementRule(pattern6988, replacement6988)

    if load_rule:
        matcher.add(pattern6988, 6988)
    pattern6989 = Pattern(Integral(S(1)/(a_ + v_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement6989(b, a, v, x):
        rubi.append(6989)
        return Dist(S(1)/(S(2)*a), Int(Together(S(1)/(-v/Rt(-a/b, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*a), Int(Together(S(1)/(v/Rt(-a/b, S(2)) + S(1))), x), x)
    rule6989 = ReplacementRule(pattern6989, replacement6989)

    if load_rule:
        matcher.add(pattern6989, 6989)
    pattern6990 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1479, cons744)
    def replacement6990(v, x, b, a, n):
        rubi.append(6990)
        return Dist(S(2)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x)
    rule6990 = ReplacementRule(pattern6990, replacement6990)

    if load_rule:
        matcher.add(pattern6990, 6990)
    pattern6991 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1482, cons165)
    def replacement6991(v, x, b, a, n):
        rubi.append(6991)
        return Dist(S(1)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x)
    rule6991 = ReplacementRule(pattern6991, replacement6991)

    if load_rule:
        matcher.add(pattern6991, 6991)
    pattern6992 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons148, cons2035)
    def replacement6992(v, u, x, b, a, n):
        rubi.append(6992)
        return Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x)
    rule6992 = ReplacementRule(pattern6992, replacement6992)

    if load_rule:
        matcher.add(pattern6992, 6992)

    pattern6993 = Pattern(Integral(u_, x_), cons_with_6993)
    def replacement6993(x, u):

        v = NormalizeIntegrand(u, x)
        rubi.append(6993)
        return Int(v, x)
    rule6993 = ReplacementRule(pattern6993, replacement6993)

    if load_rule:
        matcher.add(pattern6993, 6993)

    pattern6994 = Pattern(Integral(u_, x_), cons_with_6994)
    def replacement6994(x, u):

        v = ExpandIntegrand(u, x)
        rubi.append(6994)
        return Int(v, x)
    rule6994 = ReplacementRule(pattern6994, replacement6994)

    if load_rule:
        matcher.add(pattern6994, 6994)
    pattern6995 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons2036, cons1676, cons1255, cons2037)
    def replacement6995(d, q, p, c, u, m, x, b, a, n):
        rubi.append(6995)
        return Dist(x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q, Int(u*x**(m*p), x), x)
    rule6995 = ReplacementRule(pattern6995, replacement6995)

    if load_rule:
        matcher.add(pattern6995, 6995)
    pattern6996 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons347)
    def replacement6996(p, c, u, n2, x, b, a, n):
        rubi.append(6996)
        return Dist((S(4)*c)**(-p + S(1)/2)*sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule6996 = ReplacementRule(pattern6996, replacement6996)

    if load_rule:
        matcher.add(pattern6996, 6996)

    pattern6997 = Pattern(Integral(u_, x_), cons_with_6997)
    def replacement6997(x, u):

        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(6997)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule6997 = ReplacementRule(pattern6997, replacement6997)

    if load_rule:
        matcher.add(pattern6997, 6997)
    pattern6998 = Pattern(Integral(u_, x_))
    def replacement6998(x, u):
        rubi.append(6998)
        return Int(u, x)
    rule6998 = ReplacementRule(pattern6998, replacement6998)

    if load_rule:
        matcher.add(pattern6998, 6998)
    return matcher, [rule6931, rule6932, rule6933, rule6934, rule6935, rule6936, rule6937, rule6938, rule6939, rule6940, rule6941, rule6942, rule6943, rule6944, rule6945, rule6946, rule6947, rule6948, rule6949, rule6950, rule6951, rule6952, rule6953, rule6954, rule6955, rule6956, rule6957, rule6958, rule6959, rule6960, rule6961, rule6962, rule6963, rule6964, rule6965, rule6966, rule6967, rule6968, rule6969, rule6970, rule6971, rule6972, rule6973, rule6974, rule6975, rule6976, rule6977, rule6978, rule6979, rule6980, rule6981, rule6982, rule6983, rule6984, rule6985, rule6986, rule6987, rule6988, rule6989, rule6990, rule6991, rule6992, rule6993, rule6994, rule6995, rule6996, rule6997, rule6998, ]
