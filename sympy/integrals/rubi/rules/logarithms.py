"""
This code is automatically generated. Never edit it manually.
For details of generating the code see `rubi_parsing_guide.md` in `parsetools`.
"""

from sympy.external import import_module
matchpy = import_module("matchpy")

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
    i, ii, Pqq, Q, R, r, C, k, u = symbols('i ii Pqq Q R r C k u')
    _UseGamma = False
    ShowSteps = False
    StepCounter = None


def logarithms():
    from sympy.integrals.rubi.constraints import cons1158, cons8, cons29, cons50, cons127, cons5, cons52, cons89, cons90, cons2, cons3, cons1159, cons417, cons1160, cons1161, cons91, cons545, cons1162, cons210, cons211, cons586, cons4, cons68, cons19, cons1163, cons1164, cons1165, cons1166, cons1167, cons1168, cons1169, cons1170, cons1171, cons150, cons1172, cons1173, cons64, cons95, cons170, cons812, cons813, cons224, cons1174, cons226, cons798, cons81, cons1175, cons20, cons1176, cons1177, cons1178, cons1179, cons1180, cons1181, cons1182, cons1183, cons1184, cons1185, cons1186, cons1187, cons1188, cons1189, cons1190, cons1191, cons1192, cons799, cons1193, cons54, cons927, cons1194, cons1195, cons1196, cons1197, cons1198, cons1199, cons1200, cons1201, cons40, cons554, cons1202, cons1203, cons1204, cons27, cons654, cons1205, cons73, cons130, cons1206, cons1207, cons1208, cons1209, cons1210, cons148, cons1211, cons1212, cons13, cons165, cons1213, cons139, cons1214, cons1215, cons1216, cons1217, cons1218, cons1219, cons1220, cons1221, cons1222, cons1223, cons1224, cons72, cons1225, cons1226, cons808, cons842, cons1227, cons1228, cons70, cons1127, cons1229, cons1230, cons1231, cons1232, cons465, cons1233, cons1234, cons1235, cons1236, cons1237, cons1238, cons33, cons1101, cons1239, cons1057, cons517, cons818, cons819, cons1240, cons1241, cons1242, cons1243, cons1244, cons1245, cons1246, cons1247, cons36, cons37, cons1248, cons1249, cons1250


    pattern2009 = Pattern(Integral(log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))), x_), cons8, cons29, cons50, cons127, cons5, cons52, cons1158)
    rule2009 = ReplacementRule(pattern2009, replacement2009)

    pattern2010 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons89, cons90)
    rule2010 = ReplacementRule(pattern2010, replacement2010)

    pattern2011 = Pattern(Integral(S(1)/log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('d', S(1))), x_), cons29, cons50, cons127, cons1159)
    rule2011 = ReplacementRule(pattern2011, replacement2011)

    pattern2012 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons417)
    rule2012 = ReplacementRule(pattern2012, replacement2012)

    pattern2013 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons1160)
    rule2013 = ReplacementRule(pattern2013, replacement2013)

    pattern2014 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons1161)
    rule2014 = ReplacementRule(pattern2014, replacement2014)

    pattern2015 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons89, cons91)
    rule2015 = ReplacementRule(pattern2015, replacement2015)

    pattern2016 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons545)
    rule2016 = ReplacementRule(pattern2016, replacement2016)

    pattern2017 = Pattern(Integral(S(1)/((x_*WC('h', S(1)) + WC('g', S(0)))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1162)
    rule2017 = ReplacementRule(pattern2017, replacement2017)

    pattern2018 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons52, cons1162, cons586)
    rule2018 = ReplacementRule(pattern2018, replacement2018)

    pattern2019 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1162, cons68, cons89, cons90)
    rule2019 = ReplacementRule(pattern2019, replacement2019)

    pattern2020 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/log((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1))), x_), cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons1163, cons1162, cons1164)
    rule2020 = ReplacementRule(pattern2020, replacement2020)

    pattern2021 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_/log((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1))), x_), cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons1163, cons1162, cons1165)
    rule2021 = ReplacementRule(pattern2021, replacement2021)

    pattern2022 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1162, cons68)
    rule2022 = ReplacementRule(pattern2022, replacement2022)

    pattern2023 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1162, cons68, cons1166)
    rule2023 = ReplacementRule(pattern2023, replacement2023)

    pattern2024 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1162, cons68, cons1167)
    rule2024 = ReplacementRule(pattern2024, replacement2024)

    pattern2025 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1162, cons68, cons89, cons91)
    rule2025 = ReplacementRule(pattern2025, replacement2025)

    pattern2026 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons1162, cons68)
    rule2026 = ReplacementRule(pattern2026, replacement2026)

    pattern2027 = Pattern(Integral(log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons8, cons50, cons127, cons210, cons211, cons1168)
    rule2027 = ReplacementRule(pattern2027, replacement2027)

    pattern2028 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1))))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons211, cons1169, cons1170)
    rule2028 = ReplacementRule(pattern2028, replacement2028)

    pattern2029 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1171, cons150)
    rule2029 = ReplacementRule(pattern2029, replacement2029)

    pattern2030 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1171, cons68)
    rule2030 = ReplacementRule(pattern2030, replacement2030)

    pattern2031 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_/(x_*WC('h', S(1)) + WC('g', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1171, cons89, cons90)
    rule2031 = ReplacementRule(pattern2031, replacement2031)

    pattern2032 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons5, cons52, cons1171, cons89, cons90, cons68, cons1172, cons1173)
    rule2032 = ReplacementRule(pattern2032, replacement2032)

    pattern2033 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1171, cons64)
    rule2033 = ReplacementRule(pattern2033, replacement2033)

    pattern2034 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1171, cons95, cons91, cons170)
    rule2034 = ReplacementRule(pattern2034, replacement2034)

    pattern2035 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons52, cons1171, cons64)
    rule2035 = ReplacementRule(pattern2035, replacement2035)

    pattern2036 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((v_**p_*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons812, cons813)
    rule2036 = ReplacementRule(pattern2036, replacement2036)

    pattern2037 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons52, cons224)
    rule2037 = ReplacementRule(pattern2037, replacement2037)

    pattern2038 = Pattern(Integral(log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0))))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons8, cons50, cons127, cons210, cons211, cons226, cons798, cons1162, cons1174)
    rule2038 = ReplacementRule(pattern2038, replacement2038)

    pattern2039 = Pattern(Integral((a_ + WC('b', S(1))*log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0)))))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons211, cons226, cons798, cons1162, cons1174)
    rule2039 = ReplacementRule(pattern2039, replacement2039)

    pattern2040 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons5, cons52, cons1162, cons81)
    rule2040 = ReplacementRule(pattern2040, With2040)

    pattern2041 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons211, cons226, cons798, cons4, cons1162, cons64, cons1175)
    rule2041 = ReplacementRule(pattern2041, replacement2041)

    pattern2042 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons5, cons52, cons20, cons150, CustomConstraint(With2042))
    rule2042 = ReplacementRule(pattern2042, replacement2042)

    pattern2043 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons19, cons4, cons5, cons52, cons1176)
    rule2043 = ReplacementRule(pattern2043, replacement2043)

    pattern2044 = Pattern(Integral(log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0))))/(g_ + x_**S(2)*WC('h', S(1))), x_), cons8, cons50, cons127, cons210, cons211, cons1177, cons1178)
    rule2044 = ReplacementRule(pattern2044, replacement2044)

    pattern2045 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0)))))/(g_ + x_**S(2)*WC('h', S(1))), x_), cons8, cons50, cons127, cons210, cons211, cons1177, cons1179, cons1180)
    rule2045 = ReplacementRule(pattern2045, replacement2045)

    pattern2046 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons5, cons52, cons1181)
    rule2046 = ReplacementRule(pattern2046, replacement2046)

    pattern2047 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((e_ + x_*WC('f', S(1)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(g_ + x_**S(2)*WC('i', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons226, cons5, cons52, cons1182)
    rule2047 = ReplacementRule(pattern2047, replacement2047)

    pattern2048 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/sqrt(g_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1183)
    rule2048 = ReplacementRule(pattern2048, With2048)

    pattern2049 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(sqrt(g1_ + x_*WC('h1', S(1)))*sqrt(g2_ + x_*WC('h2', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1187, cons1188, cons1189, cons1190, cons5, cons52, cons1184, cons1185, cons1186)
    rule2049 = ReplacementRule(pattern2049, With2049)

    pattern2050 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/sqrt(g_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons1191)
    rule2050 = ReplacementRule(pattern2050, replacement2050)

    pattern2051 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(sqrt(g1_ + x_*WC('h1', S(1)))*sqrt(g2_ + x_*WC('h2', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1187, cons1188, cons1189, cons1190, cons5, cons52, cons1184)
    rule2051 = ReplacementRule(pattern2051, replacement2051)

    pattern2052 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('k', S(1)) + WC('j', S(0)))*WC('i', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons799, cons5, cons52, cons89, cons90, cons1192)
    rule2052 = ReplacementRule(pattern2052, replacement2052)

    pattern2053 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('k', S(1)) + WC('j', S(0)))**WC('m', S(1))*WC('i', S(1)) + S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons799, cons19, cons5, cons52, cons89, cons90, cons1193)
    rule2053 = ReplacementRule(pattern2053, replacement2053)

    pattern2054 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*PolyLog(r_, (x_*WC('k', S(1)) + WC('j', S(0)))**WC('m', S(1))*WC('i', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons799, cons19, cons5, cons52, cons54, cons89, cons90, cons1193)
    rule2054 = ReplacementRule(pattern2054, replacement2054)

    pattern2055 = Pattern(Integral(F_**(x_*WC('h', S(1)) + WC('g', S(0)))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))*WC('Px', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons52, cons927, cons64, cons1194)
    rule2055 = ReplacementRule(pattern2055, With2055)

    pattern2056 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((e_ + x_**m_*WC('f', S(1)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons150)
    rule2056 = ReplacementRule(pattern2056, replacement2056)

    pattern2057 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_**m_*(f_ + x_**WC('r', S(1))*WC('e', S(1))))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons1195, cons150)
    rule2057 = ReplacementRule(pattern2057, replacement2057)

    pattern2058 = Pattern(Integral(x_**WC('r1', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_**r_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons54, cons1196)
    rule2058 = ReplacementRule(pattern2058, replacement2058)

    pattern2059 = Pattern(Integral(x_**WC('r1', S(1))*(x_**r_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_**r_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons52, cons54, cons1196)
    rule2059 = ReplacementRule(pattern2059, replacement2059)

    pattern2060 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1197)
    rule2060 = ReplacementRule(pattern2060, With2060)

    pattern2061 = Pattern(Integral(log((x_**mn_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))/(x_*(d_ + x_**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1198, cons1199)
    rule2061 = ReplacementRule(pattern2061, replacement2061)

    pattern2062 = Pattern(Integral(log(x_**mn_*(x_**WC('n', S(1))*WC('a', S(1)) + WC('b', S(0)))*WC('c', S(1)))/(x_*(d_ + x_**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1198, cons1199)
    rule2062 = ReplacementRule(pattern2062, replacement2062)

    pattern2063 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons927)
    rule2063 = ReplacementRule(pattern2063, replacement2063)

    pattern2064 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons1200, cons150, CustomConstraint(With2064))
    rule2064 = ReplacementRule(pattern2064, replacement2064)

    pattern2065 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons1200, cons150, CustomConstraint(With2065))
    rule2065 = ReplacementRule(pattern2065, replacement2065)

    pattern2066 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_**S(2)*WC('g', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons52, cons4, cons1201, cons40)
    rule2066 = ReplacementRule(pattern2066, replacement2066)

    pattern2067 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((v_**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons554, cons1202)
    rule2067 = ReplacementRule(pattern2067, replacement2067)

    pattern2068 = Pattern(Integral(log(((x_**WC('n', S(1))*WC('c', S(1)))**p_*WC('b', S(1)))**q_*WC('a', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons52, cons54, cons1203)
    rule2068 = ReplacementRule(pattern2068, replacement2068)

    pattern2069 = Pattern(Integral(x_**WC('m', S(1))*log(((x_**WC('n', S(1))*WC('c', S(1)))**p_*WC('b', S(1)))**q_*WC('a', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons52, cons54, cons68, cons1204)
    rule2069 = ReplacementRule(pattern2069, replacement2069)

    pattern2070 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e1', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons654, cons27)
    rule2070 = ReplacementRule(pattern2070, replacement2070)

    pattern2071 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons654, cons1206, cons1205, cons73, cons130)
    rule2071 = ReplacementRule(pattern2071, replacement2071)

    pattern2072 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1207, cons1208)
    rule2072 = ReplacementRule(pattern2072, replacement2072)

    pattern2073 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1207, cons130)
    rule2073 = ReplacementRule(pattern2073, replacement2073)

    pattern2074 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1209, cons130)
    rule2074 = ReplacementRule(pattern2074, replacement2074)

    pattern2075 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1210)
    rule2075 = ReplacementRule(pattern2075, replacement2075)

    pattern2076 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1210, cons40, cons148)
    rule2076 = ReplacementRule(pattern2076, replacement2076)

    pattern2077 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1207)
    rule2077 = ReplacementRule(pattern2077, replacement2077)

    pattern2078 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1207)
    rule2078 = ReplacementRule(pattern2078, replacement2078)

    pattern2079 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1209)
    rule2079 = ReplacementRule(pattern2079, replacement2079)

    pattern2080 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1211, cons130)
    rule2080 = ReplacementRule(pattern2080, replacement2080)

    pattern2081 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1210, cons130)
    rule2081 = ReplacementRule(pattern2081, replacement2081)

    pattern2082 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0)))**S(3), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1211, cons1207)
    rule2082 = ReplacementRule(pattern2082, replacement2082)

    pattern2083 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0)))**S(3), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1210, cons1209)
    rule2083 = ReplacementRule(pattern2083, replacement2083)

    pattern2084 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons130, cons20, cons68)
    rule2084 = ReplacementRule(pattern2084, replacement2084)

    pattern2085 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_**n_*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1213, cons1212, cons73, cons68, cons13, cons165)
    rule2085 = ReplacementRule(pattern2085, replacement2085)

    pattern2086 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_)**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons1213, cons1212, cons73, cons68, cons13, cons165)
    rule2086 = ReplacementRule(pattern2086, replacement2086)

    pattern2087 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))/log(u_**n_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1213, cons1212, cons73, cons68)
    rule2087 = ReplacementRule(pattern2087, replacement2087)

    pattern2088 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))/log(u_), x_), cons2, cons3, cons8, cons29, cons1213, cons1212, cons73, cons68)
    rule2088 = ReplacementRule(pattern2088, replacement2088)

    pattern2089 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_**n_*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1213, cons1212, cons73, cons68, cons13, cons139)
    rule2089 = ReplacementRule(pattern2089, replacement2089)

    pattern2090 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_)**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons1213, cons1212, cons73, cons68, cons13, cons139)
    rule2090 = ReplacementRule(pattern2090, replacement2090)

    pattern2091 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1211, cons1207)
    rule2091 = ReplacementRule(pattern2091, replacement2091)

    pattern2092 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1211, cons1210, cons1214)
    rule2092 = ReplacementRule(pattern2092, replacement2092)

    pattern2093 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons654, cons1206, cons1205, cons73, cons1211, cons1210, cons130)
    rule2093 = ReplacementRule(pattern2093, replacement2093)

    pattern2094 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(f_ + x_**S(2)*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1215, cons1216)
    rule2094 = ReplacementRule(pattern2094, replacement2094)

    pattern2095 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1217)
    rule2095 = ReplacementRule(pattern2095, replacement2095)

    pattern2096 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons211, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1218)
    rule2096 = ReplacementRule(pattern2096, replacement2096)

    pattern2097 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1210, cons1209)
    rule2097 = ReplacementRule(pattern2097, replacement2097)

    pattern2098 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1219, cons1213, cons73, cons13, cons165)
    rule2098 = ReplacementRule(pattern2098, replacement2098)

    pattern2099 = Pattern(Integral(log(u_)**WC('p', S(1))*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1219, cons1213, cons73, cons13, cons165)
    rule2099 = ReplacementRule(pattern2099, replacement2099)

    pattern2100 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1219, cons1213, cons73, cons13, cons139)
    rule2100 = ReplacementRule(pattern2100, With2100)

    pattern2101 = Pattern(Integral(log(u_)**p_*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1219, cons1213, cons73, cons13, cons139)
    rule2101 = ReplacementRule(pattern2101, With2101)

    pattern2102 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1220, cons1213, cons73, cons13, cons165)
    rule2102 = ReplacementRule(pattern2102, replacement2102)

    pattern2103 = Pattern(Integral(log(u_)**WC('p', S(1))*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1220, cons1213, cons73, cons13, cons165)
    rule2103 = ReplacementRule(pattern2103, replacement2103)

    pattern2104 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1220, cons1213, cons73, cons13, cons139)
    rule2104 = ReplacementRule(pattern2104, With2104)

    pattern2105 = Pattern(Integral(log(u_)**p_*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1220, cons1213, cons73, cons13, cons139)
    rule2105 = ReplacementRule(pattern2105, With2105)

    pattern2106 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons1221, cons1213, cons73, cons13, cons148)
    rule2106 = ReplacementRule(pattern2106, replacement2106)

    pattern2107 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons52, cons1221, cons1213, cons73, cons13, cons148)
    rule2107 = ReplacementRule(pattern2107, replacement2107)

    pattern2108 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons1221, cons1213, cons73, cons13, cons139)
    rule2108 = ReplacementRule(pattern2108, replacement2108)

    pattern2109 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons52, cons1221, cons1213, cons73, cons13, cons139)
    rule2109 = ReplacementRule(pattern2109, replacement2109)

    pattern2110 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons1222, cons1213, cons73, cons13, cons148)
    rule2110 = ReplacementRule(pattern2110, replacement2110)

    pattern2111 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons52, cons1222, cons1213, cons73, cons13, cons148)
    rule2111 = ReplacementRule(pattern2111, replacement2111)

    pattern2112 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons1222, cons1213, cons73, cons13, cons139)
    rule2112 = ReplacementRule(pattern2112, replacement2112)

    pattern2113 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons52, cons1222, cons1213, cons73, cons13, cons139)
    rule2113 = ReplacementRule(pattern2113, replacement2113)

    pattern2114 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1223, cons1224)
    rule2114 = ReplacementRule(pattern2114, replacement2114)

    pattern2115 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons211, cons4, cons5, cons654, cons1206, cons1205, cons73, cons1223, cons72)
    rule2115 = ReplacementRule(pattern2115, replacement2115)

    pattern2116 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(f_ + x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons654, cons127, cons210, cons211, cons4, cons1206, cons1205)
    rule2116 = ReplacementRule(pattern2116, With2116)

    pattern2117 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(f_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons654, cons127, cons211, cons4, cons1206, cons1205)
    rule2117 = ReplacementRule(pattern2117, With2117)

    pattern2118 = Pattern(Integral(RFx_*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons654, cons1206, cons1205, cons1200, cons130, CustomConstraint(With2118))
    rule2118 = ReplacementRule(pattern2118, replacement2118)

    pattern2119 = Pattern(Integral(WC('u', S(1))*log(v_)**WC('p', S(1)), x_), cons5, cons1225, cons1226, CustomConstraint(With2119))
    rule2119 = ReplacementRule(pattern2119, replacement2119)

    pattern2120 = Pattern(Integral(log((x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons4, cons5, cons808)
    rule2120 = ReplacementRule(pattern2120, replacement2120)

    pattern2121 = Pattern(Integral(log(v_**WC('p', S(1))*WC('c', S(1))), x_), cons8, cons5, cons842, cons1227)
    rule2121 = ReplacementRule(pattern2121, replacement2121)

    pattern2122 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons1228)
    rule2122 = ReplacementRule(pattern2122, replacement2122)

    pattern2123 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons68)
    rule2123 = ReplacementRule(pattern2123, replacement2123)

    pattern2124 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(v_**WC('p', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons5, cons70, cons842, cons1127)
    rule2124 = ReplacementRule(pattern2124, replacement2124)

    pattern2125 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))*asin(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons64)
    rule2125 = ReplacementRule(pattern2125, With2125)

    pattern2126 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))/(x_**S(2)*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1229)
    rule2126 = ReplacementRule(pattern2126, With2126)

    pattern2127 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons150)
    rule2127 = ReplacementRule(pattern2127, replacement2127)

    pattern2128 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons150, cons1230)
    rule2128 = ReplacementRule(pattern2128, replacement2128)

    pattern2129 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons150, cons1231)
    rule2129 = ReplacementRule(pattern2129, replacement2129)

    pattern2130 = Pattern(Integral(u_*log(v_), x_), CustomConstraint(With2130))
    rule2130 = ReplacementRule(pattern2130, replacement2130)

    pattern2131 = Pattern(Integral(w_*(WC('a', S(0)) + WC('b', S(1))*log(u_))*log(v_), x_), cons2, cons3, cons1232, CustomConstraint(With2131))
    rule2131 = ReplacementRule(pattern2131, replacement2131)

    pattern2132 = Pattern(Integral(log((a_ + (x_*WC('e', S(1)) + WC('d', S(0)))**n_*WC('b', S(1)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons89, cons465)
    rule2132 = ReplacementRule(pattern2132, replacement2132)

    pattern2133 = Pattern(Integral(log((a_ + (x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1233)
    rule2133 = ReplacementRule(pattern2133, replacement2133)

    pattern2134 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((d_ + WC('e', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))))**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons150)
    rule2134 = ReplacementRule(pattern2134, replacement2134)

    pattern2135 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons5, cons1200, cons150)
    rule2135 = ReplacementRule(pattern2135, replacement2135)

    pattern2136 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1200, cons150)
    rule2136 = ReplacementRule(pattern2136, replacement2136)

    pattern2137 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1200, cons150, cons1234, cons68)
    rule2137 = ReplacementRule(pattern2137, replacement2137)

    pattern2138 = Pattern(Integral(log(RFx_**WC('n', S(1))*WC('c', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons8, cons29, cons50, cons4, cons1200, cons1235)
    rule2138 = ReplacementRule(pattern2138, With2138)

    pattern2139 = Pattern(Integral(log(Px_**WC('n', S(1))*WC('c', S(1)))/Qx_, x_), cons8, cons4, cons1236, cons1237)
    rule2139 = ReplacementRule(pattern2139, With2139)

    pattern2140 = Pattern(Integral(RGx_*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons5, cons1200, cons1238, cons150, CustomConstraint(With2140))
    rule2140 = ReplacementRule(pattern2140, replacement2140)

    pattern2141 = Pattern(Integral(RGx_*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons5, cons1200, cons1238, cons150, CustomConstraint(With2141))
    rule2141 = ReplacementRule(pattern2141, replacement2141)

    pattern2142 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(u_)), x_), cons2, cons3, cons1200, CustomConstraint(With2142))
    rule2142 = ReplacementRule(pattern2142, replacement2142)

    pattern2143 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log((F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1))*WC('e', S(1)) + S(1)), x_), cons1101, cons2, cons3, cons8, cons50, cons127, cons210, cons4, cons33, cons170)
    rule2143 = ReplacementRule(pattern2143, replacement2143)

    pattern2144 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log(d_ + (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1))*WC('e', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons33, cons170, cons1239)
    rule2144 = ReplacementRule(pattern2144, replacement2144)

    pattern2145 = Pattern(Integral(log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1057)
    rule2145 = ReplacementRule(pattern2145, replacement2145)

    pattern2146 = Pattern(Integral(log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons1057)
    rule2146 = ReplacementRule(pattern2146, replacement2146)

    pattern2147 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons1057, cons68, cons517)
    rule2147 = ReplacementRule(pattern2147, replacement2147)

    pattern2148 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons19, cons1057, cons68, cons517)
    rule2148 = ReplacementRule(pattern2148, replacement2148)

    pattern2149 = Pattern(Integral(WC('v', S(1))*log(sqrt(u_)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons29, cons50, cons127, cons818, cons819, cons1240)
    rule2149 = ReplacementRule(pattern2149, replacement2149)

    pattern2150 = Pattern(Integral(log(u_), x_), cons1232)
    rule2150 = ReplacementRule(pattern2150, replacement2150)

    pattern2151 = Pattern(Integral(log(u_)/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons1241, cons1242)
    rule2151 = ReplacementRule(pattern2151, replacement2151)

    pattern2152 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*log(u_), x_), cons2, cons3, cons19, cons1232, cons68)
    rule2152 = ReplacementRule(pattern2152, replacement2152)

    pattern2153 = Pattern(Integral(log(u_)/Qx_, x_), cons1243, cons1232)
    rule2153 = ReplacementRule(pattern2153, With2153)

    pattern2154 = Pattern(Integral(u_**(x_*WC('a', S(1)))*log(u_), x_), cons2, cons1232)
    rule2154 = ReplacementRule(pattern2154, replacement2154)

    pattern2155 = Pattern(Integral(v_*log(u_), x_), cons1232, CustomConstraint(With2155))
    rule2155 = ReplacementRule(pattern2155, replacement2155)

    pattern2156 = Pattern(Integral(log(v_)*log(w_), x_), cons1244, cons1245)
    rule2156 = ReplacementRule(pattern2156, replacement2156)

    pattern2157 = Pattern(Integral(u_*log(v_)*log(w_), x_), cons1244, cons1245, CustomConstraint(With2157))
    rule2157 = ReplacementRule(pattern2157, replacement2157)

    pattern2158 = Pattern(Integral(log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons4, cons5, cons1246)
    rule2158 = ReplacementRule(pattern2158, replacement2158)

    pattern2159 = Pattern(Integral(log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)))/x_, x_), cons2, cons3, cons4, cons5, cons1246)
    rule2159 = ReplacementRule(pattern2159, replacement2159)

    pattern2160 = Pattern(Integral(x_**WC('m', S(1))*log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons5, cons68)
    rule2160 = ReplacementRule(pattern2160, replacement2160)

    pattern2161 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*log(x_*WC('d', S(1)) + WC('c', S(0))))/sqrt(a_ + WC('b', S(1))*log(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons1247)
    rule2161 = ReplacementRule(pattern2161, replacement2161)

    pattern2162 = Pattern(Integral(f_**(WC('a', S(1))*log(u_)), x_), cons2, cons127, cons1248)
    rule2162 = ReplacementRule(pattern2162, replacement2162)

    pattern2163 = Pattern(Integral(u_, x_), cons1249, CustomConstraint(With2163))
    rule2163 = ReplacementRule(pattern2163, replacement2163)

    pattern2164 = Pattern(Integral(WC('u', S(1))*log(Gamma(v_)), x_))
    rule2164 = ReplacementRule(pattern2164, replacement2164)

    pattern2165 = Pattern(Integral((w_*WC('a', S(1)) + w_*WC('b', S(1))*log(v_)**WC('n', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons40)
    rule2165 = ReplacementRule(pattern2165, replacement2165)

    pattern2166 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons52, cons1250)
    rule2166 = ReplacementRule(pattern2166, replacement2166)
    return [rule2009, rule2010, rule2011, rule2012, rule2013, rule2014, rule2015, rule2016, rule2017, rule2018, rule2019, rule2020, rule2021, rule2022, rule2023, rule2024, rule2025, rule2026, rule2027, rule2028, rule2029, rule2030, rule2031, rule2032, rule2033, rule2034, rule2035, rule2036, rule2037, rule2038, rule2039, rule2040, rule2041, rule2042, rule2043, rule2044, rule2045, rule2046, rule2047, rule2048, rule2049, rule2050, rule2051, rule2052, rule2053, rule2054, rule2055, rule2056, rule2057, rule2058, rule2059, rule2060, rule2061, rule2062, rule2063, rule2064, rule2065, rule2066, rule2067, rule2068, rule2069, rule2070, rule2071, rule2072, rule2073, rule2074, rule2075, rule2076, rule2077, rule2078, rule2079, rule2080, rule2081, rule2082, rule2083, rule2084, rule2085, rule2086, rule2087, rule2088, rule2089, rule2090, rule2091, rule2092, rule2093, rule2094, rule2095, rule2096, rule2097, rule2098, rule2099, rule2100, rule2101, rule2102, rule2103, rule2104, rule2105, rule2106, rule2107, rule2108, rule2109, rule2110, rule2111, rule2112, rule2113, rule2114, rule2115, rule2116, rule2117, rule2118, rule2119, rule2120, rule2121, rule2122, rule2123, rule2124, rule2125, rule2126, rule2127, rule2128, rule2129, rule2130, rule2131, rule2132, rule2133, rule2134, rule2135, rule2136, rule2137, rule2138, rule2139, rule2140, rule2141, rule2142, rule2143, rule2144, rule2145, rule2146, rule2147, rule2148, rule2149, rule2150, rule2151, rule2152, rule2153, rule2154, rule2155, rule2156, rule2157, rule2158, rule2159, rule2160, rule2161, rule2162, rule2163, rule2164, rule2165, rule2166, ]





def replacement2009(c, d, e, f, p, q, x):
    return Simp((e + f*x)*log(c*(d*(e + f*x)**p)**q)/f, x) - Simp(p*q*x, x)


def replacement2010(a, b, c, d, e, f, n, p, q, x):
    return -Dist(b*n*p*q, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)/f, x)


def replacement2011(d, e, f, x):
    return Simp(LogIntegral(d*(e + f*x))/(d*f), x)


def replacement2012(a, b, c, d, e, f, p, q, x):
    return Simp((c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*ExpIntegralEi((a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))*exp(-a/(b*p*q))/(b*f*p*q), x)


def replacement2013(a, b, c, d, e, f, p, q, x):
    return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*Erfi(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))/Rt(b*p*q, S(2)))*Rt(b*p*q, S(2))*exp(-a/(b*p*q))/(b*f*p*q), x)


def replacement2014(a, b, c, d, e, f, p, q, x):
    return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*Erf(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))/Rt(-b*p*q, S(2)))*Rt(-b*p*q, S(2))*exp(-a/(b*p*q))/(b*f*p*q), x)


def replacement2015(a, b, c, d, e, f, n, p, q, x):
    return -Dist(S(1)/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(e + f*x)/(b*f*p*q*(n + S(1))), x)


def replacement2016(a, b, c, d, e, f, n, p, q, x):
    return Simp((c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(-(a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))**(-n)*(a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)*Gamma(n + S(1), -(a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))*exp(-a/(b*p*q))/f, x)


def replacement2017(a, b, c, d, e, f, g, h, p, q, x):
    return Simp(log(RemoveContent(a + b*log(c*(d*(e + f*x)**p)**q), x))/(b*h*p*q), x)


def replacement2018(a, b, c, d, e, f, g, h, n, p, q, x):
    return Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))/(b*h*p*q*(n + S(1))), x)


def replacement2019(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return -Dist(b*n*p*q/(m + S(1)), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*(g + h*x)**m, x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)


def replacement2020(d, e, f, g, h, m, p, x):
    return Simp((h/f)**(p + S(-1))*LogIntegral(d*(e + f*x)**p)/(d*f*p), x)


def replacement2021(d, e, f, g, h, m, p, x):
    return Dist((e + f*x)**(S(1) - p)*(g + h*x)**(p + S(-1)), Int((e + f*x)**(p + S(-1))/log(d*(e + f*x)**p), x), x)


def replacement2022(a, b, c, d, e, f, g, h, m, p, q, x):
    return Simp((c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*ExpIntegralEi((a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q), x)


def replacement2023(a, b, c, d, e, f, g, h, m, p, q, x):
    return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*Erfi(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))*Rt((m + S(1))/(b*p*q), S(2)))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q*Rt((m + S(1))/(b*p*q), S(2))), x)


def replacement2024(a, b, c, d, e, f, g, h, m, p, q, x):
    return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*Erf(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))*Rt(-(m + S(1))/(b*p*q), S(2)))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q*Rt(-(m + S(1))/(b*p*q), S(2))), x)


def replacement2025(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return -Dist((m + S(1))/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**m, x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**(m + S(1))/(b*h*p*q*(n + S(1))), x)


def replacement2026(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Simp((c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(-(a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))**(-n)*(a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))*Gamma(n + S(1), -(a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))*exp(-a*(m + S(1))/(b*p*q))/(h*(m + S(1))), x)


def replacement2027(c, e, f, g, h, x):
    return -Simp(PolyLog(S(2), -(g + h*x)*Together(c*f/h))/h, x)


def replacement2028(a, b, c, e, f, g, h, x):
    return Dist(b, Int(log(-h*(e + f*x)/(-e*h + f*g))/(g + h*x), x), x) + Simp((a + b*log(c*(e - f*g/h)))*log(g + h*x)/h, x)


def replacement2029(a, b, c, d, e, f, g, h, n, p, q, x):
    return -Dist(b*f*n*p*q/h, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*log(f*(g + h*x)/(-e*h + f*g))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*log(f*(g + h*x)/(-e*h + f*g))/h, x)


def replacement2030(a, b, c, d, e, f, g, h, m, p, q, x):
    return -Dist(b*f*p*q/(h*(m + S(1))), Int((g + h*x)**(m + S(1))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)


def replacement2031(a, b, c, d, e, f, g, h, n, p, q, x):
    return -Dist(b*f*n*p*q/(-e*h + f*g), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))/(g + h*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)/((g + h*x)*(-e*h + f*g)), x)


def replacement2032(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return -Dist(b*f*n*p*q/(h*(m + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*(g + h*x)**(m + S(1))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)


def replacement2033(a, b, c, d, e, f, g, h, m, p, q, x):
    return Int(ExpandIntegrand((g + h*x)**m/(a + b*log(c*(d*(e + f*x)**p)**q)), x), x)


def replacement2034(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return -Dist((m + S(1))/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**m, x), x) + Dist(m*(-e*h + f*g)/(b*f*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**(m + S(-1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(e + f*x)*(g + h*x)**m/(b*f*p*q*(n + S(1))), x)


def replacement2035(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Int(ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x), x)


def replacement2036(a, b, c, d, m, n, p, q, u, v, x):
    return Int((a + b*log(c*(d*ExpandToSum(v, x)**p)**q))**n*ExpandToSum(u, x)**m, x)


def replacement2037(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x)


def replacement2038(c, e, f, g, h, i, j, x):
    return Simp(f*PolyLog(S(2), f*(i + j*x)/(j*(e + f*x)))/(h*(-e*j + f*i)), x)


def replacement2039(a, b, c, e, f, g, h, i, j, x):
    return Dist(a, Int(S(1)/((g + h*x)*(i + j*x)), x), x) + Dist(b, Int(log(c/(e + f*x))/((g + h*x)*(i + j*x)), x), x)


def With2040(a, b, c, d, e, f, g, h, i, j, m, p, q, x):
    u = IntHide((i + j*x)**m/(g + h*x), x)
    return -Dist(b*h*p*q, Int(SimplifyIntegrand(u/(g + h*x), x), x), x) + Dist(a + b*log(c*(d*(e + f*x)**p)**q), u, x)


def replacement2041(a, b, c, e, f, g, h, i, j, m, n, x):
    return Dist(c**(-m)*f**(-m)/h, Subst(Int((a + b*x)**n*(-c*e*j + c*f*i + j*exp(x))**m, x), x, log(c*(e + f*x))), x)


def With2042(a, b, c, d, e, f, g, h, i, j, m, n, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, (i + j*x)**m/(g + h*x), x)
    if SumQ(u):
        return True
    return False


def replacement2042(a, b, c, d, e, f, g, h, i, j, m, n, p, q, x):

    u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, (i + j*x)**m/(g + h*x), x)
    return Int(u, x)


def replacement2043(a, b, c, d, e, f, g, h, i, j, m, n, p, q, x):
    return Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(i + j*x)**m/(g + h*x), x)


def replacement2044(c, e, f, g, h, x):
    return -Simp(f*PolyLog(S(2), (-e + f*x)/(e + f*x))/(S(2)*e*h), x)


def replacement2045(a, b, c, e, f, g, h, x):
    return Dist(b, Int(log(S(2)*e/(e + f*x))/(g + h*x**S(2)), x), x) + Dist(a + b*log(c/(S(2)*e)), Int(S(1)/(g + h*x**S(2)), x), x)


def replacement2046(a, b, c, d, e, f, g, h, i, p, q, x):
    return Dist(e*f, Int((a + b*log(c*(d*(e + f*x)**p)**q))/((e + f*x)*(e*i*x + f*g)), x), x)


def replacement2047(a, b, c, d, e, f, g, i, p, q, x):
    return Dist(e*f, Int((a + b*log(c*(d*(e + f*x)**p)**q))/((e + f*x)*(e*i*x + f*g)), x), x)


def With2048(a, b, c, d, e, f, g, h, p, q, x):
    u = IntHide(S(1)/sqrt(g + h*x**S(2)), x)
    return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Simp(u*(a + b*log(c*(d*(e + f*x)**p)**q)), x)


def With2049(a, b, c, d, e, f, g1, g2, h1, h2, p, q, x):
    u = IntHide(S(1)/sqrt(g1*g2 + h1*h2*x**S(2)), x)
    return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Simp(u*(a + b*log(c*(d*(e + f*x)**p)**q)), x)


def replacement2050(a, b, c, d, e, f, g, h, p, q, x):
    return Dist(sqrt(S(1) + h*x**S(2)/g)/sqrt(g + h*x**S(2)), Int((a + b*log(c*(d*(e + f*x)**p)**q))/sqrt(S(1) + h*x**S(2)/g), x), x)


def replacement2051(a, b, c, d, e, f, g1, g2, h1, h2, p, q, x):
    return Dist(sqrt(S(1) + h1*h2*x**S(2)/(g1*g2))/(sqrt(g1 + h1*x)*sqrt(g2 + h2*x)), Int((a + b*log(c*(d*(e + f*x)**p)**q))/sqrt(S(1) + h1*h2*x**S(2)/(g1*g2)), x), x)


def replacement2052(a, b, c, d, e, f, g, h, i, j, k, n, p, q, x):
    return Dist(b*f*n*p*q/h, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(S(2), Together(-i*(j + k*x) + S(1)))/(e + f*x), x), x) - Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(S(2), Together(-i*(j + k*x) + S(1)))/h, x)


def replacement2053(a, b, c, d, e, f, g, h, i, j, k, m, n, p, q, x):
    return Dist(b*f*n*p*q/(h*m), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(S(2), -i*(j + k*x)**m)/(e + f*x), x), x) - Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(S(2), -i*(j + k*x)**m)/(h*m), x)


def replacement2054(a, b, c, d, e, f, g, h, i, j, k, m, n, p, q, r, x):
    return -Dist(b*f*n*p*q/(h*m), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(r + S(1), i*(j + k*x)**m)/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(r + S(1), i*(j + k*x)**m)/(h*m), x)


def With2055(F, Px, a, b, c, d, e, f, g, h, m, p, q, x):
    u = IntHide(Px*F(g + h*x)**m, x)
    return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Dist(a + b*log(c*(d*(e + f*x)**p)**q), u, x)


def replacement2056(a, b, c, d, e, f, m, n, p, q, x):
    return Dist(S(1)/m, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n/x, x), x, x**m), x)


def replacement2057(a, b, c, d, e, f, m, n, p, q, r, x):
    return Dist(S(1)/m, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n/x, x), x, x**m), x)


def replacement2058(a, b, c, d, e, f, n, p, q, r, r1, x):
    return Dist(S(1)/r, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n, x), x, x**r), x)


def replacement2059(a, b, c, d, e, f, g, h, m, n, p, q, r, r1, x):
    return Dist(S(1)/r, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x), x, x**r), x)


def With2060(a, b, c, d, e, n, x):
    u = IntHide(S(1)/(d + e*x**S(2)), x)
    return -Dist(b*n, Int(u/x, x), x) + Dist(a + b*log(c*x**n), u, x)


def replacement2061(a, b, c, d, e, mn, n, x):
    return Simp(PolyLog(S(2), -Together(b*c*x**(-n)*(d + e*x**n)/d))/(d*n), x)


def replacement2062(a, b, c, d, e, mn, n, x):
    return Simp(PolyLog(S(2), -Together(b*c*x**(-n)*(d + e*x**n)/d))/(d*n), x)


def replacement2063(Px, a, b, c, d, e, f, n, p, q, x):
    return Int(ExpandIntegrand(Px*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x), x)


def With2064(RFx, a, b, c, d, e, f, n, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement2064(RFx, a, b, c, d, e, f, n, p, q, x):

    u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, RFx, x)
    return Int(u, x)


def With2065(RFx, a, b, c, d, e, f, n, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(RFx*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
    if SumQ(u):
        return True
    return False


def replacement2065(RFx, a, b, c, d, e, f, n, p, q, x):

    u = ExpandIntegrand(RFx*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
    return Int(u, x)


def replacement2066(a, b, c, d, e, f, g, n, p, q, u, x):
    return Int(u*(a + b*log(c*(S(4)**(-p)*d*g**(-p)*(f + S(2)*g*x)**(S(2)*p))**q))**n, x)


def replacement2067(a, b, c, d, n, p, q, u, v, x):
    return Int(u*(a + b*log(c*(d*ExpandToSum(v, x)**p)**q))**n, x)


def replacement2068(a, b, c, n, p, q, r, x):
    return Subst(Int(log(x**(n*p*q))**r, x), x**(n*p*q), a*(b*(c*x**n)**p)**q)


def replacement2069(a, b, c, m, n, p, q, r, x):
    return Subst(Int(x**m*log(x**(n*p*q))**r, x), x**(n*p*q), a*(b*(c*x**n)**p)**q)


def replacement2070(a, b, c, d, e, e1, n, p, u, x):
    return Dist(log(e*(b*e1/d)**n)**p, Int(u, x), x)


def replacement2071(a, b, c, d, e, e1, n, n1, n2, p, x):
    return -Dist(n*n1*p*(-a*d + b*c)/b, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/(c + d*x), x), x) + Simp((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/b, x)


def replacement2072(a, b, c, d, e, f, g, x):
    return Simp(PolyLog(S(2), Together(-a*e + c)/(c + d*x))/g, x)


def replacement2073(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(n*n1*p*(-a*d + b*c)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log((-a*d + b*c)/(b*(c + d*x)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log((-a*d + b*c)/(b*(c + d*x)))/g, x)


def replacement2074(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(n*n1*p*(-a*d + b*c)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log(-(-a*d + b*c)/(d*(a + b*x)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log(-(-a*d + b*c)/(d*(a + b*x)))/g, x)


def replacement2075(a, b, c, d, e, e1, f, g, n, n1, n2, x):
    return -Dist(n*n1*(-a*d + b*c)/g, Int(log(f + g*x)/((a + b*x)*(c + d*x)), x), x) + Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)*log(f + g*x)/g, x)


def replacement2076(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(d/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(c + d*x), x), x) - Dist((-c*g + d*f)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c + d*x)*(f + g*x)), x), x)


def replacement2077(a, b, c, d, e, f, g, x):
    return Simp(d**S(2)*LogIntegral(e*(a + b*x)/(c + d*x))/(e*g**S(2)*(-a*d + b*c)), x)


def replacement2078(a, b, c, d, e, e1, f, g, n, n1, n2, x):
    return Simp(d**S(2)*(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**(-S(1)/(n*n1))*(a + b*x)*ExpIntegralEi(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)/(n*n1))/(g**S(2)*n*n1*(c + d*x)*(-a*d + b*c)), x)


def replacement2079(a, b, c, d, e, e1, f, g, n, n1, n2, x):
    return Simp(b**S(2)*(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**(S(1)/(n*n1))*(c + d*x)*ExpIntegralEi(-log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)/(n*n1))/(g**S(2)*n*n1*(a + b*x)*(-a*d + b*c)), x)


def replacement2080(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return -Dist(n*n1*p*(-a*d + b*c)/(-a*g + b*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((c + d*x)*(f + g*x)), x), x) + Simp((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((f + g*x)*(-a*g + b*f)), x)


def replacement2081(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return -Dist(n*n1*p*(-a*d + b*c)/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((a + b*x)*(f + g*x)), x), x) + Simp((c + d*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((f + g*x)*(-c*g + d*f)), x)


def replacement2082(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(b/(-a*g + b*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(2), x), x) - Dist(g/(-a*g + b*f), Int((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(3), x), x)


def replacement2083(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(d/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(2), x), x) - Dist(g/(-c*g + d*f), Int((c + d*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(3), x), x)


def replacement2084(a, b, c, d, e, e1, f, g, m, n, n1, n2, p, x):
    return -Dist(n*n1*p*(-a*d + b*c)/(g*(m + S(1))), Int((f + g*x)**(m + S(1))*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp((f + g*x)**(m + S(1))*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(g*(m + S(1))), x)


def replacement2085(a, b, c, d, e, m, m2, n, p, u, x):
    return -Dist(n*p/(m + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(e*u**n)**(p + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(e*u**n)**p/((m + S(1))*(-a*d + b*c)), x)


def replacement2086(a, b, c, d, m, m2, p, u, x):
    return -Dist(p/(m + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(u)**(p + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(u)**p/((m + S(1))*(-a*d + b*c)), x)


def replacement2087(a, b, c, d, e, m, m2, n, u, x):
    return Simp((e*u**n)**(-(m + S(1))/n)*(a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*ExpIntegralEi((m + S(1))*log(e*u**n)/n)/(n*(-a*d + b*c)), x)


def replacement2088(a, b, c, d, m, m2, u, x):
    return Simp(u**(-m + S(-1))*(a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*ExpIntegralEi((m + S(1))*log(u))/(-a*d + b*c), x)


def replacement2089(a, b, c, d, e, m, m2, n, p, u, x):
    return -Dist((m + S(1))/(n*(p + S(1))), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(e*u**n)**(p + S(1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def replacement2090(a, b, c, d, m, m2, p, u, x):
    return -Dist((m + S(1))/(p + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(u)**(p + S(1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)


def replacement2091(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(d/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(c + d*x)**S(2), x), x)


def replacement2092(a, b, c, d, e, f, g, x):
    return Simp(PolyLog(S(2), -(f + g*x)*(a*e - c)/(f*(c + d*x)))/(-c*g + d*f), x)


def replacement2093(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(n*n1*p*(-a*d + b*c)/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log((f + g*x)*(-a*d + b*c)/((c + d*x)*(-a*g + b*f)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log((f + g*x)*(-a*d + b*c)/((c + d*x)*(-a*g + b*f)))/(-c*g + d*f), x)


def replacement2094(a, b, c, d, e, f, g, x):
    return Simp(c*PolyLog(S(2), -(c - d*x)*(a*e - c)/(c*(c + d*x)))/(S(2)*d*f), x)


def replacement2095(a, b, c, d, e, e1, f, g, h, n, n1, n2, p, x):
    return Dist(d**S(2), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c + d*x)*(-c*h + d*g + d*h*x)), x), x)


def replacement2096(a, b, c, d, e, e1, f, h, n, n1, n2, p, x):
    return -Dist(d**S(2)/h, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c - d*x)*(c + d*x)), x), x)


def replacement2097(a, b, c, d, e, e1, f, g, n, n1, n2, p, x):
    return Dist(b/(g*n*n1*(-a*d + b*c)), Subst(Int(x**p, x), x, log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)), x)


def replacement2098(a, b, c, d, e, n, p, u, v, x):
    return Dist(n*p, Int(PolyLog(S(2), Together(S(1) - v))*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(S(2), Together(S(1) - v))*log(e*u**n)**p/(-a*d + b*c), x)


def replacement2099(a, b, c, d, p, u, v, x):
    return Dist(p, Int(PolyLog(S(2), Together(S(1) - v))*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(S(2), Together(S(1) - v))*log(u)**p/(-a*d + b*c), x)


def With2100(a, b, c, d, e, n, p, u, v, x):
    f = (S(1) - v)/u
    return Dist(f/(n*(p + S(1))), Int(log(e*u**n)**(p + S(1))/((c + d*x)*(-a*f - b*f + c + d)), x), x) + Simp(log(v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def With2101(a, b, c, d, p, u, v, x):
    f = (S(1) - v)/u
    return Dist(f/(p + S(1)), Int(log(u)**(p + S(1))/((c + d*x)*(-a*f - b*f + c + d)), x), x) + Simp(log(u)**(p + S(1))*log(v)/((p + S(1))*(-a*d + b*c)), x)


def replacement2102(a, b, c, d, e, n, p, u, v, x):
    return -Dist(n*p, Int(PolyLog(S(2), Together(S(1) - v))*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(S(2), Together(S(1) - v))*log(e*u**n)**p/(-a*d + b*c), x)


def replacement2103(a, b, c, d, p, u, v, x):
    return -Dist(p, Int(PolyLog(S(2), Together(S(1) - v))*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(S(2), Together(S(1) - v))*log(u)**p/(-a*d + b*c), x)


def With2104(a, b, c, d, e, n, p, u, v, x):
    f = u*(S(1) - v)
    return -Dist(f/(n*(p + S(1))), Int(log(e*u**n)**(p + S(1))/((a + b*x)*(a + b - c*f - d*f)), x), x) + Simp(log(v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def With2105(a, b, c, d, p, u, v, x):
    f = u*(S(1) - v)
    return -Dist(f/(p + S(1)), Int(log(u)**(p + S(1))/((a + b*x)*(a + b - c*f - d*f)), x), x) + Simp(log(u)**(p + S(1))*log(v)/((p + S(1))*(-a*d + b*c)), x)


def replacement2106(a, b, c, d, e, n, p, q, u, v, x):
    return -Dist(n*p, Int(PolyLog(q + S(1), v)*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q + S(1), v)*log(e*u**n)**p/(-a*d + b*c), x)


def replacement2107(a, b, c, d, p, q, u, v, x):
    return -Dist(p, Int(PolyLog(q + S(1), v)*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q + S(1), v)*log(u)**p/(-a*d + b*c), x)


def replacement2108(a, b, c, d, e, n, p, q, u, v, x):
    return -Dist(S(1)/(n*(p + S(1))), Int(PolyLog(q + S(-1), v)*log(e*u**n)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def replacement2109(a, b, c, d, p, q, u, v, x):
    return -Dist(S(1)/(p + S(1)), Int(PolyLog(q + S(-1), v)*log(u)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)


def replacement2110(a, b, c, d, e, n, p, q, u, v, x):
    return Dist(n*p, Int(PolyLog(q + S(1), v)*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(q + S(1), v)*log(e*u**n)**p/(-a*d + b*c), x)


def replacement2111(a, b, c, d, p, q, u, v, x):
    return Dist(p, Int(PolyLog(q + S(1), v)*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(q + S(1), v)*log(u)**p/(-a*d + b*c), x)


def replacement2112(a, b, c, d, e, n, p, q, u, v, x):
    return Dist(S(1)/(n*(p + S(1))), Int(PolyLog(q + S(-1), v)*log(e*u**n)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)


def replacement2113(a, b, c, d, p, q, u, v, x):
    return Dist(S(1)/(p + S(1)), Int(PolyLog(q + S(-1), v)*log(u)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)


def replacement2114(a, b, c, d, e, e1, f, g, h, n, n1, n2, p, u, x):
    return Dist(b*d/h, Int(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((a + b*x)*(c + d*x)), x), x)


def replacement2115(a, b, c, d, e, e1, f, h, n, n1, n2, p, u, x):
    return Dist(b*d/h, Int(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((a + b*x)*(c + d*x)), x), x)


def With2116(a, b, c, d, e, e1, f, g, h, n, n1, n2, x):
    u = IntHide(S(1)/(f + g*x + h*x**S(2)), x)
    return -Dist(n*(-a*d + b*c), Int(u/((a + b*x)*(c + d*x)), x), x) + Simp(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n), x)


def With2117(a, b, c, d, e, e1, f, h, n, n1, n2, x):
    u = IntHide(S(1)/(f + h*x**S(2)), x)
    return -Dist(n*(-a*d + b*c), Int(u/((a + b*x)*(c + d*x)), x), x) + Simp(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n), x)


def With2118(RFx, a, b, c, d, e, e1, n, n1, n2, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**p, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement2118(RFx, a, b, c, d, e, e1, n, n1, n2, p, x):

    u = ExpandIntegrand(log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**p, RFx, x)
    return Int(u, x)


def With2119(p, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    lst = QuotientOfLinearsParts(v, x)
    if Not(And(OneQ(p), ZeroQ(Part(lst, S(3))))):
        return True
    return False


def replacement2119(p, u, v, x):

    lst = QuotientOfLinearsParts(v, x)
    return Int(u*log((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**p, x)


def replacement2120(a, b, c, n, p, x):
    return -Dist(b*n*p, Int(x**n/(a + b*x**n), x), x) + Simp(x*log(c*(a + b*x**n)**p), x)


def replacement2121(c, p, v, x):
    return Int(log(c*ExpandToSum(v, x)**p), x)


def replacement2122(a, b, c, d, e, f, g, n, p, x):
    return -Dist(b*e*n*p/g, Int(x**(n + S(-1))*log(f + g*x)/(d + e*x**n), x), x) + Simp((a + b*log(c*(d + e*x**n)**p))*log(f + g*x)/g, x)


def replacement2123(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(b*e*n*p/(g*(m + S(1))), Int(x**(n + S(-1))*(f + g*x)**(m + S(1))/(d + e*x**n), x), x) + Simp((a + b*log(c*(d + e*x**n)**p))*(f + g*x)**(m + S(1))/(g*(m + S(1))), x)


def replacement2124(a, b, c, m, p, u, v, x):
    return Int((a + b*log(c*ExpandToSum(v, x)**p))*ExpandToSum(u, x)**m, x)


def With2125(a, b, c, d, e, f, g, m, n, p, x):
    w = IntHide(asin(f + g*x)**m, x)
    return -Dist(b*e*n*p, Int(SimplifyIntegrand(w*x**(n + S(-1))/(d + e*x**n), x), x), x) + Dist(a + b*log(c*(d + e*x**n)**p), w, x)


def With2126(a, b, c, d, e, f, g, p, x):
    u = IntHide(S(1)/(f + g*x**S(2)), x)
    return -Dist(S(2)*b*e*p, Int(u*x/(d + e*x**S(2)), x), x) + Simp(u*(a + b*log(c*(d + e*x**S(2))**p)), x)


def replacement2127(a, b, c, d, e, n, p, x):
    return -Dist(S(2)*b*e*n*p, Int(x**S(2)*(a + b*log(c*(d + e*x**S(2))**p))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x*(a + b*log(c*(d + e*x**S(2))**p))**n, x)


def replacement2128(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/2, Subst(Int(x**(m/S(2) + S(-1)/2)*(a + b*log(c*(d + e*x)**p))**n, x), x, x**S(2)), x)


def replacement2129(a, b, c, d, e, m, n, p, x):
    return -Dist(S(2)*b*e*n*p/(m + S(1)), Int(x**(m + S(2))*(a + b*log(c*(d + e*x**S(2))**p))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*log(c*(d + e*x**S(2))**p))**n/(m + S(1)), x)


def With2130(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        w = DerivativeDivides(v, u*(S(1) - v), x)
        res = Not(FalseQ(w))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement2130(u, v, x):

    w = DerivativeDivides(v, u*(S(1) - v), x)
    return Simp(w*PolyLog(S(2), Together(S(1) - v)), x)


def With2131(a, b, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        z = DerivativeDivides(v, w*(S(1) - v), x)
        res = Not(FalseQ(z))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement2131(a, b, u, v, w, x):

    z = DerivativeDivides(v, w*(S(1) - v), x)
    return -Dist(b, Int(SimplifyIntegrand(z*D(u, x)*PolyLog(S(2), Together(S(1) - v))/u, x), x), x) + Simp(z*(a + b*log(u))*PolyLog(S(2), Together(S(1) - v)), x)


def replacement2132(a, b, c, d, e, n, p, x):
    return -Dist(b*n*p, Int(S(1)/(a*(d + e*x)**(-n) + b), x), x) + Simp((d + e*x)*log(c*(a + b*(d + e*x)**n)**p)/e, x)


def replacement2133(a, b, c, d, e, n, p, x):
    return Dist(a*n*p, Int(S(1)/(a + b*(d + e*x)**n), x), x) + Simp((d + e*x)*log(c*(a + b*(d + e*x)**n)**p)/e, x) - Simp(n*p*x, x)


def replacement2134(a, b, c, d, e, f, g, n, p, x):
    return -Dist(b*e*n*p/(d*g), Subst(Int((a + b*log(c*(d + e*x)**p))**(n + S(-1))/x, x), x, S(1)/(f + g*x)), x) + Simp((a + b*log(c*(d + e/(f + g*x))**p))**n*(d*(f + g*x) + e)/(d*g), x)


def replacement2135(RFx, a, b, c, n, p, x):
    return -Dist(b*n*p, Int(SimplifyIntegrand(x*(a + b*log(RFx**p*c))**(n + S(-1))*D(RFx, x)/RFx, x), x), x) + Simp(x*(a + b*log(RFx**p*c))**n, x)


def replacement2136(RFx, a, b, c, d, e, n, p, x):
    return -Dist(b*n*p/e, Int((a + b*log(RFx**p*c))**(n + S(-1))*D(RFx, x)*log(d + e*x)/RFx, x), x) + Simp((a + b*log(RFx**p*c))**n*log(d + e*x)/e, x)


def replacement2137(RFx, a, b, c, d, e, m, n, p, x):
    return -Dist(b*n*p/(e*(m + S(1))), Int(SimplifyIntegrand((a + b*log(RFx**p*c))**(n + S(-1))*(d + e*x)**(m + S(1))*D(RFx, x)/RFx, x), x), x) + Simp((a + b*log(RFx**p*c))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def With2138(RFx, c, d, e, n, x):
    u = IntHide(S(1)/(d + e*x**S(2)), x)
    return -Dist(n, Int(SimplifyIntegrand(u*D(RFx, x)/RFx, x), x), x) + Simp(u*log(RFx**n*c), x)


def With2139(Px, Qx, c, n, x):
    u = IntHide(S(1)/Qx, x)
    return -Dist(n, Int(SimplifyIntegrand(u*D(Px, x)/Px, x), x), x) + Simp(u*log(Px**n*c), x)


def With2140(RFx, RGx, a, b, c, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((a + b*log(RFx**p*c))**n, RGx, x)
    if SumQ(u):
        return True
    return False


def replacement2140(RFx, RGx, a, b, c, n, p, x):

    u = ExpandIntegrand((a + b*log(RFx**p*c))**n, RGx, x)
    return Int(u, x)


def With2141(RFx, RGx, a, b, c, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(RGx*(a + b*log(RFx**p*c))**n, x)
    if SumQ(u):
        return True
    return False


def replacement2141(RFx, RGx, a, b, c, n, p, x):

    u = ExpandIntegrand(RGx*(a + b*log(RFx**p*c))**n, x)
    return Int(u, x)


def With2142(RFx, a, b, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = SubstForFractionalPowerOfLinear(RFx*(a + b*log(u)), x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement2142(RFx, a, b, u, x):

    lst = SubstForFractionalPowerOfLinear(RFx*(a + b*log(u)), x)
    return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)


def replacement2143(F, a, b, c, e, f, g, m, n, x):
    return Dist(g*m/(b*c*n*log(F)), Int((f + g*x)**(m + S(-1))*PolyLog(S(2), -e*(F**(c*(a + b*x)))**n), x), x) - Simp((f + g*x)**m*PolyLog(S(2), -e*(F**(c*(a + b*x)))**n)/(b*c*n*log(F)), x)


def replacement2144(F, a, b, c, d, e, f, g, m, n, x):
    return Int((f + g*x)**m*log(S(1) + e*(F**(c*(a + b*x)))**n/d), x) - Simp((f + g*x)**(m + S(1))*log(S(1) + e*(F**(c*(a + b*x)))**n/d)/(g*(m + S(1))), x) + Simp((f + g*x)**(m + S(1))*log(d + e*(F**(c*(a + b*x)))**n)/(g*(m + S(1))), x)


def replacement2145(a, b, c, d, e, f, x):
    return Dist(f**S(2)*(-S(4)*a*c + b**S(2))/S(2), Int(x/(-f*sqrt(a + b*x + c*x**S(2))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d)) + (-b*f**S(2) + S(2)*d*e)*(a + b*x + c*x**S(2))), x), x) + Simp(x*log(d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)


def replacement2146(a, c, d, e, f, x):
    return -Dist(a*c*f**S(2), Int(x/(d*e*(a + c*x**S(2)) + f*sqrt(a + c*x**S(2))*(a*e - c*d*x)), x), x) + Simp(x*log(d + e*x + f*sqrt(a + c*x**S(2))), x)


def replacement2147(a, b, c, d, e, f, g, m, x):
    return Dist(f**S(2)*(-S(4)*a*c + b**S(2))/(S(2)*g*(m + S(1))), Int((g*x)**(m + S(1))/(-f*sqrt(a + b*x + c*x**S(2))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d)) + (-b*f**S(2) + S(2)*d*e)*(a + b*x + c*x**S(2))), x), x) + Simp((g*x)**(m + S(1))*log(d + e*x + f*sqrt(a + b*x + c*x**S(2)))/(g*(m + S(1))), x)


def replacement2148(a, c, d, e, f, g, m, x):
    return -Dist(a*c*f**S(2)/(g*(m + S(1))), Int((g*x)**(m + S(1))/(d*e*(a + c*x**S(2)) + f*sqrt(a + c*x**S(2))*(a*e - c*d*x)), x), x) + Simp((g*x)**(m + S(1))*log(d + e*x + f*sqrt(a + c*x**S(2)))/(g*(m + S(1))), x)


def replacement2149(d, e, f, u, v, x):
    return Int(v*log(d + e*x + f*sqrt(ExpandToSum(u, x))), x)


def replacement2150(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/u, x), x) + Simp(x*log(u), x)


def replacement2151(a, b, u, x):
    return -Dist(S(1)/b, Int(SimplifyIntegrand(D(u, x)*log(a + b*x)/u, x), x), x) + Simp(log(u)*log(a + b*x)/b, x)


def replacement2152(a, b, m, u, x):
    return -Dist(S(1)/(b*(m + S(1))), Int(SimplifyIntegrand((a + b*x)**(m + S(1))*D(u, x)/u, x), x), x) + Simp((a + b*x)**(m + S(1))*log(u)/(b*(m + S(1))), x)


def With2153(Qx, u, x):
    v = IntHide(S(1)/Qx, x)
    return -Int(SimplifyIntegrand(v*D(u, x)/u, x), x) + Simp(v*log(u), x)


def replacement2154(a, u, x):
    return -Int(SimplifyIntegrand(u**(a*x + S(-1))*x*D(u, x), x), x) + Simp(u**(a*x)/a, x)


def With2155(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement2155(u, v, x):

    w = IntHide(v, x)
    return Dist(log(u), w, x) - Int(SimplifyIntegrand(w*D(u, x)/u, x), x)


def replacement2156(v, w, x):
    return -Int(SimplifyIntegrand(x*D(v, x)*log(w)/v, x), x) - Int(SimplifyIntegrand(x*D(w, x)*log(v)/w, x), x) + Simp(x*log(v)*log(w), x)


def With2157(u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    z = IntHide(u, x)
    if InverseFunctionFreeQ(z, x):
        return True
    return False


def replacement2157(u, v, w, x):

    z = IntHide(u, x)
    return Dist(log(v)*log(w), z, x) - Int(SimplifyIntegrand(z*D(v, x)*log(w)/v, x), x) - Int(SimplifyIntegrand(z*D(w, x)*log(v)/w, x), x)


def replacement2158(a, b, n, p, x):
    return -Dist(n*p, Int(S(1)/log(b*x**n), x), x) + Simp(x*log(a*log(b*x**n)**p), x)


def replacement2159(a, b, n, p, x):
    return Simp((-p + log(a*log(b*x**n)**p))*log(b*x**n)/n, x)


def replacement2160(a, b, m, n, p, x):
    return -Dist(n*p/(m + S(1)), Int(x**m/log(b*x**n), x), x) + Simp(x**(m + S(1))*log(a*log(b*x**n)**p)/(m + S(1)), x)


def replacement2161(A, B, a, b, c, d, x):
    return Dist(B/b, Int(sqrt(a + b*log(c + d*x)), x), x) + Dist((A*b - B*a)/b, Int(S(1)/sqrt(a + b*log(c + d*x)), x), x)


def replacement2162(a, f, u, x):
    return Int(u**(a*log(f)), x)


def With2163(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = FunctionOfLog(u*x, x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement2163(u, x):

    lst = FunctionOfLog(u*x, x)
    return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, log(Part(lst, S(2)))), x)


def replacement2164(u, v, x):
    return Dist(-LogGamma(v) + log(Gamma(v)), Int(u, x), x) + Int(u*LogGamma(v), x)


def replacement2165(a, b, n, p, u, v, w, x):
    return Int(u*w**p*(a + b*log(v)**n)**p, x)


def replacement2166(a, b, c, d, e, f, n, p, q, u, x):
    return Int(u*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
