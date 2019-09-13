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

def miscellaneous_algebraic(rubi):
    from sympy.integrals.rubi.constraints import cons800, cons2, cons3, cons8, cons52, cons4, cons5, cons20, cons19, cons801, cons29, cons50, cons127, cons54, cons802, cons27, cons803, cons804, cons151, cons805, cons502, cons806, cons650, cons807, cons808, cons21, cons48, cons809, cons810, cons70, cons811, cons812, cons813, cons814, cons815, cons816, cons817, cons818, cons819, cons820, cons821, cons822, cons823, cons454, cons824, cons825, cons826, cons827, cons828, cons829, cons830, cons831, cons832, cons833, cons834, cons835, cons836, cons837, cons838, cons839, cons840, cons841, cons842, cons843, cons844, cons845, cons846, cons847, cons848, cons849, cons850, cons851, cons852, cons853, cons854, cons210, cons211, cons66, cons855, cons68, cons856, cons857, cons466, cons858, cons859, cons860, cons55, cons13, cons139, cons861, cons862, cons150, cons246, cons165, cons863, cons523, cons864, cons865, cons866, cons86, cons867, cons36, cons37, cons868, cons470, cons471, cons869, cons870, cons38, cons871, cons872, cons873, cons874, cons875, cons876, cons877, cons878, cons879, cons880, cons881, cons882, cons883, cons884, cons885, cons886, cons887, cons888, cons889, cons890, cons891, cons892, cons893, cons894, cons895, cons896, cons897, cons898, cons899, cons900, cons901, cons902, cons903, cons904, cons905, cons906, cons676, cons907, cons483, cons908, cons909, cons484, cons910, cons911, cons912, cons913, cons914, cons915, cons916, cons917, cons918, cons87, cons33, cons96, cons919, cons198, cons369, cons358, cons491, cons543, cons25, cons920, cons556, cons921, cons554, cons57, cons496, cons59, cons60, cons61, cons62, cons922, cons923, cons924, cons925, cons926, cons597, cons73, cons927, cons588, cons89, cons130, cons928, cons929, cons930, cons931, cons932, cons47, cons316, cons228, cons933, cons934, cons935, cons936, cons937, cons938, cons939, cons940, cons941, cons942, cons943, cons944, cons945, cons946, cons947, cons948, cons284, cons949, cons65, cons721, cons950, cons951, cons952, cons75, cons953, cons704, cons149, cons954, cons955, cons798, cons956, cons957, cons958, cons959, cons960, cons961, cons962, cons963, cons964, cons965, cons966, cons967, cons968, cons71, cons969, cons970, cons971, cons972, cons973, cons974, cons975, cons976, cons977, cons514, cons978, cons979, cons980, cons981, cons982, cons669, cons983, cons984, cons799, cons985, cons986, cons987, cons988, cons989, cons990, cons95, cons90, cons991, cons992, cons993, cons994, cons995, cons996, cons997, cons998, cons999, cons1000, cons40, cons1001, cons1002, cons1003, cons1004, cons1005, cons1006, cons1007, cons1008, cons1009, cons1010, cons1011, cons1012, cons385, cons1013, cons1014, cons1015, cons1016, cons1017, cons1018, cons1019, cons1020, cons359, cons1021, cons248, cons1022, cons1023, cons1024, cons1025, cons1026, cons1027, cons1028, cons1029, cons1030, cons1031, cons1032, cons1033, cons1034, cons1035, cons1036, cons1037, cons1038, cons1039, cons1040, cons1041, cons1042, cons1043, cons1044, cons1045, cons299, cons1046, cons1047, cons1048, cons1049, cons1050, cons707, cons384, cons1051, cons1052, cons699, cons711, cons155, cons1053, cons1054, cons1055, cons1056, cons1057, cons1058, cons1059, cons1060, cons1061, cons226, cons1062, cons517, cons1063, cons1064, cons1065, cons1066, cons1067, cons1068, cons1069, cons1070, cons1071, cons1072, cons1073, cons45, cons481, cons482, cons1074, cons1075, cons1076, cons1077, cons1078, cons1079, cons1080, cons1081, cons1082, cons1083, cons1084, cons1085, cons1086, cons1087, cons1088, cons1089, cons1090, cons1091

    pattern1476 = Pattern(Integral(((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons52, cons4, cons5, cons800)
    rule1476 = ReplacementRule(pattern1476, replacement1476)
    pattern1477 = Pattern(Integral(x_**WC('m', S(1))*((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons52, cons800, cons20)
    rule1477 = ReplacementRule(pattern1477, replacement1477)
    pattern1478 = Pattern(Integral(x_**WC('m', S(1))*((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('r', S(1))*WC('e', S(1)))**p_*((c_ + x_**WC('n', S(1))*WC('d', S(1)))**s_*WC('f', S(1)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons54, cons802, cons801)
    rule1478 = ReplacementRule(pattern1478, replacement1478)
    pattern1479 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons27)
    rule1479 = ReplacementRule(pattern1479, replacement1479)
    pattern1480 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons803, cons804)
    rule1480 = ReplacementRule(pattern1480, replacement1480)

    pattern1481 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons151, cons805)
    rule1481 = ReplacementRule(pattern1481, With1481)

    pattern1482 = Pattern(Integral(x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons151, cons502)
    rule1482 = ReplacementRule(pattern1482, With1482)

    pattern1483 = Pattern(Integral(u_**WC('r', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons806, cons151, cons805, cons650)
    rule1483 = ReplacementRule(pattern1483, With1483)

    pattern1484 = Pattern(Integral(u_**WC('r', S(1))*x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons806, cons151, cons805, cons807)
    rule1484 = ReplacementRule(pattern1484, With1484)
    pattern1485 = Pattern(Integral(((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons808)
    rule1485 = ReplacementRule(pattern1485, replacement1485)
    pattern1486 = Pattern(Integral(x_**WC('m', S(1))*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons20)
    rule1486 = ReplacementRule(pattern1486, replacement1486)
    pattern1487 = Pattern(Integral((x_*WC('d', S(1)))**m_*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons21)
    rule1487 = ReplacementRule(pattern1487, replacement1487)
    pattern1488 = Pattern(Integral(((WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons48)
    rule1488 = ReplacementRule(pattern1488, replacement1488)
    pattern1489 = Pattern(Integral(x_**WC('m', S(1))*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons48, cons20)
    rule1489 = ReplacementRule(pattern1489, replacement1489)
    pattern1490 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons21)
    rule1490 = ReplacementRule(pattern1490, replacement1490)
    pattern1491 = Pattern(Integral((x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons809, cons810)
    rule1491 = ReplacementRule(pattern1491, replacement1491)
    pattern1492 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons809, cons810, cons20)
    rule1492 = ReplacementRule(pattern1492, replacement1492)
    pattern1493 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons809, cons21, cons810)
    rule1493 = ReplacementRule(pattern1493, replacement1493)
    pattern1494 = Pattern(Integral(u_**m_, x_), cons19, cons70, cons811)
    rule1494 = ReplacementRule(pattern1494, replacement1494)
    pattern1495 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1)), x_), cons19, cons4, cons812, cons813)
    rule1495 = ReplacementRule(pattern1495, replacement1495)
    pattern1496 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), cons19, cons4, cons5, cons814, cons815)
    rule1496 = ReplacementRule(pattern1496, replacement1496)
    pattern1497 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1))*z_**WC('q', S(1)), x_), cons19, cons4, cons5, cons52, cons816, cons817)
    rule1497 = ReplacementRule(pattern1497, replacement1497)
    pattern1498 = Pattern(Integral(u_**p_, x_), cons5, cons818, cons819)
    rule1498 = ReplacementRule(pattern1498, replacement1498)
    pattern1499 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), cons19, cons5, cons70, cons820, cons821)
    rule1499 = ReplacementRule(pattern1499, replacement1499)
    pattern1500 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), cons19, cons4, cons5, cons812, cons822, cons823)
    rule1500 = ReplacementRule(pattern1500, replacement1500)
    pattern1501 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons52, cons454, cons824)
    rule1501 = ReplacementRule(pattern1501, replacement1501)
    pattern1502 = Pattern(Integral(u_**p_, x_), cons5, cons825, cons826)
    rule1502 = ReplacementRule(pattern1502, replacement1502)
    pattern1503 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('c', S(1)))**WC('m', S(1)), x_), cons8, cons19, cons5, cons825, cons826)
    rule1503 = ReplacementRule(pattern1503, replacement1503)
    pattern1504 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons52, cons827, cons828, cons829)
    rule1504 = ReplacementRule(pattern1504, replacement1504)
    pattern1505 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), cons19, cons5, cons52, cons827, cons828, cons829)
    rule1505 = ReplacementRule(pattern1505, replacement1505)
    pattern1506 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1))*w_**WC('q', S(1)), x_), cons19, cons5, cons52, cons830, cons828, cons831, cons832)
    rule1506 = ReplacementRule(pattern1506, replacement1506)
    pattern1507 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1))*z_**WC('r', S(1)), x_), cons19, cons5, cons52, cons54, cons833, cons828, cons834, cons835)
    rule1507 = ReplacementRule(pattern1507, replacement1507)
    pattern1508 = Pattern(Integral(u_**p_, x_), cons5, cons836, cons837)
    rule1508 = ReplacementRule(pattern1508, replacement1508)
    pattern1509 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), cons19, cons5, cons836, cons837)
    rule1509 = ReplacementRule(pattern1509, replacement1509)
    pattern1510 = Pattern(Integral(u_**p_, x_), cons5, cons838, cons839)
    rule1510 = ReplacementRule(pattern1510, replacement1510)
    pattern1511 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('d', S(1)))**WC('m', S(1)), x_), cons29, cons19, cons5, cons838, cons839)
    rule1511 = ReplacementRule(pattern1511, replacement1511)
    pattern1512 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), cons5, cons52, cons825, cons840, cons841)
    rule1512 = ReplacementRule(pattern1512, replacement1512)
    pattern1513 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), cons5, cons52, cons825, cons842, cons843)
    rule1513 = ReplacementRule(pattern1513, replacement1513)
    pattern1514 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), cons19, cons5, cons52, cons844, cons838, cons845)
    rule1514 = ReplacementRule(pattern1514, replacement1514)
    pattern1515 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), cons19, cons5, cons52, cons844, cons825, cons846)
    rule1515 = ReplacementRule(pattern1515, replacement1515)
    pattern1516 = Pattern(Integral(u_**p_, x_), cons5, cons847, cons848)
    rule1516 = ReplacementRule(pattern1516, replacement1516)
    pattern1517 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), cons19, cons5, cons847, cons848)
    rule1517 = ReplacementRule(pattern1517, replacement1517)
    pattern1518 = Pattern(Integral(u_**WC('p', S(1))*z_, x_), cons5, cons844, cons847, cons849, cons850)
    rule1518 = ReplacementRule(pattern1518, replacement1518)
    pattern1519 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_, x_), cons19, cons5, cons844, cons847, cons849, cons850)
    rule1519 = ReplacementRule(pattern1519, replacement1519)
    pattern1520 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons8, cons50, cons127, cons210, cons211, cons19, cons4, cons851, cons852, cons853, cons854)
    rule1520 = ReplacementRule(pattern1520, replacement1520)
    pattern1521 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons853, cons851, cons852, cons854)
    rule1521 = ReplacementRule(pattern1521, replacement1521)

    pattern1522 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons66, cons151, cons855)
    rule1522 = ReplacementRule(pattern1522, With1522)
    pattern1523 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons68, cons856, cons857)
    rule1523 = ReplacementRule(pattern1523, replacement1523)
    pattern1524 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons466, cons858)
    rule1524 = ReplacementRule(pattern1524, replacement1524)
    pattern1525 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons66, cons859)
    rule1525 = ReplacementRule(pattern1525, replacement1525)
    pattern1526 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons66, cons859)
    rule1526 = ReplacementRule(pattern1526, replacement1526)
    pattern1527 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons860, cons502)
    rule1527 = ReplacementRule(pattern1527, replacement1527)
    pattern1528 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons860, cons502)
    rule1528 = ReplacementRule(pattern1528, replacement1528)
    pattern1529 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons66, cons55, cons13, cons139)
    rule1529 = ReplacementRule(pattern1529, replacement1529)
    pattern1530 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons29, cons19, cons4, cons5, cons66, cons861)
    rule1530 = ReplacementRule(pattern1530, replacement1530)
    pattern1531 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons66, cons861, cons862)
    rule1531 = ReplacementRule(pattern1531, replacement1531)

    pattern1532 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons150, cons246, cons165, cons863)
    rule1532 = ReplacementRule(pattern1532, With1532)

    pattern1533 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons66, cons523, cons13, cons165)
    rule1533 = ReplacementRule(pattern1533, With1533)

    pattern1534 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons523, cons13, cons165)
    rule1534 = ReplacementRule(pattern1534, With1534)

    pattern1535 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons150, cons13, cons139, CustomConstraint(With1535))
    rule1535 = ReplacementRule(pattern1535, replacement1535)
    pattern1536 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons150, cons13, cons139, cons864)
    rule1536 = ReplacementRule(pattern1536, replacement1536)
    pattern1537 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons865)
    rule1537 = ReplacementRule(pattern1537, replacement1537)
    pattern1538 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons29, cons127, cons210, cons865)
    rule1538 = ReplacementRule(pattern1538, replacement1538)
    pattern1539 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons210, cons865)
    rule1539 = ReplacementRule(pattern1539, replacement1539)
    pattern1540 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons50, cons127, cons211, cons866)
    rule1540 = ReplacementRule(pattern1540, replacement1540)
    pattern1541 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons50, cons211, cons866)
    rule1541 = ReplacementRule(pattern1541, replacement1541)
    pattern1542 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons211, cons866, cons865)
    rule1542 = ReplacementRule(pattern1542, replacement1542)
    pattern1543 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons210, cons211, cons866, cons865)
    rule1543 = ReplacementRule(pattern1543, replacement1543)

    pattern1544 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons150, cons13, cons139, CustomConstraint(With1544))
    rule1544 = ReplacementRule(pattern1544, replacement1544)

    pattern1545 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons66, cons150, cons13, cons139, cons86)
    rule1545 = ReplacementRule(pattern1545, With1545)

    pattern1546 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons860, cons150, cons20, CustomConstraint(With1546))
    rule1546 = ReplacementRule(pattern1546, replacement1546)
    pattern1547 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons867)
    rule1547 = ReplacementRule(pattern1547, replacement1547)

    pattern1548 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons868, cons470)
    rule1548 = ReplacementRule(pattern1548, With1548)

    pattern1549 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons868, cons471)
    rule1549 = ReplacementRule(pattern1549, With1549)
    pattern1550 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons869, cons870)
    rule1550 = ReplacementRule(pattern1550, replacement1550)

    pattern1551 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons871)
    rule1551 = ReplacementRule(pattern1551, With1551)

    pattern1552 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons872)
    rule1552 = ReplacementRule(pattern1552, With1552)

    pattern1553 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons873)
    rule1553 = ReplacementRule(pattern1553, With1553)

    pattern1554 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons874)
    rule1554 = ReplacementRule(pattern1554, With1554)

    pattern1555 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons875)
    rule1555 = ReplacementRule(pattern1555, With1555)

    pattern1556 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons876)
    rule1556 = ReplacementRule(pattern1556, With1556)

    pattern1557 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons877)
    rule1557 = ReplacementRule(pattern1557, With1557)

    pattern1558 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons878)
    rule1558 = ReplacementRule(pattern1558, With1558)

    pattern1559 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons879)
    rule1559 = ReplacementRule(pattern1559, With1559)

    pattern1560 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons880)
    rule1560 = ReplacementRule(pattern1560, With1560)

    pattern1561 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons881)
    rule1561 = ReplacementRule(pattern1561, With1561)

    pattern1562 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons882)
    rule1562 = ReplacementRule(pattern1562, With1562)

    pattern1563 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons883)
    rule1563 = ReplacementRule(pattern1563, With1563)

    pattern1564 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons884)
    rule1564 = ReplacementRule(pattern1564, With1564)

    pattern1565 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons885)
    rule1565 = ReplacementRule(pattern1565, With1565)

    pattern1566 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons886)
    rule1566 = ReplacementRule(pattern1566, With1566)

    pattern1567 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons887)
    rule1567 = ReplacementRule(pattern1567, With1567)

    pattern1568 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons888)
    rule1568 = ReplacementRule(pattern1568, With1568)

    pattern1569 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons889)
    rule1569 = ReplacementRule(pattern1569, With1569)

    pattern1570 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons890)
    rule1570 = ReplacementRule(pattern1570, With1570)

    pattern1571 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons891)
    rule1571 = ReplacementRule(pattern1571, With1571)

    pattern1572 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons892)
    rule1572 = ReplacementRule(pattern1572, With1572)

    pattern1573 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons893)
    rule1573 = ReplacementRule(pattern1573, With1573)

    pattern1574 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons894)
    rule1574 = ReplacementRule(pattern1574, With1574)
    pattern1575 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons895)
    rule1575 = ReplacementRule(pattern1575, replacement1575)
    pattern1576 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons896)
    rule1576 = ReplacementRule(pattern1576, replacement1576)
    pattern1577 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons897)
    rule1577 = ReplacementRule(pattern1577, replacement1577)

    pattern1578 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons898)
    rule1578 = ReplacementRule(pattern1578, With1578)

    pattern1579 = Pattern(Integral(x_*(x_*WC('C', S(1)) + WC('B', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons899)
    rule1579 = ReplacementRule(pattern1579, With1579)

    pattern1580 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons900)
    rule1580 = ReplacementRule(pattern1580, With1580)

    pattern1581 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons901)
    rule1581 = ReplacementRule(pattern1581, With1581)

    pattern1582 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons902)
    rule1582 = ReplacementRule(pattern1582, With1582)

    pattern1583 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons903)
    rule1583 = ReplacementRule(pattern1583, With1583)

    pattern1584 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons868, cons904, cons905, CustomConstraint(With1584))
    rule1584 = ReplacementRule(pattern1584, replacement1584)

    pattern1585 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons904, cons905, CustomConstraint(With1585))
    rule1585 = ReplacementRule(pattern1585, replacement1585)

    pattern1586 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons904, cons905, CustomConstraint(With1586))
    rule1586 = ReplacementRule(pattern1586, replacement1586)

    pattern1587 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons868, cons904, cons906, CustomConstraint(With1587))
    rule1587 = ReplacementRule(pattern1587, replacement1587)

    pattern1588 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons37, cons38, cons904, cons906, CustomConstraint(With1588))
    rule1588 = ReplacementRule(pattern1588, replacement1588)

    pattern1589 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons36, cons38, cons904, cons906, CustomConstraint(With1589))
    rule1589 = ReplacementRule(pattern1589, replacement1589)

    pattern1590 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons19, cons66, cons676, cons907, CustomConstraint(With1590))
    rule1590 = ReplacementRule(pattern1590, replacement1590)

    pattern1591 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons66, cons676, cons907, CustomConstraint(With1591))
    rule1591 = ReplacementRule(pattern1591, replacement1591)

    pattern1592 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons483, cons908)
    rule1592 = ReplacementRule(pattern1592, With1592)

    pattern1593 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons483, cons909)
    rule1593 = ReplacementRule(pattern1593, With1593)

    pattern1594 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons484, cons910)
    rule1594 = ReplacementRule(pattern1594, With1594)

    pattern1595 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons484, cons911)
    rule1595 = ReplacementRule(pattern1595, With1595)

    pattern1596 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons912)
    rule1596 = ReplacementRule(pattern1596, With1596)

    pattern1597 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons913)
    rule1597 = ReplacementRule(pattern1597, With1597)
    pattern1598 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons914)
    rule1598 = ReplacementRule(pattern1598, replacement1598)
    pattern1599 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons915)
    rule1599 = ReplacementRule(pattern1599, replacement1599)
    pattern1600 = Pattern(Integral(Pq_/(x_*sqrt(a_ + x_**n_*WC('b', S(1)))), x_), cons2, cons3, cons66, cons150, cons916)
    rule1600 = ReplacementRule(pattern1600, replacement1600)

    pattern1601 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons66, cons676, cons917)
    rule1601 = ReplacementRule(pattern1601, With1601)

    pattern1602 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons66, cons676, cons917)
    rule1602 = ReplacementRule(pattern1602, With1602)
    pattern1603 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons66, cons150, cons918)
    rule1603 = ReplacementRule(pattern1603, replacement1603)
    pattern1604 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons19, cons66, cons87)
    rule1604 = ReplacementRule(pattern1604, replacement1604)
    pattern1605 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons66, cons87)
    rule1605 = ReplacementRule(pattern1605, replacement1605)

    pattern1606 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons66, cons150, cons33, cons96, cons919, CustomConstraint(With1606))
    rule1606 = ReplacementRule(pattern1606, replacement1606)

    pattern1607 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons66, cons150, CustomConstraint(With1607))
    rule1607 = ReplacementRule(pattern1607, replacement1607)

    pattern1608 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons66, cons150, CustomConstraint(With1608))
    rule1608 = ReplacementRule(pattern1608, replacement1608)

    pattern1609 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons66, cons198, cons20)
    rule1609 = ReplacementRule(pattern1609, With1609)

    pattern1610 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons66, cons198, cons369)
    rule1610 = ReplacementRule(pattern1610, With1610)

    pattern1611 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons66, cons198, cons358)
    rule1611 = ReplacementRule(pattern1611, With1611)

    pattern1612 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons5, cons66, cons491)
    rule1612 = ReplacementRule(pattern1612, With1612)

    pattern1613 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons66, cons491)
    rule1613 = ReplacementRule(pattern1613, With1613)
    pattern1614 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons66, cons491)
    rule1614 = ReplacementRule(pattern1614, replacement1614)
    pattern1615 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons5, cons860, cons543, cons25)
    rule1615 = ReplacementRule(pattern1615, replacement1615)
    pattern1616 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons860, cons543, cons25)
    rule1616 = ReplacementRule(pattern1616, replacement1616)
    pattern1617 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons36, cons37, cons19, cons4, cons5, cons55)
    rule1617 = ReplacementRule(pattern1617, replacement1617)
    pattern1618 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons920)
    rule1618 = ReplacementRule(pattern1618, replacement1618)
    pattern1619 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons920)
    rule1619 = ReplacementRule(pattern1619, replacement1619)
    pattern1620 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons19, cons4, cons5, cons556, cons921)
    rule1620 = ReplacementRule(pattern1620, replacement1620)
    pattern1621 = Pattern(Integral(Pq_*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons554, cons921)
    rule1621 = ReplacementRule(pattern1621, replacement1621)
    pattern1622 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons66, cons57, cons496)
    rule1622 = ReplacementRule(pattern1622, replacement1622)
    pattern1623 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons4, cons5, cons66, cons57, cons496)
    rule1623 = ReplacementRule(pattern1623, replacement1623)
    pattern1624 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons8, cons19, cons4, cons5, cons66, cons57)
    rule1624 = ReplacementRule(pattern1624, replacement1624)
    pattern1625 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons59, cons60, cons61, cons62, cons4, cons5, cons66, cons57)
    rule1625 = ReplacementRule(pattern1625, replacement1625)
    pattern1626 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons48, cons922, cons923)
    rule1626 = ReplacementRule(pattern1626, replacement1626)
    pattern1627 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons210, cons4, cons5, cons48, cons924, cons923)
    rule1627 = ReplacementRule(pattern1627, replacement1627)
    pattern1628 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons48, cons925, cons926, cons68)
    rule1628 = ReplacementRule(pattern1628, replacement1628)
    pattern1629 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons210, cons211, cons19, cons4, cons5, cons48, cons597, cons926, cons68)
    rule1629 = ReplacementRule(pattern1629, replacement1629)
    pattern1630 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons19, cons4, cons5, cons52, cons73, cons55)
    rule1630 = ReplacementRule(pattern1630, replacement1630)

    pattern1631 = Pattern(Integral(Px_**WC('q', S(1))*((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons927, cons588, cons89)
    rule1631 = ReplacementRule(pattern1631, With1631)
    pattern1632 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons48, cons860, cons55)
    rule1632 = ReplacementRule(pattern1632, replacement1632)
    pattern1633 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons48, cons66, cons130)
    rule1633 = ReplacementRule(pattern1633, replacement1633)
    pattern1634 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons48, cons66, cons130)
    rule1634 = ReplacementRule(pattern1634, replacement1634)
    pattern1635 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons928, cons929)
    rule1635 = ReplacementRule(pattern1635, replacement1635)
    pattern1636 = Pattern(Integral((d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons4, cons5, cons48, cons924, cons930)
    rule1636 = ReplacementRule(pattern1636, replacement1636)
    pattern1637 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons48, cons931, cons932, cons68)
    rule1637 = ReplacementRule(pattern1637, replacement1637)
    pattern1638 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons19, cons4, cons5, cons48, cons597, cons930, cons68)
    rule1638 = ReplacementRule(pattern1638, replacement1638)
    pattern1639 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons48, cons66, cons47, cons316)
    rule1639 = ReplacementRule(pattern1639, replacement1639)
    pattern1640 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons66, cons47, cons316)
    rule1640 = ReplacementRule(pattern1640, replacement1640)
    pattern1641 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons48, cons860, cons228, cons502)
    rule1641 = ReplacementRule(pattern1641, replacement1641)
    pattern1642 = Pattern(Integral(Pq_*(d_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons48, cons860, cons228, cons502)
    rule1642 = ReplacementRule(pattern1642, replacement1642)
    pattern1643 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons48, cons66, cons861)
    rule1643 = ReplacementRule(pattern1643, replacement1643)
    pattern1644 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons66, cons861, cons862)
    rule1644 = ReplacementRule(pattern1644, replacement1644)
    pattern1645 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons48, cons933, cons228, cons934, cons935)
    rule1645 = ReplacementRule(pattern1645, replacement1645)
    pattern1646 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons4, cons5, cons48, cons933, cons228, cons936, cons937)
    rule1646 = ReplacementRule(pattern1646, replacement1646)
    pattern1647 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons210, cons4, cons5, cons48, cons933, cons228, cons934, cons938)
    rule1647 = ReplacementRule(pattern1647, replacement1647)
    pattern1648 = Pattern(Integral((d_ + x_**WC('n3', S(1))*WC('g', S(1)))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons210, cons4, cons5, cons48, cons933, cons228, cons936, cons939)
    rule1648 = ReplacementRule(pattern1648, replacement1648)
    pattern1649 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons211, cons19, cons4, cons48, cons940, cons941, cons942, cons228, cons943, cons854)
    rule1649 = ReplacementRule(pattern1649, replacement1649)
    pattern1650 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons48, cons940, cons941, cons942, cons228, cons943, cons854)
    rule1650 = ReplacementRule(pattern1650, replacement1650)

    pattern1651 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons48, cons66, cons228, cons150, cons13, cons139, CustomConstraint(With1651))
    rule1651 = ReplacementRule(pattern1651, replacement1651)
    pattern1652 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons228, cons944)
    rule1652 = ReplacementRule(pattern1652, replacement1652)
    pattern1653 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons228, cons944)
    rule1653 = ReplacementRule(pattern1653, replacement1653)
    pattern1654 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons210, cons228, cons944)
    rule1654 = ReplacementRule(pattern1654, replacement1654)
    pattern1655 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons50, cons127, cons210, cons211, cons228, cons945, cons946)
    rule1655 = ReplacementRule(pattern1655, replacement1655)
    pattern1656 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons50, cons210, cons211, cons228, cons945, cons946)
    rule1656 = ReplacementRule(pattern1656, replacement1656)
    pattern1657 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons228, cons945, cons947)
    rule1657 = ReplacementRule(pattern1657, replacement1657)
    pattern1658 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons211, cons228, cons945, cons948)
    rule1658 = ReplacementRule(pattern1658, replacement1658)

    pattern1659 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons48, cons66, cons228, cons150, cons13, cons139, CustomConstraint(With1659))
    rule1659 = ReplacementRule(pattern1659, replacement1659)

    pattern1660 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons48, cons66, cons228, cons150, cons13, cons139, cons86, CustomConstraint(With1660))
    rule1660 = ReplacementRule(pattern1660, replacement1660)

    pattern1661 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons860, cons228, cons150, cons20, CustomConstraint(With1661))
    rule1661 = ReplacementRule(pattern1661, replacement1661)
    pattern1662 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons48, cons860, cons228, cons150, cons284)
    rule1662 = ReplacementRule(pattern1662, replacement1662)
    pattern1663 = Pattern(Integral(Pq_/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons8, cons48, cons860, cons228, cons150, cons949)
    rule1663 = ReplacementRule(pattern1663, replacement1663)

    pattern1664 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons66, cons228, cons65, CustomConstraint(With1664))
    rule1664 = ReplacementRule(pattern1664, replacement1664)

    pattern1665 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons66, cons228, cons721, cons950, CustomConstraint(With1665))
    rule1665 = ReplacementRule(pattern1665, replacement1665)

    pattern1666 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons66, cons228, cons721, cons951, CustomConstraint(With1666))
    rule1666 = ReplacementRule(pattern1666, replacement1666)

    pattern1667 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons48, cons860, cons228, cons150, CustomConstraint(With1667))
    rule1667 = ReplacementRule(pattern1667, replacement1667)

    pattern1668 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons860, cons228, cons150, CustomConstraint(With1668))
    rule1668 = ReplacementRule(pattern1668, replacement1668)

    pattern1669 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons48, cons66, cons228, cons150, cons952)
    rule1669 = ReplacementRule(pattern1669, With1669)

    pattern1670 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons66, cons228, cons150, cons952)
    rule1670 = ReplacementRule(pattern1670, With1670)
    pattern1671 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons48, cons66, cons228, cons150)
    rule1671 = ReplacementRule(pattern1671, replacement1671)
    pattern1672 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons48, cons66, cons228, cons150)
    rule1672 = ReplacementRule(pattern1672, replacement1672)

    pattern1673 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons66, cons228, cons198, cons20)
    rule1673 = ReplacementRule(pattern1673, With1673)

    pattern1674 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons48, cons66, cons228, cons198, cons369)
    rule1674 = ReplacementRule(pattern1674, With1674)

    pattern1675 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons48, cons66, cons228, cons198, cons358)
    rule1675 = ReplacementRule(pattern1675, With1675)

    pattern1676 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons48, cons66, cons228, cons491)
    rule1676 = ReplacementRule(pattern1676, With1676)

    pattern1677 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons66, cons228, cons491)
    rule1677 = ReplacementRule(pattern1677, With1677)
    pattern1678 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons48, cons66, cons228, cons491, cons75)
    rule1678 = ReplacementRule(pattern1678, replacement1678)
    pattern1679 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons48, cons66, cons228, cons491, cons953)
    rule1679 = ReplacementRule(pattern1679, replacement1679)
    pattern1680 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons48, cons66, cons228, cons491)
    rule1680 = ReplacementRule(pattern1680, replacement1680)
    pattern1681 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons48, cons860, cons228, cons543, cons25)
    rule1681 = ReplacementRule(pattern1681, replacement1681)
    pattern1682 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons48, cons860, cons228, cons543, cons25)
    rule1682 = ReplacementRule(pattern1682, replacement1682)

    pattern1683 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons48, cons66, cons228)
    rule1683 = ReplacementRule(pattern1683, With1683)

    pattern1684 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons4, cons48, cons66, cons228)
    rule1684 = ReplacementRule(pattern1684, With1684)
    pattern1685 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons48, cons66, cons704)
    rule1685 = ReplacementRule(pattern1685, replacement1685)
    pattern1686 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons66, cons704)
    rule1686 = ReplacementRule(pattern1686, replacement1686)
    pattern1687 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons48, cons920)
    rule1687 = ReplacementRule(pattern1687, replacement1687)
    pattern1688 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons48, cons920)
    rule1688 = ReplacementRule(pattern1688, replacement1688)
    pattern1689 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons48, cons556, cons921)
    rule1689 = ReplacementRule(pattern1689, replacement1689)
    pattern1690 = Pattern(Integral(Pq_*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons48, cons554, cons921)
    rule1690 = ReplacementRule(pattern1690, replacement1690)
    pattern1691 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons149, cons954, cons955)
    rule1691 = ReplacementRule(pattern1691, replacement1691)
    pattern1692 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons149, cons954, cons956, cons13, cons139)
    rule1692 = ReplacementRule(pattern1692, replacement1692)
    pattern1693 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons149, cons954, cons956, cons957)
    rule1693 = ReplacementRule(pattern1693, replacement1693)
    pattern1694 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons149, cons958, cons959, cons165, cons960)
    rule1694 = ReplacementRule(pattern1694, replacement1694)
    pattern1695 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons149, cons958, cons959, cons165, cons961)
    rule1695 = ReplacementRule(pattern1695, replacement1695)
    pattern1696 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons149, cons958, cons959, cons139, cons962)
    rule1696 = ReplacementRule(pattern1696, replacement1696)
    pattern1697 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons149, cons958, cons959, cons139)
    rule1697 = ReplacementRule(pattern1697, replacement1697)
    pattern1698 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons963, cons954, cons964)
    rule1698 = ReplacementRule(pattern1698, replacement1698)
    pattern1699 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons4, cons965)
    rule1699 = ReplacementRule(pattern1699, replacement1699)
    pattern1700 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons721, cons954, cons964)
    rule1700 = ReplacementRule(pattern1700, replacement1700)
    pattern1701 = Pattern(Integral(S(1)/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons966, cons967)
    rule1701 = ReplacementRule(pattern1701, replacement1701)
    pattern1702 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons149, cons954, cons968)
    rule1702 = ReplacementRule(pattern1702, replacement1702)
    pattern1703 = Pattern(Integral((u_**WC('j', S(1))*WC('a', S(1)) + u_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons70, cons71)
    rule1703 = ReplacementRule(pattern1703, replacement1703)
    pattern1704 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons19, cons4, cons5, cons149, cons954, cons969, cons55)
    rule1704 = ReplacementRule(pattern1704, replacement1704)
    pattern1705 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons970, cons971)
    rule1705 = ReplacementRule(pattern1705, replacement1705)
    pattern1706 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons149, cons954, cons972, cons13, cons139, cons971)
    rule1706 = ReplacementRule(pattern1706, replacement1706)
    pattern1707 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons972, cons973, cons974)
    rule1707 = ReplacementRule(pattern1707, replacement1707)
    pattern1708 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons972)
    rule1708 = ReplacementRule(pattern1708, replacement1708)
    pattern1709 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons19, cons4, cons5, cons149, cons954, cons969, cons502, cons975)
    rule1709 = ReplacementRule(pattern1709, replacement1709)
    pattern1710 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons969, cons502, cons975)
    rule1710 = ReplacementRule(pattern1710, replacement1710)
    pattern1711 = Pattern(Integral((x_*WC('c', S(1)))**m_*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons149, cons976, cons959, cons974, cons165, cons977)
    rule1711 = ReplacementRule(pattern1711, replacement1711)
    pattern1712 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons149, cons958, cons959, cons974, cons165, cons514)
    rule1712 = ReplacementRule(pattern1712, replacement1712)
    pattern1713 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons149, cons976, cons959, cons974, cons139, cons978)
    rule1713 = ReplacementRule(pattern1713, replacement1713)
    pattern1714 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons149, cons958, cons959, cons974, cons139)
    rule1714 = ReplacementRule(pattern1714, replacement1714)
    pattern1715 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons149, cons966, cons959, cons974, cons979, cons514)
    rule1715 = ReplacementRule(pattern1715, replacement1715)
    pattern1716 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons149, cons966, cons959, cons974, cons980)
    rule1716 = ReplacementRule(pattern1716, replacement1716)
    pattern1717 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons19, cons4, cons5, cons149, cons954, cons969, cons68, cons543, cons25)
    rule1717 = ReplacementRule(pattern1717, replacement1717)
    pattern1718 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons969, cons68, cons543, cons25)
    rule1718 = ReplacementRule(pattern1718, replacement1718)
    pattern1719 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons963, cons954, cons981, cons971)
    rule1719 = ReplacementRule(pattern1719, replacement1719)
    pattern1720 = Pattern(Integral(x_**WC('m', S(1))/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons798, cons4, cons982, cons954)
    rule1720 = ReplacementRule(pattern1720, replacement1720)
    pattern1721 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons721, cons954, cons981, cons971)
    rule1721 = ReplacementRule(pattern1721, replacement1721)
    pattern1722 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons669, cons954, cons981)
    rule1722 = ReplacementRule(pattern1722, replacement1722)
    pattern1723 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons149, cons954, cons968)
    rule1723 = ReplacementRule(pattern1723, replacement1723)
    pattern1724 = Pattern(Integral(u_**WC('m', S(1))*(v_**WC('j', S(1))*WC('a', S(1)) + v_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons798, cons19, cons4, cons5, cons556)
    rule1724 = ReplacementRule(pattern1724, replacement1724)
    pattern1725 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons798, cons799, cons19, cons4, cons5, cons52, cons149, cons983, cons969, cons984, cons502, cons975)
    rule1725 = ReplacementRule(pattern1725, replacement1725)
    pattern1726 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons799, cons19, cons4, cons5, cons52, cons149, cons983, cons969, cons984, cons502, cons975)
    rule1726 = ReplacementRule(pattern1726, replacement1726)
    pattern1727 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons19, cons4, cons5, cons985, cons149, cons73, cons986, cons987, cons973)
    rule1727 = ReplacementRule(pattern1727, replacement1727)
    pattern1728 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons19, cons4, cons985, cons149, cons73, cons988, cons139, cons989, cons990)
    rule1728 = ReplacementRule(pattern1728, replacement1728)
    pattern1729 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons5, cons985, cons149, cons73, cons95, cons90, cons991, cons992, cons973, cons993)
    rule1729 = ReplacementRule(pattern1729, replacement1729)
    pattern1730 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons19, cons4, cons5, cons985, cons149, cons73, cons994, cons990)
    rule1730 = ReplacementRule(pattern1730, replacement1730)
    pattern1731 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons798, cons799, cons19, cons4, cons5, cons52, cons149, cons983, cons969, cons984, cons68, cons543, cons25)
    rule1731 = ReplacementRule(pattern1731, replacement1731)
    pattern1732 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons799, cons19, cons4, cons5, cons52, cons149, cons983, cons969, cons984, cons68, cons543, cons25)
    rule1732 = ReplacementRule(pattern1732, replacement1732)
    pattern1733 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons798, cons19, cons4, cons5, cons52, cons985, cons149, cons73, cons995)
    rule1733 = ReplacementRule(pattern1733, replacement1733)

    pattern1734 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons860, cons149, cons954, cons966, cons969, cons996)
    rule1734 = ReplacementRule(pattern1734, With1734)
    pattern1735 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons19, cons4, cons5, cons860, cons149, cons954, cons969, cons502)
    rule1735 = ReplacementRule(pattern1735, replacement1735)
    pattern1736 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons4, cons5, cons860, cons149, cons954, cons969, cons502, cons33, cons997)
    rule1736 = ReplacementRule(pattern1736, replacement1736)
    pattern1737 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons860, cons149, cons954, cons969, cons502)
    rule1737 = ReplacementRule(pattern1737, replacement1737)

    pattern1738 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons5, cons860, cons149, cons998, cons20, CustomConstraint(With1738))
    rule1738 = ReplacementRule(pattern1738, replacement1738)

    pattern1739 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons66, cons149, cons999, cons1000, CustomConstraint(With1739))
    rule1739 = ReplacementRule(pattern1739, replacement1739)
    pattern1740 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons19, cons4, cons5, cons860, cons149, cons954, cons969, cons543, cons25)
    rule1740 = ReplacementRule(pattern1740, replacement1740)
    pattern1741 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons4, cons5, cons860, cons149, cons954, cons969, cons543, cons25, cons33, cons997)
    rule1741 = ReplacementRule(pattern1741, replacement1741)
    pattern1742 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons860, cons149, cons954, cons969, cons543, cons25)
    rule1742 = ReplacementRule(pattern1742, replacement1742)
    pattern1743 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons798, cons19, cons4, cons5, cons920, cons149, cons954)
    rule1743 = ReplacementRule(pattern1743, replacement1743)
    pattern1744 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons798, cons4, cons5, cons920, cons149, cons954)
    rule1744 = ReplacementRule(pattern1744, replacement1744)
    pattern1745 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons40, cons1001)
    rule1745 = ReplacementRule(pattern1745, replacement1745)
    pattern1746 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons130, cons1002)
    rule1746 = ReplacementRule(pattern1746, replacement1746)

    pattern1747 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons65, cons1002, CustomConstraint(With1747))
    rule1747 = ReplacementRule(pattern1747, replacement1747)

    pattern1748 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons65, cons1002)
    rule1748 = ReplacementRule(pattern1748, With1748)
    pattern1749 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons5, cons149, cons1001)
    rule1749 = ReplacementRule(pattern1749, replacement1749)

    pattern1750 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons5, cons149, cons1002, CustomConstraint(With1750))
    rule1750 = ReplacementRule(pattern1750, replacement1750)

    pattern1751 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons5, cons149, cons1002)
    rule1751 = ReplacementRule(pattern1751, With1751)
    pattern1752 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons40, cons1001)
    rule1752 = ReplacementRule(pattern1752, replacement1752)
    pattern1753 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons130, cons1002)
    rule1753 = ReplacementRule(pattern1753, replacement1753)

    pattern1754 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons65, cons1002, CustomConstraint(With1754))
    rule1754 = ReplacementRule(pattern1754, replacement1754)

    pattern1755 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons65, cons1002)
    rule1755 = ReplacementRule(pattern1755, With1755)
    pattern1756 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons149, cons1001)
    rule1756 = ReplacementRule(pattern1756, replacement1756)

    pattern1757 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons149, cons1002, CustomConstraint(With1757))
    rule1757 = ReplacementRule(pattern1757, replacement1757)

    pattern1758 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons149, cons1002)
    rule1758 = ReplacementRule(pattern1758, With1758)
    pattern1759 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons40, cons1003)
    rule1759 = ReplacementRule(pattern1759, replacement1759)
    pattern1760 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons130, cons1004)
    rule1760 = ReplacementRule(pattern1760, replacement1760)

    pattern1761 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons65, cons1004, CustomConstraint(With1761))
    rule1761 = ReplacementRule(pattern1761, replacement1761)

    pattern1762 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons65, cons1004)
    rule1762 = ReplacementRule(pattern1762, With1762)
    pattern1763 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons5, cons149, cons1003)
    rule1763 = ReplacementRule(pattern1763, replacement1763)

    pattern1764 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons5, cons149, cons1004, CustomConstraint(With1764))
    rule1764 = ReplacementRule(pattern1764, replacement1764)

    pattern1765 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons5, cons149, cons1004)
    rule1765 = ReplacementRule(pattern1765, With1765)
    pattern1766 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons40, cons1003)
    rule1766 = ReplacementRule(pattern1766, replacement1766)
    pattern1767 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons130, cons1004)
    rule1767 = ReplacementRule(pattern1767, replacement1767)

    pattern1768 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons65, cons1004, CustomConstraint(With1768))
    rule1768 = ReplacementRule(pattern1768, replacement1768)

    pattern1769 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons65, cons1004)
    rule1769 = ReplacementRule(pattern1769, With1769)
    pattern1770 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1003)
    rule1770 = ReplacementRule(pattern1770, replacement1770)

    pattern1771 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1004, CustomConstraint(With1771))
    rule1771 = ReplacementRule(pattern1771, replacement1771)

    pattern1772 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1004)
    rule1772 = ReplacementRule(pattern1772, With1772)
    pattern1773 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons40, cons1005, cons1006)
    rule1773 = ReplacementRule(pattern1773, replacement1773)
    pattern1774 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons40, cons1005, cons1007)
    rule1774 = ReplacementRule(pattern1774, replacement1774)

    pattern1775 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons40, cons1008, cons1006)
    rule1775 = ReplacementRule(pattern1775, With1775)
    pattern1776 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons130, cons1008, cons1007)
    rule1776 = ReplacementRule(pattern1776, replacement1776)

    pattern1777 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons65, cons1008, cons1007, CustomConstraint(With1777))
    rule1777 = ReplacementRule(pattern1777, replacement1777)
    pattern1778 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons65, cons1008, cons1007)
    rule1778 = ReplacementRule(pattern1778, replacement1778)
    pattern1779 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons149, cons1005, cons1006)
    rule1779 = ReplacementRule(pattern1779, replacement1779)

    pattern1780 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons149, cons1005, cons1007)
    rule1780 = ReplacementRule(pattern1780, With1780)

    pattern1781 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons149, cons1008, cons1006)
    rule1781 = ReplacementRule(pattern1781, With1781)

    pattern1782 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons149, cons1008, cons1007, CustomConstraint(With1782))
    rule1782 = ReplacementRule(pattern1782, replacement1782)

    pattern1783 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons149, cons1008, cons1007)
    rule1783 = ReplacementRule(pattern1783, With1783)
    pattern1784 = Pattern(Integral(u_**p_, x_), cons5, cons1009, cons1010)
    rule1784 = ReplacementRule(pattern1784, replacement1784)
    pattern1785 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons40, cons1005, cons1006)
    rule1785 = ReplacementRule(pattern1785, replacement1785)

    pattern1786 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons40, cons1005, cons1007)
    rule1786 = ReplacementRule(pattern1786, With1786)

    pattern1787 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons40, cons1008, cons1006)
    rule1787 = ReplacementRule(pattern1787, With1787)
    pattern1788 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons130, cons1008, cons1007)
    rule1788 = ReplacementRule(pattern1788, replacement1788)

    pattern1789 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons65, cons1008, cons1007, CustomConstraint(With1789))
    rule1789 = ReplacementRule(pattern1789, replacement1789)
    pattern1790 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons65, cons1008, cons1007)
    rule1790 = ReplacementRule(pattern1790, replacement1790)
    pattern1791 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1005, cons1006)
    rule1791 = ReplacementRule(pattern1791, replacement1791)

    pattern1792 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1005, cons1007)
    rule1792 = ReplacementRule(pattern1792, With1792)

    pattern1793 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1008, cons1006)
    rule1793 = ReplacementRule(pattern1793, With1793)

    pattern1794 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1008, cons1007, CustomConstraint(With1794))
    rule1794 = ReplacementRule(pattern1794, replacement1794)

    pattern1795 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons149, cons1008, cons1007)
    rule1795 = ReplacementRule(pattern1795, With1795)
    pattern1796 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), cons19, cons5, cons70, cons1011, cons1012)
    rule1796 = ReplacementRule(pattern1796, replacement1796)
    pattern1797 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons385, cons1013, cons1014)
    rule1797 = ReplacementRule(pattern1797, replacement1797)
    pattern1798 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons385, cons1013, cons1015)
    rule1798 = ReplacementRule(pattern1798, replacement1798)
    pattern1799 = Pattern(Integral((x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1016, cons1017, cons1018)
    rule1799 = ReplacementRule(pattern1799, replacement1799)

    pattern1800 = Pattern(Integral(v_**p_, x_), cons5, cons1019, cons1020, cons1017, cons1018, CustomConstraint(With1800))
    rule1800 = ReplacementRule(pattern1800, replacement1800)
    pattern1801 = Pattern(Integral(u_*(x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons806, cons1016, cons359)
    rule1801 = ReplacementRule(pattern1801, replacement1801)

    pattern1802 = Pattern(Integral(u_*v_**p_, x_), cons5, cons806, cons1019, cons1020, cons359, CustomConstraint(With1802))
    rule1802 = ReplacementRule(pattern1802, replacement1802)
    pattern1803 = Pattern(Integral((a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons1021, cons248)
    rule1803 = ReplacementRule(pattern1803, replacement1803)

    pattern1804 = Pattern(Integral(v_**p_, x_), cons5, cons1019, cons1020, cons248, CustomConstraint(With1804))
    rule1804 = ReplacementRule(pattern1804, replacement1804)

    pattern1805 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons1025, cons1022, cons1023, cons1024)
    rule1805 = ReplacementRule(pattern1805, With1805)

    pattern1806 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons1025, cons1022, cons1023, cons1024)
    rule1806 = ReplacementRule(pattern1806, With1806)

    pattern1807 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons1025, cons19, cons1022, cons1023, cons1024)
    rule1807 = ReplacementRule(pattern1807, With1807)

    pattern1808 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons1025, cons19, cons1022, cons1023, cons1024)
    rule1808 = ReplacementRule(pattern1808, With1808)

    pattern1809 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1026, cons1027, cons1028)
    rule1809 = ReplacementRule(pattern1809, With1809)

    pattern1810 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1029, cons1030, cons1031)
    rule1810 = ReplacementRule(pattern1810, With1810)

    pattern1811 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1026, cons1027, cons1032)
    rule1811 = ReplacementRule(pattern1811, With1811)

    pattern1812 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1029, cons1030, cons1033)
    rule1812 = ReplacementRule(pattern1812, With1812)
    pattern1813 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1025, cons1034, cons1035)
    rule1813 = ReplacementRule(pattern1813, replacement1813)
    pattern1814 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1025, cons1036, cons1037)
    rule1814 = ReplacementRule(pattern1814, replacement1814)
    pattern1815 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1038)
    rule1815 = ReplacementRule(pattern1815, replacement1815)
    pattern1816 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1039)
    rule1816 = ReplacementRule(pattern1816, replacement1816)
    pattern1817 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1040, cons1041)
    rule1817 = ReplacementRule(pattern1817, replacement1817)
    pattern1818 = Pattern(Integral(WC('u', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons1042, cons1043)
    rule1818 = ReplacementRule(pattern1818, replacement1818)
    pattern1819 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1042, cons1044)
    rule1819 = ReplacementRule(pattern1819, replacement1819)

    pattern1820 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons470)
    rule1820 = ReplacementRule(pattern1820, With1820)

    pattern1821 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons29, cons127, cons470)
    rule1821 = ReplacementRule(pattern1821, With1821)

    pattern1822 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons471)
    rule1822 = ReplacementRule(pattern1822, With1822)

    pattern1823 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons29, cons127, cons471)
    rule1823 = ReplacementRule(pattern1823, With1823)
    pattern1824 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule1824 = ReplacementRule(pattern1824, replacement1824)
    pattern1825 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons299)
    rule1825 = ReplacementRule(pattern1825, replacement1825)
    pattern1826 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1046, cons1047)
    rule1826 = ReplacementRule(pattern1826, replacement1826)
    pattern1827 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1046, cons1048)
    rule1827 = ReplacementRule(pattern1827, replacement1827)
    pattern1828 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))**S(2)), x_), cons2, cons8, cons29, cons50, cons1049)
    rule1828 = ReplacementRule(pattern1828, replacement1828)
    pattern1829 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1050, cons707)
    rule1829 = ReplacementRule(pattern1829, replacement1829)
    pattern1830 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons36, cons37, cons1050, cons707)
    rule1830 = ReplacementRule(pattern1830, replacement1830)
    pattern1831 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons384, cons1051)
    rule1831 = ReplacementRule(pattern1831, replacement1831)
    pattern1832 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons127, cons36, cons37, cons384, cons1051)
    rule1832 = ReplacementRule(pattern1832, replacement1832)
    pattern1833 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/((d_ + x_**S(4)*WC('f', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons127, cons36, cons37, cons384, cons1051)
    rule1833 = ReplacementRule(pattern1833, replacement1833)
    pattern1834 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1052, cons699)
    rule1834 = ReplacementRule(pattern1834, replacement1834)

    pattern1835 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1052, cons711)
    rule1835 = ReplacementRule(pattern1835, With1835)
    pattern1836 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule1836 = ReplacementRule(pattern1836, replacement1836)
    pattern1837 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*sqrt(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons1053, cons1054)
    rule1837 = ReplacementRule(pattern1837, replacement1837)
    pattern1838 = Pattern(Integral((u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons127, cons210, cons211, cons798, cons799, cons19, cons4, cons70, cons820, cons1055, cons1056)
    rule1838 = ReplacementRule(pattern1838, replacement1838)
    pattern1839 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1057, cons40)
    rule1839 = ReplacementRule(pattern1839, replacement1839)
    pattern1840 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1057, cons40)
    rule1840 = ReplacementRule(pattern1840, replacement1840)
    pattern1841 = Pattern(Integral(((u_ + sqrt(v_)*WC('f', S(1)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons127, cons210, cons211, cons4, cons70, cons820, cons821, cons1058, cons40)
    rule1841 = ReplacementRule(pattern1841, replacement1841)
    pattern1842 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons8, cons50, cons127, cons210, cons211, cons4, cons1057, cons20)
    rule1842 = ReplacementRule(pattern1842, replacement1842)
    pattern1843 = Pattern(Integral(x_**WC('p', S(1))*(g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)))**WC('n', S(1)), x_), cons2, cons8, cons50, cons127, cons210, cons226, cons4, cons1057, cons1059, cons1060, cons1061)
    rule1843 = ReplacementRule(pattern1843, replacement1843)
    pattern1844 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons4, cons1057, cons1059, cons1062, cons517, cons1061)
    rule1844 = ReplacementRule(pattern1844, replacement1844)
    pattern1845 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons226, cons4, cons1057, cons1059, cons517, cons1061)
    rule1845 = ReplacementRule(pattern1845, replacement1845)
    pattern1846 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons4, cons1057, cons1059, cons1062, cons75, cons1063)
    rule1846 = ReplacementRule(pattern1846, replacement1846)
    pattern1847 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons226, cons4, cons1057, cons1059, cons75, cons1063)
    rule1847 = ReplacementRule(pattern1847, replacement1847)
    pattern1848 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons4, cons1057, cons1059, cons1062, cons953, cons1063)
    rule1848 = ReplacementRule(pattern1848, replacement1848)
    pattern1849 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons226, cons4, cons1057, cons1059, cons953, cons1063)
    rule1849 = ReplacementRule(pattern1849, replacement1849)
    pattern1850 = Pattern(Integral(w_**WC('m', S(1))*(u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1)), x_), cons127, cons798, cons799, cons19, cons4, cons70, cons1064, cons1065, cons1066)
    rule1850 = ReplacementRule(pattern1850, replacement1850)
    pattern1851 = Pattern(Integral(S(1)/((a_ + x_**WC('n', S(1))*WC('b', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + (a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons4, cons1067)
    rule1851 = ReplacementRule(pattern1851, replacement1851)
    pattern1852 = Pattern(Integral(sqrt(a_ + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons1068)
    rule1852 = ReplacementRule(pattern1852, replacement1852)
    pattern1853 = Pattern(Integral(sqrt(x_**S(2)*WC('a', S(1)) + x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1069, cons1070)
    rule1853 = ReplacementRule(pattern1853, replacement1853)
    pattern1854 = Pattern(Integral(sqrt(x_*(x_*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))*WC('e', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1069, cons1071)
    rule1854 = ReplacementRule(pattern1854, replacement1854)
    pattern1855 = Pattern(Integral(sqrt(x_**S(2)*WC('c', S(1)) + sqrt(a_ + x_**S(4)*WC('b', S(1)))*WC('d', S(1)))/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons1072)
    rule1855 = ReplacementRule(pattern1855, replacement1855)
    pattern1856 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sqrt(x_**S(2)*WC('b', S(1)) + sqrt(a_ + x_**S(4)*WC('e', S(1))))/sqrt(a_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons1073, cons45)
    rule1856 = ReplacementRule(pattern1856, replacement1856)

    pattern1857 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons483, cons481)
    rule1857 = ReplacementRule(pattern1857, With1857)

    pattern1858 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons483, cons482)
    rule1858 = ReplacementRule(pattern1858, With1858)

    pattern1859 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons484, cons481)
    rule1859 = ReplacementRule(pattern1859, With1859)

    pattern1860 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons484, cons482)
    rule1860 = ReplacementRule(pattern1860, With1860)

    pattern1861 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons483, cons481, CustomConstraint(With1861))
    rule1861 = ReplacementRule(pattern1861, replacement1861)

    pattern1862 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons483, cons481, CustomConstraint(With1862))
    rule1862 = ReplacementRule(pattern1862, replacement1862)

    pattern1863 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons483, cons482, CustomConstraint(With1863))
    rule1863 = ReplacementRule(pattern1863, replacement1863)

    pattern1864 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons483, cons482, CustomConstraint(With1864))
    rule1864 = ReplacementRule(pattern1864, replacement1864)

    pattern1865 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons484, cons481, CustomConstraint(With1865))
    rule1865 = ReplacementRule(pattern1865, replacement1865)

    pattern1866 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons484, cons481, CustomConstraint(With1866))
    rule1866 = ReplacementRule(pattern1866, replacement1866)

    pattern1867 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons484, cons482, CustomConstraint(With1867))
    rule1867 = ReplacementRule(pattern1867, replacement1867)

    pattern1868 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons484, cons482, CustomConstraint(With1868))
    rule1868 = ReplacementRule(pattern1868, replacement1868)
    pattern1869 = Pattern(Integral(x_**WC('m', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1074, cons502)
    rule1869 = ReplacementRule(pattern1869, replacement1869)
    pattern1870 = Pattern(Integral(WC('u', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1074)
    rule1870 = ReplacementRule(pattern1870, replacement1870)
    pattern1871 = Pattern(Integral((A_ + x_**n_*WC('B', S(1)))/(a_ + x_**S(2)*WC('b', S(1)) + x_**n2_*WC('d', S(1)) + x_**n_*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons4, cons48, cons965, cons1075, cons1076)
    rule1871 = ReplacementRule(pattern1871, replacement1871)
    pattern1872 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('n', S(1))*WC('B', S(1)))/(a_ + x_**n2_*WC('d', S(1)) + x_**WC('k', S(1))*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons19, cons4, cons48, cons1077, cons1078, cons1079)
    rule1872 = ReplacementRule(pattern1872, replacement1872)
    pattern1873 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons48, cons933, cons228, cons704)
    rule1873 = ReplacementRule(pattern1873, replacement1873)
    pattern1874 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons48, cons228, cons704)
    rule1874 = ReplacementRule(pattern1874, replacement1874)
    pattern1875 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons210, cons4, cons48, cons933, cons228, cons704)
    rule1875 = ReplacementRule(pattern1875, replacement1875)
    pattern1876 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1))), x_), cons2, cons3, cons8, cons29, cons127, cons210, cons4, cons48, cons933, cons228, cons704)
    rule1876 = ReplacementRule(pattern1876, replacement1876)
    pattern1877 = Pattern(Integral((d_ + x_**n2_*WC('f', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons4, cons48, cons228, cons704)
    rule1877 = ReplacementRule(pattern1877, replacement1877)
    pattern1878 = Pattern(Integral((d_ + g_*x_**n3_)*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons210, cons4, cons48, cons933, cons228, cons704)
    rule1878 = ReplacementRule(pattern1878, replacement1878)
    pattern1879 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons48, cons933, cons704)
    rule1879 = ReplacementRule(pattern1879, replacement1879)
    pattern1880 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons4, cons48, cons704)
    rule1880 = ReplacementRule(pattern1880, replacement1880)
    pattern1881 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons210, cons4, cons48, cons933, cons704)
    rule1881 = ReplacementRule(pattern1881, replacement1881)

    pattern1882 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1080, cons1081, cons1082, cons1083, cons1084, cons1085)
    rule1882 = ReplacementRule(pattern1882, With1882)

    pattern1883 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons1086, cons1087, cons1082, cons1088)
    rule1883 = ReplacementRule(pattern1883, With1883)

    pattern1884 = Pattern(Integral(u_*v_**p_, x_), cons13, cons139, cons806, cons1019, cons1089, cons1090, cons1091, CustomConstraint(With1884))
    rule1884 = ReplacementRule(pattern1884, replacement1884)
    return [rule1476, rule1477, rule1478, rule1479, rule1480, rule1481, rule1482, rule1483, rule1484, rule1485, rule1486, rule1487, rule1488, rule1489, rule1490, rule1491, rule1492, rule1493, rule1494, rule1495, rule1496, rule1497, rule1498, rule1499, rule1500, rule1501, rule1502, rule1503, rule1504, rule1505, rule1506, rule1507, rule1508, rule1509, rule1510, rule1511, rule1512, rule1513, rule1514, rule1515, rule1516, rule1517, rule1518, rule1519, rule1520, rule1521, rule1522, rule1523, rule1524, rule1525, rule1526, rule1527, rule1528, rule1529, rule1530, rule1531, rule1532, rule1533, rule1534, rule1535, rule1536, rule1537, rule1538, rule1539, rule1540, rule1541, rule1542, rule1543, rule1544, rule1545, rule1546, rule1547, rule1548, rule1549, rule1550, rule1551, rule1552, rule1553, rule1554, rule1555, rule1556, rule1557, rule1558, rule1559, rule1560, rule1561, rule1562, rule1563, rule1564, rule1565, rule1566, rule1567, rule1568, rule1569, rule1570, rule1571, rule1572, rule1573, rule1574, rule1575, rule1576, rule1577, rule1578, rule1579, rule1580, rule1581, rule1582, rule1583, rule1584, rule1585, rule1586, rule1587, rule1588, rule1589, rule1590, rule1591, rule1592, rule1593, rule1594, rule1595, rule1596, rule1597, rule1598, rule1599, rule1600, rule1601, rule1602, rule1603, rule1604, rule1605, rule1606, rule1607, rule1608, rule1609, rule1610, rule1611, rule1612, rule1613, rule1614, rule1615, rule1616, rule1617, rule1618, rule1619, rule1620, rule1621, rule1622, rule1623, rule1624, rule1625, rule1626, rule1627, rule1628, rule1629, rule1630, rule1631, rule1632, rule1633, rule1634, rule1635, rule1636, rule1637, rule1638, rule1639, rule1640, rule1641, rule1642, rule1643, rule1644, rule1645, rule1646, rule1647, rule1648, rule1649, rule1650, rule1651, rule1652, rule1653, rule1654, rule1655, rule1656, rule1657, rule1658, rule1659, rule1660, rule1661, rule1662, rule1663, rule1664, rule1665, rule1666, rule1667, rule1668, rule1669, rule1670, rule1671, rule1672, rule1673, rule1674, rule1675, rule1676, rule1677, rule1678, rule1679, rule1680, rule1681, rule1682, rule1683, rule1684, rule1685, rule1686, rule1687, rule1688, rule1689, rule1690, rule1691, rule1692, rule1693, rule1694, rule1695, rule1696, rule1697, rule1698, rule1699, rule1700, rule1701, rule1702, rule1703, rule1704, rule1705, rule1706, rule1707, rule1708, rule1709, rule1710, rule1711, rule1712, rule1713, rule1714, rule1715, rule1716, rule1717, rule1718, rule1719, rule1720, rule1721, rule1722, rule1723, rule1724, rule1725, rule1726, rule1727, rule1728, rule1729, rule1730, rule1731, rule1732, rule1733, rule1734, rule1735, rule1736, rule1737, rule1738, rule1739, rule1740, rule1741, rule1742, rule1743, rule1744, rule1745, rule1746, rule1747, rule1748, rule1749, rule1750, rule1751, rule1752, rule1753, rule1754, rule1755, rule1756, rule1757, rule1758, rule1759, rule1760, rule1761, rule1762, rule1763, rule1764, rule1765, rule1766, rule1767, rule1768, rule1769, rule1770, rule1771, rule1772, rule1773, rule1774, rule1775, rule1776, rule1777, rule1778, rule1779, rule1780, rule1781, rule1782, rule1783, rule1784, rule1785, rule1786, rule1787, rule1788, rule1789, rule1790, rule1791, rule1792, rule1793, rule1794, rule1795, rule1796, rule1797, rule1798, rule1799, rule1800, rule1801, rule1802, rule1803, rule1804, rule1805, rule1806, rule1807, rule1808, rule1809, rule1810, rule1811, rule1812, rule1813, rule1814, rule1815, rule1816, rule1817, rule1818, rule1819, rule1820, rule1821, rule1822, rule1823, rule1824, rule1825, rule1826, rule1827, rule1828, rule1829, rule1830, rule1831, rule1832, rule1833, rule1834, rule1835, rule1836, rule1837, rule1838, rule1839, rule1840, rule1841, rule1842, rule1843, rule1844, rule1845, rule1846, rule1847, rule1848, rule1849, rule1850, rule1851, rule1852, rule1853, rule1854, rule1855, rule1856, rule1857, rule1858, rule1859, rule1860, rule1861, rule1862, rule1863, rule1864, rule1865, rule1866, rule1867, rule1868, rule1869, rule1870, rule1871, rule1872, rule1873, rule1874, rule1875, rule1876, rule1877, rule1878, rule1879, rule1880, rule1881, rule1882, rule1883, rule1884, ]



def replacement1476(c, q, n, x, a, p, b):
        # rubi.append(1476)
        return Dist(x*(c*x**n)**(-S(1)/n), Subst(Int((a + b*x**(n*q))**p, x), x, (c*x**n)**(S(1)/n)), x)
def replacement1477(c, m, q, n, x, a, p, b):
        # rubi.append(1477)
        return Dist(x**(m + S(1))*(c*x**n)**(-(m + S(1))/n), Subst(Int(x**m*(a + b*x**(n*q))**p, x), x, (c*x**n)**(S(1)/n)), x)
def replacement1478(c, m, q, n, r, x, d, f, a, p, e, s, b):
        # rubi.append(1478)
        return Dist((e*(a + b*x**n)**r)**p*(f*(c + d*x**n)**s)**q*(a + b*x**n)**(-p*r)*(c + d*x**n)**(-q*s), Int(x**m*(a + b*x**n)**(p*r)*(c + d*x**n)**(q*s), x), x)
def replacement1479(c, n, x, d, u, a, p, e, b):
        # rubi.append(1479)
        return Dist((b*e/d)**p, Int(u, x), x)
def replacement1480(c, n, x, d, u, a, p, e, b):
        # rubi.append(1480)
        return Int(u*(e*(a + b*x**n))**p*(c + d*x**n)**(-p), x)

def With1481(c, n, x, d, a, p, e, b):
        q = Denominator(p)
        # rubi.append(1481)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + S(1)/n)*(b*e - d*x**q)**(S(-1) - S(1)/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)

def With1482(c, m, n, x, d, a, p, e, b):
        q = Denominator(p)
        # rubi.append(1482)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)

def With1483(c, n, r, x, d, u, a, p, e, b):
        q = Denominator(p)
        # rubi.append(1483)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + S(1)/n)*(b*e - d*x**q)**(S(-1) - S(1)/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(S(1)/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)

def With1484(c, m, n, r, x, d, u, a, p, e, b):
        q = Denominator(p)
        # rubi.append(1484)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(S(1)/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)
def replacement1485(c, n, x, a, p, b):
        # rubi.append(1485)
        return -Dist(c, Subst(Int((a + b*x**n)**p/x**S(2), x), x, c/x), x)
def replacement1486(c, m, n, x, a, p, b):
        # rubi.append(1486)
        return -Dist(c**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x), x)
def replacement1487(c, m, n, x, d, a, p, b):
        # rubi.append(1487)
        return -Dist(c*(c/x)**m*(d*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x), x)
def replacement1488(c, n, x, d, n2, a, p, b):
        # rubi.append(1488)
        return -Dist(d, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p/x**S(2), x), x, d/x), x)
def replacement1489(c, m, n, x, d, n2, a, p, b):
        # rubi.append(1489)
        return -Dist(d**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x), x)
def replacement1490(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(1490)
        return -Dist(d*(d/x)**m*(e*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x), x)
def replacement1491(c, n, x, d, n2, a, p, b):
        # rubi.append(1491)
        return -Dist(d, Subst(Int((a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p/x**S(2), x), x, d/x), x)
def replacement1492(c, m, n, x, d, n2, a, p, b):
        # rubi.append(1492)
        return -Dist(d**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x), x)
def replacement1493(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(1493)
        return -Dist(d*(d/x)**m*(e*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x), x)
def replacement1494(x, u, m):
        # rubi.append(1494)
        return Int(ExpandToSum(u, x)**m, x)
def replacement1495(m, n, x, v, u):
        # rubi.append(1495)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n, x)
def replacement1496(m, n, x, v, u, w, p):
        # rubi.append(1496)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x)
def replacement1497(m, q, n, x, v, u, w, p, z):
        # rubi.append(1497)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p*ExpandToSum(z, x)**q, x)
def replacement1498(x, u, p):
        # rubi.append(1498)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1499(m, x, v, u, p):
        # rubi.append(1499)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x)
def replacement1500(m, n, x, v, u, w, p):
        # rubi.append(1500)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x)
def replacement1501(q, x, v, u, p):
        # rubi.append(1501)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
def replacement1502(x, u, p):
        # rubi.append(1502)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1503(c, m, x, u, p):
        # rubi.append(1503)
        return Int((c*x)**m*ExpandToSum(u, x)**p, x)
def replacement1504(q, x, v, u, p):
        # rubi.append(1504)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
def replacement1505(m, q, x, v, u, p):
        # rubi.append(1505)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
def replacement1506(m, q, x, v, u, w, p):
        # rubi.append(1506)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p*ExpandToSum(w, x)**q, x)
def replacement1507(m, q, r, x, v, u, p, z):
        # rubi.append(1507)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**r, x)
def replacement1508(x, u, p):
        # rubi.append(1508)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1509(x, p, m, u):
        # rubi.append(1509)
        return Int(x**m*ExpandToSum(u, x)**p, x)
def replacement1510(x, u, p):
        # rubi.append(1510)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1511(m, x, d, u, p):
        # rubi.append(1511)
        return Int((d*x)**m*ExpandToSum(u, x)**p, x)
def replacement1512(q, x, v, u, p):
        # rubi.append(1512)
        return Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x)
def replacement1513(q, x, v, u, p):
        # rubi.append(1513)
        return Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x)
def replacement1514(m, q, x, u, p, z):
        # rubi.append(1514)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x)
def replacement1515(m, q, x, u, p, z):
        # rubi.append(1515)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x)
def replacement1516(x, u, p):
        # rubi.append(1516)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1517(x, p, m, u):
        # rubi.append(1517)
        return Int(x**m*ExpandToSum(u, x)**p, x)
def replacement1518(x, p, u, z):
        # rubi.append(1518)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(z, x), x)
def replacement1519(m, x, u, p, z):
        # rubi.append(1519)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x), x)
def replacement1520(c, m, q, n, r, x, h, f, g, a, e):
        # rubi.append(1520)
        return -Simp((S(2)*a*g + S(4)*a*h*x**(n/S(4)) - S(2)*c*f*x**(n/S(2)))/(a*c*n*sqrt(a + c*x**n)), x)
def replacement1521(c, m, q, n, r, x, d, h, f, g, a, e):
        # rubi.append(1521)
        return Dist(x**(-m)*(d*x)**m, Int(x**m*(e + f*x**(n/S(4)) + g*x**(S(3)*n/S(4)) + h*x**n)/(a + c*x**n)**(S(3)/2), x), x)

def With1522(c, m, Pq, x, a, p, b):
        n = Denominator(p)
        # rubi.append(1522)
        return Dist(n/b, Subst(Int(x**(n*p + n + S(-1))*(-a*c/b + c*x**n/b)**m*ReplaceAll(Pq, Rule(x, -a/b + x**n/b)), x), x, (a + b*x)**(S(1)/n)), x)
def replacement1523(m, n, Pq, x, a, p, b):
        # rubi.append(1523)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*SubstFor(x**(m + S(1)), Pq, x), x), x, x**(m + S(1))), x)
def replacement1524(n, Pq, x, a, p, b):
        # rubi.append(1524)
        return Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x) + Simp((a + b*x**n)**(p + S(1))*Coeff(Pq, x, n + S(-1))/(b*n*(p + S(1))), x)
def replacement1525(c, m, n, Pq, x, a, p, b):
        # rubi.append(1525)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x)
def replacement1526(n, Pq, x, a, p, b):
        # rubi.append(1526)
        return Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x)
def replacement1527(m, n, Pq, x, a, p, b):
        # rubi.append(1527)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
def replacement1528(c, m, n, Pq, x, a, p, b):
        # rubi.append(1528)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
def replacement1529(m, n, Pq, x, a, p, b):
        # rubi.append(1529)
        return -Dist(S(1)/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*D(Pq, x), x), x) + Simp(Pq*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement1530(m, n, Pq, x, d, a, p, b):
        # rubi.append(1530)
        return Dist(S(1)/d, Int((d*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x), x)
def replacement1531(n, Pq, x, a, p, b):
        # rubi.append(1531)
        return Int(x*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x)

def With1532(m, n, Pq, x, a, p, b):
        u = IntHide(Pq*x**m, x)
        # rubi.append(1532)
        return -Dist(b*n*p, Int(x**(m + n)*(a + b*x**n)**(p + S(-1))*ExpandToSum(u*x**(-m + S(-1)), x), x), x) + Simp(u*(a + b*x**n)**p, x)

def With1533(c, m, n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        i = Symbol('i')
        # rubi.append(1533)
        return Dist(a*n*p, Int((c*x)**m*(a + b*x**n)**(p + S(-1))*Sum_doit(x**i*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q)), x), x) + Simp((c*x)**m*(a + b*x**n)**p*Sum_doit(x**(i + S(1))*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q)), x)

def With1534(n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        i = Symbol('i')
        # rubi.append(1534)
        return Dist(a*n*p, Int((a + b*x**n)**(p + S(-1))*Sum_doit(x**i*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q)), x), x) + Simp((a + b*x**n)**p*Sum_doit(x**(i + S(1))*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q)), x)

def With1535(n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        i = Symbol('i')
        if Equal(q, n + S(-1)):
            return True
        return False
def replacement1535(n, Pq, x, a, p, b):

        q = Expon(Pq, x)
        i = Symbol('i')
        # rubi.append(1535)
        return Dist(S(1)/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*Sum_doit(x**i*(i + n*(p + S(1)) + S(1))*Coeff(Pq, x, i), List(i, S(0), q + S(-1))), x), x) + Simp((a + b*x**n)**(p + S(1))*(a*Coeff(Pq, x, q) - b*x*ExpandToSum(Pq - x**q*Coeff(Pq, x, q), x))/(a*b*n*(p + S(1))), x)
def replacement1536(n, Pq, x, a, p, b):
        # rubi.append(1536)
        return Dist(S(1)/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*ExpandToSum(Pq*n*(p + S(1)) + D(Pq*x, x), x), x), x) - Simp(Pq*x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))), x)
def replacement1537(x, d, f, g, a, e, b):
        # rubi.append(1537)
        return -Simp((S(2)*a*f + S(4)*a*g*x - S(2)*b*e*x**S(2))/(S(4)*a*b*sqrt(a + b*x**S(4))), x)
def replacement1538(x, d, f, g, a, b):
        # rubi.append(1538)
        return -Simp((f + S(2)*g*x)/(S(2)*b*sqrt(a + b*x**S(4))), x)
def replacement1539(x, d, g, a, e, b):
        # rubi.append(1539)
        return -Simp(x*(S(2)*a*g - b*e*x)/(S(2)*a*b*sqrt(a + b*x**S(4))), x)
def replacement1540(x, h, f, a, e, b):
        # rubi.append(1540)
        return -Simp((f - S(2)*h*x**S(3))/(S(2)*b*sqrt(a + b*x**S(4))), x)
def replacement1541(x, h, a, e, b):
        # rubi.append(1541)
        return Simp(h*x**S(3)/(b*sqrt(a + b*x**S(4))), x)
def replacement1542(x, d, h, f, g, a, e, b):
        # rubi.append(1542)
        return -Simp((a*f - S(2)*a*h*x**S(3) - S(2)*b*d*x)/(S(2)*a*b*sqrt(a + b*x**S(4))), x)
def replacement1543(x, d, h, g, a, e, b):
        # rubi.append(1543)
        return Simp(x*(a*h*x**S(2) + b*d)/(a*b*sqrt(a + b*x**S(4))), x)

def With1544(n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        R = PolynomialRemainder(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        if GreaterEqual(q, n):
            return True
        return False
def replacement1544(n, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        R = PolynomialRemainder(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        # rubi.append(1544)
        return Dist(b**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)), Int((a + b*x**n)**(p + 1)*ExpandToSum(Q*a*n*(p + 1) + R*n*(p + 1) + D(R*x, x), x), x), x) - Simp(R*b**(-Floor((q - 1)/n) - 1)*x*(a + b*x**n)**(p + 1)/(a*n*(p + 1)), x)

def With1545(m, n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)
        R = PolynomialRemainder(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)
        # rubi.append(1545)
        return Dist(b**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)), Int(x**m*(a + b*x**n)**(p + 1)*ExpandToSum(Q*n*x**(-m)*(p + 1) + Sum_doit(x**(i - m)*(i + n*(p + 1) + 1)*Coeff(R, x, i)/a, List(i, 0, n - 1)), x), x), x) - Simp(R*b**(-Floor((q - 1)/n) - 1)*x*(a + b*x**n)**(p + 1)/(a**2*n*(p + 1)), x)

def With1546(m, n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        g = GCD(m + S(1), n)
        if Unequal(g, S(1)):
            return True
        return False
def replacement1546(m, n, Pq, x, a, p, b):

        g = GCD(m + S(1), n)
        # rubi.append(1546)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)
def replacement1547(x, B, A, a, b):
        # rubi.append(1547)
        return Dist(B**S(3)/b, Int(S(1)/(A**S(2) - A*B*x + B**S(2)*x**S(2)), x), x)

def With1548(x, B, A, a, b):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        # rubi.append(1548)
        return Dist(r/(S(3)*a*s), Int((r*(S(2)*A*s + B*r) + s*x*(-A*s + B*r))/(r**S(2) - r*s*x + s**S(2)*x**S(2)), x), x) - Dist(r*(-A*s + B*r)/(S(3)*a*s), Int(S(1)/(r + s*x), x), x)

def With1549(x, B, A, a, b):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        # rubi.append(1549)
        return -Dist(r/(S(3)*a*s), Int((r*(-S(2)*A*s + B*r) - s*x*(A*s + B*r))/(r**S(2) + r*s*x + s**S(2)*x**S(2)), x), x) + Dist(r*(A*s + B*r)/(S(3)*a*s), Int(S(1)/(r - s*x), x), x)
def replacement1550(x, B, A, a, C, b):
        # rubi.append(1550)
        return -Dist(C**S(2)/b, Int(S(1)/(B - C*x), x), x)

def With1551(x, B, A, a, C, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1551)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1552(x, B, a, C, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1552)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1553(x, A, a, C, b):
        q = a**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1553)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1554(x, B, A, a, C, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1554)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1555(x, B, a, C, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1555)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1556(x, A, a, C, b):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1556)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1557(x, B, A, a, C, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1557)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1558(x, B, a, C, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1558)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1559(x, A, a, C, b):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        # rubi.append(1559)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1560(x, B, A, a, C, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1560)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1561(x, B, a, C, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1561)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1562(x, A, a, C, b):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        # rubi.append(1562)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1563(x, B, A, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1563)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1564(x, B, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1564)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1565(x, A, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1565)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1566(x, B, A, a, C, b):
        q = Rt(a/b, S(3))
        # rubi.append(1566)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1567(x, B, a, C, b):
        q = Rt(a/b, S(3))
        # rubi.append(1567)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1568(x, A, a, C, b):
        q = Rt(a/b, S(3))
        # rubi.append(1568)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)

def With1569(x, B, A, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1569)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1570(x, B, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1570)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1571(x, A, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1571)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1572(x, B, A, a, C, b):
        q = Rt(-a/b, S(3))
        # rubi.append(1572)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1573(x, B, a, C, b):
        q = Rt(-a/b, S(3))
        # rubi.append(1573)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)

def With1574(x, A, a, C, b):
        q = Rt(-a/b, S(3))
        # rubi.append(1574)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
def replacement1575(x, B, A, a, C, b):
        # rubi.append(1575)
        return Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x) + Int((A + B*x)/(a + b*x**S(3)), x)
def replacement1576(x, B, a, C, b):
        # rubi.append(1576)
        return Dist(B, Int(x/(a + b*x**S(3)), x), x) + Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x)
def replacement1577(x, A, a, C, b):
        # rubi.append(1577)
        return Dist(A, Int(S(1)/(a + b*x**S(3)), x), x) + Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x)

def With1578(x, B, A, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1578)
        return Dist(q**S(2)/a, Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x), x)

def With1579(x, B, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1579)
        return Dist(C*q**S(3)/a, Int(x/(q**S(2) - q*x + x**S(2)), x), x)

def With1580(x, A, a, C, b):
        q = (a/b)**(S(1)/3)
        # rubi.append(1580)
        return Dist(q**S(2)/a, Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x), x)

def With1581(x, B, A, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1581)
        return Dist(q/a, Int((A*q + x*(A + B*q))/(q**S(2) + q*x + x**S(2)), x), x)

def With1582(x, B, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1582)
        return Dist(B*q**S(2)/a, Int(x/(q**S(2) + q*x + x**S(2)), x), x)

def With1583(x, A, a, C, b):
        q = (-a/b)**(S(1)/3)
        # rubi.append(1583)
        return Dist(A*q/a, Int((q + x)/(q**S(2) + q*x + x**S(2)), x), x)

def With1584(x, B, A, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (a/b)**(S(1)/3)
        if NonzeroQ(A - B*q + C*q**S(2)):
            return True
        return False
def replacement1584(x, B, A, a, C, b):

        q = (a/b)**(S(1)/3)
        # rubi.append(1584)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A + B*q - C*q**S(2)) - x*(A - B*q - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) + Dist(q*(A - B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)

def With1585(x, B, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (a/b)**(S(1)/3)
        if NonzeroQ(B*q - C*q**S(2)):
            return True
        return False
def replacement1585(x, B, a, C, b):

        q = (a/b)**(S(1)/3)
        # rubi.append(1585)
        return Dist(q/(S(3)*a), Int((q*(B*q - C*q**S(2)) + x*(B*q + S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) - Dist(q*(B*q - C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)

def With1586(x, A, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (a/b)**(S(1)/3)
        if NonzeroQ(A + C*q**S(2)):
            return True
        return False
def replacement1586(x, A, a, C, b):

        q = (a/b)**(S(1)/3)
        # rubi.append(1586)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - C*q**S(2)) - x*(A - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) + Dist(q*(A + C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)

def With1587(x, B, A, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (-a/b)**(S(1)/3)
        if NonzeroQ(A + B*q + C*q**S(2)):
            return True
        return False
def replacement1587(x, B, A, a, C, b):

        q = (-a/b)**(S(1)/3)
        # rubi.append(1587)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - B*q - C*q**S(2)) + x*(A + B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(A + B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)

def With1588(x, B, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (-a/b)**(S(1)/3)
        if NonzeroQ(B*q + C*q**S(2)):
            return True
        return False
def replacement1588(x, B, a, C, b):

        q = (-a/b)**(S(1)/3)
        # rubi.append(1588)
        return Dist(q/(S(3)*a), Int((-q*(B*q + C*q**S(2)) + x*(B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)

def With1589(x, A, a, C, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = (-a/b)**(S(1)/3)
        if NonzeroQ(A + C*q**S(2)):
            return True
        return False
def replacement1589(x, A, a, C, b):

        q = (-a/b)**(S(1)/3)
        # rubi.append(1589)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - C*q**S(2)) + x*(A - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(A + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)

def With1590(c, m, n, Pq, x, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = Sum_doit(c**(-ii)*(c*x)**(ii + m)*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        if SumQ(v):
            return True
        return False
def replacement1590(c, m, n, Pq, x, a, b):

        v = Sum_doit(c**(-ii)*(c*x)**(ii + m)*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        # rubi.append(1590)
        return Int(v, x)

def With1591(n, Pq, x, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = Sum_doit(x**ii*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        if SumQ(v):
            return True
        return False
def replacement1591(n, Pq, x, a, b):

        v = Sum_doit(x**ii*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        # rubi.append(1591)
        return Int(v, x)

def With1592(c, x, d, a, b):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        # rubi.append(1592)
        return Simp(S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(S(1) + sqrt(S(3))))), x) - Simp(S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(S(2) - sqrt(S(3)))*(r*x + s)*EllipticE(asin((r*x + s*(S(1) - sqrt(S(3))))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(r**S(2)*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))), x)

def With1593(c, x, d, a, b):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        # rubi.append(1593)
        return Dist(d/r, Int((r*x + s*(S(1) - sqrt(S(3))))/sqrt(a + b*x**S(3)), x), x) + Dist((c*r - d*s*(S(1) - sqrt(S(3))))/r, Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1594(c, x, d, a, b):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        # rubi.append(1594)
        return Simp(S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(S(1) - sqrt(S(3))))), x) + Simp(S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) - sqrt(S(3))))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticE(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(S(1) - sqrt(S(3))))), S(-7) + S(4)*sqrt(S(3)))/(r**S(2)*sqrt(-s*(r*x + s)/(r*x + s*(S(1) - sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))), x)

def With1595(c, x, d, a, b):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        # rubi.append(1595)
        return Dist(d/r, Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x), x) + Dist((c*r - d*s*(S(1) + sqrt(S(3))))/r, Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1596(c, x, d, a, b):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        # rubi.append(1596)
        return Simp(d*s**S(3)*x*(S(1) + sqrt(S(3)))*sqrt(a + b*x**S(6))/(S(2)*a*r**S(2)*(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), x) - Simp(S(3)**(S(1)/4)*d*s*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticE(acos((r*x**S(2)*(S(1) - sqrt(S(3))) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(2)*r**S(2)*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6))), x)

def With1597(c, x, d, a, b):
        q = Rt(b/a, S(3))
        # rubi.append(1597)
        return Dist(d/(S(2)*q**S(2)), Int((S(2)*q**S(2)*x**S(4) - sqrt(S(3)) + S(1))/sqrt(a + b*x**S(6)), x), x) + Dist((S(2)*c*q**S(2) - d*(S(1) - sqrt(S(3))))/(S(2)*q**S(2)), Int(S(1)/sqrt(a + b*x**S(6)), x), x)
def replacement1598(c, x, d, a, b):
        # rubi.append(1598)
        return -Simp(c*d*x**S(3)*sqrt(-(c - d*x**S(2))**S(2)/(c*d*x**S(2)))*sqrt(-d**S(2)*(a + b*x**S(8))/(b*c**S(2)*x**S(4)))*EllipticF(asin(sqrt((sqrt(S(2))*c**S(2) + S(2)*c*d*x**S(2) + sqrt(S(2))*d**S(2)*x**S(4))/(c*d*x**S(2)))/S(2)), S(-2) + S(2)*sqrt(S(2)))/(sqrt(sqrt(S(2)) + S(2))*sqrt(a + b*x**S(8))*(c - d*x**S(2))), x)
def replacement1599(c, x, d, a, b):
        # rubi.append(1599)
        return -Dist((-c*Rt(b/a, S(4)) + d)/(S(2)*Rt(b/a, S(4))), Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x) + Dist((c*Rt(b/a, S(4)) + d)/(S(2)*Rt(b/a, S(4))), Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x)
def replacement1600(n, Pq, x, a, b):
        # rubi.append(1600)
        return Dist(Coeff(Pq, x, S(0)), Int(S(1)/(x*sqrt(a + b*x**n)), x), x) + Int(ExpandToSum((Pq - Coeff(Pq, x, S(0)))/x, x)/sqrt(a + b*x**n), x)

def With1601(c, m, n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        # rubi.append(1601)
        return Int(Sum_doit(c**(-j)*(c*x)**(j + m)*(a + b*x**n)**p*Sum_doit(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)

def With1602(n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        # rubi.append(1602)
        return Int(Sum_doit(x**j*(a + b*x**n)**p*Sum_doit(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)
def replacement1603(n, Pq, x, a, p, b):
        # rubi.append(1603)
        return Dist(Coeff(Pq, x, n + S(-1)), Int(x**(n + S(-1))*(a + b*x**n)**p, x), x) + Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x)
def replacement1604(c, m, n, Pq, x, a, b):
        # rubi.append(1604)
        return Int(ExpandIntegrand(Pq*(c*x)**m/(a + b*x**n), x), x)
def replacement1605(n, Pq, x, a, b):
        # rubi.append(1605)
        return Int(ExpandIntegrand(Pq/(a + b*x**n), x), x)

def With1606(c, m, n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        Pq0 = Coeff(Pq, x, S(0))
        if NonzeroQ(Pq0):
            return True
        return False
def replacement1606(c, m, n, Pq, x, a, p, b):

        Pq0 = Coeff(Pq, x, S(0))
        # rubi.append(1606)
        return Dist(S(1)/(S(2)*a*c*(m + S(1))), Int((c*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(-S(2)*Pq0*b*x**(n + S(-1))*(m + n*(p + S(1)) + S(1)) + S(2)*a*(Pq - Pq0)*(m + S(1))/x, x), x), x) + Simp(Pq0*(c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))), x)

def With1607(c, m, n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if And(NonzeroQ(m + n*p + q + S(1)), GreaterEqual(-n + q, S(0)), Or(IntegerQ(S(2)*p), IntegerQ(p + (q + S(1))/(S(2)*n)))):
            return True
        return False
def replacement1607(c, m, n, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1607)
        return Dist(1/(b*(m + n*p + q + 1)), Int((c*x)**m*(a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(m - n + q + 1) + b*(Pq - Pqq*x**q)*(m + n*p + q + 1), x), x), x) + Simp(Pqq*c**(n - q - 1)*(c*x)**(m - n + q + 1)*(a + b*x**n)**(p + 1)/(b*(m + n*p + q + 1)), x)

def With1608(n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if And(NonzeroQ(n*p + q + S(1)), GreaterEqual(-n + q, S(0)), Or(IntegerQ(S(2)*p), IntegerQ(p + (q + S(1))/(S(2)*n)))):
            return True
        return False
def replacement1608(n, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1608)
        return Dist(1/(b*(n*p + q + 1)), Int((a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(-n + q + 1) + b*(Pq - Pqq*x**q)*(n*p + q + 1), x), x), x) + Simp(Pqq*x**(-n + q + 1)*(a + b*x**n)**(p + 1)/(b*(n*p + q + 1)), x)

def With1609(m, n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        # rubi.append(1609)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x)

def With1610(c, m, n, Pq, x, a, p, b):
        g = Denominator(m)
        q = Expon(Pq, x)
        # rubi.append(1610)
        return -Dist(g/c, Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*c**(-n)*x**(-g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/c)), x), x), x, (c*x)**(-S(1)/g)), x)

def With1611(c, m, n, Pq, x, a, p, b):
        q = Expon(Pq, x)
        # rubi.append(1611)
        return -Dist((c*x)**m*(S(1)/x)**m, Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x), x)

def With1612(m, n, Pq, x, a, p, b):
        g = Denominator(n)
        # rubi.append(1612)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)

def With1613(n, Pq, x, a, p, b):
        g = Denominator(n)
        # rubi.append(1613)
        return Dist(g, Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
def replacement1614(c, m, n, Pq, x, a, p, b):
        # rubi.append(1614)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
def replacement1615(m, n, Pq, x, a, p, b):
        # rubi.append(1615)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
def replacement1616(c, m, n, Pq, x, a, p, b):
        # rubi.append(1616)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
def replacement1617(m, n, x, B, A, a, p, b):
        # rubi.append(1617)
        return Dist(A, Int((a + b*x**n)**p, x), x) + Dist(B, Int(x**m*(a + b*x**n)**p, x), x)
def replacement1618(c, m, n, Pq, x, a, p, b):
        # rubi.append(1618)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x)
def replacement1619(n, Pq, x, a, p, b):
        # rubi.append(1619)
        return Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x)
def replacement1620(m, n, Pq, x, v, u, a, p, b):
        # rubi.append(1620)
        return Dist(u**m*v**(-m)/Coeff(v, x, S(1)), Subst(Int(x**m*(a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v), x)
def replacement1621(n, Pq, x, v, a, p, b):
        # rubi.append(1621)
        return Dist(S(1)/Coeff(v, x, S(1)), Subst(Int((a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v), x)
def replacement1622(c, m, n, Pq, x, b2, b1, a1, p, a2):
        # rubi.append(1622)
        return Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
def replacement1623(n, Pq, x, b2, b1, a1, p, a2):
        # rubi.append(1623)
        return Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
def replacement1624(c, m, n, Pq, x, b2, b1, a1, p, a2):
        # rubi.append(1624)
        return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)
def replacement1625(n, Pq, x, b2, b1, a1, p, a2):
        # rubi.append(1625)
        return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)
def replacement1626(c, n, x, d, f, g, n2, a, p, e, b):
        # rubi.append(1626)
        return Simp(e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c), x)
def replacement1627(c, n, x, d, g, n2, a, p, e, b):
        # rubi.append(1627)
        return Simp(e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c), x)
def replacement1628(c, m, n, x, d, h, f, g, n2, a, p, e, b):
        # rubi.append(1628)
        return Simp(e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))), x)
def replacement1629(c, m, n, x, d, h, g, n2, a, p, e, b):
        # rubi.append(1629)
        return Simp(e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))), x)
def replacement1630(c, m, q, n, x, d, B, A, a, p, b):
        # rubi.append(1630)
        return Dist(A, Int((a + b*x**n)**p*(c + d*x**n)**q, x), x) + Dist(B, Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)

def With1631(c, q, n, x, d, Px, a, p, b):
        k = Denominator(n)
        # rubi.append(1631)
        return Dist(k/d, Subst(Int(SimplifyIntegrand(x**(k + S(-1))*(a + b*x**(k*n))**p*ReplaceAll(Px, Rule(x, -c/d + x**k/d))**q, x), x), x, (c + d*x)**(S(1)/k)), x)
def replacement1632(c, m, n, Pq, x, n2, a, p, b):
        # rubi.append(1632)
        return Dist(S(1)/n, Subst(Int((a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
def replacement1633(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1633)
        return Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1634(c, n, Pq, x, n2, a, p, b):
        # rubi.append(1634)
        return Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1635(c, n, x, d, f, n2, a, p, e, b):
        # rubi.append(1635)
        return Simp(d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a, x)
def replacement1636(c, n, x, d, f, n2, a, p, b):
        # rubi.append(1636)
        return Simp(d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a, x)
def replacement1637(c, m, n, x, d, f, g, n2, a, p, e, b):
        # rubi.append(1637)
        return Simp(d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))), x)
def replacement1638(c, m, n, x, d, f, g, n2, a, p, b):
        # rubi.append(1638)
        return Simp(d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))), x)
def replacement1639(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1639)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int(Pq*(d*x)**m*(b + S(2)*c*x**n)**(S(2)*p), x), x)
def replacement1640(c, n, Pq, x, n2, a, p, b):
        # rubi.append(1640)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int(Pq*(b + S(2)*c*x**n)**(S(2)*p), x), x)
def replacement1641(c, m, n, Pq, x, n2, a, p, b):
        # rubi.append(1641)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
def replacement1642(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1642)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1643(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1643)
        return Dist(S(1)/d, Int((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x), x)
def replacement1644(c, n, Pq, x, n2, a, p, b):
        # rubi.append(1644)
        return Int(x*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x)
def replacement1645(c, n, x, d, n3, f, g, n2, a, p, e, b):
        # rubi.append(1645)
        return Simp(x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
def replacement1646(c, n, x, d, n3, f, g, n2, a, p, b):
        # rubi.append(1646)
        return Simp(x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
def replacement1647(c, n, x, d, n3, g, n2, a, p, e, b):
        # rubi.append(1647)
        return Simp(x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
def replacement1648(c, n, x, d, n3, g, n2, a, p, b):
        # rubi.append(1648)
        return Simp(x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
def replacement1649(c, m, q, n, r, x, h, f, g, n2, a, e, s, b):
        # rubi.append(1649)
        return -Simp((S(2)*c*x**n*(-b*g + S(2)*c*f) + S(2)*c*(-S(2)*a*g + b*f) + S(2)*h*x**(n/S(2))*(-S(4)*a*c + b**S(2)))/(c*n*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**n + c*x**(S(2)*n))), x)
def replacement1650(c, m, q, n, r, x, d, h, f, g, n2, a, e, s, b):
        # rubi.append(1650)
        return Dist(x**(-m)*(d*x)**m, Int(x**m*(e + f*x**(n/S(2)) + g*x**(S(3)*n/S(2)) + h*x**(S(2)*n))/(a + b*x**n + c*x**(S(2)*n))**(S(3)/2), x), x)

def With1651(c, n, Pq, x, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        i = Symbol('i')
        if Less(q, S(2)*n):
            return True
        return False
def replacement1651(c, n, Pq, x, n2, a, p, b):

        q = Expon(Pq, x)
        i = Symbol('i')
        # rubi.append(1651)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum_doit(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i))*(i + n*(S(2)*p + S(3)) + S(1)) + x**i*(-a*b*(i + S(1))*Coeff(Pq, x, i + n) + (-S(2)*a*c*(i + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(i + n*(p + S(1)) + S(1)))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1))), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum_doit(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i)) + x**i*(-a*b*Coeff(Pq, x, i + n) + (-S(2)*a*c + b**S(2))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1652(c, x, d, f, g, a, e, b):
        # rubi.append(1652)
        return -Simp((c*x**S(2)*(-b*f + S(2)*c*e) + c*(-S(2)*a*f + b*e) + g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1653(c, x, d, f, g, a, b):
        # rubi.append(1653)
        return Simp((S(2)*a*c*f + b*c*f*x**S(2) - g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1654(c, x, d, g, a, e, b):
        # rubi.append(1654)
        return -Simp((b*c*e + S(2)*c**S(2)*e*x**S(2) + g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1655(c, x, h, f, g, a, e, b):
        # rubi.append(1655)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1656(c, x, h, g, a, e, b):
        # rubi.append(1656)
        return Simp(h*x**S(3)/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1657(c, x, d, h, f, g, a, e, b):
        # rubi.append(1657)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1658(c, x, d, h, f, a, e, b):
        # rubi.append(1658)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)

def With1659(c, n, Pq, x, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        if GreaterEqual(q, S(2)*n):
            return True
        return False
def replacement1659(c, n, Pq, x, n2, a, p, b):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        # rubi.append(1659)
        return Dist((b*c)**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)*(-4*a*c + b**2)), Int((a + b*x**n + c*x**(2*n))**(p + 1)*ExpandToSum(Q*a*n*(p + 1)*(-4*a*c + b**2) + Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i))*(i + n*(2*p + 3) + 1) + x**i*(-a*b*(i + 1)*Coeff(R, x, i + n) + (-2*a*c*(i + 2*n*(p + 1) + 1) + b**2*(i + n*(p + 1) + 1))*Coeff(R, x, i)), List(i, 0, n - 1)), x), x), x) - Simp(x*(b*c)**(-Floor((q - 1)/n) - 1)*(a + b*x**n + c*x**(2*n))**(p + 1)*Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-2*a*c + b**2)*Coeff(R, x, i)), List(i, 0, n - 1))/(a*n*(p + 1)*(-4*a*c + b**2)), x)

def With1660(c, m, n, Pq, x, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        if GreaterEqual(q, S(2)*n):
            return True
        return False
def replacement1660(c, m, n, Pq, x, n2, a, p, b):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        # rubi.append(1660)
        return Dist((b*c)**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)*(-4*a*c + b**2)), Int(x**m*(a + b*x**n + c*x**(2*n))**(p + 1)*ExpandToSum(Q*n*x**(-m)*(p + 1)*(-4*a*c + b**2) + Sum_doit(c*x**(i - m + n)*(-2*Coeff(R, x, i + n) + b*Coeff(R, x, i)/a)*(i + n*(2*p + 3) + 1) + x**(i - m)*(-b*(i + 1)*Coeff(R, x, i + n) + (-2*c*(i + 2*n*(p + 1) + 1) + b**2*(i + n*(p + 1) + 1)/a)*Coeff(R, x, i)), List(i, 0, n - 1)), x), x), x) - Simp(x*(b*c)**(-Floor((q - 1)/n) - 1)*(a + b*x**n + c*x**(2*n))**(p + 1)*Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-2*a*c + b**2)*Coeff(R, x, i)), List(i, 0, n - 1))/(a**2*n*(p + 1)*(-4*a*c + b**2)), x)

def With1661(c, m, n, Pq, x, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        g = GCD(m + S(1), n)
        if Unequal(g, S(1)):
            return True
        return False
def replacement1661(c, m, n, Pq, x, n2, a, p, b):

        g = GCD(m + S(1), n)
        # rubi.append(1661)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g) + c*x**(S(2)*n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)
def replacement1662(c, m, n, Pq, x, d, n2, a, b):
        # rubi.append(1662)
        return Int(ExpandIntegrand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)
def replacement1663(c, n, Pq, x, n2, a, b):
        # rubi.append(1663)
        return Int(ExpandIntegrand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x)

def With1664(c, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if Equal(S(2)*p + q + S(1), S(0)):
            return True
        return False
def replacement1664(c, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1664)
        return Dist(1/2, Int((a + b*x + c*x**2)**p*ExpandToSum(2*Pq - Pqq*c**p*(b + 2*c*x)*(a + b*x + c*x**2)**(-p - 1), x), x), x) + Simp(Pqq*c**p*log(a + b*x + c*x**2)/2, x)

def With1665(c, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if Equal(S(2)*p + q + S(1), S(0)):
            return True
        return False
def replacement1665(c, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1665)
        return Int((a + b*x + c*x**2)**p*ExpandToSum(Pq - Pqq*c**(p + 1/2)*(a + b*x + c*x**2)**(-p - 1/2), x), x) + Simp(Pqq*c**p*atanh((b + 2*c*x)/(2*sqrt(a + b*x + c*x**2)*Rt(c, 2))), x)

def With1666(c, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if Equal(S(2)*p + q + S(1), S(0)):
            return True
        return False
def replacement1666(c, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1666)
        return Int((a + b*x + c*x**2)**p*ExpandToSum(Pq - Pqq*(-c)**(p + 1/2)*(a + b*x + c*x**2)**(-p - 1/2), x), x) - Simp(Pqq*(-c)**p*ArcTan((b + 2*c*x)/(2*sqrt(a + b*x + c*x**2)*Rt(-c, 2))), x)

def With1667(c, m, n, Pq, x, d, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if And(GreaterEqual(q, S(2)*n), Unequal(m + S(2)*n*p + q + S(1), S(0)), Or(IntegerQ(S(2)*p), And(Equal(n, S(1)), IntegerQ(S(4)*p)), IntegerQ(p + (q + S(1))/(S(2)*n)))):
            return True
        return False
def replacement1667(c, m, n, Pq, x, d, n2, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1667)
        return Int((d*x)**m*(a + b*x**n + c*x**(2*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-2*n + q)*(m - 2*n + q + 1) + b*x**(-n + q)*(m + n*(p - 1) + q + 1))/(c*(m + 2*n*p + q + 1)), x), x) + Simp(Pqq*d**(2*n - q - 1)*(d*x)**(m - 2*n + q + 1)*(a + b*x**n + c*x**(2*n))**(p + 1)/(c*(m + 2*n*p + q + 1)), x)

def With1668(c, n, Pq, x, n2, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if And(GreaterEqual(q, S(2)*n), Unequal(S(2)*n*p + q + S(1), S(0)), Or(IntegerQ(S(2)*p), And(Equal(n, S(1)), IntegerQ(S(4)*p)), IntegerQ(p + (q + S(1))/(S(2)*n)))):
            return True
        return False
def replacement1668(c, n, Pq, x, n2, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1668)
        return Int((a + b*x**n + c*x**(2*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-2*n + q)*(-2*n + q + 1) + b*x**(-n + q)*(n*(p - 1) + q + 1))/(c*(2*n*p + q + 1)), x), x) + Simp(Pqq*x**(-2*n + q + 1)*(a + b*x**n + c*x**(2*n))**(p + 1)/(c*(2*n*p + q + 1)), x)

def With1669(c, m, n, Pq, x, d, n2, a, p, b):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        # rubi.append(1669)
        return Int(Sum_doit(d**(-j)*(d*x)**(j + m)*(a + b*x**n + c*x**(S(2)*n))**p*Sum_doit(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)

def With1670(c, n, Pq, x, n2, a, p, b):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        # rubi.append(1670)
        return Int(Sum_doit(x**j*(a + b*x**n + c*x**(S(2)*n))**p*Sum_doit(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)
def replacement1671(c, m, n, Pq, x, d, n2, a, b):
        # rubi.append(1671)
        return Int(RationalFunctionExpand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)
def replacement1672(c, n, Pq, x, n2, a, b):
        # rubi.append(1672)
        return Int(RationalFunctionExpand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x)

def With1673(c, m, n, Pq, x, n2, a, p, b):
        q = Expon(Pq, x)
        # rubi.append(1673)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x)

def With1674(c, m, n, Pq, x, d, n2, a, p, b):
        g = Denominator(m)
        q = Expon(Pq, x)
        # rubi.append(1674)
        return -Dist(g/d, Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*d**(-n)*x**(-g*n) + c*d**(-S(2)*n)*x**(-S(2)*g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/d)), x), x), x, (d*x)**(-S(1)/g)), x)

def With1675(c, m, n, Pq, x, d, n2, a, p, b):
        q = Expon(Pq, x)
        # rubi.append(1675)
        return -Dist((d*x)**m*(S(1)/x)**m, Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x), x)

def With1676(c, m, n, Pq, x, n2, a, p, b):
        g = Denominator(n)
        # rubi.append(1676)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)

def With1677(c, n, Pq, x, n2, a, p, b):
        g = Denominator(n)
        # rubi.append(1677)
        return Dist(g, Subst(Int(x**(g + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
def replacement1678(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1678)
        return Dist(d**(m + S(-1)/2)*sqrt(d*x)/sqrt(x), Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1679(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1679)
        return Dist(d**(m + S(1)/2)*sqrt(x)/sqrt(d*x), Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1680(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1680)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1681(c, m, n, Pq, x, n2, a, p, b):
        # rubi.append(1681)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
def replacement1682(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1682)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)

def With1683(c, m, n, Pq, x, d, n2, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(1683)
        return Dist(S(2)*c/q, Int(Pq*(d*x)**m/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int(Pq*(d*x)**m/(b + S(2)*c*x**n + q), x), x)

def With1684(c, n, Pq, x, n2, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(1684)
        return Dist(S(2)*c/q, Int(Pq/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int(Pq/(b + S(2)*c*x**n + q), x), x)
def replacement1685(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1685)
        return Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1686(c, n, Pq, x, n2, a, p, b):
        # rubi.append(1686)
        return Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
def replacement1687(c, m, n, Pq, x, d, n2, a, p, b):
        # rubi.append(1687)
        return Int(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x)
def replacement1688(c, n, Pq, x, n2, a, p, b):
        # rubi.append(1688)
        return Int(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x)
def replacement1689(c, m, n, Pq, x, v, u, n2, a, p, b):
        # rubi.append(1689)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v), x)
def replacement1690(c, n, Pq, x, v, n2, a, p, b):
        # rubi.append(1690)
        return Dist(S(1)/Coefficient(v, x, S(1)), Subst(Int((a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v), x)
def replacement1691(j, n, x, a, p, b):
        # rubi.append(1691)
        return Simp(x**(S(1) - n)*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
def replacement1692(j, n, x, a, p, b):
        # rubi.append(1692)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(S(1) - j)*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1693(j, n, x, a, p, b):
        # rubi.append(1693)
        return -Dist(b*(-j + n*p + n + S(1))/(a*(j*p + S(1))), Int(x**(-j + n)*(a*x**j + b*x**n)**p, x), x) + Simp(x**(S(1) - j)*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + S(1))), x)
def replacement1694(j, n, x, a, p, b):
        # rubi.append(1694)
        return -Dist(b*p*(-j + n)/(j*p + S(1)), Int(x**n*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(j*p + S(1)), x)
def replacement1695(j, n, x, a, p, b):
        # rubi.append(1695)
        return Dist(a*p*(-j + n)/(n*p + S(1)), Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(n*p + S(1)), x)
def replacement1696(j, n, x, a, p, b):
        # rubi.append(1696)
        return -Dist((j*p + j - n + S(1))/(b*(-j + n)*(p + S(1))), Int(x**(-n)*(a*x**j + b*x**n)**(p + S(1)), x), x) + Simp(x**(S(1) - n)*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
def replacement1697(j, n, x, a, p, b):
        # rubi.append(1697)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(S(1) - j)*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1698(j, n, x, a, p, b):
        # rubi.append(1698)
        return Dist(a, Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(p*(-j + n)), x)
def replacement1699(x, a, n, b):
        # rubi.append(1699)
        return Dist(S(2)/(S(2) - n), Subst(Int(S(1)/(-a*x**S(2) + S(1)), x), x, x/sqrt(a*x**S(2) + b*x**n)), x)
def replacement1700(j, n, x, a, p, b):
        # rubi.append(1700)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(S(1) - j)*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1701(j, n, x, a, b):
        # rubi.append(1701)
        return -Dist(a*(-j + S(2)*n + S(-2))/(b*(n + S(-2))), Int(x**(j - n)/sqrt(a*x**j + b*x**n), x), x) + Simp(-S(2)*x**(S(1) - n)*sqrt(a*x**j + b*x**n)/(b*(n + S(-2))), x)
def replacement1702(j, n, x, a, p, b):
        # rubi.append(1702)
        return Dist(x**(-j*FracPart(p))*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p), Int(x**(j*p)*(a + b*x**(-j + n))**p, x), x)
def replacement1703(j, n, x, u, a, p, b):
        # rubi.append(1703)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a*x**j + b*x**n)**p, x), x, u), x)
def replacement1704(m, j, n, x, a, p, b):
        # rubi.append(1704)
        return Dist(S(1)/n, Subst(Int((a*x**(j/n) + b*x)**p, x), x, x**n), x)
def replacement1705(c, m, j, n, x, a, p, b):
        # rubi.append(1705)
        return -Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1706(c, m, j, n, x, a, p, b):
        # rubi.append(1706)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1707(c, m, j, n, x, a, p, b):
        # rubi.append(1707)
        return -Dist(b*c**(j - n)*(-j + m + n*p + n + S(1))/(a*(j*p + m + S(1))), Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))), x)
def replacement1708(c, m, j, n, x, a, p, b):
        # rubi.append(1708)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1709(m, j, n, x, a, p, b):
        # rubi.append(1709)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p, x), x, x**n), x)
def replacement1710(c, m, j, n, x, a, p, b):
        # rubi.append(1710)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1711(c, m, j, n, x, a, p, b):
        # rubi.append(1711)
        return -Dist(b*c**(-n)*p*(-j + n)/(j*p + m + S(1)), Int((c*x)**(m + n)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(j*p + m + S(1))), x)
def replacement1712(c, m, j, n, x, a, p, b):
        # rubi.append(1712)
        return Dist(a*c**(-j)*p*(-j + n)/(m + n*p + S(1)), Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(m + n*p + S(1))), x)
def replacement1713(c, m, j, n, x, a, p, b):
        # rubi.append(1713)
        return -Dist(c**n*(j*p + j + m - n + S(1))/(b*(-j + n)*(p + S(1))), Int((c*x)**(m - n)*(a*x**j + b*x**n)**(p + S(1)), x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
def replacement1714(c, m, j, n, x, a, p, b):
        # rubi.append(1714)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1715(c, m, j, n, x, a, p, b):
        # rubi.append(1715)
        return -Dist(a*c**(-j + n)*(j*p + j + m - n + S(1))/(b*(m + n*p + S(1))), Int((c*x)**(j + m - n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))), x)
def replacement1716(c, m, j, n, x, a, p, b):
        # rubi.append(1716)
        return -Dist(b*c**(j - n)*(-j + m + n*p + n + S(1))/(a*(j*p + m + S(1))), Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))), x)
def replacement1717(m, j, n, x, a, p, b):
        # rubi.append(1717)
        return Dist(S(1)/(m + S(1)), Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1))), x)
def replacement1718(c, m, j, n, x, a, p, b):
        # rubi.append(1718)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1719(c, m, j, n, x, a, p, b):
        # rubi.append(1719)
        return Dist(a*c**(-j), Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*p*(-j + n)), x)
def replacement1720(m, j, n, x, a, b):
        # rubi.append(1720)
        return Dist(-S(2)/(-j + n), Subst(Int(S(1)/(-a*x**S(2) + S(1)), x), x, x**(j/S(2))/sqrt(a*x**j + b*x**n)), x)
def replacement1721(c, m, j, n, x, a, p, b):
        # rubi.append(1721)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
def replacement1722(c, m, j, n, x, a, p, b):
        # rubi.append(1722)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1723(c, m, j, n, x, a, p, b):
        # rubi.append(1723)
        return Dist(c**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(c*x)**FracPart(m)*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p), Int(x**(j*p + m)*(a + b*x**(-j + n))**p, x), x)
def replacement1724(m, j, n, x, v, u, a, p, b):
        # rubi.append(1724)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a*x**j + b*x**n)**p, x), x, v), x)
def replacement1725(c, m, q, j, n, x, d, k, a, p, b):
        # rubi.append(1725)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(c + d*x)**q*(a*x**(j/n) + b*x**(k/n))**p, x), x, x**n), x)
def replacement1726(c, m, q, n, j, x, d, k, a, p, e, b):
        # rubi.append(1726)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x), x)
def replacement1727(c, m, j, n, x, d, jn, a, p, e, b):
        # rubi.append(1727)
        return Simp(c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))), x)
def replacement1728(c, m, j, n, x, d, jn, a, p, e, b):
        # rubi.append(1728)
        return -Dist(e**j*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(a*b*n*(p + S(1))), Int((e*x)**(-j + m)*(a*x**j + b*x**(j + n))**(p + S(1)), x), x) - Simp(e**(j + S(-1))*(e*x)**(-j + m + S(1))*(-a*d + b*c)*(a*x**j + b*x**(j + n))**(p + S(1))/(a*b*n*(p + S(1))), x)
def replacement1729(c, m, j, n, x, d, jn, a, p, e, b):
        # rubi.append(1729)
        return Dist(e**(-n)*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(a*(j*p + m + S(1))), Int((e*x)**(m + n)*(a*x**j + b*x**(j + n))**p, x), x) + Simp(c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))), x)
def replacement1730(c, m, j, n, x, d, jn, a, p, e, b):
        # rubi.append(1730)
        return -Dist((a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(b*(m + n + p*(j + n) + S(1))), Int((e*x)**m*(a*x**j + b*x**(j + n))**p, x), x) + Simp(d*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(b*(m + n + p*(j + n) + S(1))), x)
def replacement1731(c, m, q, n, j, x, d, k, a, p, b):
        # rubi.append(1731)
        return Dist(S(1)/(m + S(1)), Subst(Int((c + d*x**(n/(m + S(1))))**q*(a*x**(j/(m + S(1))) + b*x**(k/(m + S(1))))**p, x), x, x**(m + S(1))), x)
def replacement1732(c, m, q, n, j, x, d, k, a, p, e, b):
        # rubi.append(1732)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x), x)
def replacement1733(c, m, j, n, q, x, d, jn, a, p, e, b):
        # rubi.append(1733)
        return Dist(e**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(e*x)**FracPart(m)*(a + b*x**n)**(-FracPart(p))*(a*x**j + b*x**(j + n))**FracPart(p), Int(x**(j*p + m)*(a + b*x**n)**p*(c + d*x**n)**q, x), x)

def With1734(j, n, Pq, x, a, p, b):
        d = Denominator(n)
        # rubi.append(1734)
        return Dist(d, Subst(Int(x**(d + S(-1))*(a*x**(d*j) + b*x**(d*n))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(d*n))), x), x, x**(S(1)/d)), x)
def replacement1735(m, j, n, Pq, x, a, p, b):
        # rubi.append(1735)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
def replacement1736(c, m, j, n, Pq, x, a, p, b):
        # rubi.append(1736)
        return Dist(c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m)), Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1737(c, m, j, n, Pq, x, a, p, b):
        # rubi.append(1737)
        return Dist(x**(-m)*(c*x)**m, Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)

def With1738(m, j, n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        g = GCD(m + S(1), n)
        if Unequal(g, S(1)):
            return True
        return False
def replacement1738(m, j, n, Pq, x, a, p, b):

        g = GCD(m + S(1), n)
        # rubi.append(1738)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a*x**(j/g) + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)

def With1739(c, m, j, n, Pq, x, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        if And(Greater(q, n + S(-1)), Unequal(m + n*p + q + S(1), S(0)), Or(IntegerQ(S(2)*p), IntegerQ(p + (q + S(1))/(S(2)*n)))):
            return True
        return False
def replacement1739(c, m, j, n, Pq, x, a, p, b):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        # rubi.append(1739)
        return Int((c*x)**m*(a*x**j + b*x**n)**p*ExpandToSum(Pq - Pqq*a*x**(-n + q)*(m - n + q + 1)/(b*(m + n*p + q + 1)) - Pqq*x**q, x), x) + Simp(Pqq*c**(n - q - 1)*(c*x)**(m - n + q + 1)*(a*x**j + b*x**n)**(p + 1)/(b*(m + n*p + q + 1)), x)
def replacement1740(m, j, n, Pq, x, a, p, b):
        # rubi.append(1740)
        return Dist(S(1)/(m + S(1)), Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
def replacement1741(c, m, j, n, Pq, x, a, p, b):
        # rubi.append(1741)
        return Dist(c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m)), Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1742(c, m, j, n, Pq, x, a, p, b):
        # rubi.append(1742)
        return Dist(x**(-m)*(c*x)**m, Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
def replacement1743(c, m, j, n, Pq, x, a, p, b):
        # rubi.append(1743)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a*x**j + b*x**n)**p, x), x)
def replacement1744(j, n, Pq, x, a, p, b):
        # rubi.append(1744)
        return Int(ExpandIntegrand(Pq*(a*x**j + b*x**n)**p, x), x)
def replacement1745(x, d, a, p, b):
        # rubi.append(1745)
        return Dist(S(3)**(-S(3)*p)*a**(-S(2)*p), Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x), x)
def replacement1746(x, d, a, p, b):
        # rubi.append(1746)
        return Int(ExpandToSum((a + b*x + d*x**S(3))**p, x), x)

def With1747(x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + b*x + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1747(x, d, a, p, b):

        u = Factor(a + b*x + d*x**S(3))
        # rubi.append(1747)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)

def With1748(x, d, a, p, b):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        # rubi.append(1748)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
def replacement1749(x, d, a, p, b):
        # rubi.append(1749)
        return Dist((S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p, Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x), x)

def With1750(x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1750(x, d, a, p, b):

        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        # rubi.append(1750)
        return Dist((a + b*x + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)

def With1751(x, d, a, p, b):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        # rubi.append(1751)
        return Dist((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p, Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
def replacement1752(m, x, d, f, a, p, e, b):
        # rubi.append(1752)
        return Dist(S(3)**(-S(3)*p)*a**(-S(2)*p), Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x), x)
def replacement1753(m, x, d, f, a, p, e, b):
        # rubi.append(1753)
        return Int(ExpandIntegrand((e + f*x)**m*(a + b*x + d*x**S(3))**p, x), x)

def With1754(m, x, d, f, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + b*x + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1754(m, x, d, f, a, p, e, b):

        u = Factor(a + b*x + d*x**S(3))
        # rubi.append(1754)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)

def With1755(m, x, d, f, a, p, e, b):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        # rubi.append(1755)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
def replacement1756(m, x, d, f, a, p, e, b):
        # rubi.append(1756)
        return Dist((S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p, Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x), x)

def With1757(m, x, d, f, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1757(m, x, d, f, a, p, e, b):

        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        # rubi.append(1757)
        return Dist((a + b*x + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)

def With1758(m, x, d, f, a, p, e, b):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        # rubi.append(1758)
        return Dist((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p, Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
def replacement1759(c, x, d, a, p):
        # rubi.append(1759)
        return -Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x), x)
def replacement1760(c, x, d, a, p):
        # rubi.append(1760)
        return Int(ExpandToSum((a + c*x**S(2) + d*x**S(3))**p, x), x)

def With1761(c, x, d, a, p):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + c*x**S(2) + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1761(c, x, d, a, p):

        u = Factor(a + c*x**S(2) + d*x**S(3))
        # rubi.append(1761)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)

def With1762(c, x, d, a, p):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        # rubi.append(1762)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
def replacement1763(c, x, d, a, p):
        # rubi.append(1763)
        return Dist((c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p, Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x), x)

def With1764(c, x, d, a, p):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1764(c, x, d, a, p):

        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        # rubi.append(1764)
        return Dist((a + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)

def With1765(c, x, d, a, p):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        # rubi.append(1765)
        return Dist((a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p), Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
def replacement1766(c, m, x, d, f, a, p, e):
        # rubi.append(1766)
        return -Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x), x)
def replacement1767(c, m, x, d, f, a, p, e):
        # rubi.append(1767)
        return Int(ExpandIntegrand((e + f*x)**m*(a + c*x**S(2) + d*x**S(3))**p, x), x)

def With1768(c, m, x, d, f, a, p, e):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + c*x**S(2) + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1768(c, m, x, d, f, a, p, e):

        u = Factor(a + c*x**S(2) + d*x**S(3))
        # rubi.append(1768)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)

def With1769(c, m, x, d, f, a, p, e):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        # rubi.append(1769)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
def replacement1770(c, m, x, d, f, a, p, e):
        # rubi.append(1770)
        return Dist((c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p, Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x), x)

def With1771(c, m, x, d, f, a, p, e):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1771(c, m, x, d, f, a, p, e):

        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        # rubi.append(1771)
        return Dist((a + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)

def With1772(c, m, x, d, f, a, p, e):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        # rubi.append(1772)
        return Dist((a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p), Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
def replacement1773(c, x, d, a, p, b):
        # rubi.append(1773)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + c*x)**(S(3)*p), x), x)
def replacement1774(c, x, d, a, p, b):
        # rubi.append(1774)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Subst(Int((S(3)*a*b*c - b**S(3) + c**S(3)*x**S(3))**p, x), x, c/(S(3)*d) + x), x)

def With1775(c, x, d, a, p, b):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        # rubi.append(1775)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p, x), x)
def replacement1776(c, x, d, a, p, b):
        # rubi.append(1776)
        return Int(ExpandToSum((a + b*x + c*x**S(2) + d*x**S(3))**p, x), x)

def With1777(c, x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1777(c, x, d, a, p, b):

        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        # rubi.append(1777)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)
def replacement1778(c, x, d, a, p, b):
        # rubi.append(1778)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x), x)
def replacement1779(c, x, d, a, p, b):
        # rubi.append(1779)
        return Dist((b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x)**(S(3)*p), x), x)

def With1780(c, x, d, a, p, b):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        # rubi.append(1780)
        return Dist((b + c*x - r)**(-p)*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**(-p)*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)

def With1781(c, x, d, a, p, b):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        # rubi.append(1781)
        return Dist((b + x*(c - r))**(-p)*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**(-p)*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p, x), x)

def With1782(c, x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1782(c, x, d, a, p, b):

        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        # rubi.append(1782)
        return Dist((a + b*x + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)

def With1783(c, x, d, a, p, b):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        # rubi.append(1783)
        return Dist((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**p, x), x)
def replacement1784(x, u, p):
        # rubi.append(1784)
        return Int(ExpandToSum(u, x)**p, x)
def replacement1785(c, m, x, d, f, a, p, e, b):
        # rubi.append(1785)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x), x)

def With1786(c, m, x, d, f, a, p, e, b):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        # rubi.append(1786)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)

def With1787(c, m, x, d, f, a, p, e, b):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        # rubi.append(1787)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p*(e + f*x)**m, x), x)
def replacement1788(c, m, x, d, f, a, p, e, b):
        # rubi.append(1788)
        return Int(ExpandIntegrand((e + f*x)**m*(a + b*x + c*x**S(2) + d*x**S(3))**p, x), x)

def With1789(c, m, x, d, f, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        if ProductQ(NonfreeFactors(u, x)):
            return True
        return False
def replacement1789(c, m, x, d, f, a, p, e, b):

        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        # rubi.append(1789)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)
def replacement1790(c, m, x, d, f, a, p, e, b):
        # rubi.append(1790)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x), x)
def replacement1791(c, m, x, d, f, a, p, e, b):
        # rubi.append(1791)
        return Dist((b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x), x)

def With1792(c, m, x, d, f, a, p, e, b):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        # rubi.append(1792)
        return Dist((b + c*x - r)**(-p)*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**(-p)*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)

def With1793(c, m, x, d, f, a, p, e, b):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        # rubi.append(1793)
        return Dist((b + x*(c - r))**(-p)*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**(-p)*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p*(e + f*x)**m, x), x)

def With1794(c, m, x, d, f, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        if ProductQ(u):
            return True
        return False
def replacement1794(c, m, x, d, f, a, p, e, b):

        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        # rubi.append(1794)
        return Dist((a + b*x + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)

def With1795(c, m, x, d, f, a, p, e, b):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        # rubi.append(1795)
        return Dist((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**p, x), x)
def replacement1796(m, x, v, u, p):
        # rubi.append(1796)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x)
def replacement1797(c, x, d, f, g, a, e, b):
        # rubi.append(1797)
        return Simp(a*f*ArcTan((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(a**S(2)*(S(2)*a - c), S(2))), x)
def replacement1798(c, x, d, f, g, a, e, b):
        # rubi.append(1798)
        return -Simp(a*f*atanh((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(-a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(-a**S(2)*(S(2)*a - c), S(2))), x)
def replacement1799(c, x, d, a, p, e, b):
        # rubi.append(1799)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x)

def With1800(x, p, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        if And(ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3)), NonzeroQ(d)):
            return True
        return False
def replacement1800(x, p, v):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        # rubi.append(1800)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x)
def replacement1801(c, x, d, u, a, p, e, b):
        # rubi.append(1801)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x)

def With1802(x, u, p, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        if And(ZeroQ(S(8)*b*e**S(2) - S(4)*c*d*e + d**S(3)), NonzeroQ(d)):
            return True
        return False
def replacement1802(x, u, p, v):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        # rubi.append(1802)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x)
def replacement1803(c, x, d, a, p, e, b):
        # rubi.append(1803)
        return Dist(-S(16)*a**S(2), Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, S(1)/x + b/(S(4)*a)), x)

def With1804(x, p, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        if And(NonzeroQ(a), NonzeroQ(b), ZeroQ(S(8)*a**S(2)*d - S(4)*a*b*c + b**S(3))):
            return True
        return False
def replacement1804(x, p, v):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        # rubi.append(1804)
        return Dist(-S(16)*a**S(2), Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, S(1)/x + b/(S(4)*a)), x)

def With1805(c, x, d, B, A, D, a, C, e, b):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        # rubi.append(1805)
        return -Dist(S(1)/q, Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)

def With1806(c, x, d, B, A, D, a, e, b):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        # rubi.append(1806)
        return -Dist(S(1)/q, Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)

def With1807(c, m, x, d, B, A, D, a, C, e, b):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        # rubi.append(1807)
        return -Dist(S(1)/q, Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)

def With1808(c, m, x, d, B, A, D, a, e, b):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        # rubi.append(1808)
        return -Dist(S(1)/q, Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)

def With1809(c, x, d, B, A, a, C, e, b):
        q = Rt(C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        # rubi.append(1809)
        return Simp(-S(2)*C**S(2)*atanh((-B*e + C*d + S(2)*C*e*x)/q)/q, x) + Simp(S(2)*C**S(2)*atanh(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q, x)

def With1810(c, x, d, A, a, C, e, b):
        q = Rt(C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        # rubi.append(1810)
        return Simp(-S(2)*C**S(2)*atanh(C*(d + S(2)*e*x)/q)/q, x) + Simp(S(2)*C**S(2)*atanh(C*(A*d - S(2)*C*d*x**S(2) - S(2)*C*e*x**S(3) - S(2)*x*(A*e + C*c))/(A*q))/q, x)

def With1811(c, x, d, B, A, a, C, e, b):
        q = Rt(-C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        # rubi.append(1811)
        return Simp(S(2)*C**S(2)*ArcTan((-B*e + C*d + S(2)*C*e*x)/q)/q, x) - Simp(S(2)*C**S(2)*ArcTan(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q, x)

def With1812(c, x, d, A, a, C, e, b):
        q = Rt(-C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        # rubi.append(1812)
        return Simp(S(2)*C**S(2)*ArcTan((C*d + S(2)*C*e*x)/q)/q, x) - Simp(S(2)*C**S(2)*ArcTan(-C*(-A*d + S(2)*C*d*x**S(2) + S(2)*C*e*x**S(3) + S(2)*x*(A*e + C*c))/(A*q))/q, x)
def replacement1813(c, x, d, B, A, D, a, C, e, b):
        # rubi.append(1813)
        return -Dist(S(1)/(S(4)*e), Int((-S(4)*A*e + D*b + x**S(2)*(-S(4)*C*e + S(3)*D*d) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x), x) + Simp(D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e), x)
def replacement1814(c, x, d, B, A, D, a, e, b):
        # rubi.append(1814)
        return -Dist(S(1)/(S(4)*e), Int((-S(4)*A*e + D*b + S(3)*D*d*x**S(2) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x), x) + Simp(D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e), x)
def replacement1815(c, x, d, u, f, a, e, b):
        # rubi.append(1815)
        return -Dist(a/(f*(-a*d + b*c)), Int(u*sqrt(c + d*x)/x, x), x) + Dist(c/(e*(-a*d + b*c)), Int(u*sqrt(a + b*x)/x, x), x)
def replacement1816(c, x, d, u, f, a, e, b):
        # rubi.append(1816)
        return Dist(b/(f*(-a*d + b*c)), Int(u*sqrt(c + d*x), x), x) - Dist(d/(e*(-a*d + b*c)), Int(u*sqrt(a + b*x), x), x)
def replacement1817(c, x, d, u, f, a, e, b):
        # rubi.append(1817)
        return Dist(e, Int(u*sqrt(a + b*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x), x) - Dist(f, Int(u*sqrt(c + d*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x), x)
def replacement1818(c, n, x, d, u, a, p, b):
        # rubi.append(1818)
        return Dist(S(1)/(a*c), Int(u*sqrt(a + b*x**(S(2)*n)), x), x) - Dist(b/(a*d), Int(u*x**n, x), x)
def replacement1819(c, m, n, x, d, a, p, b):
        # rubi.append(1819)
        return Dist(c, Int(x**m*sqrt(a + b*x**(S(2)*n))/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x), x) - Dist(d, Int(x**(m + n)/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x), x)

def With1820(x, d, f, a, e, b):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        # rubi.append(1820)
        return Dist(r/(S(3)*a), Int(S(1)/((r + s*x)*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(r/(S(3)*a), Int((S(2)*r - s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x), x)

def With1821(x, d, f, a, b):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        # rubi.append(1821)
        return Dist(r/(S(3)*a), Int(S(1)/(sqrt(d + f*x**S(2))*(r + s*x)), x), x) + Dist(r/(S(3)*a), Int((S(2)*r - s*x)/(sqrt(d + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x), x)

def With1822(x, d, f, a, e, b):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        # rubi.append(1822)
        return Dist(r/(S(3)*a), Int(S(1)/((r - s*x)*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(r/(S(3)*a), Int((S(2)*r + s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x), x)

def With1823(x, d, f, a, b):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        # rubi.append(1823)
        return Dist(r/(S(3)*a), Int(S(1)/(sqrt(d + f*x**S(2))*(r - s*x)), x), x) + Dist(r/(S(3)*a), Int((S(2)*r + s*x)/(sqrt(d + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x), x)
def replacement1824(c, x, d, a, e, b):
        # rubi.append(1824)
        return Dist(d, Int(S(1)/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) - Dist(e, Int(x/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
def replacement1825(c, x, d, a, e):
        # rubi.append(1825)
        return Dist(d, Int(S(1)/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x), x) - Dist(e, Int(x/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x), x)
def replacement1826(c, x, d, a, e, b):
        # rubi.append(1826)
        return -Dist(c/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Simp(e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), x)
def replacement1827(c, x, d, a, e, b):
        # rubi.append(1827)
        return -Dist(c/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((b*d*e**S(2) + S(2)*c*d**S(3))/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int(S(1)/((d + e*x)*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) - Simp(e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), x)
def replacement1828(c, x, d, a, e):
        # rubi.append(1828)
        return -Dist(c/(a*e**S(4) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist(S(2)*c*d**S(3)/(a*e**S(4) + c*d**S(4)), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x)), x), x) - Simp(e**S(3)*sqrt(a + c*x**S(4))/((d + e*x)*(a*e**S(4) + c*d**S(4))), x)
def replacement1829(c, x, d, B, A, a, e, b):
        # rubi.append(1829)
        return Dist(A, Subst(Int(S(1)/(d - x**S(2)*(-S(2)*a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1830(c, x, d, B, A, a, e):
        # rubi.append(1830)
        return Dist(A, Subst(Int(S(1)/(S(2)*a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))), x)
def replacement1831(c, x, d, B, A, f, a, e, b):
        # rubi.append(1831)
        return Dist(A, Subst(Int(S(1)/(d - x**S(2)*(-a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1832(c, x, d, B, A, f, a, e):
        # rubi.append(1832)
        return Dist(A, Subst(Int(S(1)/(a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))), x)
def replacement1833(c, x, d, B, A, f, a, b):
        # rubi.append(1833)
        return Dist(A, Subst(Int(S(1)/(-b*d*x**S(2) + d), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
def replacement1834(c, x, d, a, e, b):
        # rubi.append(1834)
        return Dist(a/d, Subst(Int(S(1)/(-S(2)*b*x**S(2) + x**S(4)*(-S(4)*a*c + b**S(2)) + S(1)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)

def With1835(c, x, d, a, e, b):
        q = sqrt(-S(4)*a*c + b**S(2))
        # rubi.append(1835)
        return Simp(sqrt(S(2))*a*sqrt(-b + q)*atanh(sqrt(S(2))*x*sqrt(-b + q)*(b + S(2)*c*x**S(2) + q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2))), x) - Simp(sqrt(S(2))*a*sqrt(b + q)*ArcTan(sqrt(S(2))*x*sqrt(b + q)*(b + S(2)*c*x**S(2) - q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2))), x)
def replacement1836(c, x, d, f, a, e, b):
        # rubi.append(1836)
        return Dist(a, Int(S(1)/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) - Dist(b, Int(x/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)
def replacement1837(c, x, d, h, f, g, a, e, b):
        # rubi.append(1837)
        return Simp(S(2)*sqrt(d + e*x + f*sqrt(a + b*x + c*x**S(2)))*(S(9)*c**S(2)*f*g*h*x**S(2) + S(3)*c**S(2)*f*h**S(2)*x**S(3) + c*f*x*(a*h**S(2) - b*g*h + S(10)*c*g**S(2)) + f*(S(2)*a*b*h**S(2) - S(3)*a*c*g*h - S(2)*b**S(2)*g*h + S(5)*b*c*g**S(2)) - (-d*h + e*g)*sqrt(a + b*x + c*x**S(2))*(-S(2)*b*h + S(5)*c*g + c*h*x))/(S(15)*c**S(2)*f*(g + h*x)), x)
def replacement1838(m, j, n, x, v, k, h, u, f, g):
        # rubi.append(1838)
        return Int((g + h*x)**m*(f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n, x)
def replacement1839(c, n, x, d, h, f, g, a, p, e, b):
        # rubi.append(1839)
        return Dist(S(2), Subst(Int((g + h*x**n)**p*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))/(b*f**S(2) - S(2)*d*e + S(2)*e*x)**S(2), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)
def replacement1840(c, n, x, d, h, f, g, a, p, e):
        # rubi.append(1840)
        return Dist(S(1)/(S(2)*e), Subst(Int((g + h*x**n)**p*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))/(d - x)**S(2), x), x, d + e*x + f*sqrt(a + c*x**S(2))), x)
def replacement1841(n, x, v, h, u, f, g, p):
        # rubi.append(1841)
        return Int((g + h*(f*sqrt(ExpandToSum(v, x)) + ExpandToSum(u, x))**n)**p, x)
def replacement1842(c, m, n, x, h, f, g, a, e):
        # rubi.append(1842)
        return Dist(S(2)**(-m + S(-1))*e**(-m + S(-1)), Subst(Int(x**(-m + n + S(-2))*(a*f**S(2) + x**S(2))*(-a*f**S(2)*h + S(2)*e*g*x + h*x**S(2))**m, x), x, e*x + f*sqrt(a + c*x**S(2))), x)
def replacement1843(c, m, n, x, f, i, g, a, p, e):
        # rubi.append(1843)
        return Dist(S(2)**(-S(2)*m - p + S(-1))*e**(-p + S(-1))*f**(-S(2)*m)*(i/c)**m, Subst(Int(x**(-S(2)*m + n - p + S(-2))*(-a*f**S(2) + x**S(2))**p*(a*f**S(2) + x**S(2))**(S(2)*m + S(1)), x), x, e*x + f*sqrt(a + c*x**S(2))), x)
def replacement1844(c, m, n, x, d, h, f, i, g, a, e, b):
        # rubi.append(1844)
        return Dist(S(2)*f**(-S(2)*m)*(i/c)**m, Subst(Int(x**n*(b*f**S(2) - S(2)*d*e + S(2)*e*x)**(-S(2)*m + S(-2))*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)
def replacement1845(c, m, n, x, d, f, i, g, a, e):
        # rubi.append(1845)
        return Dist(S(2)**(-S(2)*m + S(-1))*f**(-S(2)*m)*(i/c)**m/e, Subst(Int(x**n*(-d + x)**(-S(2)*m + S(-2))*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + c*x**S(2))), x)
def replacement1846(c, m, n, x, d, h, f, i, g, a, e, b):
        # rubi.append(1846)
        return Dist((i/c)**(m + S(-1)/2)*sqrt(g + h*x + i*x**S(2))/sqrt(a + b*x + c*x**S(2)), Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x), x)
def replacement1847(c, m, n, x, d, f, i, g, a, e):
        # rubi.append(1847)
        return Dist((i/c)**(m + S(-1)/2)*sqrt(g + i*x**S(2))/sqrt(a + c*x**S(2)), Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x), x)
def replacement1848(c, m, n, x, d, h, f, i, g, a, e, b):
        # rubi.append(1848)
        return Dist((i/c)**(m + S(1)/2)*sqrt(a + b*x + c*x**S(2))/sqrt(g + h*x + i*x**S(2)), Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x), x)
def replacement1849(c, m, n, x, d, f, i, g, a, e):
        # rubi.append(1849)
        return Dist((i/c)**(m + S(1)/2)*sqrt(a + c*x**S(2))/sqrt(g + i*x**S(2)), Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x), x)
def replacement1850(m, j, n, x, v, k, u, w, f):
        # rubi.append(1850)
        return Int((f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n*ExpandToSum(w, x)**m, x)
def replacement1851(c, n, x, d, a, p, b):
        # rubi.append(1851)
        return Dist(S(1)/a, Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*(a + b*x**n)**(S(2)/n))), x)
def replacement1852(c, x, d, a, b):
        # rubi.append(1852)
        return Simp(S(2)*a*x/sqrt(a + b*sqrt(c + d*x**S(2))), x) + Simp(S(2)*b**S(2)*d*x**S(3)/(S(3)*(a + b*sqrt(c + d*x**S(2)))**(S(3)/2)), x)
def replacement1853(c, x, d, a, b):
        # rubi.append(1853)
        return Dist(sqrt(S(2))*b/a, Subst(Int(S(1)/sqrt(S(1) + x**S(2)/a), x), x, a*x + b*sqrt(c + d*x**S(2))), x)
def replacement1854(c, x, d, a, e, b):
        # rubi.append(1854)
        return Int(sqrt(a*e*x**S(2) + b*e*x*sqrt(c + d*x**S(2)))/(x*sqrt(c + d*x**S(2))), x)
def replacement1855(c, x, d, a, b):
        # rubi.append(1855)
        return Dist(d, Subst(Int(S(1)/(-S(2)*c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*sqrt(a + b*x**S(4)))), x)
def replacement1856(c, m, x, d, a, e, b):
        # rubi.append(1856)
        return Dist(S(1)/2 - I/S(2), Int((c + d*x)**m/sqrt(sqrt(a) - I*b*x**S(2)), x), x) + Dist(S(1)/2 + I/S(2), Int((c + d*x)**m/sqrt(sqrt(a) + I*b*x**S(2)), x), x)

def With1857(c, x, d, a, b):
        q = Rt(b/a, S(3))
        # rubi.append(1857)
        return Dist(d/(-c*q + d*(S(1) + sqrt(S(3)))), Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) - Dist(q/(-c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1858(c, x, d, a, b):
        q = Rt(-b/a, S(3))
        # rubi.append(1858)
        return Dist(d/(c*q + d*(S(1) + sqrt(S(3)))), Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist(q/(c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1859(c, x, d, a, b):
        q = Rt(-b/a, S(3))
        # rubi.append(1859)
        return Dist(d/(c*q + d*(S(1) - sqrt(S(3)))), Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist(q/(c*q + d*(S(1) - sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1860(c, x, d, a, b):
        q = Rt(b/a, S(3))
        # rubi.append(1860)
        return Dist(d/(-c*q + d*(S(1) - sqrt(S(3)))), Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) - Dist(q/(-c*q + d*(S(1) - sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1861(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(b/a, S(3))
        if ZeroQ(-e*q + f*(S(1) + sqrt(S(3)))):
            return True
        return False
def replacement1861(c, x, d, f, a, e, b):

        q = Rt(b/a, S(3))
        # rubi.append(1861)
        return Dist(S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(S(2) - sqrt(S(3)))*(q*x + S(1))/(q*sqrt((q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(S(1) - sqrt(S(3))) + x*(-c*q + d*(S(1) + sqrt(S(3)))))), x), x, (-q*x + S(-1) + sqrt(S(3)))/(q*x + S(1) + sqrt(S(3)))), x)

def With1862(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(b/a, S(3))
        if NonzeroQ(-e*q + f*(S(1) + sqrt(S(3)))):
            return True
        return False
def replacement1862(c, x, d, f, a, e, b):

        q = Rt(b/a, S(3))
        # rubi.append(1862)
        return Dist((-c*f + d*e)/(-c*q + d*(S(1) + sqrt(S(3)))), Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((-e*q + f*(S(1) + sqrt(S(3))))/(-c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1863(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(-b/a, S(3))
        if ZeroQ(e*q + f*(S(1) + sqrt(S(3)))):
            return True
        return False
def replacement1863(c, x, d, f, a, e, b):

        q = Rt(-b/a, S(3))
        # rubi.append(1863)
        return Dist(-S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(S(2) - sqrt(S(3)))*(-q*x + S(1))/(q*sqrt((-q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(c*q + d*(S(1) - sqrt(S(3))) + x*(c*q + d*(S(1) + sqrt(S(3)))))), x), x, (q*x + S(-1) + sqrt(S(3)))/(-q*x + S(1) + sqrt(S(3)))), x)

def With1864(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(-b/a, S(3))
        if NonzeroQ(e*q + f*(S(1) + sqrt(S(3)))):
            return True
        return False
def replacement1864(c, x, d, f, a, e, b):

        q = Rt(-b/a, S(3))
        # rubi.append(1864)
        return Dist((-c*f + d*e)/(c*q + d*(S(1) + sqrt(S(3)))), Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((e*q + f*(S(1) + sqrt(S(3))))/(c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1865(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(-b/a, S(3))
        if ZeroQ(e*q + f*(S(1) - sqrt(S(3)))):
            return True
        return False
def replacement1865(c, x, d, f, a, e, b):

        q = Rt(-b/a, S(3))
        # rubi.append(1865)
        return Dist(S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(-q*x + S(1))/(q*sqrt(-(-q*x + S(1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(c*q + d*(S(1) + sqrt(S(3))) + x*(c*q + d*(S(1) - sqrt(S(3)))))), x), x, (-q*x + S(1) + sqrt(S(3)))/(q*x + S(-1) + sqrt(S(3)))), x)

def With1866(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(-b/a, S(3))
        if NonzeroQ(e*q + f*(S(1) - sqrt(S(3)))):
            return True
        return False
def replacement1866(c, x, d, f, a, e, b):

        q = Rt(-b/a, S(3))
        # rubi.append(1866)
        return Dist((-c*f + d*e)/(c*q + d*(S(1) - sqrt(S(3)))), Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((e*q + f*(S(1) - sqrt(S(3))))/(c*q + d*(S(1) - sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)

def With1867(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(b/a, S(3))
        if ZeroQ(-e*q + f*(S(1) - sqrt(S(3)))):
            return True
        return False
def replacement1867(c, x, d, f, a, e, b):

        q = Rt(b/a, S(3))
        # rubi.append(1867)
        return Dist(-S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(q*x + S(1))/(q*sqrt(-(q*x + S(1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(S(1) - x**S(2))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(S(1) + sqrt(S(3))) + x*(-c*q + d*(S(1) - sqrt(S(3)))))), x), x, (q*x + S(1) + sqrt(S(3)))/(-q*x + S(-1) + sqrt(S(3)))), x)

def With1868(c, x, d, f, a, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Rt(b/a, S(3))
        if NonzeroQ(-e*q + f*(S(1) - sqrt(S(3)))):
            return True
        return False
def replacement1868(c, x, d, f, a, e, b):

        q = Rt(b/a, S(3))
        # rubi.append(1868)
        return Dist((-c*f + d*e)/(-c*q + d*(S(1) - sqrt(S(3)))), Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((-e*q + f*(S(1) - sqrt(S(3))))/(-c*q + d*(S(1) - sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
def replacement1869(c, m, n, x, d, a, e, b):
        # rubi.append(1869)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)/(c + d*x + e*sqrt(a + b*x)), x), x, x**n), x)
def replacement1870(c, n, x, d, u, a, e, b):
        # rubi.append(1870)
        return Dist(c, Int(u/(-a*e**S(2) + c**S(2) + c*d*x**n), x), x) - Dist(a*e, Int(u/(sqrt(a + b*x**n)*(-a*e**S(2) + c**S(2) + c*d*x**n)), x), x)
def replacement1871(c, n, x, d, B, A, n2, a, b):
        # rubi.append(1871)
        return Dist(A**S(2)*(n + S(-1)), Subst(Int(S(1)/(A**S(2)*b*x**S(2)*(n + S(-1))**S(2) + a), x), x, x/(A*(n + S(-1)) - B*x**n)), x)
def replacement1872(c, m, n, x, d, k, B, A, n2, a, b):
        # rubi.append(1872)
        return Dist(A**S(2)*(m - n + S(1))/(m + S(1)), Subst(Int(S(1)/(A**S(2)*b*x**S(2)*(m - n + S(1))**S(2) + a), x), x, x**(m + S(1))/(A*(m - n + S(1)) + B*x**n*(m + S(1)))), x)
def replacement1873(c, n, x, d, n3, f, g, n2, a, p, e, b):
        # rubi.append(1873)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1874(c, n, x, d, f, n2, a, p, e, b):
        # rubi.append(1874)
        return -Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*e - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*d*(n*p + n + S(1)) - x**n*(-S(2)*a*c*e*(n*(S(2)*p + S(3)) + S(1)) + b*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d + x**n*(-S(2)*a*c*e + b*(a*f + c*d)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1875(c, n, x, d, n3, g, n2, a, p, e, b):
        # rubi.append(1875)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1876(c, n, x, d, n3, f, g, n2, a, p, b):
        # rubi.append(1876)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1877(c, n, x, d, f, n2, a, p, b):
        # rubi.append(1877)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) + b**S(2)*d*(n*p + n + S(1)) + b*x**n*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*(-a*f + c*d) + b**S(2)*d + b*x**n*(a*f + c*d))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1878(c, n, x, d, n3, g, n2, a, p, b):
        # rubi.append(1878)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
def replacement1879(c, n, x, d, n3, f, g, n2, a, p, e):
        # rubi.append(1879)
        return -Dist(-S(1)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*x**n*(-a*g + c*e) - S(2)*a*c*(-a*f + c*d))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), x)
def replacement1880(c, n, x, d, f, n2, a, p, e):
        # rubi.append(1880)
        return -Dist(-S(1)/(S(4)*a**S(2)*c*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c*e*x**n*(n*(S(2)*p + S(3)) + S(1)) - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*e*x**n - S(2)*a*(-a*f + c*d))/(S(4)*a**S(2)*c*n*(p + S(1))), x)
def replacement1881(c, n, x, d, n3, g, n2, a, p, e):
        # rubi.append(1881)
        return -Dist(-S(1)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c**S(2)*d - S(2)*a*c*x**n*(-a*g + c*e))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), x)

def With1882(c, x, d, f, g, a, e, b):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + f*(-S(2)*a*b*g + S(3)*c**S(2)*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((a*c*f**S(2) - f*(S(2)*a*b*g + S(3)*c**S(2)*d) + S(4)*g*(a**S(2)*g + b*c*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        # rubi.append(1882)
        return -Simp(c*ArcTan((r - S(2)*x)/q)/(g*q), x) + Simp(c*ArcTan((r + S(2)*x)/q)/(g*q), x) - Simp(c*ArcTan(x*(-a*f + S(3)*c*d)*(S(6)*a**S(2)*b*g**S(2) - S(2)*a**S(2)*c*f*g - a*b**S(2)*f*g + b*c**S(2)*d*f + c**S(2)*g*x**S(4)*(-a*f + S(3)*c*d) + c*x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) - b*c*d*g + S(3)*c**S(2)*d*f))/(g*q*(-S(2)*a**S(2)*g + b*c*d)*(S(4)*a**S(2)*g - a*b*f + b*c*d)))/(g*q), x)

def With1883(c, x, d, f, g, a, e):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((S(4)*a**S(2)*g**S(2) + a*c*f**S(2) - S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        # rubi.append(1883)
        return -Simp(c*ArcTan((r - S(2)*x)/q)/(g*q), x) + Simp(c*ArcTan((r + S(2)*x)/q)/(g*q), x) - Simp(c*ArcTan(c*x*(-a*f + S(3)*c*d)*(S(2)*a**S(2)*f*g - c*g*x**S(4)*(-a*f + S(3)*c*d) - x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f))/(S(8)*a**S(4)*g**S(3)*q))/(g*q), x)

def With1884(x, u, p, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            m = Exponent(u, x)
            n = Exponent(v, x)
            c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
            c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
            w = Apart(-c*x**(m - n)*(v*(m - n + S(1)) + x*(p + S(1))*D(v, x)) + u, x)
            res = And(Inequality(S(1), Less, n, LessEqual, m + S(1)), Less(m + n*p, S(-1)), FalseQ(DerivativeDivides(v, u, x)))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement1884(x, u, p, v):

        m = Exponent(u, x)
        n = Exponent(v, x)
        c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
        c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
        w = Apart(-c*x**(m - n)*(v*(m - n + S(1)) + x*(p + S(1))*D(v, x)) + u, x)
        # rubi.append(1884)
        return Simp(If(ZeroQ(w), c*v**(p + 1)*x**(m - n + 1), c*v**(p + 1)*x**(m - n + 1) + Int(v**p*w, x)), x)
