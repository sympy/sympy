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

def miscellaneous_algebraic(rubi, matcher, load_rule):
    from sympy.integrals.rubi.constraints import cons798, cons2, cons3, cons7, cons50, cons4, cons5, cons17, cons21, cons799, cons27, cons48, cons125, cons52, cons800, cons25, cons801, cons802, cons149, cons803, cons500, cons804, cons648, cons805, cons806, cons18, cons46, cons807, cons808, cons68, cons809, cons810, cons811, cons812, cons813, cons814, cons815, cons816, cons817, cons818, cons819, cons820, cons821, cons452, cons822, cons823, cons824, cons825, cons826, cons827, cons828, cons829, cons830, cons831, cons832, cons833, cons834, cons835, cons836, cons837, cons838, cons839, cons840, cons841, cons842, cons843, cons844, cons845, cons846, cons847, cons848, cons849, cons850, cons851, cons852, cons208, cons209, cons64, cons853, cons66, cons854, cons855, cons464, cons856, cons857, cons858, cons53, cons13, cons137, cons859, cons860, cons148, cons244, cons163, cons861, cons521, cons_with_1532, cons862, cons863, cons864, cons_with_1541, cons84, cons_with_1543, cons865, cons34, cons35, cons866, cons468, cons469, cons867, cons868, cons36, cons869, cons870, cons871, cons872, cons873, cons874, cons875, cons876, cons877, cons878, cons879, cons880, cons881, cons882, cons883, cons884, cons885, cons886, cons887, cons888, cons889, cons890, cons891, cons892, cons893, cons894, cons895, cons896, cons897, cons898, cons899, cons900, cons901, cons902, cons903, cons_with_1581, cons_with_1582, cons_with_1583, cons904, cons_with_1584, cons_with_1585, cons_with_1586, cons674, cons905, cons_with_1587, cons_with_1588, cons481, cons906, cons907, cons482, cons908, cons909, cons910, cons911, cons912, cons913, cons914, cons915, cons916, cons85, cons31, cons94, cons917, cons_with_1603, cons_with_1604, cons_with_1605, cons196, cons367, cons356, cons489, cons541, cons23, cons918, cons554, cons919, cons552, cons55, cons494, cons57, cons58, cons59, cons60, cons920, cons921, cons922, cons923, cons924, cons595, cons71, cons925, cons586, cons87, cons128, cons926, cons927, cons928, cons929, cons930, cons45, cons314, cons226, cons931, cons932, cons933, cons934, cons935, cons936, cons937, cons938, cons939, cons940, cons941, cons_with_1648, cons942, cons943, cons944, cons945, cons946, cons_with_1656, cons_with_1657, cons_with_1658, cons282, cons947, cons63, cons_with_1661, cons719, cons948, cons_with_1662, cons949, cons_with_1663, cons_with_1664, cons_with_1665, cons950, cons73, cons951, cons702, cons147, cons952, cons953, cons796, cons954, cons955, cons956, cons957, cons958, cons959, cons960, cons961, cons962, cons963, cons964, cons965, cons966, cons69, cons967, cons968, cons969, cons970, cons971, cons972, cons973, cons974, cons975, cons512, cons976, cons977, cons978, cons979, cons980, cons667, cons981, cons982, cons797, cons983, cons984, cons985, cons986, cons987, cons988, cons93, cons88, cons989, cons990, cons991, cons992, cons993, cons994, cons995, cons996, cons_with_1735, cons997, cons998, cons_with_1736, cons38, cons999, cons1000, cons_with_1744, cons_with_1747, cons_with_1751, cons_with_1754, cons1001, cons1002, cons_with_1758, cons_with_1761, cons_with_1765, cons_with_1768, cons1003, cons1004, cons1005, cons1006, cons_with_1774, cons_with_1779, cons1007, cons1008, cons_with_1786, cons_with_1791, cons1009, cons1010, cons383, cons1011, cons1012, cons1013, cons1014, cons1015, cons1016, cons1017, cons1018, cons_with_1797, cons357, cons_with_1799, cons1019, cons246, cons_with_1801, cons1020, cons1021, cons1022, cons1023, cons1024, cons1025, cons1026, cons1027, cons1028, cons1029, cons1030, cons1031, cons1032, cons1033, cons1034, cons1035, cons1036, cons1037, cons1038, cons1039, cons1040, cons1041, cons1042, cons1043, cons297, cons1044, cons1045, cons1046, cons1047, cons1048, cons705, cons382, cons1049, cons1050, cons697, cons709, cons153, cons1051, cons1052, cons1053, cons1054, cons1055, cons1056, cons1057, cons1058, cons1059, cons224, cons1060, cons515, cons1061, cons1062, cons1063, cons1064, cons1065, cons1066, cons1067, cons1068, cons1069, cons1070, cons1071, cons43, cons479, cons480, cons_with_1858, cons_with_1859, cons_with_1860, cons_with_1861, cons_with_1862, cons_with_1863, cons_with_1864, cons_with_1865, cons1072, cons1073, cons1074, cons1075, cons1076, cons1077, cons1078, cons1079, cons1080, cons1081, cons1082, cons1083, cons1084, cons1085, cons1086, cons1087, cons1088, cons1089, cons_with_1881

    pattern1473 = Pattern(Integral(((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons50, cons4, cons5, cons798)
    def replacement1473(q, p, c, x, b, a, n):
        rubi.append(1473)
        return Dist(x*(c*x**n)**(-S(1)/n), Subst(Int((a + b*x**(n*q))**p, x), x, (c*x**n)**(S(1)/n)), x)
    rule1473 = ReplacementRule(pattern1473, replacement1473)

    if load_rule:
        matcher.add(pattern1473, 1473)
    pattern1474 = Pattern(Integral(x_**WC('m', S(1))*((x_**n_*WC('c', S(1)))**q_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons50, cons798, cons17)
    def replacement1474(q, p, c, m, x, b, a, n):
        rubi.append(1474)
        return Dist(x**(m + S(1))*(c*x**n)**(-(m + S(1))/n), Subst(Int(x**m*(a + b*x**(n*q))**p, x), x, (c*x**n)**(S(1)/n)), x)
    rule1474 = ReplacementRule(pattern1474, replacement1474)

    if load_rule:
        matcher.add(pattern1474, 1474)
    pattern1475 = Pattern(Integral(x_**WC('m', S(1))*((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('r', S(1))*WC('e', S(1)))**p_*((c_ + x_**WC('n', S(1))*WC('d', S(1)))**s_*WC('f', S(1)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons52, cons800, cons799)
    def replacement1475(d, s, q, p, c, f, e, x, r, b, a, m, n):
        rubi.append(1475)
        return Dist((e*(a + b*x**n)**r)**p*(f*(c + d*x**n)**s)**q*(a + b*x**n)**(-p*r)*(c + d*x**n)**(-q*s), Int(x**m*(a + b*x**n)**(p*r)*(c + d*x**n)**(q*s), x), x)
    rule1475 = ReplacementRule(pattern1475, replacement1475)

    if load_rule:
        matcher.add(pattern1475, 1475)
    pattern1476 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons25)
    def replacement1476(d, p, c, e, u, x, b, a, n):
        rubi.append(1476)
        return Dist((b*e/d)**p, Int(u, x), x)
    rule1476 = ReplacementRule(pattern1476, replacement1476)

    if load_rule:
        matcher.add(pattern1476, 1476)
    pattern1477 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons801, cons802)
    def replacement1477(d, p, c, e, u, x, b, a, n):
        rubi.append(1477)
        return Int(u*(e*(a + b*x**n))**p*(c + d*x**n)**(-p), x)
    rule1477 = ReplacementRule(pattern1477, replacement1477)

    if load_rule:
        matcher.add(pattern1477, 1477)
    def With1478(d, p, c, e, x, b, a, n):
        q = Denominator(p)
        rubi.append(1478)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + S(1)/n)*(b*e - d*x**q)**(S(-1) - S(1)/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)
    pattern1478 = Pattern(Integral(((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons149, cons803)
    rule1478 = ReplacementRule(pattern1478, With1478)

    if load_rule:
        matcher.add(pattern1478, 1478)
    def With1479(d, p, c, e, m, x, b, a, n):
        q = Denominator(p)
        rubi.append(1479)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)
    pattern1479 = Pattern(Integral(x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons149, cons500)
    rule1479 = ReplacementRule(pattern1479, With1479)

    if load_rule:
        matcher.add(pattern1479, 1479)
    def With1480(d, p, c, e, u, x, r, b, a, n):
        q = Denominator(p)
        rubi.append(1480)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + S(1)/n)*(b*e - d*x**q)**(S(-1) - S(1)/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(S(1)/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)
    pattern1480 = Pattern(Integral(u_**WC('r', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons804, cons149, cons803, cons648)
    rule1480 = ReplacementRule(pattern1480, With1480)

    if load_rule:
        matcher.add(pattern1480, 1480)
    def With1481(d, p, c, e, u, m, x, r, b, a, n):
        q = Denominator(p)
        rubi.append(1481)
        return Dist(e*q*(-a*d + b*c)/n, Subst(Int(SimplifyIntegrand(x**(q*(p + S(1)) + S(-1))*(-a*e + c*x**q)**(S(-1) + (m + S(1))/n)*(b*e - d*x**q)**(S(-1) - (m + S(1))/n)*ReplaceAll(u, Rule(x, (-a*e + c*x**q)**(S(1)/n)*(b*e - d*x**q)**(-S(1)/n)))**r, x), x), x, (e*(a + b*x**n)/(c + d*x**n))**(S(1)/q)), x)
    pattern1481 = Pattern(Integral(u_**WC('r', S(1))*x_**WC('m', S(1))*((x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(c_ + x_**WC('n', S(1))*WC('d', S(1))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons804, cons149, cons803, cons805)
    rule1481 = ReplacementRule(pattern1481, With1481)

    if load_rule:
        matcher.add(pattern1481, 1481)
    pattern1482 = Pattern(Integral(((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons806)
    def replacement1482(p, c, x, b, a, n):
        rubi.append(1482)
        return -Dist(c, Subst(Int((a + b*x**n)**p/x**S(2), x), x, c/x), x)
    rule1482 = ReplacementRule(pattern1482, replacement1482)

    if load_rule:
        matcher.add(pattern1482, 1482)
    pattern1483 = Pattern(Integral(x_**WC('m', S(1))*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons17)
    def replacement1483(p, c, m, x, b, a, n):
        rubi.append(1483)
        return -Dist(c**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x), x)
    rule1483 = ReplacementRule(pattern1483, replacement1483)

    if load_rule:
        matcher.add(pattern1483, 1483)
    pattern1484 = Pattern(Integral((x_*WC('d', S(1)))**m_*((WC('c', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons18)
    def replacement1484(d, p, c, m, x, b, a, n):
        rubi.append(1484)
        return -Dist(c*(c/x)**m*(d*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n)**p, x), x, c/x), x)
    rule1484 = ReplacementRule(pattern1484, replacement1484)

    if load_rule:
        matcher.add(pattern1484, 1484)
    pattern1485 = Pattern(Integral(((WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons46)
    def replacement1485(d, p, n2, x, b, c, a, n):
        rubi.append(1485)
        return -Dist(d, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p/x**S(2), x), x, d/x), x)
    rule1485 = ReplacementRule(pattern1485, replacement1485)

    if load_rule:
        matcher.add(pattern1485, 1485)
    pattern1486 = Pattern(Integral(x_**WC('m', S(1))*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons46, cons17)
    def replacement1486(d, p, n2, x, b, c, a, m, n):
        rubi.append(1486)
        return -Dist(d**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x), x)
    rule1486 = ReplacementRule(pattern1486, replacement1486)

    if load_rule:
        matcher.add(pattern1486, 1486)
    pattern1487 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + (WC('d', S(1))/x_)**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons46, cons18)
    def replacement1487(d, p, e, n2, m, x, b, c, a, n):
        rubi.append(1487)
        return -Dist(d*(d/x)**m*(e*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, d/x), x)
    rule1487 = ReplacementRule(pattern1487, replacement1487)

    if load_rule:
        matcher.add(pattern1487, 1487)
    pattern1488 = Pattern(Integral((x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons807, cons808)
    def replacement1488(d, p, c, n2, x, b, a, n):
        rubi.append(1488)
        return -Dist(d, Subst(Int((a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p/x**S(2), x), x, d/x), x)
    rule1488 = ReplacementRule(pattern1488, replacement1488)

    if load_rule:
        matcher.add(pattern1488, 1488)
    pattern1489 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons807, cons808, cons17)
    def replacement1489(d, p, c, n2, x, b, a, m, n):
        rubi.append(1489)
        return -Dist(d**(m + S(1)), Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x), x)
    rule1489 = ReplacementRule(pattern1489, replacement1489)

    if load_rule:
        matcher.add(pattern1489, 1489)
    pattern1490 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)) + (WC('d', S(1))/x_)**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons807, cons18, cons808)
    def replacement1490(d, p, c, e, n2, x, b, a, m, n):
        rubi.append(1490)
        return -Dist(d*(d/x)**m*(e*x)**m, Subst(Int(x**(-m + S(-2))*(a + b*x**n + c*d**(-S(2)*n)*x**(S(2)*n))**p, x), x, d/x), x)
    rule1490 = ReplacementRule(pattern1490, replacement1490)

    if load_rule:
        matcher.add(pattern1490, 1490)
    pattern1491 = Pattern(Integral(u_**m_, x_), cons21, cons68, cons809)
    def replacement1491(m, u, x):
        rubi.append(1491)
        return Int(ExpandToSum(u, x)**m, x)
    rule1491 = ReplacementRule(pattern1491, replacement1491)

    if load_rule:
        matcher.add(pattern1491, 1491)
    pattern1492 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1)), x_), cons21, cons4, cons810, cons811)
    def replacement1492(v, u, x, m, n):
        rubi.append(1492)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n, x)
    rule1492 = ReplacementRule(pattern1492, replacement1492)

    if load_rule:
        matcher.add(pattern1492, 1492)
    pattern1493 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), cons21, cons4, cons5, cons812, cons813)
    def replacement1493(v, p, u, x, w, m, n):
        rubi.append(1493)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x)
    rule1493 = ReplacementRule(pattern1493, replacement1493)

    if load_rule:
        matcher.add(pattern1493, 1493)
    pattern1494 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1))*z_**WC('q', S(1)), x_), cons21, cons4, cons5, cons50, cons814, cons815)
    def replacement1494(q, p, v, z, u, x, w, m, n):
        rubi.append(1494)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p*ExpandToSum(z, x)**q, x)
    rule1494 = ReplacementRule(pattern1494, replacement1494)

    if load_rule:
        matcher.add(pattern1494, 1494)
    pattern1495 = Pattern(Integral(u_**p_, x_), cons5, cons816, cons817)
    def replacement1495(x, u, p):
        rubi.append(1495)
        return Int(ExpandToSum(u, x)**p, x)
    rule1495 = ReplacementRule(pattern1495, replacement1495)

    if load_rule:
        matcher.add(pattern1495, 1495)
    pattern1496 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), cons21, cons5, cons68, cons818, cons819)
    def replacement1496(v, p, u, x, m):
        rubi.append(1496)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x)
    rule1496 = ReplacementRule(pattern1496, replacement1496)

    if load_rule:
        matcher.add(pattern1496, 1496)
    pattern1497 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('n', S(1))*w_**WC('p', S(1)), x_), cons21, cons4, cons5, cons810, cons820, cons821)
    def replacement1497(v, p, u, x, w, m, n):
        rubi.append(1497)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**n*ExpandToSum(w, x)**p, x)
    rule1497 = ReplacementRule(pattern1497, replacement1497)

    if load_rule:
        matcher.add(pattern1497, 1497)
    pattern1498 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons50, cons452, cons822)
    def replacement1498(v, q, p, u, x):
        rubi.append(1498)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
    rule1498 = ReplacementRule(pattern1498, replacement1498)

    if load_rule:
        matcher.add(pattern1498, 1498)
    pattern1499 = Pattern(Integral(u_**p_, x_), cons5, cons823, cons824)
    def replacement1499(x, u, p):
        rubi.append(1499)
        return Int(ExpandToSum(u, x)**p, x)
    rule1499 = ReplacementRule(pattern1499, replacement1499)

    if load_rule:
        matcher.add(pattern1499, 1499)
    pattern1500 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('c', S(1)))**WC('m', S(1)), x_), cons7, cons21, cons5, cons823, cons824)
    def replacement1500(p, u, x, c, m):
        rubi.append(1500)
        return Int((c*x)**m*ExpandToSum(u, x)**p, x)
    rule1500 = ReplacementRule(pattern1500, replacement1500)

    if load_rule:
        matcher.add(pattern1500, 1500)
    pattern1501 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons5, cons50, cons825, cons826, cons827)
    def replacement1501(v, q, p, u, x):
        rubi.append(1501)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
    rule1501 = ReplacementRule(pattern1501, replacement1501)

    if load_rule:
        matcher.add(pattern1501, 1501)
    pattern1502 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), cons21, cons5, cons50, cons825, cons826, cons827)
    def replacement1502(v, q, p, u, x, m):
        rubi.append(1502)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q, x)
    rule1502 = ReplacementRule(pattern1502, replacement1502)

    if load_rule:
        matcher.add(pattern1502, 1502)
    pattern1503 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1))*w_**WC('q', S(1)), x_), cons21, cons5, cons50, cons828, cons826, cons829, cons830)
    def replacement1503(v, q, p, u, x, w, m):
        rubi.append(1503)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p*ExpandToSum(w, x)**q, x)
    rule1503 = ReplacementRule(pattern1503, replacement1503)

    if load_rule:
        matcher.add(pattern1503, 1503)
    pattern1504 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1))*z_**WC('r', S(1)), x_), cons21, cons5, cons50, cons52, cons831, cons826, cons832, cons833)
    def replacement1504(q, p, v, z, u, x, r, m):
        rubi.append(1504)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**r, x)
    rule1504 = ReplacementRule(pattern1504, replacement1504)

    if load_rule:
        matcher.add(pattern1504, 1504)
    pattern1505 = Pattern(Integral(u_**p_, x_), cons5, cons834, cons835)
    def replacement1505(x, u, p):
        rubi.append(1505)
        return Int(ExpandToSum(u, x)**p, x)
    rule1505 = ReplacementRule(pattern1505, replacement1505)

    if load_rule:
        matcher.add(pattern1505, 1505)
    pattern1506 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), cons21, cons5, cons834, cons835)
    def replacement1506(m, x, p, u):
        rubi.append(1506)
        return Int(x**m*ExpandToSum(u, x)**p, x)
    rule1506 = ReplacementRule(pattern1506, replacement1506)

    if load_rule:
        matcher.add(pattern1506, 1506)
    pattern1507 = Pattern(Integral(u_**p_, x_), cons5, cons836, cons837)
    def replacement1507(x, u, p):
        rubi.append(1507)
        return Int(ExpandToSum(u, x)**p, x)
    rule1507 = ReplacementRule(pattern1507, replacement1507)

    if load_rule:
        matcher.add(pattern1507, 1507)
    pattern1508 = Pattern(Integral(u_**WC('p', S(1))*(x_*WC('d', S(1)))**WC('m', S(1)), x_), cons27, cons21, cons5, cons836, cons837)
    def replacement1508(d, p, u, x, m):
        rubi.append(1508)
        return Int((d*x)**m*ExpandToSum(u, x)**p, x)
    rule1508 = ReplacementRule(pattern1508, replacement1508)

    if load_rule:
        matcher.add(pattern1508, 1508)
    pattern1509 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), cons5, cons50, cons823, cons838, cons839)
    def replacement1509(v, q, p, u, x):
        rubi.append(1509)
        return Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x)
    rule1509 = ReplacementRule(pattern1509, replacement1509)

    if load_rule:
        matcher.add(pattern1509, 1509)
    pattern1510 = Pattern(Integral(u_**WC('q', S(1))*v_**WC('p', S(1)), x_), cons5, cons50, cons823, cons840, cons841)
    def replacement1510(v, q, p, u, x):
        rubi.append(1510)
        return Int(ExpandToSum(u, x)**q*ExpandToSum(v, x)**p, x)
    rule1510 = ReplacementRule(pattern1510, replacement1510)

    if load_rule:
        matcher.add(pattern1510, 1510)
    pattern1511 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), cons21, cons5, cons50, cons842, cons836, cons843)
    def replacement1511(q, p, z, u, x, m):
        rubi.append(1511)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x)
    rule1511 = ReplacementRule(pattern1511, replacement1511)

    if load_rule:
        matcher.add(pattern1511, 1511)
    pattern1512 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_**WC('q', S(1)), x_), cons21, cons5, cons50, cons842, cons823, cons844)
    def replacement1512(q, p, z, u, x, m):
        rubi.append(1512)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x)**q, x)
    rule1512 = ReplacementRule(pattern1512, replacement1512)

    if load_rule:
        matcher.add(pattern1512, 1512)
    pattern1513 = Pattern(Integral(u_**p_, x_), cons5, cons845, cons846)
    def replacement1513(x, u, p):
        rubi.append(1513)
        return Int(ExpandToSum(u, x)**p, x)
    rule1513 = ReplacementRule(pattern1513, replacement1513)

    if load_rule:
        matcher.add(pattern1513, 1513)
    pattern1514 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1)), x_), cons21, cons5, cons845, cons846)
    def replacement1514(m, x, p, u):
        rubi.append(1514)
        return Int(x**m*ExpandToSum(u, x)**p, x)
    rule1514 = ReplacementRule(pattern1514, replacement1514)

    if load_rule:
        matcher.add(pattern1514, 1514)
    pattern1515 = Pattern(Integral(u_**WC('p', S(1))*z_, x_), cons5, cons842, cons845, cons847, cons848)
    def replacement1515(x, z, p, u):
        rubi.append(1515)
        return Int(ExpandToSum(u, x)**p*ExpandToSum(z, x), x)
    rule1515 = ReplacementRule(pattern1515, replacement1515)

    if load_rule:
        matcher.add(pattern1515, 1515)
    pattern1516 = Pattern(Integral(u_**WC('p', S(1))*x_**WC('m', S(1))*z_, x_), cons21, cons5, cons842, cons845, cons847, cons848)
    def replacement1516(z, p, u, x, m):
        rubi.append(1516)
        return Int(x**m*ExpandToSum(u, x)**p*ExpandToSum(z, x), x)
    rule1516 = ReplacementRule(pattern1516, replacement1516)

    if load_rule:
        matcher.add(pattern1516, 1516)
    pattern1517 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons7, cons48, cons125, cons208, cons209, cons21, cons4, cons849, cons850, cons851, cons852)
    def replacement1517(q, f, e, h, g, x, r, c, a, m, n):
        rubi.append(1517)
        return -Simp((S(2)*a*g + S(4)*a*h*x**(n/S(4)) - S(2)*c*f*x**(n/S(2)))/(a*c*n*sqrt(a + c*x**n)), x)
    rule1517 = ReplacementRule(pattern1517, replacement1517)

    if load_rule:
        matcher.add(pattern1517, 1517)
    pattern1518 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('n', S(1))*WC('h', S(1)) + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)))/(a_ + x_**WC('n', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons851, cons849, cons850, cons852)
    def replacement1518(d, q, f, e, h, g, x, r, c, a, m, n):
        rubi.append(1518)
        return Dist(x**(-m)*(d*x)**m, Int(x**m*(e + f*x**(n/S(4)) + g*x**(S(3)*n/S(4)) + h*x**n)/(a + c*x**n)**(S(3)/2), x), x)
    rule1518 = ReplacementRule(pattern1518, replacement1518)

    if load_rule:
        matcher.add(pattern1518, 1518)
    def With1519(Pq, p, x, b, c, a, m):
        n = Denominator(p)
        rubi.append(1519)
        return Dist(n/b, Subst(Int(x**(n*p + n + S(-1))*(-a*c/b + c*x**n/b)**m*ReplaceAll(Pq, Rule(x, -a/b + x**n/b)), x), x, (a + b*x)**(S(1)/n)), x)
    pattern1519 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons64, cons149, cons853)
    rule1519 = ReplacementRule(pattern1519, With1519)

    if load_rule:
        matcher.add(pattern1519, 1519)
    pattern1520 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons66, cons854, cons855)
    def replacement1520(Pq, p, x, b, a, m, n):
        rubi.append(1520)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*SubstFor(x**(m + S(1)), Pq, x), x), x, x**(m + S(1))), x)
    rule1520 = ReplacementRule(pattern1520, replacement1520)

    if load_rule:
        matcher.add(pattern1520, 1520)
    pattern1521 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons464, cons856)
    def replacement1521(Pq, p, x, b, a, n):
        rubi.append(1521)
        return Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x) + Simp((a + b*x**n)**(p + S(1))*Coeff(Pq, x, n + S(-1))/(b*n*(p + S(1))), x)
    rule1521 = ReplacementRule(pattern1521, replacement1521)

    if load_rule:
        matcher.add(pattern1521, 1521)
    pattern1522 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons64, cons857)
    def replacement1522(Pq, p, x, b, c, a, m, n):
        rubi.append(1522)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x)
    rule1522 = ReplacementRule(pattern1522, replacement1522)

    if load_rule:
        matcher.add(pattern1522, 1522)
    pattern1523 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons64, cons857)
    def replacement1523(Pq, p, x, b, a, n):
        rubi.append(1523)
        return Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x)
    rule1523 = ReplacementRule(pattern1523, replacement1523)

    if load_rule:
        matcher.add(pattern1523, 1523)
    pattern1524 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons858, cons500)
    def replacement1524(Pq, p, x, b, a, m, n):
        rubi.append(1524)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
    rule1524 = ReplacementRule(pattern1524, replacement1524)

    if load_rule:
        matcher.add(pattern1524, 1524)
    pattern1525 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons858, cons500)
    def replacement1525(Pq, p, c, x, b, a, m, n):
        rubi.append(1525)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
    rule1525 = ReplacementRule(pattern1525, replacement1525)

    if load_rule:
        matcher.add(pattern1525, 1525)
    pattern1526 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons64, cons53, cons13, cons137)
    def replacement1526(Pq, p, m, x, b, a, n):
        rubi.append(1526)
        return -Dist(S(1)/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*D(Pq, x), x), x) + Simp(Pq*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
    rule1526 = ReplacementRule(pattern1526, replacement1526)

    if load_rule:
        matcher.add(pattern1526, 1526)
    pattern1527 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons27, cons21, cons4, cons5, cons64, cons859)
    def replacement1527(d, Pq, p, x, b, a, m, n):
        rubi.append(1527)
        return Dist(S(1)/d, Int((d*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x), x)
    rule1527 = ReplacementRule(pattern1527, replacement1527)

    if load_rule:
        matcher.add(pattern1527, 1527)
    pattern1528 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons64, cons859, cons860)
    def replacement1528(Pq, p, x, b, a, n):
        rubi.append(1528)
        return Int(x*(a + b*x**n)**p*ExpandToSum(Pq/x, x), x)
    rule1528 = ReplacementRule(pattern1528, replacement1528)

    if load_rule:
        matcher.add(pattern1528, 1528)
    def With1529(Pq, p, x, b, a, m, n):
        u = IntHide(Pq*x**m, x)
        rubi.append(1529)
        return -Dist(b*n*p, Int(x**(m + n)*(a + b*x**n)**(p + S(-1))*ExpandToSum(u*x**(-m + S(-1)), x), x), x) + Simp(u*(a + b*x**n)**p, x)
    pattern1529 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons148, cons244, cons163, cons861)
    rule1529 = ReplacementRule(pattern1529, With1529)

    if load_rule:
        matcher.add(pattern1529, 1529)
    def With1530(Pq, p, x, b, c, a, m, n):
        q = Expon(Pq, x)
        i = Symbol('i')
        rubi.append(1530)
        return Dist(a*n*p, Int((c*x)**m*(a + b*x**n)**(p + S(-1))*Sum_doit(x**i*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q)), x), x) + Simp((c*x)**m*(a + b*x**n)**p*Sum_doit(x**(i + S(1))*Coeff(Pq, x, i)/(i + m + n*p + S(1)), List(i, S(0), q)), x)
    pattern1530 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons64, cons521, cons13, cons163)
    rule1530 = ReplacementRule(pattern1530, With1530)

    if load_rule:
        matcher.add(pattern1530, 1530)
    def With1531(Pq, p, x, b, a, n):
        q = Expon(Pq, x)
        i = Symbol('i')
        rubi.append(1531)
        return Dist(a*n*p, Int((a + b*x**n)**(p + S(-1))*Sum_doit(x**i*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q)), x), x) + Simp((a + b*x**n)**p*Sum_doit(x**(i + S(1))*Coeff(Pq, x, i)/(i + n*p + S(1)), List(i, S(0), q)), x)
    pattern1531 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons521, cons13, cons163)
    rule1531 = ReplacementRule(pattern1531, With1531)

    if load_rule:
        matcher.add(pattern1531, 1531)

    pattern1532 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons148, cons13, cons137, cons_with_1532)
    def replacement1532(Pq, p, x, b, a, n):

        q = Expon(Pq, x)
        i = Symbol('i')
        rubi.append(1532)
        return Dist(S(1)/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*Sum_doit(x**i*(i + n*(p + S(1)) + S(1))*Coeff(Pq, x, i), List(i, S(0), q + S(-1))), x), x) + Simp((a + b*x**n)**(p + S(1))*(a*Coeff(Pq, x, q) - b*x*ExpandToSum(Pq - x**q*Coeff(Pq, x, q), x))/(a*b*n*(p + S(1))), x)
    rule1532 = ReplacementRule(pattern1532, replacement1532)

    if load_rule:
        matcher.add(pattern1532, 1532)
    pattern1533 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons148, cons13, cons137, cons862)
    def replacement1533(Pq, p, x, b, a, n):
        rubi.append(1533)
        return Dist(S(1)/(a*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*ExpandToSum(Pq*n*(p + S(1)) + D(Pq*x, x), x), x), x) - Simp(Pq*x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))), x)
    rule1533 = ReplacementRule(pattern1533, replacement1533)

    if load_rule:
        matcher.add(pattern1533, 1533)
    pattern1534 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons27, cons48, cons125, cons208, cons863)
    def replacement1534(d, f, e, x, b, a, g):
        rubi.append(1534)
        return -Simp((S(2)*a*f + S(4)*a*g*x - S(2)*b*e*x**S(2))/(S(4)*a*b*sqrt(a + b*x**S(4))), x)
    rule1534 = ReplacementRule(pattern1534, replacement1534)

    if load_rule:
        matcher.add(pattern1534, 1534)
    pattern1535 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons27, cons125, cons208, cons863)
    def replacement1535(d, f, x, b, a, g):
        rubi.append(1535)
        return -Simp((f + S(2)*g*x)/(S(2)*b*sqrt(a + b*x**S(4))), x)
    rule1535 = ReplacementRule(pattern1535, replacement1535)

    if load_rule:
        matcher.add(pattern1535, 1535)
    pattern1536 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons27, cons48, cons208, cons863)
    def replacement1536(d, e, x, b, a, g):
        rubi.append(1536)
        return -Simp(x*(S(2)*a*g - b*e*x)/(S(2)*a*b*sqrt(a + b*x**S(4))), x)
    rule1536 = ReplacementRule(pattern1536, replacement1536)

    if load_rule:
        matcher.add(pattern1536, 1536)
    pattern1537 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons48, cons125, cons209, cons864)
    def replacement1537(f, e, h, x, b, a):
        rubi.append(1537)
        return -Simp((f - S(2)*h*x**S(3))/(S(2)*b*sqrt(a + b*x**S(4))), x)
    rule1537 = ReplacementRule(pattern1537, replacement1537)

    if load_rule:
        matcher.add(pattern1537, 1537)
    pattern1538 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons48, cons209, cons864)
    def replacement1538(e, h, x, b, a):
        rubi.append(1538)
        return Simp(h*x**S(3)/(b*sqrt(a + b*x**S(4))), x)
    rule1538 = ReplacementRule(pattern1538, replacement1538)

    if load_rule:
        matcher.add(pattern1538, 1538)
    pattern1539 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons27, cons48, cons125, cons208, cons209, cons864, cons863)
    def replacement1539(d, f, e, h, x, b, a, g):
        rubi.append(1539)
        return -Simp((a*f - S(2)*a*h*x**S(3) - S(2)*b*d*x)/(S(2)*a*b*sqrt(a + b*x**S(4))), x)
    rule1539 = ReplacementRule(pattern1539, replacement1539)

    if load_rule:
        matcher.add(pattern1539, 1539)
    pattern1540 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons27, cons48, cons208, cons209, cons864, cons863)
    def replacement1540(d, e, h, x, b, a, g):
        rubi.append(1540)
        return Simp(x*(a*h*x**S(2) + b*d)/(a*b*sqrt(a + b*x**S(4))), x)
    rule1540 = ReplacementRule(pattern1540, replacement1540)

    if load_rule:
        matcher.add(pattern1540, 1540)

    pattern1541 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons148, cons13, cons137, cons_with_1541)
    def replacement1541(Pq, p, x, b, a, n):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        R = PolynomialRemainder(Pq*b**(Floor((q + S(-1))/n) + S(1)), a + b*x**n, x)
        rubi.append(1541)
        return Dist(b**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)), Int((a + b*x**n)**(p + 1)*ExpandToSum(Q*a*n*(p + 1) + R*n*(p + 1) + D(R*x, x), x), x), x) - Simp(R*b**(-Floor((q - 1)/n) - 1)*x*(a + b*x**n)**(p + 1)/(a*n*(p + 1)), x)
    rule1541 = ReplacementRule(pattern1541, replacement1541)

    if load_rule:
        matcher.add(pattern1541, 1541)
    def With1542(Pq, p, x, b, a, m, n):
        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)
        R = PolynomialRemainder(Pq*a*b**(Floor((q + S(-1))/n) + S(1))*x**m, a + b*x**n, x)
        rubi.append(1542)
        return Dist(b**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)), Int(x**m*(a + b*x**n)**(p + 1)*ExpandToSum(Q*n*x**(-m)*(p + 1) + Sum_doit(x**(i - m)*(i + n*(p + 1) + 1)*Coeff(R, x, i)/a, List(i, 0, n - 1)), x), x), x) - Simp(R*b**(-Floor((q - 1)/n) - 1)*x*(a + b*x**n)**(p + 1)/(a**2*n*(p + 1)), x)
    pattern1542 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons64, cons148, cons13, cons137, cons84)
    rule1542 = ReplacementRule(pattern1542, With1542)

    if load_rule:
        matcher.add(pattern1542, 1542)

    pattern1543 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons858, cons148, cons17, cons_with_1543)
    def replacement1543(Pq, p, m, x, b, a, n):

        g = GCD(m + S(1), n)
        rubi.append(1543)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)
    rule1543 = ReplacementRule(pattern1543, replacement1543)

    if load_rule:
        matcher.add(pattern1543, 1543)
    pattern1544 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons865)
    def replacement1544(B, A, x, b, a):
        rubi.append(1544)
        return Dist(B**S(3)/b, Int(S(1)/(A**S(2) - A*B*x + B**S(2)*x**S(2)), x), x)
    rule1544 = ReplacementRule(pattern1544, replacement1544)

    if load_rule:
        matcher.add(pattern1544, 1544)
    def With1545(B, A, x, b, a):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        rubi.append(1545)
        return Dist(r/(S(3)*a*s), Int((r*(S(2)*A*s + B*r) + s*x*(-A*s + B*r))/(r**S(2) - r*s*x + s**S(2)*x**S(2)), x), x) - Dist(r*(-A*s + B*r)/(S(3)*a*s), Int(S(1)/(r + s*x), x), x)
    pattern1545 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons866, cons468)
    rule1545 = ReplacementRule(pattern1545, With1545)

    if load_rule:
        matcher.add(pattern1545, 1545)
    def With1546(B, A, x, b, a):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        rubi.append(1546)
        return -Dist(r/(S(3)*a*s), Int((r*(-S(2)*A*s + B*r) - s*x*(A*s + B*r))/(r**S(2) + r*s*x + s**S(2)*x**S(2)), x), x) + Dist(r*(A*s + B*r)/(S(3)*a*s), Int(S(1)/(r - s*x), x), x)
    pattern1546 = Pattern(Integral((A_ + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons866, cons469)
    rule1546 = ReplacementRule(pattern1546, With1546)

    if load_rule:
        matcher.add(pattern1546, 1546)
    pattern1547 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons867, cons868)
    def replacement1547(B, A, x, b, C, a):
        rubi.append(1547)
        return -Dist(C**S(2)/b, Int(S(1)/(B - C*x), x), x)
    rule1547 = ReplacementRule(pattern1547, replacement1547)

    if load_rule:
        matcher.add(pattern1547, 1547)
    def With1548(B, A, C, x, b, a):
        q = a**(S(1)/3)/b**(S(1)/3)
        rubi.append(1548)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1548 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons869)
    rule1548 = ReplacementRule(pattern1548, With1548)

    if load_rule:
        matcher.add(pattern1548, 1548)
    def With1549(B, x, b, C, a):
        q = a**(S(1)/3)/b**(S(1)/3)
        rubi.append(1549)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1549 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons870)
    rule1549 = ReplacementRule(pattern1549, With1549)

    if load_rule:
        matcher.add(pattern1549, 1549)
    def With1550(A, x, b, C, a):
        q = a**(S(1)/3)/b**(S(1)/3)
        rubi.append(1550)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1550 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons871)
    rule1550 = ReplacementRule(pattern1550, With1550)

    if load_rule:
        matcher.add(pattern1550, 1550)
    def With1551(B, A, C, x, b, a):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1551)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1551 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons872)
    rule1551 = ReplacementRule(pattern1551, With1551)

    if load_rule:
        matcher.add(pattern1551, 1551)
    def With1552(B, x, b, C, a):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1552)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1552 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons873)
    rule1552 = ReplacementRule(pattern1552, With1552)

    if load_rule:
        matcher.add(pattern1552, 1552)
    def With1553(A, x, b, C, a):
        q = (-a)**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1553)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1553 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons874)
    rule1553 = ReplacementRule(pattern1553, With1553)

    if load_rule:
        matcher.add(pattern1553, 1553)
    def With1554(B, A, C, x, b, a):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        rubi.append(1554)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1554 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons875)
    rule1554 = ReplacementRule(pattern1554, With1554)

    if load_rule:
        matcher.add(pattern1554, 1554)
    def With1555(B, x, b, C, a):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        rubi.append(1555)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1555 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons876)
    rule1555 = ReplacementRule(pattern1555, With1555)

    if load_rule:
        matcher.add(pattern1555, 1555)
    def With1556(A, x, b, C, a):
        q = (-a)**(S(1)/3)/b**(S(1)/3)
        rubi.append(1556)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1556 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons877)
    rule1556 = ReplacementRule(pattern1556, With1556)

    if load_rule:
        matcher.add(pattern1556, 1556)
    def With1557(B, A, C, x, b, a):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1557)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1557 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons878)
    rule1557 = ReplacementRule(pattern1557, With1557)

    if load_rule:
        matcher.add(pattern1557, 1557)
    def With1558(B, x, b, C, a):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1558)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1558 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons879)
    rule1558 = ReplacementRule(pattern1558, With1558)

    if load_rule:
        matcher.add(pattern1558, 1558)
    def With1559(A, x, b, C, a):
        q = a**(S(1)/3)/(-b)**(S(1)/3)
        rubi.append(1559)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1559 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons880)
    rule1559 = ReplacementRule(pattern1559, With1559)

    if load_rule:
        matcher.add(pattern1559, 1559)
    def With1560(B, A, C, x, b, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1560)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1560 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons881)
    rule1560 = ReplacementRule(pattern1560, With1560)

    if load_rule:
        matcher.add(pattern1560, 1560)
    def With1561(B, x, b, C, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1561)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1561 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons882)
    rule1561 = ReplacementRule(pattern1561, With1561)

    if load_rule:
        matcher.add(pattern1561, 1561)
    def With1562(A, x, b, C, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1562)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1562 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons883)
    rule1562 = ReplacementRule(pattern1562, With1562)

    if load_rule:
        matcher.add(pattern1562, 1562)
    def With1563(B, A, C, x, b, a):
        q = Rt(a/b, S(3))
        rubi.append(1563)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1563 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons884)
    rule1563 = ReplacementRule(pattern1563, With1563)

    if load_rule:
        matcher.add(pattern1563, 1563)
    def With1564(B, x, b, C, a):
        q = Rt(a/b, S(3))
        rubi.append(1564)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist((B + C*q)/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1564 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons885)
    rule1564 = ReplacementRule(pattern1564, With1564)

    if load_rule:
        matcher.add(pattern1564, 1564)
    def With1565(A, x, b, C, a):
        q = Rt(a/b, S(3))
        rubi.append(1565)
        return Dist(C/b, Int(S(1)/(q + x), x), x) + Dist(C*q/b, Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1565 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons886)
    rule1565 = ReplacementRule(pattern1565, With1565)

    if load_rule:
        matcher.add(pattern1565, 1565)
    def With1566(B, A, C, x, b, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1566)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1566 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons887)
    rule1566 = ReplacementRule(pattern1566, With1566)

    if load_rule:
        matcher.add(pattern1566, 1566)
    def With1567(B, x, b, C, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1567)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1567 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons888)
    rule1567 = ReplacementRule(pattern1567, With1567)

    if load_rule:
        matcher.add(pattern1567, 1567)
    def With1568(A, x, b, C, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1568)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1568 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons889)
    rule1568 = ReplacementRule(pattern1568, With1568)

    if load_rule:
        matcher.add(pattern1568, 1568)
    def With1569(B, A, C, x, b, a):
        q = Rt(-a/b, S(3))
        rubi.append(1569)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1569 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons890)
    rule1569 = ReplacementRule(pattern1569, With1569)

    if load_rule:
        matcher.add(pattern1569, 1569)
    def With1570(B, x, b, C, a):
        q = Rt(-a/b, S(3))
        rubi.append(1570)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) + Dist((B - C*q)/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1570 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons891)
    rule1570 = ReplacementRule(pattern1570, With1570)

    if load_rule:
        matcher.add(pattern1570, 1570)
    def With1571(A, x, b, C, a):
        q = Rt(-a/b, S(3))
        rubi.append(1571)
        return -Dist(C/b, Int(S(1)/(q - x), x), x) - Dist(C*q/b, Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1571 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons892)
    rule1571 = ReplacementRule(pattern1571, With1571)

    if load_rule:
        matcher.add(pattern1571, 1571)
    pattern1572 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons893)
    def replacement1572(B, A, C, x, b, a):
        rubi.append(1572)
        return Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x) + Int((A + B*x)/(a + b*x**S(3)), x)
    rule1572 = ReplacementRule(pattern1572, replacement1572)

    if load_rule:
        matcher.add(pattern1572, 1572)
    pattern1573 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons894)
    def replacement1573(B, x, b, C, a):
        rubi.append(1573)
        return Dist(B, Int(x/(a + b*x**S(3)), x), x) + Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x)
    rule1573 = ReplacementRule(pattern1573, replacement1573)

    if load_rule:
        matcher.add(pattern1573, 1573)
    pattern1574 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons895)
    def replacement1574(A, x, b, C, a):
        rubi.append(1574)
        return Dist(A, Int(S(1)/(a + b*x**S(3)), x), x) + Dist(C, Int(x**S(2)/(a + b*x**S(3)), x), x)
    rule1574 = ReplacementRule(pattern1574, replacement1574)

    if load_rule:
        matcher.add(pattern1574, 1574)
    def With1575(B, A, C, x, b, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1575)
        return Dist(q**S(2)/a, Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1575 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons896)
    rule1575 = ReplacementRule(pattern1575, With1575)

    if load_rule:
        matcher.add(pattern1575, 1575)
    def With1576(B, x, b, C, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1576)
        return Dist(C*q**S(3)/a, Int(x/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1576 = Pattern(Integral(x_*(x_*WC('C', S(1)) + WC('B', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons897)
    rule1576 = ReplacementRule(pattern1576, With1576)

    if load_rule:
        matcher.add(pattern1576, 1576)
    def With1577(A, x, b, C, a):
        q = (a/b)**(S(1)/3)
        rubi.append(1577)
        return Dist(q**S(2)/a, Int((A + C*q*x)/(q**S(2) - q*x + x**S(2)), x), x)
    pattern1577 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons898)
    rule1577 = ReplacementRule(pattern1577, With1577)

    if load_rule:
        matcher.add(pattern1577, 1577)
    def With1578(B, A, C, x, b, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1578)
        return Dist(q/a, Int((A*q + x*(A + B*q))/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1578 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons899)
    rule1578 = ReplacementRule(pattern1578, With1578)

    if load_rule:
        matcher.add(pattern1578, 1578)
    def With1579(B, x, b, C, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1579)
        return Dist(B*q**S(2)/a, Int(x/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1579 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons900)
    rule1579 = ReplacementRule(pattern1579, With1579)

    if load_rule:
        matcher.add(pattern1579, 1579)
    def With1580(A, x, b, C, a):
        q = (-a/b)**(S(1)/3)
        rubi.append(1580)
        return Dist(A*q/a, Int((q + x)/(q**S(2) + q*x + x**S(2)), x), x)
    pattern1580 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons901)
    rule1580 = ReplacementRule(pattern1580, With1580)

    if load_rule:
        matcher.add(pattern1580, 1580)

    pattern1581 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons866, cons902, cons903, cons_with_1581)
    def replacement1581(B, A, C, x, b, a):

        q = (a/b)**(S(1)/3)
        rubi.append(1581)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A + B*q - C*q**S(2)) - x*(A - B*q - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) + Dist(q*(A - B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)
    rule1581 = ReplacementRule(pattern1581, replacement1581)

    if load_rule:
        matcher.add(pattern1581, 1581)

    pattern1582 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons902, cons903, cons_with_1582)
    def replacement1582(B, x, b, C, a):

        q = (a/b)**(S(1)/3)
        rubi.append(1582)
        return Dist(q/(S(3)*a), Int((q*(B*q - C*q**S(2)) + x*(B*q + S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) - Dist(q*(B*q - C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)
    rule1582 = ReplacementRule(pattern1582, replacement1582)

    if load_rule:
        matcher.add(pattern1582, 1582)

    pattern1583 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons902, cons903, cons_with_1583)
    def replacement1583(A, x, b, C, a):

        q = (a/b)**(S(1)/3)
        rubi.append(1583)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - C*q**S(2)) - x*(A - S(2)*C*q**S(2)))/(q**S(2) - q*x + x**S(2)), x), x) + Dist(q*(A + C*q**S(2))/(S(3)*a), Int(S(1)/(q + x), x), x)
    rule1583 = ReplacementRule(pattern1583, replacement1583)

    if load_rule:
        matcher.add(pattern1583, 1583)

    pattern1584 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons35, cons36, cons866, cons902, cons904, cons_with_1584)
    def replacement1584(B, A, C, x, b, a):

        q = (-a/b)**(S(1)/3)
        rubi.append(1584)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - B*q - C*q**S(2)) + x*(A + B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(A + B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)
    rule1584 = ReplacementRule(pattern1584, replacement1584)

    if load_rule:
        matcher.add(pattern1584, 1584)

    pattern1585 = Pattern(Integral(x_*(B_ + x_*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons35, cons36, cons902, cons904, cons_with_1585)
    def replacement1585(B, x, b, C, a):

        q = (-a/b)**(S(1)/3)
        rubi.append(1585)
        return Dist(q/(S(3)*a), Int((-q*(B*q + C*q**S(2)) + x*(B*q - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(B*q + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)
    rule1585 = ReplacementRule(pattern1585, replacement1585)

    if load_rule:
        matcher.add(pattern1585, 1585)

    pattern1586 = Pattern(Integral((A_ + x_**S(2)*WC('C', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons34, cons36, cons902, cons904, cons_with_1586)
    def replacement1586(A, x, b, C, a):

        q = (-a/b)**(S(1)/3)
        rubi.append(1586)
        return Dist(q/(S(3)*a), Int((q*(S(2)*A - C*q**S(2)) + x*(A - S(2)*C*q**S(2)))/(q**S(2) + q*x + x**S(2)), x), x) + Dist(q*(A + C*q**S(2))/(S(3)*a), Int(S(1)/(q - x), x), x)
    rule1586 = ReplacementRule(pattern1586, replacement1586)

    if load_rule:
        matcher.add(pattern1586, 1586)

    pattern1587 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons64, cons674, cons905, cons_with_1587)
    def replacement1587(Pq, x, b, c, a, m, n):

        v = Sum_doit(c**(-ii)*(c*x)**(ii + m)*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        rubi.append(1587)
        return Int(v, x)
    rule1587 = ReplacementRule(pattern1587, replacement1587)

    if load_rule:
        matcher.add(pattern1587, 1587)

    pattern1588 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons64, cons674, cons905, cons_with_1588)
    def replacement1588(Pq, x, b, a, n):

        v = Sum_doit(x**ii*(x**(n/S(2))*Coeff(Pq, x, ii + n/S(2)) + Coeff(Pq, x, ii))/(a + b*x**n), List(ii, S(0), n/S(2) + S(-1)))
        rubi.append(1588)
        return Int(v, x)
    rule1588 = ReplacementRule(pattern1588, replacement1588)

    if load_rule:
        matcher.add(pattern1588, 1588)
    def With1589(d, c, x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        rubi.append(1589)
        return Simp(S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(S(1) + sqrt(S(3))))), x) - Simp(S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(r*x + s)*EllipticE(asin((r*x + s*(-sqrt(S(3)) + S(1)))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(r**S(2)*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3))), x)
    pattern1589 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons481, cons906)
    rule1589 = ReplacementRule(pattern1589, With1589)

    if load_rule:
        matcher.add(pattern1589, 1589)
    def With1590(d, c, x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        rubi.append(1590)
        return Dist(d/r, Int((r*x + s*(-sqrt(S(3)) + S(1)))/sqrt(a + b*x**S(3)), x), x) + Dist((c*r - d*s*(-sqrt(S(3)) + S(1)))/r, Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1590 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons481, cons907)
    rule1590 = ReplacementRule(pattern1590, With1590)

    if load_rule:
        matcher.add(pattern1590, 1590)
    def With1591(d, c, x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        rubi.append(1591)
        return Simp(S(2)*d*s**S(3)*sqrt(a + b*x**S(3))/(a*r**S(2)*(r*x + s*(-sqrt(S(3)) + S(1)))), x) + Simp(S(3)**(S(1)/4)*d*s*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticE(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(-sqrt(S(3)) + S(1)))), S(-7) + S(4)*sqrt(S(3)))/(r**S(2)*sqrt(-s*(r*x + s)/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(a + b*x**S(3))), x)
    pattern1591 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons482, cons908)
    rule1591 = ReplacementRule(pattern1591, With1591)

    if load_rule:
        matcher.add(pattern1591, 1591)
    def With1592(d, c, x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        rubi.append(1592)
        return Dist(d/r, Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x), x) + Dist((c*r - d*s*(S(1) + sqrt(S(3))))/r, Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1592 = Pattern(Integral((c_ + x_*WC('d', S(1)))/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons482, cons909)
    rule1592 = ReplacementRule(pattern1592, With1592)

    if load_rule:
        matcher.add(pattern1592, 1592)
    def With1593(d, c, x, b, a):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        rubi.append(1593)
        return Simp(d*s**S(3)*x*(S(1) + sqrt(S(3)))*sqrt(a + b*x**S(6))/(S(2)*a*r**S(2)*(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), x) - Simp(S(3)**(S(1)/4)*d*s*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticE(acos((r*x**S(2)*(-sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(2)*r**S(2)*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6))), x)
    pattern1593 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons910)
    rule1593 = ReplacementRule(pattern1593, With1593)

    if load_rule:
        matcher.add(pattern1593, 1593)
    def With1594(d, c, x, b, a):
        q = Rt(b/a, S(3))
        rubi.append(1594)
        return Dist(d/(S(2)*q**S(2)), Int((S(2)*q**S(2)*x**S(4) - sqrt(S(3)) + S(1))/sqrt(a + b*x**S(6)), x), x) + Dist((S(2)*c*q**S(2) - d*(-sqrt(S(3)) + S(1)))/(S(2)*q**S(2)), Int(S(1)/sqrt(a + b*x**S(6)), x), x)
    pattern1594 = Pattern(Integral((c_ + x_**S(4)*WC('d', S(1)))/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons911)
    rule1594 = ReplacementRule(pattern1594, With1594)

    if load_rule:
        matcher.add(pattern1594, 1594)
    pattern1595 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons912)
    def replacement1595(d, c, x, b, a):
        rubi.append(1595)
        return -Simp(c*d*x**S(3)*sqrt(-(c - d*x**S(2))**S(2)/(c*d*x**S(2)))*sqrt(-d**S(2)*(a + b*x**S(8))/(b*c**S(2)*x**S(4)))*EllipticF(asin(sqrt((sqrt(S(2))*c**S(2) + S(2)*c*d*x**S(2) + sqrt(S(2))*d**S(2)*x**S(4))/(c*d*x**S(2)))/S(2)), S(-2) + S(2)*sqrt(S(2)))/(sqrt(sqrt(S(2)) + S(2))*sqrt(a + b*x**S(8))*(c - d*x**S(2))), x)
    rule1595 = ReplacementRule(pattern1595, replacement1595)

    if load_rule:
        matcher.add(pattern1595, 1595)
    pattern1596 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons913)
    def replacement1596(d, c, x, b, a):
        rubi.append(1596)
        return -Dist((-c*Rt(b/a, S(4)) + d)/(S(2)*Rt(b/a, S(4))), Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x) + Dist((c*Rt(b/a, S(4)) + d)/(S(2)*Rt(b/a, S(4))), Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x), x)
    rule1596 = ReplacementRule(pattern1596, replacement1596)

    if load_rule:
        matcher.add(pattern1596, 1596)
    pattern1597 = Pattern(Integral(Pq_/(x_*sqrt(a_ + x_**n_*WC('b', S(1)))), x_), cons2, cons3, cons64, cons148, cons914)
    def replacement1597(Pq, x, b, a, n):
        rubi.append(1597)
        return Dist(Coeff(Pq, x, S(0)), Int(S(1)/(x*sqrt(a + b*x**n)), x), x) + Int(ExpandToSum((Pq - Coeff(Pq, x, S(0)))/x, x)/sqrt(a + b*x**n), x)
    rule1597 = ReplacementRule(pattern1597, replacement1597)

    if load_rule:
        matcher.add(pattern1597, 1597)
    def With1598(Pq, p, x, b, c, a, m, n):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        rubi.append(1598)
        return Int(Sum_doit(c**(-j)*(c*x)**(j + m)*(a + b*x**n)**p*Sum_doit(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)
    pattern1598 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons64, cons674, cons915)
    rule1598 = ReplacementRule(pattern1598, With1598)

    if load_rule:
        matcher.add(pattern1598, 1598)
    def With1599(Pq, p, x, b, a, n):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        rubi.append(1599)
        return Int(Sum_doit(x**j*(a + b*x**n)**p*Sum_doit(x**(k*n/S(2))*Coeff(Pq, x, j + k*n/S(2)), List(k, S(0), S(1) + S(2)*(-j + q)/n)), List(j, S(0), n/S(2) + S(-1))), x)
    pattern1599 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons64, cons674, cons915)
    rule1599 = ReplacementRule(pattern1599, With1599)

    if load_rule:
        matcher.add(pattern1599, 1599)
    pattern1600 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons64, cons148, cons916)
    def replacement1600(Pq, p, x, b, a, n):
        rubi.append(1600)
        return Dist(Coeff(Pq, x, n + S(-1)), Int(x**(n + S(-1))*(a + b*x**n)**p, x), x) + Int((a + b*x**n)**p*ExpandToSum(Pq - x**(n + S(-1))*Coeff(Pq, x, n + S(-1)), x), x)
    rule1600 = ReplacementRule(pattern1600, replacement1600)

    if load_rule:
        matcher.add(pattern1600, 1600)
    pattern1601 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons64, cons85)
    def replacement1601(Pq, x, b, c, a, m, n):
        rubi.append(1601)
        return Int(ExpandIntegrand(Pq*(c*x)**m/(a + b*x**n), x), x)
    rule1601 = ReplacementRule(pattern1601, replacement1601)

    if load_rule:
        matcher.add(pattern1601, 1601)
    pattern1602 = Pattern(Integral(Pq_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons64, cons85)
    def replacement1602(Pq, x, b, a, n):
        rubi.append(1602)
        return Int(ExpandIntegrand(Pq/(a + b*x**n), x), x)
    rule1602 = ReplacementRule(pattern1602, replacement1602)

    if load_rule:
        matcher.add(pattern1602, 1602)

    pattern1603 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons64, cons148, cons31, cons94, cons917, cons_with_1603)
    def replacement1603(Pq, p, x, b, c, a, m, n):

        Pq0 = Coeff(Pq, x, S(0))
        rubi.append(1603)
        return Dist(S(1)/(S(2)*a*c*(m + S(1))), Int((c*x)**(m + S(1))*(a + b*x**n)**p*ExpandToSum(-S(2)*Pq0*b*x**(n + S(-1))*(m + n*(p + S(1)) + S(1)) + S(2)*a*(Pq - Pq0)*(m + S(1))/x, x), x), x) + Simp(Pq0*(c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))), x)
    rule1603 = ReplacementRule(pattern1603, replacement1603)

    if load_rule:
        matcher.add(pattern1603, 1603)

    pattern1604 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons64, cons148, cons_with_1604)
    def replacement1604(Pq, p, x, b, c, a, m, n):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1604)
        return Dist(1/(b*(m + n*p + q + 1)), Int((c*x)**m*(a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(m - n + q + 1) + b*(Pq - Pqq*x**q)*(m + n*p + q + 1), x), x), x) + Simp(Pqq*c**(n - q - 1)*(c*x)**(m - n + q + 1)*(a + b*x**n)**(p + 1)/(b*(m + n*p + q + 1)), x)
    rule1604 = ReplacementRule(pattern1604, replacement1604)

    if load_rule:
        matcher.add(pattern1604, 1604)

    pattern1605 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons64, cons148, cons_with_1605)
    def replacement1605(Pq, p, x, b, a, n):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1605)
        return Dist(1/(b*(n*p + q + 1)), Int((a + b*x**n)**p*ExpandToSum(-Pqq*a*x**(-n + q)*(-n + q + 1) + b*(Pq - Pqq*x**q)*(n*p + q + 1), x), x), x) + Simp(Pqq*x**(-n + q + 1)*(a + b*x**n)**(p + 1)/(b*(n*p + q + 1)), x)
    rule1605 = ReplacementRule(pattern1605, replacement1605)

    if load_rule:
        matcher.add(pattern1605, 1605)
    def With1606(Pq, p, m, x, b, a, n):
        q = Expon(Pq, x)
        rubi.append(1606)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x)
    pattern1606 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons64, cons196, cons17)
    rule1606 = ReplacementRule(pattern1606, With1606)

    if load_rule:
        matcher.add(pattern1606, 1606)
    def With1607(Pq, p, x, b, c, a, m, n):
        g = Denominator(m)
        q = Expon(Pq, x)
        rubi.append(1607)
        return -Dist(g/c, Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*c**(-n)*x**(-g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/c)), x), x), x, (c*x)**(-S(1)/g)), x)
    pattern1607 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons64, cons196, cons367)
    rule1607 = ReplacementRule(pattern1607, With1607)

    if load_rule:
        matcher.add(pattern1607, 1607)
    def With1608(Pq, p, x, b, c, a, m, n):
        q = Expon(Pq, x)
        rubi.append(1608)
        return -Dist((c*x)**m*(S(1)/x)**m, Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x), x)
    pattern1608 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons64, cons196, cons356)
    rule1608 = ReplacementRule(pattern1608, With1608)

    if load_rule:
        matcher.add(pattern1608, 1608)
    def With1609(Pq, p, m, x, b, a, n):
        g = Denominator(n)
        rubi.append(1609)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
    pattern1609 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons5, cons64, cons489)
    rule1609 = ReplacementRule(pattern1609, With1609)

    if load_rule:
        matcher.add(pattern1609, 1609)
    def With1610(Pq, p, x, b, a, n):
        g = Denominator(n)
        rubi.append(1610)
        return Dist(g, Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
    pattern1610 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons5, cons64, cons489)
    rule1610 = ReplacementRule(pattern1610, With1610)

    if load_rule:
        matcher.add(pattern1610, 1610)
    pattern1611 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons64, cons489)
    def replacement1611(Pq, p, c, x, b, a, m, n):
        rubi.append(1611)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
    rule1611 = ReplacementRule(pattern1611, replacement1611)

    if load_rule:
        matcher.add(pattern1611, 1611)
    pattern1612 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons858, cons541, cons23)
    def replacement1612(Pq, p, m, x, b, a, n):
        rubi.append(1612)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
    rule1612 = ReplacementRule(pattern1612, replacement1612)

    if load_rule:
        matcher.add(pattern1612, 1612)
    pattern1613 = Pattern(Integral(Pq_*(c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons858, cons541, cons23)
    def replacement1613(Pq, p, c, x, b, a, m, n):
        rubi.append(1613)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(Pq*x**m*(a + b*x**n)**p, x), x)
    rule1613 = ReplacementRule(pattern1613, replacement1613)

    if load_rule:
        matcher.add(pattern1613, 1613)
    pattern1614 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons34, cons35, cons21, cons4, cons5, cons53)
    def replacement1614(B, p, A, x, b, a, m, n):
        rubi.append(1614)
        return Dist(A, Int((a + b*x**n)**p, x), x) + Dist(B, Int(x**m*(a + b*x**n)**p, x), x)
    rule1614 = ReplacementRule(pattern1614, replacement1614)

    if load_rule:
        matcher.add(pattern1614, 1614)
    pattern1615 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons918)
    def replacement1615(Pq, p, x, b, c, a, m, n):
        rubi.append(1615)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a + b*x**n)**p, x), x)
    rule1615 = ReplacementRule(pattern1615, replacement1615)

    if load_rule:
        matcher.add(pattern1615, 1615)
    pattern1616 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons5, cons918)
    def replacement1616(Pq, p, x, b, a, n):
        rubi.append(1616)
        return Int(ExpandIntegrand(Pq*(a + b*x**n)**p, x), x)
    rule1616 = ReplacementRule(pattern1616, replacement1616)

    if load_rule:
        matcher.add(pattern1616, 1616)
    pattern1617 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons21, cons4, cons5, cons554, cons919)
    def replacement1617(Pq, v, p, u, x, b, a, m, n):
        rubi.append(1617)
        return Dist(u**m*v**(-m)/Coeff(v, x, S(1)), Subst(Int(x**m*(a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v), x)
    rule1617 = ReplacementRule(pattern1617, replacement1617)

    if load_rule:
        matcher.add(pattern1617, 1617)
    pattern1618 = Pattern(Integral(Pq_*(a_ + v_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons5, cons552, cons919)
    def replacement1618(Pq, v, p, x, b, a, n):
        rubi.append(1618)
        return Dist(S(1)/Coeff(v, x, S(1)), Subst(Int((a + b*x**n)**p*SubstFor(v, Pq, x), x), x, v), x)
    rule1618 = ReplacementRule(pattern1618, replacement1618)

    if load_rule:
        matcher.add(pattern1618, 1618)
    pattern1619 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons64, cons55, cons494)
    def replacement1619(a2, Pq, p, x, b2, c, m, n, b1, a1):
        rubi.append(1619)
        return Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule1619 = ReplacementRule(pattern1619, replacement1619)

    if load_rule:
        matcher.add(pattern1619, 1619)
    pattern1620 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons4, cons5, cons64, cons55, cons494)
    def replacement1620(Pq, p, x, b2, b1, a2, n, a1):
        rubi.append(1620)
        return Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x)
    rule1620 = ReplacementRule(pattern1620, replacement1620)

    if load_rule:
        matcher.add(pattern1620, 1620)
    pattern1621 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons7, cons21, cons4, cons5, cons64, cons55)
    def replacement1621(a2, Pq, p, x, b2, c, m, n, b1, a1):
        rubi.append(1621)
        return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int(Pq*(c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)
    rule1621 = ReplacementRule(pattern1621, replacement1621)

    if load_rule:
        matcher.add(pattern1621, 1621)
    pattern1622 = Pattern(Integral(Pq_*(a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons57, cons58, cons59, cons60, cons4, cons5, cons64, cons55)
    def replacement1622(Pq, p, x, b2, b1, a2, n, a1):
        rubi.append(1622)
        return Dist((a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p)), Int(Pq*(a1*a2 + b1*b2*x**(S(2)*n))**p, x), x)
    rule1622 = ReplacementRule(pattern1622, replacement1622)

    if load_rule:
        matcher.add(pattern1622, 1622)
    pattern1623 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons46, cons920, cons921)
    def replacement1623(d, p, c, f, e, g, n2, x, b, a, n):
        rubi.append(1623)
        return Simp(e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c), x)
    rule1623 = ReplacementRule(pattern1623, replacement1623)

    if load_rule:
        matcher.add(pattern1623, 1623)
    pattern1624 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons208, cons4, cons5, cons46, cons922, cons921)
    def replacement1624(d, p, c, e, g, n2, x, b, a, n):
        rubi.append(1624)
        return Simp(e*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c), x)
    rule1624 = ReplacementRule(pattern1624, replacement1624)

    if load_rule:
        matcher.add(pattern1624, 1624)
    pattern1625 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)) + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons46, cons923, cons924, cons66)
    def replacement1625(d, p, c, f, e, h, g, n2, x, b, a, m, n):
        rubi.append(1625)
        return Simp(e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))), x)
    rule1625 = ReplacementRule(pattern1625, replacement1625)

    if load_rule:
        matcher.add(pattern1625, 1625)
    pattern1626 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('p', S(1))*(e_ + x_**WC('n2', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons208, cons209, cons21, cons4, cons5, cons46, cons595, cons924, cons66)
    def replacement1626(d, p, c, e, h, g, n2, x, b, a, m, n):
        rubi.append(1626)
        return Simp(e*(h*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(p + S(1))/(a*c*h*(m + S(1))), x)
    rule1626 = ReplacementRule(pattern1626, replacement1626)

    if load_rule:
        matcher.add(pattern1626, 1626)
    pattern1627 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons34, cons35, cons21, cons4, cons5, cons50, cons71, cons53)
    def replacement1627(B, d, q, p, c, A, x, b, a, m, n):
        rubi.append(1627)
        return Dist(A, Int((a + b*x**n)**p*(c + d*x**n)**q, x), x) + Dist(B, Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule1627 = ReplacementRule(pattern1627, replacement1627)

    if load_rule:
        matcher.add(pattern1627, 1627)
    def With1628(d, q, p, c, x, b, a, n, Px):
        k = Denominator(n)
        rubi.append(1628)
        return Dist(k/d, Subst(Int(SimplifyIntegrand(x**(k + S(-1))*(a + b*x**(k*n))**p*ReplaceAll(Px, Rule(x, -c/d + x**k/d))**q, x), x), x, (c + d*x)**(S(1)/k)), x)
    pattern1628 = Pattern(Integral(Px_**WC('q', S(1))*((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons925, cons586, cons87)
    rule1628 = ReplacementRule(pattern1628, With1628)

    if load_rule:
        matcher.add(pattern1628, 1628)
    pattern1629 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons858, cons53)
    def replacement1629(Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1629)
        return Dist(S(1)/n, Subst(Int((a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
    rule1629 = ReplacementRule(pattern1629, replacement1629)

    if load_rule:
        matcher.add(pattern1629, 1629)
    pattern1630 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons46, cons64, cons128)
    def replacement1630(d, Pq, p, c, n2, x, b, a, m, n):
        rubi.append(1630)
        return Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1630 = ReplacementRule(pattern1630, replacement1630)

    if load_rule:
        matcher.add(pattern1630, 1630)
    pattern1631 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons64, cons128)
    def replacement1631(Pq, p, c, n2, x, b, a, n):
        rubi.append(1631)
        return Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1631 = ReplacementRule(pattern1631, replacement1631)

    if load_rule:
        matcher.add(pattern1631, 1631)
    pattern1632 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons46, cons926, cons927)
    def replacement1632(d, p, c, f, e, n2, x, b, a, n):
        rubi.append(1632)
        return Simp(d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a, x)
    rule1632 = ReplacementRule(pattern1632, replacement1632)

    if load_rule:
        matcher.add(pattern1632, 1632)
    pattern1633 = Pattern(Integral((d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons4, cons5, cons46, cons922, cons928)
    def replacement1633(d, p, c, f, n2, x, b, a, n):
        rubi.append(1633)
        return Simp(d*x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/a, x)
    rule1633 = ReplacementRule(pattern1633, replacement1633)

    if load_rule:
        matcher.add(pattern1633, 1633)
    pattern1634 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n', S(1))*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons46, cons929, cons930, cons66)
    def replacement1634(d, p, c, f, e, n2, x, n, b, a, m, g):
        rubi.append(1634)
        return Simp(d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))), x)
    rule1634 = ReplacementRule(pattern1634, replacement1634)

    if load_rule:
        matcher.add(pattern1634, 1634)
    pattern1635 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons21, cons4, cons5, cons46, cons595, cons928, cons66)
    def replacement1635(d, p, c, f, n2, x, n, b, a, m, g):
        rubi.append(1635)
        return Simp(d*(g*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*g*(m + S(1))), x)
    rule1635 = ReplacementRule(pattern1635, replacement1635)

    if load_rule:
        matcher.add(pattern1635, 1635)
    pattern1636 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons46, cons64, cons45, cons314)
    def replacement1636(d, Pq, p, c, n2, x, b, a, m, n):
        rubi.append(1636)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int(Pq*(d*x)**m*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1636 = ReplacementRule(pattern1636, replacement1636)

    if load_rule:
        matcher.add(pattern1636, 1636)
    pattern1637 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons64, cons45, cons314)
    def replacement1637(Pq, p, c, n2, x, b, a, n):
        rubi.append(1637)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int(Pq*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1637 = ReplacementRule(pattern1637, replacement1637)

    if load_rule:
        matcher.add(pattern1637, 1637)
    pattern1638 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons858, cons226, cons500)
    def replacement1638(Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1638)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
    rule1638 = ReplacementRule(pattern1638, replacement1638)

    if load_rule:
        matcher.add(pattern1638, 1638)
    pattern1639 = Pattern(Integral(Pq_*(d_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons46, cons858, cons226, cons500)
    def replacement1639(d, Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1639)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1639 = ReplacementRule(pattern1639, replacement1639)

    if load_rule:
        matcher.add(pattern1639, 1639)
    pattern1640 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons46, cons64, cons859)
    def replacement1640(d, Pq, p, c, n2, x, b, a, m, n):
        rubi.append(1640)
        return Dist(S(1)/d, Int((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x), x)
    rule1640 = ReplacementRule(pattern1640, replacement1640)

    if load_rule:
        matcher.add(pattern1640, 1640)
    pattern1641 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons64, cons859, cons860)
    def replacement1641(Pq, p, c, n2, x, b, a, n):
        rubi.append(1641)
        return Int(x*(a + b*x**n + c*x**(S(2)*n))**p*ExpandToSum(Pq/x, x), x)
    rule1641 = ReplacementRule(pattern1641, replacement1641)

    if load_rule:
        matcher.add(pattern1641, 1641)
    pattern1642 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons46, cons931, cons226, cons932, cons933)
    def replacement1642(n3, d, p, f, e, n2, x, b, n, c, a, g):
        rubi.append(1642)
        return Simp(x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
    rule1642 = ReplacementRule(pattern1642, replacement1642)

    if load_rule:
        matcher.add(pattern1642, 1642)
    pattern1643 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('n2', S(1))*WC('f', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons4, cons5, cons46, cons931, cons226, cons934, cons935)
    def replacement1643(n3, d, p, f, n2, x, b, n, c, a, g):
        rubi.append(1643)
        return Simp(x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
    rule1643 = ReplacementRule(pattern1643, replacement1643)

    if load_rule:
        matcher.add(pattern1643, 1643)
    pattern1644 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)) + x_**WC('n3', S(1))*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons208, cons4, cons5, cons46, cons931, cons226, cons932, cons936)
    def replacement1644(n3, d, p, e, n2, x, b, n, c, a, g):
        rubi.append(1644)
        return Simp(x*(S(3)*a*d - x**S(2)*(-a*e + S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
    rule1644 = ReplacementRule(pattern1644, replacement1644)

    if load_rule:
        matcher.add(pattern1644, 1644)
    pattern1645 = Pattern(Integral((d_ + x_**WC('n3', S(1))*WC('g', S(1)))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons208, cons4, cons5, cons46, cons931, cons226, cons934, cons937)
    def replacement1645(n3, d, p, n2, x, b, n, c, a, g):
        rubi.append(1645)
        return Simp(x*(S(3)*a*d - x**S(2)*(S(2)*b*d*p + S(3)*b*d))*(a + b*x**S(2) + c*x**S(4))**(p + S(1))/(S(3)*a**S(2)), x)
    rule1645 = ReplacementRule(pattern1645, replacement1645)

    if load_rule:
        matcher.add(pattern1645, 1645)
    pattern1646 = Pattern(Integral(x_**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons209, cons21, cons4, cons46, cons938, cons939, cons940, cons226, cons941, cons852)
    def replacement1646(s, q, c, f, e, h, g, n2, x, r, b, a, m, n):
        rubi.append(1646)
        return -Simp((S(2)*c*x**n*(-b*g + S(2)*c*f) + S(2)*c*(-S(2)*a*g + b*f) + S(2)*h*x**(n/S(2))*(-S(4)*a*c + b**S(2)))/(c*n*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**n + c*x**(S(2)*n))), x)
    rule1646 = ReplacementRule(pattern1646, replacement1646)

    if load_rule:
        matcher.add(pattern1646, 1646)
    pattern1647 = Pattern(Integral((d_*x_)**WC('m', S(1))*(e_ + x_**WC('q', S(1))*WC('f', S(1)) + x_**WC('r', S(1))*WC('g', S(1)) + x_**WC('s', S(1))*WC('h', S(1)))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons46, cons938, cons939, cons940, cons226, cons941, cons852)
    def replacement1647(s, d, q, c, f, e, h, g, n2, x, r, b, a, m, n):
        rubi.append(1647)
        return Dist(x**(-m)*(d*x)**m, Int(x**m*(e + f*x**(n/S(2)) + g*x**(S(3)*n/S(2)) + h*x**(S(2)*n))/(a + b*x**n + c*x**(S(2)*n))**(S(3)/2), x), x)
    rule1647 = ReplacementRule(pattern1647, replacement1647)

    if load_rule:
        matcher.add(pattern1647, 1647)

    pattern1648 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons46, cons64, cons226, cons148, cons13, cons137, cons_with_1648)
    def replacement1648(Pq, p, c, n2, x, b, a, n):

        q = Expon(Pq, x)
        i = Symbol('i')
        rubi.append(1648)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum_doit(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i))*(i + n*(S(2)*p + S(3)) + S(1)) + x**i*(-a*b*(i + S(1))*Coeff(Pq, x, i + n) + (-S(2)*a*c*(i + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(i + n*(p + S(1)) + S(1)))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1))), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Sum_doit(c*x**(i + n)*(-S(2)*a*Coeff(Pq, x, i + n) + b*Coeff(Pq, x, i)) + x**i*(-a*b*Coeff(Pq, x, i + n) + (-S(2)*a*c + b**S(2))*Coeff(Pq, x, i)), List(i, S(0), n + S(-1)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1648 = ReplacementRule(pattern1648, replacement1648)

    if load_rule:
        matcher.add(pattern1648, 1648)
    pattern1649 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons226, cons942)
    def replacement1649(d, c, f, e, x, b, a, g):
        rubi.append(1649)
        return -Simp((c*x**S(2)*(-b*f + S(2)*c*e) + c*(-S(2)*a*f + b*e) + g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1649 = ReplacementRule(pattern1649, replacement1649)

    if load_rule:
        matcher.add(pattern1649, 1649)
    pattern1650 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons226, cons942)
    def replacement1650(d, c, f, x, b, a, g):
        rubi.append(1650)
        return Simp((S(2)*a*c*f + b*c*f*x**S(2) - g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1650 = ReplacementRule(pattern1650, replacement1650)

    if load_rule:
        matcher.add(pattern1650, 1650)
    pattern1651 = Pattern(Integral((d_ + x_**S(4)*WC('g', S(1)) + x_*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons208, cons226, cons942)
    def replacement1651(d, c, e, x, b, a, g):
        rubi.append(1651)
        return -Simp((b*c*e + S(2)*c**S(2)*e*x**S(2) + g*x*(-S(4)*a*c + b**S(2)))/(c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1651 = ReplacementRule(pattern1651, replacement1651)

    if load_rule:
        matcher.add(pattern1651, 1651)
    pattern1652 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons209, cons226, cons943, cons944)
    def replacement1652(c, f, e, h, x, b, a, g):
        rubi.append(1652)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1652 = ReplacementRule(pattern1652, replacement1652)

    if load_rule:
        matcher.add(pattern1652, 1652)
    pattern1653 = Pattern(Integral(x_**S(2)*(x_**S(4)*WC('h', S(1)) + x_**S(2)*WC('g', S(1)) + WC('e', S(0)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons48, cons208, cons209, cons226, cons943, cons944)
    def replacement1653(c, e, h, x, b, a, g):
        rubi.append(1653)
        return Simp(h*x**S(3)/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1653 = ReplacementRule(pattern1653, replacement1653)

    if load_rule:
        matcher.add(pattern1653, 1653)
    pattern1654 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(4)*WC('g', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons226, cons943, cons945)
    def replacement1654(d, c, f, e, h, x, b, a, g):
        rubi.append(1654)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1654 = ReplacementRule(pattern1654, replacement1654)

    if load_rule:
        matcher.add(pattern1654, 1654)
    pattern1655 = Pattern(Integral((d_ + x_**S(6)*WC('h', S(1)) + x_**S(3)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons209, cons226, cons943, cons946)
    def replacement1655(d, c, f, e, h, x, b, a):
        rubi.append(1655)
        return Simp((S(2)*a**S(2)*c*f + a*b*c*f*x**S(2) + a*h*x**S(3)*(-S(4)*a*c + b**S(2)) + c*d*x*(-S(4)*a*c + b**S(2)))/(a*c*(-S(4)*a*c + b**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1655 = ReplacementRule(pattern1655, replacement1655)

    if load_rule:
        matcher.add(pattern1655, 1655)

    pattern1656 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons46, cons64, cons226, cons148, cons13, cons137, cons_with_1656)
    def replacement1656(Pq, p, c, n2, x, b, a, n):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        rubi.append(1656)
        return Dist((b*c)**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)*(-4*a*c + b**2)), Int((a + b*x**n + c*x**(2*n))**(p + 1)*ExpandToSum(Q*a*n*(p + 1)*(-4*a*c + b**2) + Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i))*(i + n*(2*p + 3) + 1) + x**i*(-a*b*(i + 1)*Coeff(R, x, i + n) + (-2*a*c*(i + 2*n*(p + 1) + 1) + b**2*(i + n*(p + 1) + 1))*Coeff(R, x, i)), List(i, 0, n - 1)), x), x), x) - Simp(x*(b*c)**(-Floor((q - 1)/n) - 1)*(a + b*x**n + c*x**(2*n))**(p + 1)*Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-2*a*c + b**2)*Coeff(R, x, i)), List(i, 0, n - 1))/(a*n*(p + 1)*(-4*a*c + b**2)), x)
    rule1656 = ReplacementRule(pattern1656, replacement1656)

    if load_rule:
        matcher.add(pattern1656, 1656)

    pattern1657 = Pattern(Integral(Pq_*x_**m_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons46, cons64, cons226, cons148, cons13, cons137, cons84, cons_with_1657)
    def replacement1657(Pq, p, c, n2, x, b, a, m, n):

        q = Expon(Pq, x)
        Q = PolynomialQuotient(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        R = PolynomialRemainder(Pq*a*x**m*(b*c)**(Floor((q + S(-1))/n) + S(1)), a + b*x**n + c*x**(S(2)*n), x)
        rubi.append(1657)
        return Dist((b*c)**(-Floor((q - 1)/n) - 1)/(a*n*(p + 1)*(-4*a*c + b**2)), Int(x**m*(a + b*x**n + c*x**(2*n))**(p + 1)*ExpandToSum(Q*n*x**(-m)*(p + 1)*(-4*a*c + b**2) + Sum_doit(c*x**(i - m + n)*(-2*Coeff(R, x, i + n) + b*Coeff(R, x, i)/a)*(i + n*(2*p + 3) + 1) + x**(i - m)*(-b*(i + 1)*Coeff(R, x, i + n) + (-2*c*(i + 2*n*(p + 1) + 1) + b**2*(i + n*(p + 1) + 1)/a)*Coeff(R, x, i)), List(i, 0, n - 1)), x), x), x) - Simp(x*(b*c)**(-Floor((q - 1)/n) - 1)*(a + b*x**n + c*x**(2*n))**(p + 1)*Sum_doit(c*x**(i + n)*(-2*a*Coeff(R, x, i + n) + b*Coeff(R, x, i)) + x**i*(-a*b*Coeff(R, x, i + n) + (-2*a*c + b**2)*Coeff(R, x, i)), List(i, 0, n - 1))/(a**2*n*(p + 1)*(-4*a*c + b**2)), x)
    rule1657 = ReplacementRule(pattern1657, replacement1657)

    if load_rule:
        matcher.add(pattern1657, 1657)

    pattern1658 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons858, cons226, cons148, cons17, cons_with_1658)
    def replacement1658(Pq, p, c, n2, x, b, a, m, n):

        g = GCD(m + S(1), n)
        rubi.append(1658)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a + b*x**(n/g) + c*x**(S(2)*n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)
    rule1658 = ReplacementRule(pattern1658, replacement1658)

    if load_rule:
        matcher.add(pattern1658, 1658)
    pattern1659 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons46, cons858, cons226, cons148, cons282)
    def replacement1659(d, Pq, c, n2, x, b, a, m, n):
        rubi.append(1659)
        return Int(ExpandIntegrand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1659 = ReplacementRule(pattern1659, replacement1659)

    if load_rule:
        matcher.add(pattern1659, 1659)
    pattern1660 = Pattern(Integral(Pq_/(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons7, cons46, cons858, cons226, cons148, cons947)
    def replacement1660(Pq, c, n2, x, b, a, n):
        rubi.append(1660)
        return Int(ExpandIntegrand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1660 = ReplacementRule(pattern1660, replacement1660)

    if load_rule:
        matcher.add(pattern1660, 1660)

    pattern1661 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons64, cons226, cons63, cons_with_1661)
    def replacement1661(Pq, p, c, x, b, a):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1661)
        return Dist(1/2, Int((a + b*x + c*x**2)**p*ExpandToSum(2*Pq - Pqq*c**p*(b + 2*c*x)*(a + b*x + c*x**2)**(-p - 1), x), x), x) + Simp(Pqq*c**p*log(a + b*x + c*x**2)/2, x)
    rule1661 = ReplacementRule(pattern1661, replacement1661)

    if load_rule:
        matcher.add(pattern1661, 1661)

    pattern1662 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons64, cons226, cons719, cons948, cons_with_1662)
    def replacement1662(Pq, p, c, x, b, a):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1662)
        return Int((a + b*x + c*x**2)**p*ExpandToSum(Pq - Pqq*c**(p + 1/2)*(a + b*x + c*x**2)**(-p - 1/2), x), x) + Simp(Pqq*c**p*atanh((b + 2*c*x)/(2*sqrt(a + b*x + c*x**2)*Rt(c, 2))), x)
    rule1662 = ReplacementRule(pattern1662, replacement1662)

    if load_rule:
        matcher.add(pattern1662, 1662)

    pattern1663 = Pattern(Integral(Pq_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons64, cons226, cons719, cons949, cons_with_1663)
    def replacement1663(Pq, p, c, x, b, a):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1663)
        return Int((a + b*x + c*x**2)**p*ExpandToSum(Pq - Pqq*(-c)**(p + 1/2)*(a + b*x + c*x**2)**(-p - 1/2), x), x) - Simp(Pqq*(-c)**p*ArcTan((b + 2*c*x)/(2*sqrt(a + b*x + c*x**2)*Rt(-c, 2))), x)
    rule1663 = ReplacementRule(pattern1663, replacement1663)

    if load_rule:
        matcher.add(pattern1663, 1663)

    pattern1664 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons46, cons858, cons226, cons148, cons_with_1664)
    def replacement1664(d, Pq, p, c, n2, x, b, a, m, n):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1664)
        return Int((d*x)**m*(a + b*x**n + c*x**(2*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-2*n + q)*(m - 2*n + q + 1) + b*x**(-n + q)*(m + n*(p - 1) + q + 1))/(c*(m + 2*n*p + q + 1)), x), x) + Simp(Pqq*d**(2*n - q - 1)*(d*x)**(m - 2*n + q + 1)*(a + b*x**n + c*x**(2*n))**(p + 1)/(c*(m + 2*n*p + q + 1)), x)
    rule1664 = ReplacementRule(pattern1664, replacement1664)

    if load_rule:
        matcher.add(pattern1664, 1664)

    pattern1665 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons858, cons226, cons148, cons_with_1665)
    def replacement1665(Pq, p, c, n2, x, b, a, n):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1665)
        return Int((a + b*x**n + c*x**(2*n))**p*ExpandToSum(Pq - Pqq*x**q - Pqq*(a*x**(-2*n + q)*(-2*n + q + 1) + b*x**(-n + q)*(n*(p - 1) + q + 1))/(c*(2*n*p + q + 1)), x), x) + Simp(Pqq*x**(-2*n + q + 1)*(a + b*x**n + c*x**(2*n))**(p + 1)/(c*(2*n*p + q + 1)), x)
    rule1665 = ReplacementRule(pattern1665, replacement1665)

    if load_rule:
        matcher.add(pattern1665, 1665)
    def With1666(d, Pq, p, c, n2, x, b, a, m, n):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        rubi.append(1666)
        return Int(Sum_doit(d**(-j)*(d*x)**(j + m)*(a + b*x**n + c*x**(S(2)*n))**p*Sum_doit(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)
    pattern1666 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons46, cons64, cons226, cons148, cons950)
    rule1666 = ReplacementRule(pattern1666, With1666)

    if load_rule:
        matcher.add(pattern1666, 1666)
    def With1667(Pq, p, c, n2, x, b, a, n):
        q = Expon(Pq, x)
        j = Symbol('j')
        k = Symbol('k')
        rubi.append(1667)
        return Int(Sum_doit(x**j*(a + b*x**n + c*x**(S(2)*n))**p*Sum_doit(x**(k*n)*Coeff(Pq, x, j + k*n), List(k, S(0), S(1) + (-j + q)/n)), List(j, S(0), n + S(-1))), x)
    pattern1667 = Pattern(Integral(Pq_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons64, cons226, cons148, cons950)
    rule1667 = ReplacementRule(pattern1667, With1667)

    if load_rule:
        matcher.add(pattern1667, 1667)
    pattern1668 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons46, cons64, cons226, cons148)
    def replacement1668(d, Pq, c, n2, x, b, a, m, n):
        rubi.append(1668)
        return Int(RationalFunctionExpand(Pq*(d*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1668 = ReplacementRule(pattern1668, replacement1668)

    if load_rule:
        matcher.add(pattern1668, 1668)
    pattern1669 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons46, cons64, cons226, cons148)
    def replacement1669(Pq, c, n2, x, b, a, n):
        rubi.append(1669)
        return Int(RationalFunctionExpand(Pq/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1669 = ReplacementRule(pattern1669, replacement1669)

    if load_rule:
        matcher.add(pattern1669, 1669)
    def With1670(Pq, p, n2, x, b, c, a, m, n):
        q = Expon(Pq, x)
        rubi.append(1670)
        return -Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x)
    pattern1670 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons64, cons226, cons196, cons17)
    rule1670 = ReplacementRule(pattern1670, With1670)

    if load_rule:
        matcher.add(pattern1670, 1670)
    def With1671(d, Pq, p, n2, x, b, c, a, m, n):
        g = Denominator(m)
        q = Expon(Pq, x)
        rubi.append(1671)
        return -Dist(g/d, Subst(Int(x**(-g*(m + q + S(1)) + S(-1))*(a + b*d**(-n)*x**(-g*n) + c*d**(-S(2)*n)*x**(-S(2)*g*n))**p*ExpandToSum(x**(g*q)*ReplaceAll(Pq, Rule(x, x**(-g)/d)), x), x), x, (d*x)**(-S(1)/g)), x)
    pattern1671 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons46, cons64, cons226, cons196, cons367)
    rule1671 = ReplacementRule(pattern1671, With1671)

    if load_rule:
        matcher.add(pattern1671, 1671)
    def With1672(d, Pq, p, n2, x, b, c, a, m, n):
        q = Expon(Pq, x)
        rubi.append(1672)
        return -Dist((d*x)**m*(S(1)/x)**m, Subst(Int(x**(-m - q + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p*ExpandToSum(x**q*ReplaceAll(Pq, Rule(x, S(1)/x)), x), x), x, S(1)/x), x)
    pattern1672 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons46, cons64, cons226, cons196, cons356)
    rule1672 = ReplacementRule(pattern1672, With1672)

    if load_rule:
        matcher.add(pattern1672, 1672)
    def With1673(Pq, p, n2, x, b, c, a, m, n):
        g = Denominator(n)
        rubi.append(1673)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
    pattern1673 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons46, cons64, cons226, cons489)
    rule1673 = ReplacementRule(pattern1673, With1673)

    if load_rule:
        matcher.add(pattern1673, 1673)
    def With1674(Pq, p, n2, x, b, c, a, n):
        g = Denominator(n)
        rubi.append(1674)
        return Dist(g, Subst(Int(x**(g + S(-1))*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p*ReplaceAll(Pq, Rule(x, x**g)), x), x, x**(S(1)/g)), x)
    pattern1674 = Pattern(Integral(Pq_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons64, cons226, cons489)
    rule1674 = ReplacementRule(pattern1674, With1674)

    if load_rule:
        matcher.add(pattern1674, 1674)
    pattern1675 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons46, cons64, cons226, cons489, cons73)
    def replacement1675(d, Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1675)
        return Dist(d**(m + S(-1)/2)*sqrt(d*x)/sqrt(x), Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1675 = ReplacementRule(pattern1675, replacement1675)

    if load_rule:
        matcher.add(pattern1675, 1675)
    pattern1676 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons46, cons64, cons226, cons489, cons951)
    def replacement1676(d, Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1676)
        return Dist(d**(m + S(1)/2)*sqrt(x)/sqrt(d*x), Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1676 = ReplacementRule(pattern1676, replacement1676)

    if load_rule:
        matcher.add(pattern1676, 1676)
    pattern1677 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons46, cons64, cons226, cons489)
    def replacement1677(d, Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1677)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1677 = ReplacementRule(pattern1677, replacement1677)

    if load_rule:
        matcher.add(pattern1677, 1677)
    pattern1678 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons858, cons226, cons541, cons23)
    def replacement1678(Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1678)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
    rule1678 = ReplacementRule(pattern1678, replacement1678)

    if load_rule:
        matcher.add(pattern1678, 1678)
    pattern1679 = Pattern(Integral(Pq_*(d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons46, cons858, cons226, cons541, cons23)
    def replacement1679(d, Pq, p, n2, x, b, c, a, m, n):
        rubi.append(1679)
        return Dist(x**(-m)*(d*x)**m, Int(Pq*x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1679 = ReplacementRule(pattern1679, replacement1679)

    if load_rule:
        matcher.add(pattern1679, 1679)
    def With1680(d, Pq, c, n2, x, b, a, m, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1680)
        return Dist(S(2)*c/q, Int(Pq*(d*x)**m/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int(Pq*(d*x)**m/(b + S(2)*c*x**n + q), x), x)
    pattern1680 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons46, cons64, cons226)
    rule1680 = ReplacementRule(pattern1680, With1680)

    if load_rule:
        matcher.add(pattern1680, 1680)
    def With1681(Pq, c, n2, x, b, a, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1681)
        return Dist(S(2)*c/q, Int(Pq/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int(Pq/(b + S(2)*c*x**n + q), x), x)
    pattern1681 = Pattern(Integral(Pq_/(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons4, cons46, cons64, cons226)
    rule1681 = ReplacementRule(pattern1681, With1681)

    if load_rule:
        matcher.add(pattern1681, 1681)
    pattern1682 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons46, cons64, cons702)
    def replacement1682(d, Pq, p, c, n2, x, b, a, m, n):
        rubi.append(1682)
        return Int(ExpandIntegrand(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1682 = ReplacementRule(pattern1682, replacement1682)

    if load_rule:
        matcher.add(pattern1682, 1682)
    pattern1683 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons64, cons702)
    def replacement1683(Pq, p, c, n2, x, b, a, n):
        rubi.append(1683)
        return Int(ExpandIntegrand(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1683 = ReplacementRule(pattern1683, replacement1683)

    if load_rule:
        matcher.add(pattern1683, 1683)
    pattern1684 = Pattern(Integral(Pq_*(x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons46, cons918)
    def replacement1684(d, Pq, p, c, n2, x, b, a, m, n):
        rubi.append(1684)
        return Int(Pq*(d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1684 = ReplacementRule(pattern1684, replacement1684)

    if load_rule:
        matcher.add(pattern1684, 1684)
    pattern1685 = Pattern(Integral(Pq_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons918)
    def replacement1685(Pq, p, c, n2, x, b, a, n):
        rubi.append(1685)
        return Int(Pq*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1685 = ReplacementRule(pattern1685, replacement1685)

    if load_rule:
        matcher.add(pattern1685, 1685)
    pattern1686 = Pattern(Integral(Pq_*u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons554, cons919)
    def replacement1686(Pq, p, v, u, n2, x, b, c, a, m, n):
        rubi.append(1686)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v), x)
    rule1686 = ReplacementRule(pattern1686, replacement1686)

    if load_rule:
        matcher.add(pattern1686, 1686)
    pattern1687 = Pattern(Integral(Pq_*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons552, cons919)
    def replacement1687(Pq, p, v, n2, x, b, c, a, n):
        rubi.append(1687)
        return Dist(S(1)/Coefficient(v, x, S(1)), Subst(Int((a + b*x**n + c*x**(S(2)*n))**p*SubstFor(v, Pq, x), x), x, v), x)
    rule1687 = ReplacementRule(pattern1687, replacement1687)

    if load_rule:
        matcher.add(pattern1687, 1687)
    pattern1688 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons147, cons952, cons953)
    def replacement1688(p, j, x, b, a, n):
        rubi.append(1688)
        return Simp(x**(-n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
    rule1688 = ReplacementRule(pattern1688, replacement1688)

    if load_rule:
        matcher.add(pattern1688, 1688)
    pattern1689 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons147, cons952, cons954, cons13, cons137)
    def replacement1689(p, j, x, b, a, n):
        rubi.append(1689)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1689 = ReplacementRule(pattern1689, replacement1689)

    if load_rule:
        matcher.add(pattern1689, 1689)
    pattern1690 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons147, cons952, cons954, cons955)
    def replacement1690(p, j, x, b, a, n):
        rubi.append(1690)
        return -Dist(b*(-j + n*p + n + S(1))/(a*(j*p + S(1))), Int(x**(-j + n)*(a*x**j + b*x**n)**p, x), x) + Simp(x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + S(1))), x)
    rule1690 = ReplacementRule(pattern1690, replacement1690)

    if load_rule:
        matcher.add(pattern1690, 1690)
    pattern1691 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons147, cons956, cons957, cons163, cons958)
    def replacement1691(p, j, x, b, a, n):
        rubi.append(1691)
        return -Dist(b*p*(-j + n)/(j*p + S(1)), Int(x**n*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(j*p + S(1)), x)
    rule1691 = ReplacementRule(pattern1691, replacement1691)

    if load_rule:
        matcher.add(pattern1691, 1691)
    pattern1692 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons147, cons956, cons957, cons163, cons959)
    def replacement1692(p, j, x, b, a, n):
        rubi.append(1692)
        return Dist(a*p*(-j + n)/(n*p + S(1)), Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(n*p + S(1)), x)
    rule1692 = ReplacementRule(pattern1692, replacement1692)

    if load_rule:
        matcher.add(pattern1692, 1692)
    pattern1693 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons147, cons956, cons957, cons137, cons960)
    def replacement1693(p, j, x, b, a, n):
        rubi.append(1693)
        return -Dist((j*p + j - n + S(1))/(b*(-j + n)*(p + S(1))), Int(x**(-n)*(a*x**j + b*x**n)**(p + S(1)), x), x) + Simp(x**(-n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
    rule1693 = ReplacementRule(pattern1693, replacement1693)

    if load_rule:
        matcher.add(pattern1693, 1693)
    pattern1694 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons147, cons956, cons957, cons137)
    def replacement1694(p, j, x, b, a, n):
        rubi.append(1694)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1694 = ReplacementRule(pattern1694, replacement1694)

    if load_rule:
        matcher.add(pattern1694, 1694)
    pattern1695 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons961, cons952, cons962)
    def replacement1695(p, j, x, b, a, n):
        rubi.append(1695)
        return Dist(a, Int(x**j*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp(x*(a*x**j + b*x**n)**p/(p*(-j + n)), x)
    rule1695 = ReplacementRule(pattern1695, replacement1695)

    if load_rule:
        matcher.add(pattern1695, 1695)
    pattern1696 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons4, cons963)
    def replacement1696(b, a, n, x):
        rubi.append(1696)
        return Dist(S(2)/(-n + S(2)), Subst(Int(S(1)/(-a*x**S(2) + S(1)), x), x, x/sqrt(a*x**S(2) + b*x**n)), x)
    rule1696 = ReplacementRule(pattern1696, replacement1696)

    if load_rule:
        matcher.add(pattern1696, 1696)
    pattern1697 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons719, cons952, cons962)
    def replacement1697(p, j, x, b, a, n):
        rubi.append(1697)
        return Dist((-j + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int(x**(-j)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(x**(-j + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1697 = ReplacementRule(pattern1697, replacement1697)

    if load_rule:
        matcher.add(pattern1697, 1697)
    pattern1698 = Pattern(Integral(S(1)/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons964, cons965)
    def replacement1698(j, x, b, a, n):
        rubi.append(1698)
        return -Dist(a*(-j + S(2)*n + S(-2))/(b*(n + S(-2))), Int(x**(j - n)/sqrt(a*x**j + b*x**n), x), x) + Simp(-S(2)*x**(-n + S(1))*sqrt(a*x**j + b*x**n)/(b*(n + S(-2))), x)
    rule1698 = ReplacementRule(pattern1698, replacement1698)

    if load_rule:
        matcher.add(pattern1698, 1698)
    pattern1699 = Pattern(Integral((x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons147, cons952, cons966)
    def replacement1699(p, j, x, b, a, n):
        rubi.append(1699)
        return Dist(x**(-j*FracPart(p))*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p), Int(x**(j*p)*(a + b*x**(-j + n))**p, x), x)
    rule1699 = ReplacementRule(pattern1699, replacement1699)

    if load_rule:
        matcher.add(pattern1699, 1699)
    pattern1700 = Pattern(Integral((u_**WC('j', S(1))*WC('a', S(1)) + u_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons68, cons69)
    def replacement1700(p, j, u, x, b, a, n):
        rubi.append(1700)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a*x**j + b*x**n)**p, x), x, u), x)
    rule1700 = ReplacementRule(pattern1700, replacement1700)

    if load_rule:
        matcher.add(pattern1700, 1700)
    pattern1701 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons21, cons4, cons5, cons147, cons952, cons967, cons53)
    def replacement1701(p, j, m, x, b, a, n):
        rubi.append(1701)
        return Dist(S(1)/n, Subst(Int((a*x**(j/n) + b*x)**p, x), x, x**n), x)
    rule1701 = ReplacementRule(pattern1701, replacement1701)

    if load_rule:
        matcher.add(pattern1701, 1701)
    pattern1702 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons968, cons969)
    def replacement1702(p, j, m, x, b, c, a, n):
        rubi.append(1702)
        return -Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1702 = ReplacementRule(pattern1702, replacement1702)

    if load_rule:
        matcher.add(pattern1702, 1702)
    pattern1703 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons147, cons952, cons970, cons13, cons137, cons969)
    def replacement1703(p, j, m, x, b, c, a, n):
        rubi.append(1703)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1703 = ReplacementRule(pattern1703, replacement1703)

    if load_rule:
        matcher.add(pattern1703, 1703)
    pattern1704 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons970, cons971, cons972)
    def replacement1704(p, j, m, x, b, c, a, n):
        rubi.append(1704)
        return -Dist(b*c**(j - n)*(-j + m + n*p + n + S(1))/(a*(j*p + m + S(1))), Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))), x)
    rule1704 = ReplacementRule(pattern1704, replacement1704)

    if load_rule:
        matcher.add(pattern1704, 1704)
    pattern1705 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons970)
    def replacement1705(p, j, c, m, x, b, a, n):
        rubi.append(1705)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
    rule1705 = ReplacementRule(pattern1705, replacement1705)

    if load_rule:
        matcher.add(pattern1705, 1705)
    pattern1706 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons21, cons4, cons5, cons147, cons952, cons967, cons500, cons973)
    def replacement1706(p, j, m, x, b, a, n):
        rubi.append(1706)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p, x), x, x**n), x)
    rule1706 = ReplacementRule(pattern1706, replacement1706)

    if load_rule:
        matcher.add(pattern1706, 1706)
    pattern1707 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons967, cons500, cons973)
    def replacement1707(p, j, c, m, x, b, a, n):
        rubi.append(1707)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
    rule1707 = ReplacementRule(pattern1707, replacement1707)

    if load_rule:
        matcher.add(pattern1707, 1707)
    pattern1708 = Pattern(Integral((x_*WC('c', S(1)))**m_*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons147, cons974, cons957, cons972, cons163, cons975)
    def replacement1708(p, j, m, x, b, c, a, n):
        rubi.append(1708)
        return -Dist(b*c**(-n)*p*(-j + n)/(j*p + m + S(1)), Int((c*x)**(m + n)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(j*p + m + S(1))), x)
    rule1708 = ReplacementRule(pattern1708, replacement1708)

    if load_rule:
        matcher.add(pattern1708, 1708)
    pattern1709 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons147, cons956, cons957, cons972, cons163, cons512)
    def replacement1709(p, j, m, x, b, c, a, n):
        rubi.append(1709)
        return Dist(a*c**(-j)*p*(-j + n)/(m + n*p + S(1)), Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*(m + n*p + S(1))), x)
    rule1709 = ReplacementRule(pattern1709, replacement1709)

    if load_rule:
        matcher.add(pattern1709, 1709)
    pattern1710 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons147, cons974, cons957, cons972, cons137, cons976)
    def replacement1710(p, j, m, x, b, c, a, n):
        rubi.append(1710)
        return -Dist(c**n*(j*p + j + m - n + S(1))/(b*(-j + n)*(p + S(1))), Int((c*x)**(m - n)*(a*x**j + b*x**n)**(p + S(1)), x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(-j + n)*(p + S(1))), x)
    rule1710 = ReplacementRule(pattern1710, replacement1710)

    if load_rule:
        matcher.add(pattern1710, 1710)
    pattern1711 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons147, cons956, cons957, cons972, cons137)
    def replacement1711(p, j, m, x, b, c, a, n):
        rubi.append(1711)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1711 = ReplacementRule(pattern1711, replacement1711)

    if load_rule:
        matcher.add(pattern1711, 1711)
    pattern1712 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons147, cons964, cons957, cons972, cons977, cons512)
    def replacement1712(p, j, m, x, b, c, a, n):
        rubi.append(1712)
        return -Dist(a*c**(-j + n)*(j*p + j + m - n + S(1))/(b*(m + n*p + S(1))), Int((c*x)**(j + m - n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(n + S(-1))*(c*x)**(m - n + S(1))*(a*x**j + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))), x)
    rule1712 = ReplacementRule(pattern1712, replacement1712)

    if load_rule:
        matcher.add(pattern1712, 1712)
    pattern1713 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons147, cons964, cons957, cons972, cons978)
    def replacement1713(p, j, m, x, b, c, a, n):
        rubi.append(1713)
        return -Dist(b*c**(j - n)*(-j + m + n*p + n + S(1))/(a*(j*p + m + S(1))), Int((c*x)**(-j + m + n)*(a*x**j + b*x**n)**p, x), x) + Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(j*p + m + S(1))), x)
    rule1713 = ReplacementRule(pattern1713, replacement1713)

    if load_rule:
        matcher.add(pattern1713, 1713)
    pattern1714 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons21, cons4, cons5, cons147, cons952, cons967, cons66, cons541, cons23)
    def replacement1714(p, j, m, x, b, a, n):
        rubi.append(1714)
        return Dist(S(1)/(m + S(1)), Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1))), x)
    rule1714 = ReplacementRule(pattern1714, replacement1714)

    if load_rule:
        matcher.add(pattern1714, 1714)
    pattern1715 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons967, cons66, cons541, cons23)
    def replacement1715(p, j, c, m, x, b, a, n):
        rubi.append(1715)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
    rule1715 = ReplacementRule(pattern1715, replacement1715)

    if load_rule:
        matcher.add(pattern1715, 1715)
    pattern1716 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons961, cons952, cons979, cons969)
    def replacement1716(p, j, m, x, b, c, a, n):
        rubi.append(1716)
        return Dist(a*c**(-j), Int((c*x)**(j + m)*(a*x**j + b*x**n)**(p + S(-1)), x), x) + Simp((c*x)**(m + S(1))*(a*x**j + b*x**n)**p/(c*p*(-j + n)), x)
    rule1716 = ReplacementRule(pattern1716, replacement1716)

    if load_rule:
        matcher.add(pattern1716, 1716)
    pattern1717 = Pattern(Integral(x_**WC('m', S(1))/sqrt(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons796, cons4, cons980, cons952)
    def replacement1717(j, m, x, b, a, n):
        rubi.append(1717)
        return Dist(-S(2)/(-j + n), Subst(Int(S(1)/(-a*x**S(2) + S(1)), x), x, x**(j/S(2))/sqrt(a*x**j + b*x**n)), x)
    rule1717 = ReplacementRule(pattern1717, replacement1717)

    if load_rule:
        matcher.add(pattern1717, 1717)
    pattern1718 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons719, cons952, cons979, cons969)
    def replacement1718(p, j, m, x, b, c, a, n):
        rubi.append(1718)
        return Dist(c**j*(-j + m + n*p + n + S(1))/(a*(-j + n)*(p + S(1))), Int((c*x)**(-j + m)*(a*x**j + b*x**n)**(p + S(1)), x), x) - Simp(c**(j + S(-1))*(c*x)**(-j + m + S(1))*(a*x**j + b*x**n)**(p + S(1))/(a*(-j + n)*(p + S(1))), x)
    rule1718 = ReplacementRule(pattern1718, replacement1718)

    if load_rule:
        matcher.add(pattern1718, 1718)
    pattern1719 = Pattern(Integral((c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons667, cons952, cons979)
    def replacement1719(p, j, c, m, x, b, a, n):
        rubi.append(1719)
        return Dist(c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m), Int(x**m*(a*x**j + b*x**n)**p, x), x)
    rule1719 = ReplacementRule(pattern1719, replacement1719)

    if load_rule:
        matcher.add(pattern1719, 1719)
    pattern1720 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons147, cons952, cons966)
    def replacement1720(p, j, m, x, b, c, a, n):
        rubi.append(1720)
        return Dist(c**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(c*x)**FracPart(m)*(a + b*x**(-j + n))**(-FracPart(p))*(a*x**j + b*x**n)**FracPart(p), Int(x**(j*p + m)*(a + b*x**(-j + n))**p, x), x)
    rule1720 = ReplacementRule(pattern1720, replacement1720)

    if load_rule:
        matcher.add(pattern1720, 1720)
    pattern1721 = Pattern(Integral(u_**WC('m', S(1))*(v_**WC('j', S(1))*WC('a', S(1)) + v_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons796, cons21, cons4, cons5, cons554)
    def replacement1721(v, p, j, u, m, x, b, a, n):
        rubi.append(1721)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a*x**j + b*x**n)**p, x), x, v), x)
    rule1721 = ReplacementRule(pattern1721, replacement1721)

    if load_rule:
        matcher.add(pattern1721, 1721)
    pattern1722 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons796, cons797, cons21, cons4, cons5, cons50, cons147, cons981, cons967, cons982, cons500, cons973)
    def replacement1722(d, q, p, k, j, c, m, x, b, a, n):
        rubi.append(1722)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(c + d*x)**q*(a*x**(j/n) + b*x**(k/n))**p, x), x, x**n), x)
    rule1722 = ReplacementRule(pattern1722, replacement1722)

    if load_rule:
        matcher.add(pattern1722, 1722)
    pattern1723 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons797, cons21, cons4, cons5, cons50, cons147, cons981, cons967, cons982, cons500, cons973)
    def replacement1723(d, q, p, k, c, j, e, x, b, a, m, n):
        rubi.append(1723)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x), x)
    rule1723 = ReplacementRule(pattern1723, replacement1723)

    if load_rule:
        matcher.add(pattern1723, 1723)
    pattern1724 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons21, cons4, cons5, cons983, cons147, cons71, cons984, cons985, cons971)
    def replacement1724(d, p, j, c, e, jn, x, b, a, m, n):
        rubi.append(1724)
        return Simp(c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))), x)
    rule1724 = ReplacementRule(pattern1724, replacement1724)

    if load_rule:
        matcher.add(pattern1724, 1724)
    pattern1725 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons21, cons4, cons983, cons147, cons71, cons986, cons137, cons987, cons988)
    def replacement1725(d, p, j, c, e, jn, x, b, a, m, n):
        rubi.append(1725)
        return -Dist(e**j*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(a*b*n*(p + S(1))), Int((e*x)**(-j + m)*(a*x**j + b*x**(j + n))**(p + S(1)), x), x) - Simp(e**(j + S(-1))*(e*x)**(-j + m + S(1))*(-a*d + b*c)*(a*x**j + b*x**(j + n))**(p + S(1))/(a*b*n*(p + S(1))), x)
    rule1725 = ReplacementRule(pattern1725, replacement1725)

    if load_rule:
        matcher.add(pattern1725, 1725)
    pattern1726 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons5, cons983, cons147, cons71, cons93, cons88, cons989, cons990, cons971, cons991)
    def replacement1726(d, p, j, c, e, jn, x, b, a, m, n):
        rubi.append(1726)
        return Dist(e**(-n)*(a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(a*(j*p + m + S(1))), Int((e*x)**(m + n)*(a*x**j + b*x**(j + n))**p, x), x) + Simp(c*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(a*(j*p + m + S(1))), x)
    rule1726 = ReplacementRule(pattern1726, replacement1726)

    if load_rule:
        matcher.add(pattern1726, 1726)
    pattern1727 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons21, cons4, cons5, cons983, cons147, cons71, cons992, cons988)
    def replacement1727(d, p, j, c, e, jn, x, b, a, m, n):
        rubi.append(1727)
        return -Dist((a*d*(j*p + m + S(1)) - b*c*(m + n + p*(j + n) + S(1)))/(b*(m + n + p*(j + n) + S(1))), Int((e*x)**m*(a*x**j + b*x**(j + n))**p, x), x) + Simp(d*e**(j + S(-1))*(e*x)**(-j + m + S(1))*(a*x**j + b*x**(j + n))**(p + S(1))/(b*(m + n + p*(j + n) + S(1))), x)
    rule1727 = ReplacementRule(pattern1727, replacement1727)

    if load_rule:
        matcher.add(pattern1727, 1727)
    pattern1728 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons796, cons797, cons21, cons4, cons5, cons50, cons147, cons981, cons967, cons982, cons66, cons541, cons23)
    def replacement1728(d, q, p, k, c, j, x, b, a, m, n):
        rubi.append(1728)
        return Dist(S(1)/(m + S(1)), Subst(Int((c + d*x**(n/(m + S(1))))**q*(a*x**(j/(m + S(1))) + b*x**(k/(m + S(1))))**p, x), x, x**(m + S(1))), x)
    rule1728 = ReplacementRule(pattern1728, replacement1728)

    if load_rule:
        matcher.add(pattern1728, 1728)
    pattern1729 = Pattern(Integral((e_*x_)**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**j_*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons797, cons21, cons4, cons5, cons50, cons147, cons981, cons967, cons982, cons66, cons541, cons23)
    def replacement1729(d, q, p, k, c, j, e, x, b, a, m, n):
        rubi.append(1729)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(c + d*x**n)**q*(a*x**j + b*x**k)**p, x), x)
    rule1729 = ReplacementRule(pattern1729, replacement1729)

    if load_rule:
        matcher.add(pattern1729, 1729)
    pattern1730 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('jn', S(1))*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons796, cons21, cons4, cons5, cons50, cons983, cons147, cons71, cons993)
    def replacement1730(d, q, p, j, c, e, jn, x, b, a, m, n):
        rubi.append(1730)
        return Dist(e**IntPart(m)*x**(-j*FracPart(p) - FracPart(m))*(e*x)**FracPart(m)*(a + b*x**n)**(-FracPart(p))*(a*x**j + b*x**(j + n))**FracPart(p), Int(x**(j*p + m)*(a + b*x**n)**p*(c + d*x**n)**q, x), x)
    rule1730 = ReplacementRule(pattern1730, replacement1730)

    if load_rule:
        matcher.add(pattern1730, 1730)
    def With1731(Pq, p, j, x, b, a, n):
        d = Denominator(n)
        rubi.append(1731)
        return Dist(d, Subst(Int(x**(d + S(-1))*(a*x**(d*j) + b*x**(d*n))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(d*n))), x), x, x**(S(1)/d)), x)
    pattern1731 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons858, cons147, cons952, cons964, cons967, cons994)
    rule1731 = ReplacementRule(pattern1731, With1731)

    if load_rule:
        matcher.add(pattern1731, 1731)
    pattern1732 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons21, cons4, cons5, cons858, cons147, cons952, cons967, cons500)
    def replacement1732(Pq, p, j, m, x, b, a, n):
        rubi.append(1732)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a*x**(j/n) + b*x)**p*SubstFor(x**n, Pq, x), x), x, x**n), x)
    rule1732 = ReplacementRule(pattern1732, replacement1732)

    if load_rule:
        matcher.add(pattern1732, 1732)
    pattern1733 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons4, cons5, cons858, cons147, cons952, cons967, cons500, cons31, cons995)
    def replacement1733(Pq, p, j, c, m, x, b, a, n):
        rubi.append(1733)
        return Dist(c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m)), Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
    rule1733 = ReplacementRule(pattern1733, replacement1733)

    if load_rule:
        matcher.add(pattern1733, 1733)
    pattern1734 = Pattern(Integral(Pq_*(c_*x_)**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons858, cons147, cons952, cons967, cons500)
    def replacement1734(Pq, p, j, c, m, x, b, a, n):
        rubi.append(1734)
        return Dist(x**(-m)*(c*x)**m, Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
    rule1734 = ReplacementRule(pattern1734, replacement1734)

    if load_rule:
        matcher.add(pattern1734, 1734)

    pattern1735 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons5, cons858, cons147, cons996, cons17, cons_with_1735)
    def replacement1735(Pq, p, j, m, x, b, a, n):

        g = GCD(m + S(1), n)
        rubi.append(1735)
        return Dist(S(1)/g, Subst(Int(x**(S(-1) + (m + S(1))/g)*(a*x**(j/g) + b*x**(n/g))**p*ReplaceAll(Pq, Rule(x, x**(S(1)/g))), x), x, x**g), x)
    rule1735 = ReplacementRule(pattern1735, replacement1735)

    if load_rule:
        matcher.add(pattern1735, 1735)

    pattern1736 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons64, cons147, cons997, cons998, cons_with_1736)
    def replacement1736(Pq, p, j, m, x, b, c, a, n):

        q = Expon(Pq, x)
        Pqq = Coeff(Pq, x, q)
        rubi.append(1736)
        return Int((c*x)**m*(a*x**j + b*x**n)**p*ExpandToSum(Pq - Pqq*a*x**(-n + q)*(m - n + q + 1)/(b*(m + n*p + q + 1)) - Pqq*x**q, x), x) + Simp(Pqq*c**(n - q - 1)*(c*x)**(m - n + q + 1)*(a*x**j + b*x**n)**(p + 1)/(b*(m + n*p + q + 1)), x)
    rule1736 = ReplacementRule(pattern1736, replacement1736)

    if load_rule:
        matcher.add(pattern1736, 1736)
    pattern1737 = Pattern(Integral(Pq_*x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons21, cons4, cons5, cons858, cons147, cons952, cons967, cons541, cons23)
    def replacement1737(Pq, p, j, m, x, b, a, n):
        rubi.append(1737)
        return Dist(S(1)/(m + S(1)), Subst(Int((a*x**(j/(m + S(1))) + b*x**(n/(m + S(1))))**p*ReplaceAll(SubstFor(x**n, Pq, x), Rule(x, x**(n/(m + S(1))))), x), x, x**(m + S(1))), x)
    rule1737 = ReplacementRule(pattern1737, replacement1737)

    if load_rule:
        matcher.add(pattern1737, 1737)
    pattern1738 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons4, cons5, cons858, cons147, cons952, cons967, cons541, cons23, cons31, cons995)
    def replacement1738(Pq, p, j, c, m, x, b, a, n):
        rubi.append(1738)
        return Dist(c**(Quotient(m, sign(m))*sign(m))*x**(-Mod(m, sign(m)))*(c*x)**Mod(m, sign(m)), Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
    rule1738 = ReplacementRule(pattern1738, replacement1738)

    if load_rule:
        matcher.add(pattern1738, 1738)
    pattern1739 = Pattern(Integral(Pq_*(c_*x_)**m_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons858, cons147, cons952, cons967, cons541, cons23)
    def replacement1739(Pq, p, j, c, m, x, b, a, n):
        rubi.append(1739)
        return Dist(x**(-m)*(c*x)**m, Int(Pq*x**m*(a*x**j + b*x**n)**p, x), x)
    rule1739 = ReplacementRule(pattern1739, replacement1739)

    if load_rule:
        matcher.add(pattern1739, 1739)
    pattern1740 = Pattern(Integral(Pq_*(x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons796, cons21, cons4, cons5, cons918, cons147, cons952)
    def replacement1740(Pq, p, j, m, x, b, c, a, n):
        rubi.append(1740)
        return Int(ExpandIntegrand(Pq*(c*x)**m*(a*x**j + b*x**n)**p, x), x)
    rule1740 = ReplacementRule(pattern1740, replacement1740)

    if load_rule:
        matcher.add(pattern1740, 1740)
    pattern1741 = Pattern(Integral(Pq_*(x_**n_*WC('b', S(1)) + x_**WC('j', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons796, cons4, cons5, cons918, cons147, cons952)
    def replacement1741(Pq, p, j, x, b, a, n):
        rubi.append(1741)
        return Int(ExpandIntegrand(Pq*(a*x**j + b*x**n)**p, x), x)
    rule1741 = ReplacementRule(pattern1741, replacement1741)

    if load_rule:
        matcher.add(pattern1741, 1741)
    pattern1742 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons38, cons999)
    def replacement1742(d, p, x, b, a):
        rubi.append(1742)
        return Dist(S(3)**(-S(3)*p)*a**(-S(2)*p), Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x), x)
    rule1742 = ReplacementRule(pattern1742, replacement1742)

    if load_rule:
        matcher.add(pattern1742, 1742)
    pattern1743 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons128, cons1000)
    def replacement1743(d, p, x, b, a):
        rubi.append(1743)
        return Int(ExpandToSum((a + b*x + d*x**S(3))**p, x), x)
    rule1743 = ReplacementRule(pattern1743, replacement1743)

    if load_rule:
        matcher.add(pattern1743, 1743)

    pattern1744 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons63, cons1000, cons_with_1744)
    def replacement1744(d, p, x, b, a):

        u = Factor(a + b*x + d*x**S(3))
        rubi.append(1744)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1744 = ReplacementRule(pattern1744, replacement1744)

    if load_rule:
        matcher.add(pattern1744, 1744)
    def With1745(d, p, x, b, a):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        rubi.append(1745)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
    pattern1745 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons63, cons1000)
    rule1745 = ReplacementRule(pattern1745, With1745)

    if load_rule:
        matcher.add(pattern1745, 1745)
    pattern1746 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons5, cons147, cons999)
    def replacement1746(d, p, x, b, a):
        rubi.append(1746)
        return Dist((S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p, Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p), x), x)
    rule1746 = ReplacementRule(pattern1746, replacement1746)

    if load_rule:
        matcher.add(pattern1746, 1746)

    pattern1747 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons5, cons147, cons1000, cons_with_1747)
    def replacement1747(d, p, x, b, a):

        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        rubi.append(1747)
        return Dist((a + b*x + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)
    rule1747 = ReplacementRule(pattern1747, replacement1747)

    if load_rule:
        matcher.add(pattern1747, 1747)
    def With1748(d, p, x, b, a):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        rubi.append(1748)
        return Dist((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p, Int((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
    pattern1748 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons5, cons147, cons1000)
    rule1748 = ReplacementRule(pattern1748, With1748)

    if load_rule:
        matcher.add(pattern1748, 1748)
    pattern1749 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons38, cons999)
    def replacement1749(d, p, f, e, x, b, a, m):
        rubi.append(1749)
        return Dist(S(3)**(-S(3)*p)*a**(-S(2)*p), Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x), x)
    rule1749 = ReplacementRule(pattern1749, replacement1749)

    if load_rule:
        matcher.add(pattern1749, 1749)
    pattern1750 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons128, cons1000)
    def replacement1750(d, p, f, e, x, b, a, m):
        rubi.append(1750)
        return Int(ExpandIntegrand((e + f*x)**m*(a + b*x + d*x**S(3))**p, x), x)
    rule1750 = ReplacementRule(pattern1750, replacement1750)

    if load_rule:
        matcher.add(pattern1750, 1750)

    pattern1751 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons63, cons1000, cons_with_1751)
    def replacement1751(d, p, f, e, x, b, a, m):

        u = Factor(a + b*x + d*x**S(3))
        rubi.append(1751)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1751 = ReplacementRule(pattern1751, replacement1751)

    if load_rule:
        matcher.add(pattern1751, 1751)
    def With1752(d, p, f, e, x, b, a, m):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        rubi.append(1752)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
    pattern1752 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons63, cons1000)
    rule1752 = ReplacementRule(pattern1752, With1752)

    if load_rule:
        matcher.add(pattern1752, 1752)
    pattern1753 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons5, cons147, cons999)
    def replacement1753(d, p, f, e, x, b, a, m):
        rubi.append(1753)
        return Dist((S(3)*a - b*x)**(-p)*(S(3)*a + S(2)*b*x)**(-S(2)*p)*(a + b*x + d*x**S(3))**p, Int((S(3)*a - b*x)**p*(S(3)*a + S(2)*b*x)**(S(2)*p)*(e + f*x)**m, x), x)
    rule1753 = ReplacementRule(pattern1753, replacement1753)

    if load_rule:
        matcher.add(pattern1753, 1753)

    pattern1754 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons5, cons147, cons1000, cons_with_1754)
    def replacement1754(d, p, f, e, x, b, a, m):

        u = NonfreeFactors(Factor(a + b*x + d*x**S(3)), x)
        rubi.append(1754)
        return Dist((a + b*x + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)
    rule1754 = ReplacementRule(pattern1754, replacement1754)

    if load_rule:
        matcher.add(pattern1754, 1754)
    def With1755(d, p, f, e, x, b, a, m):
        r = Rt(-S(27)*a*d**S(2) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*b**S(3)*d), S(3))
        rubi.append(1755)
        return Dist((-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p)*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(a + b*x + d*x**S(3))**p, Int((e + f*x)**m*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) - sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(-S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p*(S(3)*d*x + S(2)**(S(1)/3)*(S(6)*b*d - S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p, x), x)
    pattern1755 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons27, cons48, cons125, cons21, cons5, cons147, cons1000)
    rule1755 = ReplacementRule(pattern1755, With1755)

    if load_rule:
        matcher.add(pattern1755, 1755)
    pattern1756 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons38, cons1001)
    def replacement1756(d, p, x, c, a):
        rubi.append(1756)
        return -Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x), x)
    rule1756 = ReplacementRule(pattern1756, replacement1756)

    if load_rule:
        matcher.add(pattern1756, 1756)
    pattern1757 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons128, cons1002)
    def replacement1757(d, p, x, c, a):
        rubi.append(1757)
        return Int(ExpandToSum((a + c*x**S(2) + d*x**S(3))**p, x), x)
    rule1757 = ReplacementRule(pattern1757, replacement1757)

    if load_rule:
        matcher.add(pattern1757, 1757)

    pattern1758 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons63, cons1002, cons_with_1758)
    def replacement1758(d, p, x, c, a):

        u = Factor(a + c*x**S(2) + d*x**S(3))
        rubi.append(1758)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1758 = ReplacementRule(pattern1758, replacement1758)

    if load_rule:
        matcher.add(pattern1758, 1758)
    def With1759(d, p, x, c, a):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        rubi.append(1759)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
    pattern1759 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons63, cons1002)
    rule1759 = ReplacementRule(pattern1759, With1759)

    if load_rule:
        matcher.add(pattern1759, 1759)
    pattern1760 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons5, cons147, cons1001)
    def replacement1760(d, p, x, c, a):
        rubi.append(1760)
        return Dist((c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p, Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p), x), x)
    rule1760 = ReplacementRule(pattern1760, replacement1760)

    if load_rule:
        matcher.add(pattern1760, 1760)

    pattern1761 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons5, cons147, cons1002, cons_with_1761)
    def replacement1761(d, p, x, c, a):

        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        rubi.append(1761)
        return Dist((a + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)
    rule1761 = ReplacementRule(pattern1761, replacement1761)

    if load_rule:
        matcher.add(pattern1761, 1761)
    def With1762(d, p, x, c, a):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        rubi.append(1762)
        return Dist((a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p), Int((c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
    pattern1762 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons5, cons147, cons1002)
    rule1762 = ReplacementRule(pattern1762, With1762)

    if load_rule:
        matcher.add(pattern1762, 1762)
    pattern1763 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons38, cons1001)
    def replacement1763(d, p, f, e, x, c, a, m):
        rubi.append(1763)
        return -Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x), x)
    rule1763 = ReplacementRule(pattern1763, replacement1763)

    if load_rule:
        matcher.add(pattern1763, 1763)
    pattern1764 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons128, cons1002)
    def replacement1764(d, p, f, e, x, c, a, m):
        rubi.append(1764)
        return Int(ExpandIntegrand((e + f*x)**m*(a + c*x**S(2) + d*x**S(3))**p, x), x)
    rule1764 = ReplacementRule(pattern1764, replacement1764)

    if load_rule:
        matcher.add(pattern1764, 1764)

    pattern1765 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons63, cons1002, cons_with_1765)
    def replacement1765(d, p, f, e, x, c, a, m):

        u = Factor(a + c*x**S(2) + d*x**S(3))
        rubi.append(1765)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1765 = ReplacementRule(pattern1765, replacement1765)

    if load_rule:
        matcher.add(pattern1765, 1765)
    def With1766(d, p, f, e, x, c, a, m):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        rubi.append(1766)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
    pattern1766 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons63, cons1002)
    rule1766 = ReplacementRule(pattern1766, With1766)

    if load_rule:
        matcher.add(pattern1766, 1766)
    pattern1767 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1001)
    def replacement1767(d, p, f, e, x, c, a, m):
        rubi.append(1767)
        return Dist((c - S(3)*d*x)**(-p)*(S(2)*c + S(3)*d*x)**(-S(2)*p)*(a + c*x**S(2) + d*x**S(3))**p, Int((c - S(3)*d*x)**p*(S(2)*c + S(3)*d*x)**(S(2)*p)*(e + f*x)**m, x), x)
    rule1767 = ReplacementRule(pattern1767, replacement1767)

    if load_rule:
        matcher.add(pattern1767, 1767)

    pattern1768 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1002, cons_with_1768)
    def replacement1768(d, p, f, e, x, c, a, m):

        u = NonfreeFactors(Factor(a + c*x**S(2) + d*x**S(3)), x)
        rubi.append(1768)
        return Dist((a + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)
    rule1768 = ReplacementRule(pattern1768, replacement1768)

    if load_rule:
        matcher.add(pattern1768, 1768)
    def With1769(d, p, f, e, x, c, a, m):
        r = Rt(-S(27)*a*d**S(2) - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) + S(4)*a*c**S(3)), S(3))
        rubi.append(1769)
        return Dist((a + c*x**S(2) + d*x**S(3))**p*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**(-p), Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) + sqrt(S(3))*I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) + S(2)**(S(1)/3)*r**S(2)*(S(1) - sqrt(S(3))*I))/(S(4)*r))**p, x), x)
    pattern1769 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1002)
    rule1769 = ReplacementRule(pattern1769, With1769)

    if load_rule:
        matcher.add(pattern1769, 1769)
    pattern1770 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons38, cons1003, cons1004)
    def replacement1770(d, p, c, x, b, a):
        rubi.append(1770)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + c*x)**(S(3)*p), x), x)
    rule1770 = ReplacementRule(pattern1770, replacement1770)

    if load_rule:
        matcher.add(pattern1770, 1770)
    pattern1771 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons38, cons1003, cons1005)
    def replacement1771(d, p, c, x, b, a):
        rubi.append(1771)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Subst(Int((S(3)*a*b*c - b**S(3) + c**S(3)*x**S(3))**p, x), x, c/(S(3)*d) + x), x)
    rule1771 = ReplacementRule(pattern1771, replacement1771)

    if load_rule:
        matcher.add(pattern1771, 1771)
    def With1772(d, p, c, x, b, a):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        rubi.append(1772)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p, x), x)
    pattern1772 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons38, cons1006, cons1004)
    rule1772 = ReplacementRule(pattern1772, With1772)

    if load_rule:
        matcher.add(pattern1772, 1772)
    pattern1773 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons128, cons1006, cons1005)
    def replacement1773(d, p, c, x, b, a):
        rubi.append(1773)
        return Int(ExpandToSum((a + b*x + c*x**S(2) + d*x**S(3))**p, x), x)
    rule1773 = ReplacementRule(pattern1773, replacement1773)

    if load_rule:
        matcher.add(pattern1773, 1773)

    pattern1774 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons63, cons1006, cons1005, cons_with_1774)
    def replacement1774(d, p, c, x, b, a):

        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        rubi.append(1774)
        return Dist(FreeFactors(u, x)**p, Int(DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1774 = ReplacementRule(pattern1774, replacement1774)

    if load_rule:
        matcher.add(pattern1774, 1774)
    pattern1775 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons63, cons1006, cons1005)
    def replacement1775(d, p, c, x, b, a):
        rubi.append(1775)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x), x)
    rule1775 = ReplacementRule(pattern1775, replacement1775)

    if load_rule:
        matcher.add(pattern1775, 1775)
    pattern1776 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons147, cons1003, cons1004)
    def replacement1776(d, p, c, x, b, a):
        rubi.append(1776)
        return Dist((b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x)**(S(3)*p), x), x)
    rule1776 = ReplacementRule(pattern1776, replacement1776)

    if load_rule:
        matcher.add(pattern1776, 1776)
    def With1777(d, p, c, x, b, a):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        rubi.append(1777)
        return Dist((b + c*x - r)**(-p)*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**(-p)*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)
    pattern1777 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons147, cons1003, cons1005)
    rule1777 = ReplacementRule(pattern1777, With1777)

    if load_rule:
        matcher.add(pattern1777, 1777)
    def With1778(d, p, c, x, b, a):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        rubi.append(1778)
        return Dist((b + x*(c - r))**(-p)*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**(-p)*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p, x), x)
    pattern1778 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons147, cons1006, cons1004)
    rule1778 = ReplacementRule(pattern1778, With1778)

    if load_rule:
        matcher.add(pattern1778, 1778)

    pattern1779 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons147, cons1006, cons1005, cons_with_1779)
    def replacement1779(d, p, c, x, b, a):

        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        rubi.append(1779)
        return Dist((a + b*x + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int(DistributeDegree(u, p), x), x)
    rule1779 = ReplacementRule(pattern1779, replacement1779)

    if load_rule:
        matcher.add(pattern1779, 1779)
    def With1780(d, p, c, x, b, a):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        rubi.append(1780)
        return Dist((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**p, x), x)
    pattern1780 = Pattern(Integral((x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons147, cons1006, cons1005)
    rule1780 = ReplacementRule(pattern1780, With1780)

    if load_rule:
        matcher.add(pattern1780, 1780)
    pattern1781 = Pattern(Integral(u_**p_, x_), cons5, cons1007, cons1008)
    def replacement1781(x, u, p):
        rubi.append(1781)
        return Int(ExpandToSum(u, x)**p, x)
    rule1781 = ReplacementRule(pattern1781, replacement1781)

    if load_rule:
        matcher.add(pattern1781, 1781)
    pattern1782 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons38, cons1003, cons1004)
    def replacement1782(d, p, c, f, e, x, b, a, m):
        rubi.append(1782)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x), x)
    rule1782 = ReplacementRule(pattern1782, replacement1782)

    if load_rule:
        matcher.add(pattern1782, 1782)
    def With1783(d, p, c, f, e, x, b, a, m):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        rubi.append(1783)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)
    pattern1783 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons38, cons1003, cons1005)
    rule1783 = ReplacementRule(pattern1783, With1783)

    if load_rule:
        matcher.add(pattern1783, 1783)
    def With1784(d, p, c, f, e, x, b, a, m):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        rubi.append(1784)
        return Dist(S(3)**(-p)*b**(-p)*c**(-p), Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p*(e + f*x)**m, x), x)
    pattern1784 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons38, cons1006, cons1004)
    rule1784 = ReplacementRule(pattern1784, With1784)

    if load_rule:
        matcher.add(pattern1784, 1784)
    pattern1785 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons128, cons1006, cons1005)
    def replacement1785(d, p, c, f, e, x, b, a, m):
        rubi.append(1785)
        return Int(ExpandIntegrand((e + f*x)**m*(a + b*x + c*x**S(2) + d*x**S(3))**p, x), x)
    rule1785 = ReplacementRule(pattern1785, replacement1785)

    if load_rule:
        matcher.add(pattern1785, 1785)

    pattern1786 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons63, cons1006, cons1005, cons_with_1786)
    def replacement1786(d, p, c, f, e, x, b, a, m):

        u = Factor(a + b*x + c*x**S(2) + d*x**S(3))
        rubi.append(1786)
        return Dist(FreeFactors(u, x)**p, Int((e + f*x)**m*DistributeDegree(NonfreeFactors(u, x), p), x), x)
    rule1786 = ReplacementRule(pattern1786, replacement1786)

    if load_rule:
        matcher.add(pattern1786, 1786)
    pattern1787 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons63, cons1006, cons1005)
    def replacement1787(d, p, c, f, e, x, b, a, m):
        rubi.append(1787)
        return Dist(S(3)**(-S(3)*p)*d**(-S(2)*p), Subst(Int((S(27)*a*d**S(2) - S(9)*b*c*d + S(2)*c**S(3) + S(27)*d**S(3)*x**S(3) - S(9)*d*x*(-S(3)*b*d + c**S(2)))**p, x), x, c/(S(3)*d) + x), x)
    rule1787 = ReplacementRule(pattern1787, replacement1787)

    if load_rule:
        matcher.add(pattern1787, 1787)
    pattern1788 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1003, cons1004)
    def replacement1788(d, p, c, f, e, x, b, a, m):
        rubi.append(1788)
        return Dist((b + c*x)**(-S(3)*p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + c*x)**(S(3)*p)*(e + f*x)**m, x), x)
    rule1788 = ReplacementRule(pattern1788, replacement1788)

    if load_rule:
        matcher.add(pattern1788, 1788)
    def With1789(d, p, c, f, e, x, b, a, m):
        r = Rt(-S(3)*a*b*c + b**S(3), S(3))
        rubi.append(1789)
        return Dist((b + c*x - r)**(-p)*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**(-p)*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((e + f*x)**m*(b + c*x - r)**p*(b + c*x + r*(S(1) - sqrt(S(3))*I)/S(2))**p*(b + c*x + r*(S(1) + sqrt(S(3))*I)/S(2))**p, x), x)
    pattern1789 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1003, cons1005)
    rule1789 = ReplacementRule(pattern1789, With1789)

    if load_rule:
        matcher.add(pattern1789, 1789)
    def With1790(d, p, c, f, e, x, b, a, m):
        r = Rt(-S(3)*b*c*d + c**S(3), S(3))
        rubi.append(1790)
        return Dist((b + x*(c - r))**(-p)*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**(-p)*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((b + x*(c - r))**p*(b + x*(c + r*(S(1) - sqrt(S(3))*I)/S(2)))**p*(b + x*(c + r*(S(1) + sqrt(S(3))*I)/S(2)))**p*(e + f*x)**m, x), x)
    pattern1790 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1006, cons1004)
    rule1790 = ReplacementRule(pattern1790, With1790)

    if load_rule:
        matcher.add(pattern1790, 1790)

    pattern1791 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1006, cons1005, cons_with_1791)
    def replacement1791(d, p, c, f, e, x, b, a, m):

        u = NonfreeFactors(Factor(a + b*x + c*x**S(2) + d*x**S(3)), x)
        rubi.append(1791)
        return Dist((a + b*x + c*x**S(2) + d*x**S(3))**p/DistributeDegree(u, p), Int((e + f*x)**m*DistributeDegree(u, p), x), x)
    rule1791 = ReplacementRule(pattern1791, replacement1791)

    if load_rule:
        matcher.add(pattern1791, 1791)
    def With1792(d, p, c, f, e, x, b, a, m):
        r = Rt(-S(27)*a*d**S(2) + S(9)*b*c*d - S(2)*c**S(3) + S(3)*sqrt(S(3))*d*sqrt(S(27)*a**S(2)*d**S(2) - S(18)*a*b*c*d + S(4)*a*c**S(3) + S(4)*b**S(3)*d - b**S(2)*c**S(2)), S(3))
        rubi.append(1792)
        return Dist((c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**(-p)*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**(-p)*(a + b*x + c*x**S(2) + d*x**S(3))**p, Int((e + f*x)**m*(c + S(3)*d*x - S(2)**(S(1)/3)*(-S(6)*b*d + S(2)*c**S(2) + S(2)**(S(1)/3)*r**S(2))/(S(2)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) - sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) - sqrt(S(3))*I) + S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) - I))/(S(4)*r))**p*(c + S(3)*d*x + S(2)**(S(1)/3)*(-S(6)*b*d*(S(1) + sqrt(S(3))*I) + S(2)*c**S(2)*(S(1) + sqrt(S(3))*I) - S(2)**(S(1)/3)*I*r**S(2)*(sqrt(S(3)) + I))/(S(4)*r))**p, x), x)
    pattern1792 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons147, cons1006, cons1005)
    rule1792 = ReplacementRule(pattern1792, With1792)

    if load_rule:
        matcher.add(pattern1792, 1792)
    pattern1793 = Pattern(Integral(u_**WC('m', S(1))*v_**WC('p', S(1)), x_), cons21, cons5, cons68, cons1009, cons1010)
    def replacement1793(v, p, u, x, m):
        rubi.append(1793)
        return Int(ExpandToSum(u, x)**m*ExpandToSum(v, x)**p, x)
    rule1793 = ReplacementRule(pattern1793, replacement1793)

    if load_rule:
        matcher.add(pattern1793, 1793)
    pattern1794 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons383, cons1011, cons1012)
    def replacement1794(d, c, f, e, x, b, a, g):
        rubi.append(1794)
        return Simp(a*f*ArcTan((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(a**S(2)*(S(2)*a - c), S(2))), x)
    rule1794 = ReplacementRule(pattern1794, replacement1794)

    if load_rule:
        matcher.add(pattern1794, 1794)
    pattern1795 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/((d_ + x_**S(2)*WC('d', S(1)) + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('a', S(1)) + x_**S(3)*WC('b', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons383, cons1011, cons1013)
    def replacement1795(d, c, f, e, x, b, a, g):
        rubi.append(1795)
        return -Simp(a*f*atanh((a*b*x**S(2) + a*b + x*(S(4)*a**S(2) - S(2)*a*c + b**S(2)))/(S(2)*sqrt(a*x**S(4) + a + b*x**S(3) + b*x + c*x**S(2))*Rt(-a**S(2)*(S(2)*a - c), S(2))))/(d*Rt(-a**S(2)*(S(2)*a - c), S(2))), x)
    rule1795 = ReplacementRule(pattern1795, replacement1795)

    if load_rule:
        matcher.add(pattern1795, 1795)
    pattern1796 = Pattern(Integral((x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons1014, cons1015, cons1016)
    def replacement1796(d, p, c, e, x, b, a):
        rubi.append(1796)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x)
    rule1796 = ReplacementRule(pattern1796, replacement1796)

    if load_rule:
        matcher.add(pattern1796, 1796)

    pattern1797 = Pattern(Integral(v_**p_, x_), cons5, cons1017, cons1018, cons1015, cons1016, cons_with_1797)
    def replacement1797(x, v, p):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        rubi.append(1797)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p, x), x), x, d/(S(4)*e) + x)
    rule1797 = ReplacementRule(pattern1797, replacement1797)

    if load_rule:
        matcher.add(pattern1797, 1797)
    pattern1798 = Pattern(Integral(u_*(x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons804, cons1014, cons357)
    def replacement1798(d, p, c, e, u, x, b, a):
        rubi.append(1798)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x)
    rule1798 = ReplacementRule(pattern1798, replacement1798)

    if load_rule:
        matcher.add(pattern1798, 1798)

    pattern1799 = Pattern(Integral(u_*v_**p_, x_), cons5, cons804, cons1017, cons1018, cons357, cons_with_1799)
    def replacement1799(x, v, u, p):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        rubi.append(1799)
        return Subst(Int(SimplifyIntegrand((a - b*d/(S(8)*e) + d**S(4)/(S(256)*e**S(3)) + e*x**S(4) + x**S(2)*(c - S(3)*d**S(2)/(S(8)*e)))**p*ReplaceAll(u, Rule(x, -d/(S(4)*e) + x)), x), x), x, d/(S(4)*e) + x)
    rule1799 = ReplacementRule(pattern1799, replacement1799)

    if load_rule:
        matcher.add(pattern1799, 1799)
    pattern1800 = Pattern(Integral((a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons1019, cons246)
    def replacement1800(d, p, c, e, x, b, a):
        rubi.append(1800)
        return Dist(-S(16)*a**S(2), Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, S(1)/x + b/(S(4)*a)), x)
    rule1800 = ReplacementRule(pattern1800, replacement1800)

    if load_rule:
        matcher.add(pattern1800, 1800)

    pattern1801 = Pattern(Integral(v_**p_, x_), cons5, cons1017, cons1018, cons246, cons_with_1801)
    def replacement1801(x, v, p):

        a = Coefficient(v, x, S(0))
        b = Coefficient(v, x, S(1))
        c = Coefficient(v, x, S(2))
        d = Coefficient(v, x, S(3))
        e = Coefficient(v, x, S(4))
        rubi.append(1801)
        return Dist(-S(16)*a**S(2), Subst(Int((a*(S(256)*a**S(4)*x**S(4) + S(256)*a**S(3)*e - S(64)*a**S(2)*b*d - S(32)*a**S(2)*x**S(2)*(-S(8)*a*c + S(3)*b**S(2)) + S(16)*a*b**S(2)*c - S(3)*b**S(4))/(-S(4)*a*x + b)**S(4))**p/(-S(4)*a*x + b)**S(2), x), x, S(1)/x + b/(S(4)*a)), x)
    rule1801 = ReplacementRule(pattern1801, replacement1801)

    if load_rule:
        matcher.add(pattern1801, 1801)
    def With1802(B, d, c, e, D, A, C, x, b, a):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        rubi.append(1802)
        return -Dist(S(1)/q, Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)
    pattern1802 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons1023, cons1020, cons1021, cons1022)
    rule1802 = ReplacementRule(pattern1802, With1802)

    if load_rule:
        matcher.add(pattern1802, 1802)
    def With1803(B, d, c, e, D, A, x, b, a):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        rubi.append(1803)
        return -Dist(S(1)/q, Int((A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int((A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)
    pattern1803 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons1023, cons1020, cons1021, cons1022)
    rule1803 = ReplacementRule(pattern1803, With1803)

    if load_rule:
        matcher.add(pattern1803, 1803)
    def With1804(B, d, c, e, D, A, x, b, C, a, m):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        rubi.append(1804)
        return -Dist(S(1)/q, Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a - S(2)*C*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)
    pattern1804 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons1023, cons21, cons1020, cons1021, cons1022)
    rule1804 = ReplacementRule(pattern1804, With1804)

    if load_rule:
        matcher.add(pattern1804, 1804)
    def With1805(B, d, c, e, D, A, x, b, a, m):
        q = sqrt(S(8)*a**S(2) - S(4)*a*c + b**S(2))
        rubi.append(1805)
        return -Dist(S(1)/q, Int(x**m*(A*b - A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b - D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b - q)), x), x) + Dist(S(1)/q, Int(x**m*(A*b + A*q - S(2)*B*a + S(2)*D*a + x*(S(2)*A*a + D*b + D*q))/(S(2)*a*x**S(2) + S(2)*a + x*(b + q)), x), x)
    pattern1805 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons1023, cons21, cons1020, cons1021, cons1022)
    rule1805 = ReplacementRule(pattern1805, With1805)

    if load_rule:
        matcher.add(pattern1805, 1805)
    def With1806(B, d, e, A, C, x, b, c, a):
        q = Rt(C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        rubi.append(1806)
        return Simp(-S(2)*C**S(2)*atanh((-B*e + C*d + S(2)*C*e*x)/q)/q, x) + Simp(S(2)*C**S(2)*atanh(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q, x)
    pattern1806 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons36, cons1024, cons1025, cons1026)
    rule1806 = ReplacementRule(pattern1806, With1806)

    if load_rule:
        matcher.add(pattern1806, 1806)
    def With1807(d, c, e, A, x, b, C, a):
        q = Rt(C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        rubi.append(1807)
        return Simp(-S(2)*C**S(2)*atanh(C*(d + S(2)*e*x)/q)/q, x) + Simp(S(2)*C**S(2)*atanh(C*(A*d - S(2)*C*d*x**S(2) - S(2)*C*e*x**S(3) - S(2)*x*(A*e + C*c))/(A*q))/q, x)
    pattern1807 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons36, cons1027, cons1028, cons1029)
    rule1807 = ReplacementRule(pattern1807, With1807)

    if load_rule:
        matcher.add(pattern1807, 1807)
    def With1808(B, d, e, A, C, x, b, c, a):
        q = Rt(-C*(C*(-S(4)*c*e + d**S(2)) + S(2)*e*(-S(4)*A*e + B*d)), S(2))
        rubi.append(1808)
        return Simp(S(2)*C**S(2)*ArcTan((-B*e + C*d + S(2)*C*e*x)/q)/q, x) - Simp(S(2)*C**S(2)*ArcTan(C*(S(12)*A*B*e - S(4)*A*C*d - S(3)*B**S(2)*d + S(4)*B*C*c + S(8)*C**S(2)*e*x**S(3) + S(4)*C*x**S(2)*(-B*e + S(2)*C*d) + S(4)*C*x*(S(2)*A*e - B*d + S(2)*C*c))/(q*(-S(4)*A*C + B**S(2))))/q, x)
    pattern1808 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons36, cons1024, cons1025, cons1030)
    rule1808 = ReplacementRule(pattern1808, With1808)

    if load_rule:
        matcher.add(pattern1808, 1808)
    def With1809(d, c, e, A, x, b, C, a):
        q = Rt(-C*(-S(8)*A*e**S(2) + C*(-S(4)*c*e + d**S(2))), S(2))
        rubi.append(1809)
        return Simp(S(2)*C**S(2)*ArcTan((C*d + S(2)*C*e*x)/q)/q, x) - Simp(S(2)*C**S(2)*ArcTan(-C*(-A*d + S(2)*C*d*x**S(2) + S(2)*C*e*x**S(3) + S(2)*x*(A*e + C*c))/(A*q))/q, x)
    pattern1809 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons36, cons1027, cons1028, cons1031)
    rule1809 = ReplacementRule(pattern1809, With1809)

    if load_rule:
        matcher.add(pattern1809, 1809)
    pattern1810 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons36, cons1023, cons1032, cons1033)
    def replacement1810(B, d, c, e, D, A, C, x, b, a):
        rubi.append(1810)
        return -Dist(S(1)/(S(4)*e), Int((-S(4)*A*e + D*b + x**S(2)*(-S(4)*C*e + S(3)*D*d) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x), x) + Simp(D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e), x)
    rule1810 = ReplacementRule(pattern1810, replacement1810)

    if load_rule:
        matcher.add(pattern1810, 1810)
    pattern1811 = Pattern(Integral((x_**S(3)*WC('D', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(a_ + x_**S(4)*WC('e', S(1)) + x_**S(3)*WC('d', S(1)) + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons1023, cons1034, cons1035)
    def replacement1811(B, d, c, e, D, A, x, b, a):
        rubi.append(1811)
        return -Dist(S(1)/(S(4)*e), Int((-S(4)*A*e + D*b + S(3)*D*d*x**S(2) + S(2)*x*(-S(2)*B*e + D*c))/(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4)), x), x) + Simp(D*log(a + b*x + c*x**S(2) + d*x**S(3) + e*x**S(4))/(S(4)*e), x)
    rule1811 = ReplacementRule(pattern1811, replacement1811)

    if load_rule:
        matcher.add(pattern1811, 1811)
    pattern1812 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons71, cons1036)
    def replacement1812(d, c, f, e, u, x, b, a):
        rubi.append(1812)
        return -Dist(a/(f*(-a*d + b*c)), Int(u*sqrt(c + d*x)/x, x), x) + Dist(c/(e*(-a*d + b*c)), Int(u*sqrt(a + b*x)/x, x), x)
    rule1812 = ReplacementRule(pattern1812, replacement1812)

    if load_rule:
        matcher.add(pattern1812, 1812)
    pattern1813 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons71, cons1037)
    def replacement1813(d, c, f, e, u, x, b, a):
        rubi.append(1813)
        return Dist(b/(f*(-a*d + b*c)), Int(u*sqrt(c + d*x), x), x) - Dist(d/(e*(-a*d + b*c)), Int(u*sqrt(a + b*x), x), x)
    rule1813 = ReplacementRule(pattern1813, replacement1813)

    if load_rule:
        matcher.add(pattern1813, 1813)
    pattern1814 = Pattern(Integral(u_/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons1038, cons1039)
    def replacement1814(d, c, f, e, u, x, b, a):
        rubi.append(1814)
        return Dist(e, Int(u*sqrt(a + b*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x), x) - Dist(f, Int(u*sqrt(c + d*x)/(a*e**S(2) - c*f**S(2) + x*(b*e**S(2) - d*f**S(2))), x), x)
    rule1814 = ReplacementRule(pattern1814, replacement1814)

    if load_rule:
        matcher.add(pattern1814, 1814)
    pattern1815 = Pattern(Integral(WC('u', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons1040, cons1041)
    def replacement1815(d, p, u, x, b, c, a, n):
        rubi.append(1815)
        return Dist(S(1)/(a*c), Int(u*sqrt(a + b*x**(S(2)*n)), x), x) - Dist(b/(a*d), Int(u*x**n, x), x)
    rule1815 = ReplacementRule(pattern1815, replacement1815)

    if load_rule:
        matcher.add(pattern1815, 1815)
    pattern1816 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('d', S(1)) + sqrt(x_**WC('p', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1040, cons1042)
    def replacement1816(d, p, x, b, c, a, m, n):
        rubi.append(1816)
        return Dist(c, Int(x**m*sqrt(a + b*x**(S(2)*n))/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x), x) - Dist(d, Int(x**(m + n)/(a*c**S(2) + x**(S(2)*n)*(b*c**S(2) - d**S(2))), x), x)
    rule1816 = ReplacementRule(pattern1816, replacement1816)

    if load_rule:
        matcher.add(pattern1816, 1816)
    def With1817(d, f, e, x, b, a):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        rubi.append(1817)
        return Dist(r/(S(3)*a), Int(S(1)/((r + s*x)*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(r/(S(3)*a), Int((S(2)*r - s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x), x)
    pattern1817 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons27, cons48, cons125, cons468)
    rule1817 = ReplacementRule(pattern1817, With1817)

    if load_rule:
        matcher.add(pattern1817, 1817)
    def With1818(d, f, x, b, a):
        r = Numerator(Rt(a/b, S(3)))
        s = Denominator(Rt(a/b, S(3)))
        rubi.append(1818)
        return Dist(r/(S(3)*a), Int(S(1)/(sqrt(d + f*x**S(2))*(r + s*x)), x), x) + Dist(r/(S(3)*a), Int((S(2)*r - s*x)/(sqrt(d + f*x**S(2))*(r**S(2) - r*s*x + s**S(2)*x**S(2))), x), x)
    pattern1818 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons27, cons125, cons468)
    rule1818 = ReplacementRule(pattern1818, With1818)

    if load_rule:
        matcher.add(pattern1818, 1818)
    def With1819(d, f, e, x, b, a):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        rubi.append(1819)
        return Dist(r/(S(3)*a), Int(S(1)/((r - s*x)*sqrt(d + e*x + f*x**S(2))), x), x) + Dist(r/(S(3)*a), Int((S(2)*r + s*x)/(sqrt(d + e*x + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x), x)
    pattern1819 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons27, cons48, cons125, cons469)
    rule1819 = ReplacementRule(pattern1819, With1819)

    if load_rule:
        matcher.add(pattern1819, 1819)
    def With1820(d, f, x, b, a):
        r = Numerator(Rt(-a/b, S(3)))
        s = Denominator(Rt(-a/b, S(3)))
        rubi.append(1820)
        return Dist(r/(S(3)*a), Int(S(1)/(sqrt(d + f*x**S(2))*(r - s*x)), x), x) + Dist(r/(S(3)*a), Int((S(2)*r + s*x)/(sqrt(d + f*x**S(2))*(r**S(2) + r*s*x + s**S(2)*x**S(2))), x), x)
    pattern1820 = Pattern(Integral(S(1)/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons27, cons125, cons469)
    rule1820 = ReplacementRule(pattern1820, With1820)

    if load_rule:
        matcher.add(pattern1820, 1820)
    pattern1821 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons1043)
    def replacement1821(d, c, e, x, b, a):
        rubi.append(1821)
        return Dist(d, Int(S(1)/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) - Dist(e, Int(x/((d**S(2) - e**S(2)*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
    rule1821 = ReplacementRule(pattern1821, replacement1821)

    if load_rule:
        matcher.add(pattern1821, 1821)
    pattern1822 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons297)
    def replacement1822(d, e, x, c, a):
        rubi.append(1822)
        return Dist(d, Int(S(1)/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x), x) - Dist(e, Int(x/(sqrt(a + c*x**S(4))*(d**S(2) - e**S(2)*x**S(2))), x), x)
    rule1822 = ReplacementRule(pattern1822, replacement1822)

    if load_rule:
        matcher.add(pattern1822, 1822)
    pattern1823 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons1044, cons1045)
    def replacement1823(d, c, e, x, b, a):
        rubi.append(1823)
        return -Dist(c/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Simp(e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), x)
    rule1823 = ReplacementRule(pattern1823, replacement1823)

    if load_rule:
        matcher.add(pattern1823, 1823)
    pattern1824 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons1044, cons1046)
    def replacement1824(d, c, e, x, b, a):
        rubi.append(1824)
        return -Dist(c/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((b*d*e**S(2) + S(2)*c*d**S(3))/(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4)), Int(S(1)/((d + e*x)*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) - Simp(e**S(3)*sqrt(a + b*x**S(2) + c*x**S(4))/((d + e*x)*(a*e**S(4) + b*d**S(2)*e**S(2) + c*d**S(4))), x)
    rule1824 = ReplacementRule(pattern1824, replacement1824)

    if load_rule:
        matcher.add(pattern1824, 1824)
    pattern1825 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))**S(2)), x_), cons2, cons7, cons27, cons48, cons1047)
    def replacement1825(d, e, x, c, a):
        rubi.append(1825)
        return -Dist(c/(a*e**S(4) + c*d**S(4)), Int((d**S(2) - e**S(2)*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist(S(2)*c*d**S(3)/(a*e**S(4) + c*d**S(4)), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x)), x), x) - Simp(e**S(3)*sqrt(a + c*x**S(4))/((d + e*x)*(a*e**S(4) + c*d**S(4))), x)
    rule1825 = ReplacementRule(pattern1825, replacement1825)

    if load_rule:
        matcher.add(pattern1825, 1825)
    pattern1826 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons1048, cons705)
    def replacement1826(B, d, c, e, A, x, b, a):
        rubi.append(1826)
        return Dist(A, Subst(Int(S(1)/(d - x**S(2)*(-S(2)*a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1826 = ReplacementRule(pattern1826, replacement1826)

    if load_rule:
        matcher.add(pattern1826, 1826)
    pattern1827 = Pattern(Integral((A_ + x_**S(2)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons34, cons35, cons1048, cons705)
    def replacement1827(B, d, e, A, x, c, a):
        rubi.append(1827)
        return Dist(A, Subst(Int(S(1)/(S(2)*a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))), x)
    rule1827 = ReplacementRule(pattern1827, replacement1827)

    if load_rule:
        matcher.add(pattern1827, 1827)
    pattern1828 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons382, cons1049)
    def replacement1828(B, d, c, f, e, A, x, b, a):
        rubi.append(1828)
        return Dist(A, Subst(Int(S(1)/(d - x**S(2)*(-a*e + b*d)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1828 = ReplacementRule(pattern1828, replacement1828)

    if load_rule:
        matcher.add(pattern1828, 1828)
    pattern1829 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons125, cons34, cons35, cons382, cons1049)
    def replacement1829(B, d, f, e, A, x, c, a):
        rubi.append(1829)
        return Dist(A, Subst(Int(S(1)/(a*e*x**S(2) + d), x), x, x/sqrt(a + c*x**S(4))), x)
    rule1829 = ReplacementRule(pattern1829, replacement1829)

    if load_rule:
        matcher.add(pattern1829, 1829)
    pattern1830 = Pattern(Integral((A_ + x_**S(4)*WC('B', S(1)))/((d_ + x_**S(4)*WC('f', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons125, cons34, cons35, cons382, cons1049)
    def replacement1830(B, d, c, f, A, x, b, a):
        rubi.append(1830)
        return Dist(A, Subst(Int(S(1)/(-b*d*x**S(2) + d), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1830 = ReplacementRule(pattern1830, replacement1830)

    if load_rule:
        matcher.add(pattern1830, 1830)
    pattern1831 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons1050, cons697)
    def replacement1831(d, c, e, x, b, a):
        rubi.append(1831)
        return Dist(a/d, Subst(Int(S(1)/(-S(2)*b*x**S(2) + x**S(4)*(-S(4)*a*c + b**S(2)) + S(1)), x), x, x/sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1831 = ReplacementRule(pattern1831, replacement1831)

    if load_rule:
        matcher.add(pattern1831, 1831)
    def With1832(d, c, e, x, b, a):
        q = sqrt(-S(4)*a*c + b**S(2))
        rubi.append(1832)
        return Simp(sqrt(S(2))*a*sqrt(-b + q)*atanh(sqrt(S(2))*x*sqrt(-b + q)*(b + S(2)*c*x**S(2) + q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2))), x) - Simp(sqrt(S(2))*a*sqrt(b + q)*ArcTan(sqrt(S(2))*x*sqrt(b + q)*(b + S(2)*c*x**S(2) - q)/(S(4)*sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-a*c, S(2))))/(S(4)*d*Rt(-a*c, S(2))), x)
    pattern1832 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons1050, cons709)
    rule1832 = ReplacementRule(pattern1832, With1832)

    if load_rule:
        matcher.add(pattern1832, 1832)
    pattern1833 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    def replacement1833(d, c, f, e, x, b, a):
        rubi.append(1833)
        return Dist(a, Int(S(1)/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x) - Dist(b, Int(x/((a**S(2) - b**S(2)*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x), x)
    rule1833 = ReplacementRule(pattern1833, replacement1833)

    if load_rule:
        matcher.add(pattern1833, 1833)
    pattern1834 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*sqrt(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons1051, cons1052)
    def replacement1834(d, c, f, e, h, x, b, a, g):
        rubi.append(1834)
        return Simp(S(2)*sqrt(d + e*x + f*sqrt(a + b*x + c*x**S(2)))*(S(9)*c**S(2)*f*g*h*x**S(2) + S(3)*c**S(2)*f*h**S(2)*x**S(3) + c*f*x*(a*h**S(2) - b*g*h + S(10)*c*g**S(2)) + f*(S(2)*a*b*h**S(2) - S(3)*a*c*g*h - S(2)*b**S(2)*g*h + S(5)*b*c*g**S(2)) - (-d*h + e*g)*sqrt(a + b*x + c*x**S(2))*(-S(2)*b*h + S(5)*c*g + c*h*x))/(S(15)*c**S(2)*f*(g + h*x)), x)
    rule1834 = ReplacementRule(pattern1834, replacement1834)

    if load_rule:
        matcher.add(pattern1834, 1834)
    pattern1835 = Pattern(Integral((u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons125, cons208, cons209, cons796, cons797, cons21, cons4, cons68, cons818, cons1053, cons1054)
    def replacement1835(v, k, j, f, h, g, u, x, m, n):
        rubi.append(1835)
        return Int((g + h*x)**m*(f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n, x)
    rule1835 = ReplacementRule(pattern1835, replacement1835)

    if load_rule:
        matcher.add(pattern1835, 1835)
    pattern1836 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1055, cons38)
    def replacement1836(d, p, c, f, e, h, x, n, b, a, g):
        rubi.append(1836)
        return Dist(S(2), Subst(Int((g + h*x**n)**p*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))/(b*f**S(2) - S(2)*d*e + S(2)*e*x)**S(2), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)
    rule1836 = ReplacementRule(pattern1836, replacement1836)

    if load_rule:
        matcher.add(pattern1836, 1836)
    pattern1837 = Pattern(Integral(((x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1055, cons38)
    def replacement1837(d, p, f, e, h, x, n, c, a, g):
        rubi.append(1837)
        return Dist(S(1)/(S(2)*e), Subst(Int((g + h*x**n)**p*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))/(d - x)**S(2), x), x, d + e*x + f*sqrt(a + c*x**S(2))), x)
    rule1837 = ReplacementRule(pattern1837, replacement1837)

    if load_rule:
        matcher.add(pattern1837, 1837)
    pattern1838 = Pattern(Integral(((u_ + sqrt(v_)*WC('f', S(1)))**n_*WC('h', S(1)) + WC('g', S(0)))**WC('p', S(1)), x_), cons125, cons208, cons209, cons4, cons68, cons818, cons819, cons1056, cons38)
    def replacement1838(v, p, f, h, u, x, n, g):
        rubi.append(1838)
        return Int((g + h*(f*sqrt(ExpandToSum(v, x)) + ExpandToSum(u, x))**n)**p, x)
    rule1838 = ReplacementRule(pattern1838, replacement1838)

    if load_rule:
        matcher.add(pattern1838, 1838)
    pattern1839 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons7, cons48, cons125, cons208, cons209, cons4, cons1055, cons17)
    def replacement1839(f, e, h, x, n, c, a, m, g):
        rubi.append(1839)
        return Dist(S(2)**(-m + S(-1))*e**(-m + S(-1)), Subst(Int(x**(-m + n + S(-2))*(a*f**S(2) + x**S(2))*(-a*f**S(2)*h + S(2)*e*g*x + h*x**S(2))**m, x), x, e*x + f*sqrt(a + c*x**S(2))), x)
    rule1839 = ReplacementRule(pattern1839, replacement1839)

    if load_rule:
        matcher.add(pattern1839, 1839)
    pattern1840 = Pattern(Integral(x_**WC('p', S(1))*(g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)))**WC('n', S(1)), x_), cons2, cons7, cons48, cons125, cons208, cons224, cons4, cons1055, cons1057, cons1058, cons1059)
    def replacement1840(p, f, e, i, g, x, c, a, m, n):
        rubi.append(1840)
        return Dist(S(2)**(-S(2)*m - p + S(-1))*e**(-p + S(-1))*f**(-S(2)*m)*(i/c)**m, Subst(Int(x**(-S(2)*m + n - p + S(-2))*(-a*f**S(2) + x**S(2))**p*(a*f**S(2) + x**S(2))**(S(2)*m + S(1)), x), x, e*x + f*sqrt(a + c*x**S(2))), x)
    rule1840 = ReplacementRule(pattern1840, replacement1840)

    if load_rule:
        matcher.add(pattern1840, 1840)
    pattern1841 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons4, cons1055, cons1057, cons1060, cons515, cons1059)
    def replacement1841(d, c, f, e, h, i, x, n, b, a, m, g):
        rubi.append(1841)
        return Dist(S(2)*f**(-S(2)*m)*(i/c)**m, Subst(Int(x**n*(b*f**S(2) - S(2)*d*e + S(2)*e*x)**(-S(2)*m + S(-2))*(d**S(2)*e + e*x**S(2) - f**S(2)*(-a*e + b*d) - x*(-b*f**S(2) + S(2)*d*e))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)
    rule1841 = ReplacementRule(pattern1841, replacement1841)

    if load_rule:
        matcher.add(pattern1841, 1841)
    pattern1842 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons224, cons4, cons1055, cons1057, cons515, cons1059)
    def replacement1842(d, f, e, i, g, x, c, a, m, n):
        rubi.append(1842)
        return Dist(S(2)**(-S(2)*m + S(-1))*f**(-S(2)*m)*(i/c)**m/e, Subst(Int(x**n*(-d + x)**(-S(2)*m + S(-2))*(a*f**S(2) + d**S(2) - S(2)*d*x + x**S(2))**(S(2)*m + S(1)), x), x, d + e*x + f*sqrt(a + c*x**S(2))), x)
    rule1842 = ReplacementRule(pattern1842, replacement1842)

    if load_rule:
        matcher.add(pattern1842, 1842)
    pattern1843 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons4, cons1055, cons1057, cons1060, cons73, cons1061)
    def replacement1843(d, c, f, e, h, i, x, n, b, a, m, g):
        rubi.append(1843)
        return Dist((i/c)**(m + S(-1)/2)*sqrt(g + h*x + i*x**S(2))/sqrt(a + b*x + c*x**S(2)), Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x), x)
    rule1843 = ReplacementRule(pattern1843, replacement1843)

    if load_rule:
        matcher.add(pattern1843, 1843)
    pattern1844 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons224, cons4, cons1055, cons1057, cons73, cons1061)
    def replacement1844(d, f, e, i, g, x, c, a, m, n):
        rubi.append(1844)
        return Dist((i/c)**(m + S(-1)/2)*sqrt(g + i*x**S(2))/sqrt(a + c*x**S(2)), Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x), x)
    rule1844 = ReplacementRule(pattern1844, replacement1844)

    if load_rule:
        matcher.add(pattern1844, 1844)
    pattern1845 = Pattern(Integral((x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1))*(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons4, cons1055, cons1057, cons1060, cons951, cons1061)
    def replacement1845(d, c, f, e, h, i, x, n, b, a, m, g):
        rubi.append(1845)
        return Dist((i/c)**(m + S(1)/2)*sqrt(a + b*x + c*x**S(2))/sqrt(g + h*x + i*x**S(2)), Int((a + b*x + c*x**S(2))**m*(d + e*x + f*sqrt(a + b*x + c*x**S(2)))**n, x), x)
    rule1845 = ReplacementRule(pattern1845, replacement1845)

    if load_rule:
        matcher.add(pattern1845, 1845)
    pattern1846 = Pattern(Integral((g_ + x_**S(2)*WC('i', S(1)))**WC('m', S(1))*(x_*WC('e', S(1)) + sqrt(a_ + x_**S(2)*WC('c', S(1)))*WC('f', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons224, cons4, cons1055, cons1057, cons951, cons1061)
    def replacement1846(d, f, e, i, g, x, c, a, m, n):
        rubi.append(1846)
        return Dist((i/c)**(m + S(1)/2)*sqrt(a + c*x**S(2))/sqrt(g + i*x**S(2)), Int((a + c*x**S(2))**m*(d + e*x + f*sqrt(a + c*x**S(2)))**n, x), x)
    rule1846 = ReplacementRule(pattern1846, replacement1846)

    if load_rule:
        matcher.add(pattern1846, 1846)
    pattern1847 = Pattern(Integral(w_**WC('m', S(1))*(u_ + (sqrt(v_)*WC('k', S(1)) + WC('j', S(0)))*WC('f', S(1)))**WC('n', S(1)), x_), cons125, cons796, cons797, cons21, cons4, cons68, cons1062, cons1063, cons1064)
    def replacement1847(v, k, j, f, u, x, w, m, n):
        rubi.append(1847)
        return Int((f*k*sqrt(ExpandToSum(v, x)) + ExpandToSum(f*j + u, x))**n*ExpandToSum(w, x)**m, x)
    rule1847 = ReplacementRule(pattern1847, replacement1847)

    if load_rule:
        matcher.add(pattern1847, 1847)
    pattern1848 = Pattern(Integral(S(1)/((a_ + x_**WC('n', S(1))*WC('b', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + (a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons4, cons1065)
    def replacement1848(d, p, c, x, b, a, n):
        rubi.append(1848)
        return Dist(S(1)/a, Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*(a + b*x**n)**(S(2)/n))), x)
    rule1848 = ReplacementRule(pattern1848, replacement1848)

    if load_rule:
        matcher.add(pattern1848, 1848)
    pattern1849 = Pattern(Integral(sqrt(a_ + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons1066)
    def replacement1849(d, c, x, b, a):
        rubi.append(1849)
        return Simp(S(2)*a*x/sqrt(a + b*sqrt(c + d*x**S(2))), x) + Simp(S(2)*b**S(2)*d*x**S(3)/(S(3)*(a + b*sqrt(c + d*x**S(2)))**(S(3)/2)), x)
    rule1849 = ReplacementRule(pattern1849, replacement1849)

    if load_rule:
        matcher.add(pattern1849, 1849)
    pattern1850 = Pattern(Integral(sqrt(x_**S(2)*WC('a', S(1)) + x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons1067, cons1068)
    def replacement1850(d, c, x, b, a):
        rubi.append(1850)
        return Dist(sqrt(S(2))*b/a, Subst(Int(S(1)/sqrt(S(1) + x**S(2)/a), x), x, a*x + b*sqrt(c + d*x**S(2))), x)
    rule1850 = ReplacementRule(pattern1850, replacement1850)

    if load_rule:
        matcher.add(pattern1850, 1850)
    pattern1851 = Pattern(Integral(sqrt(x_*(x_*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1)))*WC('e', S(1)))/(x_*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons1067, cons1069)
    def replacement1851(d, c, e, x, b, a):
        rubi.append(1851)
        return Int(sqrt(a*e*x**S(2) + b*e*x*sqrt(c + d*x**S(2)))/(x*sqrt(c + d*x**S(2))), x)
    rule1851 = ReplacementRule(pattern1851, replacement1851)

    if load_rule:
        matcher.add(pattern1851, 1851)
    pattern1852 = Pattern(Integral(sqrt(x_**S(2)*WC('c', S(1)) + sqrt(a_ + x_**S(4)*WC('b', S(1)))*WC('d', S(1)))/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons1070)
    def replacement1852(d, c, x, b, a):
        rubi.append(1852)
        return Dist(d, Subst(Int(S(1)/(-S(2)*c*x**S(2) + S(1)), x), x, x/sqrt(c*x**S(2) + d*sqrt(a + b*x**S(4)))), x)
    rule1852 = ReplacementRule(pattern1852, replacement1852)

    if load_rule:
        matcher.add(pattern1852, 1852)
    pattern1853 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sqrt(x_**S(2)*WC('b', S(1)) + sqrt(a_ + x_**S(4)*WC('e', S(1))))/sqrt(a_ + x_**S(4)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons1071, cons43)
    def replacement1853(d, e, x, b, c, a, m):
        rubi.append(1853)
        return Dist(S(1)/2 - I/S(2), Int((c + d*x)**m/sqrt(sqrt(a) - I*b*x**S(2)), x), x) + Dist(S(1)/2 + I/S(2), Int((c + d*x)**m/sqrt(sqrt(a) + I*b*x**S(2)), x), x)
    rule1853 = ReplacementRule(pattern1853, replacement1853)

    if load_rule:
        matcher.add(pattern1853, 1853)
    def With1854(d, c, x, b, a):
        q = Rt(b/a, S(3))
        rubi.append(1854)
        return Dist(d/(-c*q + d*(S(1) + sqrt(S(3)))), Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) - Dist(q/(-c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1854 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons481, cons479)
    rule1854 = ReplacementRule(pattern1854, With1854)

    if load_rule:
        matcher.add(pattern1854, 1854)
    def With1855(d, c, x, b, a):
        q = Rt(-b/a, S(3))
        rubi.append(1855)
        return Dist(d/(c*q + d*(S(1) + sqrt(S(3)))), Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist(q/(c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1855 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons481, cons480)
    rule1855 = ReplacementRule(pattern1855, With1855)

    if load_rule:
        matcher.add(pattern1855, 1855)
    def With1856(d, c, x, b, a):
        q = Rt(-b/a, S(3))
        rubi.append(1856)
        return Dist(d/(c*q + d*(-sqrt(S(3)) + S(1))), Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist(q/(c*q + d*(-sqrt(S(3)) + S(1))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1856 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons482, cons479)
    rule1856 = ReplacementRule(pattern1856, With1856)

    if load_rule:
        matcher.add(pattern1856, 1856)
    def With1857(d, c, x, b, a):
        q = Rt(b/a, S(3))
        rubi.append(1857)
        return Dist(d/(-c*q + d*(-sqrt(S(3)) + S(1))), Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) - Dist(q/(-c*q + d*(-sqrt(S(3)) + S(1))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    pattern1857 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons482, cons480)
    rule1857 = ReplacementRule(pattern1857, With1857)

    if load_rule:
        matcher.add(pattern1857, 1857)

    pattern1858 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons481, cons479, cons_with_1858)
    def replacement1858(d, c, f, e, x, b, a):

        q = Rt(b/a, S(3))
        rubi.append(1858)
        return Dist(S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(q*x + S(1))/(q*sqrt((q*x + S(1))/(q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(-sqrt(S(3)) + S(1)) + x*(-c*q + d*(S(1) + sqrt(S(3)))))), x), x, (-q*x + S(-1) + sqrt(S(3)))/(q*x + S(1) + sqrt(S(3)))), x)
    rule1858 = ReplacementRule(pattern1858, replacement1858)

    if load_rule:
        matcher.add(pattern1858, 1858)

    pattern1859 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons481, cons479, cons_with_1859)
    def replacement1859(d, c, f, e, x, b, a):

        q = Rt(b/a, S(3))
        rubi.append(1859)
        return Dist((-c*f + d*e)/(-c*q + d*(S(1) + sqrt(S(3)))), Int((q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((-e*q + f*(S(1) + sqrt(S(3))))/(-c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    rule1859 = ReplacementRule(pattern1859, replacement1859)

    if load_rule:
        matcher.add(pattern1859, 1859)

    pattern1860 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons481, cons480, cons_with_1860)
    def replacement1860(d, c, f, e, x, b, a):

        q = Rt(-b/a, S(3))
        rubi.append(1860)
        return Dist(-S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(-q*x + S(1))/(q*sqrt((-q*x + S(1))/(-q*x + S(1) + sqrt(S(3)))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) - S(4)*sqrt(S(3)) + S(7))*(c*q + d*(-sqrt(S(3)) + S(1)) + x*(c*q + d*(S(1) + sqrt(S(3)))))), x), x, (q*x + S(-1) + sqrt(S(3)))/(-q*x + S(1) + sqrt(S(3)))), x)
    rule1860 = ReplacementRule(pattern1860, replacement1860)

    if load_rule:
        matcher.add(pattern1860, 1860)

    pattern1861 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons481, cons480, cons_with_1861)
    def replacement1861(d, c, f, e, x, b, a):

        q = Rt(-b/a, S(3))
        rubi.append(1861)
        return Dist((-c*f + d*e)/(c*q + d*(S(1) + sqrt(S(3)))), Int((-q*x + S(1) + sqrt(S(3)))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((e*q + f*(S(1) + sqrt(S(3))))/(c*q + d*(S(1) + sqrt(S(3)))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    rule1861 = ReplacementRule(pattern1861, replacement1861)

    if load_rule:
        matcher.add(pattern1861, 1861)

    pattern1862 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons482, cons479, cons_with_1862)
    def replacement1862(d, c, f, e, x, b, a):

        q = Rt(-b/a, S(3))
        rubi.append(1862)
        return Dist(S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) + q*x + S(1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(-q*x + S(1))/(q*sqrt(-(-q*x + S(1))/(-q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(c*q + d*(S(1) + sqrt(S(3))) + x*(c*q + d*(-sqrt(S(3)) + S(1))))), x), x, (-q*x + S(1) + sqrt(S(3)))/(q*x + S(-1) + sqrt(S(3)))), x)
    rule1862 = ReplacementRule(pattern1862, replacement1862)

    if load_rule:
        matcher.add(pattern1862, 1862)

    pattern1863 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons482, cons479, cons_with_1863)
    def replacement1863(d, c, f, e, x, b, a):

        q = Rt(-b/a, S(3))
        rubi.append(1863)
        return Dist((-c*f + d*e)/(c*q + d*(-sqrt(S(3)) + S(1))), Int((-q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((e*q + f*(-sqrt(S(3)) + S(1)))/(c*q + d*(-sqrt(S(3)) + S(1))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    rule1863 = ReplacementRule(pattern1863, replacement1863)

    if load_rule:
        matcher.add(pattern1863, 1863)

    pattern1864 = Pattern(Integral((e_ + x_*WC('f', S(1)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons482, cons480, cons_with_1864)
    def replacement1864(d, c, f, e, x, b, a):

        q = Rt(b/a, S(3))
        rubi.append(1864)
        return Dist(-S(4)*S(3)**(S(1)/4)*f*sqrt((q**S(2)*x**S(2) - q*x + S(1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(sqrt(S(3)) + S(2))*(q*x + S(1))/(q*sqrt(-(q*x + S(1))/(q*x - sqrt(S(3)) + S(1))**S(2))*sqrt(a + b*x**S(3))), Subst(Int(S(1)/(sqrt(-x**S(2) + S(1))*sqrt(x**S(2) + S(4)*sqrt(S(3)) + S(7))*(-c*q + d*(S(1) + sqrt(S(3))) + x*(-c*q + d*(-sqrt(S(3)) + S(1))))), x), x, (q*x + S(1) + sqrt(S(3)))/(-q*x + S(-1) + sqrt(S(3)))), x)
    rule1864 = ReplacementRule(pattern1864, replacement1864)

    if load_rule:
        matcher.add(pattern1864, 1864)

    pattern1865 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_**S(3)*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons482, cons480, cons_with_1865)
    def replacement1865(d, c, f, e, x, b, a):

        q = Rt(b/a, S(3))
        rubi.append(1865)
        return Dist((-c*f + d*e)/(-c*q + d*(-sqrt(S(3)) + S(1))), Int((q*x - sqrt(S(3)) + S(1))/(sqrt(a + b*x**S(3))*(c + d*x)), x), x) + Dist((-e*q + f*(-sqrt(S(3)) + S(1)))/(-c*q + d*(-sqrt(S(3)) + S(1))), Int(S(1)/sqrt(a + b*x**S(3)), x), x)
    rule1865 = ReplacementRule(pattern1865, replacement1865)

    if load_rule:
        matcher.add(pattern1865, 1865)
    pattern1866 = Pattern(Integral(x_**WC('m', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons1072, cons500)
    def replacement1866(d, c, e, x, b, a, m, n):
        rubi.append(1866)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)/(c + d*x + e*sqrt(a + b*x)), x), x, x**n), x)
    rule1866 = ReplacementRule(pattern1866, replacement1866)

    if load_rule:
        matcher.add(pattern1866, 1866)
    pattern1867 = Pattern(Integral(WC('u', S(1))/(c_ + x_**n_*WC('d', S(1)) + sqrt(a_ + x_**n_*WC('b', S(1)))*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1072)
    def replacement1867(d, c, e, u, x, b, a, n):
        rubi.append(1867)
        return Dist(c, Int(u/(-a*e**S(2) + c**S(2) + c*d*x**n), x), x) - Dist(a*e, Int(u/(sqrt(a + b*x**n)*(-a*e**S(2) + c**S(2) + c*d*x**n)), x), x)
    rule1867 = ReplacementRule(pattern1867, replacement1867)

    if load_rule:
        matcher.add(pattern1867, 1867)
    pattern1868 = Pattern(Integral((A_ + x_**n_*WC('B', S(1)))/(a_ + x_**S(2)*WC('b', S(1)) + x_**n2_*WC('d', S(1)) + x_**n_*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons34, cons35, cons4, cons46, cons963, cons1073, cons1074)
    def replacement1868(B, d, c, A, n2, x, b, a, n):
        rubi.append(1868)
        return Dist(A**S(2)*(n + S(-1)), Subst(Int(S(1)/(A**S(2)*b*x**S(2)*(n + S(-1))**S(2) + a), x), x, x/(A*(n + S(-1)) - B*x**n)), x)
    rule1868 = ReplacementRule(pattern1868, replacement1868)

    if load_rule:
        matcher.add(pattern1868, 1868)
    pattern1869 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('n', S(1))*WC('B', S(1)))/(a_ + x_**n2_*WC('d', S(1)) + x_**WC('k', S(1))*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons34, cons35, cons21, cons4, cons46, cons1075, cons1076, cons1077)
    def replacement1869(B, d, k, c, A, n2, x, b, a, m, n):
        rubi.append(1869)
        return Dist(A**S(2)*(m - n + S(1))/(m + S(1)), Subst(Int(S(1)/(A**S(2)*b*x**S(2)*(m - n + S(1))**S(2) + a), x), x, x**(m + S(1))/(A*(m - n + S(1)) + B*x**n*(m + S(1)))), x)
    rule1869 = ReplacementRule(pattern1869, replacement1869)

    if load_rule:
        matcher.add(pattern1869, 1869)
    pattern1870 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons46, cons931, cons226, cons702)
    def replacement1870(d, n3, p, c, f, e, g, n2, x, b, a, n):
        rubi.append(1870)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1870 = ReplacementRule(pattern1870, replacement1870)

    if load_rule:
        matcher.add(pattern1870, 1870)
    pattern1871 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons46, cons226, cons702)
    def replacement1871(d, p, c, f, e, n2, x, b, a, n):
        rubi.append(1871)
        return -Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*e - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*d*(n*p + n + S(1)) - x**n*(-S(2)*a*c*e*(n*(S(2)*p + S(3)) + S(1)) + b*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d + x**n*(-S(2)*a*c*e + b*(a*f + c*d)))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1871 = ReplacementRule(pattern1871, replacement1871)

    if load_rule:
        matcher.add(pattern1871, 1871)
    pattern1872 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons208, cons4, cons46, cons931, cons226, cons702)
    def replacement1872(d, n3, p, c, e, g, n2, x, b, a, n):
        rubi.append(1872)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a*b*(a*g + c*e) + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(a*b**S(2)*g*(n*(p + S(2)) + S(1)) - S(2)*a*c*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*(a*g + c*e) - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(-a*b**S(2)*g - S(2)*a*c*(-a*g + c*e) + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1872 = ReplacementRule(pattern1872, replacement1872)

    if load_rule:
        matcher.add(pattern1872, 1872)
    pattern1873 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1))), x_), cons2, cons3, cons7, cons27, cons125, cons208, cons4, cons46, cons931, cons226, cons702)
    def replacement1873(d, n3, p, c, f, g, n2, x, b, a, n):
        rubi.append(1873)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c*(-a*f + c*d) + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c*(a*f + c*d)))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1873 = ReplacementRule(pattern1873, replacement1873)

    if load_rule:
        matcher.add(pattern1873, 1873)
    pattern1874 = Pattern(Integral((d_ + x_**n2_*WC('f', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons4, cons46, cons226, cons702)
    def replacement1874(d, p, c, f, n2, x, b, a, n):
        rubi.append(1874)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))) + b**S(2)*d*(n*p + n + S(1)) + b*x**n*(a*f + c*d)*(n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*(-a*f + c*d) + b**S(2)*d + b*x**n*(a*f + c*d))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1874 = ReplacementRule(pattern1874, replacement1874)

    if load_rule:
        matcher.add(pattern1874, 1874)
    pattern1875 = Pattern(Integral((d_ + g_*x_**n3_)*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons208, cons4, cons46, cons931, cons226, cons702)
    def replacement1875(d, n3, p, c, g, n2, x, b, a, n):
        rubi.append(1875)
        return -Dist(S(1)/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(a**S(2)*b*g + S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - b**S(2)*c*d*(n*p + n + S(1)) + x**n*(-S(2)*a**S(2)*c*g*(n + S(1)) + a*b**S(2)*g*(n*(p + S(2)) + S(1)) - b*c**S(2)*d*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a**S(2)*b*g - S(2)*a*c**S(2)*d + b**S(2)*c*d + x**n*(S(2)*a**S(2)*c*g - a*b**S(2)*g + b*c**S(2)*d))/(a*c*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1875 = ReplacementRule(pattern1875, replacement1875)

    if load_rule:
        matcher.add(pattern1875, 1875)
    pattern1876 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons46, cons931, cons702)
    def replacement1876(d, n3, p, f, e, g, n2, x, c, a, n):
        rubi.append(1876)
        return -Dist(-S(1)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))) - S(2)*a*c*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*x**n*(-a*g + c*e) - S(2)*a*c*(-a*f + c*d))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), x)
    rule1876 = ReplacementRule(pattern1876, replacement1876)

    if load_rule:
        matcher.add(pattern1876, 1876)
    pattern1877 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n2_*WC('f', S(1)) + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons4, cons46, cons702)
    def replacement1877(d, p, f, e, n2, x, c, a, n):
        rubi.append(1877)
        return -Dist(-S(1)/(S(4)*a**S(2)*c*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c*e*x**n*(n*(S(2)*p + S(3)) + S(1)) - S(2)*a*(a*f - c*d*(S(2)*n*(p + S(1)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c*e*x**n - S(2)*a*(-a*f + c*d))/(S(4)*a**S(2)*c*n*(p + S(1))), x)
    rule1877 = ReplacementRule(pattern1877, replacement1877)

    if load_rule:
        matcher.add(pattern1877, 1877)
    pattern1878 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + g_*x_**n3_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons208, cons4, cons46, cons931, cons702)
    def replacement1878(d, n3, p, e, g, n2, x, c, a, n):
        rubi.append(1878)
        return -Dist(-S(1)/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*Simp(S(2)*a*c**S(2)*d*(S(2)*n*(p + S(1)) + S(1)) - S(2)*a*c*x**n*(a*g*(n + S(1)) - c*e*(n*(S(2)*p + S(3)) + S(1))), x), x), x) - Simp(-x*(a + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c**S(2)*d - S(2)*a*c*x**n*(-a*g + c*e))/(S(4)*a**S(2)*c**S(2)*n*(p + S(1))), x)
    rule1878 = ReplacementRule(pattern1878, replacement1878)

    if load_rule:
        matcher.add(pattern1878, 1878)
    def With1879(d, c, f, e, x, b, a, g):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + f*(-S(2)*a*b*g + S(3)*c**S(2)*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((a*c*f**S(2) - f*(S(2)*a*b*g + S(3)*c**S(2)*d) + S(4)*g*(a**S(2)*g + b*c*d))/(c*g*(-a*f + S(3)*c*d)), S(2))
        rubi.append(1879)
        return -Simp(c*ArcTan((r - S(2)*x)/q)/(g*q), x) + Simp(c*ArcTan((r + S(2)*x)/q)/(g*q), x) - Simp(c*ArcTan(x*(-a*f + S(3)*c*d)*(S(6)*a**S(2)*b*g**S(2) - S(2)*a**S(2)*c*f*g - a*b**S(2)*f*g + b*c**S(2)*d*f + c**S(2)*g*x**S(4)*(-a*f + S(3)*c*d) + c*x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) - b*c*d*g + S(3)*c**S(2)*d*f))/(g*q*(-S(2)*a**S(2)*g + b*c*d)*(S(4)*a**S(2)*g - a*b*f + b*c*d)))/(g*q), x)
    pattern1879 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1078, cons1079, cons1080, cons1081, cons1082, cons1083)
    rule1879 = ReplacementRule(pattern1879, With1879)

    if load_rule:
        matcher.add(pattern1879, 1879)
    def With1880(d, f, e, x, c, a, g):
        q = Rt((S(12)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        r = Rt((S(4)*a**S(2)*g**S(2) + a*c*f**S(2) - S(3)*c**S(2)*d*f)/(c*g*(-a*f + S(3)*c*d)), S(2))
        rubi.append(1880)
        return -Simp(c*ArcTan((r - S(2)*x)/q)/(g*q), x) + Simp(c*ArcTan((r + S(2)*x)/q)/(g*q), x) - Simp(c*ArcTan(c*x*(-a*f + S(3)*c*d)*(S(2)*a**S(2)*f*g - c*g*x**S(4)*(-a*f + S(3)*c*d) - x**S(2)*(S(2)*a**S(2)*g**S(2) - a*c*f**S(2) + S(3)*c**S(2)*d*f))/(S(8)*a**S(4)*g**S(3)*q))/(g*q), x)
    pattern1880 = Pattern(Integral((x_**S(4)*WC('c', S(1)) + WC('a', S(0)))/(d_ + x_**S(6)*WC('g', S(1)) + x_**S(4)*WC('f', S(1)) + x_**S(2)*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons1084, cons1085, cons1080, cons1086)
    rule1880 = ReplacementRule(pattern1880, With1880)

    if load_rule:
        matcher.add(pattern1880, 1880)

    pattern1881 = Pattern(Integral(u_*v_**p_, x_), cons13, cons137, cons804, cons1017, cons1087, cons1088, cons1089, cons_with_1881)
    def replacement1881(x, v, u, p):

        m = Exponent(u, x)
        n = Exponent(v, x)
        c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
        c = Coefficient(u, x, m)/((m + n*p + S(1))*Coefficient(v, x, n))
        w = Apart(-c*x**(m - n)*(v*(m - n + S(1)) + x*(p + S(1))*D(v, x)) + u, x)
        rubi.append(1881)
        return Simp(If(ZeroQ(w), c*v**(p + 1)*x**(m - n + 1), c*v**(p + 1)*x**(m - n + 1) + Int(v**p*w, x)), x)
    rule1881 = ReplacementRule(pattern1881, replacement1881)

    if load_rule:
        matcher.add(pattern1881, 1881)
    return matcher, [rule1473, rule1474, rule1475, rule1476, rule1477, rule1478, rule1479, rule1480, rule1481, rule1482, rule1483, rule1484, rule1485, rule1486, rule1487, rule1488, rule1489, rule1490, rule1491, rule1492, rule1493, rule1494, rule1495, rule1496, rule1497, rule1498, rule1499, rule1500, rule1501, rule1502, rule1503, rule1504, rule1505, rule1506, rule1507, rule1508, rule1509, rule1510, rule1511, rule1512, rule1513, rule1514, rule1515, rule1516, rule1517, rule1518, rule1519, rule1520, rule1521, rule1522, rule1523, rule1524, rule1525, rule1526, rule1527, rule1528, rule1529, rule1530, rule1531, rule1532, rule1533, rule1534, rule1535, rule1536, rule1537, rule1538, rule1539, rule1540, rule1541, rule1542, rule1543, rule1544, rule1545, rule1546, rule1547, rule1548, rule1549, rule1550, rule1551, rule1552, rule1553, rule1554, rule1555, rule1556, rule1557, rule1558, rule1559, rule1560, rule1561, rule1562, rule1563, rule1564, rule1565, rule1566, rule1567, rule1568, rule1569, rule1570, rule1571, rule1572, rule1573, rule1574, rule1575, rule1576, rule1577, rule1578, rule1579, rule1580, rule1581, rule1582, rule1583, rule1584, rule1585, rule1586, rule1587, rule1588, rule1589, rule1590, rule1591, rule1592, rule1593, rule1594, rule1595, rule1596, rule1597, rule1598, rule1599, rule1600, rule1601, rule1602, rule1603, rule1604, rule1605, rule1606, rule1607, rule1608, rule1609, rule1610, rule1611, rule1612, rule1613, rule1614, rule1615, rule1616, rule1617, rule1618, rule1619, rule1620, rule1621, rule1622, rule1623, rule1624, rule1625, rule1626, rule1627, rule1628, rule1629, rule1630, rule1631, rule1632, rule1633, rule1634, rule1635, rule1636, rule1637, rule1638, rule1639, rule1640, rule1641, rule1642, rule1643, rule1644, rule1645, rule1646, rule1647, rule1648, rule1649, rule1650, rule1651, rule1652, rule1653, rule1654, rule1655, rule1656, rule1657, rule1658, rule1659, rule1660, rule1661, rule1662, rule1663, rule1664, rule1665, rule1666, rule1667, rule1668, rule1669, rule1670, rule1671, rule1672, rule1673, rule1674, rule1675, rule1676, rule1677, rule1678, rule1679, rule1680, rule1681, rule1682, rule1683, rule1684, rule1685, rule1686, rule1687, rule1688, rule1689, rule1690, rule1691, rule1692, rule1693, rule1694, rule1695, rule1696, rule1697, rule1698, rule1699, rule1700, rule1701, rule1702, rule1703, rule1704, rule1705, rule1706, rule1707, rule1708, rule1709, rule1710, rule1711, rule1712, rule1713, rule1714, rule1715, rule1716, rule1717, rule1718, rule1719, rule1720, rule1721, rule1722, rule1723, rule1724, rule1725, rule1726, rule1727, rule1728, rule1729, rule1730, rule1731, rule1732, rule1733, rule1734, rule1735, rule1736, rule1737, rule1738, rule1739, rule1740, rule1741, rule1742, rule1743, rule1744, rule1745, rule1746, rule1747, rule1748, rule1749, rule1750, rule1751, rule1752, rule1753, rule1754, rule1755, rule1756, rule1757, rule1758, rule1759, rule1760, rule1761, rule1762, rule1763, rule1764, rule1765, rule1766, rule1767, rule1768, rule1769, rule1770, rule1771, rule1772, rule1773, rule1774, rule1775, rule1776, rule1777, rule1778, rule1779, rule1780, rule1781, rule1782, rule1783, rule1784, rule1785, rule1786, rule1787, rule1788, rule1789, rule1790, rule1791, rule1792, rule1793, rule1794, rule1795, rule1796, rule1797, rule1798, rule1799, rule1800, rule1801, rule1802, rule1803, rule1804, rule1805, rule1806, rule1807, rule1808, rule1809, rule1810, rule1811, rule1812, rule1813, rule1814, rule1815, rule1816, rule1817, rule1818, rule1819, rule1820, rule1821, rule1822, rule1823, rule1824, rule1825, rule1826, rule1827, rule1828, rule1829, rule1830, rule1831, rule1832, rule1833, rule1834, rule1835, rule1836, rule1837, rule1838, rule1839, rule1840, rule1841, rule1842, rule1843, rule1844, rule1845, rule1846, rule1847, rule1848, rule1849, rule1850, rule1851, rule1852, rule1853, rule1854, rule1855, rule1856, rule1857, rule1858, rule1859, rule1860, rule1861, rule1862, rule1863, rule1864, rule1865, rule1866, rule1867, rule1868, rule1869, rule1870, rule1871, rule1872, rule1873, rule1874, rule1875, rule1876, rule1877, rule1878, rule1879, rule1880, rule1881, ]
