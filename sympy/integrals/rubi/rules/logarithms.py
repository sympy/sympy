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
        Zeta, ProductLog, DerivativeDivides, HypergeometricPFQ, IntHide, OneQ, exp, log
    )
    from sympy import (Integral, S, sqrt, And, Or, Integer, Float, Mod, I, Abs, simplify, Mul, Add, Pow)
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii , Pqq, Q, R, r, C = symbols('i ii Pqq Q R r C_1')
    _UseGamma = False

def logarithms(rubi):
    from sympy.integrals.rubi.constraints import cons1148, cons7, cons27, cons48, cons125, cons5, cons50, cons87, cons88, cons2, cons3, cons1149, cons415, cons1150, cons1151, cons89, cons543, cons1152, cons208, cons209, cons584, cons4, cons66, cons21, cons1153, cons1154, cons1155, cons1156, cons1157, cons1158, cons1159, cons1160, cons1161, cons148, cons1162, cons1163, cons62, cons93, cons168, cons810, cons811, cons222, cons1164, cons224, cons796, cons79, cons1165, cons17, cons1166, cons1167, cons1168, cons1169, cons1170, cons1171, cons1172, cons1173, cons1174, cons1175, cons1176, cons1177, cons1178, cons1179, cons1180, cons1181, cons1182, cons797, cons1183, cons52, cons925, cons1184, cons1185, cons1186, cons1187, cons1188, cons1189, cons1190, cons1191, cons38, cons552, cons1192, cons1193, cons1194, cons25, cons652, cons1195, cons71, cons128, cons1196, cons1197, cons1198, cons1199, cons1200, cons146, cons1201, cons1202, cons13, cons163, cons1203, cons137, cons1204, cons1205, cons1206, cons1207, cons1208, cons1209, cons1210, cons1211, cons1212, cons1213, cons1214, cons70, cons1215, cons1216, cons806, cons840, cons1217, cons1218, cons68, cons1117, cons1219, cons1220, cons1221, cons1222, cons463, cons1223, cons1224, cons1225, cons1226, cons1227, cons1228, cons31, cons1091, cons1229, cons1055, cons515, cons816, cons817, cons1230, cons1231, cons1232, cons1233, cons1234, cons1235, cons1236, cons1237, cons34, cons35, cons1238, cons1239, cons1240

    pattern1987 = Pattern(Integral(log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))), x_), cons7, cons27, cons48, cons125, cons5, cons50, cons1148)
    def replacement1987(e, q, p, f, x, d, c):
        rubi.append(1987)
        return Simp((e + f*x)*log(c*(d*(e + f*x)**p)**q)/f, x) - Simp(p*q*x, x)
    rule1987 = ReplacementRule(pattern1987, replacement1987)
    pattern1988 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons87, cons88)
    def replacement1988(e, a, n, q, p, b, f, x, d, c):
        rubi.append(1988)
        return -Dist(b*n*p*q, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)/f, x)
    rule1988 = ReplacementRule(pattern1988, replacement1988)
    pattern1989 = Pattern(Integral(S(1)/log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('d', S(1))), x_), cons27, cons48, cons125, cons1149)
    def replacement1989(e, f, x, d):
        rubi.append(1989)
        return Simp(LogIntegral(d*(e + f*x))/(d*f), x)
    rule1989 = ReplacementRule(pattern1989, replacement1989)
    pattern1990 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons415)
    def replacement1990(e, a, q, p, b, f, x, d, c):
        rubi.append(1990)
        return Simp((c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*ExpIntegralEi((a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))*exp(-a/(b*p*q))/(b*f*p*q), x)
    rule1990 = ReplacementRule(pattern1990, replacement1990)
    pattern1991 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons1150)
    def replacement1991(e, a, q, p, b, f, x, d, c):
        rubi.append(1991)
        return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*Erfi(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))/Rt(b*p*q, S(2)))*Rt(b*p*q, S(2))*exp(-a/(b*p*q))/(b*f*p*q), x)
    rule1991 = ReplacementRule(pattern1991, replacement1991)
    pattern1992 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons1151)
    def replacement1992(e, a, q, p, b, f, x, d, c):
        rubi.append(1992)
        return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(e + f*x)*Erf(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))/Rt(-b*p*q, S(2)))*Rt(-b*p*q, S(2))*exp(-a/(b*p*q))/(b*f*p*q), x)
    rule1992 = ReplacementRule(pattern1992, replacement1992)
    pattern1993 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons87, cons89)
    def replacement1993(e, a, n, q, p, b, f, x, d, c):
        rubi.append(1993)
        return -Dist(S(1)/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(e + f*x)/(b*f*p*q*(n + S(1))), x)
    rule1993 = ReplacementRule(pattern1993, replacement1993)
    pattern1994 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons543)
    def replacement1994(e, a, n, q, p, b, f, x, d, c):
        rubi.append(1994)
        return Simp((c*(d*(e + f*x)**p)**q)**(-S(1)/(p*q))*(-(a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))**(-n)*(a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)*Gamma(n + S(1), -(a + b*log(c*(d*(e + f*x)**p)**q))/(b*p*q))*exp(-a/(b*p*q))/f, x)
    rule1994 = ReplacementRule(pattern1994, replacement1994)
    pattern1995 = Pattern(Integral(S(1)/((x_*WC('h', S(1)) + WC('g', S(0)))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1152)
    def replacement1995(e, a, h, q, g, p, b, f, x, d, c):
        rubi.append(1995)
        return Simp(log(RemoveContent(a + b*log(c*(d*(e + f*x)**p)**q), x))/(b*h*p*q), x)
    rule1995 = ReplacementRule(pattern1995, replacement1995)
    pattern1996 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons50, cons1152, cons584)
    def replacement1996(e, a, n, h, q, g, p, b, f, x, d, c):
        rubi.append(1996)
        return Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))/(b*h*p*q*(n + S(1))), x)
    rule1996 = ReplacementRule(pattern1996, replacement1996)
    pattern1997 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1152, cons66, cons87, cons88)
    def replacement1997(e, a, n, h, q, m, g, p, b, f, x, d, c):
        rubi.append(1997)
        return -Dist(b*n*p*q/(m + S(1)), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*(g + h*x)**m, x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)
    rule1997 = ReplacementRule(pattern1997, replacement1997)
    pattern1998 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/log((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1))), x_), cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons1153, cons1152, cons1154)
    def replacement1998(e, h, m, g, p, f, x, d):
        rubi.append(1998)
        return Simp((h/f)**(p + S(-1))*LogIntegral(d*(e + f*x)**p)/(d*f*p), x)
    rule1998 = ReplacementRule(pattern1998, replacement1998)
    pattern1999 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_/log((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1))), x_), cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons1153, cons1152, cons1155)
    def replacement1999(e, h, m, g, p, f, x, d):
        rubi.append(1999)
        return Dist((e + f*x)**(-p + S(1))*(g + h*x)**(p + S(-1)), Int((e + f*x)**(p + S(-1))/log(d*(e + f*x)**p), x), x)
    rule1999 = ReplacementRule(pattern1999, replacement1999)
    pattern2000 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1152, cons66)
    def replacement2000(e, a, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2000)
        return Simp((c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*ExpIntegralEi((a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q), x)
    rule2000 = ReplacementRule(pattern2000, replacement2000)
    pattern2001 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1152, cons66, cons1156)
    def replacement2001(e, a, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2001)
        return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*Erfi(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))*Rt((m + S(1))/(b*p*q), S(2)))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q*Rt((m + S(1))/(b*p*q), S(2))), x)
    rule2001 = ReplacementRule(pattern2001, replacement2001)
    pattern2002 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/sqrt(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1152, cons66, cons1157)
    def replacement2002(e, a, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2002)
        return Simp(sqrt(Pi)*(c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(g + h*x)**(m + S(1))*Erf(sqrt(a + b*log(c*(d*(e + f*x)**p)**q))*Rt(-(m + S(1))/(b*p*q), S(2)))*exp(-a*(m + S(1))/(b*p*q))/(b*h*p*q*Rt(-(m + S(1))/(b*p*q), S(2))), x)
    rule2002 = ReplacementRule(pattern2002, replacement2002)
    pattern2003 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1152, cons66, cons87, cons89)
    def replacement2003(e, a, h, n, q, m, g, p, b, f, x, d, c):
        rubi.append(2003)
        return -Dist((m + S(1))/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**m, x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**(m + S(1))/(b*h*p*q*(n + S(1))), x)
    rule2003 = ReplacementRule(pattern2003, replacement2003)
    pattern2004 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons1152, cons66)
    def replacement2004(e, a, n, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2004)
        return Simp((c*(d*(e + f*x)**p)**q)**(-(m + S(1))/(p*q))*(-(a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))**(-n)*(a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))*Gamma(n + S(1), -(a + b*log(c*(d*(e + f*x)**p)**q))*(m + S(1))/(b*p*q))*exp(-a*(m + S(1))/(b*p*q))/(h*(m + S(1))), x)
    rule2004 = ReplacementRule(pattern2004, replacement2004)
    pattern2005 = Pattern(Integral(log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons7, cons48, cons125, cons208, cons209, cons1158)
    def replacement2005(e, h, g, f, x, c):
        rubi.append(2005)
        return -Simp(PolyLog(S(2), -(g + h*x)*Together(c*f/h))/h, x)
    rule2005 = ReplacementRule(pattern2005, replacement2005)
    pattern2006 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1))))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons209, cons1159, cons1160)
    def replacement2006(e, a, h, g, b, f, x, c):
        rubi.append(2006)
        return Dist(b, Int(log(-h*(e + f*x)/(-e*h + f*g))/(g + h*x), x), x) + Simp((a + b*log(c*(e - f*g/h)))*log(g + h*x)/h, x)
    rule2006 = ReplacementRule(pattern2006, replacement2006)
    pattern2007 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1161, cons148)
    def replacement2007(e, a, n, h, q, g, p, b, f, x, d, c):
        rubi.append(2007)
        return -Dist(b*f*n*p*q/h, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*log(f*(g + h*x)/(-e*h + f*g))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*log(f*(g + h*x)/(-e*h + f*g))/h, x)
    rule2007 = ReplacementRule(pattern2007, replacement2007)
    pattern2008 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1161, cons66)
    def replacement2008(e, a, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2008)
        return -Dist(b*f*p*q/(h*(m + S(1))), Int((g + h*x)**(m + S(1))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)
    rule2008 = ReplacementRule(pattern2008, replacement2008)
    pattern2009 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_/(x_*WC('h', S(1)) + WC('g', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1161, cons87, cons88)
    def replacement2009(e, a, h, n, q, g, p, b, f, x, d, c):
        rubi.append(2009)
        return -Dist(b*f*n*p*q/(-e*h + f*g), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))/(g + h*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(e + f*x)/((g + h*x)*(-e*h + f*g)), x)
    rule2009 = ReplacementRule(pattern2009, replacement2009)
    pattern2010 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons5, cons50, cons1161, cons87, cons88, cons66, cons1162, cons1163)
    def replacement2010(e, a, h, n, q, m, g, p, b, f, x, d, c):
        rubi.append(2010)
        return -Dist(b*f*n*p*q/(h*(m + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*(g + h*x)**(m + S(1))/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)
    rule2010 = ReplacementRule(pattern2010, replacement2010)
    pattern2011 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1161, cons62)
    def replacement2011(e, a, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2011)
        return Int(ExpandIntegrand((g + h*x)**m/(a + b*log(c*(d*(e + f*x)**p)**q)), x), x)
    rule2011 = ReplacementRule(pattern2011, replacement2011)
    pattern2012 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1161, cons93, cons89, cons168)
    def replacement2012(e, a, h, n, q, m, g, p, b, f, x, d, c):
        rubi.append(2012)
        return -Dist((m + S(1))/(b*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**m, x), x) + Dist(m*(-e*h + f*g)/(b*f*p*q*(n + S(1))), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(g + h*x)**(m + S(-1)), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(1))*(e + f*x)*(g + h*x)**m/(b*f*p*q*(n + S(1))), x)
    rule2012 = ReplacementRule(pattern2012, replacement2012)
    pattern2013 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons50, cons1161, cons62)
    def replacement2013(e, a, h, n, q, m, g, p, b, f, x, d, c):
        rubi.append(2013)
        return Int(ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x), x)
    rule2013 = ReplacementRule(pattern2013, replacement2013)
    pattern2014 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((v_**p_*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons810, cons811)
    def replacement2014(a, n, q, m, v, p, b, u, x, d, c):
        rubi.append(2014)
        return Int((a + b*log(c*(d*ExpandToSum(v, x)**p)**q))**n*ExpandToSum(u, x)**m, x)
    rule2014 = ReplacementRule(pattern2014, replacement2014)
    pattern2015 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons222)
    def replacement2015(e, a, n, h, q, m, g, p, b, f, x, d, c):
        rubi.append(2015)
        return Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x)
    rule2015 = ReplacementRule(pattern2015, replacement2015)
    pattern2016 = Pattern(Integral(log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0))))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons7, cons48, cons125, cons208, cons209, cons224, cons796, cons1152, cons1164)
    def replacement2016(e, h, g, j, f, x, i, c):
        rubi.append(2016)
        return Simp(f*PolyLog(S(2), f*(i + j*x)/(j*(e + f*x)))/(h*(-e*j + f*i)), x)
    rule2016 = ReplacementRule(pattern2016, replacement2016)
    pattern2017 = Pattern(Integral((a_ + WC('b', S(1))*log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0)))))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons209, cons224, cons796, cons1152, cons1164)
    def replacement2017(e, a, h, g, b, j, f, x, i, c):
        rubi.append(2017)
        return Dist(a, Int(S(1)/((g + h*x)*(i + j*x)), x), x) + Dist(b, Int(log(c/(e + f*x))/((g + h*x)*(i + j*x)), x), x)
    rule2017 = ReplacementRule(pattern2017, replacement2017)
    def With2018(e, a, h, q, i, m, g, p, b, j, f, x, d, c):
        u = IntHide((i + j*x)**m/(g + h*x), x)
        rubi.append(2018)
        return -Dist(b*h*p*q, Int(SimplifyIntegrand(u/(g + h*x), x), x), x) + Dist(a + b*log(c*(d*(e + f*x)**p)**q), u, x)
    pattern2018 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons5, cons50, cons1152, cons79)
    rule2018 = ReplacementRule(pattern2018, With2018)
    pattern2019 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((x_*WC('f', S(1)) + WC('e', S(0)))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons48, cons125, cons208, cons209, cons224, cons796, cons4, cons1152, cons62, cons1165)
    def replacement2019(e, a, n, h, m, g, b, j, f, x, i, c):
        rubi.append(2019)
        return Dist(c**(-m)*f**(-m)/h, Subst(Int((a + b*x)**n*(-c*e*j + c*f*i + j*exp(x))**m, x), x, log(c*(e + f*x))), x)
    rule2019 = ReplacementRule(pattern2019, replacement2019)
    def With2020(e, a, n, h, q, i, m, g, p, b, j, f, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, (i + j*x)**m/(g + h*x), x)
        if SumQ(u):
            return True
        return False
    pattern2020 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons5, cons50, cons17, cons148, CustomConstraint(With2020))
    def replacement2020(e, a, n, h, q, i, m, g, p, b, j, f, x, d, c):

        u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, (i + j*x)**m/(g + h*x), x)
        rubi.append(2020)
        return Int(u, x)
    rule2020 = ReplacementRule(pattern2020, replacement2020)
    pattern2021 = Pattern(Integral((x_*WC('j', S(1)) + WC('i', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons21, cons4, cons5, cons50, cons1166)
    def replacement2021(e, a, n, h, q, i, m, g, p, b, j, f, x, d, c):
        rubi.append(2021)
        return Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(i + j*x)**m/(g + h*x), x)
    rule2021 = ReplacementRule(pattern2021, replacement2021)
    pattern2022 = Pattern(Integral(log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0))))/(g_ + x_**S(2)*WC('h', S(1))), x_), cons7, cons48, cons125, cons208, cons209, cons1167, cons1168)
    def replacement2022(e, h, g, f, x, c):
        rubi.append(2022)
        return -Simp(f*PolyLog(S(2), (-e + f*x)/(e + f*x))/(S(2)*e*h), x)
    rule2022 = ReplacementRule(pattern2022, replacement2022)
    pattern2023 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(WC('c', S(1))/(x_*WC('f', S(1)) + WC('e', S(0)))))/(g_ + x_**S(2)*WC('h', S(1))), x_), cons7, cons48, cons125, cons208, cons209, cons1167, cons1169, cons1170)
    def replacement2023(e, a, h, g, b, f, x, c):
        rubi.append(2023)
        return Dist(b, Int(log(S(2)*e/(e + f*x))/(g + h*x**S(2)), x), x) + Dist(a + b*log(c/(S(2)*e)), Int(S(1)/(g + h*x**S(2)), x), x)
    rule2023 = ReplacementRule(pattern2023, replacement2023)
    pattern2024 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(x_**S(2)*WC('i', S(1)) + x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons5, cons50, cons1171)
    def replacement2024(e, a, h, q, i, g, p, b, f, x, d, c):
        rubi.append(2024)
        return Dist(e*f, Int((a + b*log(c*(d*(e + f*x)**p)**q))/((e + f*x)*(e*i*x + f*g)), x), x)
    rule2024 = ReplacementRule(pattern2024, replacement2024)
    pattern2025 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((e_ + x_*WC('f', S(1)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(g_ + x_**S(2)*WC('i', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons224, cons5, cons50, cons1172)
    def replacement2025(e, a, q, i, g, p, b, f, x, d, c):
        rubi.append(2025)
        return Dist(e*f, Int((a + b*log(c*(d*(e + f*x)**p)**q))/((e + f*x)*(e*i*x + f*g)), x), x)
    rule2025 = ReplacementRule(pattern2025, replacement2025)
    def With2026(e, a, h, q, g, p, b, f, x, d, c):
        u = IntHide(S(1)/sqrt(g + h*x**S(2)), x)
        rubi.append(2026)
        return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Simp(u*(a + b*log(c*(d*(e + f*x)**p)**q)), x)
    pattern2026 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/sqrt(g_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1173)
    rule2026 = ReplacementRule(pattern2026, With2026)
    def With2027(e, h1, a, g1, q, h2, g2, p, b, f, x, d, c):
        u = IntHide(S(1)/sqrt(g1*g2 + h1*h2*x**S(2)), x)
        rubi.append(2027)
        return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Simp(u*(a + b*log(c*(d*(e + f*x)**p)**q)), x)
    pattern2027 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(sqrt(g1_ + x_*WC('h1', S(1)))*sqrt(g2_ + x_*WC('h2', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons1177, cons1178, cons1179, cons1180, cons5, cons50, cons1174, cons1175, cons1176)
    rule2027 = ReplacementRule(pattern2027, With2027)
    pattern2028 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/sqrt(g_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons1181)
    def replacement2028(e, a, h, q, g, p, b, f, x, d, c):
        rubi.append(2028)
        return Dist(sqrt(S(1) + h*x**S(2)/g)/sqrt(g + h*x**S(2)), Int((a + b*log(c*(d*(e + f*x)**p)**q))/sqrt(S(1) + h*x**S(2)/g), x), x)
    rule2028 = ReplacementRule(pattern2028, replacement2028)
    pattern2029 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))/(sqrt(g1_ + x_*WC('h1', S(1)))*sqrt(g2_ + x_*WC('h2', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons1177, cons1178, cons1179, cons1180, cons5, cons50, cons1174)
    def replacement2029(e, h1, a, g1, q, h2, g2, p, b, f, x, d, c):
        rubi.append(2029)
        return Dist(sqrt(S(1) + h1*h2*x**S(2)/(g1*g2))/(sqrt(g1 + h1*x)*sqrt(g2 + h2*x)), Int((a + b*log(c*(d*(e + f*x)**p)**q))/sqrt(S(1) + h1*h2*x**S(2)/(g1*g2)), x), x)
    rule2029 = ReplacementRule(pattern2029, replacement2029)
    pattern2030 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('k', S(1)) + WC('j', S(0)))*WC('i', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons797, cons5, cons50, cons87, cons88, cons1182)
    def replacement2030(e, a, n, h, q, g, p, d, j, k, b, f, x, i, c):
        rubi.append(2030)
        return Dist(b*f*n*p*q/h, Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(S(2), Together(-i*(j + k*x) + S(1)))/(e + f*x), x), x) - Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(S(2), Together(-i*(j + k*x) + S(1)))/h, x)
    rule2030 = ReplacementRule(pattern2030, replacement2030)
    pattern2031 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('k', S(1)) + WC('j', S(0)))**WC('m', S(1))*WC('i', S(1)) + S(1))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons797, cons21, cons5, cons50, cons87, cons88, cons1183)
    def replacement2031(e, a, n, h, q, m, g, p, d, j, k, b, f, x, i, c):
        rubi.append(2031)
        return Dist(b*f*n*p*q/(h*m), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(S(2), -i*(j + k*x)**m)/(e + f*x), x), x) - Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(S(2), -i*(j + k*x)**m)/(h*m), x)
    rule2031 = ReplacementRule(pattern2031, replacement2031)
    pattern2032 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*PolyLog(r_, (x_*WC('k', S(1)) + WC('j', S(0)))**WC('m', S(1))*WC('i', S(1)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons797, cons21, cons5, cons50, cons52, cons87, cons88, cons1183)
    def replacement2032(e, a, n, h, q, i, m, g, p, r, b, j, k, f, x, d, c):
        rubi.append(2032)
        return -Dist(b*f*n*p*q/(h*m), Int((a + b*log(c*(d*(e + f*x)**p)**q))**(n + S(-1))*PolyLog(r + S(1), i*(j + k*x)**m)/(e + f*x), x), x) + Simp((a + b*log(c*(d*(e + f*x)**p)**q))**n*PolyLog(r + S(1), i*(j + k*x)**m)/(h*m), x)
    rule2032 = ReplacementRule(pattern2032, replacement2032)
    def With2033(e, a, h, q, m, g, p, b, F, f, x, d, Px, c):
        u = IntHide(Px*F(g + h*x)**m, x)
        rubi.append(2033)
        return -Dist(b*f*p*q, Int(SimplifyIntegrand(u/(e + f*x), x), x), x) + Dist(a + b*log(c*(d*(e + f*x)**p)**q), u, x)
    pattern2033 = Pattern(Integral(F_**(x_*WC('h', S(1)) + WC('g', S(0)))*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))*WC('Px', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons50, cons925, cons62, cons1184)
    rule2033 = ReplacementRule(pattern2033, With2033)
    pattern2034 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((e_ + x_**m_*WC('f', S(1)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons148)
    def replacement2034(e, a, n, q, m, p, b, f, x, d, c):
        rubi.append(2034)
        return Dist(S(1)/m, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n/x, x), x, x**m), x)
    rule2034 = ReplacementRule(pattern2034, replacement2034)
    pattern2035 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_**m_*(f_ + x_**WC('r', S(1))*WC('e', S(1))))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons1185, cons148)
    def replacement2035(e, a, n, q, m, p, b, r, f, x, d, c):
        rubi.append(2035)
        return Dist(S(1)/m, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n/x, x), x, x**m), x)
    rule2035 = ReplacementRule(pattern2035, replacement2035)
    pattern2036 = Pattern(Integral(x_**WC('r1', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_**r_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons52, cons1186)
    def replacement2036(e, a, n, q, p, r1, b, r, f, x, d, c):
        rubi.append(2036)
        return Dist(S(1)/r, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n, x), x, x**r), x)
    rule2036 = ReplacementRule(pattern2036, replacement2036)
    pattern2037 = Pattern(Integral(x_**WC('r1', S(1))*(x_**r_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(((x_**r_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons52, cons1186)
    def replacement2037(e, a, n, h, q, m, g, p, r1, b, r, f, x, d, c):
        rubi.append(2037)
        return Dist(S(1)/r, Subst(Int((a + b*log(c*(d*(e + f*x)**p)**q))**n*(g + h*x)**m, x), x, x**r), x)
    rule2037 = ReplacementRule(pattern2037, replacement2037)
    def With2038(e, a, n, b, x, d, c):
        u = IntHide(S(1)/(d + e*x**S(2)), x)
        rubi.append(2038)
        return -Dist(b*n, Int(u/x, x), x) + Dist(a + b*log(c*x**n), u, x)
    pattern2038 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1187)
    rule2038 = ReplacementRule(pattern2038, With2038)
    pattern2039 = Pattern(Integral(log((x_**mn_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))/(x_*(d_ + x_**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1188, cons1189)
    def replacement2039(e, a, n, mn, b, x, d, c):
        rubi.append(2039)
        return Simp(PolyLog(S(2), -Together(b*c*x**(-n)*(d + e*x**n)/d))/(d*n), x)
    rule2039 = ReplacementRule(pattern2039, replacement2039)
    pattern2040 = Pattern(Integral(log(x_**mn_*(x_**WC('n', S(1))*WC('a', S(1)) + WC('b', S(0)))*WC('c', S(1)))/(x_*(d_ + x_**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1188, cons1189)
    def replacement2040(e, a, n, mn, b, x, d, c):
        rubi.append(2040)
        return Simp(PolyLog(S(2), -Together(b*c*x**(-n)*(d + e*x**n)/d))/(d*n), x)
    rule2040 = ReplacementRule(pattern2040, replacement2040)
    pattern2041 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons925)
    def replacement2041(e, a, n, q, p, b, f, x, d, Px, c):
        rubi.append(2041)
        return Int(ExpandIntegrand(Px*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x), x)
    rule2041 = ReplacementRule(pattern2041, replacement2041)
    def With2042(e, a, n, q, RFx, p, b, f, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, RFx, x)
        if SumQ(u):
            return True
        return False
    pattern2042 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons1190, cons148, CustomConstraint(With2042))
    def replacement2042(e, a, n, q, RFx, p, b, f, x, d, c):

        u = ExpandIntegrand((a + b*log(c*(d*(e + f*x)**p)**q))**n, RFx, x)
        rubi.append(2042)
        return Int(u, x)
    rule2042 = ReplacementRule(pattern2042, replacement2042)
    def With2043(e, a, n, q, RFx, p, b, f, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand(RFx*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
        if SumQ(u):
            return True
        return False
    pattern2043 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons1190, cons148, CustomConstraint(With2043))
    def replacement2043(e, a, n, q, RFx, p, b, f, x, d, c):

        u = ExpandIntegrand(RFx*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
        rubi.append(2043)
        return Int(u, x)
    rule2043 = ReplacementRule(pattern2043, replacement2043)
    pattern2044 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_**S(2)*WC('g', S(1)) + x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons50, cons4, cons1191, cons38)
    def replacement2044(e, a, n, q, u, g, p, b, f, x, d, c):
        rubi.append(2044)
        return Int(u*(a + b*log(c*(S(4)**(-p)*d*g**(-p)*(f + S(2)*g*x)**(S(2)*p))**q))**n, x)
    rule2044 = ReplacementRule(pattern2044, replacement2044)
    pattern2045 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((v_**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons552, cons1192)
    def replacement2045(a, n, q, v, p, b, u, x, d, c):
        rubi.append(2045)
        return Int(u*(a + b*log(c*(d*ExpandToSum(v, x)**p)**q))**n, x)
    rule2045 = ReplacementRule(pattern2045, replacement2045)
    pattern2046 = Pattern(Integral(log(((x_**WC('n', S(1))*WC('c', S(1)))**p_*WC('b', S(1)))**q_*WC('a', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons50, cons52, cons1193)
    def replacement2046(a, n, p, b, r, x, q, c):
        rubi.append(2046)
        return Subst(Int(log(x**(n*p*q))**r, x), x**(n*p*q), a*(b*(c*x**n)**p)**q)
    rule2046 = ReplacementRule(pattern2046, replacement2046)
    pattern2047 = Pattern(Integral(x_**WC('m', S(1))*log(((x_**WC('n', S(1))*WC('c', S(1)))**p_*WC('b', S(1)))**q_*WC('a', S(1)))**WC('r', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons50, cons52, cons66, cons1194)
    def replacement2047(a, n, m, p, b, r, x, q, c):
        rubi.append(2047)
        return Subst(Int(x**m*log(x**(n*p*q))**r, x), x**(n*p*q), a*(b*(c*x**n)**p)**q)
    rule2047 = ReplacementRule(pattern2047, replacement2047)
    pattern2048 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e1', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons652, cons25)
    def replacement2048(e, e1, a, n, p, b, u, x, d, c):
        rubi.append(2048)
        return Dist(log(e*(b*e1/d)**n)**p, Int(u, x), x)
    rule2048 = ReplacementRule(pattern2048, replacement2048)
    pattern2049 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons652, cons1196, cons1195, cons71, cons128)
    def replacement2049(e, n1, a, n, n2, p, b, c, x, d, e1):
        rubi.append(2049)
        return -Dist(n*n1*p*(-a*d + b*c)/b, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/(c + d*x), x), x) + Simp((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/b, x)
    rule2049 = ReplacementRule(pattern2049, replacement2049)
    pattern2050 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons71, cons1197, cons1198)
    def replacement2050(e, a, g, b, f, x, d, c):
        rubi.append(2050)
        return Simp(PolyLog(S(2), Together(-a*e + c)/(c + d*x))/g, x)
    rule2050 = ReplacementRule(pattern2050, replacement2050)
    pattern2051 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1197, cons128)
    def replacement2051(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2051)
        return Dist(n*n1*p*(-a*d + b*c)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log((-a*d + b*c)/(b*(c + d*x)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log((-a*d + b*c)/(b*(c + d*x)))/g, x)
    rule2051 = ReplacementRule(pattern2051, replacement2051)
    pattern2052 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1199, cons128)
    def replacement2052(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2052)
        return Dist(n*n1*p*(-a*d + b*c)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log(-(-a*d + b*c)/(d*(a + b*x)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log(-(-a*d + b*c)/(d*(a + b*x)))/g, x)
    rule2052 = ReplacementRule(pattern2052, replacement2052)
    pattern2053 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1200)
    def replacement2053(e, n1, a, n, n2, g, b, f, c, x, d, e1):
        rubi.append(2053)
        return -Dist(n*n1*(-a*d + b*c)/g, Int(log(f + g*x)/((a + b*x)*(c + d*x)), x), x) + Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)*log(f + g*x)/g, x)
    rule2053 = ReplacementRule(pattern2053, replacement2053)
    pattern2054 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1200, cons38, cons146)
    def replacement2054(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2054)
        return Dist(d/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(c + d*x), x), x) - Dist((-c*g + d*f)/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c + d*x)*(f + g*x)), x), x)
    rule2054 = ReplacementRule(pattern2054, replacement2054)
    pattern2055 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons71, cons1197)
    def replacement2055(e, a, g, b, f, x, d, c):
        rubi.append(2055)
        return Simp(d**S(2)*LogIntegral(e*(a + b*x)/(c + d*x))/(e*g**S(2)*(-a*d + b*c)), x)
    rule2055 = ReplacementRule(pattern2055, replacement2055)
    pattern2056 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1197)
    def replacement2056(e, n1, a, n, n2, g, b, f, c, x, d, e1):
        rubi.append(2056)
        return Simp(d**S(2)*(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**(-S(1)/(n*n1))*(a + b*x)*ExpIntegralEi(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)/(n*n1))/(g**S(2)*n*n1*(c + d*x)*(-a*d + b*c)), x)
    rule2056 = ReplacementRule(pattern2056, replacement2056)
    pattern2057 = Pattern(Integral(S(1)/((x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1199)
    def replacement2057(e, n1, a, n, n2, g, b, f, c, x, d, e1):
        rubi.append(2057)
        return Simp(b**S(2)*(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**(S(1)/(n*n1))*(c + d*x)*ExpIntegralEi(-log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)/(n*n1))/(g**S(2)*n*n1*(a + b*x)*(-a*d + b*c)), x)
    rule2057 = ReplacementRule(pattern2057, replacement2057)
    pattern2058 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1201, cons128)
    def replacement2058(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2058)
        return -Dist(n*n1*p*(-a*d + b*c)/(-a*g + b*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((c + d*x)*(f + g*x)), x), x) + Simp((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((f + g*x)*(-a*g + b*f)), x)
    rule2058 = ReplacementRule(pattern2058, replacement2058)
    pattern2059 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1200, cons128)
    def replacement2059(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2059)
        return -Dist(n*n1*p*(-a*d + b*c)/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((a + b*x)*(f + g*x)), x), x) + Simp((c + d*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((f + g*x)*(-c*g + d*f)), x)
    rule2059 = ReplacementRule(pattern2059, replacement2059)
    pattern2060 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0)))**S(3), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1201, cons1197)
    def replacement2060(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2060)
        return Dist(b/(-a*g + b*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(2), x), x) - Dist(g/(-a*g + b*f), Int((a + b*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(3), x), x)
    rule2060 = ReplacementRule(pattern2060, replacement2060)
    pattern2061 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**p_/(x_*WC('g', S(1)) + WC('f', S(0)))**S(3), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1200, cons1199)
    def replacement2061(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2061)
        return Dist(d/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(2), x), x) - Dist(g/(-c*g + d*f), Int((c + d*x)*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(f + g*x)**S(3), x), x)
    rule2061 = ReplacementRule(pattern2061, replacement2061)
    pattern2062 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons128, cons17, cons66)
    def replacement2062(e, n1, a, n, m, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2062)
        return -Dist(n*n1*p*(-a*d + b*c)/(g*(m + S(1))), Int((f + g*x)**(m + S(1))*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp((f + g*x)**(m + S(1))*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(g*(m + S(1))), x)
    rule2062 = ReplacementRule(pattern2062, replacement2062)
    pattern2063 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_**n_*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1203, cons1202, cons71, cons66, cons13, cons163)
    def replacement2063(e, a, m2, n, m, p, b, u, x, d, c):
        rubi.append(2063)
        return -Dist(n*p/(m + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(e*u**n)**(p + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(e*u**n)**p/((m + S(1))*(-a*d + b*c)), x)
    rule2063 = ReplacementRule(pattern2063, replacement2063)
    pattern2064 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_)**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons1203, cons1202, cons71, cons66, cons13, cons163)
    def replacement2064(a, m2, m, p, b, u, x, d, c):
        rubi.append(2064)
        return -Dist(p/(m + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(u)**(p + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(u)**p/((m + S(1))*(-a*d + b*c)), x)
    rule2064 = ReplacementRule(pattern2064, replacement2064)
    pattern2065 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))/log(u_**n_*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1203, cons1202, cons71, cons66)
    def replacement2065(e, a, m2, n, m, b, u, x, d, c):
        rubi.append(2065)
        return Simp((e*u**n)**(-(m + S(1))/n)*(a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*ExpIntegralEi((m + S(1))*log(e*u**n)/n)/(n*(-a*d + b*c)), x)
    rule2065 = ReplacementRule(pattern2065, replacement2065)
    pattern2066 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))/log(u_), x_), cons2, cons3, cons7, cons27, cons1203, cons1202, cons71, cons66)
    def replacement2066(a, m2, m, b, u, x, d, c):
        rubi.append(2066)
        return Simp(u**(-m + S(-1))*(a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*ExpIntegralEi((m + S(1))*log(u))/(-a*d + b*c), x)
    rule2066 = ReplacementRule(pattern2066, replacement2066)
    pattern2067 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_**n_*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1203, cons1202, cons71, cons66, cons13, cons137)
    def replacement2067(e, a, m2, n, m, p, b, u, x, d, c):
        rubi.append(2067)
        return -Dist((m + S(1))/(n*(p + S(1))), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(e*u**n)**(p + S(1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)
    rule2067 = ReplacementRule(pattern2067, replacement2067)
    pattern2068 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m2', S(1))*log(u_)**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons1203, cons1202, cons71, cons66, cons13, cons137)
    def replacement2068(a, m2, m, p, b, u, x, d, c):
        rubi.append(2068)
        return -Dist((m + S(1))/(p + S(1)), Int((a + b*x)**m*(c + d*x)**(-m + S(-2))*log(u)**(p + S(1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(-m + S(-1))*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)
    rule2068 = ReplacementRule(pattern2068, replacement2068)
    pattern2069 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1201, cons1197)
    def replacement2069(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2069)
        return Dist(d/g, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/(c + d*x)**S(2), x), x)
    rule2069 = ReplacementRule(pattern2069, replacement2069)
    pattern2070 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons71, cons1201, cons1200, cons1204)
    def replacement2070(e, a, g, b, f, x, d, c):
        rubi.append(2070)
        return Simp(PolyLog(S(2), -(f + g*x)*(a*e - c)/(f*(c + d*x)))/(-c*g + d*f), x)
    rule2070 = ReplacementRule(pattern2070, replacement2070)
    pattern2071 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons652, cons1196, cons1195, cons71, cons1201, cons1200, cons128)
    def replacement2071(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2071)
        return Dist(n*n1*p*(-a*d + b*c)/(-c*g + d*f), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**(p + S(-1))*log((f + g*x)*(-a*d + b*c)/((c + d*x)*(-a*g + b*f)))/((a + b*x)*(c + d*x)), x), x) - Simp(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p*log((f + g*x)*(-a*d + b*c)/((c + d*x)*(-a*g + b*f)))/(-c*g + d*f), x)
    rule2071 = ReplacementRule(pattern2071, replacement2071)
    pattern2072 = Pattern(Integral(log((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(f_ + x_**S(2)*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1205, cons1206)
    def replacement2072(e, a, g, b, f, x, d, c):
        rubi.append(2072)
        return Simp(c*PolyLog(S(2), -(c - d*x)*(a*e - c)/(c*(c + d*x)))/(S(2)*d*f), x)
    rule2072 = ReplacementRule(pattern2072, replacement2072)
    pattern2073 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1207)
    def replacement2073(e, n1, a, n, h, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2073)
        return Dist(d**S(2), Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c + d*x)*(-c*h + d*g + d*h*x)), x), x)
    rule2073 = ReplacementRule(pattern2073, replacement2073)
    pattern2074 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons209, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1208)
    def replacement2074(e, n1, a, n, h, n2, p, b, f, c, x, d, e1):
        rubi.append(2074)
        return -Dist(d**S(2)/h, Int(log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((c - d*x)*(c + d*x)), x), x)
    rule2074 = ReplacementRule(pattern2074, replacement2074)
    pattern2075 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/((x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1200, cons1199)
    def replacement2075(e, n1, a, n, n2, g, p, b, f, c, x, d, e1):
        rubi.append(2075)
        return Dist(b/(g*n*n1*(-a*d + b*c)), Subst(Int(x**p, x), x, log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)), x)
    rule2075 = ReplacementRule(pattern2075, replacement2075)
    pattern2076 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1209, cons1203, cons71, cons13, cons163)
    def replacement2076(e, a, n, v, p, b, u, x, d, c):
        rubi.append(2076)
        return Dist(n*p, Int(PolyLog(S(2), Together(-v + S(1)))*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(S(2), Together(-v + S(1)))*log(e*u**n)**p/(-a*d + b*c), x)
    rule2076 = ReplacementRule(pattern2076, replacement2076)
    pattern2077 = Pattern(Integral(log(u_)**WC('p', S(1))*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons1209, cons1203, cons71, cons13, cons163)
    def replacement2077(a, v, p, b, u, x, d, c):
        rubi.append(2077)
        return Dist(p, Int(PolyLog(S(2), Together(-v + S(1)))*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(S(2), Together(-v + S(1)))*log(u)**p/(-a*d + b*c), x)
    rule2077 = ReplacementRule(pattern2077, replacement2077)
    def With2078(e, a, n, v, p, b, u, x, d, c):
        f = (-v + S(1))/u
        rubi.append(2078)
        return Dist(f/(n*(p + S(1))), Int(log(e*u**n)**(p + S(1))/((c + d*x)*(-a*f - b*f + c + d)), x), x) + Simp(log(v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)
    pattern2078 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1209, cons1203, cons71, cons13, cons137)
    rule2078 = ReplacementRule(pattern2078, With2078)
    def With2079(a, v, p, b, u, x, d, c):
        f = (-v + S(1))/u
        rubi.append(2079)
        return Dist(f/(p + S(1)), Int(log(u)**(p + S(1))/((c + d*x)*(-a*f - b*f + c + d)), x), x) + Simp(log(u)**(p + S(1))*log(v)/((p + S(1))*(-a*d + b*c)), x)
    pattern2079 = Pattern(Integral(log(u_)**p_*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons1209, cons1203, cons71, cons13, cons137)
    rule2079 = ReplacementRule(pattern2079, With2079)
    pattern2080 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1210, cons1203, cons71, cons13, cons163)
    def replacement2080(e, a, n, v, p, b, u, x, d, c):
        rubi.append(2080)
        return -Dist(n*p, Int(PolyLog(S(2), Together(-v + S(1)))*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(S(2), Together(-v + S(1)))*log(e*u**n)**p/(-a*d + b*c), x)
    rule2080 = ReplacementRule(pattern2080, replacement2080)
    pattern2081 = Pattern(Integral(log(u_)**WC('p', S(1))*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons1210, cons1203, cons71, cons13, cons163)
    def replacement2081(a, v, p, b, u, x, d, c):
        rubi.append(2081)
        return -Dist(p, Int(PolyLog(S(2), Together(-v + S(1)))*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(S(2), Together(-v + S(1)))*log(u)**p/(-a*d + b*c), x)
    rule2081 = ReplacementRule(pattern2081, replacement2081)
    def With2082(e, a, n, v, p, b, u, x, d, c):
        f = u*(-v + S(1))
        rubi.append(2082)
        return -Dist(f/(n*(p + S(1))), Int(log(e*u**n)**(p + S(1))/((a + b*x)*(a + b - c*f - d*f)), x), x) + Simp(log(v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)
    pattern2082 = Pattern(Integral(log(v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons1210, cons1203, cons71, cons13, cons137)
    rule2082 = ReplacementRule(pattern2082, With2082)
    def With2083(a, v, p, b, u, x, d, c):
        f = u*(-v + S(1))
        rubi.append(2083)
        return -Dist(f/(p + S(1)), Int(log(u)**(p + S(1))/((a + b*x)*(a + b - c*f - d*f)), x), x) + Simp(log(u)**(p + S(1))*log(v)/((p + S(1))*(-a*d + b*c)), x)
    pattern2083 = Pattern(Integral(log(u_)**p_*log(v_)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons1210, cons1203, cons71, cons13, cons137)
    rule2083 = ReplacementRule(pattern2083, With2083)
    pattern2084 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons1211, cons1203, cons71, cons13, cons146)
    def replacement2084(e, a, n, q, v, p, b, u, x, d, c):
        rubi.append(2084)
        return -Dist(n*p, Int(PolyLog(q + S(1), v)*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q + S(1), v)*log(e*u**n)**p/(-a*d + b*c), x)
    rule2084 = ReplacementRule(pattern2084, replacement2084)
    pattern2085 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons50, cons1211, cons1203, cons71, cons13, cons146)
    def replacement2085(a, q, v, p, b, u, x, d, c):
        rubi.append(2085)
        return -Dist(p, Int(PolyLog(q + S(1), v)*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q + S(1), v)*log(u)**p/(-a*d + b*c), x)
    rule2085 = ReplacementRule(pattern2085, replacement2085)
    pattern2086 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons1211, cons1203, cons71, cons13, cons137)
    def replacement2086(e, a, n, q, v, p, b, u, x, d, c):
        rubi.append(2086)
        return -Dist(S(1)/(n*(p + S(1))), Int(PolyLog(q + S(-1), v)*log(e*u**n)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)
    rule2086 = ReplacementRule(pattern2086, replacement2086)
    pattern2087 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons50, cons1211, cons1203, cons71, cons13, cons137)
    def replacement2087(a, q, v, p, b, u, x, d, c):
        rubi.append(2087)
        return -Dist(S(1)/(p + S(1)), Int(PolyLog(q + S(-1), v)*log(u)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)
    rule2087 = ReplacementRule(pattern2087, replacement2087)
    pattern2088 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons1212, cons1203, cons71, cons13, cons146)
    def replacement2088(e, a, n, q, v, p, b, u, x, d, c):
        rubi.append(2088)
        return Dist(n*p, Int(PolyLog(q + S(1), v)*log(e*u**n)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(q + S(1), v)*log(e*u**n)**p/(-a*d + b*c), x)
    rule2088 = ReplacementRule(pattern2088, replacement2088)
    pattern2089 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons50, cons1212, cons1203, cons71, cons13, cons146)
    def replacement2089(a, q, v, p, b, u, x, d, c):
        rubi.append(2089)
        return Dist(p, Int(PolyLog(q + S(1), v)*log(u)**(p + S(-1))/((a + b*x)*(c + d*x)), x), x) - Simp(PolyLog(q + S(1), v)*log(u)**p/(-a*d + b*c), x)
    rule2089 = ReplacementRule(pattern2089, replacement2089)
    pattern2090 = Pattern(Integral(PolyLog(q_, v_)*log(u_**n_*WC('e', S(1)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons1212, cons1203, cons71, cons13, cons137)
    def replacement2090(e, a, n, q, v, p, b, u, x, d, c):
        rubi.append(2090)
        return Dist(S(1)/(n*(p + S(1))), Int(PolyLog(q + S(-1), v)*log(e*u**n)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(e*u**n)**(p + S(1))/(n*(p + S(1))*(-a*d + b*c)), x)
    rule2090 = ReplacementRule(pattern2090, replacement2090)
    pattern2091 = Pattern(Integral(PolyLog(q_, v_)*log(u_)**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons50, cons1212, cons1203, cons71, cons13, cons137)
    def replacement2091(a, q, v, p, b, u, x, d, c):
        rubi.append(2091)
        return Dist(S(1)/(p + S(1)), Int(PolyLog(q + S(-1), v)*log(u)**(p + S(1))/((a + b*x)*(c + d*x)), x), x) + Simp(PolyLog(q, v)*log(u)**(p + S(1))/((p + S(1))*(-a*d + b*c)), x)
    rule2091 = ReplacementRule(pattern2091, replacement2091)
    pattern2092 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1213, cons1214)
    def replacement2092(e, n1, a, n, h, n2, f, g, p, b, u, c, x, d, e1):
        rubi.append(2092)
        return Dist(b*d/h, Int(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((a + b*x)*(c + d*x)), x), x)
    rule2092 = ReplacementRule(pattern2092, replacement2092)
    pattern2093 = Pattern(Integral(WC('u', S(1))*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1))/(x_**S(2)*WC('h', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons209, cons4, cons5, cons652, cons1196, cons1195, cons71, cons1213, cons70)
    def replacement2093(e, n1, a, n, h, n2, f, p, b, u, c, x, d, e1):
        rubi.append(2093)
        return Dist(b*d/h, Int(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**(-n1))**n)**p/((a + b*x)*(c + d*x)), x), x)
    rule2093 = ReplacementRule(pattern2093, replacement2093)
    def With2094(e, n1, a, n, h, n2, g, b, f, c, x, d, e1):
        u = IntHide(S(1)/(f + g*x + h*x**S(2)), x)
        rubi.append(2094)
        return -Dist(n*(-a*d + b*c), Int(u/((a + b*x)*(c + d*x)), x), x) + Simp(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n), x)
    pattern2094 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(f_ + x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons652, cons125, cons208, cons209, cons4, cons1196, cons1195)
    rule2094 = ReplacementRule(pattern2094, With2094)
    def With2095(e, n1, a, n, h, n2, b, f, c, x, d, e1):
        u = IntHide(S(1)/(f + h*x**S(2)), x)
        rubi.append(2095)
        return -Dist(n*(-a*d + b*c), Int(u/((a + b*x)*(c + d*x)), x), x) + Simp(u*log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n), x)
    pattern2095 = Pattern(Integral(log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))/(f_ + x_**S(2)*WC('h', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons652, cons125, cons209, cons4, cons1196, cons1195)
    rule2095 = ReplacementRule(pattern2095, With2095)
    def With2096(e, n1, a, n, RFx, n2, p, b, c, x, d, e1):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand(log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**p, RFx, x)
        if SumQ(u):
            return True
        return False
    pattern2096 = Pattern(Integral(RFx_*log(((x_*WC('b', S(1)) + WC('a', S(0)))**WC('n1', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**n2_*WC('e1', S(1)))**WC('n', S(1))*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons652, cons1196, cons1195, cons1190, cons128, CustomConstraint(With2096))
    def replacement2096(e, n1, a, n, RFx, n2, p, b, c, x, d, e1):

        u = ExpandIntegrand(log(e*(e1*(a + b*x)**n1*(c + d*x)**n2)**n)**p, RFx, x)
        rubi.append(2096)
        return Int(u, x)
    rule2096 = ReplacementRule(pattern2096, replacement2096)
    def With2097(p, u, x, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = QuotientOfLinearsParts(v, x)
        if Not(And(OneQ(p), ZeroQ(Part(lst, S(3))))):
            return True
        return False
    pattern2097 = Pattern(Integral(WC('u', S(1))*log(v_)**WC('p', S(1)), x_), cons5, cons1215, cons1216, CustomConstraint(With2097))
    def replacement2097(p, u, x, v):

        lst = QuotientOfLinearsParts(v, x)
        rubi.append(2097)
        return Int(u*log((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**p, x)
    rule2097 = ReplacementRule(pattern2097, replacement2097)
    pattern2098 = Pattern(Integral(log((x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons4, cons5, cons806)
    def replacement2098(a, n, p, b, x, c):
        rubi.append(2098)
        return -Dist(b*n*p, Int(x**n/(a + b*x**n), x), x) + Simp(x*log(c*(a + b*x**n)**p), x)
    rule2098 = ReplacementRule(pattern2098, replacement2098)
    pattern2099 = Pattern(Integral(log(v_**WC('p', S(1))*WC('c', S(1))), x_), cons7, cons5, cons840, cons1217)
    def replacement2099(p, x, v, c):
        rubi.append(2099)
        return Int(log(c*ExpandToSum(v, x)**p), x)
    rule2099 = ReplacementRule(pattern2099, replacement2099)
    pattern2100 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))/(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons1218)
    def replacement2100(e, a, n, g, p, b, f, x, d, c):
        rubi.append(2100)
        return -Dist(b*e*n*p/g, Int(x**(n + S(-1))*log(f + g*x)/(d + e*x**n), x), x) + Simp((a + b*log(c*(d + e*x**n)**p))*log(f + g*x)/g, x)
    rule2100 = ReplacementRule(pattern2100, replacement2100)
    pattern2101 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons66)
    def replacement2101(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(2101)
        return -Dist(b*e*n*p/(g*(m + S(1))), Int(x**(n + S(-1))*(f + g*x)**(m + S(1))/(d + e*x**n), x), x) + Simp((a + b*log(c*(d + e*x**n)**p))*(f + g*x)**(m + S(1))/(g*(m + S(1))), x)
    rule2101 = ReplacementRule(pattern2101, replacement2101)
    pattern2102 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(v_**WC('p', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons21, cons5, cons68, cons840, cons1117)
    def replacement2102(a, m, v, p, b, u, x, c):
        rubi.append(2102)
        return Int((a + b*log(c*ExpandToSum(v, x)**p))*ExpandToSum(u, x)**m, x)
    rule2102 = ReplacementRule(pattern2102, replacement2102)
    def With2103(e, a, n, m, g, p, b, f, x, d, c):
        w = IntHide(asin(f + g*x)**m, x)
        rubi.append(2103)
        return -Dist(b*e*n*p, Int(SimplifyIntegrand(w*x**(n + S(-1))/(d + e*x**n), x), x), x) + Dist(a + b*log(c*(d + e*x**n)**p), w, x)
    pattern2103 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**n_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))*asin(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons62)
    rule2103 = ReplacementRule(pattern2103, With2103)
    def With2104(e, a, g, p, b, f, x, d, c):
        u = IntHide(S(1)/(f + g*x**S(2)), x)
        rubi.append(2104)
        return -Dist(S(2)*b*e*p, Int(u*x/(d + e*x**S(2)), x), x) + Simp(u*(a + b*log(c*(d + e*x**S(2))**p)), x)
    pattern2104 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*WC('c', S(1))))/(x_**S(2)*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons1219)
    rule2104 = ReplacementRule(pattern2104, With2104)
    pattern2105 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons148)
    def replacement2105(e, a, n, p, b, x, d, c):
        rubi.append(2105)
        return -Dist(S(2)*b*e*n*p, Int(x**S(2)*(a + b*log(c*(d + e*x**S(2))**p))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x*(a + b*log(c*(d + e*x**S(2))**p))**n, x)
    rule2105 = ReplacementRule(pattern2105, replacement2105)
    pattern2106 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons148, cons1220)
    def replacement2106(e, a, n, m, p, b, x, d, c):
        rubi.append(2106)
        return Dist(S(1)/2, Subst(Int(x**(m/S(2) + S(-1)/2)*(a + b*log(c*(d + e*x)**p))**n, x), x, x**S(2)), x)
    rule2106 = ReplacementRule(pattern2106, replacement2106)
    pattern2107 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*WC('c', S(1))))**n_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons148, cons1221)
    def replacement2107(e, a, n, m, p, b, x, d, c):
        rubi.append(2107)
        return -Dist(S(2)*b*e*n*p/(m + S(1)), Int(x**(m + S(2))*(a + b*log(c*(d + e*x**S(2))**p))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*log(c*(d + e*x**S(2))**p))**n/(m + S(1)), x)
    rule2107 = ReplacementRule(pattern2107, replacement2107)
    def With2108(u, x, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        w = DerivativeDivides(v, u*(-v + S(1)), x)
        if Not(FalseQ(w)):
            return True
        return False
    pattern2108 = Pattern(Integral(u_*log(v_), x_), CustomConstraint(With2108))
    def replacement2108(u, x, v):

        w = DerivativeDivides(v, u*(-v + S(1)), x)
        rubi.append(2108)
        return Simp(w*PolyLog(S(2), Together(-v + S(1))), x)
    rule2108 = ReplacementRule(pattern2108, replacement2108)
    def With2109(a, u, v, b, w, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = DerivativeDivides(v, w*(-v + S(1)), x)
        if Not(FalseQ(z)):
            return True
        return False
    pattern2109 = Pattern(Integral(w_*(WC('a', S(0)) + WC('b', S(1))*log(u_))*log(v_), x_), cons2, cons3, cons1222, CustomConstraint(With2109))
    def replacement2109(a, u, v, b, w, x):

        z = DerivativeDivides(v, w*(-v + S(1)), x)
        rubi.append(2109)
        return -Dist(b, Int(SimplifyIntegrand(z*D(u, x)*PolyLog(S(2), Together(-v + S(1)))/u, x), x), x) + Simp(z*(a + b*log(u))*PolyLog(S(2), Together(-v + S(1))), x)
    rule2109 = ReplacementRule(pattern2109, replacement2109)
    pattern2110 = Pattern(Integral(log((a_ + (x_*WC('e', S(1)) + WC('d', S(0)))**n_*WC('b', S(1)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons87, cons463)
    def replacement2110(e, a, n, p, b, x, d, c):
        rubi.append(2110)
        return -Dist(b*n*p, Int(S(1)/(a*(d + e*x)**(-n) + b), x), x) + Simp((d + e*x)*log(c*(a + b*(d + e*x)**n)**p)/e, x)
    rule2110 = ReplacementRule(pattern2110, replacement2110)
    pattern2111 = Pattern(Integral(log((a_ + (x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons1223)
    def replacement2111(e, n, a, p, b, x, d, c):
        rubi.append(2111)
        return Dist(a*n*p, Int(S(1)/(a + b*(d + e*x)**n), x), x) + Simp((d + e*x)*log(c*(a + b*(d + e*x)**n)**p)/e, x) - Simp(n*p*x, x)
    rule2111 = ReplacementRule(pattern2111, replacement2111)
    pattern2112 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log((d_ + WC('e', S(1))/(x_*WC('g', S(1)) + WC('f', S(0))))**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons148)
    def replacement2112(e, a, n, g, p, b, f, x, d, c):
        rubi.append(2112)
        return -Dist(b*e*n*p/(d*g), Subst(Int((a + b*log(c*(d + e*x)**p))**(n + S(-1))/x, x), x, S(1)/(f + g*x)), x) + Simp((a + b*log(c*(d + e/(f + g*x))**p))**n*(d*(f + g*x) + e)/(d*g), x)
    rule2112 = ReplacementRule(pattern2112, replacement2112)
    pattern2113 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons5, cons1190, cons148)
    def replacement2113(a, n, p, b, x, RFx, c):
        rubi.append(2113)
        return -Dist(b*n*p, Int(SimplifyIntegrand(x*(a + b*log(RFx**p*c))**(n + S(-1))*D(RFx, x)/RFx, x), x), x) + Simp(x*(a + b*log(RFx**p*c))**n, x)
    rule2113 = ReplacementRule(pattern2113, replacement2113)
    pattern2114 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons1190, cons148)
    def replacement2114(e, x, a, n, p, b, RFx, d, c):
        rubi.append(2114)
        return -Dist(b*n*p/e, Int((a + b*log(RFx**p*c))**(n + S(-1))*D(RFx, x)*log(d + e*x)/RFx, x), x) + Simp((a + b*log(RFx**p*c))**n*log(d + e*x)/e, x)
    rule2114 = ReplacementRule(pattern2114, replacement2114)
    pattern2115 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons1190, cons148, cons1224, cons66)
    def replacement2115(e, x, a, n, m, p, b, RFx, d, c):
        rubi.append(2115)
        return -Dist(b*n*p/(e*(m + S(1))), Int(SimplifyIntegrand((a + b*log(RFx**p*c))**(n + S(-1))*(d + e*x)**(m + S(1))*D(RFx, x)/RFx, x), x), x) + Simp((a + b*log(RFx**p*c))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)
    rule2115 = ReplacementRule(pattern2115, replacement2115)
    def With2116(e, x, n, RFx, d, c):
        u = IntHide(S(1)/(d + e*x**S(2)), x)
        rubi.append(2116)
        return -Dist(n, Int(SimplifyIntegrand(u*D(RFx, x)/RFx, x), x), x) + Simp(u*log(RFx**n*c), x)
    pattern2116 = Pattern(Integral(log(RFx_**WC('n', S(1))*WC('c', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons7, cons27, cons48, cons4, cons1190, cons1225)
    rule2116 = ReplacementRule(pattern2116, With2116)
    def With2117(n, Qx, x, Px, c):
        u = IntHide(S(1)/Qx, x)
        rubi.append(2117)
        return -Dist(n, Int(SimplifyIntegrand(u*D(Px, x)/Px, x), x), x) + Simp(u*log(Px**n*c), x)
    pattern2117 = Pattern(Integral(log(Px_**WC('n', S(1))*WC('c', S(1)))/Qx_, x_), cons7, cons4, cons1226, cons1227)
    rule2117 = ReplacementRule(pattern2117, With2117)
    def With2118(a, n, p, b, RGx, x, RFx, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand((a + b*log(RFx**p*c))**n, RGx, x)
        if SumQ(u):
            return True
        return False
    pattern2118 = Pattern(Integral(RGx_*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons5, cons1190, cons1228, cons148, CustomConstraint(With2118))
    def replacement2118(a, n, p, b, RGx, x, RFx, c):

        u = ExpandIntegrand((a + b*log(RFx**p*c))**n, RGx, x)
        rubi.append(2118)
        return Int(u, x)
    rule2118 = ReplacementRule(pattern2118, replacement2118)
    def With2119(a, n, p, b, RGx, x, RFx, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        u = ExpandIntegrand(RGx*(a + b*log(RFx**p*c))**n, x)
        if SumQ(u):
            return True
        return False
    pattern2119 = Pattern(Integral(RGx_*(WC('a', S(0)) + WC('b', S(1))*log(RFx_**WC('p', S(1))*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons5, cons1190, cons1228, cons148, CustomConstraint(With2119))
    def replacement2119(a, n, p, b, RGx, x, RFx, c):

        u = ExpandIntegrand(RGx*(a + b*log(RFx**p*c))**n, x)
        rubi.append(2119)
        return Int(u, x)
    rule2119 = ReplacementRule(pattern2119, replacement2119)
    def With2120(a, b, u, x, RFx):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfLinear(RFx*(a + b*log(u)), x)
        try:
            res = Not(FalseQ(lst))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern2120 = Pattern(Integral(RFx_*(WC('a', S(0)) + WC('b', S(1))*log(u_)), x_), cons2, cons3, cons1190, CustomConstraint(With2120))
    def replacement2120(a, b, u, x, RFx):

        lst = SubstForFractionalPowerOfLinear(RFx*(a + b*log(u)), x)
        rubi.append(2120)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule2120 = ReplacementRule(pattern2120, replacement2120)
    pattern2121 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log((F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1))*WC('e', S(1)) + S(1)), x_), cons1091, cons2, cons3, cons7, cons48, cons125, cons208, cons4, cons31, cons168)
    def replacement2121(e, a, n, m, g, b, F, f, x, c):
        rubi.append(2121)
        return Dist(g*m/(b*c*n*log(F)), Int((f + g*x)**(m + S(-1))*PolyLog(S(2), -e*(F**(c*(a + b*x)))**n), x), x) - Simp((f + g*x)**m*PolyLog(S(2), -e*(F**(c*(a + b*x)))**n)/(b*c*n*log(F)), x)
    rule2121 = ReplacementRule(pattern2121, replacement2121)
    pattern2122 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*log(d_ + (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1))*WC('e', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons31, cons168, cons1229)
    def replacement2122(e, a, n, m, g, b, F, f, x, d, c):
        rubi.append(2122)
        return Int((f + g*x)**m*log(S(1) + e*(F**(c*(a + b*x)))**n/d), x) - Simp((f + g*x)**(m + S(1))*log(S(1) + e*(F**(c*(a + b*x)))**n/d)/(g*(m + S(1))), x) + Simp((f + g*x)**(m + S(1))*log(d + e*(F**(c*(a + b*x)))**n)/(g*(m + S(1))), x)
    rule2122 = ReplacementRule(pattern2122, replacement2122)
    pattern2123 = Pattern(Integral(log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons1055)
    def replacement2123(e, a, b, f, x, d, c):
        rubi.append(2123)
        return Dist(f**S(2)*(-S(4)*a*c + b**S(2))/S(2), Int(x/(-f*sqrt(a + b*x + c*x**S(2))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d)) + (-b*f**S(2) + S(2)*d*e)*(a + b*x + c*x**S(2))), x), x) + Simp(x*log(d + e*x + f*sqrt(a + b*x + c*x**S(2))), x)
    rule2123 = ReplacementRule(pattern2123, replacement2123)
    pattern2124 = Pattern(Integral(log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons1055)
    def replacement2124(e, a, f, x, d, c):
        rubi.append(2124)
        return -Dist(a*c*f**S(2), Int(x/(d*e*(a + c*x**S(2)) + f*sqrt(a + c*x**S(2))*(a*e - c*d*x)), x), x) + Simp(x*log(d + e*x + f*sqrt(a + c*x**S(2))), x)
    rule2124 = ReplacementRule(pattern2124, replacement2124)
    pattern2125 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons1055, cons66, cons515)
    def replacement2125(e, a, m, g, b, f, x, d, c):
        rubi.append(2125)
        return Dist(f**S(2)*(-S(4)*a*c + b**S(2))/(S(2)*g*(m + S(1))), Int((g*x)**(m + S(1))/(-f*sqrt(a + b*x + c*x**S(2))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d)) + (-b*f**S(2) + S(2)*d*e)*(a + b*x + c*x**S(2))), x), x) + Simp((g*x)**(m + S(1))*log(d + e*x + f*sqrt(a + b*x + c*x**S(2)))/(g*(m + S(1))), x)
    rule2125 = ReplacementRule(pattern2125, replacement2125)
    pattern2126 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*log(x_*WC('e', S(1)) + sqrt(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons21, cons1055, cons66, cons515)
    def replacement2126(e, a, m, g, f, x, d, c):
        rubi.append(2126)
        return -Dist(a*c*f**S(2)/(g*(m + S(1))), Int((g*x)**(m + S(1))/(d*e*(a + c*x**S(2)) + f*sqrt(a + c*x**S(2))*(a*e - c*d*x)), x), x) + Simp((g*x)**(m + S(1))*log(d + e*x + f*sqrt(a + c*x**S(2)))/(g*(m + S(1))), x)
    rule2126 = ReplacementRule(pattern2126, replacement2126)
    pattern2127 = Pattern(Integral(WC('v', S(1))*log(sqrt(u_)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons27, cons48, cons125, cons816, cons817, cons1230)
    def replacement2127(e, u, v, f, x, d):
        rubi.append(2127)
        return Int(v*log(d + e*x + f*sqrt(ExpandToSum(u, x))), x)
    rule2127 = ReplacementRule(pattern2127, replacement2127)
    pattern2128 = Pattern(Integral(log(u_), x_), cons1222)
    def replacement2128(u, x):
        rubi.append(2128)
        return -Int(SimplifyIntegrand(x*D(u, x)/u, x), x) + Simp(x*log(u), x)
    rule2128 = ReplacementRule(pattern2128, replacement2128)
    pattern2129 = Pattern(Integral(log(u_)/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons1231, cons1232)
    def replacement2129(u, a, x, b):
        rubi.append(2129)
        return -Dist(S(1)/b, Int(SimplifyIntegrand(D(u, x)*log(a + b*x)/u, x), x), x) + Simp(log(u)*log(a + b*x)/b, x)
    rule2129 = ReplacementRule(pattern2129, replacement2129)
    pattern2130 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*log(u_), x_), cons2, cons3, cons21, cons1222, cons66)
    def replacement2130(a, m, b, u, x):
        rubi.append(2130)
        return -Dist(S(1)/(b*(m + S(1))), Int(SimplifyIntegrand((a + b*x)**(m + S(1))*D(u, x)/u, x), x), x) + Simp((a + b*x)**(m + S(1))*log(u)/(b*(m + S(1))), x)
    rule2130 = ReplacementRule(pattern2130, replacement2130)
    def With2131(Qx, u, x):
        v = IntHide(S(1)/Qx, x)
        rubi.append(2131)
        return -Int(SimplifyIntegrand(v*D(u, x)/u, x), x) + Simp(v*log(u), x)
    pattern2131 = Pattern(Integral(log(u_)/Qx_, x_), cons1233, cons1222)
    rule2131 = ReplacementRule(pattern2131, With2131)
    pattern2132 = Pattern(Integral(u_**(x_*WC('a', S(1)))*log(u_), x_), cons2, cons1222)
    def replacement2132(u, a, x):
        rubi.append(2132)
        return -Int(SimplifyIntegrand(u**(a*x + S(-1))*x*D(u, x), x), x) + Simp(u**(a*x)/a, x)
    rule2132 = ReplacementRule(pattern2132, replacement2132)
    def With2133(u, x, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        w = IntHide(v, x)
        if InverseFunctionFreeQ(w, x):
            return True
        return False
    pattern2133 = Pattern(Integral(v_*log(u_), x_), cons1222, CustomConstraint(With2133))
    def replacement2133(u, x, v):

        w = IntHide(v, x)
        rubi.append(2133)
        return Dist(log(u), w, x) - Int(SimplifyIntegrand(w*D(u, x)/u, x), x)
    rule2133 = ReplacementRule(pattern2133, replacement2133)
    pattern2134 = Pattern(Integral(log(v_)*log(w_), x_), cons1234, cons1235)
    def replacement2134(w, x, v):
        rubi.append(2134)
        return -Int(SimplifyIntegrand(x*D(v, x)*log(w)/v, x), x) - Int(SimplifyIntegrand(x*D(w, x)*log(v)/w, x), x) + Simp(x*log(v)*log(w), x)
    rule2134 = ReplacementRule(pattern2134, replacement2134)
    def With2135(w, x, u, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = IntHide(u, x)
        if InverseFunctionFreeQ(z, x):
            return True
        return False
    pattern2135 = Pattern(Integral(u_*log(v_)*log(w_), x_), cons1234, cons1235, CustomConstraint(With2135))
    def replacement2135(w, x, u, v):

        z = IntHide(u, x)
        rubi.append(2135)
        return Dist(log(v)*log(w), z, x) - Int(SimplifyIntegrand(z*D(v, x)*log(w)/v, x), x) - Int(SimplifyIntegrand(z*D(w, x)*log(v)/w, x), x)
    rule2135 = ReplacementRule(pattern2135, replacement2135)
    pattern2136 = Pattern(Integral(log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons4, cons5, cons1236)
    def replacement2136(n, a, p, b, x):
        rubi.append(2136)
        return -Dist(n*p, Int(S(1)/log(b*x**n), x), x) + Simp(x*log(a*log(b*x**n)**p), x)
    rule2136 = ReplacementRule(pattern2136, replacement2136)
    pattern2137 = Pattern(Integral(log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1)))/x_, x_), cons2, cons3, cons4, cons5, cons1236)
    def replacement2137(n, a, p, b, x):
        rubi.append(2137)
        return Simp((-p + log(a*log(b*x**n)**p))*log(b*x**n)/n, x)
    rule2137 = ReplacementRule(pattern2137, replacement2137)
    pattern2138 = Pattern(Integral(x_**WC('m', S(1))*log(WC('a', S(1))*log(x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons21, cons4, cons5, cons66)
    def replacement2138(n, a, m, p, b, x):
        rubi.append(2138)
        return -Dist(n*p/(m + S(1)), Int(x**m/log(b*x**n), x), x) + Simp(x**(m + S(1))*log(a*log(b*x**n)**p)/(m + S(1)), x)
    rule2138 = ReplacementRule(pattern2138, replacement2138)
    pattern2139 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*log(x_*WC('d', S(1)) + WC('c', S(0))))/sqrt(a_ + WC('b', S(1))*log(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons34, cons35, cons1237)
    def replacement2139(a, B, b, A, x, d, c):
        rubi.append(2139)
        return Dist(B/b, Int(sqrt(a + b*log(c + d*x)), x), x) + Dist((A*b - B*a)/b, Int(S(1)/sqrt(a + b*log(c + d*x)), x), x)
    rule2139 = ReplacementRule(pattern2139, replacement2139)
    pattern2140 = Pattern(Integral(f_**(WC('a', S(1))*log(u_)), x_), cons2, cons125, cons1238)
    def replacement2140(f, a, x, u):
        rubi.append(2140)
        return Int(u**(a*log(f)), x)
    rule2140 = ReplacementRule(pattern2140, replacement2140)
    def With2141(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfLog(u*x, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern2141 = Pattern(Integral(u_, x_), cons1239, CustomConstraint(With2141))
    def replacement2141(u, x):

        lst = FunctionOfLog(u*x, x)
        rubi.append(2141)
        return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, log(Part(lst, S(2)))), x)
    rule2141 = ReplacementRule(pattern2141, replacement2141)
    pattern2142 = Pattern(Integral(WC('u', S(1))*log(Gamma(v_)), x_))
    def replacement2142(u, x, v):
        rubi.append(2142)
        return Dist(-LogGamma(v) + log(Gamma(v)), Int(u, x), x) + Int(u*LogGamma(v), x)
    rule2142 = ReplacementRule(pattern2142, replacement2142)
    pattern2143 = Pattern(Integral((w_*WC('a', S(1)) + w_*WC('b', S(1))*log(v_)**WC('n', S(1)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons38)
    def replacement2143(w, a, n, v, p, b, u, x):
        rubi.append(2143)
        return Int(u*w**p*(a + b*log(v)**n)**p, x)
    rule2143 = ReplacementRule(pattern2143, replacement2143)
    pattern2144 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*log(((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('c', S(1))))**n_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons50, cons1240)
    def replacement2144(e, a, n, q, u, p, b, f, x, d, c):
        rubi.append(2144)
        return Int(u*(a + b*log(c*(d*(e + f*x)**p)**q))**n, x)
    rule2144 = ReplacementRule(pattern2144, replacement2144)
    return [rule1987, rule1988, rule1989, rule1990, rule1991, rule1992, rule1993, rule1994, rule1995, rule1996, rule1997, rule1998, rule1999, rule2000, rule2001, rule2002, rule2003, rule2004, rule2005, rule2006, rule2007, rule2008, rule2009, rule2010, rule2011, rule2012, rule2013, rule2014, rule2015, rule2016, rule2017, rule2018, rule2019, rule2020, rule2021, rule2022, rule2023, rule2024, rule2025, rule2026, rule2027, rule2028, rule2029, rule2030, rule2031, rule2032, rule2033, rule2034, rule2035, rule2036, rule2037, rule2038, rule2039, rule2040, rule2041, rule2042, rule2043, rule2044, rule2045, rule2046, rule2047, rule2048, rule2049, rule2050, rule2051, rule2052, rule2053, rule2054, rule2055, rule2056, rule2057, rule2058, rule2059, rule2060, rule2061, rule2062, rule2063, rule2064, rule2065, rule2066, rule2067, rule2068, rule2069, rule2070, rule2071, rule2072, rule2073, rule2074, rule2075, rule2076, rule2077, rule2078, rule2079, rule2080, rule2081, rule2082, rule2083, rule2084, rule2085, rule2086, rule2087, rule2088, rule2089, rule2090, rule2091, rule2092, rule2093, rule2094, rule2095, rule2096, rule2097, rule2098, rule2099, rule2100, rule2101, rule2102, rule2103, rule2104, rule2105, rule2106, rule2107, rule2108, rule2109, rule2110, rule2111, rule2112, rule2113, rule2114, rule2115, rule2116, rule2117, rule2118, rule2119, rule2120, rule2121, rule2122, rule2123, rule2124, rule2125, rule2126, rule2127, rule2128, rule2129, rule2130, rule2131, rule2132, rule2133, rule2134, rule2135, rule2136, rule2137, rule2138, rule2139, rule2140, rule2141, rule2142, rule2143, rule2144, ]
