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

def trinomial_products(rubi, matcher, load_rule):
    from sympy.integrals.rubi.constraints import cons46, cons87, cons463, cons38, cons2, cons3, cons7, cons489, cons5, cons45, cons147, cons664, cons4, cons665, cons584, cons666, cons13, cons163, cons667, cons314, cons668, cons462, cons196, cons669, cons670, cons146, cons671, cons672, cons338, cons137, cons226, cons128, cons246, cons673, cons674, cons413, cons675, cons293, cons676, cons677, cons484, cons177, cons_with_1095, cons_with_1097, cons_with_1098, cons_with_1099, cons_with_1100, cons678, cons679, cons680, cons585, cons681, cons68, cons69, cons53, cons21, cons501, cons27, cons63, cons502, cons682, cons155, cons683, cons684, cons225, cons56, cons243, cons148, cons244, cons685, cons17, cons686, cons687, cons510, cons688, cons689, cons690, cons529, cons31, cons530, cons691, cons94, cons367, cons356, cons500, cons_with_1134, cons692, cons693, cons694, cons695, cons696, cons697, cons698, cons699, cons700, cons701, cons_with_1155, cons_with_1156, cons_with_1157, cons_with_1158, cons541, cons23, cons702, cons552, cons553, cons554, cons220, cons48, cons50, cons703, cons704, cons256, cons257, cons279, cons221, cons280, cons395, cons396, cons705, cons706, cons707, cons708, cons709, cons85, cons710, cons711, cons712, cons713, cons586, cons386, cons149, cons714, cons43, cons_with_1221, cons_with_1222, cons_with_1223, cons_with_1224, cons_with_1225, cons_with_1226, cons_with_1227, cons_with_1228, cons715, cons_with_1230, cons_with_1231, cons_with_1232, cons_with_1233, cons_with_1234, cons_with_1235, cons_with_1236, cons_with_1237, cons448, cons716, cons400, cons717, cons718, cons_with_1256, cons_with_1257, cons719, cons347, cons564, cons720, cons268, cons721, cons722, cons723, cons724, cons725, cons726, cons727, cons728, cons125, cons208, cons52, cons_with_1293, cons_with_1294, cons593, cons729, cons730, cons731, cons652, cons732, cons654, cons34, cons35, cons733, cons18, cons734, cons464, cons735, cons168, cons267, cons736, cons737, cons_with_1327, cons_with_1328, cons738, cons739, cons740, cons81, cons_with_1343, cons434, cons_with_1344, cons741, cons742, cons743, cons744, cons_with_1347, cons_with_1348, cons611, cons403, cons745, cons746, cons747, cons748, cons749, cons750, cons751, cons752, cons753, cons754, cons755, cons756, cons757, cons758, cons759, cons760, cons606, cons761, cons762, cons763, cons764, cons765, cons766, cons767, cons768, cons769, cons770, cons771, cons772, cons773, cons774, cons775, cons776, cons777, cons778, cons779, cons780, cons781, cons782, cons_with_1452, cons_with_1454, cons783, cons784, cons785, cons_with_1459, cons786, cons_with_1461, cons787, cons_with_1463, cons788, cons_with_1465, cons789, cons_with_1467, cons790, cons_with_1469, cons791, cons792, cons793, cons794, cons795, cons796, cons797

    pattern1076 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons46, cons87, cons463, cons38)
    def replacement1076(p, n2, x, b, c, a, n):
        rubi.append(1076)
        return Int(x**(S(2)*n*p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1076 = ReplacementRule(pattern1076, replacement1076)

    if load_rule:
        matcher.add(pattern1076, 1076)
    def With1077(p, n2, x, b, c, a, n):
        k = Denominator(n)
        rubi.append(1077)
        return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*x**(k*n) + c*x**(S(2)*k*n))**p, x), x, x**(S(1)/k)), x)
    pattern1077 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons489)
    rule1077 = ReplacementRule(pattern1077, With1077)

    if load_rule:
        matcher.add(pattern1077, 1077)
    pattern1078 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons664)
    def replacement1078(p, n2, x, b, c, a, n):
        rubi.append(1078)
        return Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a), x)
    rule1078 = ReplacementRule(pattern1078, replacement1078)

    if load_rule:
        matcher.add(pattern1078, 1078)
    pattern1079 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons665, cons584)
    def replacement1079(p, n2, x, b, c, a, n):
        rubi.append(1079)
        return -Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*(S(2)*p + S(1))), x) + Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*(n + S(1))), x)
    rule1079 = ReplacementRule(pattern1079, replacement1079)

    if load_rule:
        matcher.add(pattern1079, 1079)
    pattern1080 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons45, cons666, cons13, cons163, cons667)
    def replacement1080(p, n2, x, b, c, a, n):
        rubi.append(1080)
        return Dist(sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int((b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)/2), x), x)
    rule1080 = ReplacementRule(pattern1080, replacement1080)

    if load_rule:
        matcher.add(pattern1080, 1080)
    pattern1081 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons45, cons666, cons13, cons163, cons314)
    def replacement1081(p, n2, x, b, c, a, n):
        rubi.append(1081)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1081 = ReplacementRule(pattern1081, replacement1081)

    if load_rule:
        matcher.add(pattern1081, 1081)
    pattern1082 = Pattern(Integral(sqrt(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons584, cons668, cons462)
    def replacement1082(c, n2, x, b, a, n):
        rubi.append(1082)
        return Simp(x*sqrt(a + b*x**n + c*x**(S(2)*n))/(n + S(1)), x) + Simp(b*n*x*sqrt(a + b*x**n + c*x**(S(2)*n))/((b + S(2)*c*x**n)*(n + S(1))), x)
    rule1082 = ReplacementRule(pattern1082, replacement1082)

    if load_rule:
        matcher.add(pattern1082, 1082)
    pattern1083 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons45, cons147, cons196)
    def replacement1083(p, n2, x, b, c, a, n):
        rubi.append(1083)
        return -Subst(Int((a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)
    rule1083 = ReplacementRule(pattern1083, replacement1083)

    if load_rule:
        matcher.add(pattern1083, 1083)
    pattern1084 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons45, cons147, cons669, cons670, cons13, cons146)
    def replacement1084(p, n2, x, b, c, a, n):
        rubi.append(1084)
        return Dist(S(2)*a*n**S(2)*p*(S(2)*p + S(-1))/((S(2)*n*p + S(1))*(n*(S(2)*p + S(-1)) + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p + S(1)), x) + Simp(n*p*x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/((S(2)*n*p + S(1))*(n*(S(2)*p + S(-1)) + S(1))), x)
    rule1084 = ReplacementRule(pattern1084, replacement1084)

    if load_rule:
        matcher.add(pattern1084, 1084)
    pattern1085 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons45, cons147, cons671, cons672, cons338, cons137)
    def replacement1085(p, n2, x, b, c, a, n):
        rubi.append(1085)
        return Dist((S(2)*n*(p + S(1)) + S(1))*(n*(S(2)*p + S(1)) + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*n*(S(2)*p + S(1))), x) - Simp(x*(n*(S(2)*p + S(1)) + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)
    rule1085 = ReplacementRule(pattern1085, replacement1085)

    if load_rule:
        matcher.add(pattern1085, 1085)
    pattern1086 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147)
    def replacement1086(p, n2, x, b, c, a, n):
        rubi.append(1086)
        return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b/S(2) + c*x**n)**(S(2)*p), x), x)
    rule1086 = ReplacementRule(pattern1086, replacement1086)

    if load_rule:
        matcher.add(pattern1086, 1086)
    pattern1087 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons46, cons196)
    def replacement1087(p, n2, x, b, c, a, n):
        rubi.append(1087)
        return -Subst(Int((a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)
    rule1087 = ReplacementRule(pattern1087, replacement1087)

    if load_rule:
        matcher.add(pattern1087, 1087)
    pattern1088 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons226, cons128)
    def replacement1088(p, n2, x, b, c, a, n):
        rubi.append(1088)
        return Int(ExpandIntegrand((a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1088 = ReplacementRule(pattern1088, replacement1088)

    if load_rule:
        matcher.add(pattern1088, 1088)
    pattern1089 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons226, cons13, cons163, cons669, cons246, cons673)
    def replacement1089(p, n2, x, b, c, a, n):
        rubi.append(1089)
        return Dist(n*p/(S(2)*n*p + S(1)), Int((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p + S(1)), x)
    rule1089 = ReplacementRule(pattern1089, replacement1089)

    if load_rule:
        matcher.add(pattern1089, 1089)
    pattern1090 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons46, cons226, cons13, cons137, cons246, cons673)
    def replacement1090(p, n2, x, b, c, a, n):
        rubi.append(1090)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n*(n*(S(2)*p + S(3)) + S(1)) + n*(p + S(1))*(-S(4)*a*c + b**S(2))), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1090 = ReplacementRule(pattern1090, replacement1090)

    if load_rule:
        matcher.add(pattern1090, 1090)
    def With1091(c, n2, x, b, a, n):
        q = Rt(a/c, S(2))
        r = Rt(-b/c + S(2)*q, S(2))
        rubi.append(1091)
        return Dist(1/(2*c*q*r), Int((r - x**(n/2))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((r + x**(n/2))/(q + r*x**(n/2) + x**n), x), x)
    pattern1091 = Pattern(Integral(S(1)/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons46, cons226, cons674, cons413)
    rule1091 = ReplacementRule(pattern1091, With1091)

    if load_rule:
        matcher.add(pattern1091, 1091)
    def With1092(c, n2, x, b, a, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1092)
        return Dist(c/q, Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x) - Dist(c/q, Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x)
    pattern1092 = Pattern(Integral(S(1)/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons46, cons226)
    rule1092 = ReplacementRule(pattern1092, With1092)

    if load_rule:
        matcher.add(pattern1092, 1092)
    def With1093(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1093)
        return Dist(S(2)*sqrt(-c), Int(S(1)/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)
    pattern1093 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons293)
    rule1093 = ReplacementRule(pattern1093, With1093)

    if load_rule:
        matcher.add(pattern1093, 1093)
    def With1094(b, a, x, c):
        q = Rt(c/a, S(4))
        rubi.append(1094)
        return Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(2)*q*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    pattern1094 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons676, cons677)
    rule1094 = ReplacementRule(pattern1094, With1094)

    if load_rule:
        matcher.add(pattern1094, 1094)

    pattern1095 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons484, cons177, cons_with_1095)
    def replacement1095(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1095)
        return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt(-S(2)*a - x**S(2)*(b - q))*EllipticF(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*sqrt(-a)*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1095 = ReplacementRule(pattern1095, replacement1095)

    if load_rule:
        matcher.add(pattern1095, 1095)
    def With1096(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1096)
        return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*EllipticF(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*sqrt(a/(S(2)*a + x**S(2)*(b + q)))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    pattern1096 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons484, cons177)
    rule1096 = ReplacementRule(pattern1096, With1096)

    if load_rule:
        matcher.add(pattern1096, 1096)

    pattern1097 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1097)
    def replacement1097(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1097)
        return Simp(sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*(S(2)*a + x**S(2)*(b + q))*EllipticF(ArcTan(x*Rt((b + q)/(S(2)*a), S(2))), S(2)*q/(b + q))/(S(2)*a*sqrt(a + b*x**S(2) + c*x**S(4))*Rt((b + q)/(S(2)*a), S(2))), x)
    rule1097 = ReplacementRule(pattern1097, replacement1097)

    if load_rule:
        matcher.add(pattern1097, 1097)

    pattern1098 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1098)
    def replacement1098(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1098)
        return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/(S(2)*a + x**S(2)*(b - q)))*(S(2)*a + x**S(2)*(b - q))*EllipticF(ArcTan(x*Rt((b - q)/(S(2)*a), S(2))), -S(2)*q/(b - q))/(S(2)*a*sqrt(a + b*x**S(2) + c*x**S(4))*Rt((b - q)/(S(2)*a), S(2))), x)
    rule1098 = ReplacementRule(pattern1098, replacement1098)

    if load_rule:
        matcher.add(pattern1098, 1098)

    pattern1099 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1099)
    def replacement1099(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1099)
        return Simp(sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticF(asin(x*Rt(-(b + q)/(S(2)*a), S(2))), (b - q)/(b + q))/(sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-(b + q)/(S(2)*a), S(2))), x)
    rule1099 = ReplacementRule(pattern1099, replacement1099)

    if load_rule:
        matcher.add(pattern1099, 1099)

    pattern1100 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1100)
    def replacement1100(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1100)
        return Simp(sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticF(asin(x*Rt(-(b - q)/(S(2)*a), S(2))), (b + q)/(b - q))/(sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-(b - q)/(S(2)*a), S(2))), x)
    rule1100 = ReplacementRule(pattern1100, replacement1100)

    if load_rule:
        matcher.add(pattern1100, 1100)
    def With1101(b, a, x, c):
        q = Rt(c/a, S(4))
        rubi.append(1101)
        return Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(2)*q*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    pattern1101 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons678)
    rule1101 = ReplacementRule(pattern1101, With1101)

    if load_rule:
        matcher.add(pattern1101, 1101)
    def With1102(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1102)
        return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(S(1)/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)
    pattern1102 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons679)
    rule1102 = ReplacementRule(pattern1102, With1102)

    if load_rule:
        matcher.add(pattern1102, 1102)
    pattern1103 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons680)
    def replacement1103(p, n2, x, b, c, a, n):
        rubi.append(1103)
        return Dist(a**IntPart(p)*(S(2)*c*x**n/(b - Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(S(2)*c*x**n/(b + Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((S(2)*c*x**n/(b - sqrt(-S(4)*a*c + b**S(2))) + S(1))**p*(S(2)*c*x**n/(b + sqrt(-S(4)*a*c + b**S(2))) + S(1))**p, x), x)
    rule1103 = ReplacementRule(pattern1103, replacement1103)

    if load_rule:
        matcher.add(pattern1103, 1103)
    pattern1104 = Pattern(Integral((a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons585, cons38, cons681)
    def replacement1104(mn, p, x, b, c, a, n):
        rubi.append(1104)
        return Int(x**(-n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x)
    rule1104 = ReplacementRule(pattern1104, replacement1104)

    if load_rule:
        matcher.add(pattern1104, 1104)
    pattern1105 = Pattern(Integral((a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons585, cons147, cons681)
    def replacement1105(mn, p, x, b, c, a, n):
        rubi.append(1105)
        return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(-n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x), x)
    rule1105 = ReplacementRule(pattern1105, replacement1105)

    if load_rule:
        matcher.add(pattern1105, 1105)
    pattern1106 = Pattern(Integral((a_ + u_**n_*WC('b', S(1)) + u_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons68, cons69)
    def replacement1106(p, u, n2, x, b, c, a, n):
        rubi.append(1106)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, u), x)
    rule1106 = ReplacementRule(pattern1106, replacement1106)

    if load_rule:
        matcher.add(pattern1106, 1106)
    pattern1107 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons680, cons53)
    def replacement1107(p, n2, x, b, c, a, m, n):
        rubi.append(1107)
        return Dist(S(1)/n, Subst(Int((a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1107 = ReplacementRule(pattern1107, replacement1107)

    if load_rule:
        matcher.add(pattern1107, 1107)
    pattern1108 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons680, cons128, cons501)
    def replacement1108(d, p, n2, x, b, c, a, m, n):
        rubi.append(1108)
        return Int(ExpandIntegrand((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1108 = ReplacementRule(pattern1108, replacement1108)

    if load_rule:
        matcher.add(pattern1108, 1108)
    pattern1109 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons680, cons63, cons502)
    def replacement1109(p, n2, x, b, c, a, m, n):
        rubi.append(1109)
        return Int(x**(m + S(2)*n*p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1109 = ReplacementRule(pattern1109, replacement1109)

    if load_rule:
        matcher.add(pattern1109, 1109)
    pattern1110 = Pattern(Integral(sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))/x_, x_), cons2, cons3, cons7, cons4, cons680, cons45)
    def replacement1110(n2, x, b, c, a, n):
        rubi.append(1110)
        return Simp(sqrt(a + b*x**n + c*x**(S(2)*n))/n, x) + Simp(b*sqrt(a + b*x**n + c*x**(S(2)*n))*log(x)/(b + S(2)*c*x**n), x)
    rule1110 = ReplacementRule(pattern1110, replacement1110)

    if load_rule:
        matcher.add(pattern1110, 1110)
    pattern1111 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons7, cons4, cons680, cons45, cons13, cons146)
    def replacement1111(p, n2, x, b, c, a, n):
        rubi.append(1111)
        return Dist(a, Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/x, x), x) + Simp((a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p), x) + Simp((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(S(2)*n*(S(2)*p + S(-1))), x)
    rule1111 = ReplacementRule(pattern1111, replacement1111)

    if load_rule:
        matcher.add(pattern1111, 1111)
    pattern1112 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons7, cons4, cons680, cons45, cons13, cons137)
    def replacement1112(p, n2, x, b, c, a, n):
        rubi.append(1112)
        return Dist(S(1)/a, Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/x, x), x) - Simp((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*n*(p + S(1))), x) - Simp((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*n*(S(2)*p + S(1))), x)
    rule1112 = ReplacementRule(pattern1112, replacement1112)

    if load_rule:
        matcher.add(pattern1112, 1112)
    pattern1113 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons7, cons4, cons5, cons680, cons45, cons147)
    def replacement1113(p, n2, x, b, c, a, n):
        rubi.append(1113)
        return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b/S(2) + c*x**n)**(S(2)*p)/x, x), x)
    rule1113 = ReplacementRule(pattern1113, replacement1113)

    if load_rule:
        matcher.add(pattern1113, 1113)
    pattern1114 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680, cons45, cons682)
    def replacement1114(d, p, n2, x, b, c, a, m, n):
        rubi.append(1114)
        return Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*(m + S(1))), x)
    rule1114 = ReplacementRule(pattern1114, replacement1114)

    if load_rule:
        matcher.add(pattern1114, 1114)
    pattern1115 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons680, cons45, cons155)
    def replacement1115(d, n2, x, b, c, a, m, n):
        rubi.append(1115)
        return Dist(sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int((d*x)**m*(b + S(2)*c*x**n), x), x)
    rule1115 = ReplacementRule(pattern1115, replacement1115)

    if load_rule:
        matcher.add(pattern1115, 1115)
    pattern1116 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons680, cons45, cons683)
    def replacement1116(d, n2, x, b, c, a, m, n):
        rubi.append(1116)
        return Simp((d*x)**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(d*(m + n + S(1))), x) + Simp(b*n*(d*x)**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(d*(b + S(2)*c*x**n)*(m + S(1))*(m + n + S(1))), x)
    rule1116 = ReplacementRule(pattern1116, replacement1116)

    if load_rule:
        matcher.add(pattern1116, 1116)
    pattern1117 = Pattern(Integral(x_**WC('m', S(1))/sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons680, cons45, cons155)
    def replacement1117(n2, x, b, c, a, m, n):
        rubi.append(1117)
        return -Dist(b/(S(2)*a), Int(S(1)/(x*sqrt(a + b*x**n + c*x**(S(2)*n))), x), x) - Simp(x**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(a*n), x)
    rule1117 = ReplacementRule(pattern1117, replacement1117)

    if load_rule:
        matcher.add(pattern1117, 1117)
    pattern1118 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680, cons45, cons684, cons225)
    def replacement1118(d, p, n2, x, b, c, a, m, n):
        rubi.append(1118)
        return -Simp((d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*d*n*(S(2)*p + S(1))), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*d*n*(p + S(1))*(S(2)*p + S(1))), x)
    rule1118 = ReplacementRule(pattern1118, replacement1118)

    if load_rule:
        matcher.add(pattern1118, 1118)
    pattern1119 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons680, cons45, cons56, cons243)
    def replacement1119(p, n2, x, b, c, a, m, n):
        rubi.append(1119)
        return -Dist(b/(S(2)*c), Int(x**(n + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*c*n*(p + S(1))), x)
    rule1119 = ReplacementRule(pattern1119, replacement1119)

    if load_rule:
        matcher.add(pattern1119, 1119)
    pattern1120 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons45, cons148, cons244, cons146, cons685, cons246, cons17)
    def replacement1120(d, p, n2, x, b, c, a, m, n):
        rubi.append(1120)
        return -Dist(b*d**(-n)*n**S(2)*p*(S(2)*p + S(-1))/((m + S(1))*(m + S(2)*n*p + S(1))), Int((d*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(1))*(m + S(2)*n*p + S(1))), x)
    rule1120 = ReplacementRule(pattern1120, replacement1120)

    if load_rule:
        matcher.add(pattern1120, 1120)
    pattern1121 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons45, cons148, cons244, cons146, cons686, cons687, cons246, cons17)
    def replacement1121(d, p, n2, x, b, c, a, m, n):
        rubi.append(1121)
        return Dist(S(2)*c*d**(-S(2)*n)*n**S(2)*p*(S(2)*p + S(-1))/((m + S(1))*(m + n + S(1))), Int((d*x)**(m + S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*(m - n*(S(2)*p + S(-1)) + S(1))/(d*(m + S(1))*(m + n + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(1))*(m + n + S(1))), x)
    rule1121 = ReplacementRule(pattern1121, replacement1121)

    if load_rule:
        matcher.add(pattern1121, 1121)
    pattern1122 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons680, cons45, cons148, cons13, cons146, cons510, cons688, cons687, cons689, cons246)
    def replacement1122(d, p, n2, x, b, c, a, m, n):
        rubi.append(1122)
        return Dist(S(2)*a*n**S(2)*p*(S(2)*p + S(-1))/((m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), x)
    rule1122 = ReplacementRule(pattern1122, replacement1122)

    if load_rule:
        matcher.add(pattern1122, 1122)
    pattern1123 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons45, cons148, cons244, cons137, cons690, cons246)
    def replacement1123(d, p, n2, x, b, c, a, m, n):
        rubi.append(1123)
        return -Dist(d**n*(m - n + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(b*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*n*(S(2)*p + S(1))), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(b*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)
    rule1123 = ReplacementRule(pattern1123, replacement1123)

    if load_rule:
        matcher.add(pattern1123, 1123)
    pattern1124 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons45, cons148, cons244, cons137, cons529, cons246)
    def replacement1124(d, p, n2, x, b, c, a, m, n):
        rubi.append(1124)
        return Dist(d**(S(2)*n)*(m - S(2)*n + S(1))*(m - n + S(1))/(S(2)*c*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*c*n*(S(2)*p + S(1))), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m - S(2)*n*p - S(3)*n + S(1))/(S(2)*c*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)
    rule1124 = ReplacementRule(pattern1124, replacement1124)

    if load_rule:
        matcher.add(pattern1124, 1124)
    pattern1125 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons680, cons45, cons148, cons244, cons137, cons246)
    def replacement1125(d, p, n2, x, b, c, a, m, n):
        rubi.append(1125)
        return Dist((m + S(2)*n*(p + S(1)) + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp((d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*d*n*(S(2)*p + S(1))), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(S(2)*a*d*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)
    rule1125 = ReplacementRule(pattern1125, replacement1125)

    if load_rule:
        matcher.add(pattern1125, 1125)
    pattern1126 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons45, cons148, cons31, cons530, cons510, cons691)
    def replacement1126(d, p, n2, x, b, c, a, m, n):
        rubi.append(1126)
        return -Dist(b*d**n*(m - n + S(1))/(S(2)*c*(m + S(2)*n*p + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*c*(m + S(2)*n*p + S(1))), x)
    rule1126 = ReplacementRule(pattern1126, replacement1126)

    if load_rule:
        matcher.add(pattern1126, 1126)
    pattern1127 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons45, cons148, cons31, cons94, cons691)
    def replacement1127(d, p, n2, x, b, c, a, m, n):
        rubi.append(1127)
        return -Dist(S(2)*c*d**(-n)*(m + n*(S(2)*p + S(1)) + S(1))/(b*(m + S(1))), Int((d*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*(m + S(1))), x)
    rule1127 = ReplacementRule(pattern1127, replacement1127)

    if load_rule:
        matcher.add(pattern1127, 1127)
    pattern1128 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons680, cons45, cons196, cons17)
    def replacement1128(p, n2, x, b, c, a, m, n):
        rubi.append(1128)
        return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)
    rule1128 = ReplacementRule(pattern1128, replacement1128)

    if load_rule:
        matcher.add(pattern1128, 1128)
    def With1129(d, p, n2, x, b, c, a, m, n):
        k = Denominator(m)
        rubi.append(1129)
        return -Dist(k/d, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(-k*n) + c*d**(-S(2)*n)*x**(-S(2)*k*n))**p, x), x, (d*x)**(-S(1)/k)), x)
    pattern1129 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons45, cons196, cons367)
    rule1129 = ReplacementRule(pattern1129, With1129)

    if load_rule:
        matcher.add(pattern1129, 1129)
    pattern1130 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons680, cons45, cons196, cons356)
    def replacement1130(d, p, n2, x, b, c, a, m, n):
        rubi.append(1130)
        return -Dist(d**IntPart(m)*(d*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)
    rule1130 = ReplacementRule(pattern1130, replacement1130)

    if load_rule:
        matcher.add(pattern1130, 1130)
    pattern1131 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680, cons45, cons147)
    def replacement1131(d, p, n2, x, b, c, a, m, n):
        rubi.append(1131)
        return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d*x)**m*(b/S(2) + c*x**n)**(S(2)*p), x), x)
    rule1131 = ReplacementRule(pattern1131, replacement1131)

    if load_rule:
        matcher.add(pattern1131, 1131)
    pattern1132 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons680, cons226, cons500)
    def replacement1132(p, n2, x, b, c, a, m, n):
        rubi.append(1132)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1132 = ReplacementRule(pattern1132, replacement1132)

    if load_rule:
        matcher.add(pattern1132, 1132)
    pattern1133 = Pattern(Integral((d_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680, cons226, cons500)
    def replacement1133(d, p, n2, x, b, c, a, m, n):
        rubi.append(1133)
        return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1133 = ReplacementRule(pattern1133, replacement1133)

    if load_rule:
        matcher.add(pattern1133, 1133)

    pattern1134 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons680, cons226, cons148, cons17, cons_with_1134)
    def replacement1134(p, n2, x, b, c, a, m, n):

        k = GCD(m + S(1), n)
        rubi.append(1134)
        return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k) + c*x**(S(2)*n/k))**p, x), x, x**k), x)
    rule1134 = ReplacementRule(pattern1134, replacement1134)

    if load_rule:
        matcher.add(pattern1134, 1134)
    def With1135(d, p, n2, x, b, c, a, m, n):
        k = Denominator(m)
        rubi.append(1135)
        return Dist(k/d, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(k*n) + c*d**(-S(2)*n)*x**(S(2)*k*n))**p, x), x, (d*x)**(S(1)/k)), x)
    pattern1135 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons226, cons148, cons367, cons38)
    rule1135 = ReplacementRule(pattern1135, With1135)

    if load_rule:
        matcher.add(pattern1135, 1135)
    pattern1136 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons244, cons163, cons530, cons692, cons693, cons694)
    def replacement1136(d, p, n2, x, b, c, a, m, n):
        rubi.append(1136)
        return -Dist(d**n*n*p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(a*b*(m - n + S(1)) - x**n*(S(2)*a*c*(m + n*(S(2)*p + S(-1)) + S(1)) - b**S(2)*(m + n*(p + S(-1)) + S(1))), x), x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b*n*p + c*x**n*(m + n*(S(2)*p + S(-1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), x)
    rule1136 = ReplacementRule(pattern1136, replacement1136)

    if load_rule:
        matcher.add(pattern1136, 1136)
    pattern1137 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons244, cons163, cons94, cons694)
    def replacement1137(d, p, n2, x, b, c, a, m, n):
        rubi.append(1137)
        return -Dist(d**(-n)*n*p/(m + S(1)), Int((d*x)**(m + n)*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(1))), x)
    rule1137 = ReplacementRule(pattern1137, replacement1137)

    if load_rule:
        matcher.add(pattern1137, 1137)
    pattern1138 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons680, cons226, cons148, cons13, cons163, cons510, cons694)
    def replacement1138(d, p, n2, x, b, c, a, m, n):
        rubi.append(1138)
        return Dist(n*p/(m + S(2)*n*p + S(1)), Int((d*x)**m*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x)
    rule1138 = ReplacementRule(pattern1138, replacement1138)

    if load_rule:
        matcher.add(pattern1138, 1138)
    pattern1139 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons244, cons137, cons690, cons694)
    def replacement1139(d, p, n2, x, b, c, a, m, n):
        rubi.append(1139)
        return -Dist(d**n/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**(m - n)*(b*(m - n + S(1)) + S(2)*c*x**n*(m + S(2)*n*(p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1139 = ReplacementRule(pattern1139, replacement1139)

    if load_rule:
        matcher.add(pattern1139, 1139)
    pattern1140 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons244, cons137, cons529, cons694)
    def replacement1140(d, p, n2, x, b, c, a, m, n):
        rubi.append(1140)
        return Dist(d**(S(2)*n)/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**(m - S(2)*n)*(S(2)*a*(m - S(2)*n + S(1)) + b*x**n*(m + n*(S(2)*p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1140 = ReplacementRule(pattern1140, replacement1140)

    if load_rule:
        matcher.add(pattern1140, 1140)
    pattern1141 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons680, cons226, cons148, cons13, cons137, cons694)
    def replacement1141(d, p, n2, x, b, c, a, m, n):
        rubi.append(1141)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1)) + b*c*x**n*(m + S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*d*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1141 = ReplacementRule(pattern1141, replacement1141)

    if load_rule:
        matcher.add(pattern1141, 1141)
    pattern1142 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons226, cons148, cons31, cons529, cons510, cons694)
    def replacement1142(d, p, n2, x, b, c, a, m, n):
        rubi.append(1142)
        return -Dist(d**(S(2)*n)/(c*(m + S(2)*n*p + S(1))), Int((d*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*(m - S(2)*n + S(1)) + b*x**n*(m + n*(p + S(-1)) + S(1)), x), x), x) + Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(m + S(2)*n*p + S(1))), x)
    rule1142 = ReplacementRule(pattern1142, replacement1142)

    if load_rule:
        matcher.add(pattern1142, 1142)
    pattern1143 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons226, cons148, cons31, cons94, cons694)
    def replacement1143(d, p, n2, x, b, c, a, m, n):
        rubi.append(1143)
        return -Dist(d**(-n)/(a*(m + S(1))), Int((d*x)**(m + n)*(b*(m + n*(p + S(1)) + S(1)) + c*x**n*(m + S(2)*n*(p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*d*(m + S(1))), x)
    rule1143 = ReplacementRule(pattern1143, replacement1143)

    if load_rule:
        matcher.add(pattern1143, 1143)
    pattern1144 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons31, cons94)
    def replacement1144(d, n2, x, b, c, a, m, n):
        rubi.append(1144)
        return -Dist(d**(-n)/a, Int((d*x)**(m + n)*(b + c*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Simp((d*x)**(m + S(1))/(a*d*(m + S(1))), x)
    rule1144 = ReplacementRule(pattern1144, replacement1144)

    if load_rule:
        matcher.add(pattern1144, 1144)
    pattern1145 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons680, cons226, cons148, cons17, cons695)
    def replacement1145(n2, x, b, c, a, m, n):
        rubi.append(1145)
        return Int(PolynomialDivide(x**m, a + b*x**n + c*x**(S(2)*n), x), x)
    rule1145 = ReplacementRule(pattern1145, replacement1145)

    if load_rule:
        matcher.add(pattern1145, 1145)
    pattern1146 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons31, cons529)
    def replacement1146(d, n2, x, b, c, a, m, n):
        rubi.append(1146)
        return -Dist(d**(S(2)*n)/c, Int((d*x)**(m - S(2)*n)*(a + b*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))/(c*(m - S(2)*n + S(1))), x)
    rule1146 = ReplacementRule(pattern1146, replacement1146)

    if load_rule:
        matcher.add(pattern1146, 1146)
    def With1147(b, a, x, c):
        q = Rt(a/c, S(2))
        rubi.append(1147)
        return -Dist(S(1)/2, Int((q - x**S(2))/(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/2, Int((q + x**S(2))/(a + b*x**S(2) + c*x**S(4)), x), x)
    pattern1147 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons696, cons697)
    rule1147 = ReplacementRule(pattern1147, With1147)

    if load_rule:
        matcher.add(pattern1147, 1147)
    def With1148(n2, x, b, c, a, m, n):
        q = Rt(a/c, S(2))
        r = Rt(-b/c + S(2)*q, S(2))
        rubi.append(1148)
        return -Dist(1/(2*c*r), Int(x**(m - 3*n/2)*(q - r*x**(n/2))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*r), Int(x**(m - 3*n/2)*(q + r*x**(n/2))/(q + r*x**(n/2) + x**n), x), x)
    pattern1148 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons46, cons226, cons698, cons699, cons413)
    rule1148 = ReplacementRule(pattern1148, With1148)

    if load_rule:
        matcher.add(pattern1148, 1148)
    def With1149(n2, x, b, c, a, m, n):
        q = Rt(a/c, S(2))
        r = Rt(-b/c + S(2)*q, S(2))
        rubi.append(1149)
        return Dist(1/(2*c*r), Int(x**(m - n/2)/(q - r*x**(n/2) + x**n), x), x) - Dist(1/(2*c*r), Int(x**(m - n/2)/(q + r*x**(n/2) + x**n), x), x)
    pattern1149 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons46, cons226, cons698, cons700, cons413)
    rule1149 = ReplacementRule(pattern1149, With1149)

    if load_rule:
        matcher.add(pattern1149, 1149)
    def With1150(d, n2, x, b, c, a, m, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1150)
        return -Dist(d**n*(b/q + S(-1))/S(2), Int((d*x)**(m - n)/(b/S(2) + c*x**n - q/S(2)), x), x) + Dist(d**n*(b/q + S(1))/S(2), Int((d*x)**(m - n)/(b/S(2) + c*x**n + q/S(2)), x), x)
    pattern1150 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons680, cons226, cons148, cons31, cons701)
    rule1150 = ReplacementRule(pattern1150, With1150)

    if load_rule:
        matcher.add(pattern1150, 1150)
    def With1151(d, n2, x, b, c, a, m, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1151)
        return Dist(c/q, Int((d*x)**m/(b/S(2) + c*x**n - q/S(2)), x), x) - Dist(c/q, Int((d*x)**m/(b/S(2) + c*x**n + q/S(2)), x), x)
    pattern1151 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons680, cons226, cons148)
    rule1151 = ReplacementRule(pattern1151, With1151)

    if load_rule:
        matcher.add(pattern1151, 1151)
    def With1152(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1152)
        return Dist(S(2)*sqrt(-c), Int(x**S(2)/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)
    pattern1152 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons293)
    rule1152 = ReplacementRule(pattern1152, With1152)

    if load_rule:
        matcher.add(pattern1152, 1152)
    def With1153(b, a, x, c):
        q = Rt(c/a, S(2))
        rubi.append(1153)
        return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    pattern1153 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons676, cons677)
    rule1153 = ReplacementRule(pattern1153, With1153)

    if load_rule:
        matcher.add(pattern1153, 1153)
    def With1154(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1154)
        return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b - q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    pattern1154 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons484, cons177)
    rule1154 = ReplacementRule(pattern1154, With1154)

    if load_rule:
        matcher.add(pattern1154, 1154)

    pattern1155 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1155)
    def replacement1155(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1155)
        return Simp(x*(b + S(2)*c*x**S(2) + q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*(S(2)*a + x**S(2)*(b + q))*EllipticE(ArcTan(x*Rt((b + q)/(S(2)*a), S(2))), S(2)*q/(b + q))*Rt((b + q)/(S(2)*a), S(2))/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1155 = ReplacementRule(pattern1155, replacement1155)

    if load_rule:
        matcher.add(pattern1155, 1155)

    pattern1156 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1156)
    def replacement1156(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1156)
        return Simp(x*(b + S(2)*c*x**S(2) - q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(sqrt((S(2)*a + x**S(2)*(b + q))/(S(2)*a + x**S(2)*(b - q)))*(S(2)*a + x**S(2)*(b - q))*EllipticE(ArcTan(x*Rt((b - q)/(S(2)*a), S(2))), -S(2)*q/(b - q))*Rt((b - q)/(S(2)*a), S(2))/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1156 = ReplacementRule(pattern1156, replacement1156)

    if load_rule:
        matcher.add(pattern1156, 1156)

    pattern1157 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1157)
    def replacement1157(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1157)
        return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) + q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b + q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1157 = ReplacementRule(pattern1157, replacement1157)

    if load_rule:
        matcher.add(pattern1157, 1157)

    pattern1158 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons675, cons_with_1158)
    def replacement1158(b, a, x, c):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1158)
        return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b - q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1158 = ReplacementRule(pattern1158, replacement1158)

    if load_rule:
        matcher.add(pattern1158, 1158)
    def With1159(b, a, x, c):
        q = Rt(c/a, S(2))
        rubi.append(1159)
        return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    pattern1159 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons678)
    rule1159 = ReplacementRule(pattern1159, With1159)

    if load_rule:
        matcher.add(pattern1159, 1159)
    def With1160(b, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1160)
        return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(x**S(2)/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)
    pattern1160 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons226, cons679)
    rule1160 = ReplacementRule(pattern1160, With1160)

    if load_rule:
        matcher.add(pattern1160, 1160)
    pattern1161 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons5, cons680, cons226, cons196, cons17)
    def replacement1161(p, n2, x, b, c, a, m, n):
        rubi.append(1161)
        return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)
    rule1161 = ReplacementRule(pattern1161, replacement1161)

    if load_rule:
        matcher.add(pattern1161, 1161)
    def With1162(d, p, n2, x, b, c, a, m, n):
        k = Denominator(m)
        rubi.append(1162)
        return -Dist(k/d, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(-k*n) + c*d**(-S(2)*n)*x**(-S(2)*k*n))**p, x), x, (d*x)**(-S(1)/k)), x)
    pattern1162 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons5, cons680, cons226, cons196, cons367)
    rule1162 = ReplacementRule(pattern1162, With1162)

    if load_rule:
        matcher.add(pattern1162, 1162)
    pattern1163 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons680, cons226, cons196, cons356)
    def replacement1163(d, p, n2, x, b, c, a, m, n):
        rubi.append(1163)
        return -Dist(d**IntPart(m)*(d*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)
    rule1163 = ReplacementRule(pattern1163, replacement1163)

    if load_rule:
        matcher.add(pattern1163, 1163)
    def With1164(p, n2, x, b, c, a, m, n):
        k = Denominator(n)
        rubi.append(1164)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n) + c*x**(S(2)*k*n))**p, x), x, x**(S(1)/k)), x)
    pattern1164 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons5, cons680, cons226, cons489)
    rule1164 = ReplacementRule(pattern1164, With1164)

    if load_rule:
        matcher.add(pattern1164, 1164)
    pattern1165 = Pattern(Integral((d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons5, cons680, cons226, cons489)
    def replacement1165(d, p, n2, x, b, c, a, m, n):
        rubi.append(1165)
        return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1165 = ReplacementRule(pattern1165, replacement1165)

    if load_rule:
        matcher.add(pattern1165, 1165)
    pattern1166 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons680, cons226, cons541, cons23)
    def replacement1166(p, n2, x, b, c, a, m, n):
        rubi.append(1166)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p, x), x, x**(m + S(1))), x)
    rule1166 = ReplacementRule(pattern1166, replacement1166)

    if load_rule:
        matcher.add(pattern1166, 1166)
    pattern1167 = Pattern(Integral((d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680, cons226, cons541, cons23)
    def replacement1167(d, p, n2, x, b, c, a, m, n):
        rubi.append(1167)
        return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1167 = ReplacementRule(pattern1167, replacement1167)

    if load_rule:
        matcher.add(pattern1167, 1167)
    def With1168(d, n2, x, b, c, a, m, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1168)
        return Dist(S(2)*c/q, Int((d*x)**m/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int((d*x)**m/(b + S(2)*c*x**n + q), x), x)
    pattern1168 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons680, cons226)
    rule1168 = ReplacementRule(pattern1168, With1168)

    if load_rule:
        matcher.add(pattern1168, 1168)
    pattern1169 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons680, cons226, cons702)
    def replacement1169(d, p, n2, x, b, c, a, m, n):
        rubi.append(1169)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1)) + b*c*x**n*(m + S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*d*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1169 = ReplacementRule(pattern1169, replacement1169)

    if load_rule:
        matcher.add(pattern1169, 1169)
    pattern1170 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons680)
    def replacement1170(d, p, n2, x, b, c, a, m, n):
        rubi.append(1170)
        return Dist(a**IntPart(p)*(S(2)*c*x**n/(b - Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(S(2)*c*x**n/(b + Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d*x)**m*(S(2)*c*x**n/(b - sqrt(-S(4)*a*c + b**S(2))) + S(1))**p*(S(2)*c*x**n/(b + sqrt(-S(4)*a*c + b**S(2))) + S(1))**p, x), x)
    rule1170 = ReplacementRule(pattern1170, replacement1170)

    if load_rule:
        matcher.add(pattern1170, 1170)
    pattern1171 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons585, cons38, cons681)
    def replacement1171(mn, p, x, b, c, a, m, n):
        rubi.append(1171)
        return Int(x**(m - n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x)
    rule1171 = ReplacementRule(pattern1171, replacement1171)

    if load_rule:
        matcher.add(pattern1171, 1171)
    pattern1172 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons585, cons147, cons681)
    def replacement1172(mn, p, x, b, c, a, m, n):
        rubi.append(1172)
        return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(m - n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x), x)
    rule1172 = ReplacementRule(pattern1172, replacement1172)

    if load_rule:
        matcher.add(pattern1172, 1172)
    pattern1173 = Pattern(Integral((d_*x_)**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons585)
    def replacement1173(mn, d, p, x, b, c, a, m, n):
        rubi.append(1173)
        return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**(-n) + c*x**n)**p, x), x)
    rule1173 = ReplacementRule(pattern1173, replacement1173)

    if load_rule:
        matcher.add(pattern1173, 1173)
    pattern1174 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons680, cons552, cons17, cons553)
    def replacement1174(v, p, n2, m, x, b, c, a, n):
        rubi.append(1174)
        return Dist(Coefficient(v, x, S(1))**(-m + S(-1)), Subst(Int(SimplifyIntegrand((x - Coefficient(v, x, S(0)))**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x), x, v), x)
    rule1174 = ReplacementRule(pattern1174, replacement1174)

    if load_rule:
        matcher.add(pattern1174, 1174)
    pattern1175 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons680, cons554)
    def replacement1175(v, p, u, m, n2, x, b, c, a, n):
        rubi.append(1175)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, v), x)
    rule1175 = ReplacementRule(pattern1175, replacement1175)

    if load_rule:
        matcher.add(pattern1175, 1175)
    pattern1176 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons220, cons502)
    def replacement1176(d, q, p, e, n2, x, b, c, a, n):
        rubi.append(1176)
        return Int(x**(n*(S(2)*p + q))*(d*x**(-n) + e)**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1176 = ReplacementRule(pattern1176, replacement1176)

    if load_rule:
        matcher.add(pattern1176, 1176)
    pattern1177 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons220, cons502)
    def replacement1177(d, q, p, e, n2, x, c, a, n):
        rubi.append(1177)
        return Int(x**(n*(S(2)*p + q))*(a*x**(-S(2)*n) + c)**p*(d*x**(-n) + e)**q, x)
    rule1177 = ReplacementRule(pattern1177, replacement1177)

    if load_rule:
        matcher.add(pattern1177, 1177)
    pattern1178 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons46, cons196)
    def replacement1178(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1178)
        return -Subst(Int((d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)
    rule1178 = ReplacementRule(pattern1178, replacement1178)

    if load_rule:
        matcher.add(pattern1178, 1178)
    pattern1179 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons50, cons46, cons196)
    def replacement1179(d, q, p, e, n2, x, c, a, n):
        rubi.append(1179)
        return -Subst(Int((a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q/x**S(2), x), x, S(1)/x)
    rule1179 = ReplacementRule(pattern1179, replacement1179)

    if load_rule:
        matcher.add(pattern1179, 1179)
    def With1180(d, q, p, e, n2, x, b, c, a, n):
        g = Denominator(n)
        rubi.append(1180)
        return Dist(g, Subst(Int(x**(g + S(-1))*(d + e*x**(g*n))**q*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p, x), x, x**(S(1)/g)), x)
    pattern1180 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons46, cons489)
    rule1180 = ReplacementRule(pattern1180, With1180)

    if load_rule:
        matcher.add(pattern1180, 1180)
    def With1181(d, q, p, e, n2, x, c, a, n):
        g = Denominator(n)
        rubi.append(1181)
        return Dist(g, Subst(Int(x**(g + S(-1))*(a + c*x**(S(2)*g*n))**p*(d + e*x**(g*n))**q, x), x, x**(S(1)/g)), x)
    pattern1181 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons50, cons46, cons489)
    rule1181 = ReplacementRule(pattern1181, With1181)

    if load_rule:
        matcher.add(pattern1181, 1181)
    pattern1182 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons4, cons5, cons46, cons147, cons664)
    def replacement1182(d, p, c, e, n2, x, b, n):
        rubi.append(1182)
        return Dist(e/c, Int(x**(-n)*(b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) + Simp(x**(-S(2)*n*(p + S(1)))*(b*e - c*d)*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(b*c*n*(p + S(1))), x)
    rule1182 = ReplacementRule(pattern1182, replacement1182)

    if load_rule:
        matcher.add(pattern1182, 1182)
    pattern1183 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons4, cons5, cons46, cons147, cons671, cons703)
    def replacement1183(d, p, c, e, n2, x, b, n):
        rubi.append(1183)
        return Simp(e*x**(-n + S(1))*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(n*(S(2)*p + S(1)) + S(1))), x)
    rule1183 = ReplacementRule(pattern1183, replacement1183)

    if load_rule:
        matcher.add(pattern1183, 1183)
    pattern1184 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons4, cons5, cons46, cons147, cons671, cons704)
    def replacement1184(d, p, c, e, n2, x, b, n):
        rubi.append(1184)
        return -Dist((b*e*(n*p + S(1)) - c*d*(n*(S(2)*p + S(1)) + S(1)))/(c*(n*(S(2)*p + S(1)) + S(1))), Int((b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(e*x**(-n + S(1))*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(n*(S(2)*p + S(1)) + S(1))), x)
    rule1184 = ReplacementRule(pattern1184, replacement1184)

    if load_rule:
        matcher.add(pattern1184, 1184)
    pattern1185 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons147)
    def replacement1185(d, q, p, c, e, n2, x, b, n):
        rubi.append(1185)
        return Dist(x**(-n*FracPart(p))*(b + c*x**n)**(-FracPart(p))*(b*x**n + c*x**(S(2)*n))**FracPart(p), Int(x**(n*p)*(b + c*x**n)**p*(d + e*x**n)**q, x), x)
    rule1185 = ReplacementRule(pattern1185, replacement1185)

    if load_rule:
        matcher.add(pattern1185, 1185)
    pattern1186 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons45, cons147)
    def replacement1186(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1186)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p)*(d + e*x**n)**q, x), x)
    rule1186 = ReplacementRule(pattern1186, replacement1186)

    if load_rule:
        matcher.add(pattern1186, 1186)
    pattern1187 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons46, cons226, cons256, cons38)
    def replacement1187(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1187)
        return Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)
    rule1187 = ReplacementRule(pattern1187, replacement1187)

    if load_rule:
        matcher.add(pattern1187, 1187)
    pattern1188 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons4, cons50, cons46, cons257, cons38)
    def replacement1188(d, q, p, e, n2, x, c, a, n):
        rubi.append(1188)
        return Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)
    rule1188 = ReplacementRule(pattern1188, replacement1188)

    if load_rule:
        matcher.add(pattern1188, 1188)
    pattern1189 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons226, cons256, cons147)
    def replacement1189(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1189)
        return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)
    rule1189 = ReplacementRule(pattern1189, replacement1189)

    if load_rule:
        matcher.add(pattern1189, 1189)
    pattern1190 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons257, cons147)
    def replacement1190(d, q, p, e, n2, x, c, a, n):
        rubi.append(1190)
        return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)
    rule1190 = ReplacementRule(pattern1190, replacement1190)

    if load_rule:
        matcher.add(pattern1190, 1190)
    pattern1191 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons279, cons221)
    def replacement1191(d, q, c, e, n2, x, b, a, n):
        rubi.append(1191)
        return Int(ExpandIntegrand((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1191 = ReplacementRule(pattern1191, replacement1191)

    if load_rule:
        matcher.add(pattern1191, 1191)
    pattern1192 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons221)
    def replacement1192(d, q, e, n2, x, c, a, n):
        rubi.append(1192)
        return Int(ExpandIntegrand((a + c*x**(S(2)*n))*(d + e*x**n)**q, x), x)
    rule1192 = ReplacementRule(pattern1192, replacement1192)

    if load_rule:
        matcher.add(pattern1192, 1192)
    pattern1193 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons279, cons395, cons396)
    def replacement1193(d, q, c, e, n2, x, b, a, n):
        rubi.append(1193)
        return Dist(S(1)/(d*e**S(2)*n*(q + S(1))), Int((d + e*x**n)**(q + S(1))*Simp(a*e**S(2)*(n*(q + S(1)) + S(1)) - b*d*e + c*d**S(2) + c*d*e*n*x**n*(q + S(1)), x), x), x) - Simp(x*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))/(d*e**S(2)*n*(q + S(1))), x)
    rule1193 = ReplacementRule(pattern1193, replacement1193)

    if load_rule:
        matcher.add(pattern1193, 1193)
    pattern1194 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons395, cons396)
    def replacement1194(d, q, e, n2, x, c, a, n):
        rubi.append(1194)
        return Dist(S(1)/(d*e**S(2)*n*(q + S(1))), Int((d + e*x**n)**(q + S(1))*Simp(a*e**S(2)*(n*(q + S(1)) + S(1)) + c*d**S(2) + c*d*e*n*x**n*(q + S(1)), x), x), x) - Simp(x*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))/(d*e**S(2)*n*(q + S(1))), x)
    rule1194 = ReplacementRule(pattern1194, replacement1194)

    if load_rule:
        matcher.add(pattern1194, 1194)
    pattern1195 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons46, cons226, cons279)
    def replacement1195(d, q, c, e, n2, x, b, a, n):
        rubi.append(1195)
        return Dist(S(1)/(e*(n*(q + S(2)) + S(1))), Int((d + e*x**n)**q*(a*e*(n*(q + S(2)) + S(1)) - x**n*(-b*e*(n*(q + S(2)) + S(1)) + c*d*(n + S(1)))), x), x) + Simp(c*x**(n + S(1))*(d + e*x**n)**(q + S(1))/(e*(n*(q + S(2)) + S(1))), x)
    rule1195 = ReplacementRule(pattern1195, replacement1195)

    if load_rule:
        matcher.add(pattern1195, 1195)
    pattern1196 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons50, cons46, cons280)
    def replacement1196(d, q, e, n2, x, c, a, n):
        rubi.append(1196)
        return Dist(S(1)/(e*(n*(q + S(2)) + S(1))), Int((d + e*x**n)**q*(a*e*(n*(q + S(2)) + S(1)) - c*d*x**n*(n + S(1))), x), x) + Simp(c*x**(n + S(1))*(d + e*x**n)**(q + S(1))/(e*(n*(q + S(2)) + S(1))), x)
    rule1196 = ReplacementRule(pattern1196, replacement1196)

    if load_rule:
        matcher.add(pattern1196, 1196)
    def With1197(d, e, n2, x, c, a, n):
        q = Rt(S(2)*d*e, S(2))
        rubi.append(1197)
        return Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*x**n - q*x**(n/S(2))), x), x) + Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*x**n + q*x**(n/S(2))), x), x)
    pattern1197 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons46, cons705, cons674, cons706)
    rule1197 = ReplacementRule(pattern1197, With1197)

    if load_rule:
        matcher.add(pattern1197, 1197)
    def With1198(d, e, n2, x, c, a, n):
        q = Rt(-S(2)*d*e, S(2))
        rubi.append(1198)
        return Dist(d/(S(2)*a), Int((d - q*x**(n/S(2)))/(d - e*x**n - q*x**(n/S(2))), x), x) + Dist(d/(S(2)*a), Int((d + q*x**(n/S(2)))/(d - e*x**n + q*x**(n/S(2))), x), x)
    pattern1198 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons46, cons705, cons674, cons707)
    rule1198 = ReplacementRule(pattern1198, With1198)

    if load_rule:
        matcher.add(pattern1198, 1198)
    def With1199(d, e, x, c, a):
        q = Rt(a*c, S(2))
        rubi.append(1199)
        return Dist((-a*e + d*q)/(S(2)*a*c), Int((-c*x**S(2) + q)/(a + c*x**S(4)), x), x) + Dist((a*e + d*q)/(S(2)*a*c), Int((c*x**S(2) + q)/(a + c*x**S(4)), x), x)
    pattern1199 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons708, cons697)
    rule1199 = ReplacementRule(pattern1199, With1199)

    if load_rule:
        matcher.add(pattern1199, 1199)
    def With1200(d, e, n2, x, c, a, n):
        q = Rt(a/c, S(4))
        rubi.append(1200)
        return Dist(sqrt(S(2))/(S(4)*c*q**S(3)), Int((sqrt(S(2))*d*q - x**(n/S(2))*(d - e*q**S(2)))/(q**S(2) - sqrt(S(2))*q*x**(n/S(2)) + x**n), x), x) + Dist(sqrt(S(2))/(S(4)*c*q**S(3)), Int((sqrt(S(2))*d*q + x**(n/S(2))*(d - e*q**S(2)))/(q**S(2) + sqrt(S(2))*q*x**(n/S(2)) + x**n), x), x)
    pattern1200 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons46, cons280, cons708, cons674, cons697)
    rule1200 = ReplacementRule(pattern1200, With1200)

    if load_rule:
        matcher.add(pattern1200, 1200)
    def With1201(d, e, x, c, a):
        q = Rt(c/a, S(6))
        rubi.append(1201)
        return Dist(S(1)/(S(6)*a*q**S(2)), Int((S(2)*d*q**S(2) - x*(sqrt(S(3))*d*q**S(3) - e))/(q**S(2)*x**S(2) - sqrt(S(3))*q*x + S(1)), x), x) + Dist(S(1)/(S(6)*a*q**S(2)), Int((S(2)*d*q**S(2) + x*(sqrt(S(3))*d*q**S(3) + e))/(q**S(2)*x**S(2) + sqrt(S(3))*q*x + S(1)), x), x) + Dist(S(1)/(S(3)*a*q**S(2)), Int((d*q**S(2) - e*x)/(q**S(2)*x**S(2) + S(1)), x), x)
    pattern1201 = Pattern(Integral((d_ + x_**S(3)*WC('e', S(1)))/(a_ + x_**S(6)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons678)
    rule1201 = ReplacementRule(pattern1201, With1201)

    if load_rule:
        matcher.add(pattern1201, 1201)
    def With1202(d, e, n2, x, c, a, n):
        q = Rt(-a/c, S(2))
        rubi.append(1202)
        return Dist(d/S(2) - e*q/S(2), Int(S(1)/(a - c*q*x**n), x), x) + Dist(d/S(2) + e*q/S(2), Int(S(1)/(a + c*q*x**n), x), x)
    pattern1202 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons709, cons85)
    rule1202 = ReplacementRule(pattern1202, With1202)

    if load_rule:
        matcher.add(pattern1202, 1202)
    pattern1203 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons710)
    def replacement1203(d, e, n2, x, c, a, n):
        rubi.append(1203)
        return Dist(d, Int(S(1)/(a + c*x**(S(2)*n)), x), x) + Dist(e, Int(x**n/(a + c*x**(S(2)*n)), x), x)
    rule1203 = ReplacementRule(pattern1203, replacement1203)

    if load_rule:
        matcher.add(pattern1203, 1203)
    def With1204(d, c, e, n2, x, b, a, n):
        q = Rt(-b/c + S(2)*d/e, S(2))
        rubi.append(1204)
        return Dist(e**S(2)/(S(2)*c), Int(S(1)/(d - e*q*x**(n/S(2)) + e*x**n), x), x) + Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*q*x**(n/S(2)) + e*x**n), x), x)
    pattern1204 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons46, cons226, cons705, cons674, cons711)
    rule1204 = ReplacementRule(pattern1204, With1204)

    if load_rule:
        matcher.add(pattern1204, 1204)
    def With1205(d, c, e, n2, x, b, a, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1205)
        return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x)
    pattern1205 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons705, cons674, cons675)
    rule1205 = ReplacementRule(pattern1205, With1205)

    if load_rule:
        matcher.add(pattern1205, 1205)
    def With1206(d, c, e, n2, x, b, a, n):
        q = Rt(a/c, S(2))
        r = Rt(-b/c + S(2)*q, S(2))
        rubi.append(1206)
        return Dist(1/(2*c*q*r), Int((d*r - x**(n/2)*(d - e*q))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((d*r + x**(n/2)*(d - e*q))/(q + r*x**(n/2) + x**n), x), x)
    pattern1206 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons46, cons226, cons705, cons674, cons712)
    rule1206 = ReplacementRule(pattern1206, With1206)

    if load_rule:
        matcher.add(pattern1206, 1206)
    def With1207(d, c, e, n2, x, b, a, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1207)
        return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x)
    pattern1207 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons279, cons713)
    rule1207 = ReplacementRule(pattern1207, With1207)

    if load_rule:
        matcher.add(pattern1207, 1207)
    def With1208(d, c, e, n2, x, b, a, n):
        q = Rt(a/c, S(2))
        r = Rt(-b/c + S(2)*q, S(2))
        rubi.append(1208)
        return Dist(1/(2*c*q*r), Int((d*r - x**(n/2)*(d - e*q))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((d*r + x**(n/2)*(d - e*q))/(q + r*x**(n/2) + x**n), x), x)
    pattern1208 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons46, cons226, cons279, cons674, cons413)
    rule1208 = ReplacementRule(pattern1208, With1208)

    if load_rule:
        matcher.add(pattern1208, 1208)
    pattern1209 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons279, cons586)
    def replacement1209(d, q, c, e, n2, x, b, a, n):
        rubi.append(1209)
        return Int(ExpandIntegrand((d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1209 = ReplacementRule(pattern1209, replacement1209)

    if load_rule:
        matcher.add(pattern1209, 1209)
    pattern1210 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons586)
    def replacement1210(d, q, e, n2, x, c, a, n):
        rubi.append(1210)
        return Int(ExpandIntegrand((d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)
    rule1210 = ReplacementRule(pattern1210, replacement1210)

    if load_rule:
        matcher.add(pattern1210, 1210)
    pattern1211 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons279, cons386, cons395, cons396)
    def replacement1211(d, q, c, e, n2, x, b, a, n):
        rubi.append(1211)
        return Dist(e**S(2)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x**n)**q, x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x**n)**(q + S(1))*(-b*e + c*d - c*e*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1211 = ReplacementRule(pattern1211, replacement1211)

    if load_rule:
        matcher.add(pattern1211, 1211)
    pattern1212 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons280, cons386, cons395, cons396)
    def replacement1212(d, q, e, n2, x, c, a, n):
        rubi.append(1212)
        return Dist(c/(a*e**S(2) + c*d**S(2)), Int((d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(e**S(2)/(a*e**S(2) + c*d**S(2)), Int((d + e*x**n)**q, x), x)
    rule1212 = ReplacementRule(pattern1212, replacement1212)

    if load_rule:
        matcher.add(pattern1212, 1212)
    def With1213(d, q, c, e, n2, x, b, a, n):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1213)
        return Dist(S(2)*c/r, Int((d + e*x**n)**q/(b + S(2)*c*x**n - r), x), x) - Dist(S(2)*c/r, Int((d + e*x**n)**q/(b + S(2)*c*x**n + r), x), x)
    pattern1213 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons46, cons226, cons279, cons386)
    rule1213 = ReplacementRule(pattern1213, With1213)

    if load_rule:
        matcher.add(pattern1213, 1213)
    def With1214(d, q, e, n2, x, c, a, n):
        r = Rt(-a*c, S(2))
        rubi.append(1214)
        return -Dist(c/(S(2)*r), Int((d + e*x**n)**q/(-c*x**n + r), x), x) - Dist(c/(S(2)*r), Int((d + e*x**n)**q/(c*x**n + r), x), x)
    pattern1214 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons50, cons46, cons280, cons386)
    rule1214 = ReplacementRule(pattern1214, With1214)

    if load_rule:
        matcher.add(pattern1214, 1214)
    pattern1215 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons149, cons163, cons669, cons714, cons246, cons673)
    def replacement1215(d, p, c, e, n2, x, b, a, n):
        rubi.append(1215)
        return Dist(n*p/(c*(S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(-a*b*e + S(2)*a*c*d*(S(2)*n*p + n + S(1)) + x**n*(S(2)*a*c*e*(S(2)*n*p + S(1)) - b**S(2)*e*(n*p + S(1)) + b*c*d*(S(2)*n*p + n + S(1))), x), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p*(b*e*n*p + c*d*(S(2)*n*p + n + S(1)) + c*e*x**n*(S(2)*n*p + S(1)))/(c*(S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), x)
    rule1215 = ReplacementRule(pattern1215, replacement1215)

    if load_rule:
        matcher.add(pattern1215, 1215)
    pattern1216 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons149, cons163, cons669, cons714, cons246, cons673)
    def replacement1216(d, p, e, n2, x, c, a, n):
        rubi.append(1216)
        return Dist(S(2)*a*n*p/((S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), Int((a + c*x**(S(2)*n))**(p + S(-1))*(d*(S(2)*n*p + n + S(1)) + e*x**n*(S(2)*n*p + S(1))), x), x) + Simp(x*(a + c*x**(S(2)*n))**p*(d*(S(2)*n*p + n + S(1)) + e*x**n*(S(2)*n*p + S(1)))/((S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), x)
    rule1216 = ReplacementRule(pattern1216, replacement1216)

    if load_rule:
        matcher.add(pattern1216, 1216)
    pattern1217 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226, cons13, cons137, cons246, cons673)
    def replacement1217(d, p, c, e, n2, x, b, a, n):
        rubi.append(1217)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e - S(2)*a*c*d*(S(2)*n*p + S(2)*n + S(1)) + b**S(2)*d*(n*p + n + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e - S(2)*a*c*d + b**S(2)*d + c*x**n*(-S(2)*a*e + b*d))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1217 = ReplacementRule(pattern1217, replacement1217)

    if load_rule:
        matcher.add(pattern1217, 1217)
    pattern1218 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46, cons13, cons137, cons246, cons673)
    def replacement1218(d, p, e, n2, x, c, a, n):
        rubi.append(1218)
        return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*(d*(S(2)*n*p + S(2)*n + S(1)) + e*x**n*(S(2)*n*p + S(3)*n + S(1))), x), x) - Simp(x*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*n*(p + S(1))), x)
    rule1218 = ReplacementRule(pattern1218, replacement1218)

    if load_rule:
        matcher.add(pattern1218, 1218)
    def With1219(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1219)
        return Dist(S(2)*sqrt(-c), Int((d + e*x**S(2))/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)
    pattern1219 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons293)
    rule1219 = ReplacementRule(pattern1219, With1219)

    if load_rule:
        matcher.add(pattern1219, 1219)
    def With1220(d, e, x, c, a):
        q = Rt(-a*c, S(2))
        rubi.append(1220)
        return Dist(sqrt(-c), Int((d + e*x**S(2))/(sqrt(-c*x**S(2) + q)*sqrt(c*x**S(2) + q)), x), x)
    pattern1220 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons43, cons293)
    rule1220 = ReplacementRule(pattern1220, With1220)

    if load_rule:
        matcher.add(pattern1220, 1220)

    pattern1221 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons676, cons677, cons_with_1221)
    def replacement1221(d, c, e, x, b, a):

        q = Rt(c/a, S(4))
        rubi.append(1221)
        return -Simp(d*x*sqrt(a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(q*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1221 = ReplacementRule(pattern1221, replacement1221)

    if load_rule:
        matcher.add(pattern1221, 1221)

    pattern1222 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons676, cons677, cons_with_1222)
    def replacement1222(d, c, e, x, b, a):

        q = Rt(c/a, S(2))
        rubi.append(1222)
        return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1222 = ReplacementRule(pattern1222, replacement1222)

    if load_rule:
        matcher.add(pattern1222, 1222)

    pattern1223 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons484, cons177, cons_with_1223)
    def replacement1223(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1223)
        return Simp(e*x*(b + S(2)*c*x**S(2) + q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(e*q*sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*EllipticE(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*c*sqrt(a/(S(2)*a + x**S(2)*(b + q)))*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1223 = ReplacementRule(pattern1223, replacement1223)

    if load_rule:
        matcher.add(pattern1223, 1223)

    pattern1224 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons484, cons177, cons_with_1224)
    def replacement1224(d, e, x, c, a):

        q = Rt(-a*c, S(2))
        rubi.append(1224)
        return Simp(e*x*(c*x**S(2) + q)/(c*sqrt(a + c*x**S(4))), x) - Simp(sqrt(S(2))*e*q*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticE(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(c*sqrt(-a)*sqrt(a + c*x**S(4))), x)
    rule1224 = ReplacementRule(pattern1224, replacement1224)

    if load_rule:
        matcher.add(pattern1224, 1224)

    pattern1225 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons484, cons177, cons_with_1225)
    def replacement1225(d, e, x, c, a):

        q = Rt(-a*c, S(2))
        rubi.append(1225)
        return Simp(e*x*(c*x**S(2) + q)/(c*sqrt(a + c*x**S(4))), x) - Simp(sqrt(S(2))*e*q*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticE(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(c*sqrt(a/(a + q*x**S(2)))*sqrt(a + c*x**S(4))), x)
    rule1225 = ReplacementRule(pattern1225, replacement1225)

    if load_rule:
        matcher.add(pattern1225, 1225)

    pattern1226 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons484, cons177, cons_with_1226)
    def replacement1226(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1226)
        return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b - q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1226 = ReplacementRule(pattern1226, replacement1226)

    if load_rule:
        matcher.add(pattern1226, 1226)

    pattern1227 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons484, cons177, cons_with_1227)
    def replacement1227(d, e, x, c, a):

        q = Rt(-a*c, S(2))
        rubi.append(1227)
        return -Dist(e/c, Int((-c*x**S(2) + q)/sqrt(a + c*x**S(4)), x), x) + Dist((c*d + e*q)/c, Int(S(1)/sqrt(a + c*x**S(4)), x), x)
    rule1227 = ReplacementRule(pattern1227, replacement1227)

    if load_rule:
        matcher.add(pattern1227, 1227)

    pattern1228 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons_with_1228)
    def replacement1228(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1228)
        return Dist(d, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(e, Int(x**S(2)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1228 = ReplacementRule(pattern1228, replacement1228)

    if load_rule:
        matcher.add(pattern1228, 1228)
    pattern1229 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons715)
    def replacement1229(d, e, x, c, a):
        rubi.append(1229)
        return Dist(d, Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist(e, Int(x**S(2)/sqrt(a + c*x**S(4)), x), x)
    rule1229 = ReplacementRule(pattern1229, replacement1229)

    if load_rule:
        matcher.add(pattern1229, 1229)

    pattern1230 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons_with_1230)
    def replacement1230(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1230)
        return -Simp(a*e*sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticE(asin(x*Rt(-(b + q)/(S(2)*a), S(2))), (b - q)/(b + q))*Rt(-(b + q)/(S(2)*a), S(2))/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1230 = ReplacementRule(pattern1230, replacement1230)

    if load_rule:
        matcher.add(pattern1230, 1230)

    pattern1231 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons_with_1231)
    def replacement1231(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1231)
        return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) + q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b + q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1231 = ReplacementRule(pattern1231, replacement1231)

    if load_rule:
        matcher.add(pattern1231, 1231)

    pattern1232 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons_with_1232)
    def replacement1232(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1232)
        return -Simp(a*e*sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticE(asin(x*Rt(-(b - q)/(S(2)*a), S(2))), (b + q)/(b - q))*Rt(-(b - q)/(S(2)*a), S(2))/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1232 = ReplacementRule(pattern1232, replacement1232)

    if load_rule:
        matcher.add(pattern1232, 1232)

    pattern1233 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons_with_1233)
    def replacement1233(d, c, e, x, b, a):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1233)
        return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b - q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1233 = ReplacementRule(pattern1233, replacement1233)

    if load_rule:
        matcher.add(pattern1233, 1233)

    pattern1234 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons678, cons_with_1234)
    def replacement1234(d, c, e, x, b, a):

        q = Rt(c/a, S(4))
        rubi.append(1234)
        return -Simp(d*x*sqrt(a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(q*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1234 = ReplacementRule(pattern1234, replacement1234)

    if load_rule:
        matcher.add(pattern1234, 1234)

    pattern1235 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons678, cons_with_1235)
    def replacement1235(d, e, x, c, a):

        q = Rt(c/a, S(4))
        rubi.append(1235)
        return -Simp(d*x*sqrt(a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), S(1)/2)/(q*sqrt(a + c*x**S(4))), x)
    rule1235 = ReplacementRule(pattern1235, replacement1235)

    if load_rule:
        matcher.add(pattern1235, 1235)

    pattern1236 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons678, cons_with_1236)
    def replacement1236(d, c, e, x, b, a):

        q = Rt(c/a, S(2))
        rubi.append(1236)
        return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)
    rule1236 = ReplacementRule(pattern1236, replacement1236)

    if load_rule:
        matcher.add(pattern1236, 1236)

    pattern1237 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons678, cons_with_1237)
    def replacement1237(d, e, x, c, a):

        q = Rt(c/a, S(2))
        rubi.append(1237)
        return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + c*x**S(4)), x), x)
    rule1237 = ReplacementRule(pattern1237, replacement1237)

    if load_rule:
        matcher.add(pattern1237, 1237)
    pattern1238 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons679, cons257, cons43)
    def replacement1238(d, e, x, c, a):
        rubi.append(1238)
        return Dist(d/sqrt(a), Int(sqrt(S(1) + e*x**S(2)/d)/sqrt(S(1) - e*x**S(2)/d), x), x)
    rule1238 = ReplacementRule(pattern1238, replacement1238)

    if load_rule:
        matcher.add(pattern1238, 1238)
    pattern1239 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons679, cons257, cons448)
    def replacement1239(d, e, x, c, a):
        rubi.append(1239)
        return Dist(sqrt(S(1) + c*x**S(4)/a)/sqrt(a + c*x**S(4)), Int((d + e*x**S(2))/sqrt(S(1) + c*x**S(4)/a), x), x)
    rule1239 = ReplacementRule(pattern1239, replacement1239)

    if load_rule:
        matcher.add(pattern1239, 1239)
    def With1240(d, e, x, c, a):
        q = Rt(-c/a, S(2))
        rubi.append(1240)
        return Dist(e/q, Int((q*x**S(2) + S(1))/sqrt(a + c*x**S(4)), x), x) + Dist((d*q - e)/q, Int(S(1)/sqrt(a + c*x**S(4)), x), x)
    pattern1240 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons679, cons280)
    rule1240 = ReplacementRule(pattern1240, With1240)

    if load_rule:
        matcher.add(pattern1240, 1240)
    def With1241(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1241)
        return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int((d + e*x**S(2))/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)
    pattern1241 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons679)
    rule1241 = ReplacementRule(pattern1241, With1241)

    if load_rule:
        matcher.add(pattern1241, 1241)
    pattern1242 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons46, cons226)
    def replacement1242(d, p, c, e, n2, x, b, a, n):
        rubi.append(1242)
        return Int(ExpandIntegrand((d + e*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1242 = ReplacementRule(pattern1242, replacement1242)

    if load_rule:
        matcher.add(pattern1242, 1242)
    pattern1243 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons4, cons46)
    def replacement1243(d, p, e, n2, x, c, a, n):
        rubi.append(1243)
        return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p*(d + e*x**n), x), x)
    rule1243 = ReplacementRule(pattern1243, replacement1243)

    if load_rule:
        matcher.add(pattern1243, 1243)
    pattern1244 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons46, cons226, cons128, cons716, cons148, cons400)
    def replacement1244(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1244)
        return Int((d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(S(2)*n*p - n + S(1))/(e*(S(2)*n*p + n*q + S(1))) - c**p*x**(S(2)*n*p) + (a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(c**p*x**(S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(S(2)*n*p + n*q + S(1))), x)
    rule1244 = ReplacementRule(pattern1244, replacement1244)

    if load_rule:
        matcher.add(pattern1244, 1244)
    pattern1245 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons50, cons46, cons128, cons716, cons148, cons400)
    def replacement1245(d, q, p, e, n2, x, c, a, n):
        rubi.append(1245)
        return Int((d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(S(2)*n*p - n + S(1))/(e*(S(2)*n*p + n*q + S(1))) - c**p*x**(S(2)*n*p) + (a + c*x**(S(2)*n))**p, x), x) + Simp(c**p*x**(S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(S(2)*n*p + n*q + S(1))), x)
    rule1245 = ReplacementRule(pattern1245, replacement1245)

    if load_rule:
        matcher.add(pattern1245, 1245)
    pattern1246 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279)
    def replacement1246(d, c, e, x, b, a):
        rubi.append(1246)
        return -Dist(e**(S(-2)), Int((-b*e + c*d - c*e*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((a*e**S(2) - b*d*e + c*d**S(2))/e**S(2), Int(S(1)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
    rule1246 = ReplacementRule(pattern1246, replacement1246)

    if load_rule:
        matcher.add(pattern1246, 1246)
    pattern1247 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement1247(d, e, x, c, a):
        rubi.append(1247)
        return -Dist(c/e**S(2), Int((d - e*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist((a*e**S(2) + c*d**S(2))/e**S(2), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)
    rule1247 = ReplacementRule(pattern1247, replacement1247)

    if load_rule:
        matcher.add(pattern1247, 1247)
    def With1248(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1248)
        return -Dist(e**(S(-4)), Int(Simp(-c**S(2)*e**S(3)*x**S(6) + c*e**S(2)*x**S(4)*(-S(2)*b*e + c*d) - S(2)*c*(a*e**S(2) - b*d*e + c*d**S(2))**S(2)/(S(2)*c*d - e*(b + q)) - e*x**S(2)*(b**S(2)*e**S(2) + c**S(2)*d**S(2) - S(2)*c*e*(-a*e + b*d)) + (-b*e + c*d)*(S(2)*a*e**S(2) - b*d*e + c*d**S(2)), x)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((a*e**S(2) - b*d*e + c*d**S(2))**S(2)/(e**S(3)*(S(2)*c*d - e*(b + q))), Int((b + S(2)*c*x**S(2) + q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
    pattern1248 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons675)
    rule1248 = ReplacementRule(pattern1248, With1248)

    if load_rule:
        matcher.add(pattern1248, 1248)
    def With1249(d, e, x, c, a):
        q = Rt(-a*c, S(2))
        rubi.append(1249)
        return -Dist(c/e**S(4), Int(Simp(c*d*e**S(2)*x**S(4) - c*e**S(3)*x**S(6) + d*(S(2)*a*e**S(2) + c*d**S(2)) - e*x**S(2)*(S(2)*a*e**S(2) + c*d**S(2)) - (a*e**S(2) + c*d**S(2))**S(2)/(c*d - e*q), x)/sqrt(a + c*x**S(4)), x), x) - Dist((a*e**S(2) + c*d**S(2))**S(2)/(e**S(3)*(c*d - e*q)), Int((c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)
    pattern1249 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**(S(3)/2)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons715)
    rule1249 = ReplacementRule(pattern1249, With1249)

    if load_rule:
        matcher.add(pattern1249, 1249)
    pattern1250 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons717)
    def replacement1250(d, p, c, e, x, b, a):
        rubi.append(1250)
        return Dist(a, Int((a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(b, Int(x**S(2)*(a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(c, Int(x**S(4)*(a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x)
    rule1250 = ReplacementRule(pattern1250, replacement1250)

    if load_rule:
        matcher.add(pattern1250, 1250)
    pattern1251 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons717)
    def replacement1251(d, p, e, x, c, a):
        rubi.append(1251)
        return Dist(a, Int((a + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(c, Int(x**S(4)*(a + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x)
    rule1251 = ReplacementRule(pattern1251, replacement1251)

    if load_rule:
        matcher.add(pattern1251, 1251)
    def With1252(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1252)
        return Dist(S(2)*sqrt(-c), Int(S(1)/((d + e*x**S(2))*sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)
    pattern1252 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons293)
    rule1252 = ReplacementRule(pattern1252, With1252)

    if load_rule:
        matcher.add(pattern1252, 1252)
    def With1253(d, e, x, c, a):
        q = Rt(-a*c, S(2))
        rubi.append(1253)
        return Dist(sqrt(-c), Int(S(1)/((d + e*x**S(2))*sqrt(-c*x**S(2) + q)*sqrt(c*x**S(2) + q)), x), x)
    pattern1253 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons43, cons293)
    rule1253 = ReplacementRule(pattern1253, With1253)

    if load_rule:
        matcher.add(pattern1253, 1253)
    def With1254(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1254)
        return Dist(S(2)*c/(S(2)*c*d - e*(b - q)), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist(e/(S(2)*c*d - e*(b - q)), Int((b + S(2)*c*x**S(2) - q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
    pattern1254 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons675, cons718)
    rule1254 = ReplacementRule(pattern1254, With1254)

    if load_rule:
        matcher.add(pattern1254, 1254)
    def With1255(d, e, x, c, a):
        q = Rt(-a*c, S(2))
        rubi.append(1255)
        return Dist(c/(c*d + e*q), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist(e/(c*d + e*q), Int((-c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)
    pattern1255 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons715, cons718)
    rule1255 = ReplacementRule(pattern1255, With1255)

    if load_rule:
        matcher.add(pattern1255, 1255)

    pattern1256 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons678, cons_with_1256)
    def replacement1256(d, c, e, x, b, a):

        q = Rt(c/a, S(4))
        rubi.append(1256)
        return -Dist(q**S(2)/(-d*q**S(2) + e), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Simp(ArcTan(x*sqrt((a*e**S(2) - b*d*e + c*d**S(2))/(d*e))/sqrt(a + b*x**S(2) + c*x**S(4)))/(S(2)*d*sqrt((a*e**S(2) - b*d*e + c*d**S(2))/(d*e))), x) + Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(d*q**S(2) + e)*(q**S(2)*x**S(2) + S(1))*EllipticPi(-(-d*q**S(2) + e)**S(2)/(S(4)*d*e*q**S(2)), S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(4)*d*q*(-d*q**S(2) + e)*sqrt(a + b*x**S(2) + c*x**S(4))), x)
    rule1256 = ReplacementRule(pattern1256, replacement1256)

    if load_rule:
        matcher.add(pattern1256, 1256)

    pattern1257 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons280, cons678, cons_with_1257)
    def replacement1257(d, e, x, c, a):

        q = Rt(c/a, S(4))
        rubi.append(1257)
        return -Dist(q**S(2)/(-d*q**S(2) + e), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Simp(ArcTan(x*sqrt((a*e**S(2) + c*d**S(2))/(d*e))/sqrt(a + c*x**S(4)))/(S(2)*d*sqrt((a*e**S(2) + c*d**S(2))/(d*e))), x) + Simp(sqrt((a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(d*q**S(2) + e)*(q**S(2)*x**S(2) + S(1))*EllipticPi(-(-d*q**S(2) + e)**S(2)/(S(4)*d*e*q**S(2)), S(2)*ArcTan(q*x), S(1)/2)/(S(4)*d*q*sqrt(a + c*x**S(4))*(-d*q**S(2) + e)), x)
    rule1257 = ReplacementRule(pattern1257, replacement1257)

    if load_rule:
        matcher.add(pattern1257, 1257)
    def With1258(d, e, x, c, a):
        q = Rt(-c/a, S(4))
        rubi.append(1258)
        return Simp(EllipticPi(-e/(d*q**S(2)), asin(q*x), S(-1))/(sqrt(a)*d*q), x)
    pattern1258 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons679, cons43)
    rule1258 = ReplacementRule(pattern1258, With1258)

    if load_rule:
        matcher.add(pattern1258, 1258)
    pattern1259 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons679, cons448)
    def replacement1259(d, e, x, c, a):
        rubi.append(1259)
        return Dist(sqrt(S(1) + c*x**S(4)/a)/sqrt(a + c*x**S(4)), Int(S(1)/(sqrt(S(1) + c*x**S(4)/a)*(d + e*x**S(2))), x), x)
    rule1259 = ReplacementRule(pattern1259, replacement1259)

    if load_rule:
        matcher.add(pattern1259, 1259)
    def With1260(d, c, e, x, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1260)
        return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(S(1)/((d + e*x**S(2))*sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)
    pattern1260 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons679)
    rule1260 = ReplacementRule(pattern1260, With1260)

    if load_rule:
        matcher.add(pattern1260, 1260)
    pattern1261 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279, cons719)
    def replacement1261(d, p, c, e, x, b, a):
        rubi.append(1261)
        return -Dist(S(1)/(S(2)*a*(p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((a + b*x**S(2) + c*x**S(4))**(p + S(1))*Simp(-a*b*c*d*e*(S(8)*p + S(11)) + S(2)*a*c*(S(4)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(4)*p + S(5))) + b**S(3)*d*e*(S(2)*p + S(3)) - b**S(2)*(S(2)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(2)*p + S(3))) - c*e*x**S(4)*(S(4)*p + S(7))*(S(2)*a*c*e - b**S(2)*e + b*c*d) - x**S(2)*(S(4)*a*c**S(2)*d*e - b**S(3)*e**S(2)*(S(2)*p + S(3)) - S(2)*b**S(2)*c*d*e*(p + S(2)) + b*c*(a*e**S(2)*(S(8)*p + S(11)) + c*d**S(2)*(S(4)*p + S(7)))), x)/(d + e*x**S(2)), x), x) - Simp(x*(a + b*x**S(2) + c*x**S(4))**(p + S(1))*(S(3)*a*b*c*e - S(2)*a*c**S(2)*d - b**S(3)*e + b**S(2)*c*d + c*x**S(2)*(S(2)*a*c*e - b**S(2)*e + b*c*d))/(S(2)*a*(p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)
    rule1261 = ReplacementRule(pattern1261, replacement1261)

    if load_rule:
        matcher.add(pattern1261, 1261)
    pattern1262 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons280, cons719)
    def replacement1262(d, p, e, x, c, a):
        rubi.append(1262)
        return -Dist(-S(1)/(S(8)*a**S(2)*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(4))**(p + S(1))*Simp(-S(4)*a*c**S(2)*d*e*x**S(2) - S(2)*a*c**S(2)*e**S(2)*x**S(4)*(S(4)*p + S(7)) + S(2)*a*c*(S(4)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(4)*p + S(5))), x)/(d + e*x**S(2)), x), x) - Simp(-x*(a + c*x**S(4))**(p + S(1))*(-S(2)*a*c**S(2)*d + S(2)*a*c**S(2)*e*x**S(2))/(S(8)*a**S(2)*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)
    rule1262 = ReplacementRule(pattern1262, replacement1262)

    if load_rule:
        matcher.add(pattern1262, 1262)
    pattern1263 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons226, cons279)
    def replacement1263(d, c, e, x, b, a):
        rubi.append(1263)
        return -Dist(c/(S(2)*d*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((a*e**S(2) - S(2)*b*d*e + S(3)*c*d**S(2))/(S(2)*d*(a*e**S(2) - b*d*e + c*d**S(2))), Int(S(1)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) + Simp(e**S(2)*x*sqrt(a + b*x**S(2) + c*x**S(4))/(S(2)*d*(d + e*x**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)
    rule1263 = ReplacementRule(pattern1263, replacement1263)

    if load_rule:
        matcher.add(pattern1263, 1263)
    pattern1264 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**S(2)), x_), cons2, cons7, cons27, cons48, cons280)
    def replacement1264(d, e, x, c, a):
        rubi.append(1264)
        return -Dist(c/(S(2)*d*(a*e**S(2) + c*d**S(2))), Int((d + e*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist((a*e**S(2) + S(3)*c*d**S(2))/(S(2)*d*(a*e**S(2) + c*d**S(2))), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x) + Simp(e**S(2)*x*sqrt(a + c*x**S(4))/(S(2)*d*(d + e*x**S(2))*(a*e**S(2) + c*d**S(2))), x)
    rule1264 = ReplacementRule(pattern1264, replacement1264)

    if load_rule:
        matcher.add(pattern1264, 1264)
    def With1265(d, q, p, c, e, x, b, a):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1265)
        return Dist(a**IntPart(p)*(S(2)*c*x**S(2)/(b - r) + S(1))**(-FracPart(p))*(S(2)*c*x**S(2)/(b + r) + S(1))**(-FracPart(p))*(a + b*x**S(2) + c*x**S(4))**FracPart(p), Int((d + e*x**S(2))**q*(S(2)*c*x**S(2)/(b - r) + S(1))**p*(S(2)*c*x**S(2)/(b + r) + S(1))**p, x), x)
    pattern1265 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**q_*(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons50, cons226, cons279, cons347, cons564)
    rule1265 = ReplacementRule(pattern1265, With1265)

    if load_rule:
        matcher.add(pattern1265, 1265)
    def With1266(d, q, p, e, x, c, a):
        r = Rt(-a*c, S(2))
        rubi.append(1266)
        return Dist(a**IntPart(p)*(a + c*x**S(4))**FracPart(p)*(-c*x**S(2)/r + S(1))**(-FracPart(p))*(c*x**S(2)/r + S(1))**(-FracPart(p)), Int((d + e*x**S(2))**q*(-c*x**S(2)/r + S(1))**p*(c*x**S(2)/r + S(1))**p, x), x)
    pattern1266 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_*(d_ + x_**S(2)*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons50, cons280, cons347, cons564)
    rule1266 = ReplacementRule(pattern1266, With1266)

    if load_rule:
        matcher.add(pattern1266, 1266)
    pattern1267 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons720, cons43, cons268)
    def replacement1267(d, c, e, x, b, a):
        rubi.append(1267)
        return Simp(EllipticF(S(2)*asin(x*Rt(-e/d, S(2))), b*d/(S(4)*a*e))/(S(2)*sqrt(a)*sqrt(d)*Rt(-e/d, S(2))), x)
    rule1267 = ReplacementRule(pattern1267, replacement1267)

    if load_rule:
        matcher.add(pattern1267, 1267)
    pattern1268 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons720, cons721)
    def replacement1268(d, c, e, x, b, a):
        rubi.append(1268)
        return Dist(sqrt((a + b*x**S(2) + c*x**S(4))/a)*sqrt((d + e*x**S(2))/d)/(sqrt(d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), Int(S(1)/(sqrt(S(1) + e*x**S(2)/d)*sqrt(S(1) + b*x**S(2)/a + c*x**S(4)/a)), x), x)
    rule1268 = ReplacementRule(pattern1268, replacement1268)

    if load_rule:
        matcher.add(pattern1268, 1268)
    pattern1269 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons720, cons43, cons268)
    def replacement1269(d, c, e, x, b, a):
        rubi.append(1269)
        return Simp(sqrt(a)*EllipticE(S(2)*asin(x*Rt(-e/d, S(2))), b*d/(S(4)*a*e))/(S(2)*sqrt(d)*Rt(-e/d, S(2))), x)
    rule1269 = ReplacementRule(pattern1269, replacement1269)

    if load_rule:
        matcher.add(pattern1269, 1269)
    pattern1270 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons720, cons721)
    def replacement1270(d, c, e, x, b, a):
        rubi.append(1270)
        return Dist(sqrt((d + e*x**S(2))/d)*sqrt(a + b*x**S(2) + c*x**S(4))/(sqrt((a + b*x**S(2) + c*x**S(4))/a)*sqrt(d + e*x**S(2))), Int(sqrt(S(1) + b*x**S(2)/a + c*x**S(4)/a)/sqrt(S(1) + e*x**S(2)/d), x), x)
    rule1270 = ReplacementRule(pattern1270, replacement1270)

    if load_rule:
        matcher.add(pattern1270, 1270)
    pattern1271 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons226, cons279, cons722)
    def replacement1271(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1271)
        return Int(ExpandIntegrand((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1271 = ReplacementRule(pattern1271, replacement1271)

    if load_rule:
        matcher.add(pattern1271, 1271)
    pattern1272 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons280, cons722)
    def replacement1272(d, q, p, e, n2, x, c, a, n):
        rubi.append(1272)
        return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1272 = ReplacementRule(pattern1272, replacement1272)

    if load_rule:
        matcher.add(pattern1272, 1272)
    pattern1273 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons5, cons46, cons280, cons564, cons723)
    def replacement1273(d, q, p, e, n2, x, c, a, n):
        rubi.append(1273)
        return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p, (d/(d**S(2) - e**S(2)*x**(S(2)*n)) - e*x**n/(d**S(2) - e**S(2)*x**(S(2)*n)))**(-q), x), x)
    rule1273 = ReplacementRule(pattern1273, replacement1273)

    if load_rule:
        matcher.add(pattern1273, 1273)
    pattern1274 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons724)
    def replacement1274(d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1274)
        return Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1274 = ReplacementRule(pattern1274, replacement1274)

    if load_rule:
        matcher.add(pattern1274, 1274)
    pattern1275 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons724)
    def replacement1275(d, q, p, e, n2, x, c, a, n):
        rubi.append(1275)
        return Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x)
    rule1275 = ReplacementRule(pattern1275, replacement1275)

    if load_rule:
        matcher.add(pattern1275, 1275)
    pattern1276 = Pattern(Integral((d_ + u_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + u_**n2_*WC('c', S(1)) + u_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons68, cons69)
    def replacement1276(d, q, p, c, e, u, n2, x, b, a, n):
        rubi.append(1276)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x, u), x)
    rule1276 = ReplacementRule(pattern1276, replacement1276)

    if load_rule:
        matcher.add(pattern1276, 1276)
    pattern1277 = Pattern(Integral((a_ + u_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + u_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons4, cons5, cons50, cons46, cons68, cons69)
    def replacement1277(d, q, p, e, u, n2, x, c, a, n):
        rubi.append(1277)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x, u), x)
    rule1277 = ReplacementRule(pattern1277, replacement1277)

    if load_rule:
        matcher.add(pattern1277, 1277)
    pattern1278 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons680, cons585, cons586)
    def replacement1278(mn, d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1278)
        return Int(x**(-n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1278 = ReplacementRule(pattern1278, replacement1278)

    if load_rule:
        matcher.add(pattern1278, 1278)
    pattern1279 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons726, cons5, cons725, cons586)
    def replacement1279(mn, d, q, p, e, n2, x, c, a):
        rubi.append(1279)
        return Int(x**(mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x)
    rule1279 = ReplacementRule(pattern1279, replacement1279)

    if load_rule:
        matcher.add(pattern1279, 1279)
    pattern1280 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons680, cons585, cons386, cons38)
    def replacement1280(mn, d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1280)
        return Int(x**(S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1280 = ReplacementRule(pattern1280, replacement1280)

    if load_rule:
        matcher.add(pattern1280, 1280)
    pattern1281 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons726, cons50, cons725, cons386, cons38)
    def replacement1281(mn, d, q, p, e, n2, x, c, a):
        rubi.append(1281)
        return Int(x**(-S(2)*mn*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x)
    rule1281 = ReplacementRule(pattern1281, replacement1281)

    if load_rule:
        matcher.add(pattern1281, 1281)
    pattern1282 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons680, cons585, cons386, cons147, cons681)
    def replacement1282(mn, d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1282)
        return Dist(x**(n*FracPart(q))*(d + e*x**(-n))**FracPart(q)*(d*x**n + e)**(-FracPart(q)), Int(x**(-n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1282 = ReplacementRule(pattern1282, replacement1282)

    if load_rule:
        matcher.add(pattern1282, 1282)
    pattern1283 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons726, cons5, cons50, cons725, cons386, cons147, cons727)
    def replacement1283(mn, d, q, p, e, n2, x, c, a):
        rubi.append(1283)
        return Dist(x**(-mn*FracPart(q))*(d + e*x**mn)**FracPart(q)*(d*x**(-mn) + e)**(-FracPart(q)), Int(x**(mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x), x)
    rule1283 = ReplacementRule(pattern1283, replacement1283)

    if load_rule:
        matcher.add(pattern1283, 1283)
    pattern1284 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons680, cons585, cons386, cons147, cons502)
    def replacement1284(mn, d, q, p, c, e, n2, x, b, a, n):
        rubi.append(1284)
        return Dist(x**(-S(2)*n*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**(-FracPart(p)), Int(x**(S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x), x)
    rule1284 = ReplacementRule(pattern1284, replacement1284)

    if load_rule:
        matcher.add(pattern1284, 1284)
    pattern1285 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons726, cons50, cons725, cons386, cons147, cons728)
    def replacement1285(mn, d, q, p, e, n2, x, c, a):
        rubi.append(1285)
        return Dist(x**(-n2*FracPart(p))*(a + c*x**n2)**FracPart(p)*(a*x**(S(2)*mn) + c)**(-FracPart(p)), Int(x**(n2*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x), x)
    rule1285 = ReplacementRule(pattern1285, replacement1285)

    if load_rule:
        matcher.add(pattern1285, 1285)
    pattern1286 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons50, cons585, cons38)
    def replacement1286(mn, d, q, p, e, x, b, c, a, n):
        rubi.append(1286)
        return Int(x**(-n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x)
    rule1286 = ReplacementRule(pattern1286, replacement1286)

    if load_rule:
        matcher.add(pattern1286, 1286)
    pattern1287 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons50, cons585, cons147)
    def replacement1287(mn, d, q, p, e, x, b, c, a, n):
        rubi.append(1287)
        return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(-n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x), x)
    rule1287 = ReplacementRule(pattern1287, replacement1287)

    if load_rule:
        matcher.add(pattern1287, 1287)
    pattern1288 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons50, cons52, cons46, cons45, cons147)
    def replacement1288(d, q, p, c, f, e, n2, x, r, n, b, a, g):
        rubi.append(1288)
        return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p)*(d + e*x**n)**q*(f + g*x**n)**r, x), x)
    rule1288 = ReplacementRule(pattern1288, replacement1288)

    if load_rule:
        matcher.add(pattern1288, 1288)
    pattern1289 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons50, cons52, cons46, cons226, cons256, cons38)
    def replacement1289(d, q, p, c, f, e, n2, x, r, n, b, a, g):
        rubi.append(1289)
        return Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x)
    rule1289 = ReplacementRule(pattern1289, replacement1289)

    if load_rule:
        matcher.add(pattern1289, 1289)
    pattern1290 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons50, cons52, cons46, cons257, cons38)
    def replacement1290(d, q, p, f, e, n2, x, r, n, c, a, g):
        rubi.append(1290)
        return Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x)
    rule1290 = ReplacementRule(pattern1290, replacement1290)

    if load_rule:
        matcher.add(pattern1290, 1290)
    pattern1291 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons50, cons52, cons46, cons226, cons256, cons147)
    def replacement1291(d, q, p, c, f, e, n2, x, r, n, b, a, g):
        rubi.append(1291)
        return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x), x)
    rule1291 = ReplacementRule(pattern1291, replacement1291)

    if load_rule:
        matcher.add(pattern1291, 1291)
    pattern1292 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons50, cons52, cons46, cons257, cons147)
    def replacement1292(d, q, p, f, e, n2, x, r, n, c, a, g):
        rubi.append(1292)
        return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x), x)
    rule1292 = ReplacementRule(pattern1292, replacement1292)

    if load_rule:
        matcher.add(pattern1292, 1292)

    pattern1293 = Pattern(Integral((x_**S(2)*WC('g', S(1)) + WC('f', S(0)))/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons675, cons279, cons718, cons_with_1293)
    def replacement1293(d, c, f, e, x, b, a, g):

        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1293)
        return Dist((S(2)*c*f - g*(b - q))/(S(2)*c*d - e*(b - q)), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((-d*g + e*f)/(S(2)*c*d - e*(b - q)), Int((b + S(2)*c*x**S(2) - q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)
    rule1293 = ReplacementRule(pattern1293, replacement1293)

    if load_rule:
        matcher.add(pattern1293, 1293)

    pattern1294 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons7, cons27, cons48, cons125, cons208, cons715, cons280, cons718, cons_with_1294)
    def replacement1294(d, f, e, x, c, a, g):

        q = Rt(-a*c, S(2))
        rubi.append(1294)
        return Dist((c*f + g*q)/(c*d + e*q), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist((-d*g + e*f)/(c*d + e*q), Int((-c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)
    rule1294 = ReplacementRule(pattern1294, replacement1294)

    if load_rule:
        matcher.add(pattern1294, 1294)
    pattern1295 = Pattern(Integral((d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons731, cons652, cons732, cons654, cons4, cons5, cons50, cons46, cons593, cons729, cons730)
    def replacement1295(q, p, c, e2, non2, d1, d2, n2, x, e1, b, a, n):
        rubi.append(1295)
        return Int((d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1295 = ReplacementRule(pattern1295, replacement1295)

    if load_rule:
        matcher.add(pattern1295, 1295)
    pattern1296 = Pattern(Integral((d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons731, cons652, cons732, cons654, cons4, cons5, cons50, cons46, cons593, cons729)
    def replacement1296(q, p, c, e2, non2, d1, d2, n2, x, e1, b, a, n):
        rubi.append(1296)
        return Dist((d1 + e1*x**(n/S(2)))**FracPart(q)*(d2 + e2*x**(n/S(2)))**FracPart(q)*(d1*d2 + e1*e2*x**n)**(-FracPart(q)), Int((d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1296 = ReplacementRule(pattern1296, replacement1296)

    if load_rule:
        matcher.add(pattern1296, 1296)
    pattern1297 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons34, cons35, cons21, cons4, cons5, cons50, cons46, cons53)
    def replacement1297(B, d, q, p, c, e, A, n2, x, b, a, m, n):
        rubi.append(1297)
        return Dist(A, Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Dist(B, Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1297 = ReplacementRule(pattern1297, replacement1297)

    if load_rule:
        matcher.add(pattern1297, 1297)
    pattern1298 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons34, cons35, cons21, cons4, cons5, cons50, cons46, cons53)
    def replacement1298(B, d, q, p, e, A, n2, x, c, a, m, n):
        rubi.append(1298)
        return Dist(A, Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x) + Dist(B, Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1298 = ReplacementRule(pattern1298, replacement1298)

    if load_rule:
        matcher.add(pattern1298, 1298)
    pattern1299 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons733, cons500)
    def replacement1299(q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1299)
        return Dist(e**(S(1) - (m + S(1))/n)*f**m/n, Subst(Int((e*x)**(q + S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1299 = ReplacementRule(pattern1299, replacement1299)

    if load_rule:
        matcher.add(pattern1299, 1299)
    pattern1300 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons733, cons500)
    def replacement1300(q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1300)
        return Dist(e**(S(1) - (m + S(1))/n)*f**m/n, Subst(Int((e*x)**(q + S(-1) + (m + S(1))/n)*(a + c*x**S(2))**p, x), x, x**n), x)
    rule1300 = ReplacementRule(pattern1300, replacement1300)

    if load_rule:
        matcher.add(pattern1300, 1300)
    pattern1301 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons733, cons501)
    def replacement1301(q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1301)
        return Dist(e**IntPart(q)*f**m*x**(-n*FracPart(q))*(e*x**n)**FracPart(q), Int(x**(m + n*q)*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1301 = ReplacementRule(pattern1301, replacement1301)

    if load_rule:
        matcher.add(pattern1301, 1301)
    pattern1302 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons733, cons501)
    def replacement1302(q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1302)
        return Dist(e**IntPart(q)*f**m*x**(-n*FracPart(q))*(e*x**n)**FracPart(q), Int(x**(m + n*q)*(a + c*x**(S(2)*n))**p, x), x)
    rule1302 = ReplacementRule(pattern1302, replacement1302)

    if load_rule:
        matcher.add(pattern1302, 1302)
    pattern1303 = Pattern(Integral((f_*x_)**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons18)
    def replacement1303(q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1303)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1303 = ReplacementRule(pattern1303, replacement1303)

    if load_rule:
        matcher.add(pattern1303, 1303)
    pattern1304 = Pattern(Integral((f_*x_)**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons18)
    def replacement1304(q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1304)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(e*x**n)**q*(a + c*x**(S(2)*n))**p, x), x)
    rule1304 = ReplacementRule(pattern1304, replacement1304)

    if load_rule:
        matcher.add(pattern1304, 1304)
    pattern1305 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons53)
    def replacement1305(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1305)
        return Dist(S(1)/n, Subst(Int((d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1305 = ReplacementRule(pattern1305, replacement1305)

    if load_rule:
        matcher.add(pattern1305, 1305)
    pattern1306 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons53)
    def replacement1306(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1306)
        return Dist(S(1)/n, Subst(Int((a + c*x**S(2))**p*(d + e*x)**q, x), x, x**n), x)
    rule1306 = ReplacementRule(pattern1306, replacement1306)

    if load_rule:
        matcher.add(pattern1306, 1306)
    pattern1307 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons46, cons220, cons502)
    def replacement1307(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1307)
        return Int(x**(m + n*(S(2)*p + q))*(d*x**(-n) + e)**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1307 = ReplacementRule(pattern1307, replacement1307)

    if load_rule:
        matcher.add(pattern1307, 1307)
    pattern1308 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons4, cons46, cons220, cons502)
    def replacement1308(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1308)
        return Int(x**(m + n*(S(2)*p + q))*(a*x**(-S(2)*n) + c)**p*(d*x**(-n) + e)**q, x)
    rule1308 = ReplacementRule(pattern1308, replacement1308)

    if load_rule:
        matcher.add(pattern1308, 1308)
    pattern1309 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons46, cons45, cons147, cons734)
    def replacement1309(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1309)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1309 = ReplacementRule(pattern1309, replacement1309)

    if load_rule:
        matcher.add(pattern1309, 1309)
    pattern1310 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons45, cons147)
    def replacement1310(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1310)
        return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((f*x)**m*(b/S(2) + c*x**n)**(S(2)*p)*(d + e*x**n)**q, x), x)
    rule1310 = ReplacementRule(pattern1310, replacement1310)

    if load_rule:
        matcher.add(pattern1310, 1310)
    pattern1311 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons500)
    def replacement1311(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1311)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)
    rule1311 = ReplacementRule(pattern1311, replacement1311)

    if load_rule:
        matcher.add(pattern1311, 1311)
    pattern1312 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons500)
    def replacement1312(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1312)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + c*x**S(2))**p*(d + e*x)**q, x), x, x**n), x)
    rule1312 = ReplacementRule(pattern1312, replacement1312)

    if load_rule:
        matcher.add(pattern1312, 1312)
    pattern1313 = Pattern(Integral((f_*x_)**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons500)
    def replacement1313(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1313)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1313 = ReplacementRule(pattern1313, replacement1313)

    if load_rule:
        matcher.add(pattern1313, 1313)
    pattern1314 = Pattern(Integral((f_*x_)**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons500)
    def replacement1314(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1314)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1314 = ReplacementRule(pattern1314, replacement1314)

    if load_rule:
        matcher.add(pattern1314, 1314)
    pattern1315 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons50, cons46, cons226, cons256, cons38)
    def replacement1315(d, q, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1315)
        return Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)
    rule1315 = ReplacementRule(pattern1315, replacement1315)

    if load_rule:
        matcher.add(pattern1315, 1315)
    pattern1316 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons50, cons21, cons4, cons50, cons46, cons257, cons38)
    def replacement1316(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1316)
        return Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)
    rule1316 = ReplacementRule(pattern1316, replacement1316)

    if load_rule:
        matcher.add(pattern1316, 1316)
    pattern1317 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons226, cons256, cons147)
    def replacement1317(d, q, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1317)
        return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)
    rule1317 = ReplacementRule(pattern1317, replacement1317)

    if load_rule:
        matcher.add(pattern1317, 1317)
    pattern1318 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons257, cons147)
    def replacement1318(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1318)
        return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)
    rule1318 = ReplacementRule(pattern1318, replacement1318)

    if load_rule:
        matcher.add(pattern1318, 1318)
    pattern1319 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons46, cons226, cons464, cons735, cons396, cons168)
    def replacement1319(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1319)
        return Dist(e**(-S(2)*p - (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**Mod(m, n)*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p + (m - Mod(m, n))/n)*n*x**(m - Mod(m, n))*(q + S(1))*(a + b*x**n + c*x**(S(2)*n))**p - (-d)**(S(-1) + (m - Mod(m, n))/n)*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1)))*(a*e**S(2) - b*d*e + c*d**S(2))**p)/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))**p/(n*(q + S(1))), x)
    rule1319 = ReplacementRule(pattern1319, replacement1319)

    if load_rule:
        matcher.add(pattern1319, 1319)
    pattern1320 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons46, cons464, cons735, cons396, cons168)
    def replacement1320(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1320)
        return Dist(e**(-S(2)*p - (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**Mod(m, n)*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p + (m - Mod(m, n))/n)*n*x**(m - Mod(m, n))*(a + c*x**(S(2)*n))**p*(q + S(1)) - (-d)**(S(-1) + (m - Mod(m, n))/n)*(a*e**S(2) + c*d**S(2))**p*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1))))/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))**p/(n*(q + S(1))), x)
    rule1320 = ReplacementRule(pattern1320, replacement1320)

    if load_rule:
        matcher.add(pattern1320, 1320)
    pattern1321 = Pattern(Integral(x_**m_*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons46, cons226, cons464, cons735, cons396, cons267)
    def replacement1321(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1321)
        return Dist(e**(-S(2)*p)*(-d)**(S(-1) + (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**m*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p)*n*(-d)**(S(1) - (m - Mod(m, n))/n)*(q + S(1))*(a + b*x**n + c*x**(S(2)*n))**p - e**(-(m - Mod(m, n))/n)*x**(-m + Mod(m, n))*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1)))*(a*e**S(2) - b*d*e + c*d**S(2))**p)/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))**p/(n*(q + S(1))), x)
    rule1321 = ReplacementRule(pattern1321, replacement1321)

    if load_rule:
        matcher.add(pattern1321, 1321)
    pattern1322 = Pattern(Integral(x_**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons46, cons464, cons735, cons396, cons267)
    def replacement1322(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1322)
        return Dist(e**(-S(2)*p)*(-d)**(S(-1) + (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**m*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p)*n*(-d)**(S(1) - (m - Mod(m, n))/n)*(a + c*x**(S(2)*n))**p*(q + S(1)) - e**(-(m - Mod(m, n))/n)*x**(-m + Mod(m, n))*(a*e**S(2) + c*d**S(2))**p*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1))))/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))**p/(n*(q + S(1))), x)
    rule1322 = ReplacementRule(pattern1322, replacement1322)

    if load_rule:
        matcher.add(pattern1322, 1322)
    pattern1323 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons226, cons464, cons736, cons386, cons737)
    def replacement1323(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1323)
        return Dist(S(1)/(e*(m + S(2)*n*p + n*q + S(1))), Int((f*x)**m*(d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(m + S(2)*n*p - n + S(1)) + e*(-c**p*x**(S(2)*n*p) + (a + b*x**n + c*x**(S(2)*n))**p)*(m + S(2)*n*p + n*q + S(1)), x), x), x) + Simp(c**p*f**(-S(2)*n*p + n + S(-1))*(f*x)**(m + S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(m + S(2)*n*p + n*q + S(1))), x)
    rule1323 = ReplacementRule(pattern1323, replacement1323)

    if load_rule:
        matcher.add(pattern1323, 1323)
    pattern1324 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons464, cons736, cons386, cons737)
    def replacement1324(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1324)
        return Dist(S(1)/(e*(m + S(2)*n*p + n*q + S(1))), Int((f*x)**m*(d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(m + S(2)*n*p - n + S(1)) + e*(-c**p*x**(S(2)*n*p) + (a + c*x**(S(2)*n))**p)*(m + S(2)*n*p + n*q + S(1)), x), x), x) + Simp(c**p*f**(-S(2)*n*p + n + S(-1))*(f*x)**(m + S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(m + S(2)*n*p + n*q + S(1))), x)
    rule1324 = ReplacementRule(pattern1324, replacement1324)

    if load_rule:
        matcher.add(pattern1324, 1324)
    pattern1325 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons464)
    def replacement1325(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1325)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1325 = ReplacementRule(pattern1325, replacement1325)

    if load_rule:
        matcher.add(pattern1325, 1325)
    pattern1326 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons46, cons464)
    def replacement1326(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1326)
        return Int(ExpandIntegrand((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1326 = ReplacementRule(pattern1326, replacement1326)

    if load_rule:
        matcher.add(pattern1326, 1326)

    pattern1327 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons46, cons226, cons148, cons17, cons_with_1327)
    def replacement1327(d, q, p, e, n2, x, b, c, a, m, n):

        k = GCD(m + S(1), n)
        rubi.append(1327)
        return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(d + e*x**(n/k))**q*(a + b*x**(n/k) + c*x**(S(2)*n/k))**p, x), x, x**k), x)
    rule1327 = ReplacementRule(pattern1327, replacement1327)

    if load_rule:
        matcher.add(pattern1327, 1327)

    pattern1328 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons50, cons46, cons148, cons17, cons_with_1328)
    def replacement1328(d, q, p, e, n2, x, c, a, m, n):

        k = GCD(m + S(1), n)
        rubi.append(1328)
        return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + c*x**(S(2)*n/k))**p*(d + e*x**(n/k))**q, x), x, x**k), x)
    rule1328 = ReplacementRule(pattern1328, replacement1328)

    if load_rule:
        matcher.add(pattern1328, 1328)
    def With1329(d, q, p, f, e, n2, x, b, c, a, m, n):
        k = Denominator(m)
        rubi.append(1329)
        return Dist(k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))*(d + e*f**(-n)*x**(k*n))**q*(a + b*f**(-n)*x**(k*n) + c*f**(-S(2)*n)*x**(S(2)*k*n))**p, x), x, (f*x)**(S(1)/k)), x)
    pattern1329 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons46, cons226, cons148, cons367, cons38)
    rule1329 = ReplacementRule(pattern1329, With1329)

    if load_rule:
        matcher.add(pattern1329, 1329)
    def With1330(d, q, p, f, e, n2, x, c, a, m, n):
        k = Denominator(m)
        rubi.append(1330)
        return Dist(k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + c*x**(S(2)*k*n)/f)**p*(d + e*x**(k*n)/f)**q, x), x, (f*x)**(S(1)/k)), x)
    pattern1330 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons5, cons50, cons46, cons148, cons367, cons38)
    rule1330 = ReplacementRule(pattern1330, With1330)

    if load_rule:
        matcher.add(pattern1330, 1330)
    pattern1331 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons163, cons94, cons738, cons694)
    def replacement1331(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1331)
        return Dist(f**(-n)*n*p/((m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(S(2)*a*e*(m + S(1)) - b*d*(m + n*(S(2)*p + S(1)) + S(1)) + x**n*(b*e*(m + S(1)) - S(2)*c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp((f*x)**(m + S(1))*(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p/(f*(m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1331 = ReplacementRule(pattern1331, replacement1331)

    if load_rule:
        matcher.add(pattern1331, 1331)
    pattern1332 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons163, cons94, cons738, cons694)
    def replacement1332(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1332)
        return Dist(S(2)*f**(-n)*n*p/((m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**(p + S(-1))*(a*e*(m + S(1)) - c*d*x**n*(m + n*(S(2)*p + S(1)) + S(1))), x), x) + Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**p*(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(1)))/(f*(m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1332 = ReplacementRule(pattern1332, replacement1332)

    if load_rule:
        matcher.add(pattern1332, 1332)
    pattern1333 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148, cons13, cons163, cons510, cons739, cons694)
    def replacement1333(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1333)
        return Dist(n*p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(-a*b*e*(m + S(1)) + S(2)*a*c*d*(m + n*(S(2)*p + S(1)) + S(1)) + x**n*(S(2)*a*c*e*(m + S(2)*n*p + S(1)) - b**S(2)*e*(m + n*p + S(1)) + b*c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*(b*e*n*p + c*d*(m + n*(S(2)*p + S(1)) + S(1)) + c*e*x**n*(m + S(2)*n*p + S(1)))/(c*f*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1333 = ReplacementRule(pattern1333, replacement1333)

    if load_rule:
        matcher.add(pattern1333, 1333)
    pattern1334 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148, cons13, cons163, cons510, cons739, cons694)
    def replacement1334(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1334)
        return Dist(S(2)*a*n*p/((m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*Simp(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(2)*n*p + S(1)), x), x), x) + Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**p*(c*d*(m + n*(S(2)*p + S(1)) + S(1)) + c*e*x**n*(m + S(2)*n*p + S(1)))/(c*f*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1334 = ReplacementRule(pattern1334, replacement1334)

    if load_rule:
        matcher.add(pattern1334, 1334)
    pattern1335 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons137, cons530, cons694)
    def replacement1335(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1335)
        return Dist(f**n/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(x**n*(b*e - S(2)*c*d)*(m + S(2)*n*p + S(2)*n + S(1)) + (-S(2)*a*e + b*d)*(-m + n + S(-1)), x), x), x) + Simp(f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*e + b*d - x**n*(b*e - S(2)*c*d))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1335 = ReplacementRule(pattern1335, replacement1335)

    if load_rule:
        matcher.add(pattern1335, 1335)
    pattern1336 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons137, cons530, cons694)
    def replacement1336(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1336)
        return Dist(f**n/(S(2)*a*c*n*(p + S(1))), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**(p + S(1))*(a*e*(-m + n + S(-1)) + c*d*x**n*(m + S(2)*n*p + S(2)*n + S(1))), x), x) + Simp(f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(a*e - c*d*x**n)/(S(2)*a*c*n*(p + S(1))), x)
    rule1336 = ReplacementRule(pattern1336, replacement1336)

    if load_rule:
        matcher.add(pattern1336, 1336)
    pattern1337 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148, cons13, cons137, cons694)
    def replacement1337(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1337)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e*(m + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(m + n*(S(2)*p + S(3)) + S(1)) + d*(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1))), x), x), x) - Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e + c*x**n*(-S(2)*a*e + b*d) + d*(-S(2)*a*c + b**S(2)))/(a*f*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1337 = ReplacementRule(pattern1337, replacement1337)

    if load_rule:
        matcher.add(pattern1337, 1337)
    pattern1338 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148, cons13, cons137, cons694)
    def replacement1338(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1338)
        return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(1))*Simp(d*(m + S(2)*n*(p + S(1)) + S(1)) + e*x**n*(m + n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*f*n*(p + S(1))), x)
    rule1338 = ReplacementRule(pattern1338, replacement1338)

    if load_rule:
        matcher.add(pattern1338, 1338)
    pattern1339 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons46, cons226, cons148, cons31, cons530, cons739, cons694)
    def replacement1339(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1339)
        return -Dist(f**n/(c*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*e*(m - n + S(1)) + x**n*(b*e*(m + n*p + S(1)) - c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp(e*f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1339 = ReplacementRule(pattern1339, replacement1339)

    if load_rule:
        matcher.add(pattern1339, 1339)
    pattern1340 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons5, cons46, cons148, cons31, cons530, cons739, cons694)
    def replacement1340(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1340)
        return -Dist(f**n/(c*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**p*(a*e*(m - n + S(1)) - c*d*x**n*(m + n*(S(2)*p + S(1)) + S(1))), x), x) + Simp(e*f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + c*x**(S(2)*n))**(p + S(1))/(c*(m + n*(S(2)*p + S(1)) + S(1))), x)
    rule1340 = ReplacementRule(pattern1340, replacement1340)

    if load_rule:
        matcher.add(pattern1340, 1340)
    pattern1341 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons46, cons226, cons148, cons31, cons94, cons694)
    def replacement1341(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1341)
        return Dist(f**(-n)/(a*(m + S(1))), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*e*(m + S(1)) - b*d*(m + n*(p + S(1)) + S(1)) - c*d*x**n*(m + S(2)*n*(p + S(1)) + S(1)), x), x), x) + Simp(d*(f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*f*(m + S(1))), x)
    rule1341 = ReplacementRule(pattern1341, replacement1341)

    if load_rule:
        matcher.add(pattern1341, 1341)
    pattern1342 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons5, cons46, cons148, cons31, cons94, cons694)
    def replacement1342(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1342)
        return Dist(f**(-n)/(a*(m + S(1))), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**p*(a*e*(m + S(1)) - c*d*x**n*(m + S(2)*n*(p + S(1)) + S(1))), x), x) + Simp(d*(f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))/(a*f*(m + S(1))), x)
    rule1342 = ReplacementRule(pattern1342, replacement1342)

    if load_rule:
        matcher.add(pattern1342, 1342)

    pattern1343 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons696, cons740, cons81, cons697, cons_with_1343)
    def replacement1343(d, c, f, e, n2, x, b, a, m, n):

        q = Rt(a*c, S(2))
        r = Rt(-b*c + S(2)*c*q, S(2))
        rubi.append(1343)
        return Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r - x**(n/2)*(c*d - e*q), x)/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r + x**(n/2)*(c*d - e*q), x)/(c*x**n + q + r*x**(n/2)), x), x)
    rule1343 = ReplacementRule(pattern1343, replacement1343)

    if load_rule:
        matcher.add(pattern1343, 1343)

    pattern1344 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons434, cons740, cons81, cons_with_1344)
    def replacement1344(d, f, e, n2, x, c, a, m, n):

        q = Rt(a*c, S(2))
        r = Rt(S(2)*c*q, S(2))
        rubi.append(1344)
        return Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r - x**(n/2)*(c*d - e*q), x)/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r + x**(n/2)*(c*d - e*q), x)/(c*x**n + q + r*x**(n/2)), x), x)
    rule1344 = ReplacementRule(pattern1344, replacement1344)

    if load_rule:
        matcher.add(pattern1344, 1344)
    def With1345(d, c, f, e, x, b, a, m):
        r = Rt(c*(-b*e + S(2)*c*d)/e, S(2))
        rubi.append(1345)
        return Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) - r*x), x), x) + Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) + r*x), x), x)
    pattern1345 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons226, cons705, cons741, cons742)
    rule1345 = ReplacementRule(pattern1345, With1345)

    if load_rule:
        matcher.add(pattern1345, 1345)
    def With1346(d, f, e, x, c, a, m):
        r = Rt(S(2)*c**S(2)*d/e, S(2))
        rubi.append(1346)
        return Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) - r*x), x), x) + Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) + r*x), x), x)
    pattern1346 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons705, cons741)
    rule1346 = ReplacementRule(pattern1346, With1346)

    if load_rule:
        matcher.add(pattern1346, 1346)

    pattern1347 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons696, cons743, cons744, cons697, cons_with_1347)
    def replacement1347(d, c, f, e, n2, x, b, a, m, n):

        q = Rt(a*c, S(2))
        r = Rt(-b*c + S(2)*c*q, S(2))
        rubi.append(1347)
        return Dist(c/(2*q*r), Int((f*x)**m*(d*r - x**(n/2)*(c*d - e*q))/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*(d*r + x**(n/2)*(c*d - e*q))/(c*x**n + q + r*x**(n/2)), x), x)
    rule1347 = ReplacementRule(pattern1347, replacement1347)

    if load_rule:
        matcher.add(pattern1347, 1347)

    pattern1348 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons743, cons744, cons434, cons_with_1348)
    def replacement1348(d, f, e, n2, x, c, a, m, n):

        q = Rt(a*c, S(2))
        r = Rt(S(2)*c*q, S(2))
        rubi.append(1348)
        return Dist(c/(2*q*r), Int((f*x)**m*(d*r - x**(n/2)*(c*d - e*q))/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*(d*r + x**(n/2)*(c*d - e*q))/(c*x**n + q + r*x**(n/2)), x), x)
    rule1348 = ReplacementRule(pattern1348, replacement1348)

    if load_rule:
        matcher.add(pattern1348, 1348)
    def With1349(d, c, f, e, n2, x, b, a, m, n):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1349)
        return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int((f*x)**m/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int((f*x)**m/(b/S(2) + c*x**n - q/S(2)), x), x)
    pattern1349 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148)
    rule1349 = ReplacementRule(pattern1349, With1349)

    if load_rule:
        matcher.add(pattern1349, 1349)
    def With1350(d, f, e, n2, x, c, a, m, n):
        q = Rt(-a*c, S(2))
        rubi.append(1350)
        return Dist(-c*d/(S(2)*q) + e/S(2), Int((f*x)**m/(c*x**n + q), x), x) - Dist(c*d/(S(2)*q) + e/S(2), Int((f*x)**m/(-c*x**n + q), x), x)
    pattern1350 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148)
    rule1350 = ReplacementRule(pattern1350, With1350)

    if load_rule:
        matcher.add(pattern1350, 1350)
    pattern1351 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148, cons586, cons17)
    def replacement1351(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1351)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1351 = ReplacementRule(pattern1351, replacement1351)

    if load_rule:
        matcher.add(pattern1351, 1351)
    pattern1352 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148, cons586, cons17)
    def replacement1352(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1352)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)
    rule1352 = ReplacementRule(pattern1352, replacement1352)

    if load_rule:
        matcher.add(pattern1352, 1352)
    pattern1353 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148, cons586, cons18)
    def replacement1353(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1353)
        return Int(ExpandIntegrand((f*x)**m, (d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1353 = ReplacementRule(pattern1353, replacement1353)

    if load_rule:
        matcher.add(pattern1353, 1353)
    pattern1354 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148, cons586, cons18)
    def replacement1354(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1354)
        return Int(ExpandIntegrand((f*x)**m, (d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)
    rule1354 = ReplacementRule(pattern1354, replacement1354)

    if load_rule:
        matcher.add(pattern1354, 1354)
    pattern1355 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons386, cons611, cons403, cons529)
    def replacement1355(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1355)
        return Dist(f**(S(2)*n)/c**S(2), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(-1))*(-b*e + c*d + c*e*x**n), x), x) - Dist(f**(S(2)*n)/c**S(2), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(-1))*Simp(a*(-b*e + c*d) + x**n*(a*c*e - b**S(2)*e + b*c*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1355 = ReplacementRule(pattern1355, replacement1355)

    if load_rule:
        matcher.add(pattern1355, 1355)
    pattern1356 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons50, cons46, cons148, cons386, cons31, cons529)
    def replacement1356(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1356)
        return Dist(f**(S(2)*n)/c, Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x) - Dist(a*f**(S(2)*n)/c, Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)
    rule1356 = ReplacementRule(pattern1356, replacement1356)

    if load_rule:
        matcher.add(pattern1356, 1356)
    pattern1357 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons386, cons611, cons403, cons690)
    def replacement1357(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1357)
        return -Dist(f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - x**n*(-b*e + c*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Dist(e*f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1)), x), x)
    rule1357 = ReplacementRule(pattern1357, replacement1357)

    if load_rule:
        matcher.add(pattern1357, 1357)
    pattern1358 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons386, cons611, cons403, cons690)
    def replacement1358(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1358)
        return -Dist(f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x) + Dist(e*f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1)), x), x)
    rule1358 = ReplacementRule(pattern1358, replacement1358)

    if load_rule:
        matcher.add(pattern1358, 1358)
    pattern1359 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons386, cons611, cons403, cons267)
    def replacement1359(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1359)
        return Dist(d/a, Int((f*x)**m*(d + e*x**n)**(q + S(-1)), x), x) - Dist(f**(-n)/a, Int((f*x)**(m + n)*(d + e*x**n)**(q + S(-1))*Simp(-a*e + b*d + c*d*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1359 = ReplacementRule(pattern1359, replacement1359)

    if load_rule:
        matcher.add(pattern1359, 1359)
    pattern1360 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons386, cons611, cons403, cons267)
    def replacement1360(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1360)
        return Dist(d/a, Int((f*x)**m*(d + e*x**n)**(q + S(-1)), x), x) + Dist(f**(-n)/a, Int((f*x)**(m + n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x)
    rule1360 = ReplacementRule(pattern1360, replacement1360)

    if load_rule:
        matcher.add(pattern1360, 1360)
    pattern1361 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons386, cons611, cons396, cons529)
    def replacement1361(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1361)
        return -Dist(f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(1))*Simp(a*d + x**n*(-a*e + b*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x)
    rule1361 = ReplacementRule(pattern1361, replacement1361)

    if load_rule:
        matcher.add(pattern1361, 1361)
    pattern1362 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons386, cons611, cons396, cons529)
    def replacement1362(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1362)
        return -Dist(a*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x)
    rule1362 = ReplacementRule(pattern1362, replacement1362)

    if load_rule:
        matcher.add(pattern1362, 1362)
    pattern1363 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons386, cons611, cons396, cons690)
    def replacement1363(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1363)
        return Dist(f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**(q + S(1))*Simp(a*e + c*d*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x) - Dist(d*e*f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**q, x), x)
    rule1363 = ReplacementRule(pattern1363, replacement1363)

    if load_rule:
        matcher.add(pattern1363, 1363)
    pattern1364 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons386, cons611, cons396, cons690)
    def replacement1364(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1364)
        return Dist(f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**(q + S(1))*Simp(a*e + c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x) - Dist(d*e*f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**q, x), x)
    rule1364 = ReplacementRule(pattern1364, replacement1364)

    if load_rule:
        matcher.add(pattern1364, 1364)
    pattern1365 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons46, cons226, cons148, cons386, cons395, cons396)
    def replacement1365(d, q, c, f, e, n2, x, b, a, m, n):
        rubi.append(1365)
        return Dist(e**S(2)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**q, x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**(q + S(1))*Simp(-b*e + c*d - c*e*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1365 = ReplacementRule(pattern1365, replacement1365)

    if load_rule:
        matcher.add(pattern1365, 1365)
    pattern1366 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons46, cons148, cons386, cons395, cons396)
    def replacement1366(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1366)
        return Dist(c/(a*e**S(2) + c*d**S(2)), Int((f*x)**m*(d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(e**S(2)/(a*e**S(2) + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**q, x), x)
    rule1366 = ReplacementRule(pattern1366, replacement1366)

    if load_rule:
        matcher.add(pattern1366, 1366)
    pattern1367 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons50, cons4, cons46, cons226, cons148, cons386, cons17)
    def replacement1367(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1367)
        return Int(ExpandIntegrand((d + e*x**n)**q, (f*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1367 = ReplacementRule(pattern1367, replacement1367)

    if load_rule:
        matcher.add(pattern1367, 1367)
    pattern1368 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons50, cons4, cons46, cons148, cons386, cons17)
    def replacement1368(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1368)
        return Int(ExpandIntegrand((d + e*x**n)**q, (f*x)**m/(a + c*x**(S(2)*n)), x), x)
    rule1368 = ReplacementRule(pattern1368, replacement1368)

    if load_rule:
        matcher.add(pattern1368, 1368)
    pattern1369 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons50, cons4, cons46, cons226, cons148, cons386, cons18)
    def replacement1369(d, q, f, e, n2, x, b, c, a, m, n):
        rubi.append(1369)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q, S(1)/(a + b*x**n + c*x**(S(2)*n)), x), x)
    rule1369 = ReplacementRule(pattern1369, replacement1369)

    if load_rule:
        matcher.add(pattern1369, 1369)
    pattern1370 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons50, cons4, cons46, cons148, cons386, cons18)
    def replacement1370(d, q, f, e, n2, x, c, a, m, n):
        rubi.append(1370)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q, S(1)/(a + c*x**(S(2)*n)), x), x)
    rule1370 = ReplacementRule(pattern1370, replacement1370)

    if load_rule:
        matcher.add(pattern1370, 1370)
    pattern1371 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons163, cons745)
    def replacement1371(d, p, f, e, n2, m, x, b, c, a, n):
        rubi.append(1371)
        return Dist(d**(S(-2)), Int((f*x)**m*(a*d + x**n*(-a*e + b*d))*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Dist(f**(-S(2)*n)*(a*e**S(2) - b*d*e + c*d**S(2))/d**S(2), Int((f*x)**(m + S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)
    rule1371 = ReplacementRule(pattern1371, replacement1371)

    if load_rule:
        matcher.add(pattern1371, 1371)
    pattern1372 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons163, cons745)
    def replacement1372(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1372)
        return Dist(a/d**S(2), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*(d - e*x**n), x), x) + Dist(f**(-S(2)*n)*(a*e**S(2) + c*d**S(2))/d**S(2), Int((f*x)**(m + S(2)*n)*(a + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)
    rule1372 = ReplacementRule(pattern1372, replacement1372)

    if load_rule:
        matcher.add(pattern1372, 1372)
    pattern1373 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons163, cons267)
    def replacement1373(d, p, f, e, n2, m, x, b, c, a, n):
        rubi.append(1373)
        return Dist(S(1)/(d*e), Int((f*x)**m*(a*e + c*d*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) - Dist(f**(-n)*(a*e**S(2) - b*d*e + c*d**S(2))/(d*e), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)
    rule1373 = ReplacementRule(pattern1373, replacement1373)

    if load_rule:
        matcher.add(pattern1373, 1373)
    pattern1374 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons163, cons267)
    def replacement1374(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1374)
        return Dist(S(1)/(d*e), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*(a*e + c*d*x**n), x), x) - Dist(f**(-n)*(a*e**S(2) + c*d**S(2))/(d*e), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)
    rule1374 = ReplacementRule(pattern1374, replacement1374)

    if load_rule:
        matcher.add(pattern1374, 1374)
    pattern1375 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons137, cons746)
    def replacement1375(d, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1375)
        return -Dist(f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a*d + x**n*(-a*e + b*d))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)
    rule1375 = ReplacementRule(pattern1375, replacement1375)

    if load_rule:
        matcher.add(pattern1375, 1375)
    pattern1376 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons137, cons746)
    def replacement1376(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1376)
        return -Dist(a*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + c*x**(S(2)*n))**p*(d - e*x**n), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)
    rule1376 = ReplacementRule(pattern1376, replacement1376)

    if load_rule:
        matcher.add(pattern1376, 1376)
    pattern1377 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons46, cons226, cons148, cons244, cons137, cons168)
    def replacement1377(d, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1377)
        return Dist(f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(a*e + c*d*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) - Dist(d*e*f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)
    rule1377 = ReplacementRule(pattern1377, replacement1377)

    if load_rule:
        matcher.add(pattern1377, 1377)
    pattern1378 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons7, cons27, cons48, cons125, cons46, cons148, cons244, cons137, cons168)
    def replacement1378(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1378)
        return Dist(f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**p*(a*e + c*d*x**n), x), x) - Dist(d*e*f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)
    rule1378 = ReplacementRule(pattern1378, replacement1378)

    if load_rule:
        matcher.add(pattern1378, 1378)
    pattern1379 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons226, cons148, cons747)
    def replacement1379(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1379)
        return Int(ExpandIntegrand((a + b*x**n + c*x**(S(2)*n))**p, (f*x)**m*(d + e*x**n)**q, x), x)
    rule1379 = ReplacementRule(pattern1379, replacement1379)

    if load_rule:
        matcher.add(pattern1379, 1379)
    pattern1380 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons50, cons46, cons148, cons747)
    def replacement1380(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1380)
        return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p, (f*x)**m*(d + e*x**n)**q, x), x)
    rule1380 = ReplacementRule(pattern1380, replacement1380)

    if load_rule:
        matcher.add(pattern1380, 1380)
    pattern1381 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons46, cons226, cons196, cons17)
    def replacement1381(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1381)
        return -Subst(Int(x**(-m + S(-2))*(d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)
    rule1381 = ReplacementRule(pattern1381, replacement1381)

    if load_rule:
        matcher.add(pattern1381, 1381)
    pattern1382 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons5, cons50, cons46, cons196, cons17)
    def replacement1382(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1382)
        return -Subst(Int(x**(-m + S(-2))*(a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q, x), x, S(1)/x)
    rule1382 = ReplacementRule(pattern1382, replacement1382)

    if load_rule:
        matcher.add(pattern1382, 1382)
    def With1383(d, q, p, f, e, n2, x, b, c, a, m, n):
        g = Denominator(m)
        rubi.append(1383)
        return -Dist(g/f, Subst(Int(x**(-g*(m + S(1)) + S(-1))*(d + e*f**(-n)*x**(-g*n))**q*(a + b*f**(-n)*x**(-g*n) + c*f**(-S(2)*n)*x**(-S(2)*g*n))**p, x), x, (f*x)**(-S(1)/g)), x)
    pattern1383 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons50, cons46, cons226, cons196, cons367)
    rule1383 = ReplacementRule(pattern1383, With1383)

    if load_rule:
        matcher.add(pattern1383, 1383)
    def With1384(d, q, p, f, e, n2, x, c, a, m, n):
        g = Denominator(m)
        rubi.append(1384)
        return -Dist(g/f, Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + c*f**(-S(2)*n)*x**(-S(2)*g*n))**p*(d + e*f**(-n)*x**(-g*n))**q, x), x, (f*x)**(-S(1)/g)), x)
    pattern1384 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons5, cons50, cons46, cons196, cons367)
    rule1384 = ReplacementRule(pattern1384, With1384)

    if load_rule:
        matcher.add(pattern1384, 1384)
    pattern1385 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons226, cons196, cons356)
    def replacement1385(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1385)
        return -Dist(f**IntPart(m)*(f*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)
    rule1385 = ReplacementRule(pattern1385, replacement1385)

    if load_rule:
        matcher.add(pattern1385, 1385)
    pattern1386 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons196, cons356)
    def replacement1386(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1386)
        return -Dist(f**IntPart(m)*(f*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q, x), x, S(1)/x), x)
    rule1386 = ReplacementRule(pattern1386, replacement1386)

    if load_rule:
        matcher.add(pattern1386, 1386)
    def With1387(d, q, p, e, n2, x, b, c, a, m, n):
        g = Denominator(n)
        rubi.append(1387)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(d + e*x**(g*n))**q*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p, x), x, x**(S(1)/g)), x)
    pattern1387 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons5, cons50, cons46, cons226, cons489)
    rule1387 = ReplacementRule(pattern1387, With1387)

    if load_rule:
        matcher.add(pattern1387, 1387)
    def With1388(d, q, p, e, n2, x, c, a, m, n):
        g = Denominator(n)
        rubi.append(1388)
        return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + c*x**(S(2)*g*n))**p*(d + e*x**(g*n))**q, x), x, x**(S(1)/g)), x)
    pattern1388 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons5, cons50, cons46, cons489)
    rule1388 = ReplacementRule(pattern1388, With1388)

    if load_rule:
        matcher.add(pattern1388, 1388)
    pattern1389 = Pattern(Integral((f_*x_)**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons226, cons489)
    def replacement1389(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1389)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1389 = ReplacementRule(pattern1389, replacement1389)

    if load_rule:
        matcher.add(pattern1389, 1389)
    pattern1390 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons489)
    def replacement1390(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1390)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1390 = ReplacementRule(pattern1390, replacement1390)

    if load_rule:
        matcher.add(pattern1390, 1390)
    pattern1391 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons226, cons541, cons23)
    def replacement1391(d, q, p, e, n2, x, b, c, a, m, n):
        rubi.append(1391)
        return Dist(S(1)/(m + S(1)), Subst(Int((d + e*x**(n/(m + S(1))))**q*(a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p, x), x, x**(m + S(1))), x)
    rule1391 = ReplacementRule(pattern1391, replacement1391)

    if load_rule:
        matcher.add(pattern1391, 1391)
    pattern1392 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons541, cons23)
    def replacement1392(d, q, p, e, n2, x, c, a, m, n):
        rubi.append(1392)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + c*x**(S(2)*n/(m + S(1))))**p*(d + e*x**(n/(m + S(1))))**q, x), x, x**(m + S(1))), x)
    rule1392 = ReplacementRule(pattern1392, replacement1392)

    if load_rule:
        matcher.add(pattern1392, 1392)
    pattern1393 = Pattern(Integral((f_*x_)**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons226, cons541, cons23)
    def replacement1393(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1393)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1393 = ReplacementRule(pattern1393, replacement1393)

    if load_rule:
        matcher.add(pattern1393, 1393)
    pattern1394 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons5, cons50, cons46, cons541, cons23)
    def replacement1394(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1394)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1394 = ReplacementRule(pattern1394, replacement1394)

    if load_rule:
        matcher.add(pattern1394, 1394)
    def With1395(d, q, f, e, n2, x, b, c, a, m, n):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1395)
        return Dist(S(2)*c/r, Int((f*x)**m*(d + e*x**n)**q/(b + S(2)*c*x**n - r), x), x) - Dist(S(2)*c/r, Int((f*x)**m*(d + e*x**n)**q/(b + S(2)*c*x**n + r), x), x)
    pattern1395 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons50, cons46, cons226)
    rule1395 = ReplacementRule(pattern1395, With1395)

    if load_rule:
        matcher.add(pattern1395, 1395)
    def With1396(d, q, f, e, n2, x, c, a, m, n):
        r = Rt(-a*c, S(2))
        rubi.append(1396)
        return -Dist(c/(S(2)*r), Int((f*x)**m*(d + e*x**n)**q/(-c*x**n + r), x), x) - Dist(c/(S(2)*r), Int((f*x)**m*(d + e*x**n)**q/(c*x**n + r), x), x)
    pattern1396 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons50, cons46)
    rule1396 = ReplacementRule(pattern1396, With1396)

    if load_rule:
        matcher.add(pattern1396, 1396)
    pattern1397 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons46, cons226, cons702)
    def replacement1397(d, p, c, f, e, n2, x, b, a, m, n):
        rubi.append(1397)
        return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e*(m + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(m + n*(S(2)*p + S(3)) + S(1)) + d*(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1))), x), x), x) - Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e + c*x**n*(-S(2)*a*e + b*d) + d*(-S(2)*a*c + b**S(2)))/(a*f*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1397 = ReplacementRule(pattern1397, replacement1397)

    if load_rule:
        matcher.add(pattern1397, 1397)
    pattern1398 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons46, cons702)
    def replacement1398(d, p, f, e, n2, x, c, a, m, n):
        rubi.append(1398)
        return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(1))*Simp(d*(m + S(2)*n*(p + S(1)) + S(1)) + e*x**n*(m + n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*f*n*(p + S(1))), x)
    rule1398 = ReplacementRule(pattern1398, replacement1398)

    if load_rule:
        matcher.add(pattern1398, 1398)
    pattern1399 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons226, cons748)
    def replacement1399(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1399)
        return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1399 = ReplacementRule(pattern1399, replacement1399)

    if load_rule:
        matcher.add(pattern1399, 1399)
    pattern1400 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons748)
    def replacement1400(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1400)
        return Int(ExpandIntegrand((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1400 = ReplacementRule(pattern1400, replacement1400)

    if load_rule:
        matcher.add(pattern1400, 1400)
    pattern1401 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons564, cons733)
    def replacement1401(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1401)
        return Dist(f**m, Int(ExpandIntegrand(x**m*(a + c*x**(S(2)*n))**p, (d/(d**S(2) - e**S(2)*x**(S(2)*n)) - e*x**n/(d**S(2) - e**S(2)*x**(S(2)*n)))**(-q), x), x), x)
    rule1401 = ReplacementRule(pattern1401, replacement1401)

    if load_rule:
        matcher.add(pattern1401, 1401)
    pattern1402 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46, cons564, cons749)
    def replacement1402(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1402)
        return Dist(x**(-m)*(f*x)**m, Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)
    rule1402 = ReplacementRule(pattern1402, replacement1402)

    if load_rule:
        matcher.add(pattern1402, 1402)
    pattern1403 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46)
    def replacement1403(d, q, p, f, e, n2, x, b, c, a, m, n):
        rubi.append(1403)
        return Int((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1403 = ReplacementRule(pattern1403, replacement1403)

    if load_rule:
        matcher.add(pattern1403, 1403)
    pattern1404 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons46)
    def replacement1404(d, q, p, f, e, n2, x, c, a, m, n):
        rubi.append(1404)
        return Int((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x)
    rule1404 = ReplacementRule(pattern1404, replacement1404)

    if load_rule:
        matcher.add(pattern1404, 1404)
    pattern1405 = Pattern(Integral(u_**WC('m', S(1))*(d_ + v_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons46, cons554)
    def replacement1405(d, q, p, v, e, u, n2, x, b, c, a, m, n):
        rubi.append(1405)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x, v), x)
    rule1405 = ReplacementRule(pattern1405, replacement1405)

    if load_rule:
        matcher.add(pattern1405, 1405)
    pattern1406 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + v_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons4, cons5, cons46, cons554)
    def replacement1406(d, q, p, v, e, u, n2, x, c, a, m, n):
        rubi.append(1406)
        return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x, v), x)
    rule1406 = ReplacementRule(pattern1406, replacement1406)

    if load_rule:
        matcher.add(pattern1406, 1406)
    pattern1407 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons680, cons585, cons586)
    def replacement1407(mn, d, q, p, c, e, m, n2, x, b, a, n):
        rubi.append(1407)
        return Int(x**(m - n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1407 = ReplacementRule(pattern1407, replacement1407)

    if load_rule:
        matcher.add(pattern1407, 1407)
    pattern1408 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons21, cons726, cons5, cons725, cons586)
    def replacement1408(mn, d, q, p, e, n2, x, c, a, m):
        rubi.append(1408)
        return Int(x**(m + mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x)
    rule1408 = ReplacementRule(pattern1408, replacement1408)

    if load_rule:
        matcher.add(pattern1408, 1408)
    pattern1409 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons50, cons680, cons585, cons386, cons38)
    def replacement1409(mn, d, q, p, c, e, m, n2, x, b, a, n):
        rubi.append(1409)
        return Int(x**(m + S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)
    rule1409 = ReplacementRule(pattern1409, replacement1409)

    if load_rule:
        matcher.add(pattern1409, 1409)
    pattern1410 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons21, cons726, cons50, cons725, cons386, cons38)
    def replacement1410(mn, d, q, p, e, n2, x, c, a, m):
        rubi.append(1410)
        return Int(x**(m - S(2)*mn*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x)
    rule1410 = ReplacementRule(pattern1410, replacement1410)

    if load_rule:
        matcher.add(pattern1410, 1410)
    pattern1411 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons680, cons585, cons386, cons147)
    def replacement1411(mn, d, q, p, c, e, m, n2, x, b, a, n):
        rubi.append(1411)
        return Dist(x**(n*FracPart(q))*(d + e*x**(-n))**FracPart(q)*(d*x**n + e)**(-FracPart(q)), Int(x**(m - n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1411 = ReplacementRule(pattern1411, replacement1411)

    if load_rule:
        matcher.add(pattern1411, 1411)
    pattern1412 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons7, cons27, cons48, cons21, cons726, cons5, cons50, cons725, cons386, cons147)
    def replacement1412(mn, d, q, p, e, n2, x, c, a, m):
        rubi.append(1412)
        return Dist(x**(-mn*FracPart(q))*(d + e*x**mn)**FracPart(q)*(d*x**(-mn) + e)**(-FracPart(q)), Int(x**(m + mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x), x)
    rule1412 = ReplacementRule(pattern1412, replacement1412)

    if load_rule:
        matcher.add(pattern1412, 1412)
    pattern1413 = Pattern(Integral((f_*x_)**m_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons680, cons585)
    def replacement1413(mn, d, q, p, c, f, e, n2, m, x, b, a, n):
        rubi.append(1413)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**mn)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1413 = ReplacementRule(pattern1413, replacement1413)

    if load_rule:
        matcher.add(pattern1413, 1413)
    pattern1414 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons7, cons27, cons48, cons125, cons21, cons726, cons5, cons50, cons725)
    def replacement1414(mn, d, q, p, f, e, n2, x, c, a, m):
        rubi.append(1414)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**mn)**q, x), x)
    rule1414 = ReplacementRule(pattern1414, replacement1414)

    if load_rule:
        matcher.add(pattern1414, 1414)
    pattern1415 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons50, cons585, cons38)
    def replacement1415(mn, d, q, p, e, x, b, c, a, m, n):
        rubi.append(1415)
        return Int(x**(m - n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x)
    rule1415 = ReplacementRule(pattern1415, replacement1415)

    if load_rule:
        matcher.add(pattern1415, 1415)
    pattern1416 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons50, cons585, cons147)
    def replacement1416(mn, d, q, p, e, x, b, c, a, m, n):
        rubi.append(1416)
        return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(m - n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x), x)
    rule1416 = ReplacementRule(pattern1416, replacement1416)

    if load_rule:
        matcher.add(pattern1416, 1416)
    pattern1417 = Pattern(Integral((f_*x_)**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons585)
    def replacement1417(mn, d, q, p, f, e, x, b, c, a, m, n):
        rubi.append(1417)
        return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**(-n) + c*x**n)**p, x), x)
    rule1417 = ReplacementRule(pattern1417, replacement1417)

    if load_rule:
        matcher.add(pattern1417, 1417)
    pattern1418 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons731, cons652, cons732, cons654, cons125, cons4, cons5, cons50, cons46, cons593, cons729, cons730)
    def replacement1418(q, p, c, f, non2, e2, d1, d2, n2, x, e1, b, a, m, n):
        rubi.append(1418)
        return Int((f*x)**m*(d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)
    rule1418 = ReplacementRule(pattern1418, replacement1418)

    if load_rule:
        matcher.add(pattern1418, 1418)
    pattern1419 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons731, cons652, cons732, cons654, cons125, cons4, cons5, cons50, cons46, cons593, cons729)
    def replacement1419(q, p, c, f, non2, e2, d1, d2, n2, x, e1, b, a, m, n):
        rubi.append(1419)
        return Dist((d1 + e1*x**(n/S(2)))**FracPart(q)*(d2 + e2*x**(n/S(2)))**FracPart(q)*(d1*d2 + e1*e2*x**n)**(-FracPart(q)), Int((f*x)**m*(d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)
    rule1419 = ReplacementRule(pattern1419, replacement1419)

    if load_rule:
        matcher.add(pattern1419, 1419)
    pattern1420 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons750, cons751)
    def replacement1420(q, p, c, x, r, b, a, n):
        rubi.append(1420)
        return Int((x**n*(a + b + c))**p, x)
    rule1420 = ReplacementRule(pattern1420, replacement1420)

    if load_rule:
        matcher.add(pattern1420, 1420)
    pattern1421 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons50, cons752, cons753, cons38)
    def replacement1421(q, p, c, x, r, b, a, n):
        rubi.append(1421)
        return Int(x**(p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)
    rule1421 = ReplacementRule(pattern1421, replacement1421)

    if load_rule:
        matcher.add(pattern1421, 1421)
    pattern1422 = Pattern(Integral(sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons4, cons50, cons752, cons753)
    def replacement1422(q, c, x, r, b, a, n):
        rubi.append(1422)
        return Dist(x**(-q/S(2))*sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), Int(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)
    rule1422 = ReplacementRule(pattern1422, replacement1422)

    if load_rule:
        matcher.add(pattern1422, 1422)
    pattern1423 = Pattern(Integral(S(1)/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons4, cons50, cons752, cons753)
    def replacement1423(q, c, x, r, b, a, n):
        rubi.append(1423)
        return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(-q/S(2))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)
    rule1423 = ReplacementRule(pattern1423, replacement1423)

    if load_rule:
        matcher.add(pattern1423, 1423)
    pattern1424 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons50, cons752, cons753, cons147, cons226, cons13, cons163, cons754)
    def replacement1424(q, p, c, x, r, b, a, n):
        rubi.append(1424)
        return Dist(p*(n - q)/(p*(S(2)*n - q) + S(1)), Int(x**q*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(p*(S(2)*n - q) + S(1)), x)
    rule1424 = ReplacementRule(pattern1424, replacement1424)

    if load_rule:
        matcher.add(pattern1424, 1424)
    pattern1425 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons50, cons752, cons753, cons147, cons226, cons13, cons137)
    def replacement1425(q, p, c, x, r, b, a, n):
        rubi.append(1425)
        return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(-q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(b*c*x**(n - q)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1)) + (n - q)*(p + S(1))*(-S(4)*a*c + b**S(2)) + (-S(2)*a*c + b**S(2))*(p*q + S(1))), x), x) - Simp(x**(-q + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1425 = ReplacementRule(pattern1425, replacement1425)

    if load_rule:
        matcher.add(pattern1425, 1425)
    pattern1426 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons50, cons752, cons753, cons147)
    def replacement1426(q, p, c, x, r, b, a, n):
        rubi.append(1426)
        return Dist(x**(-p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**(-p)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, Int(x**(p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x), x)
    rule1426 = ReplacementRule(pattern1426, replacement1426)

    if load_rule:
        matcher.add(pattern1426, 1426)
    pattern1427 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons50, cons752)
    def replacement1427(q, p, c, x, r, b, a, n):
        rubi.append(1427)
        return Int((a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x)
    rule1427 = ReplacementRule(pattern1427, replacement1427)

    if load_rule:
        matcher.add(pattern1427, 1427)
    pattern1428 = Pattern(Integral((u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons50, cons752, cons68, cons69)
    def replacement1428(q, p, c, u, x, r, b, a, n):
        rubi.append(1428)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)
    rule1428 = ReplacementRule(pattern1428, replacement1428)

    if load_rule:
        matcher.add(pattern1428, 1428)
    pattern1429 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons755, cons751)
    def replacement1429(q, p, c, m, x, r, b, a, n):
        rubi.append(1429)
        return Int(x**m*(x**n*(a + b + c))**p, x)
    rule1429 = ReplacementRule(pattern1429, replacement1429)

    if load_rule:
        matcher.add(pattern1429, 1429)
    pattern1430 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons50, cons752, cons38, cons753)
    def replacement1430(q, p, c, m, x, r, b, a, n):
        rubi.append(1430)
        return Int(x**(m + p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)
    rule1430 = ReplacementRule(pattern1430, replacement1430)

    if load_rule:
        matcher.add(pattern1430, 1430)
    pattern1431 = Pattern(Integral(x_**WC('m', S(1))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons50, cons752, cons753, cons756)
    def replacement1431(q, c, m, x, r, b, a, n):
        rubi.append(1431)
        return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(m - q/S(2))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)
    rule1431 = ReplacementRule(pattern1431, replacement1431)

    if load_rule:
        matcher.add(pattern1431, 1431)
    pattern1432 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons4, cons757, cons758, cons759, cons226)
    def replacement1432(q, c, m, x, r, b, a, n):
        rubi.append(1432)
        return Simp(-S(2)*x**(n/S(2) + S(-1)/2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))), x)
    rule1432 = ReplacementRule(pattern1432, replacement1432)

    if load_rule:
        matcher.add(pattern1432, 1432)
    pattern1433 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons7, cons4, cons760, cons758, cons759, cons226)
    def replacement1433(q, c, m, x, r, b, a, n):
        rubi.append(1433)
        return Simp(x**(n/S(2) + S(-1)/2)*(S(4)*a + S(2)*b*x)/((-S(4)*a*c + b**S(2))*sqrt(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))), x)
    rule1433 = ReplacementRule(pattern1433, replacement1433)

    if load_rule:
        matcher.add(pattern1433, 1433)
    pattern1434 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons761)
    def replacement1434(q, p, c, m, x, r, b, a, n):
        rubi.append(1434)
        return -Dist(b/(S(2)*c), Int(x**(m + S(-1))*(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))**p, x), x) + Simp(x**(m - n)*(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))**(p + S(1))/(S(2)*c*(p + S(1))), x)
    rule1434 = ReplacementRule(pattern1434, replacement1434)

    if load_rule:
        matcher.add(pattern1434, 1434)
    pattern1435 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons163, cons762)
    def replacement1435(q, p, c, m, x, r, b, a, n):
        rubi.append(1435)
        return -Dist(p*(-S(4)*a*c + b**S(2))/(S(2)*c*(S(2)*p + S(1))), Int(x**(m + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m - n + q + S(1))*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(S(2)*c*(n - q)*(S(2)*p + S(1))), x)
    rule1435 = ReplacementRule(pattern1435, replacement1435)

    if load_rule:
        matcher.add(pattern1435, 1435)
    pattern1436 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons163, cons763, cons764, cons765)
    def replacement1436(q, p, c, m, x, r, b, a, n):
        rubi.append(1436)
        return Dist(p*(n - q)/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1))), Int(x**(m - n + S(2)*q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-a*b*(m - n + p*q + q + S(1)) + x**(n - q)*(S(2)*a*c*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1)) - b**S(2)*(m + p*q + (n - q)*(p + S(-1)) + S(1))), x), x), x) + Simp(x**(m - n + q + S(1))*(b*p*(n - q) + c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1))), x)
    rule1436 = ReplacementRule(pattern1436, replacement1436)

    if load_rule:
        matcher.add(pattern1436, 1436)
    pattern1437 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons163, cons766, cons767)
    def replacement1437(q, p, c, m, x, r, b, a, n):
        rubi.append(1437)
        return -Dist(p*(n - q)/(m + p*q + S(1)), Int(x**(m + n)*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(m + p*q + S(1)), x)
    rule1437 = ReplacementRule(pattern1437, replacement1437)

    if load_rule:
        matcher.add(pattern1437, 1437)
    pattern1438 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons163, cons768, cons764)
    def replacement1438(q, p, c, m, x, r, b, a, n):
        rubi.append(1438)
        return Dist(p*(n - q)/(m + p*(S(2)*n - q) + S(1)), Int(x**(m + q)*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(m + p*(S(2)*n - q) + S(1)), x)
    rule1438 = ReplacementRule(pattern1438, replacement1438)

    if load_rule:
        matcher.add(pattern1438, 1438)
    pattern1439 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons137, cons769)
    def replacement1439(q, p, c, m, x, r, b, a, n):
        rubi.append(1439)
        return Dist((S(2)*a*c - b**S(2)*(p + S(2)))/(a*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) - Simp(x**(m - q + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1439 = ReplacementRule(pattern1439, replacement1439)

    if load_rule:
        matcher.add(pattern1439, 1439)
    pattern1440 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons137, cons770)
    def replacement1440(q, p, c, m, x, r, b, a, n):
        rubi.append(1440)
        return Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - S(2)*n + q)*(S(2)*a*(m - S(2)*n + p*q + S(2)*q + S(1)) + b*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) - Simp(x**(m - S(2)*n + q + S(1))*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1440 = ReplacementRule(pattern1440, replacement1440)

    if load_rule:
        matcher.add(pattern1440, 1440)
    pattern1441 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons137, cons771)
    def replacement1441(q, p, c, m, x, r, b, a, n):
        rubi.append(1441)
        return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*a*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + b**S(2)*(m + p*q + (n - q)*(p + S(1)) + S(1)) + b*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(m - q + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1441 = ReplacementRule(pattern1441, replacement1441)

    if load_rule:
        matcher.add(pattern1441, 1441)
    pattern1442 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons137, cons772)
    def replacement1442(q, p, c, m, x, r, b, a, n):
        rubi.append(1442)
        return -Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - n)*(b*(m - n + p*q + q + S(1)) + S(2)*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) + Simp(x**(m - n + S(1))*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1442 = ReplacementRule(pattern1442, replacement1442)

    if load_rule:
        matcher.add(pattern1442, 1442)
    pattern1443 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons773, cons774)
    def replacement1443(q, p, c, m, x, r, b, a, n):
        rubi.append(1443)
        return -Dist(b/(S(2)*c), Int(x**(m - n + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - S(2)*n + q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*c*(n - q)*(p + S(1))), x)
    rule1443 = ReplacementRule(pattern1443, replacement1443)

    if load_rule:
        matcher.add(pattern1443, 1443)
    pattern1444 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons773, cons775)
    def replacement1444(q, p, c, m, x, r, b, a, n):
        rubi.append(1444)
        return -Dist(b/(S(2)*a), Int(x**(m + n - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) - Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*a*(n - q)*(p + S(1))), x)
    rule1444 = ReplacementRule(pattern1444, replacement1444)

    if load_rule:
        matcher.add(pattern1444, 1444)
    pattern1445 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons773, cons770)
    def replacement1445(q, p, c, m, x, r, b, a, n):
        rubi.append(1445)
        return -Dist(S(1)/(c*(m + p*q + S(2)*p*(n - q) + S(1))), Int(x**(m - S(2)*n + S(2)*q)*(a*(m - S(2)*n + p*q + S(2)*q + S(1)) + b*x**(n - q)*(m + p*q + (n - q)*(p + S(-1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - S(2)*n + q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + S(2)*p*(n - q) + S(1))), x)
    rule1445 = ReplacementRule(pattern1445, replacement1445)

    if load_rule:
        matcher.add(pattern1445, 1445)
    pattern1446 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons752, cons753, cons147, cons226, cons148, cons606, cons773, cons776)
    def replacement1446(q, p, c, m, x, r, b, a, n):
        rubi.append(1446)
        return -Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(b*(m + p*q + (n - q)*(p + S(1)) + S(1)) + c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)
    rule1446 = ReplacementRule(pattern1446, replacement1446)

    if load_rule:
        matcher.add(pattern1446, 1446)
    pattern1447 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons50, cons752, cons147, cons753)
    def replacement1447(q, p, c, m, x, r, b, a, n):
        rubi.append(1447)
        return Dist(x**(-p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**(-p)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, Int(x**(m + p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x), x)
    rule1447 = ReplacementRule(pattern1447, replacement1447)

    if load_rule:
        matcher.add(pattern1447, 1447)
    pattern1448 = Pattern(Integral(u_**WC('m', S(1))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons50, cons752, cons68, cons69)
    def replacement1448(q, p, c, u, m, x, r, b, a, n):
        rubi.append(1448)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int(x**m*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)
    rule1448 = ReplacementRule(pattern1448, replacement1448)

    if load_rule:
        matcher.add(pattern1448, 1448)
    pattern1449 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons50, cons777, cons778, cons38, cons753)
    def replacement1449(B, q, p, j, A, x, r, b, c, a, n):
        rubi.append(1449)
        return Int(x**(p*q)*(A + B*x**(n - q))*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)
    rule1449 = ReplacementRule(pattern1449, replacement1449)

    if load_rule:
        matcher.add(pattern1449, 1449)
    pattern1450 = Pattern(Integral((A_ + x_**WC('j', S(1))*WC('B', S(1)))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons50, cons779, cons752, cons753, cons780, cons781)
    def replacement1450(B, q, j, c, A, x, r, b, a, n):
        rubi.append(1450)
        return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(-q/S(2))*(A + B*x**(n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)
    rule1450 = ReplacementRule(pattern1450, replacement1450)

    if load_rule:
        matcher.add(pattern1450, 1450)
    pattern1451 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons34, cons35, cons4, cons50, cons777, cons778, cons147, cons226, cons13, cons163, cons754, cons782)
    def replacement1451(B, q, p, j, A, x, r, b, c, a, n):
        rubi.append(1451)
        return Dist(p*(n - q)/(c*(p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**q*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*(S(2)*A*a*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) - B*a*b*(p*q + S(1)) + x**(n - q)*(A*b*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*c*(p*(S(2)*n - q) + S(1)) - B*b**S(2)*(p*q + p*(n - q) + S(1)))), x), x) + Simp(x*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*(A*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*p*(n - q) + B*c*x**(n - q)*(p*(S(2)*n - q) + S(1)))/(c*(p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1451 = ReplacementRule(pattern1451, replacement1451)

    if load_rule:
        matcher.add(pattern1451, 1451)

    pattern1452 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons7, cons34, cons35, cons50, cons147, cons13, cons163, cons_with_1452)
    def replacement1452(B, q, p, j, A, x, r, c, a):

        n = q + r
        rubi.append(1452)
        return Dist(p*(n - q)/((p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**q*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*(S(2)*A*a*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*x**(n - q)*(p*(S(2)*n - q) + S(1))), x), x) + Simp(x*(A*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(p*(S(2)*n - q) + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1452 = ReplacementRule(pattern1452, replacement1452)

    if load_rule:
        matcher.add(pattern1452, 1452)
    pattern1453 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons7, cons34, cons35, cons4, cons50, cons777, cons778, cons147, cons226, cons13, cons137)
    def replacement1453(B, q, p, j, A, x, r, b, c, a, n):
        rubi.append(1453)
        return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(-q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c*(p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + A*b**S(2)*(p*q + (n - q)*(p + S(1)) + S(1)) - B*a*b*(p*q + S(1)) + c*x**(n - q)*(A*b - S(2)*B*a)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(-q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c + A*b**S(2) - B*a*b + c*x**(n - q)*(A*b - S(2)*B*a))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1453 = ReplacementRule(pattern1453, replacement1453)

    if load_rule:
        matcher.add(pattern1453, 1453)

    pattern1454 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons7, cons34, cons35, cons50, cons147, cons13, cons137, cons_with_1454)
    def replacement1454(B, q, p, j, A, x, r, c, a):

        n = q + r
        rubi.append(1454)
        return Dist(S(1)/(S(2)*a**S(2)*c*(n - q)*(p + S(1))), Int(x**(-q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(A*a*c*(p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*c*x**(n - q)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(-q + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(A*a*c + B*a*c*x**(n - q))/(S(2)*a**S(2)*c*(n - q)*(p + S(1))), x)
    rule1454 = ReplacementRule(pattern1454, replacement1454)

    if load_rule:
        matcher.add(pattern1454, 1454)
    pattern1455 = Pattern(Integral((A_ + x_**WC('j', S(1))*WC('B', S(1)))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons50, cons779, cons752)
    def replacement1455(B, q, p, j, c, A, x, r, b, a, n):
        rubi.append(1455)
        return Int((A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x)
    rule1455 = ReplacementRule(pattern1455, replacement1455)

    if load_rule:
        matcher.add(pattern1455, 1455)
    pattern1456 = Pattern(Integral((A_ + u_**WC('j', S(1))*WC('B', S(1)))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons50, cons779, cons752, cons68, cons69)
    def replacement1456(B, q, p, j, c, u, A, x, r, b, a, n):
        rubi.append(1456)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)
    rule1456 = ReplacementRule(pattern1456, replacement1456)

    if load_rule:
        matcher.add(pattern1456, 1456)
    pattern1457 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons50, cons777, cons778, cons38, cons753)
    def replacement1457(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1457)
        return Int(x**(m + p*q)*(A + B*x**(n - q))*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)
    rule1457 = ReplacementRule(pattern1457, replacement1457)

    if load_rule:
        matcher.add(pattern1457, 1457)
    pattern1458 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons163, cons783, cons784, cons785)
    def replacement1458(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1458)
        return Dist(p*(n - q)/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + n)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-A*b*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*(m + p*q + S(1)) + x**(n - q)*(-S(2)*A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*(m + p*q + S(1))), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1458 = ReplacementRule(pattern1458, replacement1458)

    if load_rule:
        matcher.add(pattern1458, 1458)

    pattern1459 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons163, cons_with_1459)
    def replacement1459(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1459)
        return Dist(S(2)*p*(n - q)/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + n)*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-A*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1459 = ReplacementRule(pattern1459, replacement1459)

    if load_rule:
        matcher.add(pattern1459, 1459)
    pattern1460 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons137, cons786)
    def replacement1460(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1460)
        return Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - n)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*Simp(x**(n - q)*(-S(2)*A*c + B*b)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + (-A*b + S(2)*B*a)*(m - n + p*q + q + S(1)), x), x), x) + Simp(x**(m - n + S(1))*(A*b - S(2)*B*a - x**(n - q)*(-S(2)*A*c + B*b))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1460 = ReplacementRule(pattern1460, replacement1460)

    if load_rule:
        matcher.add(pattern1460, 1460)

    pattern1461 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons137, cons_with_1461)
    def replacement1461(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1461)
        return -Dist(S(1)/(S(2)*a*c*(n - q)*(p + S(1))), Int(x**(m - n)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*Simp(-A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m - n + p*q + q + S(1)), x), x), x) + Simp(x**(m - n + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(-A*c*x**(n - q) + B*a)/(S(2)*a*c*(n - q)*(p + S(1))), x)
    rule1461 = ReplacementRule(pattern1461, replacement1461)

    if load_rule:
        matcher.add(pattern1461, 1461)
    pattern1462 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons163, cons787, cons764, cons785)
    def replacement1462(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1462)
        return Dist(p*(n - q)/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(S(2)*A*a*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) - B*a*b*(m + p*q + S(1)) + x**(n - q)*(A*b*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*c*(m + p*q + S(2)*p*(n - q) + S(1)) - B*b**S(2)*(m + p*q + p*(n - q) + S(1))), x), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*(A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*p*(n - q) + B*c*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)))/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1462 = ReplacementRule(pattern1462, replacement1462)

    if load_rule:
        matcher.add(pattern1462, 1462)

    pattern1463 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons163, cons_with_1463)
    def replacement1463(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1463)
        return Dist(p*(n - q)/((m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*Simp(S(2)*A*a*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1463 = ReplacementRule(pattern1463, replacement1463)

    if load_rule:
        matcher.add(pattern1463, 1463)
    pattern1464 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons137, cons788)
    def replacement1464(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1464)
        return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*Simp(-S(2)*A*a*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + A*b**S(2)*(m + p*q + (n - q)*(p + S(1)) + S(1)) - B*a*b*(m + p*q + S(1)) + c*x**(n - q)*(A*b - S(2)*B*a)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c + A*b**S(2) - B*a*b + c*x**(n - q)*(A*b - S(2)*B*a))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)
    rule1464 = ReplacementRule(pattern1464, replacement1464)

    if load_rule:
        matcher.add(pattern1464, 1464)

    pattern1465 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons137, cons_with_1465)
    def replacement1465(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1465)
        return Dist(S(1)/(S(2)*a*c*(n - q)*(p + S(1))), Int(x**(m - q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*Simp(A*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x**(m - q + S(1))*(A*c + B*c*x**(n - q))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*a*c*(n - q)*(p + S(1))), x)
    rule1465 = ReplacementRule(pattern1465, replacement1465)

    if load_rule:
        matcher.add(pattern1465, 1465)
    pattern1466 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons773, cons789, cons785)
    def replacement1466(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1466)
        return -Dist(S(1)/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m - n + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*Simp(B*a*(m - n + p*q + q + S(1)) + x**(n - q)*(-A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*(m + p*q + p*(n - q) + S(1))), x), x), x) + Simp(B*x**(m - n + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1466 = ReplacementRule(pattern1466, replacement1466)

    if load_rule:
        matcher.add(pattern1466, 1466)

    pattern1467 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons773, cons_with_1467)
    def replacement1467(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1467)
        return -Dist(S(1)/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m - n + q)*(a*x**q + c*x**(S(2)*n - q))**p*Simp(-A*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*a*(m - n + p*q + q + S(1)), x), x), x) + Simp(B*x**(m - n + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)
    rule1467 = ReplacementRule(pattern1467, replacement1467)

    if load_rule:
        matcher.add(pattern1467, 1467)
    pattern1468 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons777, cons778, cons147, cons226, cons148, cons606, cons790, cons783, cons784)
    def replacement1468(B, q, p, j, A, x, r, b, c, a, m, n):
        rubi.append(1468)
        return Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*Simp(-A*b*(m + p*q + (n - q)*(p + S(1)) + S(1)) - A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(A*x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)
    rule1468 = ReplacementRule(pattern1468, replacement1468)

    if load_rule:
        matcher.add(pattern1468, 1468)

    pattern1469 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons34, cons35, cons147, cons606, cons_with_1469)
    def replacement1469(B, q, p, j, m, A, x, r, c, a):

        n = q + r
        rubi.append(1469)
        return Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(a*x**q + c*x**(S(2)*n - q))**p*Simp(-A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(A*x**(m - q + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)
    rule1469 = ReplacementRule(pattern1469, replacement1469)

    if load_rule:
        matcher.add(pattern1469, 1469)
    pattern1470 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('j', S(1))*WC('B', S(1)))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons50, cons779, cons752, cons753, cons791, cons780, cons792)
    def replacement1470(B, q, j, c, A, x, r, b, a, m, n):
        rubi.append(1470)
        return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(m - q/S(2))*(A + B*x**(n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)
    rule1470 = ReplacementRule(pattern1470, replacement1470)

    if load_rule:
        matcher.add(pattern1470, 1470)
    pattern1471 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**q_*WC('B', S(1)))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons34, cons35, cons796, cons797, cons21, cons5, cons793, cons794, cons147, cons795)
    def replacement1471(B, q, p, k, j, c, m, A, x, b, a, n):
        rubi.append(1471)
        return Dist(x**(-j*p)*(a + b*x**(-j + k) + c*x**(-S(2)*j + S(2)*k))**(-p)*(a*x**j + b*x**k + c*x**n)**p, Int(x**(j*p + m)*(A + B*x**(-j + k))*(a + b*x**(-j + k) + c*x**(-S(2)*j + S(2)*k))**p, x), x)
    rule1471 = ReplacementRule(pattern1471, replacement1471)

    if load_rule:
        matcher.add(pattern1471, 1471)
    pattern1472 = Pattern(Integral(u_**WC('m', S(1))*(A_ + u_**WC('j', S(1))*WC('B', S(1)))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons5, cons50, cons779, cons752, cons68, cons69)
    def replacement1472(B, q, p, j, c, u, A, x, r, b, a, m, n):
        rubi.append(1472)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int(x**m*(A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)
    rule1472 = ReplacementRule(pattern1472, replacement1472)

    if load_rule:
        matcher.add(pattern1472, 1472)
    return matcher, [rule1076, rule1077, rule1078, rule1079, rule1080, rule1081, rule1082, rule1083, rule1084, rule1085, rule1086, rule1087, rule1088, rule1089, rule1090, rule1091, rule1092, rule1093, rule1094, rule1095, rule1096, rule1097, rule1098, rule1099, rule1100, rule1101, rule1102, rule1103, rule1104, rule1105, rule1106, rule1107, rule1108, rule1109, rule1110, rule1111, rule1112, rule1113, rule1114, rule1115, rule1116, rule1117, rule1118, rule1119, rule1120, rule1121, rule1122, rule1123, rule1124, rule1125, rule1126, rule1127, rule1128, rule1129, rule1130, rule1131, rule1132, rule1133, rule1134, rule1135, rule1136, rule1137, rule1138, rule1139, rule1140, rule1141, rule1142, rule1143, rule1144, rule1145, rule1146, rule1147, rule1148, rule1149, rule1150, rule1151, rule1152, rule1153, rule1154, rule1155, rule1156, rule1157, rule1158, rule1159, rule1160, rule1161, rule1162, rule1163, rule1164, rule1165, rule1166, rule1167, rule1168, rule1169, rule1170, rule1171, rule1172, rule1173, rule1174, rule1175, rule1176, rule1177, rule1178, rule1179, rule1180, rule1181, rule1182, rule1183, rule1184, rule1185, rule1186, rule1187, rule1188, rule1189, rule1190, rule1191, rule1192, rule1193, rule1194, rule1195, rule1196, rule1197, rule1198, rule1199, rule1200, rule1201, rule1202, rule1203, rule1204, rule1205, rule1206, rule1207, rule1208, rule1209, rule1210, rule1211, rule1212, rule1213, rule1214, rule1215, rule1216, rule1217, rule1218, rule1219, rule1220, rule1221, rule1222, rule1223, rule1224, rule1225, rule1226, rule1227, rule1228, rule1229, rule1230, rule1231, rule1232, rule1233, rule1234, rule1235, rule1236, rule1237, rule1238, rule1239, rule1240, rule1241, rule1242, rule1243, rule1244, rule1245, rule1246, rule1247, rule1248, rule1249, rule1250, rule1251, rule1252, rule1253, rule1254, rule1255, rule1256, rule1257, rule1258, rule1259, rule1260, rule1261, rule1262, rule1263, rule1264, rule1265, rule1266, rule1267, rule1268, rule1269, rule1270, rule1271, rule1272, rule1273, rule1274, rule1275, rule1276, rule1277, rule1278, rule1279, rule1280, rule1281, rule1282, rule1283, rule1284, rule1285, rule1286, rule1287, rule1288, rule1289, rule1290, rule1291, rule1292, rule1293, rule1294, rule1295, rule1296, rule1297, rule1298, rule1299, rule1300, rule1301, rule1302, rule1303, rule1304, rule1305, rule1306, rule1307, rule1308, rule1309, rule1310, rule1311, rule1312, rule1313, rule1314, rule1315, rule1316, rule1317, rule1318, rule1319, rule1320, rule1321, rule1322, rule1323, rule1324, rule1325, rule1326, rule1327, rule1328, rule1329, rule1330, rule1331, rule1332, rule1333, rule1334, rule1335, rule1336, rule1337, rule1338, rule1339, rule1340, rule1341, rule1342, rule1343, rule1344, rule1345, rule1346, rule1347, rule1348, rule1349, rule1350, rule1351, rule1352, rule1353, rule1354, rule1355, rule1356, rule1357, rule1358, rule1359, rule1360, rule1361, rule1362, rule1363, rule1364, rule1365, rule1366, rule1367, rule1368, rule1369, rule1370, rule1371, rule1372, rule1373, rule1374, rule1375, rule1376, rule1377, rule1378, rule1379, rule1380, rule1381, rule1382, rule1383, rule1384, rule1385, rule1386, rule1387, rule1388, rule1389, rule1390, rule1391, rule1392, rule1393, rule1394, rule1395, rule1396, rule1397, rule1398, rule1399, rule1400, rule1401, rule1402, rule1403, rule1404, rule1405, rule1406, rule1407, rule1408, rule1409, rule1410, rule1411, rule1412, rule1413, rule1414, rule1415, rule1416, rule1417, rule1418, rule1419, rule1420, rule1421, rule1422, rule1423, rule1424, rule1425, rule1426, rule1427, rule1428, rule1429, rule1430, rule1431, rule1432, rule1433, rule1434, rule1435, rule1436, rule1437, rule1438, rule1439, rule1440, rule1441, rule1442, rule1443, rule1444, rule1445, rule1446, rule1447, rule1448, rule1449, rule1450, rule1451, rule1452, rule1453, rule1454, rule1455, rule1456, rule1457, rule1458, rule1459, rule1460, rule1461, rule1462, rule1463, rule1464, rule1465, rule1466, rule1467, rule1468, rule1469, rule1470, rule1471, rule1472, ]
