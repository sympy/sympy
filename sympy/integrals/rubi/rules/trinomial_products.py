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
    from sympy.core.add import Add
    from sympy.core.mod import Mod
    from sympy.core.mul import Mul
    from sympy.core import EulerGamma
    from sympy.core.numbers import (Float, I, Integer)
    from sympy.core.power import Pow
    from sympy.core.singleton import S
    from sympy.functions.elementary.complexes import (Abs, sign)
    from sympy.functions.elementary.miscellaneous import sqrt
    from sympy.integrals.integrals import Integral
    from sympy.logic.boolalg import (And, Or)
    from sympy.simplify.simplify import simplify
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec, atan2)
    from sympy.core.numbers import pi as Pi

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii, Pqq, Q, R, r, C, k, u = symbols('i ii Pqq Q R r C k u')
    _UseGamma = False
    ShowSteps = False
    StepCounter = None


def trinomial_products():
    from sympy.integrals.rubi.constraints import cons48, cons89, cons465, cons40, cons2, cons3, cons8, cons491, cons5, cons47, cons149, cons666, cons4, cons667, cons586, cons668, cons13, cons165, cons669, cons316, cons670, cons464, cons198, cons671, cons672, cons148, cons673, cons674, cons340, cons139, cons228, cons130, cons248, cons675, cons676, cons415, cons677, cons295, cons678, cons679, cons486, cons179, cons680, cons681, cons682, cons587, cons683, cons70, cons71, cons55, cons19, cons503, cons29, cons65, cons504, cons684, cons157, cons685, cons686, cons227, cons58, cons245, cons150, cons246, cons687, cons20, cons688, cons689, cons512, cons690, cons691, cons692, cons531, cons33, cons532, cons693, cons96, cons369, cons358, cons502, cons694, cons695, cons696, cons697, cons698, cons699, cons700, cons701, cons702, cons703, cons543, cons25, cons704, cons554, cons555, cons556, cons222, cons50, cons52, cons705, cons706, cons258, cons259, cons281, cons223, cons282, cons397, cons398, cons707, cons708, cons709, cons710, cons711, cons87, cons712, cons713, cons714, cons715, cons588, cons388, cons151, cons716, cons45, cons717, cons450, cons718, cons402, cons719, cons720, cons721, cons349, cons566, cons722, cons270, cons723, cons724, cons725, cons726, cons727, cons728, cons729, cons730, cons127, cons210, cons54, cons595, cons731, cons732, cons733, cons654, cons734, cons656, cons36, cons37, cons735, cons21, cons736, cons466, cons737, cons170, cons269, cons738, cons739, cons740, cons741, cons742, cons83, cons436, cons743, cons744, cons745, cons746, cons613, cons405, cons747, cons748, cons749, cons750, cons751, cons752, cons753, cons754, cons755, cons756, cons757, cons758, cons759, cons760, cons761, cons762, cons608, cons763, cons764, cons765, cons766, cons767, cons768, cons769, cons770, cons771, cons772, cons773, cons774, cons775, cons776, cons777, cons778, cons779, cons780, cons781, cons782, cons783, cons784, cons785, cons786, cons787, cons788, cons789, cons790, cons791, cons792, cons793, cons794, cons795, cons796, cons797, cons798, cons799


    pattern1079 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons48, cons89, cons465, cons40)
    rule1079 = ReplacementRule(pattern1079, replacement1079)

    pattern1080 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons491)
    rule1080 = ReplacementRule(pattern1080, With1080)

    pattern1081 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons47, cons149, cons666)
    rule1081 = ReplacementRule(pattern1081, replacement1081)

    pattern1082 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons47, cons149, cons667, cons586)
    rule1082 = ReplacementRule(pattern1082, replacement1082)

    pattern1083 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons47, cons668, cons13, cons165, cons669)
    rule1083 = ReplacementRule(pattern1083, replacement1083)

    pattern1084 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons47, cons668, cons13, cons165, cons316)
    rule1084 = ReplacementRule(pattern1084, replacement1084)

    pattern1085 = Pattern(Integral(sqrt(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons4, cons48, cons47, cons586, cons670, cons464)
    rule1085 = ReplacementRule(pattern1085, replacement1085)

    pattern1086 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons47, cons149, cons198)
    rule1086 = ReplacementRule(pattern1086, replacement1086)

    pattern1087 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons47, cons149, cons671, cons672, cons13, cons148)
    rule1087 = ReplacementRule(pattern1087, replacement1087)

    pattern1088 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons47, cons149, cons673, cons674, cons340, cons139)
    rule1088 = ReplacementRule(pattern1088, replacement1088)

    pattern1089 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons47, cons149)
    rule1089 = ReplacementRule(pattern1089, replacement1089)

    pattern1090 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons48, cons198)
    rule1090 = ReplacementRule(pattern1090, replacement1090)

    pattern1091 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons228, cons130)
    rule1091 = ReplacementRule(pattern1091, replacement1091)

    pattern1092 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons228, cons13, cons165, cons671, cons248, cons675)
    rule1092 = ReplacementRule(pattern1092, replacement1092)

    pattern1093 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons48, cons228, cons13, cons139, cons248, cons675)
    rule1093 = ReplacementRule(pattern1093, replacement1093)

    pattern1094 = Pattern(Integral(S(1)/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons48, cons228, cons676, cons415)
    rule1094 = ReplacementRule(pattern1094, With1094)

    pattern1095 = Pattern(Integral(S(1)/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons48, cons228)
    rule1095 = ReplacementRule(pattern1095, With1095)

    pattern1096 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons295)
    rule1096 = ReplacementRule(pattern1096, With1096)

    pattern1097 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons678, cons679)
    rule1097 = ReplacementRule(pattern1097, With1097)

    pattern1098 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons486, cons179, CustomConstraint(With1098))
    rule1098 = ReplacementRule(pattern1098, replacement1098)

    pattern1099 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons486, cons179)
    rule1099 = ReplacementRule(pattern1099, With1099)

    pattern1100 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1100))
    rule1100 = ReplacementRule(pattern1100, replacement1100)

    pattern1101 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1101))
    rule1101 = ReplacementRule(pattern1101, replacement1101)

    pattern1102 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1102))
    rule1102 = ReplacementRule(pattern1102, replacement1102)

    pattern1103 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1103))
    rule1103 = ReplacementRule(pattern1103, replacement1103)

    pattern1104 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons680)
    rule1104 = ReplacementRule(pattern1104, With1104)

    pattern1105 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons681)
    rule1105 = ReplacementRule(pattern1105, With1105)

    pattern1106 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons682)
    rule1106 = ReplacementRule(pattern1106, replacement1106)

    pattern1107 = Pattern(Integral((a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons587, cons40, cons683)
    rule1107 = ReplacementRule(pattern1107, replacement1107)

    pattern1108 = Pattern(Integral((a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons587, cons149, cons683)
    rule1108 = ReplacementRule(pattern1108, replacement1108)

    pattern1109 = Pattern(Integral((a_ + u_**n_*WC('b', S(1)) + u_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons70, cons71)
    rule1109 = ReplacementRule(pattern1109, replacement1109)

    pattern1110 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons682, cons55)
    rule1110 = ReplacementRule(pattern1110, replacement1110)

    pattern1111 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons682, cons130, cons503)
    rule1111 = ReplacementRule(pattern1111, replacement1111)

    pattern1112 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons682, cons65, cons504)
    rule1112 = ReplacementRule(pattern1112, replacement1112)

    pattern1113 = Pattern(Integral(sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))/x_, x_), cons2, cons3, cons8, cons4, cons682, cons47)
    rule1113 = ReplacementRule(pattern1113, replacement1113)

    pattern1114 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons8, cons4, cons682, cons47, cons13, cons148)
    rule1114 = ReplacementRule(pattern1114, replacement1114)

    pattern1115 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons8, cons4, cons682, cons47, cons13, cons139)
    rule1115 = ReplacementRule(pattern1115, replacement1115)

    pattern1116 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_/x_, x_), cons2, cons3, cons8, cons4, cons5, cons682, cons47, cons149)
    rule1116 = ReplacementRule(pattern1116, replacement1116)

    pattern1117 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682, cons47, cons684)
    rule1117 = ReplacementRule(pattern1117, replacement1117)

    pattern1118 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons682, cons47, cons157)
    rule1118 = ReplacementRule(pattern1118, replacement1118)

    pattern1119 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons682, cons47, cons685)
    rule1119 = ReplacementRule(pattern1119, replacement1119)

    pattern1120 = Pattern(Integral(x_**WC('m', S(1))/sqrt(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons19, cons4, cons682, cons47, cons157)
    rule1120 = ReplacementRule(pattern1120, replacement1120)

    pattern1121 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682, cons47, cons686, cons227)
    rule1121 = ReplacementRule(pattern1121, replacement1121)

    pattern1122 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons682, cons47, cons58, cons245)
    rule1122 = ReplacementRule(pattern1122, replacement1122)

    pattern1123 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons47, cons150, cons246, cons148, cons687, cons248, cons20)
    rule1123 = ReplacementRule(pattern1123, replacement1123)

    pattern1124 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons47, cons150, cons246, cons148, cons688, cons689, cons248, cons20)
    rule1124 = ReplacementRule(pattern1124, replacement1124)

    pattern1125 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons682, cons47, cons150, cons13, cons148, cons512, cons690, cons689, cons691, cons248)
    rule1125 = ReplacementRule(pattern1125, replacement1125)

    pattern1126 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons47, cons150, cons246, cons139, cons692, cons248)
    rule1126 = ReplacementRule(pattern1126, replacement1126)

    pattern1127 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons47, cons150, cons246, cons139, cons531, cons248)
    rule1127 = ReplacementRule(pattern1127, replacement1127)

    pattern1128 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons682, cons47, cons150, cons246, cons139, cons248)
    rule1128 = ReplacementRule(pattern1128, replacement1128)

    pattern1129 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons47, cons150, cons33, cons532, cons512, cons693)
    rule1129 = ReplacementRule(pattern1129, replacement1129)

    pattern1130 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons47, cons150, cons33, cons96, cons693)
    rule1130 = ReplacementRule(pattern1130, replacement1130)

    pattern1131 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons682, cons47, cons198, cons20)
    rule1131 = ReplacementRule(pattern1131, replacement1131)

    pattern1132 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons47, cons198, cons369)
    rule1132 = ReplacementRule(pattern1132, With1132)

    pattern1133 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons682, cons47, cons198, cons358)
    rule1133 = ReplacementRule(pattern1133, replacement1133)

    pattern1134 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682, cons47, cons149)
    rule1134 = ReplacementRule(pattern1134, replacement1134)

    pattern1135 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons682, cons228, cons502)
    rule1135 = ReplacementRule(pattern1135, replacement1135)

    pattern1136 = Pattern(Integral((d_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682, cons228, cons502)
    rule1136 = ReplacementRule(pattern1136, replacement1136)

    pattern1137 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons682, cons228, cons150, cons20, CustomConstraint(With1137))
    rule1137 = ReplacementRule(pattern1137, replacement1137)

    pattern1138 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons228, cons150, cons369, cons40)
    rule1138 = ReplacementRule(pattern1138, With1138)

    pattern1139 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons246, cons165, cons532, cons694, cons695, cons696)
    rule1139 = ReplacementRule(pattern1139, replacement1139)

    pattern1140 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons246, cons165, cons96, cons696)
    rule1140 = ReplacementRule(pattern1140, replacement1140)

    pattern1141 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons682, cons228, cons150, cons13, cons165, cons512, cons696)
    rule1141 = ReplacementRule(pattern1141, replacement1141)

    pattern1142 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons246, cons139, cons692, cons696)
    rule1142 = ReplacementRule(pattern1142, replacement1142)

    pattern1143 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons246, cons139, cons531, cons696)
    rule1143 = ReplacementRule(pattern1143, replacement1143)

    pattern1144 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons682, cons228, cons150, cons13, cons139, cons696)
    rule1144 = ReplacementRule(pattern1144, replacement1144)

    pattern1145 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons228, cons150, cons33, cons531, cons512, cons696)
    rule1145 = ReplacementRule(pattern1145, replacement1145)

    pattern1146 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons228, cons150, cons33, cons96, cons696)
    rule1146 = ReplacementRule(pattern1146, replacement1146)

    pattern1147 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons33, cons96)
    rule1147 = ReplacementRule(pattern1147, replacement1147)

    pattern1148 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons682, cons228, cons150, cons20, cons697)
    rule1148 = ReplacementRule(pattern1148, replacement1148)

    pattern1149 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons33, cons531)
    rule1149 = ReplacementRule(pattern1149, replacement1149)

    pattern1150 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons698, cons699)
    rule1150 = ReplacementRule(pattern1150, With1150)

    pattern1151 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons48, cons228, cons700, cons701, cons415)
    rule1151 = ReplacementRule(pattern1151, With1151)

    pattern1152 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons48, cons228, cons700, cons702, cons415)
    rule1152 = ReplacementRule(pattern1152, With1152)

    pattern1153 = Pattern(Integral((x_*WC('d', S(1)))**m_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons682, cons228, cons150, cons33, cons703)
    rule1153 = ReplacementRule(pattern1153, With1153)

    pattern1154 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons682, cons228, cons150)
    rule1154 = ReplacementRule(pattern1154, With1154)

    pattern1155 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons295)
    rule1155 = ReplacementRule(pattern1155, With1155)

    pattern1156 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons678, cons679)
    rule1156 = ReplacementRule(pattern1156, With1156)

    pattern1157 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, cons486, cons179)
    rule1157 = ReplacementRule(pattern1157, With1157)

    pattern1158 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1158))
    rule1158 = ReplacementRule(pattern1158, replacement1158)

    pattern1159 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1159))
    rule1159 = ReplacementRule(pattern1159, replacement1159)

    pattern1160 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1160))
    rule1160 = ReplacementRule(pattern1160, replacement1160)

    pattern1161 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons677, CustomConstraint(With1161))
    rule1161 = ReplacementRule(pattern1161, replacement1161)

    pattern1162 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons680)
    rule1162 = ReplacementRule(pattern1162, With1162)

    pattern1163 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons228, cons681)
    rule1163 = ReplacementRule(pattern1163, With1163)

    pattern1164 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons5, cons682, cons228, cons198, cons20)
    rule1164 = ReplacementRule(pattern1164, replacement1164)

    pattern1165 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons682, cons228, cons198, cons369)
    rule1165 = ReplacementRule(pattern1165, With1165)

    pattern1166 = Pattern(Integral((x_*WC('d', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons682, cons228, cons198, cons358)
    rule1166 = ReplacementRule(pattern1166, replacement1166)

    pattern1167 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons5, cons682, cons228, cons491)
    rule1167 = ReplacementRule(pattern1167, With1167)

    pattern1168 = Pattern(Integral((d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons682, cons228, cons491)
    rule1168 = ReplacementRule(pattern1168, replacement1168)

    pattern1169 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons682, cons228, cons543, cons25)
    rule1169 = ReplacementRule(pattern1169, replacement1169)

    pattern1170 = Pattern(Integral((d_*x_)**m_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682, cons228, cons543, cons25)
    rule1170 = ReplacementRule(pattern1170, replacement1170)

    pattern1171 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons682, cons228)
    rule1171 = ReplacementRule(pattern1171, With1171)

    pattern1172 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons682, cons228, cons704)
    rule1172 = ReplacementRule(pattern1172, replacement1172)

    pattern1173 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons682)
    rule1173 = ReplacementRule(pattern1173, replacement1173)

    pattern1174 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons587, cons40, cons683)
    rule1174 = ReplacementRule(pattern1174, replacement1174)

    pattern1175 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons587, cons149, cons683)
    rule1175 = ReplacementRule(pattern1175, replacement1175)

    pattern1176 = Pattern(Integral((d_*x_)**WC('m', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons587)
    rule1176 = ReplacementRule(pattern1176, replacement1176)

    pattern1177 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons682, cons554, cons20, cons555)
    rule1177 = ReplacementRule(pattern1177, replacement1177)

    pattern1178 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons682, cons556)
    rule1178 = ReplacementRule(pattern1178, replacement1178)

    pattern1179 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons222, cons504)
    rule1179 = ReplacementRule(pattern1179, replacement1179)

    pattern1180 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons222, cons504)
    rule1180 = ReplacementRule(pattern1180, replacement1180)

    pattern1181 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons48, cons198)
    rule1181 = ReplacementRule(pattern1181, replacement1181)

    pattern1182 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons52, cons48, cons198)
    rule1182 = ReplacementRule(pattern1182, replacement1182)

    pattern1183 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons48, cons491)
    rule1183 = ReplacementRule(pattern1183, With1183)

    pattern1184 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons52, cons48, cons491)
    rule1184 = ReplacementRule(pattern1184, With1184)

    pattern1185 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons149, cons666)
    rule1185 = ReplacementRule(pattern1185, replacement1185)

    pattern1186 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons149, cons673, cons705)
    rule1186 = ReplacementRule(pattern1186, replacement1186)

    pattern1187 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons149, cons673, cons706)
    rule1187 = ReplacementRule(pattern1187, replacement1187)

    pattern1188 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons149)
    rule1188 = ReplacementRule(pattern1188, replacement1188)

    pattern1189 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons47, cons149)
    rule1189 = ReplacementRule(pattern1189, replacement1189)

    pattern1190 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons48, cons228, cons258, cons40)
    rule1190 = ReplacementRule(pattern1190, replacement1190)

    pattern1191 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons4, cons52, cons48, cons259, cons40)
    rule1191 = ReplacementRule(pattern1191, replacement1191)

    pattern1192 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons228, cons258, cons149)
    rule1192 = ReplacementRule(pattern1192, replacement1192)

    pattern1193 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons259, cons149)
    rule1193 = ReplacementRule(pattern1193, replacement1193)

    pattern1194 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons281, cons223)
    rule1194 = ReplacementRule(pattern1194, replacement1194)

    pattern1195 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons223)
    rule1195 = ReplacementRule(pattern1195, replacement1195)

    pattern1196 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons281, cons397, cons398)
    rule1196 = ReplacementRule(pattern1196, replacement1196)

    pattern1197 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons397, cons398)
    rule1197 = ReplacementRule(pattern1197, replacement1197)

    pattern1198 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons48, cons228, cons281)
    rule1198 = ReplacementRule(pattern1198, replacement1198)

    pattern1199 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons52, cons48, cons282)
    rule1199 = ReplacementRule(pattern1199, replacement1199)

    pattern1200 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons48, cons707, cons676, cons708)
    rule1200 = ReplacementRule(pattern1200, With1200)

    pattern1201 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons48, cons707, cons676, cons709)
    rule1201 = ReplacementRule(pattern1201, With1201)

    pattern1202 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons710, cons699)
    rule1202 = ReplacementRule(pattern1202, With1202)

    pattern1203 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons48, cons282, cons710, cons676, cons699)
    rule1203 = ReplacementRule(pattern1203, With1203)

    pattern1204 = Pattern(Integral((d_ + x_**S(3)*WC('e', S(1)))/(a_ + x_**S(6)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons680)
    rule1204 = ReplacementRule(pattern1204, With1204)

    pattern1205 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons711, cons87)
    rule1205 = ReplacementRule(pattern1205, With1205)

    pattern1206 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons712)
    rule1206 = ReplacementRule(pattern1206, replacement1206)

    pattern1207 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons707, cons676, cons713)
    rule1207 = ReplacementRule(pattern1207, With1207)

    pattern1208 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons707, cons676, cons677)
    rule1208 = ReplacementRule(pattern1208, With1208)

    pattern1209 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons707, cons676, cons714)
    rule1209 = ReplacementRule(pattern1209, With1209)

    pattern1210 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons281, cons715)
    rule1210 = ReplacementRule(pattern1210, With1210)

    pattern1211 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons281, cons676, cons415)
    rule1211 = ReplacementRule(pattern1211, With1211)

    pattern1212 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons281, cons588)
    rule1212 = ReplacementRule(pattern1212, replacement1212)

    pattern1213 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons588)
    rule1213 = ReplacementRule(pattern1213, replacement1213)

    pattern1214 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons281, cons388, cons397, cons398)
    rule1214 = ReplacementRule(pattern1214, replacement1214)

    pattern1215 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons282, cons388, cons397, cons398)
    rule1215 = ReplacementRule(pattern1215, replacement1215)

    pattern1216 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons48, cons228, cons281, cons388)
    rule1216 = ReplacementRule(pattern1216, With1216)

    pattern1217 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons52, cons48, cons282, cons388)
    rule1217 = ReplacementRule(pattern1217, With1217)

    pattern1218 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons151, cons165, cons671, cons716, cons248, cons675)
    rule1218 = ReplacementRule(pattern1218, replacement1218)

    pattern1219 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons151, cons165, cons671, cons716, cons248, cons675)
    rule1219 = ReplacementRule(pattern1219, replacement1219)

    pattern1220 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228, cons13, cons139, cons248, cons675)
    rule1220 = ReplacementRule(pattern1220, replacement1220)

    pattern1221 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48, cons13, cons139, cons248, cons675)
    rule1221 = ReplacementRule(pattern1221, replacement1221)

    pattern1222 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons295)
    rule1222 = ReplacementRule(pattern1222, With1222)

    pattern1223 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons45, cons295)
    rule1223 = ReplacementRule(pattern1223, With1223)

    pattern1224 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons678, cons679, CustomConstraint(With1224))
    rule1224 = ReplacementRule(pattern1224, replacement1224)

    pattern1225 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons678, cons679, CustomConstraint(With1225))
    rule1225 = ReplacementRule(pattern1225, replacement1225)

    pattern1226 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons486, cons179, CustomConstraint(With1226))
    rule1226 = ReplacementRule(pattern1226, replacement1226)

    pattern1227 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons486, cons179, CustomConstraint(With1227))
    rule1227 = ReplacementRule(pattern1227, replacement1227)

    pattern1228 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons486, cons179, CustomConstraint(With1228))
    rule1228 = ReplacementRule(pattern1228, replacement1228)

    pattern1229 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons486, cons179, CustomConstraint(With1229))
    rule1229 = ReplacementRule(pattern1229, replacement1229)

    pattern1230 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons486, cons179, CustomConstraint(With1230))
    rule1230 = ReplacementRule(pattern1230, replacement1230)

    pattern1231 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, CustomConstraint(With1231))
    rule1231 = ReplacementRule(pattern1231, replacement1231)

    pattern1232 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons717)
    rule1232 = ReplacementRule(pattern1232, replacement1232)

    pattern1233 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, CustomConstraint(With1233))
    rule1233 = ReplacementRule(pattern1233, replacement1233)

    pattern1234 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, CustomConstraint(With1234))
    rule1234 = ReplacementRule(pattern1234, replacement1234)

    pattern1235 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, CustomConstraint(With1235))
    rule1235 = ReplacementRule(pattern1235, replacement1235)

    pattern1236 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons677, CustomConstraint(With1236))
    rule1236 = ReplacementRule(pattern1236, replacement1236)

    pattern1237 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons680, CustomConstraint(With1237))
    rule1237 = ReplacementRule(pattern1237, replacement1237)

    pattern1238 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons680, CustomConstraint(With1238))
    rule1238 = ReplacementRule(pattern1238, replacement1238)

    pattern1239 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons680, CustomConstraint(With1239))
    rule1239 = ReplacementRule(pattern1239, replacement1239)

    pattern1240 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons680, CustomConstraint(With1240))
    rule1240 = ReplacementRule(pattern1240, replacement1240)

    pattern1241 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons681, cons259, cons45)
    rule1241 = ReplacementRule(pattern1241, replacement1241)

    pattern1242 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons681, cons259, cons450)
    rule1242 = ReplacementRule(pattern1242, replacement1242)

    pattern1243 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons681, cons282)
    rule1243 = ReplacementRule(pattern1243, With1243)

    pattern1244 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))/sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons681)
    rule1244 = ReplacementRule(pattern1244, With1244)

    pattern1245 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule1245 = ReplacementRule(pattern1245, replacement1245)

    pattern1246 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons4, cons48)
    rule1246 = ReplacementRule(pattern1246, replacement1246)

    pattern1247 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons48, cons228, cons130, cons718, cons150, cons402)
    rule1247 = ReplacementRule(pattern1247, replacement1247)

    pattern1248 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons52, cons48, cons130, cons718, cons150, cons402)
    rule1248 = ReplacementRule(pattern1248, replacement1248)

    pattern1249 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281)
    rule1249 = ReplacementRule(pattern1249, replacement1249)

    pattern1250 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons282)
    rule1250 = ReplacementRule(pattern1250, replacement1250)

    pattern1251 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**(S(3)/2)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons677)
    rule1251 = ReplacementRule(pattern1251, With1251)

    pattern1252 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**(S(3)/2)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons717)
    rule1252 = ReplacementRule(pattern1252, With1252)

    pattern1253 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons719)
    rule1253 = ReplacementRule(pattern1253, replacement1253)

    pattern1254 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons719)
    rule1254 = ReplacementRule(pattern1254, replacement1254)

    pattern1255 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons295)
    rule1255 = ReplacementRule(pattern1255, With1255)

    pattern1256 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons45, cons295)
    rule1256 = ReplacementRule(pattern1256, With1256)

    pattern1257 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons677, cons720)
    rule1257 = ReplacementRule(pattern1257, With1257)

    pattern1258 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons717, cons720)
    rule1258 = ReplacementRule(pattern1258, With1258)

    pattern1259 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons680, CustomConstraint(With1259))
    rule1259 = ReplacementRule(pattern1259, replacement1259)

    pattern1260 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons282, cons680, CustomConstraint(With1260))
    rule1260 = ReplacementRule(pattern1260, replacement1260)

    pattern1261 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons681, cons45)
    rule1261 = ReplacementRule(pattern1261, With1261)

    pattern1262 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons681, cons450)
    rule1262 = ReplacementRule(pattern1262, replacement1262)

    pattern1263 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons681)
    rule1263 = ReplacementRule(pattern1263, With1263)

    pattern1264 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281, cons721)
    rule1264 = ReplacementRule(pattern1264, replacement1264)

    pattern1265 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons282, cons721)
    rule1265 = ReplacementRule(pattern1265, replacement1265)

    pattern1266 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))**S(2)*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons228, cons281)
    rule1266 = ReplacementRule(pattern1266, replacement1266)

    pattern1267 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**S(2)), x_), cons2, cons8, cons29, cons50, cons282)
    rule1267 = ReplacementRule(pattern1267, replacement1267)

    pattern1268 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**q_*(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons52, cons228, cons281, cons349, cons566)
    rule1268 = ReplacementRule(pattern1268, With1268)

    pattern1269 = Pattern(Integral((a_ + x_**S(4)*WC('c', S(1)))**p_*(d_ + x_**S(2)*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons52, cons282, cons349, cons566)
    rule1269 = ReplacementRule(pattern1269, With1269)

    pattern1270 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons722, cons45, cons270)
    rule1270 = ReplacementRule(pattern1270, replacement1270)

    pattern1271 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons722, cons723)
    rule1271 = ReplacementRule(pattern1271, replacement1271)

    pattern1272 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons722, cons45, cons270)
    rule1272 = ReplacementRule(pattern1272, replacement1272)

    pattern1273 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons722, cons723)
    rule1273 = ReplacementRule(pattern1273, replacement1273)

    pattern1274 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons228, cons281, cons724)
    rule1274 = ReplacementRule(pattern1274, replacement1274)

    pattern1275 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons282, cons724)
    rule1275 = ReplacementRule(pattern1275, replacement1275)

    pattern1276 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons5, cons48, cons282, cons566, cons725)
    rule1276 = ReplacementRule(pattern1276, replacement1276)

    pattern1277 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons726)
    rule1277 = ReplacementRule(pattern1277, replacement1277)

    pattern1278 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons726)
    rule1278 = ReplacementRule(pattern1278, replacement1278)

    pattern1279 = Pattern(Integral((d_ + u_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + u_**n2_*WC('c', S(1)) + u_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons70, cons71)
    rule1279 = ReplacementRule(pattern1279, replacement1279)

    pattern1280 = Pattern(Integral((a_ + u_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + u_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons4, cons5, cons52, cons48, cons70, cons71)
    rule1280 = ReplacementRule(pattern1280, replacement1280)

    pattern1281 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons682, cons587, cons588)
    rule1281 = ReplacementRule(pattern1281, replacement1281)

    pattern1282 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons728, cons5, cons727, cons588)
    rule1282 = ReplacementRule(pattern1282, replacement1282)

    pattern1283 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons682, cons587, cons388, cons40)
    rule1283 = ReplacementRule(pattern1283, replacement1283)

    pattern1284 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons728, cons52, cons727, cons388, cons40)
    rule1284 = ReplacementRule(pattern1284, replacement1284)

    pattern1285 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons682, cons587, cons388, cons149, cons683)
    rule1285 = ReplacementRule(pattern1285, replacement1285)

    pattern1286 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons728, cons5, cons52, cons727, cons388, cons149, cons729)
    rule1286 = ReplacementRule(pattern1286, replacement1286)

    pattern1287 = Pattern(Integral((d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons682, cons587, cons388, cons149, cons504)
    rule1287 = ReplacementRule(pattern1287, replacement1287)

    pattern1288 = Pattern(Integral((a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons728, cons52, cons727, cons388, cons149, cons730)
    rule1288 = ReplacementRule(pattern1288, replacement1288)

    pattern1289 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons52, cons587, cons40)
    rule1289 = ReplacementRule(pattern1289, replacement1289)

    pattern1290 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons52, cons587, cons149)
    rule1290 = ReplacementRule(pattern1290, replacement1290)

    pattern1291 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons52, cons54, cons48, cons47, cons149)
    rule1291 = ReplacementRule(pattern1291, replacement1291)

    pattern1292 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons52, cons54, cons48, cons228, cons258, cons40)
    rule1292 = ReplacementRule(pattern1292, replacement1292)

    pattern1293 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons52, cons54, cons48, cons259, cons40)
    rule1293 = ReplacementRule(pattern1293, replacement1293)

    pattern1294 = Pattern(Integral((d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons52, cons54, cons48, cons228, cons258, cons149)
    rule1294 = ReplacementRule(pattern1294, replacement1294)

    pattern1295 = Pattern(Integral((a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(f_ + x_**n_*WC('g', S(1)))**WC('r', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons52, cons54, cons48, cons259, cons149)
    rule1295 = ReplacementRule(pattern1295, replacement1295)

    pattern1296 = Pattern(Integral((x_**S(2)*WC('g', S(1)) + WC('f', S(0)))/((d_ + x_**S(2)*WC('e', S(1)))*sqrt(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons677, cons281, cons720, CustomConstraint(With1296))
    rule1296 = ReplacementRule(pattern1296, replacement1296)

    pattern1297 = Pattern(Integral((f_ + x_**S(2)*WC('g', S(1)))/(sqrt(a_ + x_**S(4)*WC('c', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons8, cons29, cons50, cons127, cons210, cons717, cons282, cons720, CustomConstraint(With1297))
    rule1297 = ReplacementRule(pattern1297, replacement1297)

    pattern1298 = Pattern(Integral((d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons5, cons52, cons48, cons595, cons731, cons732)
    rule1298 = ReplacementRule(pattern1298, replacement1298)

    pattern1299 = Pattern(Integral((d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons5, cons52, cons48, cons595, cons731)
    rule1299 = ReplacementRule(pattern1299, replacement1299)

    pattern1300 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons19, cons4, cons5, cons52, cons48, cons55)
    rule1300 = ReplacementRule(pattern1300, replacement1300)

    pattern1301 = Pattern(Integral((A_ + x_**WC('m', S(1))*WC('B', S(1)))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons36, cons37, cons19, cons4, cons5, cons52, cons48, cons55)
    rule1301 = ReplacementRule(pattern1301, replacement1301)

    pattern1302 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons735, cons502)
    rule1302 = ReplacementRule(pattern1302, replacement1302)

    pattern1303 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons735, cons502)
    rule1303 = ReplacementRule(pattern1303, replacement1303)

    pattern1304 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons735, cons503)
    rule1304 = ReplacementRule(pattern1304, replacement1304)

    pattern1305 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons735, cons503)
    rule1305 = ReplacementRule(pattern1305, replacement1305)

    pattern1306 = Pattern(Integral((f_*x_)**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons21)
    rule1306 = ReplacementRule(pattern1306, replacement1306)

    pattern1307 = Pattern(Integral((f_*x_)**WC('m', S(1))*(x_**n_*WC('e', S(1)))**q_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons21)
    rule1307 = ReplacementRule(pattern1307, replacement1307)

    pattern1308 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons55)
    rule1308 = ReplacementRule(pattern1308, replacement1308)

    pattern1309 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons55)
    rule1309 = ReplacementRule(pattern1309, replacement1309)

    pattern1310 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons222, cons504)
    rule1310 = ReplacementRule(pattern1310, replacement1310)

    pattern1311 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons4, cons48, cons222, cons504)
    rule1311 = ReplacementRule(pattern1311, replacement1311)

    pattern1312 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons48, cons47, cons149, cons736)
    rule1312 = ReplacementRule(pattern1312, replacement1312)

    pattern1313 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons47, cons149)
    rule1313 = ReplacementRule(pattern1313, replacement1313)

    pattern1314 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons502)
    rule1314 = ReplacementRule(pattern1314, replacement1314)

    pattern1315 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons502)
    rule1315 = ReplacementRule(pattern1315, replacement1315)

    pattern1316 = Pattern(Integral((f_*x_)**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons502)
    rule1316 = ReplacementRule(pattern1316, replacement1316)

    pattern1317 = Pattern(Integral((f_*x_)**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons502)
    rule1317 = ReplacementRule(pattern1317, replacement1317)

    pattern1318 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons52, cons48, cons228, cons258, cons40)
    rule1318 = ReplacementRule(pattern1318, replacement1318)

    pattern1319 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons52, cons19, cons4, cons52, cons48, cons259, cons40)
    rule1319 = ReplacementRule(pattern1319, replacement1319)

    pattern1320 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons228, cons258, cons149)
    rule1320 = ReplacementRule(pattern1320, replacement1320)

    pattern1321 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons259, cons149)
    rule1321 = ReplacementRule(pattern1321, replacement1321)

    pattern1322 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons466, cons737, cons398, cons170)
    rule1322 = ReplacementRule(pattern1322, replacement1322)

    pattern1323 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons48, cons466, cons737, cons398, cons170)
    rule1323 = ReplacementRule(pattern1323, replacement1323)

    pattern1324 = Pattern(Integral(x_**m_*(d_ + x_**n_*WC('e', S(1)))**q_*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons466, cons737, cons398, cons269)
    rule1324 = ReplacementRule(pattern1324, replacement1324)

    pattern1325 = Pattern(Integral(x_**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons48, cons466, cons737, cons398, cons269)
    rule1325 = ReplacementRule(pattern1325, replacement1325)

    pattern1326 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons228, cons466, cons738, cons388, cons739)
    rule1326 = ReplacementRule(pattern1326, replacement1326)

    pattern1327 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons466, cons738, cons388, cons739)
    rule1327 = ReplacementRule(pattern1327, replacement1327)

    pattern1328 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons466)
    rule1328 = ReplacementRule(pattern1328, replacement1328)

    pattern1329 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons48, cons466)
    rule1329 = ReplacementRule(pattern1329, replacement1329)

    pattern1330 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons48, cons228, cons150, cons20, CustomConstraint(With1330))
    rule1330 = ReplacementRule(pattern1330, replacement1330)

    pattern1331 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons52, cons48, cons150, cons20, CustomConstraint(With1331))
    rule1331 = ReplacementRule(pattern1331, replacement1331)

    pattern1332 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons48, cons228, cons150, cons369, cons40)
    rule1332 = ReplacementRule(pattern1332, With1332)

    pattern1333 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons5, cons52, cons48, cons150, cons369, cons40)
    rule1333 = ReplacementRule(pattern1333, With1333)

    pattern1334 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons165, cons96, cons740, cons696)
    rule1334 = ReplacementRule(pattern1334, replacement1334)

    pattern1335 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons165, cons96, cons740, cons696)
    rule1335 = ReplacementRule(pattern1335, replacement1335)

    pattern1336 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150, cons13, cons165, cons512, cons741, cons696)
    rule1336 = ReplacementRule(pattern1336, replacement1336)

    pattern1337 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150, cons13, cons165, cons512, cons741, cons696)
    rule1337 = ReplacementRule(pattern1337, replacement1337)

    pattern1338 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons139, cons532, cons696)
    rule1338 = ReplacementRule(pattern1338, replacement1338)

    pattern1339 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons139, cons532, cons696)
    rule1339 = ReplacementRule(pattern1339, replacement1339)

    pattern1340 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150, cons13, cons139, cons696)
    rule1340 = ReplacementRule(pattern1340, replacement1340)

    pattern1341 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150, cons13, cons139, cons696)
    rule1341 = ReplacementRule(pattern1341, replacement1341)

    pattern1342 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons48, cons228, cons150, cons33, cons532, cons741, cons696)
    rule1342 = ReplacementRule(pattern1342, replacement1342)

    pattern1343 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons5, cons48, cons150, cons33, cons532, cons741, cons696)
    rule1343 = ReplacementRule(pattern1343, replacement1343)

    pattern1344 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons48, cons228, cons150, cons33, cons96, cons696)
    rule1344 = ReplacementRule(pattern1344, replacement1344)

    pattern1345 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons5, cons48, cons150, cons33, cons96, cons696)
    rule1345 = ReplacementRule(pattern1345, replacement1345)

    pattern1346 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons698, cons742, cons83, cons699, CustomConstraint(With1346))
    rule1346 = ReplacementRule(pattern1346, replacement1346)

    pattern1347 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons436, cons742, cons83, CustomConstraint(With1347))
    rule1347 = ReplacementRule(pattern1347, replacement1347)

    pattern1348 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1)) + x_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons228, cons707, cons743, cons744)
    rule1348 = ReplacementRule(pattern1348, With1348)

    pattern1349 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))/(a_ + x_**S(4)*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons707, cons743)
    rule1349 = ReplacementRule(pattern1349, With1349)

    pattern1350 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons698, cons745, cons746, cons699, CustomConstraint(With1350))
    rule1350 = ReplacementRule(pattern1350, replacement1350)

    pattern1351 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons745, cons746, cons436, CustomConstraint(With1351))
    rule1351 = ReplacementRule(pattern1351, replacement1351)

    pattern1352 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150)
    rule1352 = ReplacementRule(pattern1352, With1352)

    pattern1353 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150)
    rule1353 = ReplacementRule(pattern1353, With1353)

    pattern1354 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150, cons588, cons20)
    rule1354 = ReplacementRule(pattern1354, replacement1354)

    pattern1355 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150, cons588, cons20)
    rule1355 = ReplacementRule(pattern1355, replacement1355)

    pattern1356 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150, cons588, cons21)
    rule1356 = ReplacementRule(pattern1356, replacement1356)

    pattern1357 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150, cons588, cons21)
    rule1357 = ReplacementRule(pattern1357, replacement1357)

    pattern1358 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons388, cons613, cons405, cons531)
    rule1358 = ReplacementRule(pattern1358, replacement1358)

    pattern1359 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons52, cons48, cons150, cons388, cons33, cons531)
    rule1359 = ReplacementRule(pattern1359, replacement1359)

    pattern1360 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons388, cons613, cons405, cons692)
    rule1360 = ReplacementRule(pattern1360, replacement1360)

    pattern1361 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons388, cons613, cons405, cons692)
    rule1361 = ReplacementRule(pattern1361, replacement1361)

    pattern1362 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons388, cons613, cons405, cons269)
    rule1362 = ReplacementRule(pattern1362, replacement1362)

    pattern1363 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons388, cons613, cons405, cons269)
    rule1363 = ReplacementRule(pattern1363, replacement1363)

    pattern1364 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons388, cons613, cons398, cons531)
    rule1364 = ReplacementRule(pattern1364, replacement1364)

    pattern1365 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons388, cons613, cons398, cons531)
    rule1365 = ReplacementRule(pattern1365, replacement1365)

    pattern1366 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons388, cons613, cons398, cons692)
    rule1366 = ReplacementRule(pattern1366, replacement1366)

    pattern1367 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('e', S(1)) + WC('d', S(0)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons388, cons613, cons398, cons692)
    rule1367 = ReplacementRule(pattern1367, replacement1367)

    pattern1368 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons48, cons228, cons150, cons388, cons397, cons398)
    rule1368 = ReplacementRule(pattern1368, replacement1368)

    pattern1369 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n2_*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons48, cons150, cons388, cons397, cons398)
    rule1369 = ReplacementRule(pattern1369, replacement1369)

    pattern1370 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons52, cons4, cons48, cons228, cons150, cons388, cons20)
    rule1370 = ReplacementRule(pattern1370, replacement1370)

    pattern1371 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons52, cons4, cons48, cons150, cons388, cons20)
    rule1371 = ReplacementRule(pattern1371, replacement1371)

    pattern1372 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons52, cons4, cons48, cons228, cons150, cons388, cons21)
    rule1372 = ReplacementRule(pattern1372, replacement1372)

    pattern1373 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons52, cons4, cons48, cons150, cons388, cons21)
    rule1373 = ReplacementRule(pattern1373, replacement1373)

    pattern1374 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons165, cons747)
    rule1374 = ReplacementRule(pattern1374, replacement1374)

    pattern1375 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons165, cons747)
    rule1375 = ReplacementRule(pattern1375, replacement1375)

    pattern1376 = Pattern(Integral((x_*WC('f', S(1)))**m_*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons165, cons269)
    rule1376 = ReplacementRule(pattern1376, replacement1376)

    pattern1377 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons165, cons269)
    rule1377 = ReplacementRule(pattern1377, replacement1377)

    pattern1378 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons139, cons748)
    rule1378 = ReplacementRule(pattern1378, replacement1378)

    pattern1379 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons139, cons748)
    rule1379 = ReplacementRule(pattern1379, replacement1379)

    pattern1380 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons48, cons228, cons150, cons246, cons139, cons170)
    rule1380 = ReplacementRule(pattern1380, replacement1380)

    pattern1381 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_/(x_**n_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons8, cons29, cons50, cons127, cons48, cons150, cons246, cons139, cons170)
    rule1381 = ReplacementRule(pattern1381, replacement1381)

    pattern1382 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons228, cons150, cons749)
    rule1382 = ReplacementRule(pattern1382, replacement1382)

    pattern1383 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons52, cons48, cons150, cons749)
    rule1383 = ReplacementRule(pattern1383, replacement1383)

    pattern1384 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons48, cons228, cons198, cons20)
    rule1384 = ReplacementRule(pattern1384, replacement1384)

    pattern1385 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons5, cons52, cons48, cons198, cons20)
    rule1385 = ReplacementRule(pattern1385, replacement1385)

    pattern1386 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons52, cons48, cons228, cons198, cons369)
    rule1386 = ReplacementRule(pattern1386, With1386)

    pattern1387 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons5, cons52, cons48, cons198, cons369)
    rule1387 = ReplacementRule(pattern1387, With1387)

    pattern1388 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons228, cons198, cons358)
    rule1388 = ReplacementRule(pattern1388, replacement1388)

    pattern1389 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons198, cons358)
    rule1389 = ReplacementRule(pattern1389, replacement1389)

    pattern1390 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons52, cons48, cons228, cons491)
    rule1390 = ReplacementRule(pattern1390, With1390)

    pattern1391 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons5, cons52, cons48, cons491)
    rule1391 = ReplacementRule(pattern1391, With1391)

    pattern1392 = Pattern(Integral((f_*x_)**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons228, cons491)
    rule1392 = ReplacementRule(pattern1392, replacement1392)

    pattern1393 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons491)
    rule1393 = ReplacementRule(pattern1393, replacement1393)

    pattern1394 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons228, cons543, cons25)
    rule1394 = ReplacementRule(pattern1394, replacement1394)

    pattern1395 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons543, cons25)
    rule1395 = ReplacementRule(pattern1395, replacement1395)

    pattern1396 = Pattern(Integral((f_*x_)**m_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons228, cons543, cons25)
    rule1396 = ReplacementRule(pattern1396, replacement1396)

    pattern1397 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons52, cons48, cons543, cons25)
    rule1397 = ReplacementRule(pattern1397, replacement1397)

    pattern1398 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons52, cons48, cons228)
    rule1398 = ReplacementRule(pattern1398, With1398)

    pattern1399 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**q_/(a_ + x_**WC('n2', S(1))*WC('c', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons52, cons48)
    rule1399 = ReplacementRule(pattern1399, With1399)

    pattern1400 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))*(a_ + x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons48, cons228, cons704)
    rule1400 = ReplacementRule(pattern1400, replacement1400)

    pattern1401 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons48, cons704)
    rule1401 = ReplacementRule(pattern1401, replacement1401)

    pattern1402 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons228, cons750)
    rule1402 = ReplacementRule(pattern1402, replacement1402)

    pattern1403 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons750)
    rule1403 = ReplacementRule(pattern1403, replacement1403)

    pattern1404 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons566, cons735)
    rule1404 = ReplacementRule(pattern1404, replacement1404)

    pattern1405 = Pattern(Integral((x_*WC('f', S(1)))**m_*(a_ + x_**n2_*WC('c', S(1)))**p_*(d_ + x_**n_*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48, cons566, cons751)
    rule1405 = ReplacementRule(pattern1405, replacement1405)

    pattern1406 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**n_*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48)
    rule1406 = ReplacementRule(pattern1406, replacement1406)

    pattern1407 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons48)
    rule1407 = ReplacementRule(pattern1407, replacement1407)

    pattern1408 = Pattern(Integral(u_**WC('m', S(1))*(d_ + v_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + v_**n_*WC('b', S(1)) + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons48, cons556)
    rule1408 = ReplacementRule(pattern1408, replacement1408)

    pattern1409 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + v_**n_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons556)
    rule1409 = ReplacementRule(pattern1409, replacement1409)

    pattern1410 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons682, cons587, cons588)
    rule1410 = ReplacementRule(pattern1410, replacement1410)

    pattern1411 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons19, cons728, cons5, cons727, cons588)
    rule1411 = ReplacementRule(pattern1411, replacement1411)

    pattern1412 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons52, cons682, cons587, cons388, cons40)
    rule1412 = ReplacementRule(pattern1412, replacement1412)

    pattern1413 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons19, cons728, cons52, cons727, cons388, cons40)
    rule1413 = ReplacementRule(pattern1413, replacement1413)

    pattern1414 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons682, cons587, cons388, cons149)
    rule1414 = ReplacementRule(pattern1414, replacement1414)

    pattern1415 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**q_, x_), cons2, cons8, cons29, cons50, cons19, cons728, cons5, cons52, cons727, cons388, cons149)
    rule1415 = ReplacementRule(pattern1415, replacement1415)

    pattern1416 = Pattern(Integral((f_*x_)**m_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons682, cons587)
    rule1416 = ReplacementRule(pattern1416, replacement1416)

    pattern1417 = Pattern(Integral((f_*x_)**m_*(a_ + x_**WC('n2', S(1))*WC('c', S(1)))**p_*(d_ + x_**WC('mn', S(1))*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons728, cons5, cons52, cons727)
    rule1417 = ReplacementRule(pattern1417, replacement1417)

    pattern1418 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons52, cons587, cons40)
    rule1418 = ReplacementRule(pattern1418, replacement1418)

    pattern1419 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons52, cons587, cons149)
    rule1419 = ReplacementRule(pattern1419, replacement1419)

    pattern1420 = Pattern(Integral((f_*x_)**WC('m', S(1))*(d_ + x_**n_*WC('e', S(1)))**WC('q', S(1))*(a_ + x_**mn_*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons587)
    rule1420 = ReplacementRule(pattern1420, replacement1420)

    pattern1421 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons4, cons5, cons52, cons48, cons595, cons731, cons732)
    rule1421 = ReplacementRule(pattern1421, replacement1421)

    pattern1422 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_**WC('non2', S(1))*WC('e1', S(1)))**WC('q', S(1))*(d2_ + x_**WC('non2', S(1))*WC('e2', S(1)))**WC('q', S(1))*(x_**n2_*WC('c', S(1)) + x_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons4, cons5, cons52, cons48, cons595, cons731)
    rule1422 = ReplacementRule(pattern1422, replacement1422)

    pattern1423 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons752, cons753)
    rule1423 = ReplacementRule(pattern1423, replacement1423)

    pattern1424 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons52, cons754, cons755, cons40)
    rule1424 = ReplacementRule(pattern1424, replacement1424)

    pattern1425 = Pattern(Integral(sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons4, cons52, cons754, cons755)
    rule1425 = ReplacementRule(pattern1425, replacement1425)

    pattern1426 = Pattern(Integral(S(1)/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons4, cons52, cons754, cons755)
    rule1426 = ReplacementRule(pattern1426, replacement1426)

    pattern1427 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons52, cons754, cons755, cons149, cons228, cons13, cons165, cons756)
    rule1427 = ReplacementRule(pattern1427, replacement1427)

    pattern1428 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons52, cons754, cons755, cons149, cons228, cons13, cons139)
    rule1428 = ReplacementRule(pattern1428, replacement1428)

    pattern1429 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons52, cons754, cons755, cons149)
    rule1429 = ReplacementRule(pattern1429, replacement1429)

    pattern1430 = Pattern(Integral((x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons52, cons754)
    rule1430 = ReplacementRule(pattern1430, replacement1430)

    pattern1431 = Pattern(Integral((u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons52, cons754, cons70, cons71)
    rule1431 = ReplacementRule(pattern1431, replacement1431)

    pattern1432 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons757, cons753)
    rule1432 = ReplacementRule(pattern1432, replacement1432)

    pattern1433 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons52, cons754, cons40, cons755)
    rule1433 = ReplacementRule(pattern1433, replacement1433)

    pattern1434 = Pattern(Integral(x_**WC('m', S(1))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons19, cons4, cons52, cons754, cons755, cons758)
    rule1434 = ReplacementRule(pattern1434, replacement1434)

    pattern1435 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons4, cons759, cons760, cons761, cons228)
    rule1435 = ReplacementRule(pattern1435, replacement1435)

    pattern1436 = Pattern(Integral(x_**WC('m', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons4, cons762, cons760, cons761, cons228)
    rule1436 = ReplacementRule(pattern1436, replacement1436)

    pattern1437 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons763)
    rule1437 = ReplacementRule(pattern1437, replacement1437)

    pattern1438 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons165, cons764)
    rule1438 = ReplacementRule(pattern1438, replacement1438)

    pattern1439 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons165, cons765, cons766, cons767)
    rule1439 = ReplacementRule(pattern1439, replacement1439)

    pattern1440 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons165, cons768, cons769)
    rule1440 = ReplacementRule(pattern1440, replacement1440)

    pattern1441 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons165, cons770, cons766)
    rule1441 = ReplacementRule(pattern1441, replacement1441)

    pattern1442 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons139, cons771)
    rule1442 = ReplacementRule(pattern1442, replacement1442)

    pattern1443 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons139, cons772)
    rule1443 = ReplacementRule(pattern1443, replacement1443)

    pattern1444 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons139, cons773)
    rule1444 = ReplacementRule(pattern1444, replacement1444)

    pattern1445 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons139, cons774)
    rule1445 = ReplacementRule(pattern1445, replacement1445)

    pattern1446 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons775, cons776)
    rule1446 = ReplacementRule(pattern1446, replacement1446)

    pattern1447 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons775, cons777)
    rule1447 = ReplacementRule(pattern1447, replacement1447)

    pattern1448 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons775, cons772)
    rule1448 = ReplacementRule(pattern1448, replacement1448)

    pattern1449 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons754, cons755, cons149, cons228, cons150, cons608, cons775, cons778)
    rule1449 = ReplacementRule(pattern1449, replacement1449)

    pattern1450 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons52, cons754, cons149, cons755)
    rule1450 = ReplacementRule(pattern1450, replacement1450)

    pattern1451 = Pattern(Integral(u_**WC('m', S(1))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons52, cons754, cons70, cons71)
    rule1451 = ReplacementRule(pattern1451, replacement1451)

    pattern1452 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons52, cons779, cons780, cons40, cons755)
    rule1452 = ReplacementRule(pattern1452, replacement1452)

    pattern1453 = Pattern(Integral((A_ + x_**WC('j', S(1))*WC('B', S(1)))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons52, cons781, cons754, cons755, cons782, cons783)
    rule1453 = ReplacementRule(pattern1453, replacement1453)

    pattern1454 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons36, cons37, cons4, cons52, cons779, cons780, cons149, cons228, cons13, cons165, cons756, cons784)
    rule1454 = ReplacementRule(pattern1454, replacement1454)

    pattern1455 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons8, cons36, cons37, cons52, cons149, cons13, cons165, CustomConstraint(With1455))
    rule1455 = ReplacementRule(pattern1455, replacement1455)

    pattern1456 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons3, cons8, cons36, cons37, cons4, cons52, cons779, cons780, cons149, cons228, cons13, cons139)
    rule1456 = ReplacementRule(pattern1456, replacement1456)

    pattern1457 = Pattern(Integral((A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**p_, x_), cons2, cons8, cons36, cons37, cons52, cons149, cons13, cons139, CustomConstraint(With1457))
    rule1457 = ReplacementRule(pattern1457, replacement1457)

    pattern1458 = Pattern(Integral((A_ + x_**WC('j', S(1))*WC('B', S(1)))*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons5, cons52, cons781, cons754)
    rule1458 = ReplacementRule(pattern1458, replacement1458)

    pattern1459 = Pattern(Integral((A_ + u_**WC('j', S(1))*WC('B', S(1)))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons5, cons52, cons781, cons754, cons70, cons71)
    rule1459 = ReplacementRule(pattern1459, replacement1459)

    pattern1460 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons19, cons4, cons52, cons779, cons780, cons40, cons755)
    rule1460 = ReplacementRule(pattern1460, replacement1460)

    pattern1461 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons165, cons785, cons786, cons787)
    rule1461 = ReplacementRule(pattern1461, replacement1461)

    pattern1462 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, cons165, CustomConstraint(With1462))
    rule1462 = ReplacementRule(pattern1462, replacement1462)

    pattern1463 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons139, cons788)
    rule1463 = ReplacementRule(pattern1463, replacement1463)

    pattern1464 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, cons139, CustomConstraint(With1464))
    rule1464 = ReplacementRule(pattern1464, replacement1464)

    pattern1465 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons165, cons789, cons766, cons787)
    rule1465 = ReplacementRule(pattern1465, replacement1465)

    pattern1466 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, cons165, CustomConstraint(With1466))
    rule1466 = ReplacementRule(pattern1466, replacement1466)

    pattern1467 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons139, cons790)
    rule1467 = ReplacementRule(pattern1467, replacement1467)

    pattern1468 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, cons139, CustomConstraint(With1468))
    rule1468 = ReplacementRule(pattern1468, replacement1468)

    pattern1469 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons775, cons791, cons787)
    rule1469 = ReplacementRule(pattern1469, replacement1469)

    pattern1470 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, cons775, CustomConstraint(With1470))
    rule1470 = ReplacementRule(pattern1470, replacement1470)

    pattern1471 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons779, cons780, cons149, cons228, cons150, cons608, cons792, cons785, cons786)
    rule1471 = ReplacementRule(pattern1471, replacement1471)

    pattern1472 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('r', S(1))*WC('B', S(1)))*(x_**WC('j', S(1))*WC('c', S(1)) + x_**WC('q', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons36, cons37, cons149, cons608, CustomConstraint(With1472))
    rule1472 = ReplacementRule(pattern1472, replacement1472)

    pattern1473 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**WC('j', S(1))*WC('B', S(1)))/sqrt(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('q', S(1))*WC('a', S(1)) + x_**WC('r', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons36, cons37, cons19, cons4, cons52, cons781, cons754, cons755, cons793, cons782, cons794)
    rule1473 = ReplacementRule(pattern1473, replacement1473)

    pattern1474 = Pattern(Integral(x_**WC('m', S(1))*(A_ + x_**q_*WC('B', S(1)))*(x_**WC('j', S(1))*WC('a', S(1)) + x_**WC('k', S(1))*WC('b', S(1)) + x_**WC('n', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons36, cons37, cons798, cons799, cons19, cons5, cons795, cons796, cons149, cons797)
    rule1474 = ReplacementRule(pattern1474, replacement1474)

    pattern1475 = Pattern(Integral(u_**WC('m', S(1))*(A_ + u_**WC('j', S(1))*WC('B', S(1)))*(u_**WC('n', S(1))*WC('b', S(1)) + u_**WC('q', S(1))*WC('a', S(1)) + u_**WC('r', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons19, cons4, cons5, cons52, cons781, cons754, cons70, cons71)
    rule1475 = ReplacementRule(pattern1475, replacement1475)
    return [rule1079, rule1080, rule1081, rule1082, rule1083, rule1084, rule1085, rule1086, rule1087, rule1088, rule1089, rule1090, rule1091, rule1092, rule1093, rule1094, rule1095, rule1096, rule1097, rule1098, rule1099, rule1100, rule1101, rule1102, rule1103, rule1104, rule1105, rule1106, rule1107, rule1108, rule1109, rule1110, rule1111, rule1112, rule1113, rule1114, rule1115, rule1116, rule1117, rule1118, rule1119, rule1120, rule1121, rule1122, rule1123, rule1124, rule1125, rule1126, rule1127, rule1128, rule1129, rule1130, rule1131, rule1132, rule1133, rule1134, rule1135, rule1136, rule1137, rule1138, rule1139, rule1140, rule1141, rule1142, rule1143, rule1144, rule1145, rule1146, rule1147, rule1148, rule1149, rule1150, rule1151, rule1152, rule1153, rule1154, rule1155, rule1156, rule1157, rule1158, rule1159, rule1160, rule1161, rule1162, rule1163, rule1164, rule1165, rule1166, rule1167, rule1168, rule1169, rule1170, rule1171, rule1172, rule1173, rule1174, rule1175, rule1176, rule1177, rule1178, rule1179, rule1180, rule1181, rule1182, rule1183, rule1184, rule1185, rule1186, rule1187, rule1188, rule1189, rule1190, rule1191, rule1192, rule1193, rule1194, rule1195, rule1196, rule1197, rule1198, rule1199, rule1200, rule1201, rule1202, rule1203, rule1204, rule1205, rule1206, rule1207, rule1208, rule1209, rule1210, rule1211, rule1212, rule1213, rule1214, rule1215, rule1216, rule1217, rule1218, rule1219, rule1220, rule1221, rule1222, rule1223, rule1224, rule1225, rule1226, rule1227, rule1228, rule1229, rule1230, rule1231, rule1232, rule1233, rule1234, rule1235, rule1236, rule1237, rule1238, rule1239, rule1240, rule1241, rule1242, rule1243, rule1244, rule1245, rule1246, rule1247, rule1248, rule1249, rule1250, rule1251, rule1252, rule1253, rule1254, rule1255, rule1256, rule1257, rule1258, rule1259, rule1260, rule1261, rule1262, rule1263, rule1264, rule1265, rule1266, rule1267, rule1268, rule1269, rule1270, rule1271, rule1272, rule1273, rule1274, rule1275, rule1276, rule1277, rule1278, rule1279, rule1280, rule1281, rule1282, rule1283, rule1284, rule1285, rule1286, rule1287, rule1288, rule1289, rule1290, rule1291, rule1292, rule1293, rule1294, rule1295, rule1296, rule1297, rule1298, rule1299, rule1300, rule1301, rule1302, rule1303, rule1304, rule1305, rule1306, rule1307, rule1308, rule1309, rule1310, rule1311, rule1312, rule1313, rule1314, rule1315, rule1316, rule1317, rule1318, rule1319, rule1320, rule1321, rule1322, rule1323, rule1324, rule1325, rule1326, rule1327, rule1328, rule1329, rule1330, rule1331, rule1332, rule1333, rule1334, rule1335, rule1336, rule1337, rule1338, rule1339, rule1340, rule1341, rule1342, rule1343, rule1344, rule1345, rule1346, rule1347, rule1348, rule1349, rule1350, rule1351, rule1352, rule1353, rule1354, rule1355, rule1356, rule1357, rule1358, rule1359, rule1360, rule1361, rule1362, rule1363, rule1364, rule1365, rule1366, rule1367, rule1368, rule1369, rule1370, rule1371, rule1372, rule1373, rule1374, rule1375, rule1376, rule1377, rule1378, rule1379, rule1380, rule1381, rule1382, rule1383, rule1384, rule1385, rule1386, rule1387, rule1388, rule1389, rule1390, rule1391, rule1392, rule1393, rule1394, rule1395, rule1396, rule1397, rule1398, rule1399, rule1400, rule1401, rule1402, rule1403, rule1404, rule1405, rule1406, rule1407, rule1408, rule1409, rule1410, rule1411, rule1412, rule1413, rule1414, rule1415, rule1416, rule1417, rule1418, rule1419, rule1420, rule1421, rule1422, rule1423, rule1424, rule1425, rule1426, rule1427, rule1428, rule1429, rule1430, rule1431, rule1432, rule1433, rule1434, rule1435, rule1436, rule1437, rule1438, rule1439, rule1440, rule1441, rule1442, rule1443, rule1444, rule1445, rule1446, rule1447, rule1448, rule1449, rule1450, rule1451, rule1452, rule1453, rule1454, rule1455, rule1456, rule1457, rule1458, rule1459, rule1460, rule1461, rule1462, rule1463, rule1464, rule1465, rule1466, rule1467, rule1468, rule1469, rule1470, rule1471, rule1472, rule1473, rule1474, rule1475, ]





def replacement1079(a, b, c, n, n2, p, x):
    return Int(x**(S(2)*n*p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def With1080(a, b, c, n, n2, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*x**(k*n) + c*x**(S(2)*k*n))**p, x), x, x**(S(1)/k)), x)


def replacement1081(a, b, c, n, n2, p, x):
    return Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a), x)


def replacement1082(a, b, c, n, n2, p, x):
    return -Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*(S(2)*p + S(1))), x) + Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*(n + S(1))), x)


def replacement1083(a, b, c, n, n2, p, x):
    return Dist(sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int((b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)/2), x), x)


def replacement1084(a, b, c, n, n2, p, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p), x), x)


def replacement1085(a, b, c, n, n2, x):
    return Simp(x*sqrt(a + b*x**n + c*x**(S(2)*n))/(n + S(1)), x) + Simp(b*n*x*sqrt(a + b*x**n + c*x**(S(2)*n))/((b + S(2)*c*x**n)*(n + S(1))), x)


def replacement1086(a, b, c, n, n2, p, x):
    return -Subst(Int((a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)


def replacement1087(a, b, c, n, n2, p, x):
    return Dist(S(2)*a*n**S(2)*p*(S(2)*p + S(-1))/((S(2)*n*p + S(1))*(n*(S(2)*p + S(-1)) + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p + S(1)), x) + Simp(n*p*x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/((S(2)*n*p + S(1))*(n*(S(2)*p + S(-1)) + S(1))), x)


def replacement1088(a, b, c, n, n2, p, x):
    return Dist((S(2)*n*(p + S(1)) + S(1))*(n*(S(2)*p + S(1)) + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(x*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*n*(S(2)*p + S(1))), x) - Simp(x*(n*(S(2)*p + S(1)) + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)


def replacement1089(a, b, c, n, n2, p, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b/S(2) + c*x**n)**(S(2)*p), x), x)


def replacement1090(a, b, c, n, n2, p, x):
    return -Subst(Int((a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)


def replacement1091(a, b, c, n, n2, p, x):
    return Int(ExpandIntegrand((a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1092(a, b, c, n, n2, p, x):
    return Dist(n*p/(S(2)*n*p + S(1)), Int((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p + S(1)), x)


def replacement1093(a, b, c, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n*(n*(S(2)*p + S(3)) + S(1)) + n*(p + S(1))*(-S(4)*a*c + b**S(2))), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def With1094(a, b, c, n, n2, x):
    q = Rt(a/c, S(2))
    r = Rt(-b/c + S(2)*q, S(2))
    return Dist(1/(2*c*q*r), Int((r - x**(n/2))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((r + x**(n/2))/(q + r*x**(n/2) + x**n), x), x)


def With1095(a, b, c, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(c/q, Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x) - Dist(c/q, Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x)


def With1096(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*sqrt(-c), Int(S(1)/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)


def With1097(a, b, c, x):
    q = Rt(c/a, S(4))
    return Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(2)*q*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1098(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if IntegerQ(q):
        return True
    return False


def replacement1098(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt(-S(2)*a - x**S(2)*(b - q))*EllipticF(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*sqrt(-a)*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1099(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*EllipticF(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*sqrt(a/(S(2)*a + x**S(2)*(b + q)))*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1100(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(PosQ((b + q)/a), Not(And(PosQ((b - q)/a), SimplerSqrtQ((b - q)/(S(2)*a), (b + q)/(S(2)*a))))):
        return True
    return False


def replacement1100(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*(S(2)*a + x**S(2)*(b + q))*EllipticF(ArcTan(x*Rt((b + q)/(S(2)*a), S(2))), S(2)*q/(b + q))/(S(2)*a*sqrt(a + b*x**S(2) + c*x**S(4))*Rt((b + q)/(S(2)*a), S(2))), x)


def With1101(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if PosQ((b - q)/a):
        return True
    return False


def replacement1101(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt((S(2)*a + x**S(2)*(b + q))/(S(2)*a + x**S(2)*(b - q)))*(S(2)*a + x**S(2)*(b - q))*EllipticF(ArcTan(x*Rt((b - q)/(S(2)*a), S(2))), -S(2)*q/(b - q))/(S(2)*a*sqrt(a + b*x**S(2) + c*x**S(4))*Rt((b - q)/(S(2)*a), S(2))), x)


def With1102(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b + q)/a), Not(And(NegQ((b - q)/a), SimplerSqrtQ(-(b - q)/(S(2)*a), -(b + q)/(S(2)*a))))):
        return True
    return False


def replacement1102(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticF(asin(x*Rt(-(b + q)/(S(2)*a), S(2))), (b - q)/(b + q))/(sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-(b + q)/(S(2)*a), S(2))), x)


def With1103(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if NegQ((b - q)/a):
        return True
    return False


def replacement1103(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticF(asin(x*Rt(-(b - q)/(S(2)*a), S(2))), (b + q)/(b - q))/(sqrt(a + b*x**S(2) + c*x**S(4))*Rt(-(b - q)/(S(2)*a), S(2))), x)


def With1104(a, b, c, x):
    q = Rt(c/a, S(4))
    return Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(2)*q*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1105(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(S(1)/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)


def replacement1106(a, b, c, n, n2, p, x):
    return Dist(a**IntPart(p)*(S(2)*c*x**n/(b - Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(S(2)*c*x**n/(b + Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((S(2)*c*x**n/(b - sqrt(-S(4)*a*c + b**S(2))) + S(1))**p*(S(2)*c*x**n/(b + sqrt(-S(4)*a*c + b**S(2))) + S(1))**p, x), x)


def replacement1107(a, b, c, mn, n, p, x):
    return Int(x**(-n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x)


def replacement1108(a, b, c, mn, n, p, x):
    return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(-n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x), x)


def replacement1109(a, b, c, n, n2, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, u), x)


def replacement1110(a, b, c, m, n, n2, p, x):
    return Dist(S(1)/n, Subst(Int((a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1111(a, b, c, d, m, n, n2, p, x):
    return Int(ExpandIntegrand((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1112(a, b, c, m, n, n2, p, x):
    return Int(x**(m + S(2)*n*p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def replacement1113(a, b, c, n, n2, x):
    return Simp(sqrt(a + b*x**n + c*x**(S(2)*n))/n, x) + Simp(b*sqrt(a + b*x**n + c*x**(S(2)*n))*log(x)/(b + S(2)*c*x**n), x)


def replacement1114(a, b, c, n, n2, p, x):
    return Dist(a, Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/x, x), x) + Simp((a + b*x**n + c*x**(S(2)*n))**p/(S(2)*n*p), x) + Simp((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(S(2)*n*(S(2)*p + S(-1))), x)


def replacement1115(a, b, c, n, n2, p, x):
    return Dist(S(1)/a, Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/x, x), x) - Simp((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*n*(p + S(1))), x) - Simp((S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*n*(S(2)*p + S(1))), x)


def replacement1116(a, b, c, n, n2, p, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b/S(2) + c*x**n)**(S(2)*p)/x, x), x)


def replacement1117(a, b, c, d, m, n, n2, p, x):
    return Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*(m + S(1))), x)


def replacement1118(a, b, c, d, m, n, n2, x):
    return Dist(sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int((d*x)**m*(b + S(2)*c*x**n), x), x)


def replacement1119(a, b, c, d, m, n, n2, x):
    return Simp((d*x)**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(d*(m + n + S(1))), x) + Simp(b*n*(d*x)**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(d*(b + S(2)*c*x**n)*(m + S(1))*(m + n + S(1))), x)


def replacement1120(a, b, c, m, n, n2, x):
    return -Dist(b/(S(2)*a), Int(S(1)/(x*sqrt(a + b*x**n + c*x**(S(2)*n))), x), x) - Simp(x**(m + S(1))*sqrt(a + b*x**n + c*x**(S(2)*n))/(a*n), x)


def replacement1121(a, b, c, d, m, n, n2, p, x):
    return -Simp((d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*d*n*(S(2)*p + S(1))), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*a*d*n*(p + S(1))*(S(2)*p + S(1))), x)


def replacement1122(a, b, c, m, n, n2, p, x):
    return -Dist(b/(S(2)*c), Int(x**(n + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(S(2)*c*n*(p + S(1))), x)


def replacement1123(a, b, c, d, m, n, n2, p, x):
    return -Dist(b*d**(-n)*n**S(2)*p*(S(2)*p + S(-1))/((m + S(1))*(m + S(2)*n*p + S(1))), Int((d*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(1))*(m + S(2)*n*p + S(1))), x)


def replacement1124(a, b, c, d, m, n, n2, p, x):
    return Dist(S(2)*c*d**(-S(2)*n)*n**S(2)*p*(S(2)*p + S(-1))/((m + S(1))*(m + n + S(1))), Int((d*x)**(m + S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*(m - n*(S(2)*p + S(-1)) + S(1))/(d*(m + S(1))*(m + n + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(1))*(m + n + S(1))), x)


def replacement1125(a, b, c, d, m, n, n2, p, x):
    return Dist(S(2)*a*n**S(2)*p*(S(2)*p + S(-1))/((m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x) + Simp(n*p*(d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), x)


def replacement1126(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**n*(m - n + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(b*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*n*(S(2)*p + S(1))), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(b*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)


def replacement1127(a, b, c, d, m, n, n2, p, x):
    return Dist(d**(S(2)*n)*(m - S(2)*n + S(1))*(m - n + S(1))/(S(2)*c*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*c*n*(S(2)*p + S(1))), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m - S(2)*n*p - S(3)*n + S(1))/(S(2)*c*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)


def replacement1128(a, b, c, d, m, n, n2, p, x):
    return Dist((m + S(2)*n*(p + S(1)) + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(S(2)*a*n**S(2)*(p + S(1))*(S(2)*p + S(1))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp((d*x)**(m + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*a*d*n*(S(2)*p + S(1))), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))/(S(2)*a*d*n**S(2)*(p + S(1))*(S(2)*p + S(1))), x)


def replacement1129(a, b, c, d, m, n, n2, p, x):
    return -Dist(b*d**n*(m - n + S(1))/(S(2)*c*(m + S(2)*n*p + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(S(2)*c*(m + S(2)*n*p + S(1))), x)


def replacement1130(a, b, c, d, m, n, n2, p, x):
    return -Dist(S(2)*c*d**(-n)*(m + n*(S(2)*p + S(1)) + S(1))/(b*(m + S(1))), Int((d*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((d*x)**(m + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**p/(b*d*(m + S(1))), x)


def replacement1131(a, b, c, m, n, n2, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)


def With1132(a, b, c, d, m, n, n2, p, x):
    k = Denominator(m)
    return -Dist(k/d, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(-k*n) + c*d**(-S(2)*n)*x**(-S(2)*k*n))**p, x), x, (d*x)**(-S(1)/k)), x)


def replacement1133(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**IntPart(m)*(d*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)


def replacement1134(a, b, c, d, m, n, n2, p, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d*x)**m*(b/S(2) + c*x**n)**(S(2)*p), x), x)


def replacement1135(a, b, c, m, n, n2, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1136(a, b, c, d, m, n, n2, p, x):
    return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def With1137(a, b, c, m, n, n2, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement1137(a, b, c, m, n, n2, p, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k) + c*x**(S(2)*n/k))**p, x), x, x**k), x)


def With1138(a, b, c, d, m, n, n2, p, x):
    k = Denominator(m)
    return Dist(k/d, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(k*n) + c*d**(-S(2)*n)*x**(S(2)*k*n))**p, x), x, (d*x)**(S(1)/k)), x)


def replacement1139(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**n*n*p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), Int((d*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(a*b*(m - n + S(1)) - x**n*(S(2)*a*c*(m + n*(S(2)*p + S(-1)) + S(1)) - b**S(2)*(m + n*(p + S(-1)) + S(1))), x), x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b*n*p + c*x**n*(m + n*(S(2)*p + S(-1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(-1)) + S(1))), x)


def replacement1140(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**(-n)*n*p/(m + S(1)), Int((d*x)**(m + n)*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(1))), x)


def replacement1141(a, b, c, d, m, n, n2, p, x):
    return Dist(n*p/(m + S(2)*n*p + S(1)), Int((d*x)**m*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p/(d*(m + S(2)*n*p + S(1))), x)


def replacement1142(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**n/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**(m - n)*(b*(m - n + S(1)) + S(2)*c*x**n*(m + S(2)*n*(p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) + Simp(d**(n + S(-1))*(d*x)**(m - n + S(1))*(b + S(2)*c*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1143(a, b, c, d, m, n, n2, p, x):
    return Dist(d**(S(2)*n)/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**(m - S(2)*n)*(S(2)*a*(m - S(2)*n + S(1)) + b*x**n*(m + n*(S(2)*p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) - Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(S(2)*a + b*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1144(a, b, c, d, m, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1)) + b*c*x**n*(m + S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*d*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1145(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**(S(2)*n)/(c*(m + S(2)*n*p + S(1))), Int((d*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*(m - S(2)*n + S(1)) + b*x**n*(m + n*(p + S(-1)) + S(1)), x), x), x) + Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(m + S(2)*n*p + S(1))), x)


def replacement1146(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**(-n)/(a*(m + S(1))), Int((d*x)**(m + n)*(b*(m + n*(p + S(1)) + S(1)) + c*x**n*(m + S(2)*n*(p + S(1)) + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*d*(m + S(1))), x)


def replacement1147(a, b, c, d, m, n, n2, x):
    return -Dist(d**(-n)/a, Int((d*x)**(m + n)*(b + c*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Simp((d*x)**(m + S(1))/(a*d*(m + S(1))), x)


def replacement1148(a, b, c, m, n, n2, x):
    return Int(PolynomialDivide(x**m, a + b*x**n + c*x**(S(2)*n), x), x)


def replacement1149(a, b, c, d, m, n, n2, x):
    return -Dist(d**(S(2)*n)/c, Int((d*x)**(m - S(2)*n)*(a + b*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Simp(d**(S(2)*n + S(-1))*(d*x)**(m - S(2)*n + S(1))/(c*(m - S(2)*n + S(1))), x)


def With1150(a, b, c, x):
    q = Rt(a/c, S(2))
    return -Dist(S(1)/2, Int((q - x**S(2))/(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/2, Int((q + x**S(2))/(a + b*x**S(2) + c*x**S(4)), x), x)


def With1151(a, b, c, m, n, n2, x):
    q = Rt(a/c, S(2))
    r = Rt(-b/c + S(2)*q, S(2))
    return -Dist(1/(2*c*r), Int(x**(m - 3*n/2)*(q - r*x**(n/2))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*r), Int(x**(m - 3*n/2)*(q + r*x**(n/2))/(q + r*x**(n/2) + x**n), x), x)


def With1152(a, b, c, m, n, n2, x):
    q = Rt(a/c, S(2))
    r = Rt(-b/c + S(2)*q, S(2))
    return Dist(1/(2*c*r), Int(x**(m - n/2)/(q - r*x**(n/2) + x**n), x), x) - Dist(1/(2*c*r), Int(x**(m - n/2)/(q + r*x**(n/2) + x**n), x), x)


def With1153(a, b, c, d, m, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Dist(d**n*(b/q + S(-1))/S(2), Int((d*x)**(m - n)/(b/S(2) + c*x**n - q/S(2)), x), x) + Dist(d**n*(b/q + S(1))/S(2), Int((d*x)**(m - n)/(b/S(2) + c*x**n + q/S(2)), x), x)


def With1154(a, b, c, d, m, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(c/q, Int((d*x)**m/(b/S(2) + c*x**n - q/S(2)), x), x) - Dist(c/q, Int((d*x)**m/(b/S(2) + c*x**n + q/S(2)), x), x)


def With1155(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*sqrt(-c), Int(x**S(2)/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)


def With1156(a, b, c, x):
    q = Rt(c/a, S(2))
    return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1157(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b - q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1158(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(PosQ((b + q)/a), Not(And(PosQ((b - q)/a), SimplerSqrtQ((b - q)/(S(2)*a), (b + q)/(S(2)*a))))):
        return True
    return False


def replacement1158(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(x*(b + S(2)*c*x**S(2) + q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*(S(2)*a + x**S(2)*(b + q))*EllipticE(ArcTan(x*Rt((b + q)/(S(2)*a), S(2))), S(2)*q/(b + q))*Rt((b + q)/(S(2)*a), S(2))/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1159(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if PosQ((b - q)/a):
        return True
    return False


def replacement1159(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(x*(b + S(2)*c*x**S(2) - q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(sqrt((S(2)*a + x**S(2)*(b + q))/(S(2)*a + x**S(2)*(b - q)))*(S(2)*a + x**S(2)*(b - q))*EllipticE(ArcTan(x*Rt((b - q)/(S(2)*a), S(2))), -S(2)*q/(b - q))*Rt((b - q)/(S(2)*a), S(2))/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1160(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b + q)/a), Not(And(NegQ((b - q)/a), SimplerSqrtQ(-(b - q)/(S(2)*a), -(b + q)/(S(2)*a))))):
        return True
    return False


def replacement1160(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) + q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b + q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1161(a, b, c, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if NegQ((b - q)/a):
        return True
    return False


def replacement1161(a, b, c, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(1)/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((b - q)/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1162(a, b, c, x):
    q = Rt(c/a, S(2))
    return -Dist(S(1)/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(S(1)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1163(a, b, c, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(x**S(2)/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)


def replacement1164(a, b, c, m, n, n2, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)


def With1165(a, b, c, d, m, n, n2, p, x):
    k = Denominator(m)
    return -Dist(k/d, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*d**(-n)*x**(-k*n) + c*d**(-S(2)*n)*x**(-S(2)*k*n))**p, x), x, (d*x)**(-S(1)/k)), x)


def replacement1166(a, b, c, d, m, n, n2, p, x):
    return -Dist(d**IntPart(m)*(d*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)


def With1167(a, b, c, m, n, n2, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n) + c*x**(S(2)*k*n))**p, x), x, x**(S(1)/k)), x)


def replacement1168(a, b, c, d, m, n, n2, p, x):
    return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1169(a, b, c, m, n, n2, p, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p, x), x, x**(m + S(1))), x)


def replacement1170(a, b, c, d, m, n, n2, p, x):
    return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def With1171(a, b, c, d, m, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int((d*x)**m/(b + S(2)*c*x**n - q), x), x) - Dist(S(2)*c/q, Int((d*x)**m/(b + S(2)*c*x**n + q), x), x)


def replacement1172(a, b, c, d, m, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((d*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1)) + b*c*x**n*(m + S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp((d*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**n)/(a*d*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1173(a, b, c, d, m, n, n2, p, x):
    return Dist(a**IntPart(p)*(S(2)*c*x**n/(b - Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(S(2)*c*x**n/(b + Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d*x)**m*(S(2)*c*x**n/(b - sqrt(-S(4)*a*c + b**S(2))) + S(1))**p*(S(2)*c*x**n/(b + sqrt(-S(4)*a*c + b**S(2))) + S(1))**p, x), x)


def replacement1174(a, b, c, m, mn, n, p, x):
    return Int(x**(m - n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x)


def replacement1175(a, b, c, m, mn, n, p, x):
    return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(m - n*p)*(a*x**n + b + c*x**(S(2)*n))**p, x), x)


def replacement1176(a, b, c, d, m, mn, n, p, x):
    return Dist(d**IntPart(m)*x**(-FracPart(m))*(d*x)**FracPart(m), Int(x**m*(a + b*x**(-n) + c*x**n)**p, x), x)


def replacement1177(a, b, c, m, n, n2, p, v, x):
    return Dist(Coefficient(v, x, S(1))**(-m + S(-1)), Subst(Int(SimplifyIntegrand((x - Coefficient(v, x, S(0)))**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x), x, v), x)


def replacement1178(a, b, c, m, n, n2, p, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, v), x)


def replacement1179(a, b, c, d, e, n, n2, p, q, x):
    return Int(x**(n*(S(2)*p + q))*(d*x**(-n) + e)**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def replacement1180(a, c, d, e, n, n2, p, q, x):
    return Int(x**(n*(S(2)*p + q))*(a*x**(-S(2)*n) + c)**p*(d*x**(-n) + e)**q, x)


def replacement1181(a, b, c, d, e, n, n2, p, q, x):
    return -Subst(Int((d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p/x**S(2), x), x, S(1)/x)


def replacement1182(a, c, d, e, n, n2, p, q, x):
    return -Subst(Int((a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q/x**S(2), x), x, S(1)/x)


def With1183(a, b, c, d, e, n, n2, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g + S(-1))*(d + e*x**(g*n))**q*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p, x), x, x**(S(1)/g)), x)


def With1184(a, c, d, e, n, n2, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g + S(-1))*(a + c*x**(S(2)*g*n))**p*(d + e*x**(g*n))**q, x), x, x**(S(1)/g)), x)


def replacement1185(b, c, d, e, n, n2, p, x):
    return Dist(e/c, Int(x**(-n)*(b*x**n + c*x**(S(2)*n))**(p + S(1)), x), x) + Simp(x**(-S(2)*n*(p + S(1)))*(b*e - c*d)*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(b*c*n*(p + S(1))), x)


def replacement1186(b, c, d, e, n, n2, p, x):
    return Simp(e*x**(S(1) - n)*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(n*(S(2)*p + S(1)) + S(1))), x)


def replacement1187(b, c, d, e, n, n2, p, x):
    return -Dist((b*e*(n*p + S(1)) - c*d*(n*(S(2)*p + S(1)) + S(1)))/(c*(n*(S(2)*p + S(1)) + S(1))), Int((b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(e*x**(S(1) - n)*(b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(n*(S(2)*p + S(1)) + S(1))), x)


def replacement1188(b, c, d, e, n, n2, p, q, x):
    return Dist(x**(-n*FracPart(p))*(b + c*x**n)**(-FracPart(p))*(b*x**n + c*x**(S(2)*n))**FracPart(p), Int(x**(n*p)*(b + c*x**n)**p*(d + e*x**n)**q, x), x)


def replacement1189(a, b, c, d, e, n, n2, p, q, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p)*(d + e*x**n)**q, x), x)


def replacement1190(a, b, c, d, e, n, n2, p, q, x):
    return Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)


def replacement1191(a, c, d, e, n, n2, p, q, x):
    return Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)


def replacement1192(a, b, c, d, e, n, n2, p, q, x):
    return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)


def replacement1193(a, c, d, e, n, n2, p, q, x):
    return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)


def replacement1194(a, b, c, d, e, n, n2, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1195(a, c, d, e, n, n2, q, x):
    return Int(ExpandIntegrand((a + c*x**(S(2)*n))*(d + e*x**n)**q, x), x)


def replacement1196(a, b, c, d, e, n, n2, q, x):
    return Dist(S(1)/(d*e**S(2)*n*(q + S(1))), Int((d + e*x**n)**(q + S(1))*Simp(a*e**S(2)*(n*(q + S(1)) + S(1)) - b*d*e + c*d**S(2) + c*d*e*n*x**n*(q + S(1)), x), x), x) - Simp(x*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))/(d*e**S(2)*n*(q + S(1))), x)


def replacement1197(a, c, d, e, n, n2, q, x):
    return Dist(S(1)/(d*e**S(2)*n*(q + S(1))), Int((d + e*x**n)**(q + S(1))*Simp(a*e**S(2)*(n*(q + S(1)) + S(1)) + c*d**S(2) + c*d*e*n*x**n*(q + S(1)), x), x), x) - Simp(x*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))/(d*e**S(2)*n*(q + S(1))), x)


def replacement1198(a, b, c, d, e, n, n2, q, x):
    return Dist(S(1)/(e*(n*(q + S(2)) + S(1))), Int((d + e*x**n)**q*(a*e*(n*(q + S(2)) + S(1)) - x**n*(-b*e*(n*(q + S(2)) + S(1)) + c*d*(n + S(1)))), x), x) + Simp(c*x**(n + S(1))*(d + e*x**n)**(q + S(1))/(e*(n*(q + S(2)) + S(1))), x)


def replacement1199(a, c, d, e, n, n2, q, x):
    return Dist(S(1)/(e*(n*(q + S(2)) + S(1))), Int((d + e*x**n)**q*(a*e*(n*(q + S(2)) + S(1)) - c*d*x**n*(n + S(1))), x), x) + Simp(c*x**(n + S(1))*(d + e*x**n)**(q + S(1))/(e*(n*(q + S(2)) + S(1))), x)


def With1200(a, c, d, e, n, n2, x):
    q = Rt(S(2)*d*e, S(2))
    return Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*x**n - q*x**(n/S(2))), x), x) + Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*x**n + q*x**(n/S(2))), x), x)


def With1201(a, c, d, e, n, n2, x):
    q = Rt(-S(2)*d*e, S(2))
    return Dist(d/(S(2)*a), Int((d - q*x**(n/S(2)))/(d - e*x**n - q*x**(n/S(2))), x), x) + Dist(d/(S(2)*a), Int((d + q*x**(n/S(2)))/(d - e*x**n + q*x**(n/S(2))), x), x)


def With1202(a, c, d, e, x):
    q = Rt(a*c, S(2))
    return Dist((-a*e + d*q)/(S(2)*a*c), Int((-c*x**S(2) + q)/(a + c*x**S(4)), x), x) + Dist((a*e + d*q)/(S(2)*a*c), Int((c*x**S(2) + q)/(a + c*x**S(4)), x), x)


def With1203(a, c, d, e, n, n2, x):
    q = Rt(a/c, S(4))
    return Dist(sqrt(S(2))/(S(4)*c*q**S(3)), Int((sqrt(S(2))*d*q - x**(n/S(2))*(d - e*q**S(2)))/(q**S(2) - sqrt(S(2))*q*x**(n/S(2)) + x**n), x), x) + Dist(sqrt(S(2))/(S(4)*c*q**S(3)), Int((sqrt(S(2))*d*q + x**(n/S(2))*(d - e*q**S(2)))/(q**S(2) + sqrt(S(2))*q*x**(n/S(2)) + x**n), x), x)


def With1204(a, c, d, e, x):
    q = Rt(c/a, S(6))
    return Dist(S(1)/(S(6)*a*q**S(2)), Int((S(2)*d*q**S(2) - x*(sqrt(S(3))*d*q**S(3) - e))/(q**S(2)*x**S(2) - sqrt(S(3))*q*x + S(1)), x), x) + Dist(S(1)/(S(6)*a*q**S(2)), Int((S(2)*d*q**S(2) + x*(sqrt(S(3))*d*q**S(3) + e))/(q**S(2)*x**S(2) + sqrt(S(3))*q*x + S(1)), x), x) + Dist(S(1)/(S(3)*a*q**S(2)), Int((d*q**S(2) - e*x)/(q**S(2)*x**S(2) + S(1)), x), x)


def With1205(a, c, d, e, n, n2, x):
    q = Rt(-a/c, S(2))
    return Dist(d/S(2) - e*q/S(2), Int(S(1)/(a - c*q*x**n), x), x) + Dist(d/S(2) + e*q/S(2), Int(S(1)/(a + c*q*x**n), x), x)


def replacement1206(a, c, d, e, n, n2, x):
    return Dist(d, Int(S(1)/(a + c*x**(S(2)*n)), x), x) + Dist(e, Int(x**n/(a + c*x**(S(2)*n)), x), x)


def With1207(a, b, c, d, e, n, n2, x):
    q = Rt(-b/c + S(2)*d/e, S(2))
    return Dist(e**S(2)/(S(2)*c), Int(S(1)/(d - e*q*x**(n/S(2)) + e*x**n), x), x) + Dist(e**S(2)/(S(2)*c), Int(S(1)/(d + e*q*x**(n/S(2)) + e*x**n), x), x)


def With1208(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x)


def With1209(a, b, c, d, e, n, n2, x):
    q = Rt(a/c, S(2))
    r = Rt(-b/c + S(2)*q, S(2))
    return Dist(1/(2*c*q*r), Int((d*r - x**(n/2)*(d - e*q))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((d*r + x**(n/2)*(d - e*q))/(q + r*x**(n/2) + x**n), x), x)


def With1210(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int(S(1)/(b/S(2) + c*x**n - q/S(2)), x), x)


def With1211(a, b, c, d, e, n, n2, x):
    q = Rt(a/c, S(2))
    r = Rt(-b/c + S(2)*q, S(2))
    return Dist(1/(2*c*q*r), Int((d*r - x**(n/2)*(d - e*q))/(q - r*x**(n/2) + x**n), x), x) + Dist(1/(2*c*q*r), Int((d*r + x**(n/2)*(d - e*q))/(q + r*x**(n/2) + x**n), x), x)


def replacement1212(a, b, c, d, e, n, n2, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1213(a, c, d, e, n, n2, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)


def replacement1214(a, b, c, d, e, n, n2, q, x):
    return Dist(e**S(2)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x**n)**q, x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((d + e*x**n)**(q + S(1))*(-b*e + c*d - c*e*x**n)/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1215(a, c, d, e, n, n2, q, x):
    return Dist(c/(a*e**S(2) + c*d**S(2)), Int((d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(e**S(2)/(a*e**S(2) + c*d**S(2)), Int((d + e*x**n)**q, x), x)


def With1216(a, b, c, d, e, n, n2, q, x):
    r = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/r, Int((d + e*x**n)**q/(b + S(2)*c*x**n - r), x), x) - Dist(S(2)*c/r, Int((d + e*x**n)**q/(b + S(2)*c*x**n + r), x), x)


def With1217(a, c, d, e, n, n2, q, x):
    r = Rt(-a*c, S(2))
    return -Dist(c/(S(2)*r), Int((d + e*x**n)**q/(-c*x**n + r), x), x) - Dist(c/(S(2)*r), Int((d + e*x**n)**q/(c*x**n + r), x), x)


def replacement1218(a, b, c, d, e, n, n2, p, x):
    return Dist(n*p/(c*(S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(-a*b*e + S(2)*a*c*d*(S(2)*n*p + n + S(1)) + x**n*(S(2)*a*c*e*(S(2)*n*p + S(1)) - b**S(2)*e*(n*p + S(1)) + b*c*d*(S(2)*n*p + n + S(1))), x), x), x) + Simp(x*(a + b*x**n + c*x**(S(2)*n))**p*(b*e*n*p + c*d*(S(2)*n*p + n + S(1)) + c*e*x**n*(S(2)*n*p + S(1)))/(c*(S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), x)


def replacement1219(a, c, d, e, n, n2, p, x):
    return Dist(S(2)*a*n*p/((S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), Int((a + c*x**(S(2)*n))**(p + S(-1))*(d*(S(2)*n*p + n + S(1)) + e*x**n*(S(2)*n*p + S(1))), x), x) + Simp(x*(a + c*x**(S(2)*n))**p*(d*(S(2)*n*p + n + S(1)) + e*x**n*(S(2)*n*p + S(1)))/((S(2)*n*p + S(1))*(S(2)*n*p + n + S(1))), x)


def replacement1220(a, b, c, d, e, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e - S(2)*a*c*d*(S(2)*n*p + S(2)*n + S(1)) + b**S(2)*d*(n*p + n + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(S(2)*n*p + S(3)*n + S(1)), x), x), x) - Simp(x*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e - S(2)*a*c*d + b**S(2)*d + c*x**n*(-S(2)*a*e + b*d))/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1221(a, c, d, e, n, n2, p, x):
    return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((a + c*x**(S(2)*n))**(p + S(1))*(d*(S(2)*n*p + S(2)*n + S(1)) + e*x**n*(S(2)*n*p + S(3)*n + S(1))), x), x) - Simp(x*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*n*(p + S(1))), x)


def With1222(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*sqrt(-c), Int((d + e*x**S(2))/(sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)


def With1223(a, c, d, e, x):
    q = Rt(-a*c, S(2))
    return Dist(sqrt(-c), Int((d + e*x**S(2))/(sqrt(-c*x**S(2) + q)*sqrt(c*x**S(2) + q)), x), x)


def With1224(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(4))
    if ZeroQ(d*q**S(2) + e):
        return True
    return False


def replacement1224(a, b, c, d, e, x):

    q = Rt(c/a, S(4))
    return -Simp(d*x*sqrt(a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(q*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1225(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(2))
    if NonzeroQ(d*q + e):
        return True
    return False


def replacement1225(a, b, c, d, e, x):

    q = Rt(c/a, S(2))
    return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1226(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if ZeroQ(S(2)*c*d - e*(b - q)):
        return True
    return False


def replacement1226(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Simp(e*x*(b + S(2)*c*x**S(2) + q)/(S(2)*c*sqrt(a + b*x**S(2) + c*x**S(4))), x) - Simp(e*q*sqrt((S(2)*a + x**S(2)*(b + q))/q)*sqrt((S(2)*a + x**S(2)*(b - q))/(S(2)*a + x**S(2)*(b + q)))*EllipticE(asin(sqrt(S(2))*x/sqrt((S(2)*a + x**S(2)*(b + q))/q)), (b + q)/(S(2)*q))/(S(2)*c*sqrt(a/(S(2)*a + x**S(2)*(b + q)))*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1227(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-a*c, S(2))
    if And(ZeroQ(c*d + e*q), IntegerQ(q)):
        return True
    return False


def replacement1227(a, c, d, e, x):

    q = Rt(-a*c, S(2))
    return Simp(e*x*(c*x**S(2) + q)/(c*sqrt(a + c*x**S(4))), x) - Simp(sqrt(S(2))*e*q*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticE(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(c*sqrt(-a)*sqrt(a + c*x**S(4))), x)


def With1228(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-a*c, S(2))
    if ZeroQ(c*d + e*q):
        return True
    return False


def replacement1228(a, c, d, e, x):

    q = Rt(-a*c, S(2))
    return Simp(e*x*(c*x**S(2) + q)/(c*sqrt(a + c*x**S(4))), x) - Simp(sqrt(S(2))*e*q*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticE(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(c*sqrt(a/(a + q*x**S(2)))*sqrt(a + c*x**S(4))), x)


def With1229(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if NonzeroQ(S(2)*c*d - e*(b - q)):
        return True
    return False


def replacement1229(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b - q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1230(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-a*c, S(2))
    if NonzeroQ(c*d + e*q):
        return True
    return False


def replacement1230(a, c, d, e, x):

    q = Rt(-a*c, S(2))
    return -Dist(e/c, Int((-c*x**S(2) + q)/sqrt(a + c*x**S(4)), x), x) + Dist((c*d + e*q)/c, Int(S(1)/sqrt(a + c*x**S(4)), x), x)


def With1231(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if Or(PosQ((b + q)/a), PosQ((b - q)/a)):
        return True
    return False


def replacement1231(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(d, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist(e, Int(x**S(2)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def replacement1232(a, c, d, e, x):
    return Dist(d, Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist(e, Int(x**S(2)/sqrt(a + c*x**S(4)), x), x)


def With1233(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b + q)/a), ZeroQ(S(2)*c*d - e*(b + q)), Not(SimplerSqrtQ(-(b - q)/(S(2)*a), -(b + q)/(S(2)*a)))):
        return True
    return False


def replacement1233(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Simp(a*e*sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticE(asin(x*Rt(-(b + q)/(S(2)*a), S(2))), (b - q)/(b + q))*Rt(-(b + q)/(S(2)*a), S(2))/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1234(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b + q)/a), NonzeroQ(S(2)*c*d - e*(b + q)), Not(SimplerSqrtQ(-(b - q)/(S(2)*a), -(b + q)/(S(2)*a)))):
        return True
    return False


def replacement1234(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) + q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b + q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1235(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b - q)/a), ZeroQ(S(2)*c*d - e*(b - q))):
        return True
    return False


def replacement1235(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Simp(a*e*sqrt(S(1) + x**S(2)*(b - q)/(S(2)*a))*sqrt(S(1) + x**S(2)*(b + q)/(S(2)*a))*EllipticE(asin(x*Rt(-(b - q)/(S(2)*a), S(2))), (b + q)/(b - q))*Rt(-(b - q)/(S(2)*a), S(2))/(c*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1236(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if And(NegQ((b - q)/a), NonzeroQ(S(2)*c*d - e*(b - q))):
        return True
    return False


def replacement1236(a, b, c, d, e, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/(S(2)*c), Int((b + S(2)*c*x**S(2) - q)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((S(2)*c*d - e*(b - q))/(S(2)*c), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1237(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(4))
    if ZeroQ(d*q**S(2) + e):
        return True
    return False


def replacement1237(a, b, c, d, e, x):

    q = Rt(c/a, S(4))
    return -Simp(d*x*sqrt(a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(q*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1238(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(4))
    if ZeroQ(d*q**S(2) + e):
        return True
    return False


def replacement1238(a, c, d, e, x):

    q = Rt(c/a, S(4))
    return -Simp(d*x*sqrt(a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))), x) + Simp(d*sqrt((a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticE(S(2)*ArcTan(q*x), S(1)/2)/(q*sqrt(a + c*x**S(4))), x)


def With1239(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(2))
    if NonzeroQ(d*q + e):
        return True
    return False


def replacement1239(a, b, c, d, e, x):

    q = Rt(c/a, S(2))
    return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x)


def With1240(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(2))
    if NonzeroQ(d*q + e):
        return True
    return False


def replacement1240(a, c, d, e, x):

    q = Rt(c/a, S(2))
    return -Dist(e/q, Int((-q*x**S(2) + S(1))/sqrt(a + c*x**S(4)), x), x) + Dist((d*q + e)/q, Int(S(1)/sqrt(a + c*x**S(4)), x), x)


def replacement1241(a, c, d, e, x):
    return Dist(d/sqrt(a), Int(sqrt(S(1) + e*x**S(2)/d)/sqrt(S(1) - e*x**S(2)/d), x), x)


def replacement1242(a, c, d, e, x):
    return Dist(sqrt(S(1) + c*x**S(4)/a)/sqrt(a + c*x**S(4)), Int((d + e*x**S(2))/sqrt(S(1) + c*x**S(4)/a), x), x)


def With1243(a, c, d, e, x):
    q = Rt(-c/a, S(2))
    return Dist(e/q, Int((q*x**S(2) + S(1))/sqrt(a + c*x**S(4)), x), x) + Dist((d*q - e)/q, Int(S(1)/sqrt(a + c*x**S(4)), x), x)


def With1244(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int((d + e*x**S(2))/(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)


def replacement1245(a, b, c, d, e, n, n2, p, x):
    return Int(ExpandIntegrand((d + e*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1246(a, c, d, e, n, n2, p, x):
    return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p*(d + e*x**n), x), x)


def replacement1247(a, b, c, d, e, n, n2, p, q, x):
    return Int((d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(S(2)*n*p - n + S(1))/(e*(S(2)*n*p + n*q + S(1))) - c**p*x**(S(2)*n*p) + (a + b*x**n + c*x**(S(2)*n))**p, x), x) + Simp(c**p*x**(S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(S(2)*n*p + n*q + S(1))), x)


def replacement1248(a, c, d, e, n, n2, p, q, x):
    return Int((d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(S(2)*n*p - n + S(1))/(e*(S(2)*n*p + n*q + S(1))) - c**p*x**(S(2)*n*p) + (a + c*x**(S(2)*n))**p, x), x) + Simp(c**p*x**(S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(S(2)*n*p + n*q + S(1))), x)


def replacement1249(a, b, c, d, e, x):
    return -Dist(e**(S(-2)), Int((-b*e + c*d - c*e*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((a*e**S(2) - b*d*e + c*d**S(2))/e**S(2), Int(S(1)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)


def replacement1250(a, c, d, e, x):
    return -Dist(c/e**S(2), Int((d - e*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist((a*e**S(2) + c*d**S(2))/e**S(2), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)


def With1251(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Dist(e**(S(-4)), Int(Simp(-c**S(2)*e**S(3)*x**S(6) + c*e**S(2)*x**S(4)*(-S(2)*b*e + c*d) - S(2)*c*(a*e**S(2) - b*d*e + c*d**S(2))**S(2)/(S(2)*c*d - e*(b + q)) - e*x**S(2)*(b**S(2)*e**S(2) + c**S(2)*d**S(2) - S(2)*c*e*(-a*e + b*d)) + (-b*e + c*d)*(S(2)*a*e**S(2) - b*d*e + c*d**S(2)), x)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((a*e**S(2) - b*d*e + c*d**S(2))**S(2)/(e**S(3)*(S(2)*c*d - e*(b + q))), Int((b + S(2)*c*x**S(2) + q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)


def With1252(a, c, d, e, x):
    q = Rt(-a*c, S(2))
    return -Dist(c/e**S(4), Int(Simp(c*d*e**S(2)*x**S(4) - c*e**S(3)*x**S(6) + d*(S(2)*a*e**S(2) + c*d**S(2)) - e*x**S(2)*(S(2)*a*e**S(2) + c*d**S(2)) - (a*e**S(2) + c*d**S(2))**S(2)/(c*d - e*q), x)/sqrt(a + c*x**S(4)), x), x) - Dist((a*e**S(2) + c*d**S(2))**S(2)/(e**S(3)*(c*d - e*q)), Int((c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)


def replacement1253(a, b, c, d, e, p, x):
    return Dist(a, Int((a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(b, Int(x**S(2)*(a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(c, Int(x**S(4)*(a + b*x**S(2) + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x)


def replacement1254(a, c, d, e, p, x):
    return Dist(a, Int((a + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x) + Dist(c, Int(x**S(4)*(a + c*x**S(4))**(p + S(-1))/(d + e*x**S(2)), x), x)


def With1255(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*sqrt(-c), Int(S(1)/((d + e*x**S(2))*sqrt(-b - S(2)*c*x**S(2) + q)*sqrt(b + S(2)*c*x**S(2) + q)), x), x)


def With1256(a, c, d, e, x):
    q = Rt(-a*c, S(2))
    return Dist(sqrt(-c), Int(S(1)/((d + e*x**S(2))*sqrt(-c*x**S(2) + q)*sqrt(c*x**S(2) + q)), x), x)


def With1257(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/(S(2)*c*d - e*(b - q)), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist(e/(S(2)*c*d - e*(b - q)), Int((b + S(2)*c*x**S(2) - q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)


def With1258(a, c, d, e, x):
    q = Rt(-a*c, S(2))
    return Dist(c/(c*d + e*q), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist(e/(c*d + e*q), Int((-c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)


def With1259(a, b, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(4))
    if NonzeroQ(-d*q**S(2) + e):
        return True
    return False


def replacement1259(a, b, c, d, e, x):

    q = Rt(c/a, S(4))
    return -Dist(q**S(2)/(-d*q**S(2) + e), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Simp(ArcTan(x*sqrt((a*e**S(2) - b*d*e + c*d**S(2))/(d*e))/sqrt(a + b*x**S(2) + c*x**S(4)))/(S(2)*d*sqrt((a*e**S(2) - b*d*e + c*d**S(2))/(d*e))), x) + Simp(sqrt((a + b*x**S(2) + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(d*q**S(2) + e)*(q**S(2)*x**S(2) + S(1))*EllipticPi(-(-d*q**S(2) + e)**S(2)/(S(4)*d*e*q**S(2)), S(2)*ArcTan(q*x), -b*q**S(2)/(S(4)*c) + S(1)/2)/(S(4)*d*q*(-d*q**S(2) + e)*sqrt(a + b*x**S(2) + c*x**S(4))), x)


def With1260(a, c, d, e, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(c/a, S(4))
    if NonzeroQ(-d*q**S(2) + e):
        return True
    return False


def replacement1260(a, c, d, e, x):

    q = Rt(c/a, S(4))
    return -Dist(q**S(2)/(-d*q**S(2) + e), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Simp(ArcTan(x*sqrt((a*e**S(2) + c*d**S(2))/(d*e))/sqrt(a + c*x**S(4)))/(S(2)*d*sqrt((a*e**S(2) + c*d**S(2))/(d*e))), x) + Simp(sqrt((a + c*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(d*q**S(2) + e)*(q**S(2)*x**S(2) + S(1))*EllipticPi(-(-d*q**S(2) + e)**S(2)/(S(4)*d*e*q**S(2)), S(2)*ArcTan(q*x), S(1)/2)/(S(4)*d*q*sqrt(a + c*x**S(4))*(-d*q**S(2) + e)), x)


def With1261(a, c, d, e, x):
    q = Rt(-c/a, S(4))
    return Simp(EllipticPi(-e/(d*q**S(2)), asin(q*x), S(-1))/(sqrt(a)*d*q), x)


def replacement1262(a, c, d, e, x):
    return Dist(sqrt(S(1) + c*x**S(4)/a)/sqrt(a + c*x**S(4)), Int(S(1)/(sqrt(S(1) + c*x**S(4)/a)*(d + e*x**S(2))), x), x)


def With1263(a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), Int(S(1)/((d + e*x**S(2))*sqrt(S(2)*c*x**S(2)/(b - q) + S(1))*sqrt(S(2)*c*x**S(2)/(b + q) + S(1))), x), x)


def replacement1264(a, b, c, d, e, p, x):
    return -Dist(S(1)/(S(2)*a*(p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), Int((a + b*x**S(2) + c*x**S(4))**(p + S(1))*Simp(-a*b*c*d*e*(S(8)*p + S(11)) + S(2)*a*c*(S(4)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(4)*p + S(5))) + b**S(3)*d*e*(S(2)*p + S(3)) - b**S(2)*(S(2)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(2)*p + S(3))) - c*e*x**S(4)*(S(4)*p + S(7))*(S(2)*a*c*e - b**S(2)*e + b*c*d) - x**S(2)*(S(4)*a*c**S(2)*d*e - b**S(3)*e**S(2)*(S(2)*p + S(3)) - S(2)*b**S(2)*c*d*e*(p + S(2)) + b*c*(a*e**S(2)*(S(8)*p + S(11)) + c*d**S(2)*(S(4)*p + S(7)))), x)/(d + e*x**S(2)), x), x) - Simp(x*(a + b*x**S(2) + c*x**S(4))**(p + S(1))*(S(3)*a*b*c*e - S(2)*a*c**S(2)*d - b**S(3)*e + b**S(2)*c*d + c*x**S(2)*(S(2)*a*c*e - b**S(2)*e + b*c*d))/(S(2)*a*(p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement1265(a, c, d, e, p, x):
    return -Dist(-S(1)/(S(8)*a**S(2)*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), Int((a + c*x**S(4))**(p + S(1))*Simp(-S(4)*a*c**S(2)*d*e*x**S(2) - S(2)*a*c**S(2)*e**S(2)*x**S(4)*(S(4)*p + S(7)) + S(2)*a*c*(S(4)*a*e**S(2)*(p + S(1)) + c*d**S(2)*(S(4)*p + S(5))), x)/(d + e*x**S(2)), x), x) - Simp(-x*(a + c*x**S(4))**(p + S(1))*(-S(2)*a*c**S(2)*d + S(2)*a*c**S(2)*e*x**S(2))/(S(8)*a**S(2)*c*(p + S(1))*(a*e**S(2) + c*d**S(2))), x)


def replacement1266(a, b, c, d, e, x):
    return -Dist(c/(S(2)*d*(a*e**S(2) - b*d*e + c*d**S(2))), Int((d + e*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) + Dist((a*e**S(2) - S(2)*b*d*e + S(3)*c*d**S(2))/(S(2)*d*(a*e**S(2) - b*d*e + c*d**S(2))), Int(S(1)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x) + Simp(e**S(2)*x*sqrt(a + b*x**S(2) + c*x**S(4))/(S(2)*d*(d + e*x**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))), x)


def replacement1267(a, c, d, e, x):
    return -Dist(c/(S(2)*d*(a*e**S(2) + c*d**S(2))), Int((d + e*x**S(2))/sqrt(a + c*x**S(4)), x), x) + Dist((a*e**S(2) + S(3)*c*d**S(2))/(S(2)*d*(a*e**S(2) + c*d**S(2))), Int(S(1)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x) + Simp(e**S(2)*x*sqrt(a + c*x**S(4))/(S(2)*d*(d + e*x**S(2))*(a*e**S(2) + c*d**S(2))), x)


def With1268(a, b, c, d, e, p, q, x):
    r = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(a**IntPart(p)*(S(2)*c*x**S(2)/(b - r) + S(1))**(-FracPart(p))*(S(2)*c*x**S(2)/(b + r) + S(1))**(-FracPart(p))*(a + b*x**S(2) + c*x**S(4))**FracPart(p), Int((d + e*x**S(2))**q*(S(2)*c*x**S(2)/(b - r) + S(1))**p*(S(2)*c*x**S(2)/(b + r) + S(1))**p, x), x)


def With1269(a, c, d, e, p, q, x):
    r = Rt(-a*c, S(2))
    return Dist(a**IntPart(p)*(a + c*x**S(4))**FracPart(p)*(-c*x**S(2)/r + S(1))**(-FracPart(p))*(c*x**S(2)/r + S(1))**(-FracPart(p)), Int((d + e*x**S(2))**q*(-c*x**S(2)/r + S(1))**p*(c*x**S(2)/r + S(1))**p, x), x)


def replacement1270(a, b, c, d, e, x):
    return Simp(EllipticF(S(2)*asin(x*Rt(-e/d, S(2))), b*d/(S(4)*a*e))/(S(2)*sqrt(a)*sqrt(d)*Rt(-e/d, S(2))), x)


def replacement1271(a, b, c, d, e, x):
    return Dist(sqrt((a + b*x**S(2) + c*x**S(4))/a)*sqrt((d + e*x**S(2))/d)/(sqrt(d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), Int(S(1)/(sqrt(S(1) + e*x**S(2)/d)*sqrt(S(1) + b*x**S(2)/a + c*x**S(4)/a)), x), x)


def replacement1272(a, b, c, d, e, x):
    return Simp(sqrt(a)*EllipticE(S(2)*asin(x*Rt(-e/d, S(2))), b*d/(S(4)*a*e))/(S(2)*sqrt(d)*Rt(-e/d, S(2))), x)


def replacement1273(a, b, c, d, e, x):
    return Dist(sqrt((d + e*x**S(2))/d)*sqrt(a + b*x**S(2) + c*x**S(4))/(sqrt((a + b*x**S(2) + c*x**S(4))/a)*sqrt(d + e*x**S(2))), Int(sqrt(S(1) + b*x**S(2)/a + c*x**S(4)/a)/sqrt(S(1) + e*x**S(2)/d), x), x)


def replacement1274(a, b, c, d, e, n, n2, p, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1275(a, c, d, e, n, n2, p, q, x):
    return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1276(a, c, d, e, n, n2, p, q, x):
    return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p, (d/(d**S(2) - e**S(2)*x**(S(2)*n)) - e*x**n/(d**S(2) - e**S(2)*x**(S(2)*n)))**(-q), x), x)


def replacement1277(a, b, c, d, e, n, n2, p, q, x):
    return Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1278(a, c, d, e, n, n2, p, q, x):
    return Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x)


def replacement1279(a, b, c, d, e, n, n2, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x, u), x)


def replacement1280(a, c, d, e, n, n2, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x, u), x)


def replacement1281(a, b, c, d, e, mn, n, n2, p, q, x):
    return Int(x**(-n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1282(a, c, d, e, mn, n2, p, q, x):
    return Int(x**(mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x)


def replacement1283(a, b, c, d, e, mn, n, n2, p, q, x):
    return Int(x**(S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def replacement1284(a, c, d, e, mn, n2, p, q, x):
    return Int(x**(-S(2)*mn*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x)


def replacement1285(a, b, c, d, e, mn, n, n2, p, q, x):
    return Dist(x**(n*FracPart(q))*(d + e*x**(-n))**FracPart(q)*(d*x**n + e)**(-FracPart(q)), Int(x**(-n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1286(a, c, d, e, mn, n2, p, q, x):
    return Dist(x**(-mn*FracPart(q))*(d + e*x**mn)**FracPart(q)*(d*x**(-mn) + e)**(-FracPart(q)), Int(x**(mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x), x)


def replacement1287(a, b, c, d, e, mn, n, n2, p, q, x):
    return Dist(x**(-S(2)*n*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p)*(a*x**(-S(2)*n) + b*x**(-n) + c)**(-FracPart(p)), Int(x**(S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x), x)


def replacement1288(a, c, d, e, mn, n2, p, q, x):
    return Dist(x**(-n2*FracPart(p))*(a + c*x**n2)**FracPart(p)*(a*x**(S(2)*mn) + c)**(-FracPart(p)), Int(x**(n2*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x), x)


def replacement1289(a, b, c, d, e, mn, n, p, q, x):
    return Int(x**(-n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x)


def replacement1290(a, b, c, d, e, mn, n, p, q, x):
    return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(-n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x), x)


def replacement1291(a, b, c, d, e, f, g, n, n2, p, q, r, x):
    return Dist((S(4)*c)**(-IntPart(p))*(b + S(2)*c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((b + S(2)*c*x**n)**(S(2)*p)*(d + e*x**n)**q*(f + g*x**n)**r, x), x)


def replacement1292(a, b, c, d, e, f, g, n, n2, p, q, r, x):
    return Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x)


def replacement1293(a, c, d, e, f, g, n, n2, p, q, r, x):
    return Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x)


def replacement1294(a, b, c, d, e, f, g, n, n2, p, q, r, x):
    return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x), x)


def replacement1295(a, c, d, e, f, g, n, n2, p, q, r, x):
    return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((d + e*x**n)**(p + q)*(f + g*x**n)**r*(a/d + c*x**n/e)**p, x), x)


def With1296(a, b, c, d, e, f, g, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    if NonzeroQ(S(2)*c*f - g*(b - q)):
        return True
    return False


def replacement1296(a, b, c, d, e, f, g, x):

    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist((S(2)*c*f - g*(b - q))/(S(2)*c*d - e*(b - q)), Int(S(1)/sqrt(a + b*x**S(2) + c*x**S(4)), x), x) - Dist((-d*g + e*f)/(S(2)*c*d - e*(b - q)), Int((b + S(2)*c*x**S(2) - q)/((d + e*x**S(2))*sqrt(a + b*x**S(2) + c*x**S(4))), x), x)


def With1297(a, c, d, e, f, g, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(-a*c, S(2))
    if NonzeroQ(c*f + g*q):
        return True
    return False


def replacement1297(a, c, d, e, f, g, x):

    q = Rt(-a*c, S(2))
    return Dist((c*f + g*q)/(c*d + e*q), Int(S(1)/sqrt(a + c*x**S(4)), x), x) + Dist((-d*g + e*f)/(c*d + e*q), Int((-c*x**S(2) + q)/(sqrt(a + c*x**S(4))*(d + e*x**S(2))), x), x)


def replacement1298(a, b, c, d1, d2, e1, e2, n, n2, non2, p, q, x):
    return Int((d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1299(a, b, c, d1, d2, e1, e2, n, n2, non2, p, q, x):
    return Dist((d1 + e1*x**(n/S(2)))**FracPart(q)*(d2 + e2*x**(n/S(2)))**FracPart(q)*(d1*d2 + e1*e2*x**n)**(-FracPart(q)), Int((d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1300(A, B, a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(A, Int((d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Dist(B, Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1301(A, B, a, c, d, e, m, n, n2, p, q, x):
    return Dist(A, Int((a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x) + Dist(B, Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1302(a, b, c, e, f, m, n, n2, p, q, x):
    return Dist(e**(S(1) - (m + S(1))/n)*f**m/n, Subst(Int((e*x)**(q + S(-1) + (m + S(1))/n)*(a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1303(a, c, e, f, m, n, n2, p, q, x):
    return Dist(e**(S(1) - (m + S(1))/n)*f**m/n, Subst(Int((e*x)**(q + S(-1) + (m + S(1))/n)*(a + c*x**S(2))**p, x), x, x**n), x)


def replacement1304(a, b, c, e, f, m, n, n2, p, q, x):
    return Dist(e**IntPart(q)*f**m*x**(-n*FracPart(q))*(e*x**n)**FracPart(q), Int(x**(m + n*q)*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1305(a, c, e, f, m, n, n2, p, q, x):
    return Dist(e**IntPart(q)*f**m*x**(-n*FracPart(q))*(e*x**n)**FracPart(q), Int(x**(m + n*q)*(a + c*x**(S(2)*n))**p, x), x)


def replacement1306(a, b, c, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1307(a, c, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(e*x**n)**q*(a + c*x**(S(2)*n))**p, x), x)


def replacement1308(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/n, Subst(Int((d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1309(a, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/n, Subst(Int((a + c*x**S(2))**p*(d + e*x)**q, x), x, x**n), x)


def replacement1310(a, b, c, d, e, m, n, n2, p, q, x):
    return Int(x**(m + n*(S(2)*p + q))*(d*x**(-n) + e)**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def replacement1311(a, c, d, e, m, n, n2, p, q, x):
    return Int(x**(m + n*(S(2)*p + q))*(a*x**(-S(2)*n) + c)**p*(d*x**(-n) + e)**q, x)


def replacement1312(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1313(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist(c**(-IntPart(p))*(b/S(2) + c*x**n)**(-S(2)*FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((f*x)**m*(b/S(2) + c*x**n)**(S(2)*p)*(d + e*x**n)**q, x), x)


def replacement1314(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(d + e*x)**q*(a + b*x + c*x**S(2))**p, x), x, x**n), x)


def replacement1315(a, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + c*x**S(2))**p*(d + e*x)**q, x), x, x**n), x)


def replacement1316(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1317(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1318(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)


def replacement1319(a, c, d, e, f, m, n, n2, p, q, x):
    return Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x)


def replacement1320(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist((d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p))*(a + b*x**n + c*x**(S(2)*n))**FracPart(p), Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)


def replacement1321(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist((a + c*x**(S(2)*n))**FracPart(p)*(d + e*x**n)**(-FracPart(p))*(a/d + c*x**n/e)**(-FracPart(p)), Int((f*x)**m*(d + e*x**n)**(p + q)*(a/d + c*x**n/e)**p, x), x)


def replacement1322(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(e**(-S(2)*p - (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**Mod(m, n)*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p + (m - Mod(m, n))/n)*n*x**(m - Mod(m, n))*(q + S(1))*(a + b*x**n + c*x**(S(2)*n))**p - (-d)**(S(-1) + (m - Mod(m, n))/n)*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1)))*(a*e**S(2) - b*d*e + c*d**S(2))**p)/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))**p/(n*(q + S(1))), x)


def replacement1323(a, c, d, e, m, n, n2, p, q, x):
    return Dist(e**(-S(2)*p - (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**Mod(m, n)*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p + (m - Mod(m, n))/n)*n*x**(m - Mod(m, n))*(a + c*x**(S(2)*n))**p*(q + S(1)) - (-d)**(S(-1) + (m - Mod(m, n))/n)*(a*e**S(2) + c*d**S(2))**p*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1))))/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))**p/(n*(q + S(1))), x)


def replacement1324(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(e**(-S(2)*p)*(-d)**(S(-1) + (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**m*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p)*n*(-d)**(S(1) - (m - Mod(m, n))/n)*(q + S(1))*(a + b*x**n + c*x**(S(2)*n))**p - e**(-(m - Mod(m, n))/n)*x**(-m + Mod(m, n))*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1)))*(a*e**S(2) - b*d*e + c*d**S(2))**p)/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))**p/(n*(q + S(1))), x)


def replacement1325(a, c, d, e, m, n, n2, p, q, x):
    return Dist(e**(-S(2)*p)*(-d)**(S(-1) + (m - Mod(m, n))/n)/(n*(q + S(1))), Int(x**m*(d + e*x**n)**(q + S(1))*ExpandToSum(Together((e**(S(2)*p)*n*(-d)**(S(1) - (m - Mod(m, n))/n)*(a + c*x**(S(2)*n))**p*(q + S(1)) - e**(-(m - Mod(m, n))/n)*x**(-m + Mod(m, n))*(a*e**S(2) + c*d**S(2))**p*(d*(Mod(m, n) + S(1)) + e*x**n*(n*(q + S(1)) + Mod(m, n) + S(1))))/(d + e*x**n)), x), x), x) + Simp(e**(-S(2)*p - (m - Mod(m, n))/n)*x**(Mod(m, n) + S(1))*(-d)**(S(-1) + (m - Mod(m, n))/n)*(d + e*x**n)**(q + S(1))*(a*e**S(2) + c*d**S(2))**p/(n*(q + S(1))), x)


def replacement1326(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist(S(1)/(e*(m + S(2)*n*p + n*q + S(1))), Int((f*x)**m*(d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(m + S(2)*n*p - n + S(1)) + e*(-c**p*x**(S(2)*n*p) + (a + b*x**n + c*x**(S(2)*n))**p)*(m + S(2)*n*p + n*q + S(1)), x), x), x) + Simp(c**p*f**(-S(2)*n*p + n + S(-1))*(f*x)**(m + S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(m + S(2)*n*p + n*q + S(1))), x)


def replacement1327(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(S(1)/(e*(m + S(2)*n*p + n*q + S(1))), Int((f*x)**m*(d + e*x**n)**q*ExpandToSum(-c**p*d*x**(S(2)*n*p - n)*(m + S(2)*n*p - n + S(1)) + e*(-c**p*x**(S(2)*n*p) + (a + c*x**(S(2)*n))**p)*(m + S(2)*n*p + n*q + S(1)), x), x), x) + Simp(c**p*f**(-S(2)*n*p + n + S(-1))*(f*x)**(m + S(2)*n*p - n + S(1))*(d + e*x**n)**(q + S(1))/(e*(m + S(2)*n*p + n*q + S(1))), x)


def replacement1328(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1329(a, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def With1330(a, b, c, d, e, m, n, n2, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement1330(a, b, c, d, e, m, n, n2, p, q, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(d + e*x**(n/k))**q*(a + b*x**(n/k) + c*x**(S(2)*n/k))**p, x), x, x**k), x)


def With1331(a, c, d, e, m, n, n2, p, q, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    k = GCD(m + S(1), n)
    if Unequal(k, S(1)):
        return True
    return False


def replacement1331(a, c, d, e, m, n, n2, p, q, x):

    k = GCD(m + S(1), n)
    return Dist(S(1)/k, Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + c*x**(S(2)*n/k))**p*(d + e*x**(n/k))**q, x), x, x**k), x)


def With1332(a, b, c, d, e, f, m, n, n2, p, q, x):
    k = Denominator(m)
    return Dist(k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))*(d + e*f**(-n)*x**(k*n))**q*(a + b*f**(-n)*x**(k*n) + c*f**(-S(2)*n)*x**(S(2)*k*n))**p, x), x, (f*x)**(S(1)/k)), x)


def With1333(a, c, d, e, f, m, n, n2, p, q, x):
    k = Denominator(m)
    return Dist(k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + c*x**(S(2)*k*n)/f)**p*(d + e*x**(k*n)/f)**q, x), x, (f*x)**(S(1)/k)), x)


def replacement1334(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(f**(-n)*n*p/((m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(S(2)*a*e*(m + S(1)) - b*d*(m + n*(S(2)*p + S(1)) + S(1)) + x**n*(b*e*(m + S(1)) - S(2)*c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp((f*x)**(m + S(1))*(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(1)))*(a + b*x**n + c*x**(S(2)*n))**p/(f*(m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1335(a, c, d, e, f, m, n, n2, p, x):
    return Dist(S(2)*f**(-n)*n*p/((m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**(p + S(-1))*(a*e*(m + S(1)) - c*d*x**n*(m + n*(S(2)*p + S(1)) + S(1))), x), x) + Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**p*(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(1)))/(f*(m + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1336(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(n*p/(c*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))*Simp(-a*b*e*(m + S(1)) + S(2)*a*c*d*(m + n*(S(2)*p + S(1)) + S(1)) + x**n*(S(2)*a*c*e*(m + S(2)*n*p + S(1)) - b**S(2)*e*(m + n*p + S(1)) + b*c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**p*(b*e*n*p + c*d*(m + n*(S(2)*p + S(1)) + S(1)) + c*e*x**n*(m + S(2)*n*p + S(1)))/(c*f*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1337(a, c, d, e, f, m, n, n2, p, x):
    return Dist(S(2)*a*n*p/((m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*Simp(d*(m + n*(S(2)*p + S(1)) + S(1)) + e*x**n*(m + S(2)*n*p + S(1)), x), x), x) + Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**p*(c*d*(m + n*(S(2)*p + S(1)) + S(1)) + c*e*x**n*(m + S(2)*n*p + S(1)))/(c*f*(m + S(2)*n*p + S(1))*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1338(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(f**n/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(x**n*(b*e - S(2)*c*d)*(m + S(2)*n*p + S(2)*n + S(1)) + (-S(2)*a*e + b*d)*(-m + n + S(-1)), x), x), x) + Simp(f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-S(2)*a*e + b*d - x**n*(b*e - S(2)*c*d))/(n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1339(a, c, d, e, f, m, n, n2, p, x):
    return Dist(f**n/(S(2)*a*c*n*(p + S(1))), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**(p + S(1))*(a*e*(-m + n + S(-1)) + c*d*x**n*(m + S(2)*n*p + S(2)*n + S(1))), x), x) + Simp(f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(a*e - c*d*x**n)/(S(2)*a*c*n*(p + S(1))), x)


def replacement1340(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e*(m + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(m + n*(S(2)*p + S(3)) + S(1)) + d*(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1))), x), x), x) - Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e + c*x**n*(-S(2)*a*e + b*d) + d*(-S(2)*a*c + b**S(2)))/(a*f*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1341(a, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(1))*Simp(d*(m + S(2)*n*(p + S(1)) + S(1)) + e*x**n*(m + n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*f*n*(p + S(1))), x)


def replacement1342(a, b, c, d, e, f, m, n, n2, p, x):
    return -Dist(f**n/(c*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*e*(m - n + S(1)) + x**n*(b*e*(m + n*p + S(1)) - c*d*(m + n*(S(2)*p + S(1)) + S(1))), x), x), x) + Simp(e*f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(c*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1343(a, c, d, e, f, m, n, n2, p, x):
    return -Dist(f**n/(c*(m + n*(S(2)*p + S(1)) + S(1))), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**p*(a*e*(m - n + S(1)) - c*d*x**n*(m + n*(S(2)*p + S(1)) + S(1))), x), x) + Simp(e*f**(n + S(-1))*(f*x)**(m - n + S(1))*(a + c*x**(S(2)*n))**(p + S(1))/(c*(m + n*(S(2)*p + S(1)) + S(1))), x)


def replacement1344(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(f**(-n)/(a*(m + S(1))), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**p*Simp(a*e*(m + S(1)) - b*d*(m + n*(p + S(1)) + S(1)) - c*d*x**n*(m + S(2)*n*(p + S(1)) + S(1)), x), x), x) + Simp(d*(f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(a*f*(m + S(1))), x)


def replacement1345(a, c, d, e, f, m, n, n2, p, x):
    return Dist(f**(-n)/(a*(m + S(1))), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**p*(a*e*(m + S(1)) - c*d*x**n*(m + S(2)*n*(p + S(1)) + S(1))), x), x) + Simp(d*(f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))/(a*f*(m + S(1))), x)


def With1346(a, b, c, d, e, f, m, n, n2, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(a*c, S(2))
    r = Rt(-b*c + S(2)*c*q, S(2))
    if Not(NegativeQ(-b*c + S(2)*c*q)):
        return True
    return False


def replacement1346(a, b, c, d, e, f, m, n, n2, x):

    q = Rt(a*c, S(2))
    r = Rt(-b*c + S(2)*c*q, S(2))
    return Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r - x**(n/2)*(c*d - e*q), x)/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r + x**(n/2)*(c*d - e*q), x)/(c*x**n + q + r*x**(n/2)), x), x)


def With1347(a, c, d, e, f, m, n, n2, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(a*c, S(2))
    r = Rt(S(2)*c*q, S(2))
    if Not(NegativeQ(S(2)*c*q)):
        return True
    return False


def replacement1347(a, c, d, e, f, m, n, n2, x):

    q = Rt(a*c, S(2))
    r = Rt(S(2)*c*q, S(2))
    return Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r - x**(n/2)*(c*d - e*q), x)/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*Simp(d*r + x**(n/2)*(c*d - e*q), x)/(c*x**n + q + r*x**(n/2)), x), x)


def With1348(a, b, c, d, e, f, m, x):
    r = Rt(c*(-b*e + S(2)*c*d)/e, S(2))
    return Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) - r*x), x), x) + Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) + r*x), x), x)


def With1349(a, c, d, e, f, m, x):
    r = Rt(S(2)*c**S(2)*d/e, S(2))
    return Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) - r*x), x), x) + Dist(e/S(2), Int((f*x)**m/(c*d/e + c*x**S(2) + r*x), x), x)


def With1350(a, b, c, d, e, f, m, n, n2, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(a*c, S(2))
    r = Rt(-b*c + S(2)*c*q, S(2))
    if Not(NegativeQ(-b*c + S(2)*c*q)):
        return True
    return False


def replacement1350(a, b, c, d, e, f, m, n, n2, x):

    q = Rt(a*c, S(2))
    r = Rt(-b*c + S(2)*c*q, S(2))
    return Dist(c/(2*q*r), Int((f*x)**m*(d*r - x**(n/2)*(c*d - e*q))/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*(d*r + x**(n/2)*(c*d - e*q))/(c*x**n + q + r*x**(n/2)), x), x)


def With1351(a, c, d, e, f, m, n, n2, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    q = Rt(a*c, S(2))
    r = Rt(S(2)*c*q, S(2))
    if Not(NegativeQ(S(2)*c*q)):
        return True
    return False


def replacement1351(a, c, d, e, f, m, n, n2, x):

    q = Rt(a*c, S(2))
    r = Rt(S(2)*c*q, S(2))
    return Dist(c/(2*q*r), Int((f*x)**m*(d*r - x**(n/2)*(c*d - e*q))/(c*x**n + q - r*x**(n/2)), x), x) + Dist(c/(2*q*r), Int((f*x)**m*(d*r + x**(n/2)*(c*d - e*q))/(c*x**n + q + r*x**(n/2)), x), x)


def With1352(a, b, c, d, e, f, m, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(e/S(2) - (-b*e + S(2)*c*d)/(S(2)*q), Int((f*x)**m/(b/S(2) + c*x**n + q/S(2)), x), x) + Dist(e/S(2) + (-b*e + S(2)*c*d)/(S(2)*q), Int((f*x)**m/(b/S(2) + c*x**n - q/S(2)), x), x)


def With1353(a, c, d, e, f, m, n, n2, x):
    q = Rt(-a*c, S(2))
    return Dist(-c*d/(S(2)*q) + e/S(2), Int((f*x)**m/(c*x**n + q), x), x) - Dist(c*d/(S(2)*q) + e/S(2), Int((f*x)**m/(-c*x**n + q), x), x)


def replacement1354(a, b, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1355(a, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)


def replacement1356(a, b, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m, (d + e*x**n)**q/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1357(a, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m, (d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)


def replacement1358(a, b, c, d, e, f, m, n, n2, q, x):
    return Dist(f**(S(2)*n)/c**S(2), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(-1))*(-b*e + c*d + c*e*x**n), x), x) - Dist(f**(S(2)*n)/c**S(2), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(-1))*Simp(a*(-b*e + c*d) + x**n*(a*c*e - b**S(2)*e + b*c*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1359(a, c, d, e, f, m, n, n2, q, x):
    return Dist(f**(S(2)*n)/c, Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x) - Dist(a*f**(S(2)*n)/c, Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q/(a + c*x**(S(2)*n)), x), x)


def replacement1360(a, b, c, d, e, f, m, n, n2, q, x):
    return -Dist(f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - x**n*(-b*e + c*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Dist(e*f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1)), x), x)


def replacement1361(a, c, d, e, f, m, n, n2, q, x):
    return -Dist(f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x) + Dist(e*f**n/c, Int((f*x)**(m - n)*(d + e*x**n)**(q + S(-1)), x), x)


def replacement1362(a, b, c, d, e, f, m, n, n2, q, x):
    return Dist(d/a, Int((f*x)**m*(d + e*x**n)**(q + S(-1)), x), x) - Dist(f**(-n)/a, Int((f*x)**(m + n)*(d + e*x**n)**(q + S(-1))*Simp(-a*e + b*d + c*d*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1363(a, c, d, e, f, m, n, n2, q, x):
    return Dist(d/a, Int((f*x)**m*(d + e*x**n)**(q + S(-1)), x), x) + Dist(f**(-n)/a, Int((f*x)**(m + n)*(d + e*x**n)**(q + S(-1))*Simp(a*e - c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x)


def replacement1364(a, b, c, d, e, f, m, n, n2, q, x):
    return -Dist(f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**(q + S(1))*Simp(a*d + x**n*(-a*e + b*d), x)/(a + b*x**n + c*x**(S(2)*n)), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x)


def replacement1365(a, c, d, e, f, m, n, n2, q, x):
    return -Dist(a*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(d + e*x**n)**q, x), x)


def replacement1366(a, b, c, d, e, f, m, n, n2, q, x):
    return Dist(f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**(q + S(1))*Simp(a*e + c*d*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x) - Dist(d*e*f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**q, x), x)


def replacement1367(a, c, d, e, f, m, n, n2, q, x):
    return Dist(f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**(q + S(1))*Simp(a*e + c*d*x**n, x)/(a + c*x**(S(2)*n)), x), x) - Dist(d*e*f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(d + e*x**n)**q, x), x)


def replacement1368(a, b, c, d, e, f, m, n, n2, q, x):
    return Dist(e**S(2)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**q, x), x) + Dist(S(1)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**(q + S(1))*Simp(-b*e + c*d - c*e*x**n, x)/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1369(a, c, d, e, f, m, n, n2, q, x):
    return Dist(c/(a*e**S(2) + c*d**S(2)), Int((f*x)**m*(d - e*x**n)*(d + e*x**n)**(q + S(1))/(a + c*x**(S(2)*n)), x), x) + Dist(e**S(2)/(a*e**S(2) + c*d**S(2)), Int((f*x)**m*(d + e*x**n)**q, x), x)


def replacement1370(a, b, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q, (f*x)**m/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1371(a, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((d + e*x**n)**q, (f*x)**m/(a + c*x**(S(2)*n)), x), x)


def replacement1372(a, b, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q, S(1)/(a + b*x**n + c*x**(S(2)*n)), x), x)


def replacement1373(a, c, d, e, f, m, n, n2, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q, S(1)/(a + c*x**(S(2)*n)), x), x)


def replacement1374(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(d**(S(-2)), Int((f*x)**m*(a*d + x**n*(-a*e + b*d))*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) + Dist(f**(-S(2)*n)*(a*e**S(2) - b*d*e + c*d**S(2))/d**S(2), Int((f*x)**(m + S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)


def replacement1375(a, c, d, e, f, m, n, n2, p, x):
    return Dist(a/d**S(2), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*(d - e*x**n), x), x) + Dist(f**(-S(2)*n)*(a*e**S(2) + c*d**S(2))/d**S(2), Int((f*x)**(m + S(2)*n)*(a + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)


def replacement1376(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(d*e), Int((f*x)**m*(a*e + c*d*x**n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1)), x), x) - Dist(f**(-n)*(a*e**S(2) - b*d*e + c*d**S(2))/(d*e), Int((f*x)**(m + n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)


def replacement1377(a, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(d*e), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(-1))*(a*e + c*d*x**n), x), x) - Dist(f**(-n)*(a*e**S(2) + c*d**S(2))/(d*e), Int((f*x)**(m + n)*(a + c*x**(S(2)*n))**(p + S(-1))/(d + e*x**n), x), x)


def replacement1378(a, b, c, d, e, f, m, n, n2, p, x):
    return -Dist(f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a*d + x**n*(-a*e + b*d))*(a + b*x**n + c*x**(S(2)*n))**p, x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)


def replacement1379(a, c, d, e, f, m, n, n2, p, x):
    return -Dist(a*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + c*x**(S(2)*n))**p*(d - e*x**n), x), x) + Dist(d**S(2)*f**(S(2)*n)/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - S(2)*n)*(a + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)


def replacement1380(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(a*e + c*d*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x) - Dist(d*e*f**n/(a*e**S(2) - b*d*e + c*d**S(2)), Int((f*x)**(m - n)*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)


def replacement1381(a, c, d, e, f, m, n, n2, p, x):
    return Dist(f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**p*(a*e + c*d*x**n), x), x) - Dist(d*e*f**n/(a*e**S(2) + c*d**S(2)), Int((f*x)**(m - n)*(a + c*x**(S(2)*n))**(p + S(1))/(d + e*x**n), x), x)


def replacement1382(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((a + b*x**n + c*x**(S(2)*n))**p, (f*x)**m*(d + e*x**n)**q, x), x)


def replacement1383(a, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((a + c*x**(S(2)*n))**p, (f*x)**m*(d + e*x**n)**q, x), x)


def replacement1384(a, b, c, d, e, m, n, n2, p, q, x):
    return -Subst(Int(x**(-m + S(-2))*(d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x)


def replacement1385(a, c, d, e, m, n, n2, p, q, x):
    return -Subst(Int(x**(-m + S(-2))*(a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q, x), x, S(1)/x)


def With1386(a, b, c, d, e, f, m, n, n2, p, q, x):
    g = Denominator(m)
    return -Dist(g/f, Subst(Int(x**(-g*(m + S(1)) + S(-1))*(d + e*f**(-n)*x**(-g*n))**q*(a + b*f**(-n)*x**(-g*n) + c*f**(-S(2)*n)*x**(-S(2)*g*n))**p, x), x, (f*x)**(-S(1)/g)), x)


def With1387(a, c, d, e, f, m, n, n2, p, q, x):
    g = Denominator(m)
    return -Dist(g/f, Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + c*f**(-S(2)*n)*x**(-S(2)*g*n))**p*(d + e*f**(-n)*x**(-g*n))**q, x), x, (f*x)**(-S(1)/g)), x)


def replacement1388(a, b, c, d, e, f, m, n, n2, p, q, x):
    return -Dist(f**IntPart(m)*(f*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(d + e*x**(-n))**q*(a + b*x**(-n) + c*x**(-S(2)*n))**p, x), x, S(1)/x), x)


def replacement1389(a, c, d, e, f, m, n, n2, p, q, x):
    return -Dist(f**IntPart(m)*(f*x)**FracPart(m)*(S(1)/x)**FracPart(m), Subst(Int(x**(-m + S(-2))*(a + c*x**(-S(2)*n))**p*(d + e*x**(-n))**q, x), x, S(1)/x), x)


def With1390(a, b, c, d, e, m, n, n2, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(d + e*x**(g*n))**q*(a + b*x**(g*n) + c*x**(S(2)*g*n))**p, x), x, x**(S(1)/g)), x)


def With1391(a, c, d, e, m, n, n2, p, q, x):
    g = Denominator(n)
    return Dist(g, Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + c*x**(S(2)*g*n))**p*(d + e*x**(g*n))**q, x), x, x**(S(1)/g)), x)


def replacement1392(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1393(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1394(a, b, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((d + e*x**(n/(m + S(1))))**q*(a + b*x**(n/(m + S(1))) + c*x**(S(2)*n/(m + S(1))))**p, x), x, x**(m + S(1))), x)


def replacement1395(a, c, d, e, m, n, n2, p, q, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + c*x**(S(2)*n/(m + S(1))))**p*(d + e*x**(n/(m + S(1))))**q, x), x, x**(m + S(1))), x)


def replacement1396(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1397(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def With1398(a, b, c, d, e, f, m, n, n2, q, x):
    r = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/r, Int((f*x)**m*(d + e*x**n)**q/(b + S(2)*c*x**n - r), x), x) - Dist(S(2)*c/r, Int((f*x)**m*(d + e*x**n)**q/(b + S(2)*c*x**n + r), x), x)


def With1399(a, c, d, e, f, m, n, n2, q, x):
    r = Rt(-a*c, S(2))
    return -Dist(c/(S(2)*r), Int((f*x)**m*(d + e*x**n)**q/(-c*x**n + r), x), x) - Dist(c/(S(2)*r), Int((f*x)**m*(d + e*x**n)**q/(c*x**n + r), x), x)


def replacement1400(a, b, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(a*n*(p + S(1))*(-S(4)*a*c + b**S(2))), Int((f*x)**m*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*Simp(-a*b*e*(m + S(1)) + c*x**n*(-S(2)*a*e + b*d)*(m + n*(S(2)*p + S(3)) + S(1)) + d*(-S(2)*a*c*(m + S(2)*n*(p + S(1)) + S(1)) + b**S(2)*(m + n*(p + S(1)) + S(1))), x), x), x) - Simp((f*x)**(m + S(1))*(a + b*x**n + c*x**(S(2)*n))**(p + S(1))*(-a*b*e + c*x**n*(-S(2)*a*e + b*d) + d*(-S(2)*a*c + b**S(2)))/(a*f*n*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1401(a, c, d, e, f, m, n, n2, p, x):
    return Dist(S(1)/(S(2)*a*n*(p + S(1))), Int((f*x)**m*(a + c*x**(S(2)*n))**(p + S(1))*Simp(d*(m + S(2)*n*(p + S(1)) + S(1)) + e*x**n*(m + n*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp((f*x)**(m + S(1))*(a + c*x**(S(2)*n))**(p + S(1))*(d + e*x**n)/(S(2)*a*f*n*(p + S(1))), x)


def replacement1402(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1403(a, c, d, e, f, m, n, n2, p, q, x):
    return Int(ExpandIntegrand((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1404(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(f**m, Int(ExpandIntegrand(x**m*(a + c*x**(S(2)*n))**p, (d/(d**S(2) - e**S(2)*x**(S(2)*n)) - e*x**n/(d**S(2) - e**S(2)*x**(S(2)*n)))**(-q), x), x), x)


def replacement1405(a, c, d, e, f, m, n, n2, p, q, x):
    return Dist(x**(-m)*(f*x)**m, Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x)


def replacement1406(a, b, c, d, e, f, m, n, n2, p, q, x):
    return Int((f*x)**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1407(a, c, d, e, f, m, n, n2, p, q, x):
    return Int((f*x)**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x)


def replacement1408(a, b, c, d, e, m, n, n2, p, q, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(d + e*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x, v), x)


def replacement1409(a, c, d, e, m, n, n2, p, q, u, v, x):
    return Dist(u**m*v**(-m)/Coefficient(v, x, S(1)), Subst(Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**n)**q, x), x, v), x)


def replacement1410(a, b, c, d, e, m, mn, n, n2, p, q, x):
    return Int(x**(m - n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1411(a, c, d, e, m, mn, n2, p, q, x):
    return Int(x**(m + mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x)


def replacement1412(a, b, c, d, e, m, mn, n, n2, p, q, x):
    return Int(x**(m + S(2)*n*p)*(d + e*x**(-n))**q*(a*x**(-S(2)*n) + b*x**(-n) + c)**p, x)


def replacement1413(a, c, d, e, m, mn, n2, p, q, x):
    return Int(x**(m - S(2)*mn*p)*(d + e*x**mn)**q*(a*x**(S(2)*mn) + c)**p, x)


def replacement1414(a, b, c, d, e, m, mn, n, n2, p, q, x):
    return Dist(x**(n*FracPart(q))*(d + e*x**(-n))**FracPart(q)*(d*x**n + e)**(-FracPart(q)), Int(x**(m - n*q)*(d*x**n + e)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1415(a, c, d, e, m, mn, n2, p, q, x):
    return Dist(x**(-mn*FracPart(q))*(d + e*x**mn)**FracPart(q)*(d*x**(-mn) + e)**(-FracPart(q)), Int(x**(m + mn*q)*(a + c*x**n2)**p*(d*x**(-mn) + e)**q, x), x)


def replacement1416(a, b, c, d, e, f, m, mn, n, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**mn)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1417(a, c, d, e, f, m, mn, n2, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(a + c*x**(S(2)*n))**p*(d + e*x**mn)**q, x), x)


def replacement1418(a, b, c, d, e, m, mn, n, p, q, x):
    return Int(x**(m - n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x)


def replacement1419(a, b, c, d, e, m, mn, n, p, q, x):
    return Dist(x**(n*FracPart(p))*(a + b*x**(-n) + c*x**n)**FracPart(p)*(a*x**n + b + c*x**(S(2)*n))**(-FracPart(p)), Int(x**(m - n*p)*(d + e*x**n)**q*(a*x**n + b + c*x**(S(2)*n))**p, x), x)


def replacement1420(a, b, c, d, e, f, m, mn, n, p, q, x):
    return Dist(f**IntPart(m)*x**(-FracPart(m))*(f*x)**FracPart(m), Int(x**m*(d + e*x**n)**q*(a + b*x**(-n) + c*x**n)**p, x), x)


def replacement1421(a, b, c, d1, d2, e1, e2, f, m, n, n2, non2, p, q, x):
    return Int((f*x)**m*(d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x)


def replacement1422(a, b, c, d1, d2, e1, e2, f, m, n, n2, non2, p, q, x):
    return Dist((d1 + e1*x**(n/S(2)))**FracPart(q)*(d2 + e2*x**(n/S(2)))**FracPart(q)*(d1*d2 + e1*e2*x**n)**(-FracPart(q)), Int((f*x)**m*(d1*d2 + e1*e2*x**n)**q*(a + b*x**n + c*x**(S(2)*n))**p, x), x)


def replacement1423(a, b, c, n, p, q, r, x):
    return Int((x**n*(a + b + c))**p, x)


def replacement1424(a, b, c, n, p, q, r, x):
    return Int(x**(p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)


def replacement1425(a, b, c, n, q, r, x):
    return Dist(x**(-q/S(2))*sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), Int(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)


def replacement1426(a, b, c, n, q, r, x):
    return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(-q/S(2))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)


def replacement1427(a, b, c, n, p, q, r, x):
    return Dist(p*(n - q)/(p*(S(2)*n - q) + S(1)), Int(x**q*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(p*(S(2)*n - q) + S(1)), x)


def replacement1428(a, b, c, n, p, q, r, x):
    return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(-q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(b*c*x**(n - q)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1)) + (n - q)*(p + S(1))*(-S(4)*a*c + b**S(2)) + (-S(2)*a*c + b**S(2))*(p*q + S(1))), x), x) - Simp(x**(S(1) - q)*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1429(a, b, c, n, p, q, r, x):
    return Dist(x**(-p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**(-p)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, Int(x**(p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x), x)


def replacement1430(a, b, c, n, p, q, r, x):
    return Int((a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x)


def replacement1431(a, b, c, n, p, q, r, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)


def replacement1432(a, b, c, m, n, p, q, r, x):
    return Int(x**m*(x**n*(a + b + c))**p, x)


def replacement1433(a, b, c, m, n, p, q, r, x):
    return Int(x**(m + p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)


def replacement1434(a, b, c, m, n, q, r, x):
    return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(m - q/S(2))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)


def replacement1435(a, b, c, m, n, q, r, x):
    return Simp(-S(2)*x**(n/S(2) + S(-1)/2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))), x)


def replacement1436(a, b, c, m, n, q, r, x):
    return Simp(x**(n/S(2) + S(-1)/2)*(S(4)*a + S(2)*b*x)/((-S(4)*a*c + b**S(2))*sqrt(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))), x)


def replacement1437(a, b, c, m, n, p, q, r, x):
    return -Dist(b/(S(2)*c), Int(x**(m + S(-1))*(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))**p, x), x) + Simp(x**(m - n)*(a*x**(n + S(-1)) + b*x**n + c*x**(n + S(1)))**(p + S(1))/(S(2)*c*(p + S(1))), x)


def replacement1438(a, b, c, m, n, p, q, r, x):
    return -Dist(p*(-S(4)*a*c + b**S(2))/(S(2)*c*(S(2)*p + S(1))), Int(x**(m + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m - n + q + S(1))*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(S(2)*c*(n - q)*(S(2)*p + S(1))), x)


def replacement1439(a, b, c, m, n, p, q, r, x):
    return Dist(p*(n - q)/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1))), Int(x**(m - n + S(2)*q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-a*b*(m - n + p*q + q + S(1)) + x**(n - q)*(S(2)*a*c*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1)) - b**S(2)*(m + p*q + (n - q)*(p + S(-1)) + S(1))), x), x), x) + Simp(x**(m - n + q + S(1))*(b*p*(n - q) + c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(-1)) + S(1))), x)


def replacement1440(a, b, c, m, n, p, q, r, x):
    return -Dist(p*(n - q)/(m + p*q + S(1)), Int(x**(m + n)*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(m + p*q + S(1)), x)


def replacement1441(a, b, c, m, n, p, q, r, x):
    return Dist(p*(n - q)/(m + p*(S(2)*n - q) + S(1)), Int(x**(m + q)*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1)), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/(m + p*(S(2)*n - q) + S(1)), x)


def replacement1442(a, b, c, m, n, p, q, r, x):
    return Dist((S(2)*a*c - b**S(2)*(p + S(2)))/(a*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) - Simp(x**(m - q + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1443(a, b, c, m, n, p, q, r, x):
    return Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - S(2)*n + q)*(S(2)*a*(m - S(2)*n + p*q + S(2)*q + S(1)) + b*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) - Simp(x**(m - S(2)*n + q + S(1))*(S(2)*a + b*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1444(a, b, c, m, n, p, q, r, x):
    return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*a*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + b**S(2)*(m + p*q + (n - q)*(p + S(1)) + S(1)) + b*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(m - q + S(1))*(-S(2)*a*c + b**S(2) + b*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1445(a, b, c, m, n, p, q, r, x):
    return -Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - n)*(b*(m - n + p*q + q + S(1)) + S(2)*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1)), x), x) + Simp(x**(m - n + S(1))*(b + S(2)*c*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def replacement1446(a, b, c, m, n, p, q, r, x):
    return -Dist(b/(S(2)*c), Int(x**(m - n + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - S(2)*n + q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*c*(n - q)*(p + S(1))), x)


def replacement1447(a, b, c, m, n, p, q, r, x):
    return -Dist(b/(S(2)*a), Int(x**(m + n - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) - Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*a*(n - q)*(p + S(1))), x)


def replacement1448(a, b, c, m, n, p, q, r, x):
    return -Dist(S(1)/(c*(m + p*q + S(2)*p*(n - q) + S(1))), Int(x**(m - S(2)*n + S(2)*q)*(a*(m - S(2)*n + p*q + S(2)*q + S(1)) + b*x**(n - q)*(m + p*q + (n - q)*(p + S(-1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - S(2)*n + q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + S(2)*p*(n - q) + S(1))), x)


def replacement1449(a, b, c, m, n, p, q, r, x):
    return -Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(b*(m + p*q + (n - q)*(p + S(1)) + S(1)) + c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x) + Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)


def replacement1450(a, b, c, m, n, p, q, r, x):
    return Dist(x**(-p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**(-p)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, Int(x**(m + p*q)*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x), x)


def replacement1451(a, b, c, m, n, p, q, r, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int(x**m*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)


def replacement1452(A, B, a, b, c, j, n, p, q, r, x):
    return Int(x**(p*q)*(A + B*x**(n - q))*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)


def replacement1453(A, B, a, b, c, j, n, q, r, x):
    return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(-q/S(2))*(A + B*x**(n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)


def replacement1454(A, B, a, b, c, j, n, p, q, r, x):
    return Dist(p*(n - q)/(c*(p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**q*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*(S(2)*A*a*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) - B*a*b*(p*q + S(1)) + x**(n - q)*(A*b*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*c*(p*(S(2)*n - q) + S(1)) - B*b**S(2)*(p*q + p*(n - q) + S(1)))), x), x) + Simp(x*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*(A*c*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*p*(n - q) + B*c*x**(n - q)*(p*(S(2)*n - q) + S(1)))/(c*(p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def With1455(A, B, a, c, j, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), NonzeroQ(p*(S(2)*n - q) + S(1)), NonzeroQ(p*q + (n - q)*(S(2)*p + S(1)) + S(1))):
        return True
    return False


def replacement1455(A, B, a, c, j, p, q, r, x):

    n = q + r
    return Dist(p*(n - q)/((p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**q*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*(S(2)*A*a*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*x**(n - q)*(p*(S(2)*n - q) + S(1))), x), x) + Simp(x*(A*(p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(p*(S(2)*n - q) + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((p*(S(2)*n - q) + S(1))*(p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def replacement1456(A, B, a, b, c, j, n, p, q, r, x):
    return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(-q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c*(p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + A*b**S(2)*(p*q + (n - q)*(p + S(1)) + S(1)) - B*a*b*(p*q + S(1)) + c*x**(n - q)*(A*b - S(2)*B*a)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(S(1) - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c + A*b**S(2) - B*a*b + c*x**(n - q)*(A*b - S(2)*B*a))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def With1457(A, B, a, c, j, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if ZeroQ(j - S(2)*n + q):
        return True
    return False


def replacement1457(A, B, a, c, j, p, q, r, x):

    n = q + r
    return Dist(S(1)/(S(2)*a**S(2)*c*(n - q)*(p + S(1))), Int(x**(-q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(A*a*c*(p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*c*x**(n - q)*(p*q + (n - q)*(S(2)*p + S(3)) + S(1))), x), x) - Simp(x**(S(1) - q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(A*a*c + B*a*c*x**(n - q))/(S(2)*a**S(2)*c*(n - q)*(p + S(1))), x)


def replacement1458(A, B, a, b, c, j, n, p, q, r, x):
    return Int((A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x)


def replacement1459(A, B, a, b, c, j, n, p, q, r, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)


def replacement1460(A, B, a, b, c, j, m, n, p, q, r, x):
    return Int(x**(m + p*q)*(A + B*x**(n - q))*(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))**p, x)


def replacement1461(A, B, a, b, c, j, m, n, p, q, r, x):
    return Dist(p*(n - q)/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + n)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-A*b*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*(m + p*q + S(1)) + x**(n - q)*(-S(2)*A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*(m + p*q + S(1))), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(1)))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def With1462(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), LessEqual(m + p*q, -n + q), Unequal(m + p*q + S(1), S(0)), Unequal(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1), S(0))):
        return True
    return False


def replacement1462(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return Dist(S(2)*p*(n - q)/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + n)*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*Simp(-A*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((m + p*q + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def replacement1463(A, B, a, b, c, j, m, n, p, q, r, x):
    return Dist(S(1)/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - n)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*Simp(x**(n - q)*(-S(2)*A*c + B*b)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + (-A*b + S(2)*B*a)*(m - n + p*q + q + S(1)), x), x), x) + Simp(x**(m - n + S(1))*(A*b - S(2)*B*a - x**(n - q)*(-S(2)*A*c + B*b))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/((n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def With1464(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), Greater(m + p*q, n - q + S(-1))):
        return True
    return False


def replacement1464(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return -Dist(S(1)/(S(2)*a*c*(n - q)*(p + S(1))), Int(x**(m - n)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*Simp(-A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m - n + p*q + q + S(1)), x), x), x) + Simp(x**(m - n + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*(-A*c*x**(n - q) + B*a)/(S(2)*a*c*(n - q)*(p + S(1))), x)


def replacement1465(A, B, a, b, c, j, m, n, p, q, r, x):
    return Dist(p*(n - q)/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(-1))*Simp(S(2)*A*a*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) - B*a*b*(m + p*q + S(1)) + x**(n - q)*(A*b*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*c*(m + p*q + S(2)*p*(n - q) + S(1)) - B*b**S(2)*(m + p*q + p*(n - q) + S(1))), x), x), x) + Simp(x**(m + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*(A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*p*(n - q) + B*c*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)))/(c*(m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def With1466(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), Greater(m + p*q, -n + q), Unequal(m + p*q + S(2)*p*(n - q) + S(1), S(0)), Unequal(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1), S(0)), Unequal(m + S(1), n)):
        return True
    return False


def replacement1466(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return Dist(p*(n - q)/((m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m + q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(-1))*Simp(S(2)*A*a*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + S(2)*B*a*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)), x), x), x) + Simp(x**(m + S(1))*(A*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*x**(n - q)*(m + p*q + S(2)*p*(n - q) + S(1)))*(a*x**q + c*x**(S(2)*n - q))**p/((m + p*(S(2)*n - q) + S(1))*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def replacement1467(A, B, a, b, c, j, m, n, p, q, r, x):
    return Dist(S(1)/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), Int(x**(m - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*Simp(-S(2)*A*a*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + A*b**S(2)*(m + p*q + (n - q)*(p + S(1)) + S(1)) - B*a*b*(m + p*q + S(1)) + c*x**(n - q)*(A*b - S(2)*B*a)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))*(-S(2)*A*a*c + A*b**S(2) - B*a*b + c*x**(n - q)*(A*b - S(2)*B*a))/(a*(n - q)*(p + S(1))*(-S(4)*a*c + b**S(2))), x)


def With1468(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), Less(m + p*q, n - q + S(-1))):
        return True
    return False


def replacement1468(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return Dist(S(1)/(S(2)*a*c*(n - q)*(p + S(1))), Int(x**(m - q)*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))*Simp(A*c*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(3)) + S(1)), x), x), x) - Simp(x**(m - q + S(1))*(A*c + B*c*x**(n - q))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(S(2)*a*c*(n - q)*(p + S(1))), x)


def replacement1469(A, B, a, b, c, j, m, n, p, q, r, x):
    return -Dist(S(1)/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m - n + q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*Simp(B*a*(m - n + p*q + q + S(1)) + x**(n - q)*(-A*c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*b*(m + p*q + p*(n - q) + S(1))), x), x), x) + Simp(B*x**(m - n + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def With1470(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), GreaterEqual(m + p*q, n - q + S(-1)), Unequal(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1), S(0))):
        return True
    return False


def replacement1470(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return -Dist(S(1)/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), Int(x**(m - n + q)*(a*x**q + c*x**(S(2)*n - q))**p*Simp(-A*c*x**(n - q)*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1)) + B*a*(m - n + p*q + q + S(1)), x), x), x) + Simp(B*x**(m - n + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(c*(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1))), x)


def replacement1471(A, B, a, b, c, j, m, n, p, q, r, x):
    return Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p*Simp(-A*b*(m + p*q + (n - q)*(p + S(1)) + S(1)) - A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(A*x**(m - q + S(1))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)


def With1472(A, B, a, c, j, m, p, q, r, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    n = q + r
    if And(ZeroQ(j - S(2)*n + q), PositiveIntegerQ(n), Or(Inequality(S(-1), LessEqual, p, Less, S(0)), Equal(m + p*q + (n - q)*(S(2)*p + S(1)) + S(1), S(0))), LessEqual(m + p*q, -n + q), Unequal(m + p*q + S(1), S(0))):
        return True
    return False


def replacement1472(A, B, a, c, j, m, p, q, r, x):

    n = q + r
    return Dist(S(1)/(a*(m + p*q + S(1))), Int(x**(m + n - q)*(a*x**q + c*x**(S(2)*n - q))**p*Simp(-A*c*x**(n - q)*(m + p*q + S(2)*(n - q)*(p + S(1)) + S(1)) + B*a*(m + p*q + S(1)), x), x), x) + Simp(A*x**(m - q + S(1))*(a*x**q + c*x**(S(2)*n - q))**(p + S(1))/(a*(m + p*q + S(1))), x)


def replacement1473(A, B, a, b, c, j, m, n, q, r, x):
    return Dist(x**(q/S(2))*sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q))/sqrt(a*x**q + b*x**n + c*x**(S(2)*n - q)), Int(x**(m - q/S(2))*(A + B*x**(n - q))/sqrt(a + b*x**(n - q) + c*x**(S(2)*n - S(2)*q)), x), x)


def replacement1474(A, B, a, b, c, j, k, m, n, p, q, x):
    return Dist(x**(-j*p)*(a + b*x**(-j + k) + c*x**(-S(2)*j + S(2)*k))**(-p)*(a*x**j + b*x**k + c*x**n)**p, Int(x**(j*p + m)*(A + B*x**(-j + k))*(a + b*x**(-j + k) + c*x**(-S(2)*j + S(2)*k))**p, x), x)


def replacement1475(A, B, a, b, c, j, m, n, p, q, r, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int(x**m*(A + B*x**(n - q))*(a*x**q + b*x**n + c*x**(S(2)*n - q))**p, x), x, u), x)
