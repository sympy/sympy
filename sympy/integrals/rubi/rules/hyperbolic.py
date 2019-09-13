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

def hyperbolic(rubi):
    from sympy.integrals.rubi.constraints import cons33, cons170, cons8, cons29, cons50, cons127, cons96, cons1118, cons178, cons1490, cons19, cons89, cons167, cons3, cons95, cons168, cons87, cons1491, cons247, cons249, cons91, cons1444, cons150, cons1861, cons2, cons1492, cons1441, cons810, cons1493, cons1494, cons1456, cons1442, cons64, cons1269, cons1495, cons812, cons813, cons4, cons1362, cons130, cons40, cons139, cons746, cons65, cons1496, cons198, cons1497, cons5, cons55, cons13, cons598, cons1498, cons20, cons1499, cons378, cons148, cons491, cons1500, cons70, cons71, cons825, cons826, cons1501, cons1503, cons1257, cons1504, cons58, cons152, cons1505, cons1506, cons685, cons369, cons1507, cons358, cons68, cons856, cons25, cons1508, cons56, cons14, cons820, cons1133, cons1134, cons1135, cons1509, cons821, cons530, cons1267, cons1512, cons21, cons1573, cons1574, cons1575, cons1576, cons1577, cons1578, cons1579, cons1580, cons1581, cons1045, cons1582, cons1646, cons1738, cons1647, cons586, cons466, cons1685, cons1410, cons1686, cons1687, cons1862, cons1863, cons1690, cons814, cons815, cons557, cons1864, cons1865, cons27, cons1693, cons1866, cons1101, cons1867, cons1868, cons1397, cons1869, cons1695, cons1870, cons965, cons1871, cons1872, cons210, cons1702, cons1013, cons1553, cons1703, cons1704, cons211, cons226, cons1701, cons1873, cons1705, cons1706, cons1874, cons1708, cons1709, cons1875, cons1876, cons1877, cons1878, cons1879, cons1880, cons1881, cons1882, cons1883, cons1884, cons1885, cons1886, cons1887, cons1888, cons1722, cons1889, cons1890, cons1891, cons1725, cons1892, cons165, cons340, cons164, cons629, cons73, cons1727, cons1728, cons1729, cons90, cons1730, cons1458, cons465, cons1731, cons1480, cons1732, cons1893, cons36, cons37, cons1476, cons1483, cons1735

    pattern5646 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule5646 = ReplacementRule(pattern5646, replacement5646)
    pattern5647 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule5647 = ReplacementRule(pattern5647, replacement5647)
    pattern5648 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons96)
    rule5648 = ReplacementRule(pattern5648, replacement5648)
    pattern5649 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons96)
    rule5649 = ReplacementRule(pattern5649, replacement5649)
    pattern5650 = Pattern(Integral(sinh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule5650 = ReplacementRule(pattern5650, replacement5650)
    pattern5651 = Pattern(Integral(cosh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule5651 = ReplacementRule(pattern5651, replacement5651)
    pattern5652 = Pattern(Integral(sinh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule5652 = ReplacementRule(pattern5652, replacement5652)
    pattern5653 = Pattern(Integral(cosh(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule5653 = ReplacementRule(pattern5653, replacement5653)
    pattern5654 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons19, cons1490)
    rule5654 = ReplacementRule(pattern5654, replacement5654)
    pattern5655 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons19, cons1490)
    rule5655 = ReplacementRule(pattern5655, replacement5655)
    pattern5656 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167)
    rule5656 = ReplacementRule(pattern5656, replacement5656)
    pattern5657 = Pattern(Integral((WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167)
    rule5657 = ReplacementRule(pattern5657, replacement5657)
    pattern5658 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons168)
    rule5658 = ReplacementRule(pattern5658, replacement5658)
    pattern5659 = Pattern(Integral((WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons168)
    rule5659 = ReplacementRule(pattern5659, replacement5659)
    pattern5660 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons1491)
    rule5660 = ReplacementRule(pattern5660, replacement5660)
    pattern5661 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cosh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons1491)
    rule5661 = ReplacementRule(pattern5661, replacement5661)
    pattern5662 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sinh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons33, cons247)
    rule5662 = ReplacementRule(pattern5662, replacement5662)
    pattern5663 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cosh(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons33, cons247)
    rule5663 = ReplacementRule(pattern5663, replacement5663)
    pattern5664 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons249)
    rule5664 = ReplacementRule(pattern5664, replacement5664)
    pattern5665 = Pattern(Integral((WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons249)
    rule5665 = ReplacementRule(pattern5665, replacement5665)
    pattern5666 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91, cons1444)
    rule5666 = ReplacementRule(pattern5666, replacement5666)
    pattern5667 = Pattern(Integral((WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91, cons1444)
    rule5667 = ReplacementRule(pattern5667, replacement5667)
    pattern5668 = Pattern(Integral((WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons1444, cons168)
    rule5668 = ReplacementRule(pattern5668, replacement5668)
    pattern5669 = Pattern(Integral((WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons1444, cons168)
    rule5669 = ReplacementRule(pattern5669, replacement5669)
    pattern5670 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons150, cons1861)
    rule5670 = ReplacementRule(pattern5670, replacement5670)
    pattern5671 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons150, cons1492)
    rule5671 = ReplacementRule(pattern5671, replacement5671)
    pattern5672 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons87)
    rule5672 = ReplacementRule(pattern5672, replacement5672)
    pattern5673 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons810, cons1493)
    rule5673 = ReplacementRule(pattern5673, replacement5673)
    pattern5674 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1494, cons87)
    rule5674 = ReplacementRule(pattern5674, replacement5674)
    pattern5675 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1456, cons87)
    rule5675 = ReplacementRule(pattern5675, replacement5675)
    pattern5676 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1494, cons810, cons1493)
    rule5676 = ReplacementRule(pattern5676, replacement5676)
    pattern5677 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1456, cons810, cons1493)
    rule5677 = ReplacementRule(pattern5677, replacement5677)
    pattern5678 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons64)
    rule5678 = ReplacementRule(pattern5678, replacement5678)
    pattern5679 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule5679 = ReplacementRule(pattern5679, replacement5679)
    pattern5680 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons64)
    rule5680 = ReplacementRule(pattern5680, replacement5680)
    pattern5681 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule5681 = ReplacementRule(pattern5681, replacement5681)
    pattern5682 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons1495, cons64)
    rule5682 = ReplacementRule(pattern5682, replacement5682)
    pattern5683 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons1495, cons64)
    rule5683 = ReplacementRule(pattern5683, replacement5683)
    pattern5684 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule5684 = ReplacementRule(pattern5684, replacement5684)
    pattern5685 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule5685 = ReplacementRule(pattern5685, replacement5685)
    pattern5686 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5686 = ReplacementRule(pattern5686, replacement5686)
    pattern5687 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5687 = ReplacementRule(pattern5687, replacement5687)
    pattern5688 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule5688 = ReplacementRule(pattern5688, replacement5688)
    pattern5689 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule5689 = ReplacementRule(pattern5689, replacement5689)
    pattern5690 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons139, cons746)
    rule5690 = ReplacementRule(pattern5690, replacement5690)
    pattern5691 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons139, cons746)
    rule5691 = ReplacementRule(pattern5691, replacement5691)
    pattern5692 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons150, cons1496)
    rule5692 = ReplacementRule(pattern5692, replacement5692)
    pattern5693 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons150, cons1496)
    rule5693 = ReplacementRule(pattern5693, replacement5693)
    pattern5694 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons198)
    rule5694 = ReplacementRule(pattern5694, replacement5694)
    pattern5695 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons198)
    rule5695 = ReplacementRule(pattern5695, replacement5695)
    pattern5696 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5696 = ReplacementRule(pattern5696, replacement5696)
    pattern5697 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5697 = ReplacementRule(pattern5697, replacement5697)
    pattern5698 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule5698 = ReplacementRule(pattern5698, replacement5698)
    pattern5699 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule5699 = ReplacementRule(pattern5699, replacement5699)
    pattern5700 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons55, cons13, cons139, cons598)
    rule5700 = ReplacementRule(pattern5700, replacement5700)
    pattern5701 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons55, cons13, cons139, cons598)
    rule5701 = ReplacementRule(pattern5701, replacement5701)
    pattern5702 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons33, cons139, cons1498)
    rule5702 = ReplacementRule(pattern5702, replacement5702)
    pattern5703 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons33, cons139, cons1498)
    rule5703 = ReplacementRule(pattern5703, replacement5703)
    pattern5704 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons20, cons150, cons1496)
    rule5704 = ReplacementRule(pattern5704, replacement5704)
    pattern5705 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons20, cons150, cons1496)
    rule5705 = ReplacementRule(pattern5705, replacement5705)
    pattern5706 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons65, cons198)
    rule5706 = ReplacementRule(pattern5706, replacement5706)
    pattern5707 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons65, cons198)
    rule5707 = ReplacementRule(pattern5707, replacement5707)
    pattern5708 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5708 = ReplacementRule(pattern5708, replacement5708)
    pattern5709 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5709 = ReplacementRule(pattern5709, replacement5709)
    pattern5710 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons87, cons167)
    rule5710 = ReplacementRule(pattern5710, replacement5710)
    pattern5711 = Pattern(Integral(cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons87, cons167)
    rule5711 = ReplacementRule(pattern5711, replacement5711)
    pattern5712 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons378, cons167, cons148)
    rule5712 = ReplacementRule(pattern5712, replacement5712)
    pattern5713 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons378, cons167, cons148)
    rule5713 = ReplacementRule(pattern5713, replacement5713)
    pattern5714 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198)
    rule5714 = ReplacementRule(pattern5714, replacement5714)
    pattern5715 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198)
    rule5715 = ReplacementRule(pattern5715, replacement5715)

    pattern5716 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons491)
    rule5716 = ReplacementRule(pattern5716, With5716)

    pattern5717 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons491)
    rule5717 = ReplacementRule(pattern5717, With5717)
    pattern5718 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons4, cons1500)
    rule5718 = ReplacementRule(pattern5718, replacement5718)
    pattern5719 = Pattern(Integral(cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons4, cons1500)
    rule5719 = ReplacementRule(pattern5719, replacement5719)
    pattern5720 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule5720 = ReplacementRule(pattern5720, replacement5720)
    pattern5721 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule5721 = ReplacementRule(pattern5721, replacement5721)
    pattern5722 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons40, cons70, cons71)
    rule5722 = ReplacementRule(pattern5722, replacement5722)
    pattern5723 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons40, cons70, cons71)
    rule5723 = ReplacementRule(pattern5723, replacement5723)
    pattern5724 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70)
    rule5724 = ReplacementRule(pattern5724, replacement5724)
    pattern5725 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70)
    rule5725 = ReplacementRule(pattern5725, replacement5725)
    pattern5726 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sinh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5726 = ReplacementRule(pattern5726, replacement5726)
    pattern5727 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cosh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5727 = ReplacementRule(pattern5727, replacement5727)
    pattern5728 = Pattern(Integral(sinh(x_**n_*WC('d', S(1)))/x_, x_), cons29, cons4, cons1501)
    rule5728 = ReplacementRule(pattern5728, replacement5728)
    pattern5729 = Pattern(Integral(cosh(x_**n_*WC('d', S(1)))/x_, x_), cons29, cons4, cons1501)
    rule5729 = ReplacementRule(pattern5729, replacement5729)
    pattern5730 = Pattern(Integral(sinh(c_ + x_**n_*WC('d', S(1)))/x_, x_), cons8, cons29, cons4, cons1500)
    rule5730 = ReplacementRule(pattern5730, replacement5730)
    pattern5731 = Pattern(Integral(cosh(c_ + x_**n_*WC('d', S(1)))/x_, x_), cons8, cons29, cons4, cons1500)
    rule5731 = ReplacementRule(pattern5731, replacement5731)

    pattern5732 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, CustomConstraint(With5732))
    rule5732 = ReplacementRule(pattern5732, replacement5732)

    pattern5733 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, CustomConstraint(With5733))
    rule5733 = ReplacementRule(pattern5733, replacement5733)

    pattern5734 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, CustomConstraint(With5734))
    rule5734 = ReplacementRule(pattern5734, replacement5734)

    pattern5735 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, CustomConstraint(With5735))
    rule5735 = ReplacementRule(pattern5735, replacement5735)
    pattern5736 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons1503)
    rule5736 = ReplacementRule(pattern5736, replacement5736)
    pattern5737 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons1503)
    rule5737 = ReplacementRule(pattern5737, replacement5737)
    pattern5738 = Pattern(Integral((x_*WC('e', S(1)))**m_*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons96)
    rule5738 = ReplacementRule(pattern5738, replacement5738)
    pattern5739 = Pattern(Integral((x_*WC('e', S(1)))**m_*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons96)
    rule5739 = ReplacementRule(pattern5739, replacement5739)
    pattern5740 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons150)
    rule5740 = ReplacementRule(pattern5740, replacement5740)
    pattern5741 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons150)
    rule5741 = ReplacementRule(pattern5741, replacement5741)
    pattern5742 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons378, cons1257, cons148, cons1504)
    rule5742 = ReplacementRule(pattern5742, replacement5742)
    pattern5743 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons378, cons1257, cons148, cons1504)
    rule5743 = ReplacementRule(pattern5743, replacement5743)
    pattern5744 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons148)
    rule5744 = ReplacementRule(pattern5744, replacement5744)
    pattern5745 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons148)
    rule5745 = ReplacementRule(pattern5745, replacement5745)
    pattern5746 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1505)
    rule5746 = ReplacementRule(pattern5746, replacement5746)
    pattern5747 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1505)
    rule5747 = ReplacementRule(pattern5747, replacement5747)
    pattern5748 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1506, cons685)
    rule5748 = ReplacementRule(pattern5748, replacement5748)
    pattern5749 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1506, cons685)
    rule5749 = ReplacementRule(pattern5749, replacement5749)

    pattern5750 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150, cons369)
    rule5750 = ReplacementRule(pattern5750, With5750)

    pattern5751 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150, cons369)
    rule5751 = ReplacementRule(pattern5751, With5751)
    pattern5752 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons148)
    rule5752 = ReplacementRule(pattern5752, replacement5752)
    pattern5753 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons148)
    rule5753 = ReplacementRule(pattern5753, replacement5753)
    pattern5754 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons139, cons1507)
    rule5754 = ReplacementRule(pattern5754, replacement5754)
    pattern5755 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons139, cons1507)
    rule5755 = ReplacementRule(pattern5755, replacement5755)
    pattern5756 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons139, cons1507, cons1505)
    rule5756 = ReplacementRule(pattern5756, replacement5756)
    pattern5757 = Pattern(Integral(x_**WC('m', S(1))*cosh(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons139, cons1507, cons1505)
    rule5757 = ReplacementRule(pattern5757, replacement5757)
    pattern5758 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198, cons20)
    rule5758 = ReplacementRule(pattern5758, replacement5758)
    pattern5759 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198, cons20)
    rule5759 = ReplacementRule(pattern5759, replacement5759)

    pattern5760 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons198, cons369)
    rule5760 = ReplacementRule(pattern5760, With5760)

    pattern5761 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons198, cons369)
    rule5761 = ReplacementRule(pattern5761, With5761)
    pattern5762 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons198, cons358)
    rule5762 = ReplacementRule(pattern5762, replacement5762)
    pattern5763 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons198, cons358)
    rule5763 = ReplacementRule(pattern5763, replacement5763)

    pattern5764 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons40, cons491)
    rule5764 = ReplacementRule(pattern5764, With5764)

    pattern5765 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons40, cons491)
    rule5765 = ReplacementRule(pattern5765, With5765)
    pattern5766 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons491)
    rule5766 = ReplacementRule(pattern5766, replacement5766)
    pattern5767 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons491)
    rule5767 = ReplacementRule(pattern5767, replacement5767)
    pattern5768 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, cons68, cons856, cons25)
    rule5768 = ReplacementRule(pattern5768, replacement5768)
    pattern5769 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, cons68, cons856, cons25)
    rule5769 = ReplacementRule(pattern5769, replacement5769)
    pattern5770 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons68, cons856, cons25)
    rule5770 = ReplacementRule(pattern5770, replacement5770)
    pattern5771 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons68, cons856, cons25)
    rule5771 = ReplacementRule(pattern5771, replacement5771)
    pattern5772 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons4, cons1508)
    rule5772 = ReplacementRule(pattern5772, replacement5772)
    pattern5773 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons4, cons1508)
    rule5773 = ReplacementRule(pattern5773, replacement5773)
    pattern5774 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule5774 = ReplacementRule(pattern5774, replacement5774)
    pattern5775 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule5775 = ReplacementRule(pattern5775, replacement5775)
    pattern5776 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71, cons20)
    rule5776 = ReplacementRule(pattern5776, replacement5776)
    pattern5777 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71, cons20)
    rule5777 = ReplacementRule(pattern5777, replacement5777)
    pattern5778 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons70)
    rule5778 = ReplacementRule(pattern5778, replacement5778)
    pattern5779 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons70)
    rule5779 = ReplacementRule(pattern5779, replacement5779)
    pattern5780 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sinh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5780 = ReplacementRule(pattern5780, replacement5780)
    pattern5781 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cosh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5781 = ReplacementRule(pattern5781, replacement5781)
    pattern5782 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons4, cons5, cons55, cons56)
    rule5782 = ReplacementRule(pattern5782, replacement5782)
    pattern5783 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons55, cons56)
    rule5783 = ReplacementRule(pattern5783, replacement5783)
    pattern5784 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons95, cons1503, cons56)
    rule5784 = ReplacementRule(pattern5784, replacement5784)
    pattern5785 = Pattern(Integral(x_**WC('m', S(1))*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons5, cons95, cons1503, cons56)
    rule5785 = ReplacementRule(pattern5785, replacement5785)
    pattern5786 = Pattern(Integral(sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons14)
    rule5786 = ReplacementRule(pattern5786, replacement5786)
    pattern5787 = Pattern(Integral(cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons14)
    rule5787 = ReplacementRule(pattern5787, replacement5787)
    pattern5788 = Pattern(Integral(sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons87, cons167)
    rule5788 = ReplacementRule(pattern5788, replacement5788)
    pattern5789 = Pattern(Integral(cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons87, cons167)
    rule5789 = ReplacementRule(pattern5789, replacement5789)
    pattern5790 = Pattern(Integral(sinh(v_)**WC('n', S(1)), x_), cons150, cons820, cons1133)
    rule5790 = ReplacementRule(pattern5790, replacement5790)
    pattern5791 = Pattern(Integral(cosh(v_)**WC('n', S(1)), x_), cons150, cons820, cons1133)
    rule5791 = ReplacementRule(pattern5791, replacement5791)
    pattern5792 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1134)
    rule5792 = ReplacementRule(pattern5792, replacement5792)
    pattern5793 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1134)
    rule5793 = ReplacementRule(pattern5793, replacement5793)
    pattern5794 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1135)
    rule5794 = ReplacementRule(pattern5794, replacement5794)
    pattern5795 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1135)
    rule5795 = ReplacementRule(pattern5795, replacement5795)
    pattern5796 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1134)
    rule5796 = ReplacementRule(pattern5796, replacement5796)
    pattern5797 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1134)
    rule5797 = ReplacementRule(pattern5797, replacement5797)
    pattern5798 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1135)
    rule5798 = ReplacementRule(pattern5798, replacement5798)
    pattern5799 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1135)
    rule5799 = ReplacementRule(pattern5799, replacement5799)
    pattern5800 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1134)
    rule5800 = ReplacementRule(pattern5800, replacement5800)
    pattern5801 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1134)
    rule5801 = ReplacementRule(pattern5801, replacement5801)
    pattern5802 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1135)
    rule5802 = ReplacementRule(pattern5802, replacement5802)
    pattern5803 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1135)
    rule5803 = ReplacementRule(pattern5803, replacement5803)
    pattern5804 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule5804 = ReplacementRule(pattern5804, replacement5804)
    pattern5805 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule5805 = ReplacementRule(pattern5805, replacement5805)
    pattern5806 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sinh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons87, cons167)
    rule5806 = ReplacementRule(pattern5806, replacement5806)
    pattern5807 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cosh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons87, cons167)
    rule5807 = ReplacementRule(pattern5807, replacement5807)
    pattern5808 = Pattern(Integral(u_**WC('m', S(1))*sinh(v_)**WC('n', S(1)), x_), cons19, cons150, cons70, cons820, cons821)
    rule5808 = ReplacementRule(pattern5808, replacement5808)
    pattern5809 = Pattern(Integral(u_**WC('m', S(1))*cosh(v_)**WC('n', S(1)), x_), cons19, cons150, cons70, cons820, cons821)
    rule5809 = ReplacementRule(pattern5809, replacement5809)
    pattern5810 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5810 = ReplacementRule(pattern5810, replacement5810)
    pattern5811 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5811 = ReplacementRule(pattern5811, replacement5811)
    pattern5812 = Pattern(Integral((WC('c', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons95, cons167, cons170)
    rule5812 = ReplacementRule(pattern5812, replacement5812)
    pattern5813 = Pattern(Integral((WC('c', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons95, cons167, cons170)
    rule5813 = ReplacementRule(pattern5813, replacement5813)
    pattern5814 = Pattern(Integral((WC('c', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons95, cons91, cons170)
    rule5814 = ReplacementRule(pattern5814, replacement5814)
    pattern5815 = Pattern(Integral((WC('c', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons50, cons127, cons95, cons91, cons170)
    rule5815 = ReplacementRule(pattern5815, replacement5815)
    pattern5816 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule5816 = ReplacementRule(pattern5816, replacement5816)
    pattern5817 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule5817 = ReplacementRule(pattern5817, replacement5817)
    pattern5818 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons170)
    rule5818 = ReplacementRule(pattern5818, replacement5818)
    pattern5819 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons170)
    rule5819 = ReplacementRule(pattern5819, replacement5819)
    pattern5820 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267)
    rule5820 = ReplacementRule(pattern5820, replacement5820)
    pattern5821 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267)
    rule5821 = ReplacementRule(pattern5821, replacement5821)
    pattern5822 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons96, cons1512)
    rule5822 = ReplacementRule(pattern5822, replacement5822)
    pattern5823 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons96, cons1512)
    rule5823 = ReplacementRule(pattern5823, replacement5823)
    pattern5824 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267)
    rule5824 = ReplacementRule(pattern5824, replacement5824)
    pattern5825 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267)
    rule5825 = ReplacementRule(pattern5825, replacement5825)
    pattern5826 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons21)
    rule5826 = ReplacementRule(pattern5826, replacement5826)
    pattern5827 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons21)
    rule5827 = ReplacementRule(pattern5827, replacement5827)
    pattern5828 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons1573)
    rule5828 = ReplacementRule(pattern5828, replacement5828)
    pattern5829 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons1573)
    rule5829 = ReplacementRule(pattern5829, replacement5829)
    pattern5830 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons198)
    rule5830 = ReplacementRule(pattern5830, replacement5830)
    pattern5831 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons198)
    rule5831 = ReplacementRule(pattern5831, replacement5831)

    pattern5832 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons1574, cons33, cons170)
    rule5832 = ReplacementRule(pattern5832, With5832)

    pattern5833 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons1574, cons33, cons170)
    rule5833 = ReplacementRule(pattern5833, With5833)
    pattern5834 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule5834 = ReplacementRule(pattern5834, replacement5834)
    pattern5835 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule5835 = ReplacementRule(pattern5835, replacement5835)
    pattern5836 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269)
    rule5836 = ReplacementRule(pattern5836, replacement5836)
    pattern5837 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269)
    rule5837 = ReplacementRule(pattern5837, replacement5837)
    pattern5838 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons198, cons64)
    rule5838 = ReplacementRule(pattern5838, replacement5838)
    pattern5839 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons198, cons64)
    rule5839 = ReplacementRule(pattern5839, replacement5839)
    pattern5840 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule5840 = ReplacementRule(pattern5840, replacement5840)
    pattern5841 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule5841 = ReplacementRule(pattern5841, replacement5841)
    pattern5842 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5842 = ReplacementRule(pattern5842, replacement5842)
    pattern5843 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5843 = ReplacementRule(pattern5843, replacement5843)
    pattern5844 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule5844 = ReplacementRule(pattern5844, replacement5844)
    pattern5845 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule5845 = ReplacementRule(pattern5845, replacement5845)
    pattern5846 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5846 = ReplacementRule(pattern5846, replacement5846)
    pattern5847 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5847 = ReplacementRule(pattern5847, replacement5847)
    pattern5848 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule5848 = ReplacementRule(pattern5848, replacement5848)
    pattern5849 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tanh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule5849 = ReplacementRule(pattern5849, replacement5849)
    pattern5850 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tanh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5850 = ReplacementRule(pattern5850, replacement5850)
    pattern5851 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tanh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5851 = ReplacementRule(pattern5851, replacement5851)
    pattern5852 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule5852 = ReplacementRule(pattern5852, replacement5852)
    pattern5853 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule5853 = ReplacementRule(pattern5853, replacement5853)
    pattern5854 = Pattern(Integral(x_**WC('m', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons8, cons29, cons19, cons4, cons1577)
    rule5854 = ReplacementRule(pattern5854, replacement5854)
    pattern5855 = Pattern(Integral(x_**WC('m', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons8, cons29, cons19, cons4, cons1577)
    rule5855 = ReplacementRule(pattern5855, replacement5855)
    pattern5856 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule5856 = ReplacementRule(pattern5856, replacement5856)
    pattern5857 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule5857 = ReplacementRule(pattern5857, replacement5857)
    pattern5858 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5858 = ReplacementRule(pattern5858, replacement5858)
    pattern5859 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5859 = ReplacementRule(pattern5859, replacement5859)
    pattern5860 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tanh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5860 = ReplacementRule(pattern5860, replacement5860)
    pattern5861 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tanh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5861 = ReplacementRule(pattern5861, replacement5861)
    pattern5862 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*tanh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1580)
    rule5862 = ReplacementRule(pattern5862, replacement5862)
    pattern5863 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*(S(1)/tanh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('q', S(1)), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1580)
    rule5863 = ReplacementRule(pattern5863, replacement5863)
    pattern5864 = Pattern(Integral(tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5864 = ReplacementRule(pattern5864, replacement5864)
    pattern5865 = Pattern(Integral((S(1)/tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5865 = ReplacementRule(pattern5865, replacement5865)
    pattern5866 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5866 = ReplacementRule(pattern5866, replacement5866)
    pattern5867 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5867 = ReplacementRule(pattern5867, replacement5867)
    pattern5868 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule5868 = ReplacementRule(pattern5868, replacement5868)
    pattern5869 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(S(1)/tanh(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule5869 = ReplacementRule(pattern5869, replacement5869)
    pattern5870 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule5870 = ReplacementRule(pattern5870, replacement5870)
    pattern5871 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule5871 = ReplacementRule(pattern5871, replacement5871)
    pattern5872 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/cosh(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons33, cons170)
    rule5872 = ReplacementRule(pattern5872, replacement5872)
    pattern5873 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/sinh(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons33, cons170)
    rule5873 = ReplacementRule(pattern5873, replacement5873)
    pattern5874 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons89, cons167, cons1646)
    rule5874 = ReplacementRule(pattern5874, replacement5874)
    pattern5875 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons89, cons167, cons1646)
    rule5875 = ReplacementRule(pattern5875, replacement5875)
    pattern5876 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons95, cons167, cons1646, cons168)
    rule5876 = ReplacementRule(pattern5876, replacement5876)
    pattern5877 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons95, cons167, cons1646, cons168)
    rule5877 = ReplacementRule(pattern5877, replacement5877)
    pattern5878 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons89, cons91)
    rule5878 = ReplacementRule(pattern5878, replacement5878)
    pattern5879 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons89, cons91)
    rule5879 = ReplacementRule(pattern5879, replacement5879)
    pattern5880 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons95, cons91, cons168)
    rule5880 = ReplacementRule(pattern5880, replacement5880)
    pattern5881 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons95, cons91, cons168)
    rule5881 = ReplacementRule(pattern5881, replacement5881)
    pattern5882 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons25)
    rule5882 = ReplacementRule(pattern5882, replacement5882)
    pattern5883 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons25)
    rule5883 = ReplacementRule(pattern5883, replacement5883)
    pattern5884 = Pattern(Integral((a_ + WC('b', S(1))/cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule5884 = ReplacementRule(pattern5884, replacement5884)
    pattern5885 = Pattern(Integral((a_ + WC('b', S(1))/sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule5885 = ReplacementRule(pattern5885, replacement5885)
    pattern5886 = Pattern(Integral((a_ + WC('b', S(1))/cosh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons198, cons64)
    rule5886 = ReplacementRule(pattern5886, replacement5886)
    pattern5887 = Pattern(Integral((a_ + WC('b', S(1))/sinh(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons198, cons64)
    rule5887 = ReplacementRule(pattern5887, replacement5887)
    pattern5888 = Pattern(Integral(u_**WC('m', S(1))*(S(1)/cosh(v_))**WC('n', S(1)), x_), cons19, cons4, cons812, cons813)
    rule5888 = ReplacementRule(pattern5888, replacement5888)
    pattern5889 = Pattern(Integral(u_**WC('m', S(1))*(S(1)/sinh(v_))**WC('n', S(1)), x_), cons19, cons4, cons812, cons813)
    rule5889 = ReplacementRule(pattern5889, replacement5889)
    pattern5890 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule5890 = ReplacementRule(pattern5890, replacement5890)
    pattern5891 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule5891 = ReplacementRule(pattern5891, replacement5891)
    pattern5892 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule5892 = ReplacementRule(pattern5892, replacement5892)
    pattern5893 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule5893 = ReplacementRule(pattern5893, replacement5893)
    pattern5894 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5894 = ReplacementRule(pattern5894, replacement5894)
    pattern5895 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule5895 = ReplacementRule(pattern5895, replacement5895)
    pattern5896 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cosh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule5896 = ReplacementRule(pattern5896, replacement5896)
    pattern5897 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sinh(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule5897 = ReplacementRule(pattern5897, replacement5897)
    pattern5898 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cosh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5898 = ReplacementRule(pattern5898, replacement5898)
    pattern5899 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sinh(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule5899 = ReplacementRule(pattern5899, replacement5899)
    pattern5900 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule5900 = ReplacementRule(pattern5900, replacement5900)
    pattern5901 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule5901 = ReplacementRule(pattern5901, replacement5901)
    pattern5902 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule5902 = ReplacementRule(pattern5902, replacement5902)
    pattern5903 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule5903 = ReplacementRule(pattern5903, replacement5903)
    pattern5904 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cosh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5904 = ReplacementRule(pattern5904, replacement5904)
    pattern5905 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sinh(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5905 = ReplacementRule(pattern5905, replacement5905)
    pattern5906 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cosh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5906 = ReplacementRule(pattern5906, replacement5906)
    pattern5907 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sinh(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule5907 = ReplacementRule(pattern5907, replacement5907)
    pattern5908 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**p_*sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1647)
    rule5908 = ReplacementRule(pattern5908, replacement5908)
    pattern5909 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**p_*cosh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1647)
    rule5909 = ReplacementRule(pattern5909, replacement5909)
    pattern5910 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule5910 = ReplacementRule(pattern5910, replacement5910)
    pattern5911 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule5911 = ReplacementRule(pattern5911, replacement5911)
    pattern5912 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule5912 = ReplacementRule(pattern5912, replacement5912)
    pattern5913 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule5913 = ReplacementRule(pattern5913, replacement5913)
    pattern5914 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule5914 = ReplacementRule(pattern5914, replacement5914)
    pattern5915 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1685, cons33, cons170)
    rule5915 = ReplacementRule(pattern5915, replacement5915)
    pattern5916 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1685, cons33, cons170)
    rule5916 = ReplacementRule(pattern5916, replacement5916)
    pattern5917 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/cosh(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule5917 = ReplacementRule(pattern5917, replacement5917)
    pattern5918 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))/sinh(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule5918 = ReplacementRule(pattern5918, replacement5918)
    pattern5919 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**p_/cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons1410)
    rule5919 = ReplacementRule(pattern5919, replacement5919)
    pattern5920 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1410)
    rule5920 = ReplacementRule(pattern5920, replacement5920)
    pattern5921 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**p_/sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons1410)
    rule5921 = ReplacementRule(pattern5921, replacement5921)
    pattern5922 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1410)
    rule5922 = ReplacementRule(pattern5922, replacement5922)

    pattern5923 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons64, cons1686)
    rule5923 = ReplacementRule(pattern5923, With5923)

    pattern5924 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons64, cons1686)
    rule5924 = ReplacementRule(pattern5924, With5924)
    pattern5925 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons33, cons87)
    rule5925 = ReplacementRule(pattern5925, replacement5925)

    pattern5926 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sinh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/cosh(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons378, cons33, cons170, cons1687)
    rule5926 = ReplacementRule(pattern5926, With5926)
    pattern5927 = Pattern(Integral(F_**v_*G_**w_*u_**WC('m', S(1)), x_), cons19, cons4, cons5, cons1862, cons1863, cons1690, cons814, cons815)
    rule5927 = ReplacementRule(pattern5927, replacement5927)
    pattern5928 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5928 = ReplacementRule(pattern5928, replacement5928)
    pattern5929 = Pattern(Integral((a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5929 = ReplacementRule(pattern5929, replacement5929)
    pattern5930 = Pattern(Integral((a_ + WC('b', S(1))*tanh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/cosh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5930 = ReplacementRule(pattern5930, replacement5930)
    pattern5931 = Pattern(Integral((a_ + WC('b', S(1))/tanh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/sinh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5931 = ReplacementRule(pattern5931, replacement5931)
    pattern5932 = Pattern(Integral((a_ + WC('b', S(1))/cosh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*tanh(x_*WC('d', S(1)) + WC('c', S(0)))/cosh(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5932 = ReplacementRule(pattern5932, replacement5932)
    pattern5933 = Pattern(Integral((a_ + WC('b', S(1))/sinh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/(sinh(x_*WC('d', S(1)) + WC('c', S(0)))*tanh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule5933 = ReplacementRule(pattern5933, replacement5933)
    pattern5934 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons557, cons20)
    rule5934 = ReplacementRule(pattern5934, replacement5934)
    pattern5935 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons557, cons20)
    rule5935 = ReplacementRule(pattern5935, replacement5935)
    pattern5936 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons557)
    rule5936 = ReplacementRule(pattern5936, replacement5936)
    pattern5937 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1864, cons1865, cons557, cons27, cons1693)
    rule5937 = ReplacementRule(pattern5937, replacement5937)
    pattern5938 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1866)
    rule5938 = ReplacementRule(pattern5938, replacement5938)
    pattern5939 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cosh(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1866)
    rule5939 = ReplacementRule(pattern5939, replacement5939)
    pattern5940 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1867, cons89, cons167)
    rule5940 = ReplacementRule(pattern5940, replacement5940)
    pattern5941 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1867, cons89, cons167)
    rule5941 = ReplacementRule(pattern5941, replacement5941)
    pattern5942 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1868, cons586, cons1397)
    rule5942 = ReplacementRule(pattern5942, replacement5942)
    pattern5943 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1868, cons586, cons1397)
    rule5943 = ReplacementRule(pattern5943, replacement5943)
    pattern5944 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1869, cons89, cons91, cons1444)
    rule5944 = ReplacementRule(pattern5944, replacement5944)
    pattern5945 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1869, cons89, cons91, cons1444)
    rule5945 = ReplacementRule(pattern5945, replacement5945)
    pattern5946 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons25)
    rule5946 = ReplacementRule(pattern5946, replacement5946)
    pattern5947 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons25)
    rule5947 = ReplacementRule(pattern5947, replacement5947)
    pattern5948 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule5948 = ReplacementRule(pattern5948, replacement5948)
    pattern5949 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/tanh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule5949 = ReplacementRule(pattern5949, replacement5949)
    pattern5950 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cosh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1695, cons89, cons91)
    rule5950 = ReplacementRule(pattern5950, replacement5950)
    pattern5951 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sinh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1695, cons89, cons91)
    rule5951 = ReplacementRule(pattern5951, replacement5951)
    pattern5952 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cosh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1870, cons1504, cons965)
    rule5952 = ReplacementRule(pattern5952, replacement5952)
    pattern5953 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sinh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1870, cons1504, cons965)
    rule5953 = ReplacementRule(pattern5953, replacement5953)
    pattern5954 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cosh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1871, cons89, cons167, cons1646)
    rule5954 = ReplacementRule(pattern5954, replacement5954)
    pattern5955 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sinh(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1871, cons89, cons167, cons1646)
    rule5955 = ReplacementRule(pattern5955, replacement5955)
    pattern5956 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule5956 = ReplacementRule(pattern5956, replacement5956)
    pattern5957 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule5957 = ReplacementRule(pattern5957, replacement5957)
    pattern5958 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons25)
    rule5958 = ReplacementRule(pattern5958, replacement5958)
    pattern5959 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons25)
    rule5959 = ReplacementRule(pattern5959, replacement5959)
    pattern5960 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1872, cons198)
    rule5960 = ReplacementRule(pattern5960, replacement5960)
    pattern5961 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1702, cons198)
    rule5961 = ReplacementRule(pattern5961, replacement5961)
    pattern5962 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1013, cons198)
    rule5962 = ReplacementRule(pattern5962, replacement5962)
    pattern5963 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1872, cons152, cons1553)
    rule5963 = ReplacementRule(pattern5963, replacement5963)
    pattern5964 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1702, cons152, cons1553)
    rule5964 = ReplacementRule(pattern5964, replacement5964)
    pattern5965 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1013, cons152, cons1553)
    rule5965 = ReplacementRule(pattern5965, replacement5965)
    pattern5966 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons1872, cons1703, cons1704)
    rule5966 = ReplacementRule(pattern5966, replacement5966)
    pattern5967 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons1701, cons1873, cons1705)
    rule5967 = ReplacementRule(pattern5967, replacement5967)
    pattern5968 = Pattern(Integral(F_**(u_*WC('c', S(1)))*G_**v_, x_), cons1101, cons8, cons4, cons1863, cons812, cons813)
    rule5968 = ReplacementRule(pattern5968, replacement5968)

    pattern5969 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons33, cons170, cons150)
    rule5969 = ReplacementRule(pattern5969, With5969)

    pattern5970 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons33, cons170, cons150)
    rule5970 = ReplacementRule(pattern5970, With5970)
    pattern5971 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cosh(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons530)
    rule5971 = ReplacementRule(pattern5971, replacement5971)
    pattern5972 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('p', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cosh(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1706)
    rule5972 = ReplacementRule(pattern5972, replacement5972)
    pattern5973 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*G_**(x_*WC('e', S(1)) + WC('d', S(0)))*H_**(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons530, cons1863, cons1874)
    rule5973 = ReplacementRule(pattern5973, replacement5973)
    pattern5974 = Pattern(Integral(F_**u_*sinh(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons150)
    rule5974 = ReplacementRule(pattern5974, replacement5974)
    pattern5975 = Pattern(Integral(F_**u_*cosh(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons150)
    rule5975 = ReplacementRule(pattern5975, replacement5975)
    pattern5976 = Pattern(Integral(F_**u_*sinh(v_)**WC('m', S(1))*cosh(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons530)
    rule5976 = ReplacementRule(pattern5976, replacement5976)
    pattern5977 = Pattern(Integral(sinh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons8, cons1875)
    rule5977 = ReplacementRule(pattern5977, replacement5977)
    pattern5978 = Pattern(Integral(cosh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons8, cons1875)
    rule5978 = ReplacementRule(pattern5978, replacement5978)
    pattern5979 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1876, cons56)
    rule5979 = ReplacementRule(pattern5979, replacement5979)
    pattern5980 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1876, cons56)
    rule5980 = ReplacementRule(pattern5980, replacement5980)
    pattern5981 = Pattern(Integral(sqrt(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))), x_), cons2, cons3, cons8, cons4, cons1877)
    rule5981 = ReplacementRule(pattern5981, replacement5981)
    pattern5982 = Pattern(Integral(sqrt(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))), x_), cons2, cons3, cons8, cons4, cons1877)
    rule5982 = ReplacementRule(pattern5982, replacement5982)
    pattern5983 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons130, cons1878)
    rule5983 = ReplacementRule(pattern5983, replacement5983)
    pattern5984 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons130, cons1878)
    rule5984 = ReplacementRule(pattern5984, replacement5984)
    pattern5985 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1879)
    rule5985 = ReplacementRule(pattern5985, replacement5985)
    pattern5986 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1879)
    rule5986 = ReplacementRule(pattern5986, replacement5986)
    pattern5987 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1880)
    rule5987 = ReplacementRule(pattern5987, replacement5987)
    pattern5988 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1880)
    rule5988 = ReplacementRule(pattern5988, replacement5988)
    pattern5989 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1507, cons1881)
    rule5989 = ReplacementRule(pattern5989, replacement5989)
    pattern5990 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1507, cons1881)
    rule5990 = ReplacementRule(pattern5990, replacement5990)
    pattern5991 = Pattern(Integral(sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1880)
    rule5991 = ReplacementRule(pattern5991, replacement5991)
    pattern5992 = Pattern(Integral(cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1880)
    rule5992 = ReplacementRule(pattern5992, replacement5992)
    pattern5993 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1882, cons56, cons68)
    rule5993 = ReplacementRule(pattern5993, replacement5993)
    pattern5994 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1882, cons56, cons68)
    rule5994 = ReplacementRule(pattern5994, replacement5994)
    pattern5995 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons130, cons1883)
    rule5995 = ReplacementRule(pattern5995, replacement5995)
    pattern5996 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons130, cons1883)
    rule5996 = ReplacementRule(pattern5996, replacement5996)
    pattern5997 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1884, cons68)
    rule5997 = ReplacementRule(pattern5997, replacement5997)
    pattern5998 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1884, cons68)
    rule5998 = ReplacementRule(pattern5998, replacement5998)
    pattern5999 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons1885, cons13, cons148, cons68)
    rule5999 = ReplacementRule(pattern5999, replacement5999)
    pattern6000 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons1885, cons13, cons148, cons68)
    rule6000 = ReplacementRule(pattern6000, replacement6000)
    pattern6001 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons1886, cons13, cons139, cons1507, cons68)
    rule6001 = ReplacementRule(pattern6001, replacement6001)
    pattern6002 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons1886, cons13, cons139, cons1507, cons68)
    rule6002 = ReplacementRule(pattern6002, replacement6002)
    pattern6003 = Pattern(Integral(x_**WC('m', S(1))*sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1885)
    rule6003 = ReplacementRule(pattern6003, replacement6003)
    pattern6004 = Pattern(Integral(x_**WC('m', S(1))*cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1885)
    rule6004 = ReplacementRule(pattern6004, replacement6004)
    pattern6005 = Pattern(Integral((S(1)/cosh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons8, cons1875)
    rule6005 = ReplacementRule(pattern6005, replacement6005)
    pattern6006 = Pattern(Integral((S(1)/sinh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons8, cons1875)
    rule6006 = ReplacementRule(pattern6006, replacement6006)
    pattern6007 = Pattern(Integral(S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1887)
    rule6007 = ReplacementRule(pattern6007, replacement6007)
    pattern6008 = Pattern(Integral(S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1887)
    rule6008 = ReplacementRule(pattern6008, replacement6008)
    pattern6009 = Pattern(Integral((S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1888, cons1647)
    rule6009 = ReplacementRule(pattern6009, replacement6009)
    pattern6010 = Pattern(Integral((S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1888, cons1647)
    rule6010 = ReplacementRule(pattern6010, replacement6010)
    pattern6011 = Pattern(Integral((S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1722, cons1889)
    rule6011 = ReplacementRule(pattern6011, replacement6011)
    pattern6012 = Pattern(Integral((S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1722, cons1889)
    rule6012 = ReplacementRule(pattern6012, replacement6012)
    pattern6013 = Pattern(Integral((S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1880)
    rule6013 = ReplacementRule(pattern6013, replacement6013)
    pattern6014 = Pattern(Integral((S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1880)
    rule6014 = ReplacementRule(pattern6014, replacement6014)
    pattern6015 = Pattern(Integral((S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons1880)
    rule6015 = ReplacementRule(pattern6015, replacement6015)
    pattern6016 = Pattern(Integral((S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons1880)
    rule6016 = ReplacementRule(pattern6016, replacement6016)
    pattern6017 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons8, cons1890)
    rule6017 = ReplacementRule(pattern6017, replacement6017)
    pattern6018 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons8, cons1890)
    rule6018 = ReplacementRule(pattern6018, replacement6018)
    pattern6019 = Pattern(Integral(x_**WC('m', S(1))/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1891)
    rule6019 = ReplacementRule(pattern6019, replacement6019)
    pattern6020 = Pattern(Integral(x_**WC('m', S(1))/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1891)
    rule6020 = ReplacementRule(pattern6020, replacement6020)
    pattern6021 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1725, cons68, cons1647)
    rule6021 = ReplacementRule(pattern6021, replacement6021)
    pattern6022 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1725, cons68, cons1647)
    rule6022 = ReplacementRule(pattern6022, replacement6022)
    pattern6023 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1722, cons1892)
    rule6023 = ReplacementRule(pattern6023, replacement6023)
    pattern6024 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1722, cons1892)
    rule6024 = ReplacementRule(pattern6024, replacement6024)
    pattern6025 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1885)
    rule6025 = ReplacementRule(pattern6025, replacement6025)
    pattern6026 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1885)
    rule6026 = ReplacementRule(pattern6026, replacement6026)
    pattern6027 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cosh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1885)
    rule6027 = ReplacementRule(pattern6027, replacement6027)
    pattern6028 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sinh(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1885)
    rule6028 = ReplacementRule(pattern6028, replacement6028)
    pattern6029 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons165)
    rule6029 = ReplacementRule(pattern6029, replacement6029)
    pattern6030 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cosh(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons165)
    rule6030 = ReplacementRule(pattern6030, replacement6030)
    pattern6031 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons340, cons165)
    rule6031 = ReplacementRule(pattern6031, replacement6031)
    pattern6032 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cosh(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons340, cons165)
    rule6032 = ReplacementRule(pattern6032, replacement6032)
    pattern6033 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons55, cons13, cons165)
    rule6033 = ReplacementRule(pattern6033, replacement6033)
    pattern6034 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*cosh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons55, cons13, cons165)
    rule6034 = ReplacementRule(pattern6034, replacement6034)
    pattern6035 = Pattern(Integral(x_**m_*log(x_*WC('b', S(1)))**WC('p', S(1))*sinh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons164, cons165, cons629)
    rule6035 = ReplacementRule(pattern6035, replacement6035)
    pattern6036 = Pattern(Integral(x_**m_*log(x_*WC('b', S(1)))**WC('p', S(1))*cosh(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons164, cons165, cons629)
    rule6036 = ReplacementRule(pattern6036, replacement6036)
    pattern6037 = Pattern(Integral(sinh(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons8, cons29, cons150)
    rule6037 = ReplacementRule(pattern6037, replacement6037)
    pattern6038 = Pattern(Integral(cosh(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons8, cons29, cons150)
    rule6038 = ReplacementRule(pattern6038, replacement6038)
    pattern6039 = Pattern(Integral(sinh((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150, cons73)
    rule6039 = ReplacementRule(pattern6039, replacement6039)
    pattern6040 = Pattern(Integral(cosh((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150, cons73)
    rule6040 = ReplacementRule(pattern6040, replacement6040)

    pattern6041 = Pattern(Integral(sinh(u_)**WC('n', S(1)), x_), cons150, cons1727)
    rule6041 = ReplacementRule(pattern6041, With6041)

    pattern6042 = Pattern(Integral(cosh(u_)**WC('n', S(1)), x_), cons150, cons1727)
    rule6042 = ReplacementRule(pattern6042, With6042)
    pattern6043 = Pattern(Integral(WC('u', S(1))*sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), cons1690)
    rule6043 = ReplacementRule(pattern6043, replacement6043)
    pattern6044 = Pattern(Integral(WC('u', S(1))*cosh(v_)**WC('p', S(1))*cosh(w_)**WC('q', S(1)), x_), cons1690)
    rule6044 = ReplacementRule(pattern6044, replacement6044)
    pattern6045 = Pattern(Integral(sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), cons557, cons1728)
    rule6045 = ReplacementRule(pattern6045, replacement6045)
    pattern6046 = Pattern(Integral(cosh(v_)**WC('p', S(1))*cosh(w_)**WC('q', S(1)), x_), cons557, cons1728)
    rule6046 = ReplacementRule(pattern6046, replacement6046)
    pattern6047 = Pattern(Integral(x_**WC('m', S(1))*sinh(v_)**WC('p', S(1))*sinh(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule6047 = ReplacementRule(pattern6047, replacement6047)
    pattern6048 = Pattern(Integral(x_**WC('m', S(1))*cosh(v_)**WC('p', S(1))*cosh(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule6048 = ReplacementRule(pattern6048, replacement6048)
    pattern6049 = Pattern(Integral(WC('u', S(1))*sinh(v_)**WC('p', S(1))*cosh(w_)**WC('p', S(1)), x_), cons1690, cons40)
    rule6049 = ReplacementRule(pattern6049, replacement6049)
    pattern6050 = Pattern(Integral(sinh(v_)**WC('p', S(1))*cosh(w_)**WC('q', S(1)), x_), cons557, cons1728)
    rule6050 = ReplacementRule(pattern6050, replacement6050)
    pattern6051 = Pattern(Integral(x_**WC('m', S(1))*sinh(v_)**WC('p', S(1))*cosh(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule6051 = ReplacementRule(pattern6051, replacement6051)
    pattern6052 = Pattern(Integral(sinh(v_)*tanh(w_)**WC('n', S(1)), x_), cons89, cons90, cons1730)
    rule6052 = ReplacementRule(pattern6052, replacement6052)
    pattern6053 = Pattern(Integral((S(1)/tanh(w_))**WC('n', S(1))*cosh(v_), x_), cons89, cons90, cons1730)
    rule6053 = ReplacementRule(pattern6053, replacement6053)
    pattern6054 = Pattern(Integral((S(1)/tanh(w_))**WC('n', S(1))*sinh(v_), x_), cons89, cons90, cons1730)
    rule6054 = ReplacementRule(pattern6054, replacement6054)
    pattern6055 = Pattern(Integral(cosh(v_)*tanh(w_)**WC('n', S(1)), x_), cons89, cons90, cons1730)
    rule6055 = ReplacementRule(pattern6055, replacement6055)
    pattern6056 = Pattern(Integral((S(1)/cosh(w_))**WC('n', S(1))*sinh(v_), x_), cons89, cons90, cons1730)
    rule6056 = ReplacementRule(pattern6056, replacement6056)
    pattern6057 = Pattern(Integral((S(1)/sinh(w_))**WC('n', S(1))*cosh(v_), x_), cons89, cons90, cons1730)
    rule6057 = ReplacementRule(pattern6057, replacement6057)
    pattern6058 = Pattern(Integral((S(1)/sinh(w_))**WC('n', S(1))*sinh(v_), x_), cons89, cons90, cons1730)
    rule6058 = ReplacementRule(pattern6058, replacement6058)
    pattern6059 = Pattern(Integral((S(1)/cosh(w_))**WC('n', S(1))*cosh(v_), x_), cons89, cons90, cons1730)
    rule6059 = ReplacementRule(pattern6059, replacement6059)
    pattern6060 = Pattern(Integral((a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule6060 = ReplacementRule(pattern6060, replacement6060)
    pattern6061 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons1458, cons152, cons170, cons465, cons1731)
    rule6061 = ReplacementRule(pattern6061, replacement6061)
    pattern6062 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons1480, cons152, cons170, cons465, cons1731)
    rule6062 = ReplacementRule(pattern6062, replacement6062)
    pattern6063 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons13)
    rule6063 = ReplacementRule(pattern6063, replacement6063)
    pattern6064 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons13)
    rule6064 = ReplacementRule(pattern6064, replacement6064)
    pattern6065 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule6065 = ReplacementRule(pattern6065, replacement6065)
    pattern6066 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((b_ + WC('c', S(1))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons64)
    rule6066 = ReplacementRule(pattern6066, replacement6066)
    pattern6067 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((WC('a', S(1))/cosh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(0)) + WC('c', S(1))*tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*cosh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule6067 = ReplacementRule(pattern6067, replacement6067)
    pattern6068 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((c_ + WC('b', S(1))/tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons64)
    rule6068 = ReplacementRule(pattern6068, replacement6068)
    pattern6069 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((WC('a', S(1))/sinh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(1))/tanh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(0)))*sinh(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule6069 = ReplacementRule(pattern6069, replacement6069)
    pattern6070 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64)
    rule6070 = ReplacementRule(pattern6070, replacement6070)
    pattern6071 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64)
    rule6071 = ReplacementRule(pattern6071, replacement6071)
    pattern6072 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1441)
    rule6072 = ReplacementRule(pattern6072, replacement6072)
    pattern6073 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1267)
    rule6073 = ReplacementRule(pattern6073, replacement6073)
    pattern6074 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1442)
    rule6074 = ReplacementRule(pattern6074, replacement6074)
    pattern6075 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1269)
    rule6075 = ReplacementRule(pattern6075, replacement6075)
    pattern6076 = Pattern(Integral((A_ + WC('B', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1893)
    rule6076 = ReplacementRule(pattern6076, replacement6076)
    pattern6077 = Pattern(Integral((A_ + WC('B', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1476)
    rule6077 = ReplacementRule(pattern6077, replacement6077)
    pattern6078 = Pattern(Integral((a_ + WC('b', S(1))*tanh(v_))**WC('n', S(1))*(S(1)/cosh(v_))**WC('m', S(1)), x_), cons2, cons3, cons152, cons1553, cons1483)
    rule6078 = ReplacementRule(pattern6078, replacement6078)
    pattern6079 = Pattern(Integral((a_ + WC('b', S(1))/tanh(v_))**WC('n', S(1))*(S(1)/sinh(v_))**WC('m', S(1)), x_), cons2, cons3, cons152, cons1553, cons1483)
    rule6079 = ReplacementRule(pattern6079, replacement6079)
    pattern6080 = Pattern(Integral(WC('u', S(1))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*sinh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons530)
    rule6080 = ReplacementRule(pattern6080, replacement6080)
    pattern6081 = Pattern(Integral(WC('u', S(1))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*cosh(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons530)
    rule6081 = ReplacementRule(pattern6081, replacement6081)
    pattern6082 = Pattern(Integral(S(1)/(cosh(c_ + x_*WC('d', S(1)))*cosh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule6082 = ReplacementRule(pattern6082, replacement6082)
    pattern6083 = Pattern(Integral(S(1)/(sinh(c_ + x_*WC('d', S(1)))*sinh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule6083 = ReplacementRule(pattern6083, replacement6083)
    pattern6084 = Pattern(Integral(tanh(c_ + x_*WC('d', S(1)))*tanh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule6084 = ReplacementRule(pattern6084, replacement6084)
    pattern6085 = Pattern(Integral(S(1)/(tanh(c_ + x_*WC('d', S(1)))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule6085 = ReplacementRule(pattern6085, replacement6085)
    pattern6086 = Pattern(Integral((WC('a', S(1))*cosh(v_) + WC('b', S(1))*sinh(v_))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1267)
    rule6086 = ReplacementRule(pattern6086, replacement6086)
    return [rule5646, rule5647, rule5648, rule5649, rule5650, rule5651, rule5652, rule5653, rule5654, rule5655, rule5656, rule5657, rule5658, rule5659, rule5660, rule5661, rule5662, rule5663, rule5664, rule5665, rule5666, rule5667, rule5668, rule5669, rule5670, rule5671, rule5672, rule5673, rule5674, rule5675, rule5676, rule5677, rule5678, rule5679, rule5680, rule5681, rule5682, rule5683, rule5684, rule5685, rule5686, rule5687, rule5688, rule5689, rule5690, rule5691, rule5692, rule5693, rule5694, rule5695, rule5696, rule5697, rule5698, rule5699, rule5700, rule5701, rule5702, rule5703, rule5704, rule5705, rule5706, rule5707, rule5708, rule5709, rule5710, rule5711, rule5712, rule5713, rule5714, rule5715, rule5716, rule5717, rule5718, rule5719, rule5720, rule5721, rule5722, rule5723, rule5724, rule5725, rule5726, rule5727, rule5728, rule5729, rule5730, rule5731, rule5732, rule5733, rule5734, rule5735, rule5736, rule5737, rule5738, rule5739, rule5740, rule5741, rule5742, rule5743, rule5744, rule5745, rule5746, rule5747, rule5748, rule5749, rule5750, rule5751, rule5752, rule5753, rule5754, rule5755, rule5756, rule5757, rule5758, rule5759, rule5760, rule5761, rule5762, rule5763, rule5764, rule5765, rule5766, rule5767, rule5768, rule5769, rule5770, rule5771, rule5772, rule5773, rule5774, rule5775, rule5776, rule5777, rule5778, rule5779, rule5780, rule5781, rule5782, rule5783, rule5784, rule5785, rule5786, rule5787, rule5788, rule5789, rule5790, rule5791, rule5792, rule5793, rule5794, rule5795, rule5796, rule5797, rule5798, rule5799, rule5800, rule5801, rule5802, rule5803, rule5804, rule5805, rule5806, rule5807, rule5808, rule5809, rule5810, rule5811, rule5812, rule5813, rule5814, rule5815, rule5816, rule5817, rule5818, rule5819, rule5820, rule5821, rule5822, rule5823, rule5824, rule5825, rule5826, rule5827, rule5828, rule5829, rule5830, rule5831, rule5832, rule5833, rule5834, rule5835, rule5836, rule5837, rule5838, rule5839, rule5840, rule5841, rule5842, rule5843, rule5844, rule5845, rule5846, rule5847, rule5848, rule5849, rule5850, rule5851, rule5852, rule5853, rule5854, rule5855, rule5856, rule5857, rule5858, rule5859, rule5860, rule5861, rule5862, rule5863, rule5864, rule5865, rule5866, rule5867, rule5868, rule5869, rule5870, rule5871, rule5872, rule5873, rule5874, rule5875, rule5876, rule5877, rule5878, rule5879, rule5880, rule5881, rule5882, rule5883, rule5884, rule5885, rule5886, rule5887, rule5888, rule5889, rule5890, rule5891, rule5892, rule5893, rule5894, rule5895, rule5896, rule5897, rule5898, rule5899, rule5900, rule5901, rule5902, rule5903, rule5904, rule5905, rule5906, rule5907, rule5908, rule5909, rule5910, rule5911, rule5912, rule5913, rule5914, rule5915, rule5916, rule5917, rule5918, rule5919, rule5920, rule5921, rule5922, rule5923, rule5924, rule5925, rule5926, rule5927, rule5928, rule5929, rule5930, rule5931, rule5932, rule5933, rule5934, rule5935, rule5936, rule5937, rule5938, rule5939, rule5940, rule5941, rule5942, rule5943, rule5944, rule5945, rule5946, rule5947, rule5948, rule5949, rule5950, rule5951, rule5952, rule5953, rule5954, rule5955, rule5956, rule5957, rule5958, rule5959, rule5960, rule5961, rule5962, rule5963, rule5964, rule5965, rule5966, rule5967, rule5968, rule5969, rule5970, rule5971, rule5972, rule5973, rule5974, rule5975, rule5976, rule5977, rule5978, rule5979, rule5980, rule5981, rule5982, rule5983, rule5984, rule5985, rule5986, rule5987, rule5988, rule5989, rule5990, rule5991, rule5992, rule5993, rule5994, rule5995, rule5996, rule5997, rule5998, rule5999, rule6000, rule6001, rule6002, rule6003, rule6004, rule6005, rule6006, rule6007, rule6008, rule6009, rule6010, rule6011, rule6012, rule6013, rule6014, rule6015, rule6016, rule6017, rule6018, rule6019, rule6020, rule6021, rule6022, rule6023, rule6024, rule6025, rule6026, rule6027, rule6028, rule6029, rule6030, rule6031, rule6032, rule6033, rule6034, rule6035, rule6036, rule6037, rule6038, rule6039, rule6040, rule6041, rule6042, rule6043, rule6044, rule6045, rule6046, rule6047, rule6048, rule6049, rule6050, rule6051, rule6052, rule6053, rule6054, rule6055, rule6056, rule6057, rule6058, rule6059, rule6060, rule6061, rule6062, rule6063, rule6064, rule6065, rule6066, rule6067, rule6068, rule6069, rule6070, rule6071, rule6072, rule6073, rule6074, rule6075, rule6076, rule6077, rule6078, rule6079, rule6080, rule6081, rule6082, rule6083, rule6084, rule6085, rule6086, ]



def replacement5646(c, m, x, d, f, e):
        # rubi.append(5646)
        return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*cosh(e + f*x), x), x) + Simp((c + d*x)**m*cosh(e + f*x)/f, x)
def replacement5647(c, m, x, d, f, e):
        # rubi.append(5647)
        return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*sinh(e + f*x), x), x) + Simp((c + d*x)**m*sinh(e + f*x)/f, x)
def replacement5648(c, m, x, d, f, e):
        # rubi.append(5648)
        return -Dist(f/(d*(m + S(1))), Int((c + d*x)**(m + S(1))*cosh(e + f*x), x), x) + Simp((c + d*x)**(m + S(1))*sinh(e + f*x)/(d*(m + S(1))), x)
def replacement5649(c, m, x, d, f, e):
        # rubi.append(5649)
        return -Dist(f/(d*(m + S(1))), Int((c + d*x)**(m + S(1))*sinh(e + f*x), x), x) + Simp((c + d*x)**(m + S(1))*cosh(e + f*x)/(d*(m + S(1))), x)
def replacement5650(c, x, d, f, e):
        # rubi.append(5650)
        return Simp(SinhIntegral(e + f*x)/d, x)
def replacement5651(c, x, d, f, e):
        # rubi.append(5651)
        return Simp(CoshIntegral(e + f*x)/d, x)
def replacement5652(c, x, d, f, e):
        # rubi.append(5652)
        return Dist(sinh((-c*f + d*e)/d), Int(cosh(c*f/d + f*x)/(c + d*x), x), x) + Dist(cosh((-c*f + d*e)/d), Int(sinh(c*f/d + f*x)/(c + d*x), x), x)
def replacement5653(c, x, d, f, e):
        # rubi.append(5653)
        return Dist(sinh((-c*f + d*e)/d), Int(sinh(c*f/d + f*x)/(c + d*x), x), x) + Dist(cosh((-c*f + d*e)/d), Int(cosh(c*f/d + f*x)/(c + d*x), x), x)
def replacement5654(c, m, x, d, f, e):
        # rubi.append(5654)
        return -Dist(S(1)/2, Int((c + d*x)**m*exp(-e - f*x), x), x) + Dist(S(1)/2, Int((c + d*x)**m*exp(e + f*x), x), x)
def replacement5655(c, m, x, d, f, e):
        # rubi.append(5655)
        return Dist(S(1)/2, Int((c + d*x)**m*exp(-e - f*x), x), x) + Dist(S(1)/2, Int((c + d*x)**m*exp(e + f*x), x), x)
def replacement5656(c, n, x, d, f, e, b):
        # rubi.append(5656)
        return -Dist(b**S(2)*(n + S(-1))/n, Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x), x), x) - Simp(d*(b*sinh(e + f*x))**n/(f**S(2)*n**S(2)), x) + Simp(b*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)*cosh(e + f*x)/(f*n), x)
def replacement5657(c, n, x, d, f, e, b):
        # rubi.append(5657)
        return Dist(b**S(2)*(n + S(-1))/n, Int((b*cosh(e + f*x))**(n + S(-2))*(c + d*x), x), x) - Simp(d*(b*cosh(e + f*x))**n/(f**S(2)*n**S(2)), x) + Simp(b*(b*cosh(e + f*x))**(n + S(-1))*(c + d*x)*sinh(e + f*x)/(f*n), x)
def replacement5658(c, m, n, x, d, f, e, b):
        # rubi.append(5658)
        return -Dist(b**S(2)*(n + S(-1))/n, Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) + Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b*sinh(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) + Simp(b*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)**m*cosh(e + f*x)/(f*n), x) - Simp(d*m*(b*sinh(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)
def replacement5659(c, m, n, x, d, f, e, b):
        # rubi.append(5659)
        return Dist(b**S(2)*(n + S(-1))/n, Int((b*cosh(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) + Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b*cosh(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) + Simp(b*(b*cosh(e + f*x))**(n + S(-1))*(c + d*x)**m*sinh(e + f*x)/(f*n), x) - Simp(d*m*(b*cosh(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)
def replacement5660(c, m, n, x, d, f, e):
        # rubi.append(5660)
        return Int(ExpandTrigReduce((c + d*x)**m, sinh(e + f*x)**n, x), x)
def replacement5661(c, m, n, x, d, f, e):
        # rubi.append(5661)
        return Int(ExpandTrigReduce((c + d*x)**m, cosh(e + f*x)**n, x), x)
def replacement5662(c, m, n, x, d, f, e):
        # rubi.append(5662)
        return -Dist(f*n/(d*(m + S(1))), Int(ExpandTrigReduce((c + d*x)**(m + S(1)), sinh(e + f*x)**(n + S(-1))*cosh(e + f*x), x), x), x) + Simp((c + d*x)**(m + S(1))*sinh(e + f*x)**n/(d*(m + S(1))), x)
def replacement5663(c, m, n, x, d, f, e):
        # rubi.append(5663)
        return -Dist(f*n/(d*(m + S(1))), Int(ExpandTrigReduce((c + d*x)**(m + S(1)), sinh(e + f*x)*cosh(e + f*x)**(n + S(-1)), x), x), x) + Simp((c + d*x)**(m + S(1))*cosh(e + f*x)**n/(d*(m + S(1))), x)
def replacement5664(c, m, n, x, d, f, e, b):
        # rubi.append(5664)
        return Dist(f**S(2)*n**S(2)/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*sinh(e + f*x))**n*(c + d*x)**(m + S(2)), x), x) + Dist(b**S(2)*f**S(2)*n*(n + S(-1))/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*sinh(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x), x) + Simp((b*sinh(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x) - Simp(b*f*n*(b*sinh(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*cosh(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement5665(c, m, n, x, d, f, e, b):
        # rubi.append(5665)
        return Dist(f**S(2)*n**S(2)/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*cosh(e + f*x))**n*(c + d*x)**(m + S(2)), x), x) - Dist(b**S(2)*f**S(2)*n*(n + S(-1))/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*cosh(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x), x) + Simp((b*cosh(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x) - Simp(b*f*n*(b*cosh(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*sinh(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement5666(c, n, x, d, f, e, b):
        # rubi.append(5666)
        return -Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x), x), x) - Simp(d*(b*sinh(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x) + Simp((b*sinh(e + f*x))**(n + S(1))*(c + d*x)*cosh(e + f*x)/(b*f*(n + S(1))), x)
def replacement5667(c, n, x, d, f, e, b):
        # rubi.append(5667)
        return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*cosh(e + f*x))**(n + S(2))*(c + d*x), x), x) + Simp(d*(b*cosh(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x) - Simp((b*cosh(e + f*x))**(n + S(1))*(c + d*x)*sinh(e + f*x)/(b*f*(n + S(1))), x)
def replacement5668(c, m, n, x, d, f, e, b):
        # rubi.append(5668)
        return -Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), Int((b*sinh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x), x) + Simp((b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*cosh(e + f*x)/(b*f*(n + S(1))), x) - Simp(d*m*(b*sinh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5669(c, m, n, x, d, f, e, b):
        # rubi.append(5669)
        return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*cosh(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) - Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), Int((b*cosh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x), x) - Simp((b*cosh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x)/(b*f*(n + S(1))), x) + Simp(d*m*(b*cosh(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5670(c, m, n, x, d, f, a, e, b):
        # rubi.append(5670)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b*sinh(e + f*x))**n, x), x)
def replacement5671(c, m, n, x, d, f, a, e, b):
        # rubi.append(5671)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b*cosh(e + f*x))**n, x), x)
def replacement5672(c, m, n, x, d, f, a, e, b):
        # rubi.append(5672)
        return Dist((S(2)*a)**n, Int((c + d*x)**m*cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5673(c, m, n, x, d, f, a, e, b):
        # rubi.append(5673)
        return Dist((S(2)*a)**IntPart(n)*(a + b*sinh(e + f*x))**FracPart(n)*cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*cosh(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5674(c, m, n, x, d, f, a, e, b):
        # rubi.append(5674)
        return Dist((S(2)*a)**n, Int((c + d*x)**m*cosh(e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5675(c, m, n, x, d, f, a, e, b):
        # rubi.append(5675)
        return Dist((-S(2)*a)**n, Int((c + d*x)**m*sinh(e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5676(c, m, n, x, d, f, a, e, b):
        # rubi.append(5676)
        return Dist((S(2)*a)**IntPart(n)*(a + b*cosh(e + f*x))**FracPart(n)*cosh(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*cosh(e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5677(c, m, n, x, d, f, a, e, b):
        # rubi.append(5677)
        return Dist((-S(2)*a)**IntPart(n)*(a + b*cosh(e + f*x))**FracPart(n)*sinh(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*sinh(e/S(2) + f*x/S(2))**(S(2)*n), x), x)
def replacement5678(c, m, x, d, f, a, e, b):
        # rubi.append(5678)
        return Dist(S(-2), Int((c + d*x)**m*exp(e + f*x)/(-S(2)*a*exp(e + f*x) - b*exp(S(2)*e + S(2)*f*x) + b), x), x)
def replacement5679(c, m, x, d, f, a, e, b):
        # rubi.append(5679)
        return Dist(S(2), Int((c + d*x)**m*exp(e + f*x)/(S(2)*a*exp(e + f*x) + b*exp(S(2)*e + S(2)*f*x) + b), x), x)
def replacement5680(c, m, x, d, f, a, e, b):
        # rubi.append(5680)
        return Dist(a/(a**S(2) + b**S(2)), Int((c + d*x)**m/(a + b*sinh(e + f*x)), x), x) + Dist(b*d*m/(f*(a**S(2) + b**S(2))), Int((c + d*x)**(m + S(-1))*cosh(e + f*x)/(a + b*sinh(e + f*x)), x), x) - Simp(b*(c + d*x)**m*cosh(e + f*x)/(f*(a + b*sinh(e + f*x))*(a**S(2) + b**S(2))), x)
def replacement5681(c, m, x, d, f, a, e, b):
        # rubi.append(5681)
        return Dist(a/(a**S(2) - b**S(2)), Int((c + d*x)**m/(a + b*cosh(e + f*x)), x), x) + Dist(b*d*m/(f*(a**S(2) - b**S(2))), Int((c + d*x)**(m + S(-1))*sinh(e + f*x)/(a + b*cosh(e + f*x)), x), x) - Simp(b*(c + d*x)**m*sinh(e + f*x)/(f*(a + b*cosh(e + f*x))*(a**S(2) - b**S(2))), x)
def replacement5682(c, m, n, x, d, f, a, e, b):
        # rubi.append(5682)
        return Dist(a/(a**S(2) + b**S(2)), Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m, x), x) - Dist(b*(n + S(2))/((a**S(2) + b**S(2))*(n + S(1))), Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x), x), x) - Dist(b*d*m/(f*(a**S(2) + b**S(2))*(n + S(1))), Int((a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*cosh(e + f*x), x), x) + Simp(b*(a + b*sinh(e + f*x))**(n + S(1))*(c + d*x)**m*cosh(e + f*x)/(f*(a**S(2) + b**S(2))*(n + S(1))), x)
def replacement5683(c, m, n, x, d, f, a, e, b):
        # rubi.append(5683)
        return Dist(a/(a**S(2) - b**S(2)), Int((a + b*cosh(e + f*x))**(n + S(1))*(c + d*x)**m, x), x) - Dist(b*(n + S(2))/((a**S(2) - b**S(2))*(n + S(1))), Int((a + b*cosh(e + f*x))**(n + S(1))*(c + d*x)**m*cosh(e + f*x), x), x) - Dist(b*d*m/(f*(a**S(2) - b**S(2))*(n + S(1))), Int((a + b*cosh(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*sinh(e + f*x), x), x) + Simp(b*(a + b*cosh(e + f*x))**(n + S(1))*(c + d*x)**m*sinh(e + f*x)/(f*(a**S(2) - b**S(2))*(n + S(1))), x)
def replacement5684(m, n, x, v, u, a, b):
        # rubi.append(5684)
        return Int((a + b*sinh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5685(m, n, x, v, u, a, b):
        # rubi.append(5685)
        return Int((a + b*cosh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5686(c, m, n, x, d, f, a, e, b):
        # rubi.append(5686)
        return Int((a + b*sinh(e + f*x))**n*(c + d*x)**m, x)
def replacement5687(c, m, n, x, d, f, a, e, b):
        # rubi.append(5687)
        return Int((a + b*cosh(e + f*x))**n*(c + d*x)**m, x)
def replacement5688(c, n, x, d, a, p, b):
        # rubi.append(5688)
        return Int(ExpandIntegrand(sinh(c + d*x), (a + b*x**n)**p, x), x)
def replacement5689(c, n, x, d, a, p, b):
        # rubi.append(5689)
        return Int(ExpandIntegrand(cosh(c + d*x), (a + b*x**n)**p, x), x)
def replacement5690(c, n, x, d, a, p, b):
        # rubi.append(5690)
        return -Dist(d/(b*n*(p + S(1))), Int(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*cosh(c + d*x), x), x) - Dist((S(1) - n)/(b*n*(p + S(1))), Int(x**(-n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x), x) + Simp(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5691(c, n, x, d, a, p, b):
        # rubi.append(5691)
        return -Dist(d/(b*n*(p + S(1))), Int(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x), x) - Dist((S(1) - n)/(b*n*(p + S(1))), Int(x**(-n)*(a + b*x**n)**(p + S(1))*cosh(c + d*x), x), x) + Simp(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*cosh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5692(c, n, x, d, a, p, b):
        # rubi.append(5692)
        return Int(ExpandIntegrand(sinh(c + d*x), (a + b*x**n)**p, x), x)
def replacement5693(c, n, x, d, a, p, b):
        # rubi.append(5693)
        return Int(ExpandIntegrand(cosh(c + d*x), (a + b*x**n)**p, x), x)
def replacement5694(c, n, x, d, a, p, b):
        # rubi.append(5694)
        return Int(x**(n*p)*(a*x**(-n) + b)**p*sinh(c + d*x), x)
def replacement5695(c, n, x, d, a, p, b):
        # rubi.append(5695)
        return Int(x**(n*p)*(a*x**(-n) + b)**p*cosh(c + d*x), x)
def replacement5696(c, n, x, d, a, p, b):
        # rubi.append(5696)
        return Int((a + b*x**n)**p*sinh(c + d*x), x)
def replacement5697(c, n, x, d, a, p, b):
        # rubi.append(5697)
        return Int((a + b*x**n)**p*cosh(c + d*x), x)
def replacement5698(c, m, n, x, d, a, p, e, b):
        # rubi.append(5698)
        return Int(ExpandIntegrand(sinh(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x)
def replacement5699(c, m, n, x, d, a, p, e, b):
        # rubi.append(5699)
        return Int(ExpandIntegrand(cosh(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x)
def replacement5700(c, m, n, x, d, a, p, e, b):
        # rubi.append(5700)
        return -Dist(d*e**m/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*cosh(c + d*x), x), x) + Simp(e**m*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5701(c, m, n, x, d, a, p, e, b):
        # rubi.append(5701)
        return -Dist(d*e**m/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*sinh(c + d*x), x), x) + Simp(e**m*(a + b*x**n)**(p + S(1))*cosh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5702(c, m, n, x, d, a, p, b):
        # rubi.append(5702)
        return -Dist(d/(b*n*(p + S(1))), Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*cosh(c + d*x), x), x) - Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x), x) + Simp(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5703(c, m, n, x, d, a, p, b):
        # rubi.append(5703)
        return -Dist(d/(b*n*(p + S(1))), Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sinh(c + d*x), x), x) - Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*(a + b*x**n)**(p + S(1))*cosh(c + d*x), x), x) + Simp(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*cosh(c + d*x)/(b*n*(p + S(1))), x)
def replacement5704(c, m, n, x, d, a, p, b):
        # rubi.append(5704)
        return Int(ExpandIntegrand(sinh(c + d*x), x**m*(a + b*x**n)**p, x), x)
def replacement5705(c, m, n, x, d, a, p, b):
        # rubi.append(5705)
        return Int(ExpandIntegrand(cosh(c + d*x), x**m*(a + b*x**n)**p, x), x)
def replacement5706(c, m, n, x, d, a, p, b):
        # rubi.append(5706)
        return Int(x**(m + n*p)*(a*x**(-n) + b)**p*sinh(c + d*x), x)
def replacement5707(c, m, n, x, d, a, p, b):
        # rubi.append(5707)
        return Int(x**(m + n*p)*(a*x**(-n) + b)**p*cosh(c + d*x), x)
def replacement5708(c, m, n, x, d, a, p, e, b):
        # rubi.append(5708)
        return Int((e*x)**m*(a + b*x**n)**p*sinh(c + d*x), x)
def replacement5709(c, m, n, x, d, a, p, e, b):
        # rubi.append(5709)
        return Int((e*x)**m*(a + b*x**n)**p*cosh(c + d*x), x)
def replacement5710(c, d, x, n):
        # rubi.append(5710)
        return -Dist(S(1)/2, Int(exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int(exp(c + d*x**n), x), x)
def replacement5711(c, d, x, n):
        # rubi.append(5711)
        return Dist(S(1)/2, Int(exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int(exp(c + d*x**n), x), x)
def replacement5712(c, n, x, d, a, p, b):
        # rubi.append(5712)
        return Int(ExpandTrigReduce((a + b*sinh(c + d*x**n))**p, x), x)
def replacement5713(c, n, x, d, a, p, b):
        # rubi.append(5713)
        return Int(ExpandTrigReduce((a + b*cosh(c + d*x**n))**p, x), x)
def replacement5714(c, n, x, d, a, p, b):
        # rubi.append(5714)
        return -Subst(Int((a + b*sinh(c + d*x**(-n)))**p/x**S(2), x), x, S(1)/x)
def replacement5715(c, n, x, d, a, p, b):
        # rubi.append(5715)
        return -Subst(Int((a + b*cosh(c + d*x**(-n)))**p/x**S(2), x), x, S(1)/x)

def With5716(c, n, x, d, a, p, b):
        k = Denominator(n)
        # rubi.append(5716)
        return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*sinh(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)

def With5717(c, n, x, d, a, p, b):
        k = Denominator(n)
        # rubi.append(5717)
        return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*cosh(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)
def replacement5718(c, d, x, n):
        # rubi.append(5718)
        return -Dist(S(1)/2, Int(exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int(exp(c + d*x**n), x), x)
def replacement5719(c, d, x, n):
        # rubi.append(5719)
        return Dist(S(1)/2, Int(exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int(exp(c + d*x**n), x), x)
def replacement5720(c, n, x, d, a, p, b):
        # rubi.append(5720)
        return Int(ExpandTrigReduce((a + b*sinh(c + d*x**n))**p, x), x)
def replacement5721(c, n, x, d, a, p, b):
        # rubi.append(5721)
        return Int(ExpandTrigReduce((a + b*cosh(c + d*x**n))**p, x), x)
def replacement5722(c, n, x, d, u, a, p, b):
        # rubi.append(5722)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*sinh(c + d*x**n))**p, x), x, u), x)
def replacement5723(c, n, x, d, u, a, p, b):
        # rubi.append(5723)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*cosh(c + d*x**n))**p, x), x, u), x)
def replacement5724(c, n, x, d, u, a, p, b):
        # rubi.append(5724)
        return Int((a + b*sinh(c + d*u**n))**p, x)
def replacement5725(c, n, x, d, u, a, p, b):
        # rubi.append(5725)
        return Int((a + b*cosh(c + d*u**n))**p, x)
def replacement5726(x, u, a, p, b):
        # rubi.append(5726)
        return Int((a + b*sinh(ExpandToSum(u, x)))**p, x)
def replacement5727(x, u, a, p, b):
        # rubi.append(5727)
        return Int((a + b*cosh(ExpandToSum(u, x)))**p, x)
def replacement5728(x, d, n):
        # rubi.append(5728)
        return Simp(SinhIntegral(d*x**n)/n, x)
def replacement5729(x, d, n):
        # rubi.append(5729)
        return Simp(CoshIntegral(d*x**n)/n, x)
def replacement5730(x, d, c, n):
        # rubi.append(5730)
        return Dist(sinh(c), Int(cosh(d*x**n)/x, x), x) + Dist(cosh(c), Int(sinh(d*x**n)/x, x), x)
def replacement5731(c, d, x, n):
        # rubi.append(5731)
        return Dist(sinh(c), Int(sinh(d*x**n)/x, x), x) + Dist(cosh(c), Int(cosh(d*x**n)/x, x), x)

def With5732(c, m, n, x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        mn = (m + S(1))/n
        if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
            return True
        return False
def replacement5732(c, m, n, x, d, a, p, b):

        mn = (m + S(1))/n
        # rubi.append(5732)
        return Dist(S(1)/n, Subst(Int(x**(mn + S(-1))*(a + b*sinh(c + d*x))**p, x), x, x**n), x)

def With5733(c, m, n, x, d, a, p, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        mn = (m + S(1))/n
        if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
            return True
        return False
def replacement5733(c, m, n, x, d, a, p, b):

        mn = (m + S(1))/n
        # rubi.append(5733)
        return Dist(S(1)/n, Subst(Int(x**(mn + S(-1))*(a + b*cosh(c + d*x))**p, x), x, x**n), x)

def With5734(c, m, n, x, d, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        mn = (m + S(1))/n
        if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
            return True
        return False
def replacement5734(c, m, n, x, d, a, p, e, b):

        mn = (m + S(1))/n
        # rubi.append(5734)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sinh(c + d*x**n))**p, x), x)

def With5735(c, m, n, x, d, a, p, e, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        mn = (m + S(1))/n
        if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
            return True
        return False
def replacement5735(c, m, n, x, d, a, p, e, b):

        mn = (m + S(1))/n
        # rubi.append(5735)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cosh(c + d*x**n))**p, x), x)
def replacement5736(c, m, n, x, d, e):
        # rubi.append(5736)
        return -Dist(e**n*(m - n + S(1))/(d*n), Int((e*x)**(m - n)*cosh(c + d*x**n), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*cosh(c + d*x**n)/(d*n), x)
def replacement5737(c, m, n, x, d, e):
        # rubi.append(5737)
        return -Dist(e**n*(m - n + S(1))/(d*n), Int((e*x)**(m - n)*sinh(c + d*x**n), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*sinh(c + d*x**n)/(d*n), x)
def replacement5738(c, m, n, x, d, e):
        # rubi.append(5738)
        return -Dist(d*e**(-n)*n/(m + S(1)), Int((e*x)**(m + n)*cosh(c + d*x**n), x), x) + Simp((e*x)**(m + S(1))*sinh(c + d*x**n)/(e*(m + S(1))), x)
def replacement5739(c, m, n, x, d, e):
        # rubi.append(5739)
        return -Dist(d*e**(-n)*n/(m + S(1)), Int((e*x)**(m + n)*sinh(c + d*x**n), x), x) + Simp((e*x)**(m + S(1))*cosh(c + d*x**n)/(e*(m + S(1))), x)
def replacement5740(c, m, n, x, d, e):
        # rubi.append(5740)
        return -Dist(S(1)/2, Int((e*x)**m*exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(c + d*x**n), x), x)
def replacement5741(c, m, n, x, d, e):
        # rubi.append(5741)
        return Dist(S(1)/2, Int((e*x)**m*exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(c + d*x**n), x), x)
def replacement5742(m, n, x, a, p, b):
        # rubi.append(5742)
        return Dist(b*n*p/(n + S(-1)), Int(sinh(a + b*x**n)**(p + S(-1))*cosh(a + b*x**n), x), x) - Simp(x**(S(1) - n)*sinh(a + b*x**n)**p/(n + S(-1)), x)
def replacement5743(m, n, x, a, p, b):
        # rubi.append(5743)
        return Dist(b*n*p/(n + S(-1)), Int(sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(-1)), x), x) - Simp(x**(S(1) - n)*cosh(a + b*x**n)**p/(n + S(-1)), x)
def replacement5744(m, n, x, a, p, b):
        # rubi.append(5744)
        return -Dist((p + S(-1))/p, Int(x**m*sinh(a + b*x**n)**(p + S(-2)), x), x) - Simp(sinh(a + b*x**n)**p/(b**S(2)*n*p**S(2)), x) + Simp(x**n*sinh(a + b*x**n)**(p + S(-1))*cosh(a + b*x**n)/(b*n*p), x)
def replacement5745(m, n, x, a, p, b):
        # rubi.append(5745)
        return Dist((p + S(-1))/p, Int(x**m*cosh(a + b*x**n)**(p + S(-2)), x), x) - Simp(cosh(a + b*x**n)**p/(b**S(2)*n*p**S(2)), x) + Simp(x**n*sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(-1))/(b*n*p), x)
def replacement5746(m, n, x, a, p, b):
        # rubi.append(5746)
        return -Dist((p + S(-1))/p, Int(x**m*sinh(a + b*x**n)**(p + S(-2)), x), x) + Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*p**S(2)), Int(x**(m - S(2)*n)*sinh(a + b*x**n)**p, x), x) - Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*sinh(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)), x) + Simp(x**(m - n + S(1))*sinh(a + b*x**n)**(p + S(-1))*cosh(a + b*x**n)/(b*n*p), x)
def replacement5747(m, n, x, a, p, b):
        # rubi.append(5747)
        return Dist((p + S(-1))/p, Int(x**m*cosh(a + b*x**n)**(p + S(-2)), x), x) + Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*p**S(2)), Int(x**(m - S(2)*n)*cosh(a + b*x**n)**p, x), x) - Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*cosh(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)), x) + Simp(x**(m - n + S(1))*sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(-1))/(b*n*p), x)
def replacement5748(m, n, x, a, p, b):
        # rubi.append(5748)
        return Dist(b**S(2)*n**S(2)*p**S(2)/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*sinh(a + b*x**n)**p, x), x) + Dist(b**S(2)*n**S(2)*p*(p + S(-1))/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*sinh(a + b*x**n)**(p + S(-2)), x), x) + Simp(x**(m + S(1))*sinh(a + b*x**n)**p/(m + S(1)), x) - Simp(b*n*p*x**(m + n + S(1))*sinh(a + b*x**n)**(p + S(-1))*cosh(a + b*x**n)/((m + S(1))*(m + n + S(1))), x)
def replacement5749(m, n, x, a, p, b):
        # rubi.append(5749)
        return Dist(b**S(2)*n**S(2)*p**S(2)/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*cosh(a + b*x**n)**p, x), x) - Dist(b**S(2)*n**S(2)*p*(p + S(-1))/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*cosh(a + b*x**n)**(p + S(-2)), x), x) + Simp(x**(m + S(1))*cosh(a + b*x**n)**p/(m + S(1)), x) - Simp(b*n*p*x**(m + n + S(1))*sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(-1))/((m + S(1))*(m + n + S(1))), x)

def With5750(c, m, n, x, d, a, p, e, b):
        k = Denominator(m)
        # rubi.append(5750)
        return Dist(k/e, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(S(1)/k)), x)

def With5751(c, m, n, x, d, a, p, e, b):
        k = Denominator(m)
        # rubi.append(5751)
        return Dist(k/e, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*cosh(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(S(1)/k)), x)
def replacement5752(c, m, n, x, d, a, p, e, b):
        # rubi.append(5752)
        return Int(ExpandTrigReduce((e*x)**m, (a + b*sinh(c + d*x**n))**p, x), x)
def replacement5753(c, m, n, x, d, a, p, e, b):
        # rubi.append(5753)
        return Int(ExpandTrigReduce((e*x)**m, (a + b*cosh(c + d*x**n))**p, x), x)
def replacement5754(m, n, x, a, p, b):
        # rubi.append(5754)
        return -Dist((p + S(2))/(p + S(1)), Int(x**m*sinh(a + b*x**n)**(p + S(2)), x), x) - Simp(sinh(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))), x) + Simp(x**n*sinh(a + b*x**n)**(p + S(1))*cosh(a + b*x**n)/(b*n*(p + S(1))), x)
def replacement5755(m, n, x, a, p, b):
        # rubi.append(5755)
        return Dist((p + S(2))/(p + S(1)), Int(x**m*cosh(a + b*x**n)**(p + S(2)), x), x) + Simp(cosh(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))), x) - Simp(x**n*sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement5756(m, n, x, a, p, b):
        # rubi.append(5756)
        return -Dist((p + S(2))/(p + S(1)), Int(x**m*sinh(a + b*x**n)**(p + S(2)), x), x) + Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**(m - S(2)*n)*sinh(a + b*x**n)**(p + S(2)), x), x) + Simp(x**(m - n + S(1))*sinh(a + b*x**n)**(p + S(1))*cosh(a + b*x**n)/(b*n*(p + S(1))), x) - Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*sinh(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement5757(m, n, x, a, p, b):
        # rubi.append(5757)
        return Dist((p + S(2))/(p + S(1)), Int(x**m*cosh(a + b*x**n)**(p + S(2)), x), x) - Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**(m - S(2)*n)*cosh(a + b*x**n)**(p + S(2)), x), x) - Simp(x**(m - n + S(1))*sinh(a + b*x**n)*cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x) + Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*cosh(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement5758(c, m, n, x, d, a, p, b):
        # rubi.append(5758)
        return -Subst(Int(x**(-m + S(-2))*(a + b*sinh(c + d*x**(-n)))**p, x), x, S(1)/x)
def replacement5759(c, m, n, x, d, a, p, b):
        # rubi.append(5759)
        return -Subst(Int(x**(-m + S(-2))*(a + b*cosh(c + d*x**(-n)))**p, x), x, S(1)/x)

def With5760(c, m, n, x, d, a, p, e, b):
        k = Denominator(m)
        # rubi.append(5760)
        return -Dist(k/e, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k)), x)

def With5761(c, m, n, x, d, a, p, e, b):
        k = Denominator(m)
        # rubi.append(5761)
        return -Dist(k/e, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*cosh(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k)), x)
def replacement5762(c, m, n, x, d, a, p, e, b):
        # rubi.append(5762)
        return -Dist((e*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*sinh(c + d*x**(-n)))**p, x), x, S(1)/x), x)
def replacement5763(c, m, n, x, d, a, p, e, b):
        # rubi.append(5763)
        return -Dist((e*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*cosh(c + d*x**(-n)))**p, x), x, S(1)/x), x)

def With5764(c, m, n, x, d, a, p, b):
        k = Denominator(n)
        # rubi.append(5764)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sinh(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)

def With5765(c, m, n, x, d, a, p, b):
        k = Denominator(n)
        # rubi.append(5765)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*cosh(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)
def replacement5766(c, m, n, x, d, a, p, e, b):
        # rubi.append(5766)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sinh(c + d*x**n))**p, x), x)
def replacement5767(c, m, n, x, d, a, p, e, b):
        # rubi.append(5767)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cosh(c + d*x**n))**p, x), x)
def replacement5768(c, m, n, x, d, a, p, b):
        # rubi.append(5768)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*sinh(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1))), x)
def replacement5769(c, m, n, x, d, a, p, b):
        # rubi.append(5769)
        return Dist(S(1)/(m + S(1)), Subst(Int((a + b*cosh(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1))), x)
def replacement5770(c, m, n, x, d, a, p, e, b):
        # rubi.append(5770)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sinh(c + d*x**n))**p, x), x)
def replacement5771(c, m, n, x, d, a, p, e, b):
        # rubi.append(5771)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cosh(c + d*x**n))**p, x), x)
def replacement5772(c, m, n, x, d, e):
        # rubi.append(5772)
        return -Dist(S(1)/2, Int((e*x)**m*exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(c + d*x**n), x), x)
def replacement5773(c, m, n, x, d, e):
        # rubi.append(5773)
        return Dist(S(1)/2, Int((e*x)**m*exp(-c - d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(c + d*x**n), x), x)
def replacement5774(c, m, n, x, d, a, p, e, b):
        # rubi.append(5774)
        return Int(ExpandTrigReduce((e*x)**m, (a + b*sinh(c + d*x**n))**p, x), x)
def replacement5775(c, m, n, x, d, a, p, e, b):
        # rubi.append(5775)
        return Int(ExpandTrigReduce((e*x)**m, (a + b*cosh(c + d*x**n))**p, x), x)
def replacement5776(c, m, n, x, d, u, a, p, b):
        # rubi.append(5776)
        return Dist(Coefficient(u, x, S(1))**(-m + S(-1)), Subst(Int((a + b*sinh(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u), x)
def replacement5777(c, m, n, x, d, u, a, p, b):
        # rubi.append(5777)
        return Dist(Coefficient(u, x, S(1))**(-m + S(-1)), Subst(Int((a + b*cosh(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u), x)
def replacement5778(c, m, n, x, d, u, a, p, e, b):
        # rubi.append(5778)
        return Int((e*x)**m*(a + b*sinh(c + d*u**n))**p, x)
def replacement5779(c, m, n, x, d, u, a, p, e, b):
        # rubi.append(5779)
        return Int((e*x)**m*(a + b*cosh(c + d*u**n))**p, x)
def replacement5780(m, x, u, a, p, e, b):
        # rubi.append(5780)
        return Int((e*x)**m*(a + b*sinh(ExpandToSum(u, x)))**p, x)
def replacement5781(m, x, u, a, p, e, b):
        # rubi.append(5781)
        return Int((e*x)**m*(a + b*cosh(ExpandToSum(u, x)))**p, x)
def replacement5782(m, n, x, a, p, b):
        # rubi.append(5782)
        return Simp(sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement5783(m, n, x, a, p, b):
        # rubi.append(5783)
        return Simp(cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement5784(m, n, x, a, p, b):
        # rubi.append(5784)
        return -Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*sinh(a + b*x**n)**(p + S(1)), x), x) + Simp(x**(m - n + S(1))*sinh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement5785(m, n, x, a, p, b):
        # rubi.append(5785)
        return -Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*cosh(a + b*x**n)**(p + S(1)), x), x) + Simp(x**(m - n + S(1))*cosh(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)
def replacement5786(c, a, x, b):
        # rubi.append(5786)
        return -Dist(S(1)/2, Int(exp(-a - b*x - c*x**S(2)), x), x) + Dist(S(1)/2, Int(exp(a + b*x + c*x**S(2)), x), x)
def replacement5787(c, a, x, b):
        # rubi.append(5787)
        return Dist(S(1)/2, Int(exp(-a - b*x - c*x**S(2)), x), x) + Dist(S(1)/2, Int(exp(a + b*x + c*x**S(2)), x), x)
def replacement5788(c, n, x, a, b):
        # rubi.append(5788)
        return Int(ExpandTrigReduce(sinh(a + b*x + c*x**S(2))**n, x), x)
def replacement5789(c, n, x, a, b):
        # rubi.append(5789)
        return Int(ExpandTrigReduce(cosh(a + b*x + c*x**S(2))**n, x), x)
def replacement5790(x, v, n):
        # rubi.append(5790)
        return Int(sinh(ExpandToSum(v, x))**n, x)
def replacement5791(x, v, n):
        # rubi.append(5791)
        return Int(cosh(ExpandToSum(v, x))**n, x)
def replacement5792(c, x, d, a, e, b):
        # rubi.append(5792)
        return Simp(e*cosh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5793(c, x, d, a, e, b):
        # rubi.append(5793)
        return Simp(e*sinh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5794(c, x, d, a, e, b):
        # rubi.append(5794)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(sinh(a + b*x + c*x**S(2)), x), x) + Simp(e*cosh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5795(c, x, d, a, e, b):
        # rubi.append(5795)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(cosh(a + b*x + c*x**S(2)), x), x) + Simp(e*sinh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5796(c, m, x, d, a, e, b):
        # rubi.append(5796)
        return -Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*cosh(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*cosh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5797(c, m, x, d, a, e, b):
        # rubi.append(5797)
        return -Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*sinh(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5798(c, m, x, d, a, e, b):
        # rubi.append(5798)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int((d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2)), x), x) - Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*cosh(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*cosh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5799(c, m, x, d, a, e, b):
        # rubi.append(5799)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int((d + e*x)**(m + S(-1))*cosh(a + b*x + c*x**S(2)), x), x) - Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*sinh(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*sinh(a + b*x + c*x**S(2))/(S(2)*c), x)
def replacement5800(c, m, x, d, a, e, b):
        # rubi.append(5800)
        return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*cosh(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2))/(e*(m + S(1))), x)
def replacement5801(c, m, x, d, a, e, b):
        # rubi.append(5801)
        return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*sinh(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*cosh(a + b*x + c*x**S(2))/(e*(m + S(1))), x)
def replacement5802(c, m, x, d, a, e, b):
        # rubi.append(5802)
        return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*cosh(a + b*x + c*x**S(2)), x), x) - Dist((b*e - S(2)*c*d)/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(1))*cosh(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2))/(e*(m + S(1))), x)
def replacement5803(c, m, x, d, a, e, b):
        # rubi.append(5803)
        return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*sinh(a + b*x + c*x**S(2)), x), x) - Dist((b*e - S(2)*c*d)/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(1))*sinh(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*cosh(a + b*x + c*x**S(2))/(e*(m + S(1))), x)
def replacement5804(c, m, x, d, a, e, b):
        # rubi.append(5804)
        return Int((d + e*x)**m*sinh(a + b*x + c*x**S(2)), x)
def replacement5805(c, m, x, d, a, e, b):
        # rubi.append(5805)
        return Int((d + e*x)**m*cosh(a + b*x + c*x**S(2)), x)
def replacement5806(c, m, n, x, d, a, e, b):
        # rubi.append(5806)
        return Int(ExpandTrigReduce((d + e*x)**m, sinh(a + b*x + c*x**S(2))**n, x), x)
def replacement5807(c, m, n, x, d, a, e, b):
        # rubi.append(5807)
        return Int(ExpandTrigReduce((d + e*x)**m, cosh(a + b*x + c*x**S(2))**n, x), x)
def replacement5808(m, n, x, v, u):
        # rubi.append(5808)
        return Int(ExpandToSum(u, x)**m*sinh(ExpandToSum(v, x))**n, x)
def replacement5809(m, n, x, v, u):
        # rubi.append(5809)
        return Int(ExpandToSum(u, x)**m*cosh(ExpandToSum(v, x))**n, x)
def replacement5810(m, x, f, a, e, b):
        # rubi.append(5810)
        return Dist(S(2), Int((a + b*x)**m*exp(S(2)*e + S(2)*f*x)/(exp(S(2)*e + S(2)*f*x) + S(1)), x), x) - Simp((a + b*x)**(m + S(1))/(b*(m + S(1))), x)
def replacement5811(m, x, f, a, e, b):
        # rubi.append(5811)
        return -Dist(S(2), Int((a + b*x)**m*exp(S(2)*e + S(2)*f*x)/(S(1) - exp(S(2)*e + S(2)*f*x)), x), x) - Simp((a + b*x)**(m + S(1))/(b*(m + S(1))), x)
def replacement5812(c, m, n, x, f, a, e, b):
        # rubi.append(5812)
        return Dist(c**S(2), Int((c*tanh(e + f*x))**(n + S(-2))*(a + b*x)**m, x), x) + Dist(b*c*m/(f*(n + S(-1))), Int((c*tanh(e + f*x))**(n + S(-1))*(a + b*x)**(m + S(-1)), x), x) - Simp(c*(c*tanh(e + f*x))**(n + S(-1))*(a + b*x)**m/(f*(n + S(-1))), x)
def replacement5813(c, m, n, x, f, a, e, b):
        # rubi.append(5813)
        return Dist(c**S(2), Int((c/tanh(e + f*x))**(n + S(-2))*(a + b*x)**m, x), x) + Dist(b*c*m/(f*(n + S(-1))), Int((c/tanh(e + f*x))**(n + S(-1))*(a + b*x)**(m + S(-1)), x), x) - Simp(c*(c/tanh(e + f*x))**(n + S(-1))*(a + b*x)**m/(f*(n + S(-1))), x)
def replacement5814(c, m, n, x, f, a, e, b):
        # rubi.append(5814)
        return Dist(c**(S(-2)), Int((c*tanh(e + f*x))**(n + S(2))*(a + b*x)**m, x), x) - Dist(b*m/(c*f*(n + S(1))), Int((c*tanh(e + f*x))**(n + S(1))*(a + b*x)**(m + S(-1)), x), x) + Simp((c*tanh(e + f*x))**(n + S(1))*(a + b*x)**m/(c*f*(n + S(1))), x)
def replacement5815(c, m, n, x, f, a, e, b):
        # rubi.append(5815)
        return Dist(c**(S(-2)), Int((c/tanh(e + f*x))**(n + S(2))*(a + b*x)**m, x), x) - Dist(b*m/(c*f*(n + S(1))), Int((c/tanh(e + f*x))**(n + S(1))*(a + b*x)**(m + S(-1)), x), x) + Simp((c/tanh(e + f*x))**(n + S(1))*(a + b*x)**m/(c*f*(n + S(1))), x)
def replacement5816(c, m, n, x, d, f, a, e, b):
        # rubi.append(5816)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b*tanh(e + f*x))**n, x), x)
def replacement5817(c, m, n, x, d, f, a, e, b):
        # rubi.append(5817)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b/tanh(e + f*x))**n, x), x)
def replacement5818(c, m, x, d, f, a, e, b):
        # rubi.append(5818)
        return Dist(a*d*m/(S(2)*b*f), Int((c + d*x)**(m + S(-1))/(a + b*tanh(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x) - Simp(a*(c + d*x)**m/(S(2)*b*f*(a + b*tanh(e + f*x))), x)
def replacement5819(c, m, x, d, f, a, e, b):
        # rubi.append(5819)
        return Dist(a*d*m/(S(2)*b*f), Int((c + d*x)**(m + S(-1))/(a + b/tanh(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x) - Simp(a*(c + d*x)**m/(S(2)*b*f*(a + b/tanh(e + f*x))), x)
def replacement5820(c, x, d, f, a, e, b):
        # rubi.append(5820)
        return Dist(f/(a*d), Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Dist(f/(b*d), Int(cosh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Simp(S(1)/(d*(a + b*tanh(e + f*x))*(c + d*x)), x)
def replacement5821(c, x, d, f, a, e, b):
        # rubi.append(5821)
        return -Dist(f/(a*d), Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(f/(b*d), Int(cosh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Simp(S(1)/(d*(a + b/tanh(e + f*x))*(c + d*x)), x)
def replacement5822(c, m, x, d, f, a, e, b):
        # rubi.append(5822)
        return Dist(S(2)*b*f/(a*d*(m + S(1))), Int((c + d*x)**(m + S(1))/(a + b*tanh(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + b*tanh(e + f*x))*(m + S(1))), x) - Simp(f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement5823(c, m, x, d, f, a, e, b):
        # rubi.append(5823)
        return Dist(S(2)*b*f/(a*d*(m + S(1))), Int((c + d*x)**(m + S(1))/(a + b/tanh(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + b/tanh(e + f*x))*(m + S(1))), x) - Simp(f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement5824(c, x, d, f, a, e, b):
        # rubi.append(5824)
        return Dist(S(1)/(S(2)*a), Int(cosh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Dist(S(1)/(S(2)*b), Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Simp(log(c + d*x)/(S(2)*a*d), x)
def replacement5825(c, x, d, f, a, e, b):
        # rubi.append(5825)
        return -Dist(S(1)/(S(2)*a), Int(cosh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(S(1)/(S(2)*b), Int(sinh(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Simp(log(c + d*x)/(S(2)*a*d), x)
def replacement5826(c, m, x, d, f, a, e, b):
        # rubi.append(5826)
        return Dist(S(1)/(S(2)*a), Int((c + d*x)**m*exp(-S(2)*a*(e + f*x)/b), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x)
def replacement5827(c, m, x, d, f, a, e, b):
        # rubi.append(5827)
        return -Dist(S(1)/(S(2)*a), Int((c + d*x)**m*exp(-S(2)*a*(e + f*x)/b), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x)
def replacement5828(c, m, n, x, d, f, a, e, b):
        # rubi.append(5828)
        return Int(ExpandIntegrand((c + d*x)**m, (-sinh(S(2)*e + S(2)*f*x)/(S(2)*b) + cosh(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x)
def replacement5829(c, m, n, x, d, f, a, e, b):
        # rubi.append(5829)
        return Int(ExpandIntegrand((c + d*x)**m, (sinh(S(2)*e + S(2)*f*x)/(S(2)*b) - cosh(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x)
def replacement5830(c, m, n, x, d, f, a, e, b):
        # rubi.append(5830)
        return Int(ExpandIntegrand((c + d*x)**m, (S(1)/(S(2)*a) + exp(-S(2)*a*(e + f*x)/b)/(S(2)*a))**(-n), x), x)
def replacement5831(c, m, n, x, d, f, a, e, b):
        # rubi.append(5831)
        return Int(ExpandIntegrand((c + d*x)**m, (S(1)/(S(2)*a) - exp(-S(2)*a*(e + f*x)/b)/(S(2)*a))**(-n), x), x)

def With5832(c, m, n, x, d, f, a, e, b):
        u = IntHide((a + b*tanh(e + f*x))**n, x)
        # rubi.append(5832)
        return -Dist(d*m, Int(Dist((c + d*x)**(m + S(-1)), u, x), x), x) + Dist((c + d*x)**m, u, x)

def With5833(c, m, n, x, d, f, a, e, b):
        u = IntHide((a + b/tanh(e + f*x))**n, x)
        # rubi.append(5833)
        return -Dist(d*m, Int(Dist((c + d*x)**(m + S(-1)), u, x), x), x) + Dist((c + d*x)**m, u, x)
def replacement5834(c, m, x, d, f, a, e, b):
        # rubi.append(5834)
        return -Dist(S(2)*b, Int((c + d*x)**m/(a**S(2) - b**S(2) + (a - b)**S(2)*exp(-S(2)*e - S(2)*f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a - b)*(m + S(1))), x)
def replacement5835(c, m, x, d, f, a, e, b):
        # rubi.append(5835)
        return Dist(S(2)*b, Int((c + d*x)**m/(a**S(2) - b**S(2) - (a + b)**S(2)*exp(S(2)*e + S(2)*f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + b)*(m + S(1))), x)
def replacement5836(c, x, d, f, a, e, b):
        # rubi.append(5836)
        return -Dist(S(1)/(f*(a**S(2) - b**S(2))), Int((-S(2)*a*c*f - S(2)*a*d*f*x + b*d)/(a + b*tanh(e + f*x)), x), x) - Simp((c + d*x)**S(2)/(S(2)*d*(a**S(2) - b**S(2))), x) + Simp(b*(c + d*x)/(f*(a + b*tanh(e + f*x))*(a**S(2) - b**S(2))), x)
def replacement5837(c, x, d, f, a, e, b):
        # rubi.append(5837)
        return -Dist(S(1)/(f*(a**S(2) - b**S(2))), Int((-S(2)*a*c*f - S(2)*a*d*f*x + b*d)/(a + b/tanh(e + f*x)), x), x) - Simp((c + d*x)**S(2)/(S(2)*d*(a**S(2) - b**S(2))), x) + Simp(b*(c + d*x)/(f*(a + b/tanh(e + f*x))*(a**S(2) - b**S(2))), x)
def replacement5838(c, m, n, x, d, f, a, e, b):
        # rubi.append(5838)
        return Int(ExpandIntegrand((c + d*x)**m, (-S(2)*b/(a**S(2) - b**S(2) + (a - b)**S(2)*exp(-S(2)*e - S(2)*f*x)) + S(1)/(a - b))**(-n), x), x)
def replacement5839(c, m, n, x, d, f, a, e, b):
        # rubi.append(5839)
        return Int(ExpandIntegrand((c + d*x)**m, (S(2)*b/(a**S(2) - b**S(2) - (a + b)**S(2)*exp(S(2)*e + S(2)*f*x)) + S(1)/(a + b))**(-n), x), x)
def replacement5840(m, n, x, v, u, a, b):
        # rubi.append(5840)
        return Int((a + b*tanh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5841(m, n, x, v, u, a, b):
        # rubi.append(5841)
        return Int((a + b/tanh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5842(c, m, n, x, d, f, a, e, b):
        # rubi.append(5842)
        return Int((a + b*tanh(e + f*x))**n*(c + d*x)**m, x)
def replacement5843(c, m, n, x, d, f, a, e, b):
        # rubi.append(5843)
        return Int((a + b/tanh(e + f*x))**n*(c + d*x)**m, x)
def replacement5844(c, n, x, d, a, p, b):
        # rubi.append(5844)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b*tanh(c + d*x))**p, x), x, x**n), x)
def replacement5845(c, n, x, d, a, p, b):
        # rubi.append(5845)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/tanh(c + d*x))**p, x), x, x**n), x)
def replacement5846(c, n, x, d, a, p, b):
        # rubi.append(5846)
        return Int((a + b*tanh(c + d*x**n))**p, x)
def replacement5847(c, n, x, d, a, p, b):
        # rubi.append(5847)
        return Int((a + b/tanh(c + d*x**n))**p, x)
def replacement5848(c, n, x, d, u, a, p, b):
        # rubi.append(5848)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*tanh(c + d*x**n))**p, x), x, u), x)
def replacement5849(c, n, x, d, u, a, p, b):
        # rubi.append(5849)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/tanh(c + d*x**n))**p, x), x, u), x)
def replacement5850(x, u, a, p, b):
        # rubi.append(5850)
        return Int((a + b*tanh(ExpandToSum(u, x)))**p, x)
def replacement5851(x, u, a, p, b):
        # rubi.append(5851)
        return Int((a + b/tanh(ExpandToSum(u, x)))**p, x)
def replacement5852(c, m, n, x, d, a, p, b):
        # rubi.append(5852)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*tanh(c + d*x))**p, x), x, x**n), x)
def replacement5853(c, m, n, x, d, a, p, b):
        # rubi.append(5853)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/tanh(c + d*x))**p, x), x, x**n), x)
def replacement5854(c, m, n, x, d):
        # rubi.append(5854)
        return Dist((m - n + S(1))/(d*n), Int(x**(m - n)*tanh(c + d*x**n), x), x) + Int(x**m, x) - Simp(x**(m - n + S(1))*tanh(c + d*x**n)/(d*n), x)
def replacement5855(c, m, n, x, d):
        # rubi.append(5855)
        return Dist((m - n + S(1))/(d*n), Int(x**(m - n)/tanh(c + d*x**n), x), x) + Int(x**m, x) - Simp(x**(m - n + S(1))/(d*n*tanh(c + d*x**n)), x)
def replacement5856(c, m, n, x, d, a, p, b):
        # rubi.append(5856)
        return Int(x**m*(a + b*tanh(c + d*x**n))**p, x)
def replacement5857(c, m, n, x, d, a, p, b):
        # rubi.append(5857)
        return Int(x**m*(a + b/tanh(c + d*x**n))**p, x)
def replacement5858(c, m, n, x, d, a, p, e, b):
        # rubi.append(5858)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*tanh(c + d*x**n))**p, x), x)
def replacement5859(c, m, n, x, d, a, p, e, b):
        # rubi.append(5859)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/tanh(c + d*x**n))**p, x), x)
def replacement5860(m, x, u, a, p, e, b):
        # rubi.append(5860)
        return Int((e*x)**m*(a + b*tanh(ExpandToSum(u, x)))**p, x)
def replacement5861(m, x, u, a, p, e, b):
        # rubi.append(5861)
        return Int((e*x)**m*(a + b/tanh(ExpandToSum(u, x)))**p, x)
def replacement5862(m, q, n, x, a, p, b):
        # rubi.append(5862)
        return Dist((m - n + S(1))/(b*n*p), Int(x**(m - n)*(S(1)/cosh(a + b*x**n))**p, x), x) - Simp(x**(m - n + S(1))*(S(1)/cosh(a + b*x**n))**p/(b*n*p), x)
def replacement5863(m, q, n, x, a, p, b):
        # rubi.append(5863)
        return Dist((m - n + S(1))/(b*n*p), Int(x**(m - n)*(S(1)/sinh(a + b*x**n))**p, x), x) - Simp(x**(m - n + S(1))*(S(1)/sinh(a + b*x**n))**p/(b*n*p), x)
def replacement5864(c, n, x, a, b):
        # rubi.append(5864)
        return Int(tanh(a + b*x + c*x**S(2))**n, x)
def replacement5865(c, n, x, a, b):
        # rubi.append(5865)
        return Int((S(1)/tanh(a + b*x + c*x**S(2)))**n, x)
def replacement5866(c, x, d, a, e, b):
        # rubi.append(5866)
        return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(tanh(a + b*x + c*x**S(2)), x), x) + Simp(e*log(cosh(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement5867(c, x, d, a, e, b):
        # rubi.append(5867)
        return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(S(1)/tanh(a + b*x + c*x**S(2)), x), x) + Simp(e*log(sinh(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement5868(c, m, n, x, d, a, e, b):
        # rubi.append(5868)
        return Int((d + e*x)**m*tanh(a + b*x + c*x**S(2))**n, x)
def replacement5869(c, m, n, x, d, a, e, b):
        # rubi.append(5869)
        return Int((d + e*x)**m*(S(1)/tanh(a + b*x + c*x**S(2)))**n, x)
def replacement5870(c, m, x, d, a, b):
        # rubi.append(5870)
        return -Dist(I*d*m/b, Int((c + d*x)**(m + S(-1))*log(-I*exp(a + b*x) + S(1)), x), x) + Dist(I*d*m/b, Int((c + d*x)**(m + S(-1))*log(I*exp(a + b*x) + S(1)), x), x) + Simp(S(2)*(c + d*x)**m*ArcTan(exp(a + b*x))/b, x)
def replacement5871(c, m, x, d, a, b):
        # rubi.append(5871)
        return -Dist(d*m/b, Int((c + d*x)**(m + S(-1))*log(S(1) - exp(a + b*x)), x), x) + Dist(d*m/b, Int((c + d*x)**(m + S(-1))*log(exp(a + b*x) + S(1)), x), x) + Simp(-S(2)*(c + d*x)**m*atanh(exp(a + b*x))/b, x)
def replacement5872(c, m, x, d, a, b):
        # rubi.append(5872)
        return -Dist(d*m/b, Int((c + d*x)**(m + S(-1))*tanh(a + b*x), x), x) + Simp((c + d*x)**m*tanh(a + b*x)/b, x)
def replacement5873(c, m, x, d, a, b):
        # rubi.append(5873)
        return Dist(d*m/b, Int((c + d*x)**(m + S(-1))/tanh(a + b*x), x), x) - Simp((c + d*x)**m/(b*tanh(a + b*x)), x)
def replacement5874(c, n, x, d, a, b):
        # rubi.append(5874)
        return Dist((n + S(-2))/(n + S(-1)), Int((c + d*x)*(S(1)/cosh(a + b*x))**(n + S(-2)), x), x) + Simp(d*(S(1)/cosh(a + b*x))**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))), x) + Simp((c + d*x)*(S(1)/cosh(a + b*x))**(n + S(-2))*tanh(a + b*x)/(b*(n + S(-1))), x)
def replacement5875(c, n, x, d, a, b):
        # rubi.append(5875)
        return -Dist((n + S(-2))/(n + S(-1)), Int((c + d*x)*(S(1)/sinh(a + b*x))**(n + S(-2)), x), x) - Simp(d*(S(1)/sinh(a + b*x))**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))), x) - Simp((c + d*x)*(S(1)/sinh(a + b*x))**(n + S(-2))/(b*(n + S(-1))*tanh(a + b*x)), x)
def replacement5876(c, m, n, x, d, a, b):
        # rubi.append(5876)
        return Dist((n + S(-2))/(n + S(-1)), Int((c + d*x)**m*(S(1)/cosh(a + b*x))**(n + S(-2)), x), x) - Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*(n + S(-2))*(n + S(-1))), Int((c + d*x)**(m + S(-2))*(S(1)/cosh(a + b*x))**(n + S(-2)), x), x) + Simp((c + d*x)**m*(S(1)/cosh(a + b*x))**(n + S(-2))*tanh(a + b*x)/(b*(n + S(-1))), x) + Simp(d*m*(c + d*x)**(m + S(-1))*(S(1)/cosh(a + b*x))**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5877(c, m, n, x, d, a, b):
        # rubi.append(5877)
        return -Dist((n + S(-2))/(n + S(-1)), Int((c + d*x)**m*(S(1)/sinh(a + b*x))**(n + S(-2)), x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*(n + S(-2))*(n + S(-1))), Int((c + d*x)**(m + S(-2))*(S(1)/sinh(a + b*x))**(n + S(-2)), x), x) - Simp((c + d*x)**m*(S(1)/sinh(a + b*x))**(n + S(-2))/(b*(n + S(-1))*tanh(a + b*x)), x) - Simp(d*m*(c + d*x)**(m + S(-1))*(S(1)/sinh(a + b*x))**(n + S(-2))/(b**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5878(c, n, x, d, a, b):
        # rubi.append(5878)
        return Dist((n + S(1))/n, Int((c + d*x)*(S(1)/cosh(a + b*x))**(n + S(2)), x), x) - Simp(d*(S(1)/cosh(a + b*x))**n/(b**S(2)*n**S(2)), x) - Simp((c + d*x)*(S(1)/cosh(a + b*x))**(n + S(1))*sinh(a + b*x)/(b*n), x)
def replacement5879(c, n, x, d, a, b):
        # rubi.append(5879)
        return -Dist((n + S(1))/n, Int((c + d*x)*(S(1)/sinh(a + b*x))**(n + S(2)), x), x) - Simp(d*(S(1)/sinh(a + b*x))**n/(b**S(2)*n**S(2)), x) - Simp((c + d*x)*(S(1)/sinh(a + b*x))**(n + S(1))*cosh(a + b*x)/(b*n), x)
def replacement5880(c, m, n, x, d, a, b):
        # rubi.append(5880)
        return Dist((n + S(1))/n, Int((c + d*x)**m*(S(1)/cosh(a + b*x))**(n + S(2)), x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*n**S(2)), Int((c + d*x)**(m + S(-2))*(S(1)/cosh(a + b*x))**n, x), x) - Simp((c + d*x)**m*(S(1)/cosh(a + b*x))**(n + S(1))*sinh(a + b*x)/(b*n), x) - Simp(d*m*(c + d*x)**(m + S(-1))*(S(1)/cosh(a + b*x))**n/(b**S(2)*n**S(2)), x)
def replacement5881(c, m, n, x, d, a, b):
        # rubi.append(5881)
        return -Dist((n + S(1))/n, Int((c + d*x)**m*(S(1)/sinh(a + b*x))**(n + S(2)), x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*n**S(2)), Int((c + d*x)**(m + S(-2))*(S(1)/sinh(a + b*x))**n, x), x) - Simp((c + d*x)**m*(S(1)/sinh(a + b*x))**(n + S(1))*cosh(a + b*x)/(b*n), x) - Simp(d*m*(c + d*x)**(m + S(-1))*(S(1)/sinh(a + b*x))**n/(b**S(2)*n**S(2)), x)
def replacement5882(c, m, n, x, d, a, b):
        # rubi.append(5882)
        return Dist((S(1)/cosh(a + b*x))**n*cosh(a + b*x)**n, Int((c + d*x)**m*cosh(a + b*x)**(-n), x), x)
def replacement5883(c, m, n, x, d, a, b):
        # rubi.append(5883)
        return Dist((S(1)/sinh(a + b*x))**n*sinh(a + b*x)**n, Int((c + d*x)**m*sinh(a + b*x)**(-n), x), x)
def replacement5884(c, m, n, x, d, f, a, e, b):
        # rubi.append(5884)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b/cosh(e + f*x))**n, x), x)
def replacement5885(c, m, n, x, d, f, a, e, b):
        # rubi.append(5885)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b/sinh(e + f*x))**n, x), x)
def replacement5886(c, m, n, x, d, f, a, e, b):
        # rubi.append(5886)
        return Int(ExpandIntegrand((c + d*x)**m, (a*cosh(e + f*x) + b)**n*cosh(e + f*x)**(-n), x), x)
def replacement5887(c, m, n, x, d, f, a, e, b):
        # rubi.append(5887)
        return Int(ExpandIntegrand((c + d*x)**m, (a*sinh(e + f*x) + b)**n*sinh(e + f*x)**(-n), x), x)
def replacement5888(m, n, x, v, u):
        # rubi.append(5888)
        return Int((S(1)/cosh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5889(m, n, x, v, u):
        # rubi.append(5889)
        return Int((S(1)/sinh(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement5890(c, m, n, x, d, a, b):
        # rubi.append(5890)
        return Int((c + d*x)**m*(S(1)/cosh(a + b*x))**n, x)
def replacement5891(c, m, n, x, d, a, b):
        # rubi.append(5891)
        return Int((c + d*x)**m*(S(1)/sinh(a + b*x))**n, x)
def replacement5892(c, n, x, d, a, p, b):
        # rubi.append(5892)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/cosh(c + d*x))**p, x), x, x**n), x)
def replacement5893(c, n, x, d, a, p, b):
        # rubi.append(5893)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/sinh(c + d*x))**p, x), x, x**n), x)
def replacement5894(c, n, x, d, a, p, b):
        # rubi.append(5894)
        return Int((a + b/cosh(c + d*x**n))**p, x)
def replacement5895(c, n, x, d, a, p, b):
        # rubi.append(5895)
        return Int((a + b/sinh(c + d*x**n))**p, x)
def replacement5896(c, n, x, d, u, a, p, b):
        # rubi.append(5896)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/cosh(c + d*x**n))**p, x), x, u), x)
def replacement5897(c, n, x, d, u, a, p, b):
        # rubi.append(5897)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/sinh(c + d*x**n))**p, x), x, u), x)
def replacement5898(x, u, a, p, b):
        # rubi.append(5898)
        return Int((a + b/cosh(ExpandToSum(u, x)))**p, x)
def replacement5899(x, u, a, p, b):
        # rubi.append(5899)
        return Int((a + b/sinh(ExpandToSum(u, x)))**p, x)
def replacement5900(c, m, n, x, d, a, p, b):
        # rubi.append(5900)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/cosh(c + d*x))**p, x), x, x**n), x)
def replacement5901(c, m, n, x, d, a, p, b):
        # rubi.append(5901)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/sinh(c + d*x))**p, x), x, x**n), x)
def replacement5902(c, m, n, x, d, a, p, b):
        # rubi.append(5902)
        return Int(x**m*(a + b/cosh(c + d*x**n))**p, x)
def replacement5903(c, m, n, x, d, a, p, b):
        # rubi.append(5903)
        return Int(x**m*(a + b/sinh(c + d*x**n))**p, x)
def replacement5904(c, m, n, x, d, a, p, e, b):
        # rubi.append(5904)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/cosh(c + d*x**n))**p, x), x)
def replacement5905(c, m, n, x, d, a, p, e, b):
        # rubi.append(5905)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/sinh(c + d*x**n))**p, x), x)
def replacement5906(m, x, u, a, p, e, b):
        # rubi.append(5906)
        return Int((e*x)**m*(a + b/cosh(ExpandToSum(u, x)))**p, x)
def replacement5907(m, x, u, a, p, e, b):
        # rubi.append(5907)
        return Int((e*x)**m*(a + b/sinh(ExpandToSum(u, x)))**p, x)
def replacement5908(m, n, x, a, p, b):
        # rubi.append(5908)
        return Dist((m - n + S(1))/(b*n*(p + S(-1))), Int(x**(m - n)*(S(1)/cosh(a + b*x**n))**(p + S(-1)), x), x) - Simp(x**(m - n + S(1))*(S(1)/cosh(a + b*x**n))**(p + S(-1))/(b*n*(p + S(-1))), x)
def replacement5909(m, n, x, a, p, b):
        # rubi.append(5909)
        return Dist((m - n + S(1))/(b*n*(p + S(-1))), Int(x**(m - n)*(S(1)/sinh(a + b*x**n))**(p + S(-1)), x), x) - Simp(x**(m - n + S(1))*(S(1)/sinh(a + b*x**n))**(p + S(-1))/(b*n*(p + S(-1))), x)
def replacement5910(c, m, n, x, d, a, b):
        # rubi.append(5910)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*sinh(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*sinh(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement5911(c, m, n, x, d, a, b):
        # rubi.append(5911)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*cosh(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*cosh(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement5912(c, m, n, x, d, a, p, b):
        # rubi.append(5912)
        return Int(ExpandTrigReduce((c + d*x)**m, sinh(a + b*x)**n*cosh(a + b*x)**p, x), x)
def replacement5913(c, m, n, x, d, a, p, b):
        # rubi.append(5913)
        return Int((c + d*x)**m*sinh(a + b*x)**n*tanh(a + b*x)**(p + S(-2)), x) - Int((c + d*x)**m*sinh(a + b*x)**(n + S(-2))*tanh(a + b*x)**p, x)
def replacement5914(c, m, n, x, d, a, p, b):
        # rubi.append(5914)
        return Int((c + d*x)**m*(S(1)/tanh(a + b*x))**p*cosh(a + b*x)**(n + S(-2)), x) + Int((c + d*x)**m*(S(1)/tanh(a + b*x))**(p + S(-2))*cosh(a + b*x)**n, x)
def replacement5915(c, m, n, x, d, a, p, b):
        # rubi.append(5915)
        return Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*(S(1)/cosh(a + b*x))**n, x), x) - Simp((c + d*x)**m*(S(1)/cosh(a + b*x))**n/(b*n), x)
def replacement5916(c, m, n, x, d, a, p, b):
        # rubi.append(5916)
        return Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*(S(1)/sinh(a + b*x))**n, x), x) - Simp((c + d*x)**m*(S(1)/sinh(a + b*x))**n/(b*n), x)
def replacement5917(c, m, n, x, d, a, b):
        # rubi.append(5917)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*tanh(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*tanh(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement5918(c, m, n, x, d, a, b):
        # rubi.append(5918)
        return Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*(S(1)/tanh(a + b*x))**(n + S(1)), x), x) - Simp((c + d*x)**m*(S(1)/tanh(a + b*x))**(n + S(1))/(b*(n + S(1))), x)
def replacement5919(c, m, x, d, a, p, b):
        # rubi.append(5919)
        return -Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))/cosh(a + b*x)**S(3), x) + Int((c + d*x)**m*tanh(a + b*x)**(p + S(-2))/cosh(a + b*x), x)
def replacement5920(c, m, n, x, d, a, p, b):
        # rubi.append(5920)
        return Int((c + d*x)**m*(S(1)/cosh(a + b*x))**n*tanh(a + b*x)**(p + S(-2)), x) - Int((c + d*x)**m*(S(1)/cosh(a + b*x))**(n + S(2))*tanh(a + b*x)**(p + S(-2)), x)
def replacement5921(c, m, x, d, a, p, b):
        # rubi.append(5921)
        return Int((c + d*x)**m*(S(1)/tanh(a + b*x))**(p + S(-2))/sinh(a + b*x)**S(3), x) + Int((c + d*x)**m*(S(1)/tanh(a + b*x))**(p + S(-2))/sinh(a + b*x), x)
def replacement5922(c, m, n, x, d, a, p, b):
        # rubi.append(5922)
        return Int((c + d*x)**m*(S(1)/sinh(a + b*x))**n*(S(1)/tanh(a + b*x))**(p + S(-2)), x) + Int((c + d*x)**m*(S(1)/sinh(a + b*x))**(n + S(2))*(S(1)/tanh(a + b*x))**(p + S(-2)), x)

def With5923(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/cosh(a + b*x))**n*tanh(a + b*x)**p, x)
        # rubi.append(5923)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)

def With5924(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/sinh(a + b*x))**n*(S(1)/tanh(a + b*x))**p, x)
        # rubi.append(5924)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
def replacement5925(c, m, n, x, d, a, b):
        # rubi.append(5925)
        return Dist(S(2)**n, Int((c + d*x)**m*(S(1)/sinh(S(2)*a + S(2)*b*x))**n, x), x)

def With5926(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/sinh(a + b*x))**n*(S(1)/cosh(a + b*x))**p, x)
        # rubi.append(5926)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
def replacement5927(G, m, n, F, x, v, u, w, p):
        # rubi.append(5927)
        return Int(ExpandToSum(u, x)**m*F(ExpandToSum(v, x))**n*G(ExpandToSum(v, x))**p, x)
def replacement5928(c, m, n, x, d, f, a, e, b):
        # rubi.append(5928)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*sinh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*sinh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5929(c, m, n, x, d, f, a, e, b):
        # rubi.append(5929)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*cosh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*cosh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5930(c, m, n, x, d, f, a, e, b):
        # rubi.append(5930)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*tanh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*tanh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5931(c, m, n, x, d, f, a, e, b):
        # rubi.append(5931)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b/tanh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b/tanh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5932(c, m, n, x, d, f, a, e, b):
        # rubi.append(5932)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b/cosh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b/cosh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5933(c, m, n, x, d, f, a, e, b):
        # rubi.append(5933)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b/sinh(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b/sinh(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement5934(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(5934)
        return Int(ExpandTrigReduce((e + f*x)**m, sinh(a + b*x)**p*sinh(c + d*x)**q, x), x)
def replacement5935(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(5935)
        return Int(ExpandTrigReduce((e + f*x)**m, cosh(a + b*x)**p*cosh(c + d*x)**q, x), x)
def replacement5936(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(5936)
        return Int(ExpandTrigReduce((e + f*x)**m, sinh(a + b*x)**p*cosh(c + d*x)**q, x), x)
def replacement5937(c, G, m, q, F, x, d, f, a, p, e, b):
        # rubi.append(5937)
        return Int(ExpandTrigExpand((e + f*x)**m*G(c + d*x)**q, F, c + d*x, p, b/d, x), x)
def replacement5938(c, F, x, d, a, e, b):
        # rubi.append(5938)
        return Simp(F**(c*(a + b*x))*e*cosh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
def replacement5939(c, F, x, d, a, e, b):
        # rubi.append(5939)
        return Simp(F**(c*(a + b*x))*e*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
def replacement5940(c, n, F, x, d, a, e, b):
        # rubi.append(5940)
        return -Dist(e**S(2)*n*(n + S(-1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(-2)), x), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**n/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) + Simp(F**(c*(a + b*x))*e*n*sinh(d + e*x)**(n + S(-1))*cosh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement5941(c, n, F, x, d, a, e, b):
        # rubi.append(5941)
        return Dist(e**S(2)*n*(n + S(-1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*cosh(d + e*x)**(n + S(-2)), x), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)**n/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) + Simp(F**(c*(a + b*x))*e*n*sinh(d + e*x)*cosh(d + e*x)**(n + S(-1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement5942(c, n, F, x, d, a, e, b):
        # rubi.append(5942)
        return Simp(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(1))*cosh(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5943(c, n, F, x, d, a, e, b):
        # rubi.append(5943)
        return -Simp(F**(c*(a + b*x))*sinh(d + e*x)*cosh(d + e*x)**(n + S(1))/(e*(n + S(1))), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5944(c, n, F, x, d, a, e, b):
        # rubi.append(5944)
        return -Dist((-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*sinh(d + e*x)**(n + S(1))*cosh(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5945(c, n, F, x, d, a, e, b):
        # rubi.append(5945)
        return Dist((-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*cosh(d + e*x)**(n + S(2)), x), x) - Simp(F**(c*(a + b*x))*sinh(d + e*x)*cosh(d + e*x)**(n + S(1))/(e*(n + S(1))), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement5946(c, n, F, x, d, a, e, b):
        # rubi.append(5946)
        return Dist((exp(S(2)*d + S(2)*e*x) + S(-1))**(-n)*exp(n*(d + e*x))*sinh(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**n*exp(-n*(d + e*x)), x), x)
def replacement5947(c, n, F, x, d, a, e, b):
        # rubi.append(5947)
        return Dist((exp(S(2)*d + S(2)*e*x) + S(1))**(-n)*exp(n*(d + e*x))*cosh(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(1))**n*exp(-n*(d + e*x)), x), x)
def replacement5948(c, n, F, x, d, a, e, b):
        # rubi.append(5948)
        return Int(ExpandIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**n*(exp(S(2)*d + S(2)*e*x) + S(1))**(-n), x), x)
def replacement5949(c, n, F, x, d, a, e, b):
        # rubi.append(5949)
        return Int(ExpandIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(-1))**(-n)*(exp(S(2)*d + S(2)*e*x) + S(1))**n, x), x)
def replacement5950(c, n, F, x, d, a, e, b):
        # rubi.append(5950)
        return Dist(e**S(2)*n*(n + S(1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*(S(1)/cosh(d + e*x))**(n + S(2)), x), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/cosh(d + e*x))**n*log(F)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*(S(1)/cosh(d + e*x))**(n + S(1))*sinh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement5951(c, n, F, x, d, a, e, b):
        # rubi.append(5951)
        return -Dist(e**S(2)*n*(n + S(1))/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*(S(1)/sinh(d + e*x))**(n + S(2)), x), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/sinh(d + e*x))**n*log(F)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*(S(1)/sinh(d + e*x))**(n + S(1))*cosh(d + e*x)/(-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement5952(c, n, F, x, d, a, e, b):
        # rubi.append(5952)
        return Simp(F**(c*(a + b*x))*(S(1)/cosh(d + e*x))**(n + S(-1))*sinh(d + e*x)/(e*(n + S(-1))), x) + Simp(F**(c*(a + b*x))*b*c*(S(1)/cosh(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5953(c, n, F, x, d, a, e, b):
        # rubi.append(5953)
        return -Simp(F**(c*(a + b*x))*(S(1)/sinh(d + e*x))**(n + S(-1))*cosh(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/sinh(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5954(c, n, F, x, d, a, e, b):
        # rubi.append(5954)
        return Dist((-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*(S(1)/cosh(d + e*x))**(n + S(-2)), x), x) + Simp(F**(c*(a + b*x))*(S(1)/cosh(d + e*x))**(n + S(-1))*sinh(d + e*x)/(e*(n + S(-1))), x) + Simp(F**(c*(a + b*x))*b*c*(S(1)/cosh(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5955(c, n, F, x, d, a, e, b):
        # rubi.append(5955)
        return -Dist((-b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*(S(1)/sinh(d + e*x))**(n + S(-2)), x), x) - Simp(F**(c*(a + b*x))*(S(1)/sinh(d + e*x))**(n + S(-1))*cosh(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/sinh(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement5956(c, n, F, x, d, a, e, b):
        # rubi.append(5956)
        return Simp(S(2)**n*F**(c*(a + b*x))*Hypergeometric2F1(n, b*c*log(F)/(S(2)*e) + n/S(2), b*c*log(F)/(S(2)*e) + n/S(2) + S(1), -exp(S(2)*d + S(2)*e*x))*exp(n*(d + e*x))/(b*c*log(F) + e*n), x)
def replacement5957(c, n, F, x, d, a, e, b):
        # rubi.append(5957)
        return Simp((S(-2))**n*F**(c*(a + b*x))*Hypergeometric2F1(n, b*c*log(F)/(S(2)*e) + n/S(2), b*c*log(F)/(S(2)*e) + n/S(2) + S(1), exp(S(2)*d + S(2)*e*x))*exp(n*(d + e*x))/(b*c*log(F) + e*n), x)
def replacement5958(c, n, F, x, d, a, e, b):
        # rubi.append(5958)
        return Dist((exp(S(2)*d + S(2)*e*x) + S(1))**n*(S(1)/cosh(d + e*x))**n*exp(-n*(d + e*x)), Int(SimplifyIntegrand(F**(c*(a + b*x))*(exp(S(2)*d + S(2)*e*x) + S(1))**(-n)*exp(n*(d + e*x)), x), x), x)
def replacement5959(c, n, F, x, d, a, e, b):
        # rubi.append(5959)
        return Dist((S(1) - exp(-S(2)*d - S(2)*e*x))**n*(S(1)/sinh(d + e*x))**n*exp(n*(d + e*x)), Int(SimplifyIntegrand(F**(c*(a + b*x))*(S(1) - exp(-S(2)*d - S(2)*e*x))**(-n)*exp(-n*(d + e*x)), x), x), x)
def replacement5960(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(5960)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*cosh(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement5961(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(5961)
        return Dist(S(2)**n*g**n, Int(F**(c*(a + b*x))*cosh(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement5962(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(5962)
        return Dist(S(2)**n*g**n, Int(F**(c*(a + b*x))*sinh(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement5963(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(5963)
        return Dist(g**n, Int(F**(c*(a + b*x))*tanh(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**m, x), x)
def replacement5964(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(5964)
        return Dist(g**n, Int(F**(c*(a + b*x))*tanh(d/S(2) + e*x/S(2))**m, x), x)
def replacement5965(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(5965)
        return Dist(g**n, Int(F**(c*(a + b*x))*(S(1)/tanh(d/S(2) + e*x/S(2)))**m, x), x)
def replacement5966(c, F, x, d, h, f, i, g, a, e, b):
        # rubi.append(5966)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*cosh(d + e*x)/(f + g*sinh(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*cosh(d + e*x))/(f + g*sinh(d + e*x)), x)
def replacement5967(c, F, x, d, h, f, g, i, a, e, b):
        # rubi.append(5967)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*sinh(d + e*x)/(f + g*cosh(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*sinh(d + e*x))/(f + g*cosh(d + e*x)), x)
def replacement5968(c, G, n, F, x, v, u):
        # rubi.append(5968)
        return Int(F**(c*ExpandToSum(u, x))*G(ExpandToSum(v, x))**n, x)

def With5969(c, m, n, F, x, d, a, e, b):
        u = IntHide(F**(c*(a + b*x))*sinh(d + e*x)**n, x)
        # rubi.append(5969)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Simp(u*x**m, x)

def With5970(c, m, n, F, x, d, a, e, b):
        u = IntHide(F**(c*(a + b*x))*cosh(d + e*x)**n, x)
        # rubi.append(5970)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Simp(u*x**m, x)
def replacement5971(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(5971)
        return Int(ExpandTrigReduce(F**(c*(a + b*x)), sinh(d + e*x)**m*cosh(f + g*x)**n, x), x)
def replacement5972(c, m, n, F, x, d, f, g, a, p, e, b):
        # rubi.append(5972)
        return Int(ExpandTrigReduce(F**(c*(a + b*x))*x**p, sinh(d + e*x)**m*cosh(f + g*x)**n, x), x)
def replacement5973(c, G, m, n, F, x, d, H, a, e, b):
        # rubi.append(5973)
        return Int(ExpandTrigToExp(F**(c*(a + b*x)), G(d + e*x)**m*H(d + e*x)**n, x), x)
def replacement5974(n, F, x, v, u):
        # rubi.append(5974)
        return Int(ExpandTrigToExp(F**u, sinh(v)**n, x), x)
def replacement5975(n, F, x, v, u):
        # rubi.append(5975)
        return Int(ExpandTrigToExp(F**u, cosh(v)**n, x), x)
def replacement5976(m, n, F, x, v, u):
        # rubi.append(5976)
        return Int(ExpandTrigToExp(F**u, sinh(v)**m*cosh(v)**n, x), x)
def replacement5977(c, n, x, p, b):
        # rubi.append(5977)
        return Int(((c*x**n)**b/S(2) - (c*x**n)**(-b)/S(2))**p, x)
def replacement5978(c, n, x, p, b):
        # rubi.append(5978)
        return Int(((c*x**n)**b/S(2) + (c*x**n)**(-b)/S(2))**p, x)
def replacement5979(c, n, x, a, p, b):
        # rubi.append(5979)
        return -Simp(x*(p + S(2))*sinh(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) + Simp(x*sinh(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tanh(a + b*log(c*x**n))), x)
def replacement5980(c, n, x, a, p, b):
        # rubi.append(5980)
        return Simp(x*(p + S(2))*cosh(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) - Simp(x*cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement5981(c, n, x, a, b):
        # rubi.append(5981)
        return Dist(x*sqrt(sinh(a + b*log(c*x**n)))/sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(-1)), Int(sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(-1))/x, x), x)
def replacement5982(c, n, x, a, b):
        # rubi.append(5982)
        return Dist(x*sqrt(cosh(a + b*log(c*x**n)))/sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(1)), Int(sqrt((c*x**n)**(S(4)/n)*exp(S(2)*a) + S(1))/x, x), x)
def replacement5983(c, n, x, a, p, b):
        # rubi.append(5983)
        return Int(ExpandIntegrand(((c*x**n)**(S(1)/(n*p))*exp(a*b*n*p)/(S(2)*b*n*p) - (c*x**n)**(-S(1)/(n*p))*exp(-a*b*n*p)/(S(2)*b*n*p))**p, x), x)
def replacement5984(c, n, x, a, p, b):
        # rubi.append(5984)
        return Int(ExpandIntegrand(((c*x**n)**(S(1)/(n*p))*exp(a*b*n*p)/S(2) + (c*x**n)**(-S(1)/(n*p))*exp(-a*b*n*p)/S(2))**p, x), x)
def replacement5985(c, n, x, a, b):
        # rubi.append(5985)
        return -Simp(x*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)), x) + Simp(b*n*x*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)), x)
def replacement5986(c, n, x, a, b):
        # rubi.append(5986)
        return -Simp(x*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)), x) + Simp(b*n*x*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(-1)), x)
def replacement5987(c, n, x, a, p, b):
        # rubi.append(5987)
        return -Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), Int(sinh(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x*sinh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x) + Simp(b*n*p*x*sinh(a + b*log(c*x**n))**(p + S(-1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x)
def replacement5988(c, n, x, a, p, b):
        # rubi.append(5988)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), Int(cosh(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x*cosh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x) + Simp(b*n*p*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x)
def replacement5989(c, n, x, a, p, b):
        # rubi.append(5989)
        return -Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(sinh(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x*sinh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) + Simp(x*sinh(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tanh(a + b*log(c*x**n))), x)
def replacement5990(c, n, x, a, p, b):
        # rubi.append(5990)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(-1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(cosh(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x*cosh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) - Simp(x*cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement5991(c, n, x, a, p, b):
        # rubi.append(5991)
        return Simp(x*(S(2) - S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) - (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, -(b*n*p + S(1))/(S(2)*b*n), S(1) - (b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + S(1)), x)
def replacement5992(c, n, x, a, p, b):
        # rubi.append(5992)
        return Simp(x*(S(2) + S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) + (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, -(b*n*p + S(1))/(S(2)*b*n), S(1) - (b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + S(1)), x)
def replacement5993(c, m, n, x, a, p, b):
        # rubi.append(5993)
        return -Simp(x**(m + S(1))*(p + S(2))*sinh(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) + Simp(x**(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tanh(a + b*log(c*x**n))), x)
def replacement5994(c, m, n, x, a, p, b):
        # rubi.append(5994)
        return Simp(x**(m + S(1))*(p + S(2))*cosh(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) - Simp(x**(m + S(1))*cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement5995(c, m, n, x, a, p, b):
        # rubi.append(5995)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*((c*x**n)**((m + S(1))/(n*p))*(m + S(1))*exp(a*b*n*p/(m + S(1)))/(b*n*p) - (c*x**n)**(-(m + S(1))/(n*p))*(m + S(1))*exp(-a*b*n*p/(m + S(1)))/(b*n*p))**p, x), x), x)
def replacement5996(c, m, n, x, a, p, b):
        # rubi.append(5996)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*((c*x**n)**((m + S(1))/(n*p))*exp(a*b*n*p/(m + S(1))) + (c*x**n)**(-(m + S(1))/(n*p))*exp(-a*b*n*p/(m + S(1))))**p, x), x), x)
def replacement5997(c, m, n, x, a, b):
        # rubi.append(5997)
        return -Simp(x**(m + S(1))*(m + S(1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)), x) + Simp(b*n*x**(m + S(1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)), x)
def replacement5998(c, m, n, x, a, b):
        # rubi.append(5998)
        return -Simp(x**(m + S(1))*(m + S(1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)), x) + Simp(b*n*x**(m + S(1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2) - (m + S(1))**S(2)), x)
def replacement5999(c, m, n, x, a, p, b):
        # rubi.append(5999)
        return -Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), Int(x**m*sinh(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x**(m + S(1))*(m + S(1))*sinh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(-1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x)
def replacement6000(c, m, n, x, a, p, b):
        # rubi.append(6000)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), Int(x**m*cosh(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x**(m + S(1))*(m + S(1))*cosh(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x)
def replacement6001(c, m, n, x, a, p, b):
        # rubi.append(6001)
        return -Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*sinh(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x**(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tanh(a + b*log(c*x**n))), x) - Simp(x**(m + S(1))*(m + S(1))*sinh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement6002(c, m, n, x, a, p, b):
        # rubi.append(6002)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) - (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*cosh(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x**(m + S(1))*cosh(a + b*log(c*x**n))**(p + S(2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(1))), x) + Simp(x**(m + S(1))*(m + S(1))*cosh(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement6003(c, m, n, x, a, p, b):
        # rubi.append(6003)
        return Simp(x**(m + S(1))*(S(2) - S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) - (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, -(b*n*p + m + S(1))/(S(2)*b*n), S(1) - (b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + m + S(1)), x)
def replacement6004(c, m, n, x, a, p, b):
        # rubi.append(6004)
        return Simp(x**(m + S(1))*(S(2) + S(2)*(c*x**n)**(-S(2)*b)*exp(-S(2)*a))**(-p)*((c*x**n)**b*exp(a) + (c*x**n)**(-b)*exp(-a))**p*Hypergeometric2F1(-p, -(b*n*p + m + S(1))/(S(2)*b*n), S(1) - (b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(-S(2)*b)*exp(-S(2)*a))/(b*n*p + m + S(1)), x)
def replacement6005(c, n, x, p, b):
        # rubi.append(6005)
        return Dist(S(2)**p, Int(((c*x**n)**b/((c*x**n)**(S(2)*b) + S(1)))**p, x), x)
def replacement6006(c, n, x, p, b):
        # rubi.append(6006)
        return Dist(S(2)**p, Int(((c*x**n)**b/((c*x**n)**(S(2)*b) + S(-1)))**p, x), x)
def replacement6007(c, n, x, a, b):
        # rubi.append(6007)
        return Dist(S(2)*exp(-a*b*n), Int((c*x**n)**(S(1)/n)/((c*x**n)**(S(2)/n) + exp(-S(2)*a*b*n)), x), x)
def replacement6008(c, n, x, a, b):
        # rubi.append(6008)
        return Dist(-S(2)*b*n*exp(-a*b*n), Int((c*x**n)**(S(1)/n)/(-(c*x**n)**(S(2)/n) + exp(-S(2)*a*b*n)), x), x)
def replacement6009(c, n, x, a, p, b):
        # rubi.append(6009)
        return Simp(x*(p + S(-2))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))/(p + S(-1)), x) + Simp(x*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement6010(c, n, x, a, p, b):
        # rubi.append(6010)
        return -Simp(x*(p + S(-2))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(p + S(-1)), x) - Simp(x*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tanh(a + b*log(c*x**n))), x)
def replacement6011(c, n, x, a, p, b):
        # rubi.append(6011)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int((S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2)), x), x) + Simp(x*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) + Simp(x*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement6012(c, n, x, a, p, b):
        # rubi.append(6012)
        return -Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(-1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int((S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2)), x), x) - Simp(x*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) - Simp(x*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tanh(a + b*log(c*x**n))), x)
def replacement6013(c, n, x, a, p, b):
        # rubi.append(6013)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), Int((S(1)/cosh(a + b*log(c*x**n)))**(p + S(2)), x), x) - Simp(x*(S(1)/cosh(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x) - Simp(b*n*p*x*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x)
def replacement6014(c, n, x, a, p, b):
        # rubi.append(6014)
        return -Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), Int((S(1)/sinh(a + b*log(c*x**n)))**(p + S(2)), x), x) - Simp(x*(S(1)/sinh(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x) - Simp(b*n*p*x*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(-1)), x)
def replacement6015(c, n, x, a, p, b):
        # rubi.append(6015)
        return Simp(S(2)**p*x*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1)))**p*((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + S(1))/(S(2)*b*n), S(1) + (b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + S(1)), x)
def replacement6016(c, n, x, a, p, b):
        # rubi.append(6016)
        return Simp(x*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(-1)))**p*(-S(2)*(c*x**n)**(S(2)*b)*exp(S(2)*a) + S(2))**p*Hypergeometric2F1(p, (b*n*p + S(1))/(S(2)*b*n), S(1) + (b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + S(1)), x)
def replacement6017(c, m, n, x, p, b):
        # rubi.append(6017)
        return Dist(S(2)**p, Int(x**m*((c*x**n)**b/((c*x**n)**(S(2)*b) + S(1)))**p, x), x)
def replacement6018(c, m, n, x, p, b):
        # rubi.append(6018)
        return Dist(S(2)**p, Int(x**m*((c*x**n)**b/((c*x**n)**(S(2)*b) + S(-1)))**p, x), x)
def replacement6019(c, m, n, x, a, b):
        # rubi.append(6019)
        return Dist(S(2)*exp(-a*b*n/(m + S(1))), Int(x**m*(c*x**n)**((m + S(1))/n)/((c*x**n)**(S(2)*(m + S(1))/n) + exp(-S(2)*a*b*n/(m + S(1)))), x), x)
def replacement6020(c, m, n, x, a, b):
        # rubi.append(6020)
        return Dist(-S(2)*b*n*exp(-a*b*n/(m + S(1)))/(m + S(1)), Int(x**m*(c*x**n)**((m + S(1))/n)/(-(c*x**n)**(S(2)*(m + S(1))/n) + exp(-S(2)*a*b*n/(m + S(1)))), x), x)
def replacement6021(c, m, n, x, a, p, b):
        # rubi.append(6021)
        return Simp(x**(m + S(1))*(p + S(-2))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) + Simp(x**(m + S(1))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement6022(c, m, n, x, a, p, b):
        # rubi.append(6022)
        return -Simp(x**(m + S(1))*(p + S(-2))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) - Simp(x**(m + S(1))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tanh(a + b*log(c*x**n))), x)
def replacement6023(c, m, n, x, a, p, b):
        # rubi.append(6023)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))*tanh(a + b*log(c*x**n))/(b*n*(p + S(-1))), x) + Simp(x**(m + S(1))*(m + S(1))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
def replacement6024(c, m, n, x, a, p, b):
        # rubi.append(6024)
        return -Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) - (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2)), x), x) - Simp(x**(m + S(1))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tanh(a + b*log(c*x**n))), x) - Simp(x**(m + S(1))*(m + S(1))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
def replacement6025(c, m, n, x, a, p, b):
        # rubi.append(6025)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), Int(x**m*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(2)), x), x) - Simp(x**(m + S(1))*(m + S(1))*(S(1)/cosh(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*(S(1)/cosh(a + b*log(c*x**n)))**(p + S(1))*sinh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x)
def replacement6026(c, m, n, x, a, p, b):
        # rubi.append(6026)
        return -Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), Int(x**m*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(2)), x), x) - Simp(x**(m + S(1))*(m + S(1))*(S(1)/sinh(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*(S(1)/sinh(a + b*log(c*x**n)))**(p + S(1))*cosh(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) - (m + S(1))**S(2)), x)
def replacement6027(c, m, n, x, a, p, b):
        # rubi.append(6027)
        return Simp(S(2)**p*x**(m + S(1))*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1)))**p*((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + m + S(1))/(S(2)*b*n), S(1) + (b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + m + S(1)), x)
def replacement6028(c, m, n, x, a, p, b):
        # rubi.append(6028)
        return Simp(S(2)**p*x**(m + S(1))*((c*x**n)**b*exp(a)/((c*x**n)**(S(2)*b)*exp(S(2)*a) + S(-1)))**p*(-(c*x**n)**(S(2)*b)*exp(S(2)*a) + S(1))**p*Hypergeometric2F1(p, (b*n*p + m + S(1))/(S(2)*b*n), S(1) + (b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*b)*exp(S(2)*a))/(b*n*p + m + S(1)), x)
def replacement6029(x, a, p, b):
        # rubi.append(6029)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*sinh(a*x*log(b*x)**p), x), x) + Simp(cosh(a*x*log(b*x)**p)/a, x)
def replacement6030(x, a, p, b):
        # rubi.append(6030)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*cosh(a*x*log(b*x)**p), x), x) + Simp(sinh(a*x*log(b*x)**p)/a, x)
def replacement6031(n, x, a, p, b):
        # rubi.append(6031)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x), x) + Dist((n + S(-1))/(a*n), Int(x**(-n)*cosh(a*x**n*log(b*x)**p), x), x) + Simp(x**(S(1) - n)*cosh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6032(n, x, a, p, b):
        # rubi.append(6032)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*cosh(a*x**n*log(b*x)**p), x), x) + Dist((n + S(-1))/(a*n), Int(x**(-n)*sinh(a*x**n*log(b*x)**p), x), x) + Simp(x**(S(1) - n)*sinh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6033(m, n, x, a, p, b):
        # rubi.append(6033)
        return -Dist(p/n, Int(x**(n + S(-1))*log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x), x) - Simp(cosh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6034(m, n, x, a, p, b):
        # rubi.append(6034)
        return -Dist(p/n, Int(x**(n + S(-1))*log(b*x)**(p + S(-1))*cosh(a*x**n*log(b*x)**p), x), x) + Simp(sinh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6035(m, n, x, a, p, b):
        # rubi.append(6035)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*sinh(a*x**n*log(b*x)**p), x), x) - Dist((m - n + S(1))/(a*n), Int(x**(m - n)*cosh(a*x**n*log(b*x)**p), x), x) + Simp(x**(m - n + S(1))*cosh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6036(m, n, x, a, p, b):
        # rubi.append(6036)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*cosh(a*x**n*log(b*x)**p), x), x) - Dist((m - n + S(1))/(a*n), Int(x**(m - n)*sinh(a*x**n*log(b*x)**p), x), x) + Simp(x**(m - n + S(1))*sinh(a*x**n*log(b*x)**p)/(a*n), x)
def replacement6037(c, n, x, d, a):
        # rubi.append(6037)
        return -Dist(S(1)/d, Subst(Int(sinh(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement6038(c, n, x, d, a):
        # rubi.append(6038)
        return -Dist(S(1)/d, Subst(Int(cosh(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement6039(c, n, x, d, a, e, b):
        # rubi.append(6039)
        return -Dist(S(1)/d, Subst(Int(sinh(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement6040(c, n, x, d, a, e, b):
        # rubi.append(6040)
        return -Dist(S(1)/d, Subst(Int(cosh(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)

def With6041(x, u, n):
        lst = QuotientOfLinearsParts(u, x)
        # rubi.append(6041)
        return Int(sinh((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)

def With6042(x, u, n):
        lst = QuotientOfLinearsParts(u, x)
        # rubi.append(6042)
        return Int(cosh((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
def replacement6043(q, x, v, u, w, p):
        # rubi.append(6043)
        return Int(u*sinh(v)**(p + q), x)
def replacement6044(q, x, v, u, w, p):
        # rubi.append(6044)
        return Int(u*cosh(v)**(p + q), x)
def replacement6045(q, x, v, w, p):
        # rubi.append(6045)
        return Int(ExpandTrigReduce(sinh(v)**p*sinh(w)**q, x), x)
def replacement6046(q, x, v, w, p):
        # rubi.append(6046)
        return Int(ExpandTrigReduce(cosh(v)**p*cosh(w)**q, x), x)
def replacement6047(m, q, x, v, w, p):
        # rubi.append(6047)
        return Int(ExpandTrigReduce(x**m, sinh(v)**p*sinh(w)**q, x), x)
def replacement6048(m, q, x, v, w, p):
        # rubi.append(6048)
        return Int(ExpandTrigReduce(x**m, cosh(v)**p*cosh(w)**q, x), x)
def replacement6049(x, v, u, w, p):
        # rubi.append(6049)
        return Dist(S(2)**(-p), Int(u*sinh(S(2)*v)**p, x), x)
def replacement6050(q, x, v, w, p):
        # rubi.append(6050)
        return Int(ExpandTrigReduce(sinh(v)**p*cosh(w)**q, x), x)
def replacement6051(m, q, x, v, w, p):
        # rubi.append(6051)
        return Int(ExpandTrigReduce(x**m, sinh(v)**p*cosh(w)**q, x), x)
def replacement6052(x, w, v, n):
        # rubi.append(6052)
        return -Dist(cosh(v - w), Int(tanh(w)**(n + S(-1))/cosh(w), x), x) + Int(cosh(v)*tanh(w)**(n + S(-1)), x)
def replacement6053(x, w, v, n):
        # rubi.append(6053)
        return Dist(cosh(v - w), Int((S(1)/tanh(w))**(n + S(-1))/sinh(w), x), x) + Int((S(1)/tanh(w))**(n + S(-1))*sinh(v), x)
def replacement6054(x, w, v, n):
        # rubi.append(6054)
        return Dist(sinh(v - w), Int((S(1)/tanh(w))**(n + S(-1))/sinh(w), x), x) + Int((S(1)/tanh(w))**(n + S(-1))*cosh(v), x)
def replacement6055(x, w, v, n):
        # rubi.append(6055)
        return -Dist(sinh(v - w), Int(tanh(w)**(n + S(-1))/cosh(w), x), x) + Int(sinh(v)*tanh(w)**(n + S(-1)), x)
def replacement6056(x, w, v, n):
        # rubi.append(6056)
        return Dist(sinh(v - w), Int((S(1)/cosh(w))**(n + S(-1)), x), x) + Dist(cosh(v - w), Int((S(1)/cosh(w))**(n + S(-1))*tanh(w), x), x)
def replacement6057(x, w, v, n):
        # rubi.append(6057)
        return Dist(sinh(v - w), Int((S(1)/sinh(w))**(n + S(-1)), x), x) + Dist(cosh(v - w), Int((S(1)/sinh(w))**(n + S(-1))/tanh(w), x), x)
def replacement6058(x, w, v, n):
        # rubi.append(6058)
        return Dist(sinh(v - w), Int((S(1)/sinh(w))**(n + S(-1))/tanh(w), x), x) + Dist(cosh(v - w), Int((S(1)/sinh(w))**(n + S(-1)), x), x)
def replacement6059(x, w, v, n):
        # rubi.append(6059)
        return Dist(sinh(v - w), Int((S(1)/cosh(w))**(n + S(-1))*tanh(w), x), x) + Dist(cosh(v - w), Int((S(1)/cosh(w))**(n + S(-1)), x), x)
def replacement6060(c, m, n, x, d, f, a, e, b):
        # rubi.append(6060)
        return Int((a + b*sinh(S(2)*c + S(2)*d*x)/S(2))**n*(e + f*x)**m, x)
def replacement6061(c, m, n, x, d, a, b):
        # rubi.append(6061)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a + b*cosh(S(2)*c + S(2)*d*x) - b)**n, x), x)
def replacement6062(c, m, n, x, d, a, b):
        # rubi.append(6062)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a + b*cosh(S(2)*c + S(2)*d*x) + b)**n, x), x)
def replacement6063(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(6063)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*sinh(a + b*x**n)**p, x), x, c + d*x), x)
def replacement6064(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(6064)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*cosh(a + b*x**n)**p, x), x, c + d*x), x)
def replacement6065(c, m, x, d, f, g, a, e, b):
        # rubi.append(6065)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*cosh(S(2)*d + S(2)*e*x)), x), x)
def replacement6066(c, m, x, d, f, g, e, b):
        # rubi.append(6066)
        return Dist(S(2), Int((f + g*x)**m/(b - c + (b + c)*cosh(S(2)*d + S(2)*e*x)), x), x)
def replacement6067(c, m, x, d, f, g, a, e, b):
        # rubi.append(6067)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*cosh(S(2)*d + S(2)*e*x)), x), x)
def replacement6068(c, m, x, d, f, g, e, b):
        # rubi.append(6068)
        return Dist(S(2), Int((f + g*x)**m/(b - c + (b + c)*cosh(S(2)*d + S(2)*e*x)), x), x)
def replacement6069(c, m, x, d, f, g, a, e, b):
        # rubi.append(6069)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b - c + (b + c)*cosh(S(2)*d + S(2)*e*x)), x), x)
def replacement6070(c, m, x, d, f, a, e, b):
        # rubi.append(6070)
        return Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) - Rt(a**S(2) + b**S(2), S(2))), x) + Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) + Rt(a**S(2) + b**S(2), S(2))), x) - Simp((e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement6071(c, m, x, d, f, a, e, b):
        # rubi.append(6071)
        return Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) - Rt(a**S(2) - b**S(2), S(2))), x) + Int((e + f*x)**m*exp(c + d*x)/(a + b*exp(c + d*x) + Rt(a**S(2) - b**S(2), S(2))), x) - Simp((e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement6072(c, m, n, x, d, f, a, e, b):
        # rubi.append(6072)
        return Dist(S(1)/a, Int((e + f*x)**m*cosh(c + d*x)**(n + S(-2)), x), x) + Dist(S(1)/b, Int((e + f*x)**m*sinh(c + d*x)*cosh(c + d*x)**(n + S(-2)), x), x)
def replacement6073(c, m, n, x, d, f, a, e, b):
        # rubi.append(6073)
        return Dist(S(1)/a, Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2)), x), x) + Dist(S(1)/b, Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2))*cosh(c + d*x), x), x)
def replacement6074(c, m, n, x, d, f, a, e, b):
        # rubi.append(6074)
        return Dist(S(1)/b, Int((e + f*x)**m*sinh(c + d*x)*cosh(c + d*x)**(n + S(-2)), x), x) - Dist(a/b**S(2), Int((e + f*x)**m*cosh(c + d*x)**(n + S(-2)), x), x) + Dist((a**S(2) + b**S(2))/b**S(2), Int((e + f*x)**m*cosh(c + d*x)**(n + S(-2))/(a + b*sinh(c + d*x)), x), x)
def replacement6075(c, m, n, x, d, f, a, e, b):
        # rubi.append(6075)
        return Dist(S(1)/b, Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2))*cosh(c + d*x), x), x) - Dist(a/b**S(2), Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2)), x), x) + Dist((a**S(2) - b**S(2))/b**S(2), Int((e + f*x)**m*sinh(c + d*x)**(n + S(-2))/(a + b*cosh(c + d*x)), x), x)
def replacement6076(c, x, d, B, A, f, a, e, b):
        # rubi.append(6076)
        return -Dist(B*f/(a*d), Int(cosh(c + d*x)/(a + b*sinh(c + d*x)), x), x) + Simp(B*(e + f*x)*cosh(c + d*x)/(a*d*(a + b*sinh(c + d*x))), x)
def replacement6077(c, x, d, B, A, f, a, e, b):
        # rubi.append(6077)
        return -Dist(B*f/(a*d), Int(sinh(c + d*x)/(a + b*cosh(c + d*x)), x), x) + Simp(B*(e + f*x)*sinh(c + d*x)/(a*d*(a + b*cosh(c + d*x))), x)
def replacement6078(m, n, x, v, a, b):
        # rubi.append(6078)
        return Int((a*cosh(v) + b*sinh(v))**n, x)
def replacement6079(m, n, x, v, a, b):
        # rubi.append(6079)
        return Int((a*sinh(v) + b*cosh(v))**n, x)
def replacement6080(c, m, n, x, d, u, a, b):
        # rubi.append(6080)
        return Int(ExpandTrigReduce(u, sinh(a + b*x)**m*sinh(c + d*x)**n, x), x)
def replacement6081(c, m, n, x, d, u, a, b):
        # rubi.append(6081)
        return Int(ExpandTrigReduce(u, cosh(a + b*x)**m*cosh(c + d*x)**n, x), x)
def replacement6082(c, x, d, a, b):
        # rubi.append(6082)
        return Dist(S(1)/sinh((-a*d + b*c)/b), Int(tanh(c + d*x), x), x) - Dist(S(1)/sinh((-a*d + b*c)/d), Int(tanh(a + b*x), x), x)
def replacement6083(c, x, d, a, b):
        # rubi.append(6083)
        return Dist(S(1)/sinh((-a*d + b*c)/b), Int(S(1)/tanh(a + b*x), x), x) - Dist(S(1)/sinh((-a*d + b*c)/d), Int(S(1)/tanh(c + d*x), x), x)
def replacement6084(c, x, d, a, b):
        # rubi.append(6084)
        return -Dist(b*cosh((-a*d + b*c)/d)/d, Int(S(1)/(cosh(a + b*x)*cosh(c + d*x)), x), x) + Simp(b*x/d, x)
def replacement6085(c, x, d, a, b):
        # rubi.append(6085)
        return Dist(cosh((-a*d + b*c)/d), Int(S(1)/(sinh(a + b*x)*sinh(c + d*x)), x), x) + Simp(b*x/d, x)
def replacement6086(n, x, v, u, a, b):
        # rubi.append(6086)
        return Int(u*(a*exp(a*v/b))**n, x)
