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

def miscellaneous_trig(rubi):
    from sympy.integrals.rubi.constraints import cons1648, cons21, cons2, cons3, cons8, cons29, cons19, cons4, cons36, cons37, cons38, cons1649, cons1650, cons1651, cons1652, cons1653, cons1654, cons1655, cons27, cons1656, cons1410, cons210, cons40, cons50, cons127, cons149, cons345, cons5, cons242, cons246, cons1335, cons139, cons1657, cons1290, cons168, cons321, cons1658, cons33, cons251, cons96, cons255, cons13, cons165, cons248, cons1280, cons1659, cons1660, cons1661, cons172, cons1662, cons1663, cons95, cons91, cons1664, cons164, cons90, cons1665, cons1666, cons87, cons130, cons1481, cons746, cons1484, cons1667, cons25, cons1668, cons1669, cons1670, cons1671, cons1249, cons1672, cons1673, cons1674, cons1675, cons557, cons1676, cons630, cons10, cons1677, cons1678, cons1679, cons68, cons1232, cons378, cons51, cons52, cons53, cons54, cons1680, cons1441, cons1681, cons1682, cons1683, cons1684, cons64, cons586, cons466, cons1685, cons170, cons1686, cons1687, cons1688, cons1689, cons1690, cons814, cons815, cons20, cons1691, cons1692, cons1693, cons1694, cons1101, cons1695, cons89, cons167, cons1696, cons1697, cons1397, cons1698, cons1444, cons1699, cons1504, cons965, cons1700, cons1646, cons1701, cons198, cons1702, cons1013, cons152, cons1553, cons1703, cons1704, cons211, cons226, cons1705, cons812, cons813, cons150, cons530, cons1706, cons1707, cons1708, cons1709, cons1710, cons56, cons1711, cons1712, cons148, cons1713, cons1507, cons1714, cons1715, cons1716, cons1717, cons1718, cons1719, cons1720, cons1721, cons1647, cons1722, cons1723, cons1724, cons1725, cons1726, cons340, cons55, cons629, cons73, cons1727, cons1728, cons1729, cons1730, cons1362, cons1480, cons465, cons1731, cons1732, cons1733, cons1734, cons1267, cons1269, cons1476, cons1483, cons1735

    pattern4688 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1648, cons21)
    rule4688 = ReplacementRule(pattern4688, replacement4688)
    pattern4689 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1648, cons21)
    rule4689 = ReplacementRule(pattern4689, replacement4689)
    pattern4690 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1648, cons21)
    rule4690 = ReplacementRule(pattern4690, replacement4690)
    pattern4691 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1648, cons21)
    rule4691 = ReplacementRule(pattern4691, replacement4691)
    pattern4692 = Pattern(Integral(u_*(WC('c', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1648)
    rule4692 = ReplacementRule(pattern4692, replacement4692)
    pattern4693 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1648)
    rule4693 = ReplacementRule(pattern4693, replacement4693)
    pattern4694 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1648)
    rule4694 = ReplacementRule(pattern4694, replacement4694)
    pattern4695 = Pattern(Integral(u_*(WC('c', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1648)
    rule4695 = ReplacementRule(pattern4695, replacement4695)
    pattern4696 = Pattern(Integral(u_*(WC('c', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1648)
    rule4696 = ReplacementRule(pattern4696, replacement4696)
    pattern4697 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1648)
    rule4697 = ReplacementRule(pattern4697, replacement4697)
    pattern4698 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1648)
    rule4698 = ReplacementRule(pattern4698, replacement4698)
    pattern4699 = Pattern(Integral(u_*(A_ + WC('B', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1648)
    rule4699 = ReplacementRule(pattern4699, replacement4699)
    pattern4700 = Pattern(Integral(u_*(A_ + WC('B', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1648)
    rule4700 = ReplacementRule(pattern4700, replacement4700)
    pattern4701 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1648)
    rule4701 = ReplacementRule(pattern4701, replacement4701)
    pattern4702 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1648)
    rule4702 = ReplacementRule(pattern4702, replacement4702)
    pattern4703 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1648)
    rule4703 = ReplacementRule(pattern4703, replacement4703)
    pattern4704 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1648)
    rule4704 = ReplacementRule(pattern4704, replacement4704)
    pattern4705 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1648)
    rule4705 = ReplacementRule(pattern4705, replacement4705)
    pattern4706 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1648)
    rule4706 = ReplacementRule(pattern4706, replacement4706)
    pattern4707 = Pattern(Integral(u_*(A_ + WC('C', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1648)
    rule4707 = ReplacementRule(pattern4707, replacement4707)
    pattern4708 = Pattern(Integral(u_*(A_ + WC('C', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1648)
    rule4708 = ReplacementRule(pattern4708, replacement4708)
    pattern4709 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons38, cons1649)
    rule4709 = ReplacementRule(pattern4709, replacement4709)
    pattern4710 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons38, cons1649)
    rule4710 = ReplacementRule(pattern4710, replacement4710)
    pattern4711 = Pattern(Integral(u_*(WC('A', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4711 = ReplacementRule(pattern4711, replacement4711)
    pattern4712 = Pattern(Integral(u_*(WC('A', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4712 = ReplacementRule(pattern4712, replacement4712)
    pattern4713 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1652)
    rule4713 = ReplacementRule(pattern4713, replacement4713)
    pattern4714 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1653)
    rule4714 = ReplacementRule(pattern4714, replacement4714)
    pattern4715 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1652)
    rule4715 = ReplacementRule(pattern4715, replacement4715)
    pattern4716 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1653)
    rule4716 = ReplacementRule(pattern4716, replacement4716)
    pattern4717 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1652)
    rule4717 = ReplacementRule(pattern4717, replacement4717)
    pattern4718 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1653)
    rule4718 = ReplacementRule(pattern4718, replacement4718)
    pattern4719 = Pattern(Integral(u_*(A_ + WC('B', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1652)
    rule4719 = ReplacementRule(pattern4719, replacement4719)
    pattern4720 = Pattern(Integral(u_*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1653)
    rule4720 = ReplacementRule(pattern4720, replacement4720)
    pattern4721 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1652)
    rule4721 = ReplacementRule(pattern4721, replacement4721)
    pattern4722 = Pattern(Integral((WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1653)
    rule4722 = ReplacementRule(pattern4722, replacement4722)
    pattern4723 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1652)
    rule4723 = ReplacementRule(pattern4723, replacement4723)
    pattern4724 = Pattern(Integral((WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1653)
    rule4724 = ReplacementRule(pattern4724, replacement4724)
    pattern4725 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1652)
    rule4725 = ReplacementRule(pattern4725, replacement4725)
    pattern4726 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1653)
    rule4726 = ReplacementRule(pattern4726, replacement4726)
    pattern4727 = Pattern(Integral(u_*(A_ + WC('C', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1652)
    rule4727 = ReplacementRule(pattern4727, replacement4727)
    pattern4728 = Pattern(Integral(u_*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1653)
    rule4728 = ReplacementRule(pattern4728, replacement4728)
    pattern4729 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons38, cons1649)
    rule4729 = ReplacementRule(pattern4729, replacement4729)
    pattern4730 = Pattern(Integral(u_*(WC('A', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4730 = ReplacementRule(pattern4730, replacement4730)
    pattern4731 = Pattern(Integral(u_*((S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**n1_*WC('B', S(1)) + (S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**n2_*WC('C', S(1)) + (S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*WC('A', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4731 = ReplacementRule(pattern4731, replacement4731)
    pattern4732 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654)
    rule4732 = ReplacementRule(pattern4732, replacement4732)
    pattern4733 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654)
    rule4733 = ReplacementRule(pattern4733, replacement4733)
    pattern4734 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654, cons21)
    rule4734 = ReplacementRule(pattern4734, replacement4734)
    pattern4735 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654, cons21)
    rule4735 = ReplacementRule(pattern4735, replacement4735)
    pattern4736 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654, cons21)
    rule4736 = ReplacementRule(pattern4736, replacement4736)
    pattern4737 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1654, cons21)
    rule4737 = ReplacementRule(pattern4737, replacement4737)
    pattern4738 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1654)
    rule4738 = ReplacementRule(pattern4738, replacement4738)
    pattern4739 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1654)
    rule4739 = ReplacementRule(pattern4739, replacement4739)
    pattern4740 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1654)
    rule4740 = ReplacementRule(pattern4740, replacement4740)
    pattern4741 = Pattern(Integral(u_*(WC('c', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons19, cons21, cons1654)
    rule4741 = ReplacementRule(pattern4741, replacement4741)
    pattern4742 = Pattern(Integral(u_*(WC('c', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1654)
    rule4742 = ReplacementRule(pattern4742, replacement4742)
    pattern4743 = Pattern(Integral(u_*(WC('c', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons1654)
    rule4743 = ReplacementRule(pattern4743, replacement4743)
    pattern4744 = Pattern(Integral(u_*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1654)
    rule4744 = ReplacementRule(pattern4744, replacement4744)
    pattern4745 = Pattern(Integral(u_*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons36, cons37, cons1654)
    rule4745 = ReplacementRule(pattern4745, replacement4745)
    pattern4746 = Pattern(Integral((WC('c', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1654)
    rule4746 = ReplacementRule(pattern4746, replacement4746)
    pattern4747 = Pattern(Integral((WC('c', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons38, cons4, cons1654)
    rule4747 = ReplacementRule(pattern4747, replacement4747)
    pattern4748 = Pattern(Integral((WC('c', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1654)
    rule4748 = ReplacementRule(pattern4748, replacement4748)
    pattern4749 = Pattern(Integral((WC('c', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons38, cons4, cons1654)
    rule4749 = ReplacementRule(pattern4749, replacement4749)
    pattern4750 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1654)
    rule4750 = ReplacementRule(pattern4750, replacement4750)
    pattern4751 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons37, cons38, cons1654)
    rule4751 = ReplacementRule(pattern4751, replacement4751)
    pattern4752 = Pattern(Integral(u_*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1654)
    rule4752 = ReplacementRule(pattern4752, replacement4752)
    pattern4753 = Pattern(Integral(u_*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons36, cons38, cons1654)
    rule4753 = ReplacementRule(pattern4753, replacement4753)
    pattern4754 = Pattern(Integral(u_*((S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**n1_*WC('B', S(1)) + (S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**n2_*WC('C', S(1)) + (S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*WC('A', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4754 = ReplacementRule(pattern4754, replacement4754)
    pattern4755 = Pattern(Integral(u_*((S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**n1_*WC('B', S(1)) + (S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**n2_*WC('C', S(1)) + (S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*WC('A', S(1))), x_), cons2, cons3, cons36, cons37, cons38, cons4, cons1650, cons1651)
    rule4755 = ReplacementRule(pattern4755, replacement4755)
    pattern4756 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1655)
    rule4756 = ReplacementRule(pattern4756, replacement4756)
    pattern4757 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1655)
    rule4757 = ReplacementRule(pattern4757, replacement4757)
    pattern4758 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1655)
    rule4758 = ReplacementRule(pattern4758, replacement4758)
    pattern4759 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons1410)
    rule4759 = ReplacementRule(pattern4759, replacement4759)
    pattern4760 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons1410)
    rule4760 = ReplacementRule(pattern4760, replacement4760)
    pattern4761 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons27, cons1656, cons40)
    rule4761 = ReplacementRule(pattern4761, replacement4761)
    pattern4762 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons4, cons27, cons1656, cons40)
    rule4762 = ReplacementRule(pattern4762, replacement4762)
    pattern4763 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons5, cons27, cons1656, cons149, cons345)
    rule4763 = ReplacementRule(pattern4763, replacement4763)
    pattern4764 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons5, cons27, cons1656, cons149, cons345)
    rule4764 = ReplacementRule(pattern4764, replacement4764)
    pattern4765 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons5, cons27, cons1656, cons149, cons242)
    rule4765 = ReplacementRule(pattern4765, replacement4765)
    pattern4766 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons5, cons27, cons1656, cons149, cons242)
    rule4766 = ReplacementRule(pattern4766, replacement4766)
    pattern4767 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons27, cons1656, cons149, cons246, cons1335, cons139, cons1657, cons1290)
    rule4767 = ReplacementRule(pattern4767, replacement4767)
    pattern4768 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons27, cons1656, cons149, cons246, cons1335, cons139, cons1657, cons1290)
    rule4768 = ReplacementRule(pattern4768, replacement4768)
    pattern4769 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons27, cons1656, cons149, cons246, cons168, cons139, cons321, cons1658, cons1290)
    rule4769 = ReplacementRule(pattern4769, replacement4769)
    pattern4770 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons27, cons1656, cons149, cons246, cons168, cons139, cons321, cons1658, cons1290)
    rule4770 = ReplacementRule(pattern4770, replacement4770)
    pattern4771 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons5, cons27, cons1656, cons149, cons33, cons168, cons251, cons1290)
    rule4771 = ReplacementRule(pattern4771, replacement4771)
    pattern4772 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons5, cons27, cons1656, cons149, cons33, cons168, cons251, cons1290)
    rule4772 = ReplacementRule(pattern4772, replacement4772)
    pattern4773 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons5, cons27, cons1656, cons149, cons33, cons96, cons321, cons255, cons1290)
    rule4773 = ReplacementRule(pattern4773, replacement4773)
    pattern4774 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons5, cons27, cons1656, cons149, cons33, cons96, cons321, cons255, cons1290)
    rule4774 = ReplacementRule(pattern4774, replacement4774)
    pattern4775 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons149, cons13, cons165, cons248)
    rule4775 = ReplacementRule(pattern4775, replacement4775)
    pattern4776 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons149, cons13, cons165, cons248)
    rule4776 = ReplacementRule(pattern4776, replacement4776)
    pattern4777 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons149, cons13, cons139, cons248)
    rule4777 = ReplacementRule(pattern4777, replacement4777)
    pattern4778 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons149, cons13, cons139, cons248)
    rule4778 = ReplacementRule(pattern4778, replacement4778)
    pattern4779 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons27, cons1656)
    rule4779 = ReplacementRule(pattern4779, replacement4779)
    pattern4780 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons27, cons1656)
    rule4780 = ReplacementRule(pattern4780, replacement4780)
    pattern4781 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons5, cons27, cons1656, cons149, cons248)
    rule4781 = ReplacementRule(pattern4781, replacement4781)
    pattern4782 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons210, cons5, cons27, cons1656, cons149, cons248)
    rule4782 = ReplacementRule(pattern4782, replacement4782)
    pattern4783 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons210, cons19, cons5, cons27, cons1656, cons149)
    rule4783 = ReplacementRule(pattern4783, replacement4783)
    pattern4784 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons127, cons210, cons4, cons5, cons27, cons1656, cons149)
    rule4784 = ReplacementRule(pattern4784, replacement4784)
    pattern4785 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons210, cons27, cons1656, cons1410)
    rule4785 = ReplacementRule(pattern4785, replacement4785)
    pattern4786 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons27, cons1656, cons40)
    rule4786 = ReplacementRule(pattern4786, replacement4786)
    pattern4787 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons27, cons1656, cons149, cons1280)
    rule4787 = ReplacementRule(pattern4787, replacement4787)
    pattern4788 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons27, cons1656, cons149, cons1280)
    rule4788 = ReplacementRule(pattern4788, replacement4788)
    pattern4789 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons27, cons1656, cons149, cons1659, cons255)
    rule4789 = ReplacementRule(pattern4789, replacement4789)
    pattern4790 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons27, cons1656, cons149, cons246, cons1660, cons139, cons1661, cons172)
    rule4790 = ReplacementRule(pattern4790, replacement4790)
    pattern4791 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons27, cons1656, cons149, cons246, cons1660, cons139, cons1661, cons172)
    rule4791 = ReplacementRule(pattern4791, replacement4791)
    pattern4792 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons27, cons1656, cons149, cons246, cons168, cons139, cons1662, cons1661, cons172, cons1663)
    rule4792 = ReplacementRule(pattern4792, replacement4792)
    pattern4793 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons27, cons1656, cons149, cons246, cons168, cons139, cons1662, cons1661, cons172, cons1663)
    rule4793 = ReplacementRule(pattern4793, replacement4793)
    pattern4794 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons27, cons1656, cons149, cons95, cons168, cons91, cons1661, cons172)
    rule4794 = ReplacementRule(pattern4794, replacement4794)
    pattern4795 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons27, cons1656, cons149, cons95, cons168, cons91, cons1661, cons172)
    rule4795 = ReplacementRule(pattern4795, replacement4795)
    pattern4796 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons27, cons1656, cons149, cons33, cons168, cons1664, cons172)
    rule4796 = ReplacementRule(pattern4796, replacement4796)
    pattern4797 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons27, cons1656, cons149, cons33, cons168, cons1664, cons172)
    rule4797 = ReplacementRule(pattern4797, replacement4797)
    pattern4798 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons27, cons1656, cons149, cons164, cons96, cons90, cons165, cons1664, cons172)
    rule4798 = ReplacementRule(pattern4798, replacement4798)
    pattern4799 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons27, cons1656, cons149, cons164, cons96, cons90, cons165, cons1664, cons172)
    rule4799 = ReplacementRule(pattern4799, replacement4799)
    pattern4800 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons27, cons1656, cons149, cons164, cons96, cons90, cons139, cons1662, cons255, cons172)
    rule4800 = ReplacementRule(pattern4800, replacement4800)
    pattern4801 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons27, cons1656, cons149, cons164, cons96, cons90, cons139, cons1662, cons255, cons172)
    rule4801 = ReplacementRule(pattern4801, replacement4801)
    pattern4802 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons27, cons1656, cons149, cons33, cons96, cons1662, cons255, cons172)
    rule4802 = ReplacementRule(pattern4802, replacement4802)
    pattern4803 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons27, cons1656, cons149, cons33, cons96, cons1662, cons255, cons172)
    rule4803 = ReplacementRule(pattern4803, replacement4803)
    pattern4804 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons27, cons1656, cons149)
    rule4804 = ReplacementRule(pattern4804, replacement4804)
    pattern4805 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons27, cons1665)
    rule4805 = ReplacementRule(pattern4805, replacement4805)
    pattern4806 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_)**p_, x_), cons2, cons3, cons8, cons29, cons1666, cons87, cons130)
    rule4806 = ReplacementRule(pattern4806, replacement4806)
    pattern4807 = Pattern(Integral(S(1)/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons8, cons29, cons1666, cons1481, cons746)
    rule4807 = ReplacementRule(pattern4807, replacement4807)
    pattern4808 = Pattern(Integral(S(1)/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons8, cons29, cons1666, cons1484, cons746)
    rule4808 = ReplacementRule(pattern4808, replacement4808)
    pattern4809 = Pattern(Integral(G_**(x_*WC('d', S(1)) + WC('c', S(0)))/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons8, cons29, cons19, cons1667, cons87, cons746)
    rule4809 = ReplacementRule(pattern4809, replacement4809)

    pattern4810 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_, x_), cons2, cons8, cons29, cons4, cons5, cons1666, cons25, cons40)
    rule4810 = ReplacementRule(pattern4810, With4810)

    pattern4811 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1666, cons25, cons149)
    rule4811 = ReplacementRule(pattern4811, With4811)

    pattern4812 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1668, CustomConstraint(With4812))
    rule4812 = ReplacementRule(pattern4812, replacement4812)

    pattern4813 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1669, CustomConstraint(With4813))
    rule4813 = ReplacementRule(pattern4813, replacement4813)

    pattern4814 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1670, CustomConstraint(With4814))
    rule4814 = ReplacementRule(pattern4814, replacement4814)

    pattern4815 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1671, CustomConstraint(With4815))
    rule4815 = ReplacementRule(pattern4815, replacement4815)

    pattern4816 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1249, cons1672, CustomConstraint(With4816))
    rule4816 = ReplacementRule(pattern4816, replacement4816)

    pattern4817 = Pattern(Integral(u_/cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), cons2, cons3, cons8, cons1249, CustomConstraint(With4817))
    rule4817 = ReplacementRule(pattern4817, replacement4817)

    pattern4818 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1249, cons1673, CustomConstraint(With4818))
    rule4818 = ReplacementRule(pattern4818, replacement4818)

    pattern4819 = Pattern(Integral(u_/sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), cons2, cons3, cons8, cons1249, CustomConstraint(With4819))
    rule4819 = ReplacementRule(pattern4819, replacement4819)

    pattern4820 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons87, cons1670, CustomConstraint(With4820))
    rule4820 = ReplacementRule(pattern4820, replacement4820)

    pattern4821 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons87, cons1671, CustomConstraint(With4821))
    rule4821 = ReplacementRule(pattern4821, replacement4821)

    pattern4822 = Pattern(Integral(u_, x_), CustomConstraint(With4822))
    rule4822 = ReplacementRule(pattern4822, replacement4822)

    pattern4823 = Pattern(Integral(u_, x_), CustomConstraint(With4823))
    rule4823 = ReplacementRule(pattern4823, replacement4823)
    pattern4824 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1674, cons1675, cons557)
    rule4824 = ReplacementRule(pattern4824, replacement4824)
    pattern4825 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*H_**(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1674, cons1675, cons1676, cons630)
    rule4825 = ReplacementRule(pattern4825, replacement4825)

    pattern4826 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1668, CustomConstraint(With4826))
    rule4826 = ReplacementRule(pattern4826, replacement4826)

    pattern4827 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1669, CustomConstraint(With4827))
    rule4827 = ReplacementRule(pattern4827, replacement4827)

    pattern4828 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1670, CustomConstraint(With4828))
    rule4828 = ReplacementRule(pattern4828, replacement4828)

    pattern4829 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons8, cons1671, CustomConstraint(With4829))
    rule4829 = ReplacementRule(pattern4829, replacement4829)

    pattern4830 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1668, CustomConstraint(With4830))
    rule4830 = ReplacementRule(pattern4830, replacement4830)

    pattern4831 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1672, CustomConstraint(With4831))
    rule4831 = ReplacementRule(pattern4831, replacement4831)

    pattern4832 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1669, CustomConstraint(With4832))
    rule4832 = ReplacementRule(pattern4832, replacement4832)

    pattern4833 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1673, CustomConstraint(With4833))
    rule4833 = ReplacementRule(pattern4833, replacement4833)

    pattern4834 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1670, CustomConstraint(With4834))
    rule4834 = ReplacementRule(pattern4834, replacement4834)

    pattern4835 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons8, cons1484, cons1249, cons1671, CustomConstraint(With4835))
    rule4835 = ReplacementRule(pattern4835, replacement4835)

    pattern4836 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), cons2, cons3, cons8, cons29, cons10, cons1484, cons1249, cons1668, CustomConstraint(With4836))
    rule4836 = ReplacementRule(pattern4836, replacement4836)

    pattern4837 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), cons2, cons3, cons8, cons29, cons10, cons1484, cons1249, cons1669, CustomConstraint(With4837))
    rule4837 = ReplacementRule(pattern4837, replacement4837)

    pattern4838 = Pattern(Integral(u_, x_), CustomConstraint(With4838))
    rule4838 = ReplacementRule(pattern4838, replacement4838)

    pattern4839 = Pattern(Integral(u_, x_), CustomConstraint(With4839))
    rule4839 = ReplacementRule(pattern4839, replacement4839)
    pattern4840 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1677)
    rule4840 = ReplacementRule(pattern4840, replacement4840)
    pattern4841 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1678)
    rule4841 = ReplacementRule(pattern4841, replacement4841)
    pattern4842 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1678)
    rule4842 = ReplacementRule(pattern4842, replacement4842)

    pattern4843 = Pattern(Integral(u_/y_, x_), cons1679, CustomConstraint(With4843))
    rule4843 = ReplacementRule(pattern4843, replacement4843)

    pattern4844 = Pattern(Integral(u_/(w_*y_), x_), cons1679, CustomConstraint(With4844))
    rule4844 = ReplacementRule(pattern4844, replacement4844)

    pattern4845 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons19, cons68, cons1679, CustomConstraint(With4845))
    rule4845 = ReplacementRule(pattern4845, replacement4845)

    pattern4846 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons19, cons4, cons68, cons1679, CustomConstraint(With4846))
    rule4846 = ReplacementRule(pattern4846, replacement4846)

    pattern4847 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_*WC('u', S(1)), x_), cons2, cons8, cons29, cons4, cons5, cons1666, cons25, cons40)
    rule4847 = ReplacementRule(pattern4847, With4847)

    pattern4848 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1666, cons25, cons149)
    rule4848 = ReplacementRule(pattern4848, With4848)

    pattern4849 = Pattern(Integral(u_, x_), cons1232, CustomConstraint(With4849))
    rule4849 = ReplacementRule(pattern4849, replacement4849)
    pattern4850 = Pattern(Integral(((S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)))**p_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons378)
    rule4850 = ReplacementRule(pattern4850, replacement4850)
    pattern4851 = Pattern(Integral(((S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons378)
    rule4851 = ReplacementRule(pattern4851, replacement4851)
    pattern4852 = Pattern(Integral(u_*(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*a_ + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons52, cons1666, cons87, cons51)
    rule4852 = ReplacementRule(pattern4852, replacement4852)
    pattern4853 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*a_ + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons54, cons1666, cons87, cons51, cons53)
    rule4853 = ReplacementRule(pattern4853, replacement4853)
    pattern4854 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)) + a_)**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons52, cons1666, cons87, cons1680)
    rule4854 = ReplacementRule(pattern4854, replacement4854)
    pattern4855 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1441)
    rule4855 = ReplacementRule(pattern4855, replacement4855)
    pattern4856 = Pattern(Integral(u_, x_), cons1681)
    rule4856 = ReplacementRule(pattern4856, replacement4856)

    pattern4857 = Pattern(Integral((a_*v_)**p_*WC('u', S(1)), x_), cons2, cons5, cons149, cons1682)
    rule4857 = ReplacementRule(pattern4857, With4857)

    pattern4858 = Pattern(Integral((v_**m_)**p_*WC('u', S(1)), x_), cons19, cons5, cons149, cons1682)
    rule4858 = ReplacementRule(pattern4858, With4858)

    pattern4859 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1)))**p_*WC('u', S(1)), x_), cons19, cons4, cons5, cons149, cons1683)
    rule4859 = ReplacementRule(pattern4859, With4859)

    pattern4860 = Pattern(Integral(u_, x_), cons1679, CustomConstraint(With4860))
    rule4860 = ReplacementRule(pattern4860, replacement4860)

    pattern4861 = Pattern(Integral(u_, x_), cons1232, cons1684, CustomConstraint(With4861))
    rule4861 = ReplacementRule(pattern4861, replacement4861)

    pattern4862 = Pattern(Integral(u_, x_), cons1679)
    rule4862 = ReplacementRule(pattern4862, With4862)
    pattern4863 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule4863 = ReplacementRule(pattern4863, replacement4863)
    pattern4864 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule4864 = ReplacementRule(pattern4864, replacement4864)
    pattern4865 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule4865 = ReplacementRule(pattern4865, replacement4865)
    pattern4866 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule4866 = ReplacementRule(pattern4866, replacement4866)
    pattern4867 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons466)
    rule4867 = ReplacementRule(pattern4867, replacement4867)
    pattern4868 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1685, cons33, cons170)
    rule4868 = ReplacementRule(pattern4868, replacement4868)
    pattern4869 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1685, cons33, cons170)
    rule4869 = ReplacementRule(pattern4869, replacement4869)
    pattern4870 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule4870 = ReplacementRule(pattern4870, replacement4870)
    pattern4871 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))/sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons4, cons64, cons586)
    rule4871 = ReplacementRule(pattern4871, replacement4871)
    pattern4872 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_/cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons1410)
    rule4872 = ReplacementRule(pattern4872, replacement4872)
    pattern4873 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1410)
    rule4873 = ReplacementRule(pattern4873, replacement4873)
    pattern4874 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**p_/sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons1410)
    rule4874 = ReplacementRule(pattern4874, replacement4874)
    pattern4875 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1410)
    rule4875 = ReplacementRule(pattern4875, replacement4875)

    pattern4876 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons64, cons1686)
    rule4876 = ReplacementRule(pattern4876, With4876)

    pattern4877 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons64, cons1686)
    rule4877 = ReplacementRule(pattern4877, With4877)
    pattern4878 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons33, cons87)
    rule4878 = ReplacementRule(pattern4878, replacement4878)

    pattern4879 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(S(1)/sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(S(1)/cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons378, cons33, cons170, cons1687)
    rule4879 = ReplacementRule(pattern4879, With4879)
    pattern4880 = Pattern(Integral(F_**v_*G_**w_*u_**WC('m', S(1)), x_), cons19, cons4, cons5, cons1688, cons1689, cons1690, cons814, cons815)
    rule4880 = ReplacementRule(pattern4880, replacement4880)
    pattern4881 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4881 = ReplacementRule(pattern4881, replacement4881)
    pattern4882 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4882 = ReplacementRule(pattern4882, replacement4882)
    pattern4883 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4883 = ReplacementRule(pattern4883, replacement4883)
    pattern4884 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4884 = ReplacementRule(pattern4884, replacement4884)
    pattern4885 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))/cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4885 = ReplacementRule(pattern4885, replacement4885)
    pattern4886 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))/(sin(x_*WC('d', S(1)) + WC('c', S(0)))*tan(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons586)
    rule4886 = ReplacementRule(pattern4886, replacement4886)
    pattern4887 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons557, cons20)
    rule4887 = ReplacementRule(pattern4887, replacement4887)
    pattern4888 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons557, cons20)
    rule4888 = ReplacementRule(pattern4888, replacement4888)
    pattern4889 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons557)
    rule4889 = ReplacementRule(pattern4889, replacement4889)
    pattern4890 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1691, cons1692, cons557, cons27, cons1693)
    rule4890 = ReplacementRule(pattern4890, replacement4890)
    pattern4891 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1694)
    rule4891 = ReplacementRule(pattern4891, replacement4891)
    pattern4892 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1694)
    rule4892 = ReplacementRule(pattern4892, replacement4892)
    pattern4893 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1695, cons89, cons167)
    rule4893 = ReplacementRule(pattern4893, replacement4893)
    pattern4894 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1696, cons33, cons168)
    rule4894 = ReplacementRule(pattern4894, replacement4894)
    pattern4895 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1697, cons586, cons1397)
    rule4895 = ReplacementRule(pattern4895, replacement4895)
    pattern4896 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1697, cons586, cons1397)
    rule4896 = ReplacementRule(pattern4896, replacement4896)
    pattern4897 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1698, cons89, cons91, cons1444)
    rule4897 = ReplacementRule(pattern4897, replacement4897)
    pattern4898 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1698, cons89, cons91, cons1444)
    rule4898 = ReplacementRule(pattern4898, replacement4898)
    pattern4899 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons25)
    rule4899 = ReplacementRule(pattern4899, replacement4899)
    pattern4900 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons25)
    rule4900 = ReplacementRule(pattern4900, replacement4900)
    pattern4901 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule4901 = ReplacementRule(pattern4901, replacement4901)
    pattern4902 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule4902 = ReplacementRule(pattern4902, replacement4902)
    pattern4903 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1695, cons89, cons91)
    rule4903 = ReplacementRule(pattern4903, replacement4903)
    pattern4904 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1695, cons89, cons91)
    rule4904 = ReplacementRule(pattern4904, replacement4904)
    pattern4905 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1699, cons1504, cons965)
    rule4905 = ReplacementRule(pattern4905, replacement4905)
    pattern4906 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons1699, cons1504, cons965)
    rule4906 = ReplacementRule(pattern4906, replacement4906)
    pattern4907 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1700, cons89, cons167, cons1646)
    rule4907 = ReplacementRule(pattern4907, replacement4907)
    pattern4908 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1700, cons89, cons167, cons1646)
    rule4908 = ReplacementRule(pattern4908, replacement4908)
    pattern4909 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule4909 = ReplacementRule(pattern4909, replacement4909)
    pattern4910 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons87)
    rule4910 = ReplacementRule(pattern4910, replacement4910)
    pattern4911 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons25)
    rule4911 = ReplacementRule(pattern4911, replacement4911)
    pattern4912 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons25)
    rule4912 = ReplacementRule(pattern4912, replacement4912)
    pattern4913 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1701, cons198)
    rule4913 = ReplacementRule(pattern4913, replacement4913)
    pattern4914 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1702, cons198)
    rule4914 = ReplacementRule(pattern4914, replacement4914)
    pattern4915 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1013, cons198)
    rule4915 = ReplacementRule(pattern4915, replacement4915)
    pattern4916 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1701, cons152, cons1553)
    rule4916 = ReplacementRule(pattern4916, replacement4916)
    pattern4917 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1702, cons152, cons1553)
    rule4917 = ReplacementRule(pattern4917, replacement4917)
    pattern4918 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1013, cons152, cons1553)
    rule4918 = ReplacementRule(pattern4918, replacement4918)
    pattern4919 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons1701, cons1703, cons1704)
    rule4919 = ReplacementRule(pattern4919, replacement4919)
    pattern4920 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons1701, cons1703, cons1705)
    rule4920 = ReplacementRule(pattern4920, replacement4920)
    pattern4921 = Pattern(Integral(F_**(u_*WC('c', S(1)))*G_**v_, x_), cons1101, cons8, cons4, cons1689, cons812, cons813)
    rule4921 = ReplacementRule(pattern4921, replacement4921)

    pattern4922 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons33, cons170, cons150)
    rule4922 = ReplacementRule(pattern4922, With4922)

    pattern4923 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons33, cons170, cons150)
    rule4923 = ReplacementRule(pattern4923, With4923)
    pattern4924 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons530)
    rule4924 = ReplacementRule(pattern4924, replacement4924)
    pattern4925 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1706)
    rule4925 = ReplacementRule(pattern4925, replacement4925)
    pattern4926 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*G_**(x_*WC('e', S(1)) + WC('d', S(0)))*H_**(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons530, cons1689, cons1707)
    rule4926 = ReplacementRule(pattern4926, replacement4926)
    pattern4927 = Pattern(Integral(F_**u_*sin(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons150)
    rule4927 = ReplacementRule(pattern4927, replacement4927)
    pattern4928 = Pattern(Integral(F_**u_*cos(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons150)
    rule4928 = ReplacementRule(pattern4928, replacement4928)
    pattern4929 = Pattern(Integral(F_**u_*sin(v_)**WC('m', S(1))*cos(v_)**WC('n', S(1)), x_), cons1101, cons1708, cons1709, cons530)
    rule4929 = ReplacementRule(pattern4929, replacement4929)
    pattern4930 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1710, cons56)
    rule4930 = ReplacementRule(pattern4930, replacement4930)
    pattern4931 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1710, cons56)
    rule4931 = ReplacementRule(pattern4931, replacement4931)
    pattern4932 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons130, cons1711)
    rule4932 = ReplacementRule(pattern4932, replacement4932)
    pattern4933 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons130, cons1711)
    rule4933 = ReplacementRule(pattern4933, replacement4933)
    pattern4934 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1712)
    rule4934 = ReplacementRule(pattern4934, replacement4934)
    pattern4935 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1712)
    rule4935 = ReplacementRule(pattern4935, replacement4935)
    pattern4936 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1713)
    rule4936 = ReplacementRule(pattern4936, replacement4936)
    pattern4937 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1713)
    rule4937 = ReplacementRule(pattern4937, replacement4937)
    pattern4938 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1507, cons1714)
    rule4938 = ReplacementRule(pattern4938, replacement4938)
    pattern4939 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1507, cons1714)
    rule4939 = ReplacementRule(pattern4939, replacement4939)
    pattern4940 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1713)
    rule4940 = ReplacementRule(pattern4940, replacement4940)
    pattern4941 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1713)
    rule4941 = ReplacementRule(pattern4941, replacement4941)
    pattern4942 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1715, cons56, cons68)
    rule4942 = ReplacementRule(pattern4942, replacement4942)
    pattern4943 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1715, cons56, cons68)
    rule4943 = ReplacementRule(pattern4943, replacement4943)
    pattern4944 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons130, cons1716)
    rule4944 = ReplacementRule(pattern4944, replacement4944)
    pattern4945 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons130, cons1716)
    rule4945 = ReplacementRule(pattern4945, replacement4945)
    pattern4946 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1717)
    rule4946 = ReplacementRule(pattern4946, replacement4946)
    pattern4947 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1717)
    rule4947 = ReplacementRule(pattern4947, replacement4947)
    pattern4948 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1718)
    rule4948 = ReplacementRule(pattern4948, replacement4948)
    pattern4949 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1718)
    rule4949 = ReplacementRule(pattern4949, replacement4949)
    pattern4950 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1507, cons1719)
    rule4950 = ReplacementRule(pattern4950, replacement4950)
    pattern4951 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1507, cons1719)
    rule4951 = ReplacementRule(pattern4951, replacement4951)
    pattern4952 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1718)
    rule4952 = ReplacementRule(pattern4952, replacement4952)
    pattern4953 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1718)
    rule4953 = ReplacementRule(pattern4953, replacement4953)
    pattern4954 = Pattern(Integral(S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1720)
    rule4954 = ReplacementRule(pattern4954, replacement4954)
    pattern4955 = Pattern(Integral(S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1720)
    rule4955 = ReplacementRule(pattern4955, replacement4955)
    pattern4956 = Pattern(Integral((S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1721, cons1647)
    rule4956 = ReplacementRule(pattern4956, replacement4956)
    pattern4957 = Pattern(Integral((S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons1721, cons1647)
    rule4957 = ReplacementRule(pattern4957, replacement4957)
    pattern4958 = Pattern(Integral((S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1722, cons1723)
    rule4958 = ReplacementRule(pattern4958, replacement4958)
    pattern4959 = Pattern(Integral((S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons148, cons1722, cons1723)
    rule4959 = ReplacementRule(pattern4959, replacement4959)
    pattern4960 = Pattern(Integral((S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1713)
    rule4960 = ReplacementRule(pattern4960, replacement4960)
    pattern4961 = Pattern(Integral((S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons4, cons13, cons139, cons1713)
    rule4961 = ReplacementRule(pattern4961, replacement4961)
    pattern4962 = Pattern(Integral((S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons1713)
    rule4962 = ReplacementRule(pattern4962, replacement4962)
    pattern4963 = Pattern(Integral((S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons1713)
    rule4963 = ReplacementRule(pattern4963, replacement4963)
    pattern4964 = Pattern(Integral(x_**WC('m', S(1))/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1724)
    rule4964 = ReplacementRule(pattern4964, replacement4964)
    pattern4965 = Pattern(Integral(x_**WC('m', S(1))/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons4, cons1724)
    rule4965 = ReplacementRule(pattern4965, replacement4965)
    pattern4966 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1725, cons68, cons1647)
    rule4966 = ReplacementRule(pattern4966, replacement4966)
    pattern4967 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1725, cons68, cons1647)
    rule4967 = ReplacementRule(pattern4967, replacement4967)
    pattern4968 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1722, cons1726)
    rule4968 = ReplacementRule(pattern4968, replacement4968)
    pattern4969 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons148, cons1722, cons1726)
    rule4969 = ReplacementRule(pattern4969, replacement4969)
    pattern4970 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1718)
    rule4970 = ReplacementRule(pattern4970, replacement4970)
    pattern4971 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**p_, x_), cons2, cons3, cons8, cons19, cons4, cons13, cons139, cons1718)
    rule4971 = ReplacementRule(pattern4971, replacement4971)
    pattern4972 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1718)
    rule4972 = ReplacementRule(pattern4972, replacement4972)
    pattern4973 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons1718)
    rule4973 = ReplacementRule(pattern4973, replacement4973)
    pattern4974 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons165)
    rule4974 = ReplacementRule(pattern4974, replacement4974)
    pattern4975 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons165)
    rule4975 = ReplacementRule(pattern4975, replacement4975)
    pattern4976 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons340, cons165)
    rule4976 = ReplacementRule(pattern4976, replacement4976)
    pattern4977 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons340, cons165)
    rule4977 = ReplacementRule(pattern4977, replacement4977)
    pattern4978 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons55, cons13, cons165)
    rule4978 = ReplacementRule(pattern4978, replacement4978)
    pattern4979 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons55, cons13, cons165)
    rule4979 = ReplacementRule(pattern4979, replacement4979)
    pattern4980 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons13, cons165, cons629)
    rule4980 = ReplacementRule(pattern4980, replacement4980)
    pattern4981 = Pattern(Integral(x_**m_*log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons19, cons4, cons13, cons165, cons629)
    rule4981 = ReplacementRule(pattern4981, replacement4981)
    pattern4982 = Pattern(Integral(sin(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons8, cons29, cons150)
    rule4982 = ReplacementRule(pattern4982, replacement4982)
    pattern4983 = Pattern(Integral(cos(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons8, cons29, cons150)
    rule4983 = ReplacementRule(pattern4983, replacement4983)
    pattern4984 = Pattern(Integral(sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150, cons73)
    rule4984 = ReplacementRule(pattern4984, replacement4984)
    pattern4985 = Pattern(Integral(cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150, cons73)
    rule4985 = ReplacementRule(pattern4985, replacement4985)

    pattern4986 = Pattern(Integral(sin(u_)**WC('n', S(1)), x_), cons150, cons1727)
    rule4986 = ReplacementRule(pattern4986, With4986)

    pattern4987 = Pattern(Integral(cos(u_)**WC('n', S(1)), x_), cons150, cons1727)
    rule4987 = ReplacementRule(pattern4987, With4987)
    pattern4988 = Pattern(Integral(WC('u', S(1))*sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1690)
    rule4988 = ReplacementRule(pattern4988, replacement4988)
    pattern4989 = Pattern(Integral(WC('u', S(1))*cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1690)
    rule4989 = ReplacementRule(pattern4989, replacement4989)
    pattern4990 = Pattern(Integral(sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1728, cons557)
    rule4990 = ReplacementRule(pattern4990, replacement4990)
    pattern4991 = Pattern(Integral(cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1728, cons557)
    rule4991 = ReplacementRule(pattern4991, replacement4991)
    pattern4992 = Pattern(Integral(x_**WC('m', S(1))*sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule4992 = ReplacementRule(pattern4992, replacement4992)
    pattern4993 = Pattern(Integral(x_**WC('m', S(1))*cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule4993 = ReplacementRule(pattern4993, replacement4993)
    pattern4994 = Pattern(Integral(WC('u', S(1))*sin(v_)**WC('p', S(1))*cos(w_)**WC('p', S(1)), x_), cons1690, cons40)
    rule4994 = ReplacementRule(pattern4994, replacement4994)
    pattern4995 = Pattern(Integral(sin(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons557, cons1728)
    rule4995 = ReplacementRule(pattern4995, replacement4995)
    pattern4996 = Pattern(Integral(x_**WC('m', S(1))*sin(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1729, cons1728)
    rule4996 = ReplacementRule(pattern4996, replacement4996)
    pattern4997 = Pattern(Integral(sin(v_)*tan(w_)**WC('n', S(1)), x_), cons89, cons90, cons1730)
    rule4997 = ReplacementRule(pattern4997, replacement4997)
    pattern4998 = Pattern(Integral((S(1)/tan(w_))**WC('n', S(1))*cos(v_), x_), cons89, cons90, cons1730)
    rule4998 = ReplacementRule(pattern4998, replacement4998)
    pattern4999 = Pattern(Integral((S(1)/tan(w_))**WC('n', S(1))*sin(v_), x_), cons89, cons90, cons1730)
    rule4999 = ReplacementRule(pattern4999, replacement4999)
    pattern5000 = Pattern(Integral(cos(v_)*tan(w_)**WC('n', S(1)), x_), cons89, cons90, cons1730)
    rule5000 = ReplacementRule(pattern5000, replacement5000)
    pattern5001 = Pattern(Integral((S(1)/cos(w_))**WC('n', S(1))*sin(v_), x_), cons89, cons90, cons1730)
    rule5001 = ReplacementRule(pattern5001, replacement5001)
    pattern5002 = Pattern(Integral((S(1)/sin(w_))**WC('n', S(1))*cos(v_), x_), cons89, cons90, cons1730)
    rule5002 = ReplacementRule(pattern5002, replacement5002)
    pattern5003 = Pattern(Integral((S(1)/sin(w_))**WC('n', S(1))*sin(v_), x_), cons89, cons90, cons1730)
    rule5003 = ReplacementRule(pattern5003, replacement5003)
    pattern5004 = Pattern(Integral((S(1)/cos(w_))**WC('n', S(1))*cos(v_), x_), cons89, cons90, cons1730)
    rule5004 = ReplacementRule(pattern5004, replacement5004)
    pattern5005 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5005 = ReplacementRule(pattern5005, replacement5005)
    pattern5006 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons1480, cons152, cons170, cons465, cons1731)
    rule5006 = ReplacementRule(pattern5006, replacement5006)
    pattern5007 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons1480, cons152, cons170, cons465, cons1731)
    rule5007 = ReplacementRule(pattern5007, replacement5007)
    pattern5008 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons13)
    rule5008 = ReplacementRule(pattern5008, replacement5008)
    pattern5009 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons64, cons13)
    rule5009 = ReplacementRule(pattern5009, replacement5009)
    pattern5010 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule5010 = ReplacementRule(pattern5010, replacement5010)
    pattern5011 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((b_ + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons64)
    rule5011 = ReplacementRule(pattern5011, replacement5011)
    pattern5012 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((WC('a', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(0)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule5012 = ReplacementRule(pattern5012, replacement5012)
    pattern5013 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((c_ + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons210, cons64)
    rule5013 = ReplacementRule(pattern5013, replacement5013)
    pattern5014 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/((WC('a', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(0)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons64, cons1480, cons1732)
    rule5014 = ReplacementRule(pattern5014, replacement5014)
    pattern5015 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1733)
    rule5015 = ReplacementRule(pattern5015, replacement5015)
    pattern5016 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1733)
    rule5016 = ReplacementRule(pattern5016, replacement5016)
    pattern5017 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1734)
    rule5017 = ReplacementRule(pattern5017, replacement5017)
    pattern5018 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1734)
    rule5018 = ReplacementRule(pattern5018, replacement5018)
    pattern5019 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1267)
    rule5019 = ReplacementRule(pattern5019, replacement5019)
    pattern5020 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1267)
    rule5020 = ReplacementRule(pattern5020, replacement5020)
    pattern5021 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1269)
    rule5021 = ReplacementRule(pattern5021, replacement5021)
    pattern5022 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons87, cons167, cons1269)
    rule5022 = ReplacementRule(pattern5022, replacement5022)
    pattern5023 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1476)
    rule5023 = ReplacementRule(pattern5023, replacement5023)
    pattern5024 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1476)
    rule5024 = ReplacementRule(pattern5024, replacement5024)
    pattern5025 = Pattern(Integral((a_ + WC('b', S(1))*tan(v_))**WC('n', S(1))*(S(1)/cos(v_))**WC('m', S(1)), x_), cons2, cons3, cons152, cons1553, cons1483)
    rule5025 = ReplacementRule(pattern5025, replacement5025)
    pattern5026 = Pattern(Integral((a_ + WC('b', S(1))/tan(v_))**WC('n', S(1))*(S(1)/sin(v_))**WC('m', S(1)), x_), cons2, cons3, cons152, cons1553, cons1483)
    rule5026 = ReplacementRule(pattern5026, replacement5026)
    pattern5027 = Pattern(Integral(WC('u', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons530)
    rule5027 = ReplacementRule(pattern5027, replacement5027)
    pattern5028 = Pattern(Integral(WC('u', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons530)
    rule5028 = ReplacementRule(pattern5028, replacement5028)
    pattern5029 = Pattern(Integral(S(1)/(cos(c_ + x_*WC('d', S(1)))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule5029 = ReplacementRule(pattern5029, replacement5029)
    pattern5030 = Pattern(Integral(S(1)/(sin(c_ + x_*WC('d', S(1)))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule5030 = ReplacementRule(pattern5030, replacement5030)
    pattern5031 = Pattern(Integral(tan(c_ + x_*WC('d', S(1)))*tan(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule5031 = ReplacementRule(pattern5031, replacement5031)
    pattern5032 = Pattern(Integral(S(1)/(tan(c_ + x_*WC('d', S(1)))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1735, cons73)
    rule5032 = ReplacementRule(pattern5032, replacement5032)
    pattern5033 = Pattern(Integral((WC('a', S(1))*cos(v_) + WC('b', S(1))*sin(v_))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1441)
    rule5033 = ReplacementRule(pattern5033, replacement5033)
    return [rule4688, rule4689, rule4690, rule4691, rule4692, rule4693, rule4694, rule4695, rule4696, rule4697, rule4698, rule4699, rule4700, rule4701, rule4702, rule4703, rule4704, rule4705, rule4706, rule4707, rule4708, rule4709, rule4710, rule4711, rule4712, rule4713, rule4714, rule4715, rule4716, rule4717, rule4718, rule4719, rule4720, rule4721, rule4722, rule4723, rule4724, rule4725, rule4726, rule4727, rule4728, rule4729, rule4730, rule4731, rule4732, rule4733, rule4734, rule4735, rule4736, rule4737, rule4738, rule4739, rule4740, rule4741, rule4742, rule4743, rule4744, rule4745, rule4746, rule4747, rule4748, rule4749, rule4750, rule4751, rule4752, rule4753, rule4754, rule4755, rule4756, rule4757, rule4758, rule4759, rule4760, rule4761, rule4762, rule4763, rule4764, rule4765, rule4766, rule4767, rule4768, rule4769, rule4770, rule4771, rule4772, rule4773, rule4774, rule4775, rule4776, rule4777, rule4778, rule4779, rule4780, rule4781, rule4782, rule4783, rule4784, rule4785, rule4786, rule4787, rule4788, rule4789, rule4790, rule4791, rule4792, rule4793, rule4794, rule4795, rule4796, rule4797, rule4798, rule4799, rule4800, rule4801, rule4802, rule4803, rule4804, rule4805, rule4806, rule4807, rule4808, rule4809, rule4810, rule4811, rule4812, rule4813, rule4814, rule4815, rule4816, rule4817, rule4818, rule4819, rule4820, rule4821, rule4822, rule4823, rule4824, rule4825, rule4826, rule4827, rule4828, rule4829, rule4830, rule4831, rule4832, rule4833, rule4834, rule4835, rule4836, rule4837, rule4838, rule4839, rule4840, rule4841, rule4842, rule4843, rule4844, rule4845, rule4846, rule4847, rule4848, rule4849, rule4850, rule4851, rule4852, rule4853, rule4854, rule4855, rule4856, rule4857, rule4858, rule4859, rule4860, rule4861, rule4862, rule4863, rule4864, rule4865, rule4866, rule4867, rule4868, rule4869, rule4870, rule4871, rule4872, rule4873, rule4874, rule4875, rule4876, rule4877, rule4878, rule4879, rule4880, rule4881, rule4882, rule4883, rule4884, rule4885, rule4886, rule4887, rule4888, rule4889, rule4890, rule4891, rule4892, rule4893, rule4894, rule4895, rule4896, rule4897, rule4898, rule4899, rule4900, rule4901, rule4902, rule4903, rule4904, rule4905, rule4906, rule4907, rule4908, rule4909, rule4910, rule4911, rule4912, rule4913, rule4914, rule4915, rule4916, rule4917, rule4918, rule4919, rule4920, rule4921, rule4922, rule4923, rule4924, rule4925, rule4926, rule4927, rule4928, rule4929, rule4930, rule4931, rule4932, rule4933, rule4934, rule4935, rule4936, rule4937, rule4938, rule4939, rule4940, rule4941, rule4942, rule4943, rule4944, rule4945, rule4946, rule4947, rule4948, rule4949, rule4950, rule4951, rule4952, rule4953, rule4954, rule4955, rule4956, rule4957, rule4958, rule4959, rule4960, rule4961, rule4962, rule4963, rule4964, rule4965, rule4966, rule4967, rule4968, rule4969, rule4970, rule4971, rule4972, rule4973, rule4974, rule4975, rule4976, rule4977, rule4978, rule4979, rule4980, rule4981, rule4982, rule4983, rule4984, rule4985, rule4986, rule4987, rule4988, rule4989, rule4990, rule4991, rule4992, rule4993, rule4994, rule4995, rule4996, rule4997, rule4998, rule4999, rule5000, rule5001, rule5002, rule5003, rule5004, rule5005, rule5006, rule5007, rule5008, rule5009, rule5010, rule5011, rule5012, rule5013, rule5014, rule5015, rule5016, rule5017, rule5018, rule5019, rule5020, rule5021, rule5022, rule5023, rule5024, rule5025, rule5026, rule5027, rule5028, rule5029, rule5030, rule5031, rule5032, rule5033, ]



def replacement4688(c, m, n, x, d, u, a, b):
        # rubi.append(4688)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**(m + n)*(d*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4689(c, m, n, x, d, u, a, b):
        # rubi.append(4689)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4690(c, m, n, x, d, u, a, b):
        # rubi.append(4690)
        return Dist((c/tan(a + b*x))**m*(d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m), Int((d*sin(a + b*x))**(-m + n)*(d*cos(a + b*x))**m*ActivateTrig(u), x), x)
def replacement4691(c, m, n, x, d, u, a, b):
        # rubi.append(4691)
        return Dist((c/tan(a + b*x))**m*(d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m), Int((d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**(m + n)*ActivateTrig(u), x), x)
def replacement4692(c, m, n, x, d, u, a, b):
        # rubi.append(4692)
        return Dist((c/sin(a + b*x))**m*(d*sin(a + b*x))**m, Int((d*sin(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4693(c, m, x, u, a, b):
        # rubi.append(4693)
        return Dist((c*sin(a + b*x))**(-m)*(c*cos(a + b*x))**m*(c*tan(a + b*x))**m, Int((c*sin(a + b*x))**m*(c*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4694(c, m, x, u, a, b):
        # rubi.append(4694)
        return Dist((c*sin(a + b*x))**m*(c*cos(a + b*x))**(-m)*(c/tan(a + b*x))**m, Int((c*sin(a + b*x))**(-m)*(c*cos(a + b*x))**m*ActivateTrig(u), x), x)
def replacement4695(c, m, x, u, a, b):
        # rubi.append(4695)
        return Dist((c/cos(a + b*x))**m*(c*cos(a + b*x))**m, Int((c*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4696(c, m, x, u, a, b):
        # rubi.append(4696)
        return Dist((c/sin(a + b*x))**m*(c*sin(a + b*x))**m, Int((c*sin(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4697(c, n, x, B, A, u, a, b):
        # rubi.append(4697)
        return Dist(c, Int((c*sin(a + b*x))**(n + S(-1))*(A*sin(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4698(c, n, x, B, A, u, a, b):
        # rubi.append(4698)
        return Dist(c, Int((c*cos(a + b*x))**(n + S(-1))*(A*cos(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4699(x, B, A, u, a, b):
        # rubi.append(4699)
        return Int((A*sin(a + b*x) + B)*ActivateTrig(u)/sin(a + b*x), x)
def replacement4700(x, B, A, u, a, b):
        # rubi.append(4700)
        return Int((A*cos(a + b*x) + B)*ActivateTrig(u)/cos(a + b*x), x)
def replacement4701(c, n, x, B, A, u, a, C, b):
        # rubi.append(4701)
        return Dist(c**S(2), Int((c*sin(a + b*x))**(n + S(-2))*(A*sin(a + b*x)**S(2) + B*sin(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4702(c, n, x, B, A, u, a, C, b):
        # rubi.append(4702)
        return Dist(c**S(2), Int((c*cos(a + b*x))**(n + S(-2))*(A*cos(a + b*x)**S(2) + B*cos(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4703(c, n, x, A, u, a, C, b):
        # rubi.append(4703)
        return Dist(c**S(2), Int((c*sin(a + b*x))**(n + S(-2))*(A*sin(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4704(c, n, x, A, u, a, C, b):
        # rubi.append(4704)
        return Dist(c**S(2), Int((c*cos(a + b*x))**(n + S(-2))*(A*cos(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4705(x, B, A, u, a, C, b):
        # rubi.append(4705)
        return Int((A*sin(a + b*x)**S(2) + B*sin(a + b*x) + C)*ActivateTrig(u)/sin(a + b*x)**S(2), x)
def replacement4706(x, B, A, u, a, C, b):
        # rubi.append(4706)
        return Int((A*cos(a + b*x)**S(2) + B*cos(a + b*x) + C)*ActivateTrig(u)/cos(a + b*x)**S(2), x)
def replacement4707(x, A, u, a, C, b):
        # rubi.append(4707)
        return Int((A*sin(a + b*x)**S(2) + C)*ActivateTrig(u)/sin(a + b*x)**S(2), x)
def replacement4708(x, A, u, a, C, b):
        # rubi.append(4708)
        return Int((A*cos(a + b*x)**S(2) + C)*ActivateTrig(u)/cos(a + b*x)**S(2), x)
def replacement4709(x, B, A, u, a, C, b):
        # rubi.append(4709)
        return Int((A*sin(a + b*x) + B*sin(a + b*x)**S(2) + C)*ActivateTrig(u)/sin(a + b*x), x)
def replacement4710(x, B, A, u, a, C, b):
        # rubi.append(4710)
        return Int((A*cos(a + b*x) + B*cos(a + b*x)**S(2) + C)*ActivateTrig(u)/cos(a + b*x), x)
def replacement4711(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4711)
        return Int((A + B*sin(a + b*x) + C*sin(a + b*x)**S(2))*ActivateTrig(u)*sin(a + b*x)**n, x)
def replacement4712(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4712)
        return Int((A + B*cos(a + b*x) + C*cos(a + b*x)**S(2))*ActivateTrig(u)*cos(a + b*x)**n, x)
def replacement4713(c, m, n, x, d, u, a, b):
        # rubi.append(4713)
        return Dist((c/tan(a + b*x))**m*(d*tan(a + b*x))**m, Int((d*tan(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4714(c, m, n, x, d, u, a, b):
        # rubi.append(4714)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4715(c, m, x, u, a, b):
        # rubi.append(4715)
        return Dist((c/tan(a + b*x))**m*(c*tan(a + b*x))**m, Int((c*tan(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4716(c, m, x, u, a, b):
        # rubi.append(4716)
        return Dist((c/tan(a + b*x))**m*(c*tan(a + b*x))**m, Int((c/tan(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4717(c, n, x, B, A, u, a, b):
        # rubi.append(4717)
        return Dist(c, Int((c*tan(a + b*x))**(n + S(-1))*(A*tan(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4718(c, n, x, B, A, u, a, b):
        # rubi.append(4718)
        return Dist(c, Int((c/tan(a + b*x))**(n + S(-1))*(A/tan(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4719(x, B, A, u, a, b):
        # rubi.append(4719)
        return Int((A*tan(a + b*x) + B)*ActivateTrig(u)/tan(a + b*x), x)
def replacement4720(x, B, A, u, a, b):
        # rubi.append(4720)
        return Int((A/tan(a + b*x) + B)*ActivateTrig(u)*tan(a + b*x), x)
def replacement4721(c, n, x, B, A, u, a, C, b):
        # rubi.append(4721)
        return Dist(c**S(2), Int((c*tan(a + b*x))**(n + S(-2))*(A*tan(a + b*x)**S(2) + B*tan(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4722(c, n, x, B, A, u, a, C, b):
        # rubi.append(4722)
        return Dist(c**S(2), Int((c/tan(a + b*x))**(n + S(-2))*(A/tan(a + b*x)**S(2) + B/tan(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4723(c, n, x, A, u, a, C, b):
        # rubi.append(4723)
        return Dist(c**S(2), Int((c*tan(a + b*x))**(n + S(-2))*(A*tan(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4724(c, n, x, A, u, a, C, b):
        # rubi.append(4724)
        return Dist(c**S(2), Int((c/tan(a + b*x))**(n + S(-2))*(A/tan(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4725(x, B, A, u, a, C, b):
        # rubi.append(4725)
        return Int((A*tan(a + b*x)**S(2) + B*tan(a + b*x) + C)*ActivateTrig(u)/tan(a + b*x)**S(2), x)
def replacement4726(x, B, A, u, a, C, b):
        # rubi.append(4726)
        return Int((A/tan(a + b*x)**S(2) + B/tan(a + b*x) + C)*ActivateTrig(u)*tan(a + b*x)**S(2), x)
def replacement4727(x, A, u, a, C, b):
        # rubi.append(4727)
        return Int((A*tan(a + b*x)**S(2) + C)*ActivateTrig(u)/tan(a + b*x)**S(2), x)
def replacement4728(x, A, u, a, C, b):
        # rubi.append(4728)
        return Int((A/tan(a + b*x)**S(2) + C)*ActivateTrig(u)*tan(a + b*x)**S(2), x)
def replacement4729(x, B, A, u, a, C, b):
        # rubi.append(4729)
        return Int((A*tan(a + b*x) + B*tan(a + b*x)**S(2) + C)*ActivateTrig(u)/tan(a + b*x), x)
def replacement4730(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4730)
        return Int((A + B*tan(a + b*x) + C*tan(a + b*x)**S(2))*ActivateTrig(u)*tan(a + b*x)**n, x)
def replacement4731(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4731)
        return Int((A + B/tan(a + b*x) + C/tan(a + b*x)**S(2))*(S(1)/tan(a + b*x))**n*ActivateTrig(u), x)
def replacement4732(c, m, n, x, d, u, a, b):
        # rubi.append(4732)
        return Dist((c*sin(a + b*x))**m*(d/sin(a + b*x))**m, Int((d/sin(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4733(c, m, n, x, d, u, a, b):
        # rubi.append(4733)
        return Dist((c*cos(a + b*x))**m*(d/cos(a + b*x))**m, Int((d/cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4734(c, m, n, x, d, u, a, b):
        # rubi.append(4734)
        return Dist((c*tan(a + b*x))**m*(d/sin(a + b*x))**m*(d/cos(a + b*x))**(-m), Int((d/sin(a + b*x))**(-m)*(d/cos(a + b*x))**(m + n)*ActivateTrig(u), x), x)
def replacement4735(c, m, n, x, d, u, a, b):
        # rubi.append(4735)
        return Dist((c*tan(a + b*x))**m*(d/sin(a + b*x))**m*(d/cos(a + b*x))**(-m), Int((d/sin(a + b*x))**(-m + n)*(d/cos(a + b*x))**m*ActivateTrig(u), x), x)
def replacement4736(c, m, n, x, d, u, a, b):
        # rubi.append(4736)
        return Dist((c/tan(a + b*x))**m*(d/sin(a + b*x))**(-m)*(d/cos(a + b*x))**m, Int((d/sin(a + b*x))**m*(d/cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
def replacement4737(c, m, n, x, d, u, a, b):
        # rubi.append(4737)
        return Dist((c/tan(a + b*x))**m*(d/sin(a + b*x))**(-m)*(d/cos(a + b*x))**m, Int((d/sin(a + b*x))**(m + n)*(d/cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4738(c, m, x, u, a, b):
        # rubi.append(4738)
        return Dist((c/sin(a + b*x))**m*(c*sin(a + b*x))**m, Int((c/sin(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4739(c, m, x, u, a, b):
        # rubi.append(4739)
        return Dist((c/cos(a + b*x))**m*(c*cos(a + b*x))**m, Int((c/cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4740(c, m, x, u, a, b):
        # rubi.append(4740)
        return Dist((c/sin(a + b*x))**m*(c/cos(a + b*x))**(-m)*(c*tan(a + b*x))**m, Int((c/sin(a + b*x))**(-m)*(c/cos(a + b*x))**m*ActivateTrig(u), x), x)
def replacement4741(c, m, x, u, a, b):
        # rubi.append(4741)
        return Dist((c/sin(a + b*x))**(-m)*(c/cos(a + b*x))**m*(c/tan(a + b*x))**m, Int((c/sin(a + b*x))**m*(c/cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
def replacement4742(c, n, x, B, A, u, a, b):
        # rubi.append(4742)
        return Dist(c, Int((c/cos(a + b*x))**(n + S(-1))*(A/cos(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4743(c, n, x, B, A, u, a, b):
        # rubi.append(4743)
        return Dist(c, Int((c/sin(a + b*x))**(n + S(-1))*(A/sin(a + b*x) + B)*ActivateTrig(u), x), x)
def replacement4744(x, B, A, u, a, b):
        # rubi.append(4744)
        return Int((A/cos(a + b*x) + B)*ActivateTrig(u)*cos(a + b*x), x)
def replacement4745(x, B, A, u, a, b):
        # rubi.append(4745)
        return Int((A/sin(a + b*x) + B)*ActivateTrig(u)*sin(a + b*x), x)
def replacement4746(c, n, x, B, A, u, a, C, b):
        # rubi.append(4746)
        return Dist(c**S(2), Int((c/cos(a + b*x))**(n + S(-2))*(A/cos(a + b*x)**S(2) + B/cos(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4747(c, n, x, B, A, u, a, C, b):
        # rubi.append(4747)
        return Dist(c**S(2), Int((c/sin(a + b*x))**(n + S(-2))*(A/sin(a + b*x)**S(2) + B/sin(a + b*x) + C)*ActivateTrig(u), x), x)
def replacement4748(c, n, x, A, u, a, C, b):
        # rubi.append(4748)
        return Dist(c**S(2), Int((c/cos(a + b*x))**(n + S(-2))*(A/cos(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4749(c, n, x, A, u, a, C, b):
        # rubi.append(4749)
        return Dist(c**S(2), Int((c/sin(a + b*x))**(n + S(-2))*(A/sin(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
def replacement4750(x, B, A, u, a, C, b):
        # rubi.append(4750)
        return Int((A/cos(a + b*x)**S(2) + B/cos(a + b*x) + C)*ActivateTrig(u)*cos(a + b*x)**S(2), x)
def replacement4751(x, B, A, u, a, C, b):
        # rubi.append(4751)
        return Int((A/sin(a + b*x)**S(2) + B/sin(a + b*x) + C)*ActivateTrig(u)*sin(a + b*x)**S(2), x)
def replacement4752(x, A, u, a, C, b):
        # rubi.append(4752)
        return Int((A/cos(a + b*x)**S(2) + C)*ActivateTrig(u)*cos(a + b*x)**S(2), x)
def replacement4753(x, A, u, a, C, b):
        # rubi.append(4753)
        return Int((A/sin(a + b*x)**S(2) + C)*ActivateTrig(u)*sin(a + b*x)**S(2), x)
def replacement4754(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4754)
        return Int((A + B/cos(a + b*x) + C/cos(a + b*x)**S(2))*(S(1)/cos(a + b*x))**n*ActivateTrig(u), x)
def replacement4755(n, x, n1, B, A, u, n2, a, C, b):
        # rubi.append(4755)
        return Int((A + B/sin(a + b*x) + C/sin(a + b*x)**S(2))*(S(1)/sin(a + b*x))**n*ActivateTrig(u), x)
def replacement4756(c, x, d, a, b):
        # rubi.append(4756)
        return Simp(sin(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) - Simp(sin(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
def replacement4757(c, x, d, a, b):
        # rubi.append(4757)
        return Simp(sin(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) + Simp(sin(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
def replacement4758(c, x, d, a, b):
        # rubi.append(4758)
        return -Simp(cos(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) - Simp(cos(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
def replacement4759(c, x, d, g, a, p, b):
        # rubi.append(4759)
        return Dist(S(1)/2, Int((g*sin(c + d*x))**p, x), x) + Dist(S(1)/2, Int((g*sin(c + d*x))**p*cos(c + d*x), x), x)
def replacement4760(c, x, d, g, a, p, b):
        # rubi.append(4760)
        return Dist(S(1)/2, Int((g*sin(c + d*x))**p, x), x) - Dist(S(1)/2, Int((g*sin(c + d*x))**p*cos(c + d*x), x), x)
def replacement4761(c, m, x, d, a, p, e, b):
        # rubi.append(4761)
        return Dist(S(2)**p*e**(-p), Int((e*cos(a + b*x))**(m + p)*sin(a + b*x)**p, x), x)
def replacement4762(c, n, x, d, f, a, p, b):
        # rubi.append(4762)
        return Dist(S(2)**p*f**(-p), Int((f*sin(a + b*x))**(n + p)*cos(a + b*x)**p, x), x)
def replacement4763(c, m, x, d, g, a, p, e, b):
        # rubi.append(4763)
        return Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4764(c, m, x, d, g, a, p, e, b):
        # rubi.append(4764)
        return -Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4765(c, m, x, d, g, a, p, e, b):
        # rubi.append(4765)
        return -Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(b*g*m), x)
def replacement4766(c, m, x, d, g, a, p, e, b):
        # rubi.append(4766)
        return Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(b*g*m), x)
def replacement4767(c, m, x, d, g, a, p, e, b):
        # rubi.append(4767)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(p + S(1))), Int((e*cos(a + b*x))**(m + S(-4))*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4768(c, m, x, d, g, a, p, e, b):
        # rubi.append(4768)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(p + S(1))), Int((e*sin(a + b*x))**(m + S(-4))*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4769(c, m, x, d, g, a, p, e, b):
        # rubi.append(4769)
        return Dist(e**S(2)*(m + S(2)*p + S(2))/(S(4)*g**S(2)*(p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4770(c, m, x, d, g, a, p, e, b):
        # rubi.append(4770)
        return Dist(e**S(2)*(m + S(2)*p + S(2))/(S(4)*g**S(2)*(p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
def replacement4771(c, m, x, d, g, a, p, e, b):
        # rubi.append(4771)
        return Dist(e**S(2)*(m + p + S(-1))/(m + S(2)*p), Int((e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**p, x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)), x)
def replacement4772(c, m, x, d, g, a, p, e, b):
        # rubi.append(4772)
        return Dist(e**S(2)*(m + p + S(-1))/(m + S(2)*p), Int((e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**p, x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)), x)
def replacement4773(c, m, x, d, g, a, p, e, b):
        # rubi.append(4773)
        return Dist((m + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(2))*(g*sin(c + d*x))**p, x), x) - Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))), x)
def replacement4774(c, m, x, d, g, a, p, e, b):
        # rubi.append(4774)
        return Dist((m + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(2))*(g*sin(c + d*x))**p, x), x) + Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))), x)
def replacement4775(c, x, d, g, a, p, b):
        # rubi.append(4775)
        return Dist(S(2)*g*p/(S(2)*p + S(1)), Int((g*sin(c + d*x))**(p + S(-1))*sin(a + b*x), x), x) + Simp(S(2)*(g*sin(c + d*x))**p*sin(a + b*x)/(d*(S(2)*p + S(1))), x)
def replacement4776(c, x, d, g, a, p, b):
        # rubi.append(4776)
        return Dist(S(2)*g*p/(S(2)*p + S(1)), Int((g*sin(c + d*x))**(p + S(-1))*cos(a + b*x), x), x) + Simp(-S(2)*(g*sin(c + d*x))**p*cos(a + b*x)/(d*(S(2)*p + S(1))), x)
def replacement4777(c, x, d, g, a, p, b):
        # rubi.append(4777)
        return Dist((S(2)*p + S(3))/(S(2)*g*(p + S(1))), Int((g*sin(c + d*x))**(p + S(1))*sin(a + b*x), x), x) + Simp((g*sin(c + d*x))**(p + S(1))*cos(a + b*x)/(S(2)*b*g*(p + S(1))), x)
def replacement4778(c, x, d, g, a, p, b):
        # rubi.append(4778)
        return Dist((S(2)*p + S(3))/(S(2)*g*(p + S(1))), Int((g*sin(c + d*x))**(p + S(1))*cos(a + b*x), x), x) - Simp((g*sin(c + d*x))**(p + S(1))*sin(a + b*x)/(S(2)*b*g*(p + S(1))), x)
def replacement4779(c, x, d, a, b):
        # rubi.append(4779)
        return Simp(log(sin(a + b*x) + sqrt(sin(c + d*x)) + cos(a + b*x))/d, x) - Simp(-asin(sin(a + b*x) - cos(a + b*x))/d, x)
def replacement4780(c, x, d, a, b):
        # rubi.append(4780)
        return -Simp(log(sin(a + b*x) + sqrt(sin(c + d*x)) + cos(a + b*x))/d, x) - Simp(-asin(sin(a + b*x) - cos(a + b*x))/d, x)
def replacement4781(c, x, d, g, a, p, b):
        # rubi.append(4781)
        return Dist(S(2)*g, Int((g*sin(c + d*x))**(p + S(-1))*sin(a + b*x), x), x)
def replacement4782(c, x, d, g, a, p, b):
        # rubi.append(4782)
        return Dist(S(2)*g, Int((g*sin(c + d*x))**(p + S(-1))*cos(a + b*x), x), x)
def replacement4783(c, m, x, d, g, a, p, e, b):
        # rubi.append(4783)
        return Dist((e*cos(a + b*x))**(-p)*(g*sin(c + d*x))**p*sin(a + b*x)**(-p), Int((e*cos(a + b*x))**(m + p)*sin(a + b*x)**p, x), x)
def replacement4784(c, n, x, d, f, g, a, p, b):
        # rubi.append(4784)
        return Dist((f*sin(a + b*x))**(-p)*(g*sin(c + d*x))**p*cos(a + b*x)**(-p), Int((f*sin(a + b*x))**(n + p)*cos(a + b*x)**p, x), x)
def replacement4785(c, x, d, g, a, p, b):
        # rubi.append(4785)
        return Dist(S(1)/4, Int((g*sin(c + d*x))**p, x), x) - Dist(S(1)/4, Int((g*sin(c + d*x))**p*cos(c + d*x)**S(2), x), x)
def replacement4786(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(4786)
        return Dist(S(2)**p*e**(-p)*f**(-p), Int((e*cos(a + b*x))**(m + p)*(f*sin(a + b*x))**(n + p), x), x)
def replacement4787(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4787)
        return Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
def replacement4788(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4788)
        return -Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
def replacement4789(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4789)
        return -Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
def replacement4790(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4790)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-4))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
def replacement4791(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4791)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-4))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
def replacement4792(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4792)
        return Dist(e**S(2)*(m + n + S(2)*p + S(2))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp((e*cos(a + b*x))**m*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
def replacement4793(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4793)
        return Dist(e**S(2)*(m + n + S(2)*p + S(2))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp((e*sin(a + b*x))**m*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
def replacement4794(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4794)
        return Dist(e**S(2)*(m + p + S(-1))/(f**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**(n + S(2))*(g*sin(c + d*x))**p, x), x) + Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
def replacement4795(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4795)
        return Dist(e**S(2)*(m + p + S(-1))/(f**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**(n + S(2))*(g*sin(c + d*x))**p, x), x) - Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
def replacement4796(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4796)
        return Dist(e**S(2)*(m + p + S(-1))/(m + n + S(2)*p), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**p, x), x) + Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(m + n + S(2)*p)), x)
def replacement4797(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4797)
        return Dist(e**S(2)*(m + p + S(-1))/(m + n + S(2)*p), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**p, x), x) - Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(m + n + S(2)*p)), x)
def replacement4798(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4798)
        return Dist(S(2)*f*g*(n + p + S(-1))/(e*(m + n + S(2)*p)), Int((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(-1)), x), x) - Simp(f*(e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**p/(b*e*(m + n + S(2)*p)), x)
def replacement4799(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4799)
        return Dist(S(2)*f*g*(n + p + S(-1))/(e*(m + n + S(2)*p)), Int((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(-1)), x), x) + Simp(f*(e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**p/(b*e*(m + n + S(2)*p)), x)
def replacement4800(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4800)
        return Dist(f*(m + n + S(2)*p + S(2))/(S(2)*e*g*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(1)), x), x) - Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
def replacement4801(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4801)
        return Dist(f*(m + n + S(2)*p + S(2))/(S(2)*e*g*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(1)), x), x) + Simp((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
def replacement4802(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4802)
        return Dist((m + n + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**p, x), x) - Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
def replacement4803(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4803)
        return Dist((m + n + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**p, x), x) + Simp((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
def replacement4804(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(4804)
        return Dist((e*cos(a + b*x))**(-p)*(f*sin(a + b*x))**(-p)*(g*sin(c + d*x))**p, Int((e*cos(a + b*x))**(m + p)*(f*sin(a + b*x))**(n + p), x), x)
def replacement4805(c, m, x, d, a, e, b):
        # rubi.append(4805)
        return -Simp((e*cos(a + b*x))**(m + S(1))*(m + S(2))*cos((a + b*x)*(m + S(1)))/(d*e*(m + S(1))), x)
def replacement4806(c, n, F, x, d, a, p, b):
        # rubi.append(4806)
        return Int((a + b*F(c + d*x)**n)**p, x)
def replacement4807(c, n, F, x, d, a, b):
        # rubi.append(4807)
        return Dist(S(2)/(a*n), Sum_doit(Int(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*F(c + d*x)**S(2)/Rt(-a/b, n/S(2))), x), List(k, S(1), n/S(2))), x)
def replacement4808(c, n, F, x, d, a, b):
        # rubi.append(4808)
        return Int(ExpandTrig(S(1)/(a + b*F(c + d*x)**n), x), x)
def replacement4809(c, G, m, n, F, x, d, a, b):
        # rubi.append(4809)
        return Int(ExpandTrig(G(c + d*x)**m, S(1)/(a + b*F(c + d*x)**n), x), x)

def With4810(c, n, F, x, d, a, p):
        v = ActivateTrig(F(c + d*x))
        # rubi.append(4810)
        return Dist(a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*NonfreeFactors(v, x)**(-p*FracPart(n)), Int(NonfreeFactors(v, x)**(n*p), x), x)

def With4811(c, n, F, x, d, a, p, b):
        v = ActivateTrig(F(c + d*x))
        # rubi.append(4811)
        return Dist(a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n)), Int((b*v)**(n*p), x), x)

def With4812(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4812(c, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4812)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4813(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4813(c, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4813)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4814(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4814(c, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4814)
        return Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4815(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4815(c, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4815)
        return -Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4816(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if FunctionOfQ(tan(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4816(c, F, x, u, a, b):

        d = FreeFactors(tan(c*(a + b*x)), x)
        # rubi.append(4816)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)

def With4817(c, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if FunctionOfQ(tan(c*(a + b*x))/d, u, x, True):
            return True
        return False
def replacement4817(c, x, u, a, b):

        d = FreeFactors(tan(c*(a + b*x)), x)
        # rubi.append(4817)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)

def With4818(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        if FunctionOfQ(S(1)/(d*tan(c*(a + b*x))), u, x, True):
            return True
        return False
def replacement4818(c, F, x, u, a, b):

        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        # rubi.append(4818)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), S(1)/(d*tan(c*(a + b*x))), u, x), x), x, S(1)/(d*tan(c*(a + b*x)))), x)

def With4819(c, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        if FunctionOfQ(S(1)/(d*tan(c*(a + b*x))), u, x, True):
            return True
        return False
def replacement4819(c, x, u, a, b):

        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        # rubi.append(4819)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), S(1)/(d*tan(c*(a + b*x))), u, x), x), x, S(1)/(d*tan(c*(a + b*x)))), x)

def With4820(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if And(FunctionOfQ(tan(c*(a + b*x))/d, u, x, True), TryPureTanSubst((S(1)/tan(c*(a + b*x)))**n*ActivateTrig(u), x)):
            return True
        return False
def replacement4820(c, n, F, x, u, a, b):

        d = FreeFactors(tan(c*(a + b*x)), x)
        # rubi.append(4820)
        return Dist(d**(S(1) - n)/(b*c), Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)

def With4821(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        if And(FunctionOfQ(S(1)/(d*tan(c*(a + b*x))), u, x, True), TryPureTanSubst(ActivateTrig(u)*tan(c*(a + b*x))**n, x)):
            return True
        return False
def replacement4821(c, n, F, x, u, a, b):

        d = FreeFactors(S(1)/tan(c*(a + b*x)), x)
        # rubi.append(4821)
        return -Dist(d**(S(1) - n)/(b*c), Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), S(1)/(d*tan(c*(a + b*x))), u, x), x), x, S(1)/(d*tan(c*(a + b*x)))), x)

def With4822(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            v = FunctionOfTrig(u, x)
            d = FreeFactors(S(1)/tan(v), x)
            res = And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(S(1)/tan(v), x), u, x, True), TryPureTanSubst(ActivateTrig(u), x))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4822(x, u):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(S(1)/tan(v), x)
        # rubi.append(4822)
        return Simp(With(List(Set(d, FreeFactors(S(1)/tan(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1)/(d**S(2)*x**S(2) + S(1)), S(1)/(d*tan(v)), u, x), x), x, S(1)/(d*tan(v))), x)), x)

def With4823(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            v = FunctionOfTrig(u, x)
            d = FreeFactors(tan(v), x)
            res = And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(tan(v), x), u, x, True), TryPureTanSubst(ActivateTrig(u), x))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4823(x, u):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        # rubi.append(4823)
        return Simp(With(List(Set(d, FreeFactors(tan(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1)/(d**S(2)*x**S(2) + S(1)), tan(v)/d, u, x), x), x, tan(v)/d), x)), x)
def replacement4824(c, G, q, F, x, d, a, p, b):
        # rubi.append(4824)
        return Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q), x), x)
def replacement4825(c, G, q, F, r, x, d, H, f, a, p, e, b):
        # rubi.append(4825)
        return Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q*H(e + f*x)**r), x), x)

def With4826(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4826(c, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4826)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4827(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4827(c, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4827)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4828(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4828(c, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4828)
        return Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4829(c, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4829(c, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4829)
        return -Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4830(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4830(c, n, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4830)
        return Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4831(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4831(c, n, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4831)
        return Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4832(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4832(c, n, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4832)
        return -Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4833(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4833(c, n, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4833)
        return -Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4834(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4834(c, n, F, x, u, a, b):

        d = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4834)
        return Dist(d**(S(1) - n)/(b*c), Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)

def With4835(c, n, F, x, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
def replacement4835(c, n, F, x, u, a, b):

        d = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4835)
        return -Dist(d**(S(1) - n)/(b*c), Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)

def With4836(c, n, F, x, d, v, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        e = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/e, u, x):
            return True
        return False
def replacement4836(c, n, F, x, d, v, u, a, b):

        e = FreeFactors(sin(c*(a + b*x)), x)
        # rubi.append(4836)
        return Dist(d, Int(ActivateTrig(u)*cos(c*(a + b*x))**n, x), x) + Int(ActivateTrig(u*v), x)

def With4837(c, n, F, x, d, v, u, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        e = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/e, u, x):
            return True
        return False
def replacement4837(c, n, F, x, d, v, u, a, b):

        e = FreeFactors(cos(c*(a + b*x)), x)
        # rubi.append(4837)
        return Dist(d, Int(ActivateTrig(u)*sin(c*(a + b*x))**n, x), x) + Int(ActivateTrig(u*v), x)

def With4838(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            v = FunctionOfTrig(u, x)
            d = FreeFactors(sin(v), x)
            res = And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(sin(v), x), u/cos(v), x))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4838(x, u):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(sin(v), x)
        # rubi.append(4838)
        return Simp(With(List(Set(d, FreeFactors(sin(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), sin(v)/d, u/cos(v), x), x), x, sin(v)/d), x)), x)

def With4839(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            v = FunctionOfTrig(u, x)
            d = FreeFactors(cos(v), x)
            res = And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(cos(v), x), u/sin(v), x))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4839(x, u):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(cos(v), x)
        # rubi.append(4839)
        return Simp(With(List(Set(d, FreeFactors(cos(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), cos(v)/d, u/sin(v), x), x), x, cos(v)/d), x)), x)
def replacement4840(c, x, d, u, a, p, e, b):
        # rubi.append(4840)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)
def replacement4841(c, x, d, u, a, p, e, b):
        # rubi.append(4841)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)
def replacement4842(c, x, d, u, a, p, e, b):
        # rubi.append(4842)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)

def With4843(x, u, y):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
            res = Not(FalseQ(q))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4843(x, u, y):

        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        # rubi.append(4843)
        return Simp(q*log(RemoveContent(ActivateTrig(y), x)), x)

def With4844(x, u, w, y):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            q = DerivativeDivides(ActivateTrig(w*y), ActivateTrig(u), x)
            res = Not(FalseQ(q))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4844(x, u, w, y):

        q = DerivativeDivides(ActivateTrig(w*y), ActivateTrig(u), x)
        # rubi.append(4844)
        return Simp(q*log(RemoveContent(ActivateTrig(w*y), x)), x)

def With4845(x, u, m, y):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
            res = Not(FalseQ(q))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4845(x, u, m, y):

        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        # rubi.append(4845)
        return Simp(q*ActivateTrig(y**(m + S(1)))/(m + S(1)), x)

def With4846(m, n, x, u, y, z):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            q = DerivativeDivides(ActivateTrig(y*z), ActivateTrig(u*z**(-m + n)), x)
            res = Not(FalseQ(q))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4846(m, n, x, u, y, z):

        q = DerivativeDivides(ActivateTrig(y*z), ActivateTrig(u*z**(-m + n)), x)
        # rubi.append(4846)
        return Simp(q*ActivateTrig(y**(m + S(1))*z**(m + S(1)))/(m + S(1)), x)

def With4847(c, n, F, x, d, u, a, p):
        v = ActivateTrig(F(c + d*x))
        # rubi.append(4847)
        return Dist(a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*NonfreeFactors(v, x)**(-p*FracPart(n)), Int(ActivateTrig(u)*NonfreeFactors(v, x)**(n*p), x), x)

def With4848(c, n, F, x, d, u, a, p, b):
        v = ActivateTrig(F(c + d*x))
        # rubi.append(4848)
        return Dist(a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n)), Int((b*v)**(n*p)*ActivateTrig(u), x), x)

def With4849(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            v = FunctionOfTrig(u, x)
            d = FreeFactors(tan(v), x)
            res = And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(tan(v), x), u, x))
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4849(x, u):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        # rubi.append(4849)
        return Dist(d/Coefficient(v, x, 1), Subst(Int(SubstFor(1/(d**2*x**2 + 1), tan(v)/d, u, x), x), x, tan(v)/d), x)
def replacement4850(c, n, x, d, u, a, p, b):
        # rubi.append(4850)
        return Int((a*sin(c + d*x)**n + b)**p*(S(1)/cos(c + d*x))**(n*p)*ActivateTrig(u), x)
def replacement4851(c, n, x, d, u, a, p, b):
        # rubi.append(4851)
        return Int((a*cos(c + d*x)**n + b)**p*(S(1)/sin(c + d*x))**(n*p)*ActivateTrig(u), x)
def replacement4852(c, q, n, F, x, d, u, a, p, b):
        # rubi.append(4852)
        return Int(ActivateTrig(u*(a + b*F(c + d*x)**(-p + q))**n*F(c + d*x)**(n*p)), x)
def replacement4853(c, q, n, F, r, x, d, u, a, p, e, b):
        # rubi.append(4853)
        return Int(ActivateTrig(u*(a + b*F(d + e*x)**(-p + q) + c*F(d + e*x)**(-p + r))**n*F(d + e*x)**(n*p)), x)
def replacement4854(c, q, n, F, x, d, u, a, p, e, b):
        # rubi.append(4854)
        return Int(ActivateTrig(u*(a*F(d + e*x)**(-p) + b + c*F(d + e*x)**(-p + q))**n*F(d + e*x)**(n*p)), x)
def replacement4855(c, n, x, d, u, a, b):
        # rubi.append(4855)
        return Int((a*exp(-a*(c + d*x)/b))**n*ActivateTrig(u), x)
def replacement4856(x, u):
        # rubi.append(4856)
        return Int(TrigSimplify(u), x)

def With4857(x, v, u, a, p):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        # rubi.append(4857)
        return Dist(a**IntPart(p)*vv**(-FracPart(p))*(a*vv)**FracPart(p), Int(uu*vv**p, x), x)

def With4858(m, x, v, u, p):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        # rubi.append(4858)
        return Dist(vv**(-m*FracPart(p))*(vv**m)**FracPart(p), Int(uu*vv**(m*p), x), x)

def With4859(m, n, x, v, u, w, p):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        ww = ActivateTrig(w)
        # rubi.append(4859)
        return Dist(vv**(-m*FracPart(p))*ww**(-n*FracPart(p))*(vv**m*ww**n)**FracPart(p), Int(uu*vv**(m*p)*ww**(n*p), x), x)

def With4860(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = ExpandTrig(u, x)
        if SumQ(v):
            return True
        return False
def replacement4860(x, u):

        v = ExpandTrig(u, x)
        # rubi.append(4860)
        return Int(v, x)

def With4861(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        try:
            w = With(List(Set(ShowSteps, False), Set(StepCounter, Null)), Int(SubstFor(S(1)/(x**S(2)*FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x)**S(2) + S(1)), tan(FunctionOfTrig(u, x)/S(2))/FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x), u, x), x))
            v = FunctionOfTrig(u, x)
            d = FreeFactors(tan(v/S(2)), x)
            res = FreeQ(w, Int)
        except (TypeError, AttributeError):
            return False
        if res:
            return True
        return False
def replacement4861(x, u):

        w = With(List(Set(ShowSteps, False), Set(StepCounter, Null)), Int(SubstFor(S(1)/(x**S(2)*FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x)**S(2) + S(1)), tan(FunctionOfTrig(u, x)/S(2))/FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x), u, x), x))
        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v/S(2)), x)
        # rubi.append(4861)
        return Simp(Dist(2*d/Coefficient(v, x, 1), Subst(Int(SubstFor(1/(d**2*x**2 + 1), tan(v/2)/d, u, x), x), x, tan(v/2)/d), x), x)

def With4862(x, u):
        v = ActivateTrig(u)
        # rubi.append(4862)
        return Int(v, x)
def replacement4863(c, m, n, x, d, a, b):
        # rubi.append(4863)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*sin(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*sin(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement4864(c, m, n, x, d, a, b):
        # rubi.append(4864)
        return Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*cos(a + b*x)**(n + S(1)), x), x) - Simp((c + d*x)**m*cos(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement4865(c, m, n, x, d, a, p, b):
        # rubi.append(4865)
        return Int(ExpandTrigReduce((c + d*x)**m, sin(a + b*x)**n*cos(a + b*x)**p, x), x)
def replacement4866(c, m, n, x, d, a, p, b):
        # rubi.append(4866)
        return -Int((c + d*x)**m*sin(a + b*x)**n*tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*sin(a + b*x)**(n + S(-2))*tan(a + b*x)**p, x)
def replacement4867(c, m, n, x, d, a, p, b):
        # rubi.append(4867)
        return Int((c + d*x)**m*(S(1)/tan(a + b*x))**p*cos(a + b*x)**(n + S(-2)), x) - Int((c + d*x)**m*(S(1)/tan(a + b*x))**(p + S(-2))*cos(a + b*x)**n, x)
def replacement4868(c, m, n, x, d, a, p, b):
        # rubi.append(4868)
        return -Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*(S(1)/cos(a + b*x))**n, x), x) + Simp((c + d*x)**m*(S(1)/cos(a + b*x))**n/(b*n), x)
def replacement4869(c, m, n, x, d, a, p, b):
        # rubi.append(4869)
        return Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*(S(1)/sin(a + b*x))**n, x), x) - Simp((c + d*x)**m*(S(1)/sin(a + b*x))**n/(b*n), x)
def replacement4870(c, m, n, x, d, a, b):
        # rubi.append(4870)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*tan(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*tan(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
def replacement4871(c, m, n, x, d, a, b):
        # rubi.append(4871)
        return Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*(S(1)/tan(a + b*x))**(n + S(1)), x), x) - Simp((c + d*x)**m*(S(1)/tan(a + b*x))**(n + S(1))/(b*(n + S(1))), x)
def replacement4872(c, m, x, d, a, p, b):
        # rubi.append(4872)
        return Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))/cos(a + b*x)**S(3), x) - Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))/cos(a + b*x), x)
def replacement4873(c, m, n, x, d, a, p, b):
        # rubi.append(4873)
        return -Int((c + d*x)**m*(S(1)/cos(a + b*x))**n*tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*(S(1)/cos(a + b*x))**(n + S(2))*tan(a + b*x)**(p + S(-2)), x)
def replacement4874(c, m, x, d, a, p, b):
        # rubi.append(4874)
        return Int((c + d*x)**m*(S(1)/tan(a + b*x))**(p + S(-2))/sin(a + b*x)**S(3), x) - Int((c + d*x)**m*(S(1)/tan(a + b*x))**(p + S(-2))/sin(a + b*x), x)
def replacement4875(c, m, n, x, d, a, p, b):
        # rubi.append(4875)
        return -Int((c + d*x)**m*(S(1)/sin(a + b*x))**n*(S(1)/tan(a + b*x))**(p + S(-2)), x) + Int((c + d*x)**m*(S(1)/sin(a + b*x))**(n + S(2))*(S(1)/tan(a + b*x))**(p + S(-2)), x)

def With4876(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/cos(a + b*x))**n*tan(a + b*x)**p, x)
        # rubi.append(4876)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)

def With4877(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/sin(a + b*x))**n*(S(1)/tan(a + b*x))**p, x)
        # rubi.append(4877)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
def replacement4878(c, m, n, x, d, a, b):
        # rubi.append(4878)
        return Dist(S(2)**n, Int((c + d*x)**m*(S(1)/sin(S(2)*a + S(2)*b*x))**n, x), x)

def With4879(c, m, n, x, d, a, p, b):
        u = IntHide((S(1)/sin(a + b*x))**n*(S(1)/cos(a + b*x))**p, x)
        # rubi.append(4879)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
def replacement4880(G, m, n, F, x, v, u, w, p):
        # rubi.append(4880)
        return Int(ExpandToSum(u, x)**m*F(ExpandToSum(v, x))**n*G(ExpandToSum(v, x))**p, x)
def replacement4881(c, m, n, x, d, f, a, e, b):
        # rubi.append(4881)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*sin(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*sin(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4882(c, m, n, x, d, f, a, e, b):
        # rubi.append(4882)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b*cos(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b*cos(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4883(c, m, n, x, d, f, a, e, b):
        # rubi.append(4883)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*tan(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*tan(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4884(c, m, n, x, d, f, a, e, b):
        # rubi.append(4884)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b/tan(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b/tan(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4885(c, m, n, x, d, f, a, e, b):
        # rubi.append(4885)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b/cos(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b/cos(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4886(c, m, n, x, d, f, a, e, b):
        # rubi.append(4886)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b/sin(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b/sin(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
def replacement4887(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(4887)
        return Int(ExpandTrigReduce((e + f*x)**m, sin(a + b*x)**p*sin(c + d*x)**q, x), x)
def replacement4888(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(4888)
        return Int(ExpandTrigReduce((e + f*x)**m, cos(a + b*x)**p*cos(c + d*x)**q, x), x)
def replacement4889(c, m, q, x, d, f, a, p, e, b):
        # rubi.append(4889)
        return Int(ExpandTrigReduce((e + f*x)**m, sin(a + b*x)**p*cos(c + d*x)**q, x), x)
def replacement4890(c, G, m, q, F, x, d, f, a, p, e, b):
        # rubi.append(4890)
        return Int(ExpandTrigExpand((e + f*x)**m*G(c + d*x)**q, F, c + d*x, p, b/d, x), x)
def replacement4891(c, F, x, d, a, e, b):
        # rubi.append(4891)
        return -Simp(F**(c*(a + b*x))*e*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
def replacement4892(c, F, x, d, a, e, b):
        # rubi.append(4892)
        return Simp(F**(c*(a + b*x))*e*sin(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
def replacement4893(c, n, F, x, d, a, e, b):
        # rubi.append(4893)
        return Dist(e**S(2)*n*(n + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*sin(d + e*x)**(n + S(-2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**n/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*sin(d + e*x)**(n + S(-1))*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement4894(c, m, F, x, d, a, e, b):
        # rubi.append(4894)
        return Dist(e**S(2)*m*(m + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), Int(F**(c*(a + b*x))*cos(d + e*x)**(m + S(-2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**m/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), x) + Simp(F**(c*(a + b*x))*e*m*sin(d + e*x)*cos(d + e*x)**(m + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), x)
def replacement4895(c, n, F, x, d, a, e, b):
        # rubi.append(4895)
        return Simp(F**(c*(a + b*x))*sin(d + e*x)**(n + S(1))*cos(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement4896(c, n, F, x, d, a, e, b):
        # rubi.append(4896)
        return -Simp(F**(c*(a + b*x))*sin(d + e*x)*cos(d + e*x)**(n + S(1))/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement4897(c, n, F, x, d, a, e, b):
        # rubi.append(4897)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*sin(d + e*x)**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*sin(d + e*x)**(n + S(1))*cos(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement4898(c, n, F, x, d, a, e, b):
        # rubi.append(4898)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*cos(d + e*x)**(n + S(2)), x), x) - Simp(F**(c*(a + b*x))*sin(d + e*x)*cos(d + e*x)**(n + S(1))/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
def replacement4899(c, n, F, x, d, a, e, b):
        # rubi.append(4899)
        return Dist((exp(S(2)*I*(d + e*x)) + S(-1))**(-n)*exp(I*n*(d + e*x))*sin(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(-1))**n*exp(-I*n*(d + e*x)), x), x)
def replacement4900(c, n, F, x, d, a, e, b):
        # rubi.append(4900)
        return Dist((exp(S(2)*I*(d + e*x)) + S(1))**(-n)*exp(I*n*(d + e*x))*cos(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(1))**n*exp(-I*n*(d + e*x)), x), x)
def replacement4901(c, n, F, x, d, a, e, b):
        # rubi.append(4901)
        return Dist(I**n, Int(ExpandIntegrand(F**(c*(a + b*x))*(S(1) - exp(S(2)*I*(d + e*x)))**n*(exp(S(2)*I*(d + e*x)) + S(1))**(-n), x), x), x)
def replacement4902(c, n, F, x, d, a, e, b):
        # rubi.append(4902)
        return Dist((-I)**n, Int(ExpandIntegrand(F**(c*(a + b*x))*(S(1) - exp(S(2)*I*(d + e*x)))**(-n)*(exp(S(2)*I*(d + e*x)) + S(1))**n, x), x), x)
def replacement4903(c, n, F, x, d, a, e, b):
        # rubi.append(4903)
        return Dist(e**S(2)*n*(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*(S(1)/cos(d + e*x))**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*b*c*(S(1)/cos(d + e*x))**n*log(F)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*(S(1)/cos(d + e*x))**(n + S(1))*sin(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement4904(c, n, F, x, d, a, e, b):
        # rubi.append(4904)
        return Dist(e**S(2)*n*(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*(S(1)/sin(d + e*x))**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*b*c*(S(1)/sin(d + e*x))**n*log(F)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) + Simp(F**(c*(a + b*x))*e*n*(S(1)/sin(d + e*x))**(n + S(1))*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
def replacement4905(c, n, F, x, d, a, e, b):
        # rubi.append(4905)
        return Simp(F**(c*(a + b*x))*(S(1)/cos(d + e*x))**(n + S(-1))*sin(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/cos(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement4906(c, n, F, x, d, a, e, b):
        # rubi.append(4906)
        return Simp(F**(c*(a + b*x))*(S(1)/sin(d + e*x))**(n + S(-1))*cos(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/sin(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement4907(c, n, F, x, d, a, e, b):
        # rubi.append(4907)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*(S(1)/cos(d + e*x))**(n + S(-2)), x), x) + Simp(F**(c*(a + b*x))*(S(1)/cos(d + e*x))**(n + S(-1))*sin(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/cos(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement4908(c, n, F, x, d, a, e, b):
        # rubi.append(4908)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*(S(1)/sin(d + e*x))**(n + S(-2)), x), x) - Simp(F**(c*(a + b*x))*(S(1)/sin(d + e*x))**(n + S(-1))*cos(d + e*x)/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*(S(1)/sin(d + e*x))**(n + S(-2))*log(F)/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
def replacement4909(c, n, F, x, d, a, e, b):
        # rubi.append(4909)
        return Simp(S(2)**n*F**(c*(a + b*x))*Hypergeometric2F1(n, -I*b*c*log(F)/(S(2)*e) + n/S(2), -I*b*c*log(F)/(S(2)*e) + n/S(2) + S(1), -exp(S(2)*I*(d + e*x)))*exp(I*n*(d + e*x))/(b*c*log(F) + I*e*n), x)
def replacement4910(c, n, F, x, d, a, e, b):
        # rubi.append(4910)
        return Simp(F**(c*(a + b*x))*(-S(2)*I)**n*Hypergeometric2F1(n, -I*b*c*log(F)/(S(2)*e) + n/S(2), -I*b*c*log(F)/(S(2)*e) + n/S(2) + S(1), exp(S(2)*I*(d + e*x)))*exp(I*n*(d + e*x))/(b*c*log(F) + I*e*n), x)
def replacement4911(c, n, F, x, d, a, e, b):
        # rubi.append(4911)
        return Dist((exp(S(2)*I*(d + e*x)) + S(1))**n*(S(1)/cos(d + e*x))**n*exp(-I*n*(d + e*x)), Int(SimplifyIntegrand(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(1))**(-n)*exp(I*n*(d + e*x)), x), x), x)
def replacement4912(c, n, F, x, d, a, e, b):
        # rubi.append(4912)
        return Dist((S(1) - exp(-S(2)*I*(d + e*x)))**n*(S(1)/sin(d + e*x))**n*exp(I*n*(d + e*x)), Int(SimplifyIntegrand(F**(c*(a + b*x))*(S(1) - exp(-S(2)*I*(d + e*x)))**(-n)*exp(-I*n*(d + e*x)), x), x), x)
def replacement4913(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(4913)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*cos(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement4914(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(4914)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*cos(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement4915(c, n, F, x, d, f, g, a, e, b):
        # rubi.append(4915)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*sin(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
def replacement4916(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(4916)
        return Dist(g**n, Int(F**(c*(a + b*x))*(-tan(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2)))**m, x), x)
def replacement4917(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(4917)
        return Dist(f**n, Int(F**(c*(a + b*x))*tan(d/S(2) + e*x/S(2))**m, x), x)
def replacement4918(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(4918)
        return Dist(f**n, Int(F**(c*(a + b*x))*(S(1)/tan(d/S(2) + e*x/S(2)))**m, x), x)
def replacement4919(c, F, x, d, h, f, i, g, a, e, b):
        # rubi.append(4919)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*cos(d + e*x)/(f + g*sin(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*cos(d + e*x))/(f + g*sin(d + e*x)), x)
def replacement4920(c, F, x, d, h, f, g, i, a, e, b):
        # rubi.append(4920)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*sin(d + e*x)/(f + g*cos(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*sin(d + e*x))/(f + g*cos(d + e*x)), x)
def replacement4921(c, G, n, F, x, v, u):
        # rubi.append(4921)
        return Int(F**(c*ExpandToSum(u, x))*G(ExpandToSum(v, x))**n, x)

def With4922(c, m, n, F, x, d, a, e, b):
        u = IntHide(F**(c*(a + b*x))*sin(d + e*x)**n, x)
        # rubi.append(4922)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)

def With4923(c, m, n, F, x, d, a, e, b):
        u = IntHide(F**(c*(a + b*x))*cos(d + e*x)**n, x)
        # rubi.append(4923)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)
def replacement4924(c, m, n, F, x, d, f, g, a, e, b):
        # rubi.append(4924)
        return Int(ExpandTrigReduce(F**(c*(a + b*x)), sin(d + e*x)**m*cos(f + g*x)**n, x), x)
def replacement4925(c, m, n, F, x, d, f, g, a, p, e, b):
        # rubi.append(4925)
        return Int(ExpandTrigReduce(F**(c*(a + b*x))*x**p, sin(d + e*x)**m*cos(f + g*x)**n, x), x)
def replacement4926(c, G, m, n, F, x, d, H, a, e, b):
        # rubi.append(4926)
        return Int(ExpandTrigToExp(F**(c*(a + b*x)), G(d + e*x)**m*H(d + e*x)**n, x), x)
def replacement4927(n, F, x, v, u):
        # rubi.append(4927)
        return Int(ExpandTrigToExp(F**u, sin(v)**n, x), x)
def replacement4928(n, F, x, v, u):
        # rubi.append(4928)
        return Int(ExpandTrigToExp(F**u, cos(v)**n, x), x)
def replacement4929(m, n, F, x, v, u):
        # rubi.append(4929)
        return Int(ExpandTrigToExp(F**u, sin(v)**m*cos(v)**n, x), x)
def replacement4930(c, n, x, a, p, b):
        # rubi.append(4930)
        return Simp(x*(p + S(2))*sin(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) + Simp(x*sin(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tan(a + b*log(c*x**n))), x)
def replacement4931(c, n, x, a, p, b):
        # rubi.append(4931)
        return Simp(x*(p + S(2))*cos(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement4932(c, n, x, a, p, b):
        # rubi.append(4932)
        return Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/(S(2)*b*n*p) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/(S(2)*b*n*p))**p, x), x)
def replacement4933(c, n, x, a, p, b):
        # rubi.append(4933)
        return Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/S(2) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/S(2))**p, x), x)
def replacement4934(c, n, x, a, b):
        # rubi.append(4934)
        return Simp(x*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x) - Simp(b*n*x*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x)
def replacement4935(c, n, x, a, b):
        # rubi.append(4935)
        return Simp(x*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x) + Simp(b*n*x*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x)
def replacement4936(c, n, x, a, p, b):
        # rubi.append(4936)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(sin(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x*sin(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) - Simp(b*n*p*x*sin(a + b*log(c*x**n))**(p + S(-1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
def replacement4937(c, n, x, a, p, b):
        # rubi.append(4937)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(cos(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x*cos(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) + Simp(b*n*p*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
def replacement4938(c, n, x, a, p, b):
        # rubi.append(4938)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(sin(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x*sin(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) + Simp(x*sin(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tan(a + b*log(c*x**n))), x)
def replacement4939(c, n, x, a, p, b):
        # rubi.append(4939)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(cos(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement4940(c, n, x, a, p, b):
        # rubi.append(4940)
        return Simp(x*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*(-I*(c*x**n)**(I*b)*exp(I*a) + I*(c*x**n)**(-I*b)*exp(-I*a))**p*Hypergeometric2F1(-p, -I*(-I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + S(1)), x)
def replacement4941(c, n, x, a, p, b):
        # rubi.append(4941)
        return Simp(x*((c*x**n)**(I*b)*exp(I*a) + (c*x**n)**(-I*b)*exp(-I*a))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*Hypergeometric2F1(-p, -I*(-I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + S(1)), x)
def replacement4942(c, m, n, x, a, p, b):
        # rubi.append(4942)
        return Simp(x**(m + S(1))*(p + S(2))*sin(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) + Simp(x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tan(a + b*log(c*x**n))), x)
def replacement4943(c, m, n, x, a, p, b):
        # rubi.append(4943)
        return Simp(x**(m + S(1))*(p + S(2))*cos(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) - Simp(x**(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
def replacement4944(c, m, n, x, a, p, b):
        # rubi.append(4944)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*(-(c*x**n)**((m + S(1))/(n*p))*(m + S(1))*exp(-a*b*n*p/(m + S(1)))/(b*n*p) + (c*x**n)**(-(m + S(1))/(n*p))*(m + S(1))*exp(a*b*n*p/(m + S(1)))/(b*n*p))**p, x), x), x)
def replacement4945(c, m, n, x, a, p, b):
        # rubi.append(4945)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*(-(c*x**n)**((m + S(1))/(n*p))*exp(-a*b*n*p/(m + S(1))) + (c*x**n)**(-(m + S(1))/(n*p))*exp(a*b*n*p/(m + S(1))))**p, x), x), x)
def replacement4946(c, m, n, x, a, b):
        # rubi.append(4946)
        return Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*x**(m + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x)
def replacement4947(c, m, n, x, a, b):
        # rubi.append(4947)
        return Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*x**(m + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x)
def replacement4948(c, m, n, x, a, p, b):
        # rubi.append(4948)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*sin(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(-1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
def replacement4949(c, m, n, x, a, p, b):
        # rubi.append(4949)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*cos(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
def replacement4950(c, m, n, x, a, p, b):
        # rubi.append(4950)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*sin(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))/(b*n*(p + S(1))*tan(a + b*log(c*x**n))), x) - Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement4951(c, m, n, x, a, p, b):
        # rubi.append(4951)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*cos(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x**(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x) - Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
def replacement4952(c, m, n, x, a, p, b):
        # rubi.append(4952)
        return Simp(x**(m + S(1))*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*(-I*(c*x**n)**(I*b)*exp(I*a) + I*(c*x**n)**(-I*b)*exp(-I*a))**p*Hypergeometric2F1(-p, -I*(-I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + m + S(1)), x)
def replacement4953(c, m, n, x, a, p, b):
        # rubi.append(4953)
        return Simp(x**(m + S(1))*((c*x**n)**(I*b)*exp(I*a) + (c*x**n)**(-I*b)*exp(-I*a))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*Hypergeometric2F1(-p, -I*(-I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + m + S(1)), x)
def replacement4954(c, n, x, a, b):
        # rubi.append(4954)
        return Dist(S(2)*exp(a*b*n), Int((c*x**n)**(S(1)/n)/((c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x), x)
def replacement4955(c, n, x, a, b):
        # rubi.append(4955)
        return Dist(S(2)*b*n*exp(a*b*n), Int((c*x**n)**(S(1)/n)/(-(c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x), x)
def replacement4956(c, n, x, a, p, b):
        # rubi.append(4956)
        return Simp(x*(p + S(-2))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))/(p + S(-1)), x) + Simp(x*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))*tan(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement4957(c, n, x, a, p, b):
        # rubi.append(4957)
        return Simp(x*(p + S(-2))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(p + S(-1)), x) - Simp(x*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tan(a + b*log(c*x**n))), x)
def replacement4958(c, n, x, a, p, b):
        # rubi.append(4958)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int((S(1)/cos(a + b*log(c*x**n)))**(p + S(-2)), x), x) - Simp(x*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) + Simp(x*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))*tan(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement4959(c, n, x, a, p, b):
        # rubi.append(4959)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int((S(1)/sin(a + b*log(c*x**n)))**(p + S(-2)), x), x) - Simp(x*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) - Simp(x*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tan(a + b*log(c*x**n))), x)
def replacement4960(c, n, x, a, p, b):
        # rubi.append(4960)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int((S(1)/cos(a + b*log(c*x**n)))**(p + S(2)), x), x) + Simp(x*(S(1)/cos(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) - Simp(b*n*p*x*(S(1)/cos(a + b*log(c*x**n)))**(p + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
def replacement4961(c, n, x, a, p, b):
        # rubi.append(4961)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int((S(1)/sin(a + b*log(c*x**n)))**(p + S(2)), x), x) + Simp(x*(S(1)/sin(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) + Simp(b*n*p*x*(S(1)/sin(a + b*log(c*x**n)))**(p + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
def replacement4962(c, n, x, a, p, b):
        # rubi.append(4962)
        return Simp(x*((c*x**n)**(I*b)*exp(I*a)/((c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + S(1)), x)
def replacement4963(c, n, x, a, p, b):
        # rubi.append(4963)
        return Simp(x*(-I*(c*x**n)**(I*b)*exp(I*a)/(-(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + S(1)), x)
def replacement4964(c, m, n, x, a, b):
        # rubi.append(4964)
        return Dist(S(2)*exp(a*b*n/(m + S(1))), Int(x**m*(c*x**n)**((m + S(1))/n)/((c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x), x)
def replacement4965(c, m, n, x, a, b):
        # rubi.append(4965)
        return Dist(S(2)*b*n*exp(a*b*n/(m + S(1)))/(m + S(1)), Int(x**m*(c*x**n)**((m + S(1))/n)/(-(c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x), x)
def replacement4966(c, m, n, x, a, p, b):
        # rubi.append(4966)
        return Simp(x**(m + S(1))*(p + S(-2))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) + Simp(x**(m + S(1))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))*tan(a + b*log(c*x**n))/(b*n*(p + S(-1))), x)
def replacement4967(c, m, n, x, a, p, b):
        # rubi.append(4967)
        return Simp(x**(m + S(1))*(p + S(-2))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) - Simp(x**(m + S(1))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tan(a + b*log(c*x**n))), x)
def replacement4968(c, m, n, x, a, p, b):
        # rubi.append(4968)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))*tan(a + b*log(c*x**n))/(b*n*(p + S(-1))), x) - Simp(x**(m + S(1))*(m + S(1))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
def replacement4969(c, m, n, x, a, p, b):
        # rubi.append(4969)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2)), x), x) - Simp(x**(m + S(1))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b*n*(p + S(-1))*tan(a + b*log(c*x**n))), x) - Simp(x**(m + S(1))*(m + S(1))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
def replacement4970(c, m, n, x, a, p, b):
        # rubi.append(4970)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*(S(1)/cos(a + b*log(c*x**n)))**(p + S(2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*(S(1)/cos(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*(S(1)/cos(a + b*log(c*x**n)))**(p + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
def replacement4971(c, m, n, x, a, p, b):
        # rubi.append(4971)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*(S(1)/sin(a + b*log(c*x**n)))**(p + S(2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*(S(1)/sin(a + b*log(c*x**n)))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*(S(1)/sin(a + b*log(c*x**n)))**(p + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
def replacement4972(c, m, n, x, a, p, b):
        # rubi.append(4972)
        return Simp(x**(m + S(1))*((c*x**n)**(I*b)*exp(I*a)/((c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + m + S(1)), x)
def replacement4973(c, m, n, x, a, p, b):
        # rubi.append(4973)
        return Simp(x**(m + S(1))*(-I*(c*x**n)**(I*b)*exp(I*a)/(-(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + m + S(1)), x)
def replacement4974(x, a, p, b):
        # rubi.append(4974)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*sin(a*x*log(b*x)**p), x), x) - Simp(cos(a*x*log(b*x)**p)/a, x)
def replacement4975(x, a, p, b):
        # rubi.append(4975)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*cos(a*x*log(b*x)**p), x), x) + Simp(sin(a*x*log(b*x)**p)/a, x)
def replacement4976(n, x, a, p, b):
        # rubi.append(4976)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) - Dist((n + S(-1))/(a*n), Int(x**(-n)*cos(a*x**n*log(b*x)**p), x), x) - Simp(x**(S(1) - n)*cos(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4977(n, x, a, p, b):
        # rubi.append(4977)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) + Dist((n + S(-1))/(a*n), Int(x**(-n)*sin(a*x**n*log(b*x)**p), x), x) + Simp(x**(S(1) - n)*sin(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4978(m, n, x, a, p, b):
        # rubi.append(4978)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) - Simp(cos(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4979(m, n, x, a, p, b):
        # rubi.append(4979)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) + Simp(sin(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4980(m, n, x, a, p, b):
        # rubi.append(4980)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) + Dist((m - n + S(1))/(a*n), Int(x**(m - n)*cos(a*x**n*log(b*x)**p), x), x) - Simp(x**(m - n + S(1))*cos(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4981(m, n, x, a, p, b):
        # rubi.append(4981)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) - Dist((m - n + S(1))/(a*n), Int(x**(m - n)*sin(a*x**n*log(b*x)**p), x), x) + Simp(x**(m - n + S(1))*sin(a*x**n*log(b*x)**p)/(a*n), x)
def replacement4982(c, n, x, d, a):
        # rubi.append(4982)
        return -Dist(S(1)/d, Subst(Int(sin(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement4983(c, n, x, d, a):
        # rubi.append(4983)
        return -Dist(S(1)/d, Subst(Int(cos(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement4984(c, n, x, d, a, e, b):
        # rubi.append(4984)
        return -Dist(S(1)/d, Subst(Int(sin(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
def replacement4985(c, n, x, d, a, e, b):
        # rubi.append(4985)
        return -Dist(S(1)/d, Subst(Int(cos(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)

def With4986(x, u, n):
        lst = QuotientOfLinearsParts(u, x)
        # rubi.append(4986)
        return Int(sin((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)

def With4987(x, u, n):
        lst = QuotientOfLinearsParts(u, x)
        # rubi.append(4987)
        return Int(cos((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
def replacement4988(q, x, v, u, w, p):
        # rubi.append(4988)
        return Int(u*sin(v)**(p + q), x)
def replacement4989(q, x, v, u, w, p):
        # rubi.append(4989)
        return Int(u*cos(v)**(p + q), x)
def replacement4990(q, x, v, w, p):
        # rubi.append(4990)
        return Int(ExpandTrigReduce(sin(v)**p*sin(w)**q, x), x)
def replacement4991(q, x, v, w, p):
        # rubi.append(4991)
        return Int(ExpandTrigReduce(cos(v)**p*cos(w)**q, x), x)
def replacement4992(m, q, x, v, w, p):
        # rubi.append(4992)
        return Int(ExpandTrigReduce(x**m, sin(v)**p*sin(w)**q, x), x)
def replacement4993(m, q, x, v, w, p):
        # rubi.append(4993)
        return Int(ExpandTrigReduce(x**m, cos(v)**p*cos(w)**q, x), x)
def replacement4994(x, v, u, w, p):
        # rubi.append(4994)
        return Dist(S(2)**(-p), Int(u*sin(S(2)*v)**p, x), x)
def replacement4995(q, x, v, w, p):
        # rubi.append(4995)
        return Int(ExpandTrigReduce(sin(v)**p*cos(w)**q, x), x)
def replacement4996(m, q, x, v, w, p):
        # rubi.append(4996)
        return Int(ExpandTrigReduce(x**m, sin(v)**p*cos(w)**q, x), x)
def replacement4997(x, w, v, n):
        # rubi.append(4997)
        return Dist(cos(v - w), Int(tan(w)**(n + S(-1))/cos(w), x), x) - Int(cos(v)*tan(w)**(n + S(-1)), x)
def replacement4998(x, w, v, n):
        # rubi.append(4998)
        return Dist(cos(v - w), Int((S(1)/tan(w))**(n + S(-1))/sin(w), x), x) - Int((S(1)/tan(w))**(n + S(-1))*sin(v), x)
def replacement4999(x, w, v, n):
        # rubi.append(4999)
        return Dist(sin(v - w), Int((S(1)/tan(w))**(n + S(-1))/sin(w), x), x) + Int((S(1)/tan(w))**(n + S(-1))*cos(v), x)
def replacement5000(x, w, v, n):
        # rubi.append(5000)
        return -Dist(sin(v - w), Int(tan(w)**(n + S(-1))/cos(w), x), x) + Int(sin(v)*tan(w)**(n + S(-1)), x)
def replacement5001(x, w, v, n):
        # rubi.append(5001)
        return Dist(sin(v - w), Int((S(1)/cos(w))**(n + S(-1)), x), x) + Dist(cos(v - w), Int((S(1)/cos(w))**(n + S(-1))*tan(w), x), x)
def replacement5002(x, w, v, n):
        # rubi.append(5002)
        return -Dist(sin(v - w), Int((S(1)/sin(w))**(n + S(-1)), x), x) + Dist(cos(v - w), Int((S(1)/sin(w))**(n + S(-1))/tan(w), x), x)
def replacement5003(x, w, v, n):
        # rubi.append(5003)
        return Dist(sin(v - w), Int((S(1)/sin(w))**(n + S(-1))/tan(w), x), x) + Dist(cos(v - w), Int((S(1)/sin(w))**(n + S(-1)), x), x)
def replacement5004(x, w, v, n):
        # rubi.append(5004)
        return -Dist(sin(v - w), Int((S(1)/cos(w))**(n + S(-1))*tan(w), x), x) + Dist(cos(v - w), Int((S(1)/cos(w))**(n + S(-1)), x), x)
def replacement5005(c, m, n, x, d, f, a, e, b):
        # rubi.append(5005)
        return Int((a + b*sin(S(2)*c + S(2)*d*x)/S(2))**n*(e + f*x)**m, x)
def replacement5006(c, m, n, x, d, a, b):
        # rubi.append(5006)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a - b*cos(S(2)*c + S(2)*d*x) + b)**n, x), x)
def replacement5007(c, m, n, x, d, a, b):
        # rubi.append(5007)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a + b*cos(S(2)*c + S(2)*d*x) + b)**n, x), x)
def replacement5008(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(5008)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*sin(a + b*x**n)**p, x), x, c + d*x), x)
def replacement5009(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(5009)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*cos(a + b*x**n)**p, x), x, c + d*x), x)
def replacement5010(c, m, x, d, f, g, a, e, b):
        # rubi.append(5010)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
def replacement5011(c, m, x, d, f, g, e, b):
        # rubi.append(5011)
        return Dist(S(2), Int((f + g*x)**m/(b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
def replacement5012(c, m, x, d, f, g, a, e, b):
        # rubi.append(5012)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
def replacement5013(c, m, x, d, f, g, e, b):
        # rubi.append(5013)
        return Dist(S(2), Int((f + g*x)**m/(b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
def replacement5014(c, m, x, d, f, g, a, e, b):
        # rubi.append(5014)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
def replacement5015(c, m, x, d, f, a, e, b):
        # rubi.append(5015)
        return Int((e + f*x)**m*exp(I*(c + d*x))/(a - I*b*exp(I*(c + d*x)) - Rt(a**S(2) - b**S(2), S(2))), x) + Int((e + f*x)**m*exp(I*(c + d*x))/(a - I*b*exp(I*(c + d*x)) + Rt(a**S(2) - b**S(2), S(2))), x) - Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement5016(c, m, x, d, f, a, e, b):
        # rubi.append(5016)
        return -Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(a + b*exp(I*(c + d*x)) - Rt(a**S(2) - b**S(2), S(2))), x), x) - Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(a + b*exp(I*(c + d*x)) + Rt(a**S(2) - b**S(2), S(2))), x), x) + Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement5017(c, m, x, d, f, a, e, b):
        # rubi.append(5017)
        return Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + b*exp(I*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x), x) + Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + b*exp(I*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x), x) - Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement5018(c, m, x, d, f, a, e, b):
        # rubi.append(5018)
        return Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + I*b*exp(I*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x) + Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + I*b*exp(I*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x) + Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
def replacement5019(c, m, n, x, d, f, a, e, b):
        # rubi.append(5019)
        return Dist(S(1)/a, Int((e + f*x)**m*cos(c + d*x)**(n + S(-2)), x), x) - Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)*cos(c + d*x)**(n + S(-2)), x), x)
def replacement5020(c, m, n, x, d, f, a, e, b):
        # rubi.append(5020)
        return Dist(S(1)/a, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2)), x), x) - Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))*cos(c + d*x), x), x)
def replacement5021(c, m, n, x, d, f, a, e, b):
        # rubi.append(5021)
        return -Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)*cos(c + d*x)**(n + S(-2)), x), x) + Dist(a/b**S(2), Int((e + f*x)**m*cos(c + d*x)**(n + S(-2)), x), x) - Dist((a**S(2) - b**S(2))/b**S(2), Int((e + f*x)**m*cos(c + d*x)**(n + S(-2))/(a + b*sin(c + d*x)), x), x)
def replacement5022(c, m, n, x, d, f, a, e, b):
        # rubi.append(5022)
        return -Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))*cos(c + d*x), x), x) + Dist(a/b**S(2), Int((e + f*x)**m*sin(c + d*x)**(n + S(-2)), x), x) - Dist((a**S(2) - b**S(2))/b**S(2), Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))/(a + b*cos(c + d*x)), x), x)
def replacement5023(c, x, d, B, A, f, a, e, b):
        # rubi.append(5023)
        return Dist(B*f/(a*d), Int(cos(c + d*x)/(a + b*sin(c + d*x)), x), x) - Simp(B*(e + f*x)*cos(c + d*x)/(a*d*(a + b*sin(c + d*x))), x)
def replacement5024(c, x, d, B, A, f, a, e, b):
        # rubi.append(5024)
        return -Dist(B*f/(a*d), Int(sin(c + d*x)/(a + b*cos(c + d*x)), x), x) + Simp(B*(e + f*x)*sin(c + d*x)/(a*d*(a + b*cos(c + d*x))), x)
def replacement5025(m, n, x, v, a, b):
        # rubi.append(5025)
        return Int((a*cos(v) + b*sin(v))**n, x)
def replacement5026(m, n, x, v, a, b):
        # rubi.append(5026)
        return Int((a*sin(v) + b*cos(v))**n, x)
def replacement5027(c, m, n, x, d, u, a, b):
        # rubi.append(5027)
        return Int(ExpandTrigReduce(u, sin(a + b*x)**m*sin(c + d*x)**n, x), x)
def replacement5028(c, m, n, x, d, u, a, b):
        # rubi.append(5028)
        return Int(ExpandTrigReduce(u, cos(a + b*x)**m*cos(c + d*x)**n, x), x)
def replacement5029(c, x, d, a, b):
        # rubi.append(5029)
        return Dist(S(1)/sin((-a*d + b*c)/b), Int(tan(c + d*x), x), x) - Dist(S(1)/sin((-a*d + b*c)/d), Int(tan(a + b*x), x), x)
def replacement5030(c, x, d, a, b):
        # rubi.append(5030)
        return Dist(S(1)/sin((-a*d + b*c)/b), Int(S(1)/tan(a + b*x), x), x) - Dist(S(1)/sin((-a*d + b*c)/d), Int(S(1)/tan(c + d*x), x), x)
def replacement5031(c, x, d, a, b):
        # rubi.append(5031)
        return Dist(b*cos((-a*d + b*c)/d)/d, Int(S(1)/(cos(a + b*x)*cos(c + d*x)), x), x) - Simp(b*x/d, x)
def replacement5032(c, x, d, a, b):
        # rubi.append(5032)
        return Dist(cos((-a*d + b*c)/d), Int(S(1)/(sin(a + b*x)*sin(c + d*x)), x), x) - Simp(b*x/d, x)
def replacement5033(n, x, v, u, a, b):
        # rubi.append(5033)
        return Int(u*(a*exp(-a*v/b))**n, x)
