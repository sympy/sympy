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

def miscellaneous_trig(rubi):
    from sympy.integrals.rubi.constraints import cons1638, cons18, cons2, cons3, cons7, cons27, cons21, cons4, cons34, cons35, cons36, cons1639, cons1640, cons1641, cons1642, cons1643, cons1644, cons1645, cons25, cons1646, cons1400, cons208, cons38, cons48, cons125, cons147, cons343, cons5, cons240, cons244, cons1325, cons137, cons1647, cons1280, cons166, cons319, cons1648, cons31, cons249, cons94, cons253, cons13, cons163, cons246, cons1270, cons1649, cons1650, cons1651, cons170, cons1652, cons1653, cons93, cons89, cons1654, cons162, cons88, cons1655, cons1656, cons85, cons128, cons1471, cons744, cons1474, cons1657, cons23, cons1658, cons1659, cons1660, cons1661, cons1239, cons1662, cons1663, cons1664, cons1665, cons555, cons1666, cons628, cons10, cons1667, cons1668, cons1669, cons66, cons1222, cons376, cons49, cons50, cons51, cons52, cons1670, cons1431, cons1671, cons1672, cons1673, cons1674, cons62, cons584, cons464, cons1675, cons168, cons1676, cons1677, cons1678, cons1679, cons1680, cons812, cons813, cons17, cons1681, cons1682, cons1683, cons1684, cons1091, cons1685, cons87, cons165, cons1686, cons1687, cons1387, cons1688, cons1434, cons1689, cons1494, cons963, cons1690, cons1636, cons1691, cons196, cons1692, cons1011, cons150, cons1543, cons1693, cons1694, cons209, cons224, cons1695, cons810, cons811, cons148, cons528, cons1696, cons1697, cons1698, cons1699, cons1700, cons54, cons1701, cons1702, cons146, cons1703, cons1497, cons1704, cons1705, cons1706, cons1707, cons1708, cons1709, cons1710, cons1711, cons1637, cons1712, cons1713, cons1714, cons1715, cons1716, cons338, cons53, cons627, cons71, cons1717, cons1718, cons1719, cons1720, cons1352, cons1470, cons463, cons1721, cons1722, cons1723, cons1724, cons1257, cons1259, cons1466, cons1473, cons1725

    pattern4666 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1638, cons18)
    def replacement4666(a, n, m, b, u, x, d, c):
        rubi.append(4666)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**(m + n)*(d*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4666 = ReplacementRule(pattern4666, replacement4666)
    pattern4667 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1638, cons18)
    def replacement4667(a, n, m, b, u, x, d, c):
        rubi.append(4667)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4667 = ReplacementRule(pattern4667, replacement4667)
    pattern4668 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1638, cons18)
    def replacement4668(a, n, m, b, u, x, d, c):
        rubi.append(4668)
        return Dist((c*cot(a + b*x))**m*(d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m), Int((d*sin(a + b*x))**(-m + n)*(d*cos(a + b*x))**m*ActivateTrig(u), x), x)
    rule4668 = ReplacementRule(pattern4668, replacement4668)
    pattern4669 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1638, cons18)
    def replacement4669(a, n, m, b, u, x, d, c):
        rubi.append(4669)
        return Dist((c*cot(a + b*x))**m*(d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m), Int((d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**(m + n)*ActivateTrig(u), x), x)
    rule4669 = ReplacementRule(pattern4669, replacement4669)
    pattern4670 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1638)
    def replacement4670(a, n, m, b, u, x, d, c):
        rubi.append(4670)
        return Dist((c*csc(a + b*x))**m*(d*sin(a + b*x))**m, Int((d*sin(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4670 = ReplacementRule(pattern4670, replacement4670)
    pattern4671 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1638)
    def replacement4671(a, m, b, u, x, c):
        rubi.append(4671)
        return Dist((c*sin(a + b*x))**(-m)*(c*cos(a + b*x))**m*(c*tan(a + b*x))**m, Int((c*sin(a + b*x))**m*(c*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4671 = ReplacementRule(pattern4671, replacement4671)
    pattern4672 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1638)
    def replacement4672(a, m, b, u, x, c):
        rubi.append(4672)
        return Dist((c*sin(a + b*x))**m*(c*cos(a + b*x))**(-m)*(c*cot(a + b*x))**m, Int((c*sin(a + b*x))**(-m)*(c*cos(a + b*x))**m*ActivateTrig(u), x), x)
    rule4672 = ReplacementRule(pattern4672, replacement4672)
    pattern4673 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1638)
    def replacement4673(a, m, b, u, x, c):
        rubi.append(4673)
        return Dist((c*cos(a + b*x))**m*(c*sec(a + b*x))**m, Int((c*cos(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4673 = ReplacementRule(pattern4673, replacement4673)
    pattern4674 = Pattern(Integral(u_*(WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1638)
    def replacement4674(a, m, b, u, x, c):
        rubi.append(4674)
        return Dist((c*sin(a + b*x))**m*(c*csc(a + b*x))**m, Int((c*sin(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4674 = ReplacementRule(pattern4674, replacement4674)
    pattern4675 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1638)
    def replacement4675(a, n, B, b, A, u, x, c):
        rubi.append(4675)
        return Dist(c, Int((c*sin(a + b*x))**(n + S(-1))*(A*sin(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4675 = ReplacementRule(pattern4675, replacement4675)
    pattern4676 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1638)
    def replacement4676(a, n, B, b, A, u, x, c):
        rubi.append(4676)
        return Dist(c, Int((c*cos(a + b*x))**(n + S(-1))*(A*cos(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4676 = ReplacementRule(pattern4676, replacement4676)
    pattern4677 = Pattern(Integral(u_*(A_ + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1638)
    def replacement4677(a, B, b, A, u, x):
        rubi.append(4677)
        return Int((A*sin(a + b*x) + B)*ActivateTrig(u)/sin(a + b*x), x)
    rule4677 = ReplacementRule(pattern4677, replacement4677)
    pattern4678 = Pattern(Integral(u_*(A_ + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1638)
    def replacement4678(a, B, b, A, u, x):
        rubi.append(4678)
        return Int((A*cos(a + b*x) + B)*ActivateTrig(u)/cos(a + b*x), x)
    rule4678 = ReplacementRule(pattern4678, replacement4678)
    pattern4679 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1638)
    def replacement4679(a, n, B, b, A, u, x, C, c):
        rubi.append(4679)
        return Dist(c**S(2), Int((c*sin(a + b*x))**(n + S(-2))*(A*sin(a + b*x)**S(2) + B*sin(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4679 = ReplacementRule(pattern4679, replacement4679)
    pattern4680 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1638)
    def replacement4680(a, n, B, b, A, u, x, C, c):
        rubi.append(4680)
        return Dist(c**S(2), Int((c*cos(a + b*x))**(n + S(-2))*(A*cos(a + b*x)**S(2) + B*cos(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4680 = ReplacementRule(pattern4680, replacement4680)
    pattern4681 = Pattern(Integral((WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1638)
    def replacement4681(a, n, A, b, C, u, x, c):
        rubi.append(4681)
        return Dist(c**S(2), Int((c*sin(a + b*x))**(n + S(-2))*(A*sin(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4681 = ReplacementRule(pattern4681, replacement4681)
    pattern4682 = Pattern(Integral((WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1638)
    def replacement4682(a, n, A, b, C, u, x, c):
        rubi.append(4682)
        return Dist(c**S(2), Int((c*cos(a + b*x))**(n + S(-2))*(A*cos(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4682 = ReplacementRule(pattern4682, replacement4682)
    pattern4683 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1638)
    def replacement4683(a, B, b, A, u, x, C):
        rubi.append(4683)
        return Int((A*sin(a + b*x)**S(2) + B*sin(a + b*x) + C)*ActivateTrig(u)/sin(a + b*x)**S(2), x)
    rule4683 = ReplacementRule(pattern4683, replacement4683)
    pattern4684 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1638)
    def replacement4684(a, B, A, b, C, u, x):
        rubi.append(4684)
        return Int((A*cos(a + b*x)**S(2) + B*cos(a + b*x) + C)*ActivateTrig(u)/cos(a + b*x)**S(2), x)
    rule4684 = ReplacementRule(pattern4684, replacement4684)
    pattern4685 = Pattern(Integral(u_*(A_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1638)
    def replacement4685(a, A, b, C, u, x):
        rubi.append(4685)
        return Int((A*sin(a + b*x)**S(2) + C)*ActivateTrig(u)/sin(a + b*x)**S(2), x)
    rule4685 = ReplacementRule(pattern4685, replacement4685)
    pattern4686 = Pattern(Integral(u_*(A_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1638)
    def replacement4686(a, A, b, C, u, x):
        rubi.append(4686)
        return Int((A*cos(a + b*x)**S(2) + C)*ActivateTrig(u)/cos(a + b*x)**S(2), x)
    rule4686 = ReplacementRule(pattern4686, replacement4686)
    pattern4687 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons36, cons1639)
    def replacement4687(a, B, A, b, C, u, x):
        rubi.append(4687)
        return Int((A*sin(a + b*x) + B*sin(a + b*x)**S(2) + C)*ActivateTrig(u)/sin(a + b*x), x)
    rule4687 = ReplacementRule(pattern4687, replacement4687)
    pattern4688 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons36, cons1639)
    def replacement4688(a, B, b, A, u, x, C):
        rubi.append(4688)
        return Int((A*cos(a + b*x) + B*cos(a + b*x)**S(2) + C)*ActivateTrig(u)/cos(a + b*x), x)
    rule4688 = ReplacementRule(pattern4688, replacement4688)
    pattern4689 = Pattern(Integral(u_*(WC('A', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4689(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4689)
        return Int((A + B*sin(a + b*x) + C*sin(a + b*x)**S(2))*ActivateTrig(u)*sin(a + b*x)**n, x)
    rule4689 = ReplacementRule(pattern4689, replacement4689)
    pattern4690 = Pattern(Integral(u_*(WC('A', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4690(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4690)
        return Int((A + B*cos(a + b*x) + C*cos(a + b*x)**S(2))*ActivateTrig(u)*cos(a + b*x)**n, x)
    rule4690 = ReplacementRule(pattern4690, replacement4690)
    pattern4691 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1642)
    def replacement4691(a, n, m, b, u, x, d, c):
        rubi.append(4691)
        return Dist((c*cot(a + b*x))**m*(d*tan(a + b*x))**m, Int((d*tan(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4691 = ReplacementRule(pattern4691, replacement4691)
    pattern4692 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1643)
    def replacement4692(a, n, m, b, u, x, d, c):
        rubi.append(4692)
        return Dist((c*tan(a + b*x))**m*(d*sin(a + b*x))**(-m)*(d*cos(a + b*x))**m, Int((d*sin(a + b*x))**m*(d*cos(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4692 = ReplacementRule(pattern4692, replacement4692)
    pattern4693 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1642)
    def replacement4693(a, m, b, u, x, c):
        rubi.append(4693)
        return Dist((c*tan(a + b*x))**m*(c*cot(a + b*x))**m, Int((c*tan(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4693 = ReplacementRule(pattern4693, replacement4693)
    pattern4694 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1643)
    def replacement4694(a, m, b, u, x, c):
        rubi.append(4694)
        return Dist((c*tan(a + b*x))**m*(c*cot(a + b*x))**m, Int((c*cot(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4694 = ReplacementRule(pattern4694, replacement4694)
    pattern4695 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1642)
    def replacement4695(a, n, B, b, A, u, x, c):
        rubi.append(4695)
        return Dist(c, Int((c*tan(a + b*x))**(n + S(-1))*(A*tan(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4695 = ReplacementRule(pattern4695, replacement4695)
    pattern4696 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1643)
    def replacement4696(a, n, B, b, A, u, x, c):
        rubi.append(4696)
        return Dist(c, Int((c*cot(a + b*x))**(n + S(-1))*(A*cot(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4696 = ReplacementRule(pattern4696, replacement4696)
    pattern4697 = Pattern(Integral(u_*(A_ + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1642)
    def replacement4697(a, B, b, A, u, x):
        rubi.append(4697)
        return Int((A*tan(a + b*x) + B)*ActivateTrig(u)/tan(a + b*x), x)
    rule4697 = ReplacementRule(pattern4697, replacement4697)
    pattern4698 = Pattern(Integral(u_*(A_ + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1643)
    def replacement4698(a, B, b, A, u, x):
        rubi.append(4698)
        return Int((A*cot(a + b*x) + B)*ActivateTrig(u)/cot(a + b*x), x)
    rule4698 = ReplacementRule(pattern4698, replacement4698)
    pattern4699 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1642)
    def replacement4699(a, n, B, b, A, u, x, C, c):
        rubi.append(4699)
        return Dist(c**S(2), Int((c*tan(a + b*x))**(n + S(-2))*(A*tan(a + b*x)**S(2) + B*tan(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4699 = ReplacementRule(pattern4699, replacement4699)
    pattern4700 = Pattern(Integral((WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1643)
    def replacement4700(a, n, B, b, A, u, x, C, c):
        rubi.append(4700)
        return Dist(c**S(2), Int((c*cot(a + b*x))**(n + S(-2))*(A*cot(a + b*x)**S(2) + B*cot(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4700 = ReplacementRule(pattern4700, replacement4700)
    pattern4701 = Pattern(Integral((WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1642)
    def replacement4701(a, n, A, b, C, u, x, c):
        rubi.append(4701)
        return Dist(c**S(2), Int((c*tan(a + b*x))**(n + S(-2))*(A*tan(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4701 = ReplacementRule(pattern4701, replacement4701)
    pattern4702 = Pattern(Integral((WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1643)
    def replacement4702(a, n, A, b, C, u, x, c):
        rubi.append(4702)
        return Dist(c**S(2), Int((c*cot(a + b*x))**(n + S(-2))*(A*cot(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4702 = ReplacementRule(pattern4702, replacement4702)
    pattern4703 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1642)
    def replacement4703(a, B, b, A, u, x, C):
        rubi.append(4703)
        return Int((A*tan(a + b*x)**S(2) + B*tan(a + b*x) + C)*ActivateTrig(u)/tan(a + b*x)**S(2), x)
    rule4703 = ReplacementRule(pattern4703, replacement4703)
    pattern4704 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1643)
    def replacement4704(a, B, A, b, C, u, x):
        rubi.append(4704)
        return Int((A*cot(a + b*x)**S(2) + B*cot(a + b*x) + C)*ActivateTrig(u)/cot(a + b*x)**S(2), x)
    rule4704 = ReplacementRule(pattern4704, replacement4704)
    pattern4705 = Pattern(Integral(u_*(A_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1642)
    def replacement4705(a, A, b, C, u, x):
        rubi.append(4705)
        return Int((A*tan(a + b*x)**S(2) + C)*ActivateTrig(u)/tan(a + b*x)**S(2), x)
    rule4705 = ReplacementRule(pattern4705, replacement4705)
    pattern4706 = Pattern(Integral(u_*(A_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1643)
    def replacement4706(a, A, b, C, u, x):
        rubi.append(4706)
        return Int((A*cot(a + b*x)**S(2) + C)*ActivateTrig(u)/cot(a + b*x)**S(2), x)
    rule4706 = ReplacementRule(pattern4706, replacement4706)
    pattern4707 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons36, cons1639)
    def replacement4707(a, B, A, b, C, u, x):
        rubi.append(4707)
        return Int((A*tan(a + b*x) + B*tan(a + b*x)**S(2) + C)*ActivateTrig(u)/tan(a + b*x), x)
    rule4707 = ReplacementRule(pattern4707, replacement4707)
    pattern4708 = Pattern(Integral(u_*(WC('A', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4708(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4708)
        return Int((A + B*tan(a + b*x) + C*tan(a + b*x)**S(2))*ActivateTrig(u)*tan(a + b*x)**n, x)
    rule4708 = ReplacementRule(pattern4708, replacement4708)
    pattern4709 = Pattern(Integral(u_*(WC('A', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4709(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4709)
        return Int((A + B*cot(a + b*x) + C*cot(a + b*x)**S(2))*ActivateTrig(u)*cot(a + b*x)**n, x)
    rule4709 = ReplacementRule(pattern4709, replacement4709)
    pattern4710 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644)
    def replacement4710(a, n, m, b, u, x, d, c):
        rubi.append(4710)
        return Dist((c*sin(a + b*x))**m*(d*csc(a + b*x))**m, Int((d*csc(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4710 = ReplacementRule(pattern4710, replacement4710)
    pattern4711 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644)
    def replacement4711(a, n, m, b, u, x, d, c):
        rubi.append(4711)
        return Dist((c*cos(a + b*x))**m*(d*sec(a + b*x))**m, Int((d*sec(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4711 = ReplacementRule(pattern4711, replacement4711)
    pattern4712 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644, cons18)
    def replacement4712(a, n, m, b, u, x, d, c):
        rubi.append(4712)
        return Dist((c*tan(a + b*x))**m*(d*csc(a + b*x))**m*(d*sec(a + b*x))**(-m), Int((d*csc(a + b*x))**(-m)*(d*sec(a + b*x))**(m + n)*ActivateTrig(u), x), x)
    rule4712 = ReplacementRule(pattern4712, replacement4712)
    pattern4713 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644, cons18)
    def replacement4713(a, n, m, b, u, x, d, c):
        rubi.append(4713)
        return Dist((c*tan(a + b*x))**m*(d*csc(a + b*x))**m*(d*sec(a + b*x))**(-m), Int((d*csc(a + b*x))**(-m + n)*(d*sec(a + b*x))**m*ActivateTrig(u), x), x)
    rule4713 = ReplacementRule(pattern4713, replacement4713)
    pattern4714 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644, cons18)
    def replacement4714(a, n, m, b, u, x, d, c):
        rubi.append(4714)
        return Dist((c*cot(a + b*x))**m*(d*csc(a + b*x))**(-m)*(d*sec(a + b*x))**m, Int((d*csc(a + b*x))**m*(d*sec(a + b*x))**(-m + n)*ActivateTrig(u), x), x)
    rule4714 = ReplacementRule(pattern4714, replacement4714)
    pattern4715 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('d', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1644, cons18)
    def replacement4715(a, n, m, b, u, x, d, c):
        rubi.append(4715)
        return Dist((c*cot(a + b*x))**m*(d*csc(a + b*x))**(-m)*(d*sec(a + b*x))**m, Int((d*csc(a + b*x))**(m + n)*(d*sec(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4715 = ReplacementRule(pattern4715, replacement4715)
    pattern4716 = Pattern(Integral(u_*(WC('c', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1644)
    def replacement4716(a, m, b, u, x, c):
        rubi.append(4716)
        return Dist((c*sin(a + b*x))**m*(c*csc(a + b*x))**m, Int((c*csc(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4716 = ReplacementRule(pattern4716, replacement4716)
    pattern4717 = Pattern(Integral(u_*(WC('c', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1644)
    def replacement4717(a, m, b, u, x, c):
        rubi.append(4717)
        return Dist((c*cos(a + b*x))**m*(c*sec(a + b*x))**m, Int((c*sec(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4717 = ReplacementRule(pattern4717, replacement4717)
    pattern4718 = Pattern(Integral(u_*(WC('c', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1644)
    def replacement4718(a, m, b, u, x, c):
        rubi.append(4718)
        return Dist((c*tan(a + b*x))**m*(c*csc(a + b*x))**m*(c*sec(a + b*x))**(-m), Int((c*csc(a + b*x))**(-m)*(c*sec(a + b*x))**m*ActivateTrig(u), x), x)
    rule4718 = ReplacementRule(pattern4718, replacement4718)
    pattern4719 = Pattern(Integral(u_*(WC('c', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons7, cons21, cons18, cons1644)
    def replacement4719(a, m, b, u, x, c):
        rubi.append(4719)
        return Dist((c*cot(a + b*x))**m*(c*csc(a + b*x))**(-m)*(c*sec(a + b*x))**m, Int((c*csc(a + b*x))**m*(c*sec(a + b*x))**(-m)*ActivateTrig(u), x), x)
    rule4719 = ReplacementRule(pattern4719, replacement4719)
    pattern4720 = Pattern(Integral(u_*(WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1644)
    def replacement4720(a, n, B, b, A, u, x, c):
        rubi.append(4720)
        return Dist(c, Int((c*sec(a + b*x))**(n + S(-1))*(A*sec(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4720 = ReplacementRule(pattern4720, replacement4720)
    pattern4721 = Pattern(Integral(u_*(WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons1644)
    def replacement4721(a, n, B, b, A, u, x, c):
        rubi.append(4721)
        return Dist(c, Int((c*csc(a + b*x))**(n + S(-1))*(A*csc(a + b*x) + B)*ActivateTrig(u), x), x)
    rule4721 = ReplacementRule(pattern4721, replacement4721)
    pattern4722 = Pattern(Integral(u_*(A_ + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1644)
    def replacement4722(a, B, b, A, u, x):
        rubi.append(4722)
        return Int((A*sec(a + b*x) + B)*ActivateTrig(u)/sec(a + b*x), x)
    rule4722 = ReplacementRule(pattern4722, replacement4722)
    pattern4723 = Pattern(Integral(u_*(A_ + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons34, cons35, cons1644)
    def replacement4723(a, B, b, A, u, x):
        rubi.append(4723)
        return Int((A*csc(a + b*x) + B)*ActivateTrig(u)/csc(a + b*x), x)
    rule4723 = ReplacementRule(pattern4723, replacement4723)
    pattern4724 = Pattern(Integral((WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1644)
    def replacement4724(a, n, B, b, A, u, x, C, c):
        rubi.append(4724)
        return Dist(c**S(2), Int((c*sec(a + b*x))**(n + S(-2))*(A*sec(a + b*x)**S(2) + B*sec(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4724 = ReplacementRule(pattern4724, replacement4724)
    pattern4725 = Pattern(Integral((WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons36, cons4, cons1644)
    def replacement4725(a, n, B, b, A, u, x, C, c):
        rubi.append(4725)
        return Dist(c**S(2), Int((c*csc(a + b*x))**(n + S(-2))*(A*csc(a + b*x)**S(2) + B*csc(a + b*x) + C)*ActivateTrig(u), x), x)
    rule4725 = ReplacementRule(pattern4725, replacement4725)
    pattern4726 = Pattern(Integral((WC('c', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1644)
    def replacement4726(a, n, A, b, C, u, x, c):
        rubi.append(4726)
        return Dist(c**S(2), Int((c*sec(a + b*x))**(n + S(-2))*(A*sec(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4726 = ReplacementRule(pattern4726, replacement4726)
    pattern4727 = Pattern(Integral((WC('c', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons36, cons4, cons1644)
    def replacement4727(a, n, A, b, C, u, x, c):
        rubi.append(4727)
        return Dist(c**S(2), Int((c*csc(a + b*x))**(n + S(-2))*(A*csc(a + b*x)**S(2) + C)*ActivateTrig(u), x), x)
    rule4727 = ReplacementRule(pattern4727, replacement4727)
    pattern4728 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1644)
    def replacement4728(a, B, b, A, u, x, C):
        rubi.append(4728)
        return Int((A*sec(a + b*x)**S(2) + B*sec(a + b*x) + C)*ActivateTrig(u)/sec(a + b*x)**S(2), x)
    rule4728 = ReplacementRule(pattern4728, replacement4728)
    pattern4729 = Pattern(Integral(u_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))) + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons35, cons36, cons1644)
    def replacement4729(a, B, A, b, C, u, x):
        rubi.append(4729)
        return Int((A*csc(a + b*x)**S(2) + B*csc(a + b*x) + C)*ActivateTrig(u)/csc(a + b*x)**S(2), x)
    rule4729 = ReplacementRule(pattern4729, replacement4729)
    pattern4730 = Pattern(Integral(u_*(A_ + WC('C', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1644)
    def replacement4730(a, A, b, C, u, x):
        rubi.append(4730)
        return Int((A*sec(a + b*x)**S(2) + C)*ActivateTrig(u)/sec(a + b*x)**S(2), x)
    rule4730 = ReplacementRule(pattern4730, replacement4730)
    pattern4731 = Pattern(Integral(u_*(A_ + WC('C', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)), x_), cons2, cons3, cons34, cons36, cons1644)
    def replacement4731(a, A, b, C, u, x):
        rubi.append(4731)
        return Int((A*csc(a + b*x)**S(2) + C)*ActivateTrig(u)/csc(a + b*x)**S(2), x)
    rule4731 = ReplacementRule(pattern4731, replacement4731)
    pattern4732 = Pattern(Integral(u_*(WC('A', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4732(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4732)
        return Int((A + B*sec(a + b*x) + C*sec(a + b*x)**S(2))*ActivateTrig(u)*sec(a + b*x)**n, x)
    rule4732 = ReplacementRule(pattern4732, replacement4732)
    pattern4733 = Pattern(Integral(u_*(WC('A', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)) + WC('B', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**n1_ + WC('C', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**n2_), x_), cons2, cons3, cons34, cons35, cons36, cons4, cons1640, cons1641)
    def replacement4733(n1, a, n, B, n2, b, A, u, x, C):
        rubi.append(4733)
        return Int((A + B*csc(a + b*x) + C*csc(a + b*x)**S(2))*ActivateTrig(u)*csc(a + b*x)**n, x)
    rule4733 = ReplacementRule(pattern4733, replacement4733)
    pattern4734 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1645)
    def replacement4734(a, b, x, d, c):
        rubi.append(4734)
        return Simp(sin(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) - Simp(sin(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
    rule4734 = ReplacementRule(pattern4734, replacement4734)
    pattern4735 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1645)
    def replacement4735(a, b, x, d, c):
        rubi.append(4735)
        return Simp(sin(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) + Simp(sin(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
    rule4735 = ReplacementRule(pattern4735, replacement4735)
    pattern4736 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1645)
    def replacement4736(a, b, x, d, c):
        rubi.append(4736)
        return -Simp(cos(a - c + x*(b - d))/(S(2)*b - S(2)*d), x) - Simp(cos(a + c + x*(b + d))/(S(2)*b + S(2)*d), x)
    rule4736 = ReplacementRule(pattern4736, replacement4736)
    pattern4737 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons1400)
    def replacement4737(a, g, p, b, x, d, c):
        rubi.append(4737)
        return Dist(S(1)/2, Int((g*sin(c + d*x))**p, x), x) + Dist(S(1)/2, Int((g*sin(c + d*x))**p*cos(c + d*x), x), x)
    rule4737 = ReplacementRule(pattern4737, replacement4737)
    pattern4738 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons1400)
    def replacement4738(a, g, p, b, x, d, c):
        rubi.append(4738)
        return Dist(S(1)/2, Int((g*sin(c + d*x))**p, x), x) - Dist(S(1)/2, Int((g*sin(c + d*x))**p*cos(c + d*x), x), x)
    rule4738 = ReplacementRule(pattern4738, replacement4738)
    pattern4739 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons25, cons1646, cons38)
    def replacement4739(e, a, m, p, b, x, d, c):
        rubi.append(4739)
        return Dist(S(2)**p*e**(-p), Int((e*cos(a + b*x))**(m + p)*sin(a + b*x)**p, x), x)
    rule4739 = ReplacementRule(pattern4739, replacement4739)
    pattern4740 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons4, cons25, cons1646, cons38)
    def replacement4740(a, n, p, b, f, x, d, c):
        rubi.append(4740)
        return Dist(S(2)**p*f**(-p), Int((f*sin(a + b*x))**(n + p)*cos(a + b*x)**p, x), x)
    rule4740 = ReplacementRule(pattern4740, replacement4740)
    pattern4741 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons5, cons25, cons1646, cons147, cons343)
    def replacement4741(e, a, m, g, p, b, x, d, c):
        rubi.append(4741)
        return Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4741 = ReplacementRule(pattern4741, replacement4741)
    pattern4742 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons5, cons25, cons1646, cons147, cons343)
    def replacement4742(e, a, m, g, p, b, x, d, c):
        rubi.append(4742)
        return -Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4742 = ReplacementRule(pattern4742, replacement4742)
    pattern4743 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons5, cons25, cons1646, cons147, cons240)
    def replacement4743(e, a, m, g, p, b, x, d, c):
        rubi.append(4743)
        return -Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(b*g*m), x)
    rule4743 = ReplacementRule(pattern4743, replacement4743)
    pattern4744 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons5, cons25, cons1646, cons147, cons240)
    def replacement4744(e, a, m, g, p, b, x, d, c):
        rubi.append(4744)
        return Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(b*g*m), x)
    rule4744 = ReplacementRule(pattern4744, replacement4744)
    pattern4745 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons25, cons1646, cons147, cons244, cons1325, cons137, cons1647, cons1280)
    def replacement4745(e, a, m, g, p, b, x, d, c):
        rubi.append(4745)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(p + S(1))), Int((e*cos(a + b*x))**(m + S(-4))*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4745 = ReplacementRule(pattern4745, replacement4745)
    pattern4746 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons25, cons1646, cons147, cons244, cons1325, cons137, cons1647, cons1280)
    def replacement4746(e, a, m, g, p, b, x, d, c):
        rubi.append(4746)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(p + S(1))), Int((e*sin(a + b*x))**(m + S(-4))*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4746 = ReplacementRule(pattern4746, replacement4746)
    pattern4747 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons25, cons1646, cons147, cons244, cons166, cons137, cons319, cons1648, cons1280)
    def replacement4747(e, a, m, g, p, b, x, d, c):
        rubi.append(4747)
        return Dist(e**S(2)*(m + S(2)*p + S(2))/(S(4)*g**S(2)*(p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4747 = ReplacementRule(pattern4747, replacement4747)
    pattern4748 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons25, cons1646, cons147, cons244, cons166, cons137, cons319, cons1648, cons1280)
    def replacement4748(e, a, m, g, p, b, x, d, c):
        rubi.append(4748)
        return Dist(e**S(2)*(m + S(2)*p + S(2))/(S(4)*g**S(2)*(p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(p + S(1))), x)
    rule4748 = ReplacementRule(pattern4748, replacement4748)
    pattern4749 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons5, cons25, cons1646, cons147, cons31, cons166, cons249, cons1280)
    def replacement4749(e, a, m, g, p, b, x, d, c):
        rubi.append(4749)
        return Dist(e**S(2)*(m + p + S(-1))/(m + S(2)*p), Int((e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**p, x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)), x)
    rule4749 = ReplacementRule(pattern4749, replacement4749)
    pattern4750 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons5, cons25, cons1646, cons147, cons31, cons166, cons249, cons1280)
    def replacement4750(e, a, m, g, p, b, x, d, c):
        rubi.append(4750)
        return Dist(e**S(2)*(m + p + S(-1))/(m + S(2)*p), Int((e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**p, x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + S(2)*p)), x)
    rule4750 = ReplacementRule(pattern4750, replacement4750)
    pattern4751 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons5, cons25, cons1646, cons147, cons31, cons94, cons319, cons253, cons1280)
    def replacement4751(e, a, m, g, p, b, x, d, c):
        rubi.append(4751)
        return Dist((m + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(2))*(g*sin(c + d*x))**p, x), x) - Simp((e*cos(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))), x)
    rule4751 = ReplacementRule(pattern4751, replacement4751)
    pattern4752 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons5, cons25, cons1646, cons147, cons31, cons94, cons319, cons253, cons1280)
    def replacement4752(e, a, m, g, p, b, x, d, c):
        rubi.append(4752)
        return Dist((m + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(2))*(g*sin(c + d*x))**p, x), x) + Simp((e*sin(a + b*x))**m*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(m + p + S(1))), x)
    rule4752 = ReplacementRule(pattern4752, replacement4752)
    pattern4753 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons147, cons13, cons163, cons246)
    def replacement4753(a, g, p, b, x, d, c):
        rubi.append(4753)
        return Dist(S(2)*g*p/(S(2)*p + S(1)), Int((g*sin(c + d*x))**(p + S(-1))*sin(a + b*x), x), x) + Simp(S(2)*(g*sin(c + d*x))**p*sin(a + b*x)/(d*(S(2)*p + S(1))), x)
    rule4753 = ReplacementRule(pattern4753, replacement4753)
    pattern4754 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons147, cons13, cons163, cons246)
    def replacement4754(a, g, p, b, x, d, c):
        rubi.append(4754)
        return Dist(S(2)*g*p/(S(2)*p + S(1)), Int((g*sin(c + d*x))**(p + S(-1))*cos(a + b*x), x), x) + Simp(-S(2)*(g*sin(c + d*x))**p*cos(a + b*x)/(d*(S(2)*p + S(1))), x)
    rule4754 = ReplacementRule(pattern4754, replacement4754)
    pattern4755 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons147, cons13, cons137, cons246)
    def replacement4755(a, g, p, b, x, d, c):
        rubi.append(4755)
        return Dist((S(2)*p + S(3))/(S(2)*g*(p + S(1))), Int((g*sin(c + d*x))**(p + S(1))*sin(a + b*x), x), x) + Simp((g*sin(c + d*x))**(p + S(1))*cos(a + b*x)/(S(2)*b*g*(p + S(1))), x)
    rule4755 = ReplacementRule(pattern4755, replacement4755)
    pattern4756 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons147, cons13, cons137, cons246)
    def replacement4756(a, g, p, b, x, d, c):
        rubi.append(4756)
        return Dist((S(2)*p + S(3))/(S(2)*g*(p + S(1))), Int((g*sin(c + d*x))**(p + S(1))*cos(a + b*x), x), x) - Simp((g*sin(c + d*x))**(p + S(1))*sin(a + b*x)/(S(2)*b*g*(p + S(1))), x)
    rule4756 = ReplacementRule(pattern4756, replacement4756)
    pattern4757 = Pattern(Integral(cos(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons25, cons1646)
    def replacement4757(a, b, x, d, c):
        rubi.append(4757)
        return Simp(log(sin(a + b*x) + sqrt(sin(c + d*x)) + cos(a + b*x))/d, x) - Simp(-asin(sin(a + b*x) - cos(a + b*x))/d, x)
    rule4757 = ReplacementRule(pattern4757, replacement4757)
    pattern4758 = Pattern(Integral(sin(x_*WC('b', S(1)) + WC('a', S(0)))/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons25, cons1646)
    def replacement4758(a, b, x, d, c):
        rubi.append(4758)
        return -Simp(log(sin(a + b*x) + sqrt(sin(c + d*x)) + cos(a + b*x))/d, x) - Simp(-asin(sin(a + b*x) - cos(a + b*x))/d, x)
    rule4758 = ReplacementRule(pattern4758, replacement4758)
    pattern4759 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons5, cons25, cons1646, cons147, cons246)
    def replacement4759(a, g, p, b, x, d, c):
        rubi.append(4759)
        return Dist(S(2)*g, Int((g*sin(c + d*x))**(p + S(-1))*sin(a + b*x), x), x)
    rule4759 = ReplacementRule(pattern4759, replacement4759)
    pattern4760 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_/sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons208, cons5, cons25, cons1646, cons147, cons246)
    def replacement4760(a, g, p, b, x, d, c):
        rubi.append(4760)
        return Dist(S(2)*g, Int((g*sin(c + d*x))**(p + S(-1))*cos(a + b*x), x), x)
    rule4760 = ReplacementRule(pattern4760, replacement4760)
    pattern4761 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons208, cons21, cons5, cons25, cons1646, cons147)
    def replacement4761(e, a, m, g, p, b, x, d, c):
        rubi.append(4761)
        return Dist((e*cos(a + b*x))**(-p)*(g*sin(c + d*x))**p*sin(a + b*x)**(-p), Int((e*cos(a + b*x))**(m + p)*sin(a + b*x)**p, x), x)
    rule4761 = ReplacementRule(pattern4761, replacement4761)
    pattern4762 = Pattern(Integral((WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons125, cons208, cons4, cons5, cons25, cons1646, cons147)
    def replacement4762(a, n, g, p, b, f, x, d, c):
        rubi.append(4762)
        return Dist((f*sin(a + b*x))**(-p)*(g*sin(c + d*x))**p*cos(a + b*x)**(-p), Int((f*sin(a + b*x))**(n + p)*cos(a + b*x)**p, x), x)
    rule4762 = ReplacementRule(pattern4762, replacement4762)
    pattern4763 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*cos(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons208, cons25, cons1646, cons1400)
    def replacement4763(a, g, p, b, x, d, c):
        rubi.append(4763)
        return Dist(S(1)/4, Int((g*sin(c + d*x))**p, x), x) - Dist(S(1)/4, Int((g*sin(c + d*x))**p*cos(c + d*x)**S(2), x), x)
    rule4763 = ReplacementRule(pattern4763, replacement4763)
    pattern4764 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons25, cons1646, cons38)
    def replacement4764(e, a, n, m, p, b, f, x, d, c):
        rubi.append(4764)
        return Dist(S(2)**p*e**(-p)*f**(-p), Int((e*cos(a + b*x))**(m + p)*(f*sin(a + b*x))**(n + p), x), x)
    rule4764 = ReplacementRule(pattern4764, replacement4764)
    pattern4765 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons25, cons1646, cons147, cons1270)
    def replacement4765(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4765)
        return Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
    rule4765 = ReplacementRule(pattern4765, replacement4765)
    pattern4766 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons25, cons1646, cons147, cons1270)
    def replacement4766(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4766)
        return -Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
    rule4766 = ReplacementRule(pattern4766, replacement4766)
    pattern4767 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons25, cons1646, cons147, cons1649, cons253)
    def replacement4767(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4767)
        return -Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
    rule4767 = ReplacementRule(pattern4767, replacement4767)
    pattern4768 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons25, cons1646, cons147, cons244, cons1650, cons137, cons1651, cons170)
    def replacement4768(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4768)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-4))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp(e**S(2)*(e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
    rule4768 = ReplacementRule(pattern4768, replacement4768)
    pattern4769 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons25, cons1646, cons147, cons244, cons1650, cons137, cons1651, cons170)
    def replacement4769(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4769)
        return Dist(e**S(4)*(m + p + S(-1))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-4))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp(e**S(2)*(e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
    rule4769 = ReplacementRule(pattern4769, replacement4769)
    pattern4770 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons25, cons1646, cons147, cons244, cons166, cons137, cons1652, cons1651, cons170, cons1653)
    def replacement4770(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4770)
        return Dist(e**S(2)*(m + n + S(2)*p + S(2))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) + Simp((e*cos(a + b*x))**m*(f*sin(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
    rule4770 = ReplacementRule(pattern4770, replacement4770)
    pattern4771 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons25, cons1646, cons147, cons244, cons166, cons137, cons1652, cons1651, cons170, cons1653)
    def replacement4771(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4771)
        return Dist(e**S(2)*(m + n + S(2)*p + S(2))/(S(4)*g**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(2)), x), x) - Simp((e*sin(a + b*x))**m*(f*cos(a + b*x))**n*(g*sin(c + d*x))**(p + S(1))/(S(2)*b*g*(n + p + S(1))), x)
    rule4771 = ReplacementRule(pattern4771, replacement4771)
    pattern4772 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons25, cons1646, cons147, cons93, cons166, cons89, cons1651, cons170)
    def replacement4772(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4772)
        return Dist(e**S(2)*(m + p + S(-1))/(f**S(2)*(n + p + S(1))), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**(n + S(2))*(g*sin(c + d*x))**p, x), x) + Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
    rule4772 = ReplacementRule(pattern4772, replacement4772)
    pattern4773 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons5, cons25, cons1646, cons147, cons93, cons166, cons89, cons1651, cons170)
    def replacement4773(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4773)
        return Dist(e**S(2)*(m + p + S(-1))/(f**S(2)*(n + p + S(1))), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**(n + S(2))*(g*sin(c + d*x))**p, x), x) - Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(n + p + S(1))), x)
    rule4773 = ReplacementRule(pattern4773, replacement4773)
    pattern4774 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons25, cons1646, cons147, cons31, cons166, cons1654, cons170)
    def replacement4774(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4774)
        return Dist(e**S(2)*(m + p + S(-1))/(m + n + S(2)*p), Int((e*cos(a + b*x))**(m + S(-2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**p, x), x) + Simp(e*(e*cos(a + b*x))**(m + S(-1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(m + n + S(2)*p)), x)
    rule4774 = ReplacementRule(pattern4774, replacement4774)
    pattern4775 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons25, cons1646, cons147, cons31, cons166, cons1654, cons170)
    def replacement4775(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4775)
        return Dist(e**S(2)*(m + p + S(-1))/(m + n + S(2)*p), Int((e*sin(a + b*x))**(m + S(-2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**p, x), x) - Simp(e*(e*sin(a + b*x))**(m + S(-1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*f*(m + n + S(2)*p)), x)
    rule4775 = ReplacementRule(pattern4775, replacement4775)
    pattern4776 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons25, cons1646, cons147, cons162, cons94, cons88, cons163, cons1654, cons170)
    def replacement4776(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4776)
        return Dist(S(2)*f*g*(n + p + S(-1))/(e*(m + n + S(2)*p)), Int((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(-1)), x), x) - Simp(f*(e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**p/(b*e*(m + n + S(2)*p)), x)
    rule4776 = ReplacementRule(pattern4776, replacement4776)
    pattern4777 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons25, cons1646, cons147, cons162, cons94, cons88, cons163, cons1654, cons170)
    def replacement4777(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4777)
        return Dist(S(2)*f*g*(n + p + S(-1))/(e*(m + n + S(2)*p)), Int((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(-1)), x), x) + Simp(f*(e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**p/(b*e*(m + n + S(2)*p)), x)
    rule4777 = ReplacementRule(pattern4777, replacement4777)
    pattern4778 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons25, cons1646, cons147, cons162, cons94, cons88, cons137, cons1652, cons253, cons170)
    def replacement4778(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4778)
        return Dist(f*(m + n + S(2)*p + S(2))/(S(2)*e*g*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(1)), x), x) - Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
    rule4778 = ReplacementRule(pattern4778, replacement4778)
    pattern4779 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons25, cons1646, cons147, cons162, cons94, cons88, cons137, cons1652, cons253, cons170)
    def replacement4779(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4779)
        return Dist(f*(m + n + S(2)*p + S(2))/(S(2)*e*g*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(-1))*(g*sin(c + d*x))**(p + S(1)), x), x) + Simp((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
    rule4779 = ReplacementRule(pattern4779, replacement4779)
    pattern4780 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons25, cons1646, cons147, cons31, cons94, cons1652, cons253, cons170)
    def replacement4780(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4780)
        return Dist((m + n + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*cos(a + b*x))**(m + S(2))*(f*sin(a + b*x))**n*(g*sin(c + d*x))**p, x), x) - Simp((e*cos(a + b*x))**(m + S(1))*(f*sin(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
    rule4780 = ReplacementRule(pattern4780, replacement4780)
    pattern4781 = Pattern(Integral((WC('e', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**m_*(WC('f', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons5, cons25, cons1646, cons147, cons31, cons94, cons1652, cons253, cons170)
    def replacement4781(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4781)
        return Dist((m + n + S(2)*p + S(2))/(e**S(2)*(m + p + S(1))), Int((e*sin(a + b*x))**(m + S(2))*(f*cos(a + b*x))**n*(g*sin(c + d*x))**p, x), x) + Simp((e*sin(a + b*x))**(m + S(1))*(f*cos(a + b*x))**(n + S(1))*(g*sin(c + d*x))**p/(b*e*f*(m + p + S(1))), x)
    rule4781 = ReplacementRule(pattern4781, replacement4781)
    pattern4782 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*(WC('f', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons25, cons1646, cons147)
    def replacement4782(e, a, n, m, g, p, b, f, x, d, c):
        rubi.append(4782)
        return Dist((e*cos(a + b*x))**(-p)*(f*sin(a + b*x))**(-p)*(g*sin(c + d*x))**p, Int((e*cos(a + b*x))**(m + p)*(f*sin(a + b*x))**(n + p), x), x)
    rule4782 = ReplacementRule(pattern4782, replacement4782)
    pattern4783 = Pattern(Integral((WC('e', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons21, cons25, cons1655)
    def replacement4783(e, a, m, b, x, d, c):
        rubi.append(4783)
        return -Simp((e*cos(a + b*x))**(m + S(1))*(m + S(2))*cos((a + b*x)*(m + S(1)))/(d*e*(m + S(1))), x)
    rule4783 = ReplacementRule(pattern4783, replacement4783)
    pattern4784 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_)**p_, x_), cons2, cons3, cons7, cons27, cons1656, cons85, cons128)
    def replacement4784(a, n, p, F, b, x, d, c):
        rubi.append(4784)
        return Int((a + b*F(c + d*x)**n)**p, x)
    rule4784 = ReplacementRule(pattern4784, replacement4784)
    pattern4785 = Pattern(Integral(S(1)/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons7, cons27, cons1656, cons1471, cons744)
    def replacement4785(a, n, F, b, x, d, c):
        rubi.append(4785)
        return Dist(S(2)/(a*n), Sum_doit(Int(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*F(c + d*x)**S(2)/Rt(-a/b, n/S(2))), x), List(k, S(1), n/S(2))), x)
    rule4785 = ReplacementRule(pattern4785, replacement4785)
    pattern4786 = Pattern(Integral(S(1)/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons7, cons27, cons1656, cons1474, cons744)
    def replacement4786(a, n, F, b, x, d, c):
        rubi.append(4786)
        return Int(ExpandTrig(S(1)/(a + b*F(c + d*x)**n), x), x)
    rule4786 = ReplacementRule(pattern4786, replacement4786)
    pattern4787 = Pattern(Integral(G_**(x_*WC('d', S(1)) + WC('c', S(0)))/(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + a_), x_), cons2, cons3, cons7, cons27, cons21, cons1657, cons85, cons744)
    def replacement4787(G, a, n, m, b, F, x, d, c):
        rubi.append(4787)
        return Int(ExpandTrig(G(c + d*x)**m, S(1)/(a + b*F(c + d*x)**n), x), x)
    rule4787 = ReplacementRule(pattern4787, replacement4787)
    def With4788(a, n, p, F, x, d, c):
        v = ActivateTrig(F(c + d*x))
        rubi.append(4788)
        return Dist(a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*NonfreeFactors(v, x)**(-p*FracPart(n)), Int(NonfreeFactors(v, x)**(n*p), x), x)
    pattern4788 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_, x_), cons2, cons7, cons27, cons4, cons5, cons1656, cons23, cons38)
    rule4788 = ReplacementRule(pattern4788, With4788)
    def With4789(a, n, p, b, F, x, d, c):
        v = ActivateTrig(F(c + d*x))
        rubi.append(4789)
        return Dist(a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n)), Int((b*v)**(n*p), x), x)
    pattern4789 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons1656, cons23, cons147)
    rule4789 = ReplacementRule(pattern4789, With4789)
    def With4790(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4790 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1658, CustomConstraint(With4790))
    def replacement4790(a, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4790)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4790 = ReplacementRule(pattern4790, replacement4790)
    def With4791(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4791 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1659, CustomConstraint(With4791))
    def replacement4791(a, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4791)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4791 = ReplacementRule(pattern4791, replacement4791)
    def With4792(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4792 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1660, CustomConstraint(With4792))
    def replacement4792(a, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4792)
        return Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4792 = ReplacementRule(pattern4792, replacement4792)
    def With4793(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4793 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1661, CustomConstraint(With4793))
    def replacement4793(a, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4793)
        return -Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4793 = ReplacementRule(pattern4793, replacement4793)
    def With4794(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if FunctionOfQ(tan(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4794 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1239, cons1662, CustomConstraint(With4794))
    def replacement4794(a, F, b, u, x, c):

        d = FreeFactors(tan(c*(a + b*x)), x)
        rubi.append(4794)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)
    rule4794 = ReplacementRule(pattern4794, replacement4794)
    def With4795(a, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if FunctionOfQ(tan(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4795 = Pattern(Integral(u_/cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), cons2, cons3, cons7, cons1239, CustomConstraint(With4795))
    def replacement4795(a, b, u, x, c):

        d = FreeFactors(tan(c*(a + b*x)), x)
        rubi.append(4795)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)
    rule4795 = ReplacementRule(pattern4795, replacement4795)
    def With4796(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cot(c*(a + b*x)), x)
        if FunctionOfQ(cot(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4796 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1239, cons1663, CustomConstraint(With4796))
    def replacement4796(a, F, b, u, x, c):

        d = FreeFactors(cot(c*(a + b*x)), x)
        rubi.append(4796)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cot(c*(a + b*x))/d, u, x), x), x, cot(c*(a + b*x))/d), x)
    rule4796 = ReplacementRule(pattern4796, replacement4796)
    def With4797(a, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cot(c*(a + b*x)), x)
        if FunctionOfQ(cot(c*(a + b*x))/d, u, x, True):
            return True
        return False
    pattern4797 = Pattern(Integral(u_/sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**S(2), x_), cons2, cons3, cons7, cons1239, CustomConstraint(With4797))
    def replacement4797(a, b, u, x, c):

        d = FreeFactors(cot(c*(a + b*x)), x)
        rubi.append(4797)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cot(c*(a + b*x))/d, u, x), x), x, cot(c*(a + b*x))/d), x)
    rule4797 = ReplacementRule(pattern4797, replacement4797)
    def With4798(a, n, b, F, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(tan(c*(a + b*x)), x)
        if And(FunctionOfQ(tan(c*(a + b*x))/d, u, x, True), TryPureTanSubst(ActivateTrig(u)*cot(c*(a + b*x))**n, x)):
            return True
        return False
    pattern4798 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons85, cons1660, CustomConstraint(With4798))
    def replacement4798(a, n, b, F, u, x, c):

        d = FreeFactors(tan(c*(a + b*x)), x)
        rubi.append(4798)
        return Dist(d**(-n + S(1))/(b*c), Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), tan(c*(a + b*x))/d, u, x), x), x, tan(c*(a + b*x))/d), x)
    rule4798 = ReplacementRule(pattern4798, replacement4798)
    def With4799(a, n, b, F, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cot(c*(a + b*x)), x)
        if And(FunctionOfQ(cot(c*(a + b*x))/d, u, x, True), TryPureTanSubst(ActivateTrig(u)*tan(c*(a + b*x))**n, x)):
            return True
        return False
    pattern4799 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons85, cons1661, CustomConstraint(With4799))
    def replacement4799(a, n, b, F, u, x, c):

        d = FreeFactors(cot(c*(a + b*x)), x)
        rubi.append(4799)
        return -Dist(d**(-n + S(1))/(b*c), Subst(Int(SubstFor(x**(-n)/(d**S(2)*x**S(2) + S(1)), cot(c*(a + b*x))/d, u, x), x), x, cot(c*(a + b*x))/d), x)
    rule4799 = ReplacementRule(pattern4799, replacement4799)
    def With4800(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = FunctionOfTrig(u, x)
        d = FreeFactors(cot(v), x)
        if And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(cot(v), x), u, x, True), TryPureTanSubst(ActivateTrig(u), x)):
            return True
        return False
    pattern4800 = Pattern(Integral(u_, x_), CustomConstraint(With4800))
    def replacement4800(u, x):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(cot(v), x)
        rubi.append(4800)
        return Simp(With(List(Set(d, FreeFactors(cot(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1)/(d**S(2)*x**S(2) + S(1)), cot(v)/d, u, x), x), x, cot(v)/d), x)), x)
    rule4800 = ReplacementRule(pattern4800, replacement4800)
    def With4801(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        if And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(tan(v), x), u, x, True), TryPureTanSubst(ActivateTrig(u), x)):
            return True
        return False
    pattern4801 = Pattern(Integral(u_, x_), CustomConstraint(With4801))
    def replacement4801(u, x):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        rubi.append(4801)
        return Simp(With(List(Set(d, FreeFactors(tan(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1)/(d**S(2)*x**S(2) + S(1)), tan(v)/d, u, x), x), x, tan(v)/d), x)), x)
    rule4801 = ReplacementRule(pattern4801, replacement4801)
    pattern4802 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1664, cons1665, cons555)
    def replacement4802(G, a, q, p, b, F, x, d, c):
        rubi.append(4802)
        return Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q), x), x)
    rule4802 = ReplacementRule(pattern4802, replacement4802)
    pattern4803 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*H_**(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons1664, cons1665, cons1666, cons628)
    def replacement4803(H, e, G, a, q, p, b, r, F, f, x, d, c):
        rubi.append(4803)
        return Int(ExpandTrigReduce(ActivateTrig(F(a + b*x)**p*G(c + d*x)**q*H(e + f*x)**r), x), x)
    rule4803 = ReplacementRule(pattern4803, replacement4803)
    def With4804(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4804 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1658, CustomConstraint(With4804))
    def replacement4804(a, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4804)
        return Dist(d/(b*c), Subst(Int(SubstFor(S(1), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4804 = ReplacementRule(pattern4804, replacement4804)
    def With4805(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4805 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1659, CustomConstraint(With4805))
    def replacement4805(a, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4805)
        return -Dist(d/(b*c), Subst(Int(SubstFor(S(1), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4805 = ReplacementRule(pattern4805, replacement4805)
    def With4806(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4806 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1660, CustomConstraint(With4806))
    def replacement4806(a, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4806)
        return Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4806 = ReplacementRule(pattern4806, replacement4806)
    def With4807(a, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4807 = Pattern(Integral(F_*u_*(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)), x_), cons2, cons3, cons7, cons1661, CustomConstraint(With4807))
    def replacement4807(a, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4807)
        return -Dist(S(1)/(b*c), Subst(Int(SubstFor(S(1)/x, cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4807 = ReplacementRule(pattern4807, replacement4807)
    def With4808(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4808 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1658, CustomConstraint(With4808))
    def replacement4808(a, n, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4808)
        return Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4808 = ReplacementRule(pattern4808, replacement4808)
    def With4809(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4809 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1662, CustomConstraint(With4809))
    def replacement4809(a, n, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4809)
        return Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4809 = ReplacementRule(pattern4809, replacement4809)
    def With4810(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4810 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1659, CustomConstraint(With4810))
    def replacement4810(a, n, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4810)
        return -Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4810 = ReplacementRule(pattern4810, replacement4810)
    def With4811(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4811 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1663, CustomConstraint(With4811))
    def replacement4811(a, n, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4811)
        return -Dist(d/(b*c), Subst(Int(SubstFor((-d**S(2)*x**S(2) + S(1))**(-n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4811 = ReplacementRule(pattern4811, replacement4811)
    def With4812(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4812 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1660, CustomConstraint(With4812))
    def replacement4812(a, n, F, b, u, x, c):

        d = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4812)
        return Dist(d**(-n + S(1))/(b*c), Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), sin(c*(a + b*x))/d, u, x), x), x, sin(c*(a + b*x))/d), x)
    rule4812 = ReplacementRule(pattern4812, replacement4812)
    def With4813(a, n, F, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        d = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/d, u, x):
            return True
        return False
    pattern4813 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*u_, x_), cons2, cons3, cons7, cons1474, cons1239, cons1661, CustomConstraint(With4813))
    def replacement4813(a, n, F, b, u, x, c):

        d = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4813)
        return -Dist(d**(-n + S(1))/(b*c), Subst(Int(SubstFor(x**(-n)*(-d**S(2)*x**S(2) + S(1))**(n/S(2) + S(-1)/2), cos(c*(a + b*x))/d, u, x), x), x, cos(c*(a + b*x))/d), x)
    rule4813 = ReplacementRule(pattern4813, replacement4813)
    def With4814(a, n, v, b, F, u, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        e = FreeFactors(sin(c*(a + b*x)), x)
        if FunctionOfQ(sin(c*(a + b*x))/e, u, x):
            return True
        return False
    pattern4814 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), cons2, cons3, cons7, cons27, cons10, cons1474, cons1239, cons1658, CustomConstraint(With4814))
    def replacement4814(a, n, v, b, F, u, x, d, c):

        e = FreeFactors(sin(c*(a + b*x)), x)
        rubi.append(4814)
        return Dist(d, Int(ActivateTrig(u)*cos(c*(a + b*x))**n, x), x) + Int(ActivateTrig(u*v), x)
    rule4814 = ReplacementRule(pattern4814, replacement4814)
    def With4815(a, n, v, b, F, u, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        e = FreeFactors(cos(c*(a + b*x)), x)
        if FunctionOfQ(cos(c*(a + b*x))/e, u, x):
            return True
        return False
    pattern4815 = Pattern(Integral(u_*(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*WC('d', S(1)) + v_), x_), cons2, cons3, cons7, cons27, cons10, cons1474, cons1239, cons1659, CustomConstraint(With4815))
    def replacement4815(a, n, v, b, F, u, x, d, c):

        e = FreeFactors(cos(c*(a + b*x)), x)
        rubi.append(4815)
        return Dist(d, Int(ActivateTrig(u)*sin(c*(a + b*x))**n, x), x) + Int(ActivateTrig(u*v), x)
    rule4815 = ReplacementRule(pattern4815, replacement4815)
    def With4816(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = FunctionOfTrig(u, x)
        d = FreeFactors(sin(v), x)
        if And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(sin(v), x), u/cos(v), x)):
            return True
        return False
    pattern4816 = Pattern(Integral(u_, x_), CustomConstraint(With4816))
    def replacement4816(u, x):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(sin(v), x)
        rubi.append(4816)
        return Simp(With(List(Set(d, FreeFactors(sin(v), x))), Dist(d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), sin(v)/d, u/cos(v), x), x), x, sin(v)/d), x)), x)
    rule4816 = ReplacementRule(pattern4816, replacement4816)
    def With4817(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = FunctionOfTrig(u, x)
        d = FreeFactors(cos(v), x)
        if And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(cos(v), x), u/sin(v), x)):
            return True
        return False
    pattern4817 = Pattern(Integral(u_, x_), CustomConstraint(With4817))
    def replacement4817(u, x):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(cos(v), x)
        rubi.append(4817)
        return Simp(With(List(Set(d, FreeFactors(cos(v), x))), Dist(-d/Coefficient(v, x, S(1)), Subst(Int(SubstFor(S(1), cos(v)/d, u/sin(v), x), x), x, cos(v)/d), x)), x)
    rule4817 = ReplacementRule(pattern4817, replacement4817)
    pattern4818 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons1667)
    def replacement4818(e, a, p, b, u, x, d, c):
        rubi.append(4818)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)
    rule4818 = ReplacementRule(pattern4818, replacement4818)
    pattern4819 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons1668)
    def replacement4819(e, a, p, b, u, x, d, c):
        rubi.append(4819)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)
    rule4819 = ReplacementRule(pattern4819, replacement4819)
    pattern4820 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons1668)
    def replacement4820(e, a, p, b, u, x, d, c):
        rubi.append(4820)
        return Dist((a + c)**p, Int(ActivateTrig(u), x), x)
    rule4820 = ReplacementRule(pattern4820, replacement4820)
    def With4821(u, y, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern4821 = Pattern(Integral(u_/y_, x_), cons1669, CustomConstraint(With4821))
    def replacement4821(u, y, x):

        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        rubi.append(4821)
        return Simp(q*log(RemoveContent(ActivateTrig(y), x)), x)
    rule4821 = ReplacementRule(pattern4821, replacement4821)
    def With4822(w, u, y, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(ActivateTrig(w*y), ActivateTrig(u), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern4822 = Pattern(Integral(u_/(w_*y_), x_), cons1669, CustomConstraint(With4822))
    def replacement4822(w, u, y, x):

        q = DerivativeDivides(ActivateTrig(w*y), ActivateTrig(u), x)
        rubi.append(4822)
        return Simp(q*log(RemoveContent(ActivateTrig(w*y), x)), x)
    rule4822 = ReplacementRule(pattern4822, replacement4822)
    def With4823(x, u, y, m):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern4823 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons21, cons66, cons1669, CustomConstraint(With4823))
    def replacement4823(x, u, y, m):

        q = DerivativeDivides(ActivateTrig(y), ActivateTrig(u), x)
        rubi.append(4823)
        return Simp(q*ActivateTrig(y**(m + S(1)))/(m + S(1)), x)
    rule4823 = ReplacementRule(pattern4823, replacement4823)
    def With4824(n, y, m, z, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(ActivateTrig(y*z), ActivateTrig(u*z**(-m + n)), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern4824 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons21, cons4, cons66, cons1669, CustomConstraint(With4824))
    def replacement4824(n, y, m, z, u, x):

        q = DerivativeDivides(ActivateTrig(y*z), ActivateTrig(u*z**(-m + n)), x)
        rubi.append(4824)
        return Simp(q*ActivateTrig(y**(m + S(1))*z**(m + S(1)))/(m + S(1)), x)
    rule4824 = ReplacementRule(pattern4824, replacement4824)
    def With4825(a, n, p, F, u, x, d, c):
        v = ActivateTrig(F(c + d*x))
        rubi.append(4825)
        return Dist(a**IntPart(n)*(a*v**p)**FracPart(n)*(v/NonfreeFactors(v, x))**(p*IntPart(n))*NonfreeFactors(v, x)**(-p*FracPart(n)), Int(ActivateTrig(u)*NonfreeFactors(v, x)**(n*p), x), x)
    pattern4825 = Pattern(Integral((F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1)))**n_*WC('u', S(1)), x_), cons2, cons7, cons27, cons4, cons5, cons1656, cons23, cons38)
    rule4825 = ReplacementRule(pattern4825, With4825)
    def With4826(a, n, p, b, F, u, x, d, c):
        v = ActivateTrig(F(c + d*x))
        rubi.append(4826)
        return Dist(a**IntPart(n)*(a*(b*v)**p)**FracPart(n)*(b*v)**(-p*FracPart(n)), Int((b*v)**(n*p)*ActivateTrig(u), x), x)
    pattern4826 = Pattern(Integral(((F_*(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**p_*WC('a', S(1)))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons1656, cons23, cons147)
    rule4826 = ReplacementRule(pattern4826, With4826)
    def With4827(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        if And(Not(FalseQ(v)), FunctionOfQ(NonfreeFactors(tan(v), x), u, x)):
            return True
        return False
    pattern4827 = Pattern(Integral(u_, x_), cons1222, CustomConstraint(With4827))
    def replacement4827(u, x):

        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v), x)
        rubi.append(4827)
        return Dist(d/Coefficient(v, x, 1), Subst(Int(SubstFor(1/(d**2*x**2 + 1), tan(v)/d, u, x), x), x, tan(v)/d), x)
    rule4827 = ReplacementRule(pattern4827, replacement4827)
    pattern4828 = Pattern(Integral((WC('a', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)) + WC('b', S(1))*sec(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons376)
    def replacement4828(n, a, p, b, u, x, d, c):
        rubi.append(4828)
        return Int((a*sin(c + d*x)**n + b)**p*ActivateTrig(u)*sec(c + d*x)**(n*p), x)
    rule4828 = ReplacementRule(pattern4828, replacement4828)
    pattern4829 = Pattern(Integral((WC('a', S(1))*cot(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)) + WC('b', S(1))*csc(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons376)
    def replacement4829(n, a, p, b, u, x, d, c):
        rubi.append(4829)
        return Int((a*cos(c + d*x)**n + b)**p*ActivateTrig(u)*csc(c + d*x)**(n*p), x)
    rule4829 = ReplacementRule(pattern4829, replacement4829)
    pattern4830 = Pattern(Integral(u_*(F_**(x_*WC('d', S(1)) + WC('c', S(0)))*a_ + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons1656, cons85, cons49)
    def replacement4830(n, a, p, d, b, F, u, x, q, c):
        rubi.append(4830)
        return Int(ActivateTrig(u*(a + b*F(c + d*x)**(-p + q))**n*F(c + d*x)**(n*p)), x)
    rule4830 = ReplacementRule(pattern4830, replacement4830)
    pattern4831 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*a_ + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons52, cons1656, cons85, cons49, cons51)
    def replacement4831(e, n, a, p, d, b, r, F, u, x, q, c):
        rubi.append(4831)
        return Int(ActivateTrig(u*(a + b*F(d + e*x)**(-p + q) + c*F(d + e*x)**(-p + r))**n*F(d + e*x)**(n*p)), x)
    rule4831 = ReplacementRule(pattern4831, replacement4831)
    pattern4832 = Pattern(Integral(u_*(F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1)) + F_**(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1)) + a_)**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons5, cons50, cons1656, cons85, cons1670)
    def replacement4832(e, n, a, q, p, b, F, u, x, d, c):
        rubi.append(4832)
        return Int(ActivateTrig(u*(a*F(d + e*x)**(-p) + b + c*F(d + e*x)**(-p + q))**n*F(d + e*x)**(n*p)), x)
    rule4832 = ReplacementRule(pattern4832, replacement4832)
    pattern4833 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons1431)
    def replacement4833(a, n, b, u, x, d, c):
        rubi.append(4833)
        return Int((a*exp(-a*(c + d*x)/b))**n*ActivateTrig(u), x)
    rule4833 = ReplacementRule(pattern4833, replacement4833)
    pattern4834 = Pattern(Integral(u_, x_), cons1671)
    def replacement4834(u, x):
        rubi.append(4834)
        return Int(TrigSimplify(u), x)
    rule4834 = ReplacementRule(pattern4834, replacement4834)
    def With4835(a, v, p, u, x):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        rubi.append(4835)
        return Dist(a**IntPart(p)*vv**(-FracPart(p))*(a*vv)**FracPart(p), Int(uu*vv**p, x), x)
    pattern4835 = Pattern(Integral((a_*v_)**p_*WC('u', S(1)), x_), cons2, cons5, cons147, cons1672)
    rule4835 = ReplacementRule(pattern4835, With4835)
    def With4836(m, v, p, u, x):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        rubi.append(4836)
        return Dist(vv**(-m*FracPart(p))*(vv**m)**FracPart(p), Int(uu*vv**(m*p), x), x)
    pattern4836 = Pattern(Integral((v_**m_)**p_*WC('u', S(1)), x_), cons21, cons5, cons147, cons1672)
    rule4836 = ReplacementRule(pattern4836, With4836)
    def With4837(n, u, m, v, p, w, x):
        uu = ActivateTrig(u)
        vv = ActivateTrig(v)
        ww = ActivateTrig(w)
        rubi.append(4837)
        return Dist(vv**(-m*FracPart(p))*ww**(-n*FracPart(p))*(vv**m*ww**n)**FracPart(p), Int(uu*vv**(m*p)*ww**(n*p), x), x)
    pattern4837 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1)))**p_*WC('u', S(1)), x_), cons21, cons4, cons5, cons147, cons1673)
    rule4837 = ReplacementRule(pattern4837, With4837)
    def With4838(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = ExpandTrig(u, x)
        if SumQ(v):
            return True
        return False
    pattern4838 = Pattern(Integral(u_, x_), cons1669, CustomConstraint(With4838))
    def replacement4838(u, x):

        v = ExpandTrig(u, x)
        rubi.append(4838)
        return Int(v, x)
    rule4838 = ReplacementRule(pattern4838, replacement4838)
    def With4839(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        w = With(List(Set(ShowSteps, False), Set(StepCounter, Null)), Int(SubstFor(S(1)/(x**S(2)*FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x)**S(2) + S(1)), tan(FunctionOfTrig(u, x)/S(2))/FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x), u, x), x))
        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v/S(2)), x)
        if FreeQ(w, Int):
            return True
        return False
    pattern4839 = Pattern(Integral(u_, x_), cons1222, cons1674, CustomConstraint(With4839))
    def replacement4839(u, x):

        w = With(List(Set(ShowSteps, False), Set(StepCounter, Null)), Int(SubstFor(S(1)/(x**S(2)*FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x)**S(2) + S(1)), tan(FunctionOfTrig(u, x)/S(2))/FreeFactors(tan(FunctionOfTrig(u, x)/S(2)), x), u, x), x))
        v = FunctionOfTrig(u, x)
        d = FreeFactors(tan(v/S(2)), x)
        rubi.append(4839)
        return Simp(CompoundExpression(Set(d, FreeFactors(tan(v/2), x)), Dist(2*d/Coefficient(v, x, 1), Subst(Int(SubstFor(1/(d**2*x**2 + 1), tan(v/2)/d, u, x), x), x, tan(v/2)/d), x)), x)
    rule4839 = ReplacementRule(pattern4839, replacement4839)
    def With4840(u, x):
        v = ActivateTrig(u)
        rubi.append(4840)
        return Int(v, x)
    pattern4840 = Pattern(Integral(u_, x_), cons1669)
    rule4840 = ReplacementRule(pattern4840, With4840)
    pattern4841 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons4, cons62, cons584)
    def replacement4841(a, n, m, b, x, d, c):
        rubi.append(4841)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*sin(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*sin(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
    rule4841 = ReplacementRule(pattern4841, replacement4841)
    pattern4842 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons62, cons584)
    def replacement4842(a, n, m, b, x, d, c):
        rubi.append(4842)
        return Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*cos(a + b*x)**(n + S(1)), x), x) - Simp((c + d*x)**m*cos(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
    rule4842 = ReplacementRule(pattern4842, replacement4842)
    pattern4843 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons464)
    def replacement4843(a, n, m, p, b, x, d, c):
        rubi.append(4843)
        return Int(ExpandTrigReduce((c + d*x)**m, sin(a + b*x)**n*cos(a + b*x)**p, x), x)
    rule4843 = ReplacementRule(pattern4843, replacement4843)
    pattern4844 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons464)
    def replacement4844(a, n, m, p, b, x, d, c):
        rubi.append(4844)
        return -Int((c + d*x)**m*sin(a + b*x)**n*tan(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*sin(a + b*x)**(n + S(-2))*tan(a + b*x)**p, x)
    rule4844 = ReplacementRule(pattern4844, replacement4844)
    pattern4845 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons464)
    def replacement4845(a, n, m, p, b, x, d, c):
        rubi.append(4845)
        return -Int((c + d*x)**m*cos(a + b*x)**n*cot(a + b*x)**(p + S(-2)), x) + Int((c + d*x)**m*cos(a + b*x)**(n + S(-2))*cot(a + b*x)**p, x)
    rule4845 = ReplacementRule(pattern4845, replacement4845)
    pattern4846 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons1675, cons31, cons168)
    def replacement4846(a, n, m, p, b, x, d, c):
        rubi.append(4846)
        return -Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*sec(a + b*x)**n, x), x) + Simp((c + d*x)**m*sec(a + b*x)**n/(b*n), x)
    rule4846 = ReplacementRule(pattern4846, replacement4846)
    pattern4847 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons1675, cons31, cons168)
    def replacement4847(a, n, m, p, b, x, d, c):
        rubi.append(4847)
        return Dist(d*m/(b*n), Int((c + d*x)**(m + S(-1))*csc(a + b*x)**n, x), x) - Simp((c + d*x)**m*csc(a + b*x)**n/(b*n), x)
    rule4847 = ReplacementRule(pattern4847, replacement4847)
    pattern4848 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons4, cons62, cons584)
    def replacement4848(a, n, m, b, x, d, c):
        rubi.append(4848)
        return -Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*tan(a + b*x)**(n + S(1)), x), x) + Simp((c + d*x)**m*tan(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
    rule4848 = ReplacementRule(pattern4848, replacement4848)
    pattern4849 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons4, cons62, cons584)
    def replacement4849(a, n, m, b, x, d, c):
        rubi.append(4849)
        return Dist(d*m/(b*(n + S(1))), Int((c + d*x)**(m + S(-1))*cot(a + b*x)**(n + S(1)), x), x) - Simp((c + d*x)**m*cot(a + b*x)**(n + S(1))/(b*(n + S(1))), x)
    rule4849 = ReplacementRule(pattern4849, replacement4849)
    pattern4850 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_*sec(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons21, cons1400)
    def replacement4850(a, m, p, b, x, d, c):
        rubi.append(4850)
        return -Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))*sec(a + b*x), x) + Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))*sec(a + b*x)**S(3), x)
    rule4850 = ReplacementRule(pattern4850, replacement4850)
    pattern4851 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**p_*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1400)
    def replacement4851(a, n, m, p, b, x, d, c):
        rubi.append(4851)
        return -Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))*sec(a + b*x)**n, x) + Int((c + d*x)**m*tan(a + b*x)**(p + S(-2))*sec(a + b*x)**(n + S(2)), x)
    rule4851 = ReplacementRule(pattern4851, replacement4851)
    pattern4852 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**p_*csc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons21, cons1400)
    def replacement4852(a, m, p, b, x, d, c):
        rubi.append(4852)
        return -Int((c + d*x)**m*cot(a + b*x)**(p + S(-2))*csc(a + b*x), x) + Int((c + d*x)**m*cot(a + b*x)**(p + S(-2))*csc(a + b*x)**S(3), x)
    rule4852 = ReplacementRule(pattern4852, replacement4852)
    pattern4853 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**p_*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1400)
    def replacement4853(a, n, m, p, b, x, d, c):
        rubi.append(4853)
        return -Int((c + d*x)**m*cot(a + b*x)**(p + S(-2))*csc(a + b*x)**n, x) + Int((c + d*x)**m*cot(a + b*x)**(p + S(-2))*csc(a + b*x)**(n + S(2)), x)
    rule4853 = ReplacementRule(pattern4853, replacement4853)
    def With4854(a, n, m, p, b, x, d, c):
        u = IntHide(tan(a + b*x)**p*sec(a + b*x)**n, x)
        rubi.append(4854)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
    pattern4854 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons62, cons1676)
    rule4854 = ReplacementRule(pattern4854, With4854)
    def With4855(a, n, m, p, b, x, d, c):
        u = IntHide(cot(a + b*x)**p*csc(a + b*x)**n, x)
        rubi.append(4855)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
    pattern4855 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cot(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons62, cons1676)
    rule4855 = ReplacementRule(pattern4855, With4855)
    pattern4856 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons31, cons85)
    def replacement4856(a, n, m, b, x, d, c):
        rubi.append(4856)
        return Dist(S(2)**n, Int((c + d*x)**m*csc(S(2)*a + S(2)*b*x)**n, x), x)
    rule4856 = ReplacementRule(pattern4856, replacement4856)
    def With4857(a, n, m, p, b, x, d, c):
        u = IntHide(csc(a + b*x)**n*sec(a + b*x)**p, x)
        rubi.append(4857)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
    pattern4857 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*csc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*sec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons376, cons31, cons168, cons1677)
    rule4857 = ReplacementRule(pattern4857, With4857)
    pattern4858 = Pattern(Integral(F_**v_*G_**w_*u_**WC('m', S(1)), x_), cons21, cons4, cons5, cons1678, cons1679, cons1680, cons812, cons813)
    def replacement4858(w, G, n, m, v, p, F, u, x):
        rubi.append(4858)
        return Int(ExpandToSum(u, x)**m*F(ExpandToSum(v, x))**n*G(ExpandToSum(v, x))**p, x)
    rule4858 = ReplacementRule(pattern4858, replacement4858)
    pattern4859 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4859(e, n, a, m, b, f, x, d, c):
        rubi.append(4859)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*sin(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*sin(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4859 = ReplacementRule(pattern4859, replacement4859)
    pattern4860 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4860(e, n, a, m, b, f, x, d, c):
        rubi.append(4860)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b*cos(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b*cos(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4860 = ReplacementRule(pattern4860, replacement4860)
    pattern4861 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sec(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4861(e, n, a, m, b, f, x, d, c):
        rubi.append(4861)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*tan(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*tan(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4861 = ReplacementRule(pattern4861, replacement4861)
    pattern4862 = Pattern(Integral((a_ + WC('b', S(1))*cot(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*csc(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4862(e, n, a, m, b, f, x, d, c):
        rubi.append(4862)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b*cot(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b*cot(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4862 = ReplacementRule(pattern4862, replacement4862)
    pattern4863 = Pattern(Integral((a_ + WC('b', S(1))*sec(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))*sec(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4863(e, n, a, m, b, f, x, d, c):
        rubi.append(4863)
        return -Dist(f*m/(b*d*(n + S(1))), Int((a + b*sec(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) + Simp((a + b*sec(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4863 = ReplacementRule(pattern4863, replacement4863)
    pattern4864 = Pattern(Integral((a_ + WC('b', S(1))*csc(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cot(x_*WC('d', S(1)) + WC('c', S(0)))*csc(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons584)
    def replacement4864(e, n, a, m, b, f, x, d, c):
        rubi.append(4864)
        return Dist(f*m/(b*d*(n + S(1))), Int((a + b*csc(c + d*x))**(n + S(1))*(e + f*x)**(m + S(-1)), x), x) - Simp((a + b*csc(c + d*x))**(n + S(1))*(e + f*x)**m/(b*d*(n + S(1))), x)
    rule4864 = ReplacementRule(pattern4864, replacement4864)
    pattern4865 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons555, cons17)
    def replacement4865(e, a, q, m, p, b, f, x, d, c):
        rubi.append(4865)
        return Int(ExpandTrigReduce((e + f*x)**m, sin(a + b*x)**p*sin(c + d*x)**q, x), x)
    rule4865 = ReplacementRule(pattern4865, replacement4865)
    pattern4866 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons555, cons17)
    def replacement4866(e, a, q, m, p, b, f, x, d, c):
        rubi.append(4866)
        return Int(ExpandTrigReduce((e + f*x)**m, cos(a + b*x)**p*cos(c + d*x)**q, x), x)
    rule4866 = ReplacementRule(pattern4866, replacement4866)
    pattern4867 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons555)
    def replacement4867(e, a, q, m, p, b, f, x, d, c):
        rubi.append(4867)
        return Int(ExpandTrigReduce((e + f*x)**m, sin(a + b*x)**p*cos(c + d*x)**q, x), x)
    rule4867 = ReplacementRule(pattern4867, replacement4867)
    pattern4868 = Pattern(Integral(F_**(x_*WC('b', S(1)) + WC('a', S(0)))*G_**(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons1681, cons1682, cons555, cons25, cons1683)
    def replacement4868(e, G, a, q, m, p, b, F, f, x, d, c):
        rubi.append(4868)
        return Int(ExpandTrigExpand((e + f*x)**m*G(c + d*x)**q, F, c + d*x, p, b/d, x), x)
    rule4868 = ReplacementRule(pattern4868, replacement4868)
    pattern4869 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1684)
    def replacement4869(e, a, b, F, x, d, c):
        rubi.append(4869)
        return -Simp(F**(c*(a + b*x))*e*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
    rule4869 = ReplacementRule(pattern4869, replacement4869)
    pattern4870 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1684)
    def replacement4870(e, a, b, F, x, d, c):
        rubi.append(4870)
        return Simp(F**(c*(a + b*x))*e*sin(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)), x)
    rule4870 = ReplacementRule(pattern4870, replacement4870)
    pattern4871 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1685, cons87, cons165)
    def replacement4871(e, a, n, b, F, x, d, c):
        rubi.append(4871)
        return Dist(e**S(2)*n*(n + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*sin(d + e*x)**(n + S(-2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**n/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*sin(d + e*x)**(n + S(-1))*cos(d + e*x)/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
    rule4871 = ReplacementRule(pattern4871, replacement4871)
    pattern4872 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1686, cons31, cons166)
    def replacement4872(e, a, m, b, F, x, d, c):
        rubi.append(4872)
        return Dist(e**S(2)*m*(m + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), Int(F**(c*(a + b*x))*cos(d + e*x)**(m + S(-2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**m/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), x) + Simp(F**(c*(a + b*x))*e*m*sin(d + e*x)*cos(d + e*x)**(m + S(-1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*m**S(2)), x)
    rule4872 = ReplacementRule(pattern4872, replacement4872)
    pattern4873 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons1687, cons584, cons1387)
    def replacement4873(e, a, n, b, F, x, d, c):
        rubi.append(4873)
        return Simp(F**(c*(a + b*x))*sin(d + e*x)**(n + S(1))*cos(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
    rule4873 = ReplacementRule(pattern4873, replacement4873)
    pattern4874 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons1687, cons584, cons1387)
    def replacement4874(e, a, n, b, F, x, d, c):
        rubi.append(4874)
        return -Simp(F**(c*(a + b*x))*sin(d + e*x)*cos(d + e*x)**(n + S(1))/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
    rule4874 = ReplacementRule(pattern4874, replacement4874)
    pattern4875 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1688, cons87, cons89, cons1434)
    def replacement4875(e, a, n, b, F, x, d, c):
        rubi.append(4875)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*sin(d + e*x)**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*sin(d + e*x)**(n + S(1))*cos(d + e*x)/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sin(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
    rule4875 = ReplacementRule(pattern4875, replacement4875)
    pattern4876 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1688, cons87, cons89, cons1434)
    def replacement4876(e, a, n, b, F, x, d, c):
        rubi.append(4876)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(2))**S(2))/(e**S(2)*(n + S(1))*(n + S(2))), Int(F**(c*(a + b*x))*cos(d + e*x)**(n + S(2)), x), x) - Simp(F**(c*(a + b*x))*sin(d + e*x)*cos(d + e*x)**(n + S(1))/(e*(n + S(1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*cos(d + e*x)**(n + S(2))/(e**S(2)*(n + S(1))*(n + S(2))), x)
    rule4876 = ReplacementRule(pattern4876, replacement4876)
    pattern4877 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons23)
    def replacement4877(e, a, n, b, F, x, d, c):
        rubi.append(4877)
        return Dist((exp(S(2)*I*(d + e*x)) + S(-1))**(-n)*exp(I*n*(d + e*x))*sin(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(-1))**n*exp(-I*n*(d + e*x)), x), x)
    rule4877 = ReplacementRule(pattern4877, replacement4877)
    pattern4878 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons23)
    def replacement4878(e, a, n, b, F, x, d, c):
        rubi.append(4878)
        return Dist((exp(S(2)*I*(d + e*x)) + S(1))**(-n)*exp(I*n*(d + e*x))*cos(d + e*x)**n, Int(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(1))**n*exp(-I*n*(d + e*x)), x), x)
    rule4878 = ReplacementRule(pattern4878, replacement4878)
    pattern4879 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons85)
    def replacement4879(e, a, n, b, F, x, d, c):
        rubi.append(4879)
        return Dist(I**n, Int(ExpandIntegrand(F**(c*(a + b*x))*(-exp(S(2)*I*(d + e*x)) + S(1))**n*(exp(S(2)*I*(d + e*x)) + S(1))**(-n), x), x), x)
    rule4879 = ReplacementRule(pattern4879, replacement4879)
    pattern4880 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*cot(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons85)
    def replacement4880(e, n, a, b, F, x, d, c):
        rubi.append(4880)
        return Dist((-I)**n, Int(ExpandIntegrand(F**(c*(a + b*x))*(-exp(S(2)*I*(d + e*x)) + S(1))**(-n)*(exp(S(2)*I*(d + e*x)) + S(1))**n, x), x), x)
    rule4880 = ReplacementRule(pattern4880, replacement4880)
    pattern4881 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1685, cons87, cons89)
    def replacement4881(e, a, n, b, F, x, d, c):
        rubi.append(4881)
        return Dist(e**S(2)*n*(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*sec(d + e*x)**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*sec(d + e*x)**n/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) - Simp(F**(c*(a + b*x))*e*n*sin(d + e*x)*sec(d + e*x)**(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
    rule4881 = ReplacementRule(pattern4881, replacement4881)
    pattern4882 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1685, cons87, cons89)
    def replacement4882(e, a, n, b, F, x, d, c):
        rubi.append(4882)
        return Dist(e**S(2)*n*(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), Int(F**(c*(a + b*x))*csc(d + e*x)**(n + S(2)), x), x) + Simp(F**(c*(a + b*x))*b*c*log(F)*csc(d + e*x)**n/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x) + Simp(F**(c*(a + b*x))*e*n*cos(d + e*x)*csc(d + e*x)**(n + S(1))/(b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*n**S(2)), x)
    rule4882 = ReplacementRule(pattern4882, replacement4882)
    pattern4883 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons1689, cons1494, cons963)
    def replacement4883(e, a, n, b, F, x, d, c):
        rubi.append(4883)
        return Simp(F**(c*(a + b*x))*sin(d + e*x)*sec(d + e*x)**(n + S(-1))/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sec(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
    rule4883 = ReplacementRule(pattern4883, replacement4883)
    pattern4884 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons1689, cons1494, cons963)
    def replacement4884(e, a, n, b, F, x, d, c):
        rubi.append(4884)
        return Simp(F**(c*(a + b*x))*cos(d + e*x)*csc(d + e*x)**(n + S(-1))/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*csc(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
    rule4884 = ReplacementRule(pattern4884, replacement4884)
    pattern4885 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1690, cons87, cons165, cons1636)
    def replacement4885(e, a, n, b, F, x, d, c):
        rubi.append(4885)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*sec(d + e*x)**(n + S(-2)), x), x) + Simp(F**(c*(a + b*x))*sin(d + e*x)*sec(d + e*x)**(n + S(-1))/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*sec(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
    rule4885 = ReplacementRule(pattern4885, replacement4885)
    pattern4886 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1690, cons87, cons165, cons1636)
    def replacement4886(e, a, n, b, F, x, d, c):
        rubi.append(4886)
        return Dist((b**S(2)*c**S(2)*log(F)**S(2) + e**S(2)*(n + S(-2))**S(2))/(e**S(2)*(n + S(-2))*(n + S(-1))), Int(F**(c*(a + b*x))*csc(d + e*x)**(n + S(-2)), x), x) - Simp(F**(c*(a + b*x))*cos(d + e*x)*csc(d + e*x)**(n + S(-1))/(e*(n + S(-1))), x) - Simp(F**(c*(a + b*x))*b*c*log(F)*csc(d + e*x)**(n + S(-2))/(e**S(2)*(n + S(-2))*(n + S(-1))), x)
    rule4886 = ReplacementRule(pattern4886, replacement4886)
    pattern4887 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons85)
    def replacement4887(e, a, n, b, F, x, d, c):
        rubi.append(4887)
        return Simp(S(2)**n*F**(c*(a + b*x))*Hypergeometric2F1(n, -I*b*c*log(F)/(S(2)*e) + n/S(2), -I*b*c*log(F)/(S(2)*e) + n/S(2) + S(1), -exp(S(2)*I*(d + e*x)))*exp(I*n*(d + e*x))/(b*c*log(F) + I*e*n), x)
    rule4887 = ReplacementRule(pattern4887, replacement4887)
    pattern4888 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons85)
    def replacement4888(e, n, a, b, F, x, d, c):
        rubi.append(4888)
        return Simp(F**(c*(a + b*x))*(-S(2)*I)**n*Hypergeometric2F1(n, -I*b*c*log(F)/(S(2)*e) + n/S(2), -I*b*c*log(F)/(S(2)*e) + n/S(2) + S(1), exp(S(2)*I*(d + e*x)))*exp(I*n*(d + e*x))/(b*c*log(F) + I*e*n), x)
    rule4888 = ReplacementRule(pattern4888, replacement4888)
    pattern4889 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons23)
    def replacement4889(e, a, n, b, F, x, d, c):
        rubi.append(4889)
        return Dist((exp(S(2)*I*(d + e*x)) + S(1))**n*exp(-I*n*(d + e*x))*sec(d + e*x)**n, Int(SimplifyIntegrand(F**(c*(a + b*x))*(exp(S(2)*I*(d + e*x)) + S(1))**(-n)*exp(I*n*(d + e*x)), x), x), x)
    rule4889 = ReplacementRule(pattern4889, replacement4889)
    pattern4890 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons23)
    def replacement4890(e, n, a, b, F, x, d, c):
        rubi.append(4890)
        return Dist((S(1) - exp(-S(2)*I*(d + e*x)))**n*exp(I*n*(d + e*x))*csc(d + e*x)**n, Int(SimplifyIntegrand(F**(c*(a + b*x))*(S(1) - exp(-S(2)*I*(d + e*x)))**(-n)*exp(-I*n*(d + e*x)), x), x), x)
    rule4890 = ReplacementRule(pattern4890, replacement4890)
    pattern4891 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1691, cons196)
    def replacement4891(e, a, n, g, b, F, f, x, d, c):
        rubi.append(4891)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*cos(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2))**(S(2)*n), x), x)
    rule4891 = ReplacementRule(pattern4891, replacement4891)
    pattern4892 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1692, cons196)
    def replacement4892(e, n, a, g, b, F, f, x, d, c):
        rubi.append(4892)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*cos(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
    rule4892 = ReplacementRule(pattern4892, replacement4892)
    pattern4893 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1011, cons196)
    def replacement4893(e, n, a, g, b, F, f, x, d, c):
        rubi.append(4893)
        return Dist(S(2)**n*f**n, Int(F**(c*(a + b*x))*sin(d/S(2) + e*x/S(2))**(S(2)*n), x), x)
    rule4893 = ReplacementRule(pattern4893, replacement4893)
    pattern4894 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1691, cons150, cons1543)
    def replacement4894(e, a, n, m, g, b, F, f, x, d, c):
        rubi.append(4894)
        return Dist(g**n, Int(F**(c*(a + b*x))*(-tan(-Pi*f/(S(4)*g) + d/S(2) + e*x/S(2)))**m, x), x)
    rule4894 = ReplacementRule(pattern4894, replacement4894)
    pattern4895 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1692, cons150, cons1543)
    def replacement4895(e, n, a, m, g, b, F, f, x, d, c):
        rubi.append(4895)
        return Dist(f**n, Int(F**(c*(a + b*x))*tan(d/S(2) + e*x/S(2))**m, x), x)
    rule4895 = ReplacementRule(pattern4895, replacement4895)
    pattern4896 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1011, cons150, cons1543)
    def replacement4896(e, n, a, m, g, b, F, f, x, d, c):
        rubi.append(4896)
        return Dist(f**n, Int(F**(c*(a + b*x))*cot(d/S(2) + e*x/S(2))**m, x), x)
    rule4896 = ReplacementRule(pattern4896, replacement4896)
    pattern4897 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons1691, cons1693, cons1694)
    def replacement4897(e, a, h, i, g, b, F, f, x, d, c):
        rubi.append(4897)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*cos(d + e*x)/(f + g*sin(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*cos(d + e*x))/(f + g*sin(d + e*x)), x)
    rule4897 = ReplacementRule(pattern4897, replacement4897)
    pattern4898 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(h_ + WC('i', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(f_ + WC('g', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons1691, cons1693, cons1695)
    def replacement4898(e, a, h, i, g, b, F, f, x, d, c):
        rubi.append(4898)
        return Dist(S(2)*i, Int(F**(c*(a + b*x))*sin(d + e*x)/(f + g*cos(d + e*x)), x), x) + Int(F**(c*(a + b*x))*(h - i*sin(d + e*x))/(f + g*cos(d + e*x)), x)
    rule4898 = ReplacementRule(pattern4898, replacement4898)
    pattern4899 = Pattern(Integral(F_**(u_*WC('c', S(1)))*G_**v_, x_), cons1091, cons7, cons4, cons1679, cons810, cons811)
    def replacement4899(G, n, v, F, u, x, c):
        rubi.append(4899)
        return Int(F**(c*ExpandToSum(u, x))*G(ExpandToSum(v, x))**n, x)
    rule4899 = ReplacementRule(pattern4899, replacement4899)
    def With4900(e, a, n, m, b, F, x, d, c):
        u = IntHide(F**(c*(a + b*x))*sin(d + e*x)**n, x)
        rubi.append(4900)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)
    pattern4900 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons31, cons168, cons148)
    rule4900 = ReplacementRule(pattern4900, With4900)
    def With4901(e, n, a, m, b, F, x, d, c):
        u = IntHide(F**(c*(a + b*x))*cos(d + e*x)**n, x)
        rubi.append(4901)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)
    pattern4901 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons31, cons168, cons148)
    rule4901 = ReplacementRule(pattern4901, With4901)
    pattern4902 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons528)
    def replacement4902(e, n, a, m, g, b, F, f, x, d, c):
        rubi.append(4902)
        return Int(ExpandTrigReduce(F**(c*(a + b*x)), sin(d + e*x)**m*cos(f + g*x)**n, x), x)
    rule4902 = ReplacementRule(pattern4902, replacement4902)
    pattern4903 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons1696)
    def replacement4903(e, n, a, m, g, p, b, F, f, x, d, c):
        rubi.append(4903)
        return Int(ExpandTrigReduce(F**(c*(a + b*x))*x**p, sin(d + e*x)**m*cos(f + g*x)**n, x), x)
    rule4903 = ReplacementRule(pattern4903, replacement4903)
    pattern4904 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*G_**(x_*WC('e', S(1)) + WC('d', S(0)))*H_**(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons528, cons1679, cons1697)
    def replacement4904(H, e, G, a, n, m, b, F, x, d, c):
        rubi.append(4904)
        return Int(ExpandTrigToExp(F**(c*(a + b*x)), G(d + e*x)**m*H(d + e*x)**n, x), x)
    rule4904 = ReplacementRule(pattern4904, replacement4904)
    pattern4905 = Pattern(Integral(F_**u_*sin(v_)**WC('n', S(1)), x_), cons1091, cons1698, cons1699, cons148)
    def replacement4905(n, v, F, u, x):
        rubi.append(4905)
        return Int(ExpandTrigToExp(F**u, sin(v)**n, x), x)
    rule4905 = ReplacementRule(pattern4905, replacement4905)
    pattern4906 = Pattern(Integral(F_**u_*cos(v_)**WC('n', S(1)), x_), cons1091, cons1698, cons1699, cons148)
    def replacement4906(n, v, F, u, x):
        rubi.append(4906)
        return Int(ExpandTrigToExp(F**u, cos(v)**n, x), x)
    rule4906 = ReplacementRule(pattern4906, replacement4906)
    pattern4907 = Pattern(Integral(F_**u_*sin(v_)**WC('m', S(1))*cos(v_)**WC('n', S(1)), x_), cons1091, cons1698, cons1699, cons528)
    def replacement4907(n, m, v, F, u, x):
        rubi.append(4907)
        return Int(ExpandTrigToExp(F**u, sin(v)**m*cos(v)**n, x), x)
    rule4907 = ReplacementRule(pattern4907, replacement4907)
    pattern4908 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1700, cons54)
    def replacement4908(a, n, p, b, x, c):
        rubi.append(4908)
        return Simp(x*(p + S(2))*sin(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) + Simp(x*sin(a + b*log(c*x**n))**(p + S(2))*cot(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4908 = ReplacementRule(pattern4908, replacement4908)
    pattern4909 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1700, cons54)
    def replacement4909(a, n, p, b, x, c):
        rubi.append(4909)
        return Simp(x*(p + S(2))*cos(a + b*log(c*x**n))**(p + S(2))/(p + S(1)), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4909 = ReplacementRule(pattern4909, replacement4909)
    pattern4910 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons128, cons1701)
    def replacement4910(a, n, p, b, x, c):
        rubi.append(4910)
        return Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/(S(2)*b*n*p) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/(S(2)*b*n*p))**p, x), x)
    rule4910 = ReplacementRule(pattern4910, replacement4910)
    pattern4911 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons128, cons1701)
    def replacement4911(a, n, p, b, x, c):
        rubi.append(4911)
        return Int(ExpandIntegrand((-(c*x**n)**(S(1)/(n*p))*exp(-a*b*n*p)/S(2) + (c*x**n)**(-S(1)/(n*p))*exp(a*b*n*p)/S(2))**p, x), x)
    rule4911 = ReplacementRule(pattern4911, replacement4911)
    pattern4912 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons4, cons1702)
    def replacement4912(a, n, b, x, c):
        rubi.append(4912)
        return Simp(x*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x) - Simp(b*n*x*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x)
    rule4912 = ReplacementRule(pattern4912, replacement4912)
    pattern4913 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons4, cons1702)
    def replacement4913(a, n, b, x, c):
        rubi.append(4913)
        return Simp(x*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x) + Simp(b*n*x*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + S(1)), x)
    rule4913 = ReplacementRule(pattern4913, replacement4913)
    pattern4914 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons146, cons1703)
    def replacement4914(a, n, p, b, x, c):
        rubi.append(4914)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(sin(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x*sin(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) - Simp(b*n*p*x*sin(a + b*log(c*x**n))**(p + S(-1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
    rule4914 = ReplacementRule(pattern4914, replacement4914)
    pattern4915 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons146, cons1703)
    def replacement4915(a, n, p, b, x, c):
        rubi.append(4915)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(cos(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x*cos(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) + Simp(b*n*p*x*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
    rule4915 = ReplacementRule(pattern4915, replacement4915)
    pattern4916 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons137, cons1497, cons1704)
    def replacement4916(a, n, p, b, x, c):
        rubi.append(4916)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(sin(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x*sin(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) + Simp(x*sin(a + b*log(c*x**n))**(p + S(2))*cot(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4916 = ReplacementRule(pattern4916, replacement4916)
    pattern4917 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons137, cons1497, cons1704)
    def replacement4917(a, n, p, b, x, c):
        rubi.append(4917)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(cos(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x) - Simp(x*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4917 = ReplacementRule(pattern4917, replacement4917)
    pattern4918 = Pattern(Integral(sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1703)
    def replacement4918(a, n, p, b, x, c):
        rubi.append(4918)
        return Simp(x*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*(-I*(c*x**n)**(I*b)*exp(I*a) + I*(c*x**n)**(-I*b)*exp(-I*a))**p*Hypergeometric2F1(-p, -I*(-I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + S(1)), x)
    rule4918 = ReplacementRule(pattern4918, replacement4918)
    pattern4919 = Pattern(Integral(cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1703)
    def replacement4919(a, n, p, b, x, c):
        rubi.append(4919)
        return Simp(x*((c*x**n)**(I*b)*exp(I*a) + (c*x**n)**(-I*b)*exp(-I*a))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*Hypergeometric2F1(-p, -I*(-I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + S(1)), x)
    rule4919 = ReplacementRule(pattern4919, replacement4919)
    pattern4920 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1705, cons54, cons66)
    def replacement4920(n, a, m, p, b, x, c):
        rubi.append(4920)
        return Simp(x**(m + S(1))*(p + S(2))*sin(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) + Simp(x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))*cot(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4920 = ReplacementRule(pattern4920, replacement4920)
    pattern4921 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1705, cons54, cons66)
    def replacement4921(a, n, m, p, b, x, c):
        rubi.append(4921)
        return Simp(x**(m + S(1))*(p + S(2))*cos(a + b*log(c*x**n))**(p + S(2))/((m + S(1))*(p + S(1))), x) - Simp(x**(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x)
    rule4921 = ReplacementRule(pattern4921, replacement4921)
    pattern4922 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons128, cons1706)
    def replacement4922(n, a, m, p, b, x, c):
        rubi.append(4922)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*(-(c*x**n)**((m + S(1))/(n*p))*(m + S(1))*exp(-a*b*n*p/(m + S(1)))/(b*n*p) + (c*x**n)**(-(m + S(1))/(n*p))*(m + S(1))*exp(a*b*n*p/(m + S(1)))/(b*n*p))**p, x), x), x)
    rule4922 = ReplacementRule(pattern4922, replacement4922)
    pattern4923 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons128, cons1706)
    def replacement4923(a, n, m, p, b, x, c):
        rubi.append(4923)
        return Dist(S(2)**(-p), Int(ExpandIntegrand(x**m*(-(c*x**n)**((m + S(1))/(n*p))*exp(-a*b*n*p/(m + S(1))) + (c*x**n)**(-(m + S(1))/(n*p))*exp(a*b*n*p/(m + S(1))))**p, x), x), x)
    rule4923 = ReplacementRule(pattern4923, replacement4923)
    pattern4924 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons21, cons4, cons1707)
    def replacement4924(n, a, m, b, x, c):
        rubi.append(4924)
        return Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*x**(m + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x)
    rule4924 = ReplacementRule(pattern4924, replacement4924)
    pattern4925 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons21, cons4, cons1707)
    def replacement4925(a, n, m, b, x, c):
        rubi.append(4925)
        return Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*x**(m + S(1))*sin(a + b*log(c*x**n))/(b**S(2)*n**S(2) + (m + S(1))**S(2)), x)
    rule4925 = ReplacementRule(pattern4925, replacement4925)
    pattern4926 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons146, cons1708)
    def replacement4926(n, a, m, p, b, x, c):
        rubi.append(4926)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*sin(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(-1))*cos(a + b*log(c*x**n))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
    rule4926 = ReplacementRule(pattern4926, replacement4926)
    pattern4927 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons146, cons1708)
    def replacement4927(a, n, m, p, b, x, c):
        rubi.append(4927)
        return Dist(b**S(2)*n**S(2)*p*(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*cos(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*sin(a + b*log(c*x**n))*cos(a + b*log(c*x**n))**(p + S(-1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
    rule4927 = ReplacementRule(pattern4927, replacement4927)
    pattern4928 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons137, cons1497, cons1709)
    def replacement4928(n, a, m, p, b, x, c):
        rubi.append(4928)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*sin(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x**(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))*cot(a + b*log(c*x**n))/(b*n*(p + S(1))), x) - Simp(x**(m + S(1))*(m + S(1))*sin(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
    rule4928 = ReplacementRule(pattern4928, replacement4928)
    pattern4929 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons137, cons1497, cons1709)
    def replacement4929(a, n, m, p, b, x, c):
        rubi.append(4929)
        return Dist((b**S(2)*n**S(2)*(p + S(2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**m*cos(a + b*log(c*x**n))**(p + S(2)), x), x) - Simp(x**(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))*tan(a + b*log(c*x**n))/(b*n*(p + S(1))), x) - Simp(x**(m + S(1))*(m + S(1))*cos(a + b*log(c*x**n))**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)
    rule4929 = ReplacementRule(pattern4929, replacement4929)
    pattern4930 = Pattern(Integral(x_**WC('m', S(1))*sin(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1708)
    def replacement4930(n, a, m, p, b, x, c):
        rubi.append(4930)
        return Simp(x**(m + S(1))*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*(-I*(c*x**n)**(I*b)*exp(I*a) + I*(c*x**n)**(-I*b)*exp(-I*a))**p*Hypergeometric2F1(-p, -I*(-I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + m + S(1)), x)
    rule4930 = ReplacementRule(pattern4930, replacement4930)
    pattern4931 = Pattern(Integral(x_**WC('m', S(1))*cos(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1708)
    def replacement4931(a, n, m, p, b, x, c):
        rubi.append(4931)
        return Simp(x**(m + S(1))*((c*x**n)**(I*b)*exp(I*a) + (c*x**n)**(-I*b)*exp(-I*a))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**(-p)*Hypergeometric2F1(-p, -I*(-I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(-I*b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(-I*b*n*p + m + S(1)), x)
    rule4931 = ReplacementRule(pattern4931, replacement4931)
    pattern4932 = Pattern(Integral(sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons4, cons1710)
    def replacement4932(a, n, b, x, c):
        rubi.append(4932)
        return Dist(S(2)*exp(a*b*n), Int((c*x**n)**(S(1)/n)/((c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x), x)
    rule4932 = ReplacementRule(pattern4932, replacement4932)
    pattern4933 = Pattern(Integral(csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons4, cons1710)
    def replacement4933(a, n, b, x, c):
        rubi.append(4933)
        return Dist(S(2)*b*n*exp(a*b*n), Int((c*x**n)**(S(1)/n)/(-(c*x**n)**(S(2)/n) + exp(S(2)*a*b*n)), x), x)
    rule4933 = ReplacementRule(pattern4933, replacement4933)
    pattern4934 = Pattern(Integral(sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1711, cons1637)
    def replacement4934(a, n, p, b, x, c):
        rubi.append(4934)
        return Simp(x*(p + S(-2))*sec(a + b*log(c*x**n))**(p + S(-2))/(p + S(-1)), x) + Simp(x*tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4934 = ReplacementRule(pattern4934, replacement4934)
    pattern4935 = Pattern(Integral(csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons1711, cons1637)
    def replacement4935(a, n, p, b, x, c):
        rubi.append(4935)
        return Simp(x*(p + S(-2))*csc(a + b*log(c*x**n))**(p + S(-2))/(p + S(-1)), x) - Simp(x*cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4935 = ReplacementRule(pattern4935, replacement4935)
    pattern4936 = Pattern(Integral(sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons146, cons1712, cons1713)
    def replacement4936(a, n, p, b, x, c):
        rubi.append(4936)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(sec(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x*sec(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) + Simp(x*tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4936 = ReplacementRule(pattern4936, replacement4936)
    pattern4937 = Pattern(Integral(csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons146, cons1712, cons1713)
    def replacement4937(a, n, p, b, x, c):
        rubi.append(4937)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + S(1))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(csc(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x*csc(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x) - Simp(x*cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4937 = ReplacementRule(pattern4937, replacement4937)
    pattern4938 = Pattern(Integral(sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons137, cons1703)
    def replacement4938(a, n, p, b, x, c):
        rubi.append(4938)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(sec(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x*sec(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) - Simp(b*n*p*x*sin(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
    rule4938 = ReplacementRule(pattern4938, replacement4938)
    pattern4939 = Pattern(Integral(csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons4, cons13, cons137, cons1703)
    def replacement4939(a, n, p, b, x, c):
        rubi.append(4939)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), Int(csc(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x*csc(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + S(1)), x) + Simp(b*n*p*x*cos(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + S(1)), x)
    rule4939 = ReplacementRule(pattern4939, replacement4939)
    pattern4940 = Pattern(Integral(sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons1703)
    def replacement4940(a, n, p, b, x, c):
        rubi.append(4940)
        return Simp(x*((c*x**n)**(I*b)*exp(I*a)/((c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + S(1)), x)
    rule4940 = ReplacementRule(pattern4940, replacement4940)
    pattern4941 = Pattern(Integral(csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons1703)
    def replacement4941(a, n, p, b, x, c):
        rubi.append(4941)
        return Simp(x*(-I*(c*x**n)**(I*b)*exp(I*a)/(-(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + S(1)), x)
    rule4941 = ReplacementRule(pattern4941, replacement4941)
    pattern4942 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons21, cons4, cons1714)
    def replacement4942(n, a, m, b, x, c):
        rubi.append(4942)
        return Dist(S(2)*exp(a*b*n/(m + S(1))), Int(x**m*(c*x**n)**((m + S(1))/n)/((c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x), x)
    rule4942 = ReplacementRule(pattern4942, replacement4942)
    pattern4943 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1)))), x_), cons2, cons3, cons7, cons21, cons4, cons1714)
    def replacement4943(a, n, m, b, x, c):
        rubi.append(4943)
        return Dist(S(2)*b*n*exp(a*b*n/(m + S(1)))/(m + S(1)), Int(x**m*(c*x**n)**((m + S(1))/n)/(-(c*x**n)**(S(2)*(m + S(1))/n) + exp(S(2)*a*b*n/(m + S(1)))), x), x)
    rule4943 = ReplacementRule(pattern4943, replacement4943)
    pattern4944 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1715, cons66, cons1637)
    def replacement4944(n, a, m, p, b, x, c):
        rubi.append(4944)
        return Simp(x**(m + S(1))*(p + S(-2))*sec(a + b*log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) + Simp(x**(m + S(1))*tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4944 = ReplacementRule(pattern4944, replacement4944)
    pattern4945 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1715, cons66, cons1637)
    def replacement4945(a, n, m, p, b, x, c):
        rubi.append(4945)
        return Simp(x**(m + S(1))*(p + S(-2))*csc(a + b*log(c*x**n))**(p + S(-2))/((m + S(1))*(p + S(-1))), x) - Simp(x**(m + S(1))*cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x)
    rule4945 = ReplacementRule(pattern4945, replacement4945)
    pattern4946 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons146, cons1712, cons1716)
    def replacement4946(n, a, m, p, b, x, c):
        rubi.append(4946)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*sec(a + b*log(c*x**n))**(p + S(-2)), x), x) + Simp(x**(m + S(1))*tan(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x) - Simp(x**(m + S(1))*(m + S(1))*sec(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
    rule4946 = ReplacementRule(pattern4946, replacement4946)
    pattern4947 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons146, cons1712, cons1716)
    def replacement4947(a, n, m, p, b, x, c):
        rubi.append(4947)
        return Dist((b**S(2)*n**S(2)*(p + S(-2))**S(2) + (m + S(1))**S(2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), Int(x**m*csc(a + b*log(c*x**n))**(p + S(-2)), x), x) - Simp(x**(m + S(1))*cot(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(-2))/(b*n*(p + S(-1))), x) - Simp(x**(m + S(1))*(m + S(1))*csc(a + b*log(c*x**n))**(p + S(-2))/(b**S(2)*n**S(2)*(p + S(-2))*(p + S(-1))), x)
    rule4947 = ReplacementRule(pattern4947, replacement4947)
    pattern4948 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons137, cons1708)
    def replacement4948(n, a, m, p, b, x, c):
        rubi.append(4948)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*sec(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*sec(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) - Simp(b*n*p*x**(m + S(1))*sin(a + b*log(c*x**n))*sec(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
    rule4948 = ReplacementRule(pattern4948, replacement4948)
    pattern4949 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**p_, x_), cons2, cons3, cons7, cons21, cons4, cons13, cons137, cons1708)
    def replacement4949(a, n, m, p, b, x, c):
        rubi.append(4949)
        return Dist(b**S(2)*n**S(2)*p*(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), Int(x**m*csc(a + b*log(c*x**n))**(p + S(2)), x), x) + Simp(x**(m + S(1))*(m + S(1))*csc(a + b*log(c*x**n))**p/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x) + Simp(b*n*p*x**(m + S(1))*cos(a + b*log(c*x**n))*csc(a + b*log(c*x**n))**(p + S(1))/(b**S(2)*n**S(2)*p**S(2) + (m + S(1))**S(2)), x)
    rule4949 = ReplacementRule(pattern4949, replacement4949)
    pattern4950 = Pattern(Integral(x_**WC('m', S(1))*sec(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1708)
    def replacement4950(n, a, m, p, b, x, c):
        rubi.append(4950)
        return Simp(x**(m + S(1))*((c*x**n)**(I*b)*exp(I*a)/((c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + m + S(1))/(S(2)*b*n), -(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + m + S(1)), x)
    rule4950 = ReplacementRule(pattern4950, replacement4950)
    pattern4951 = Pattern(Integral(x_**WC('m', S(1))*csc(WC('a', S(0)) + WC('b', S(1))*log(x_**WC('n', S(1))*WC('c', S(1))))**WC('p', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons1708)
    def replacement4951(a, n, m, p, b, x, c):
        rubi.append(4951)
        return Simp(x**(m + S(1))*(-I*(c*x**n)**(I*b)*exp(I*a)/(-(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(1)))**p*(-S(2)*(c*x**n)**(S(2)*I*b)*exp(S(2)*I*a) + S(2))**p*Hypergeometric2F1(p, -I*(I*b*n*p + m + S(1))/(S(2)*b*n), S(1) - I*(I*b*n*p + m + S(1))/(S(2)*b*n), (c*x**n)**(S(2)*I*b)*exp(S(2)*I*a))/(I*b*n*p + m + S(1)), x)
    rule4951 = ReplacementRule(pattern4951, replacement4951)
    pattern4952 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons163)
    def replacement4952(p, a, x, b):
        rubi.append(4952)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*sin(a*x*log(b*x)**p), x), x) - Simp(cos(a*x*log(b*x)**p)/a, x)
    rule4952 = ReplacementRule(pattern4952, replacement4952)
    pattern4953 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons13, cons163)
    def replacement4953(p, a, x, b):
        rubi.append(4953)
        return -Dist(p, Int(log(b*x)**(p + S(-1))*cos(a*x*log(b*x)**p), x), x) + Simp(sin(a*x*log(b*x)**p)/a, x)
    rule4953 = ReplacementRule(pattern4953, replacement4953)
    pattern4954 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons338, cons163)
    def replacement4954(a, n, p, b, x):
        rubi.append(4954)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) - Dist((n + S(-1))/(a*n), Int(x**(-n)*cos(a*x**n*log(b*x)**p), x), x) - Simp(x**(-n + S(1))*cos(a*x**n*log(b*x)**p)/(a*n), x)
    rule4954 = ReplacementRule(pattern4954, replacement4954)
    pattern4955 = Pattern(Integral(log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**n_*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons338, cons163)
    def replacement4955(a, n, p, b, x):
        rubi.append(4955)
        return -Dist(p/n, Int(log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) + Dist((n + S(-1))/(a*n), Int(x**(-n)*sin(a*x**n*log(b*x)**p), x), x) + Simp(x**(-n + S(1))*sin(a*x**n*log(b*x)**p)/(a*n), x)
    rule4955 = ReplacementRule(pattern4955, replacement4955)
    pattern4956 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons21, cons4, cons53, cons13, cons163)
    def replacement4956(a, n, m, p, b, x):
        rubi.append(4956)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) - Simp(cos(a*x**n*log(b*x)**p)/(a*n), x)
    rule4956 = ReplacementRule(pattern4956, replacement4956)
    pattern4957 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons21, cons4, cons53, cons13, cons163)
    def replacement4957(a, n, m, p, b, x):
        rubi.append(4957)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) + Simp(sin(a*x**n*log(b*x)**p)/(a*n), x)
    rule4957 = ReplacementRule(pattern4957, replacement4957)
    pattern4958 = Pattern(Integral(x_**WC('m', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))*sin(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons21, cons4, cons13, cons163, cons627)
    def replacement4958(a, n, m, p, b, x):
        rubi.append(4958)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*sin(a*x**n*log(b*x)**p), x), x) + Dist((m - n + S(1))/(a*n), Int(x**(m - n)*cos(a*x**n*log(b*x)**p), x), x) - Simp(x**(m - n + S(1))*cos(a*x**n*log(b*x)**p)/(a*n), x)
    rule4958 = ReplacementRule(pattern4958, replacement4958)
    pattern4959 = Pattern(Integral(x_**m_*log(x_*WC('b', S(1)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('a', S(1))*log(x_*WC('b', S(1)))**WC('p', S(1))), x_), cons2, cons3, cons21, cons4, cons13, cons163, cons627)
    def replacement4959(a, n, m, p, b, x):
        rubi.append(4959)
        return -Dist(p/n, Int(x**m*log(b*x)**(p + S(-1))*cos(a*x**n*log(b*x)**p), x), x) - Dist((m - n + S(1))/(a*n), Int(x**(m - n)*sin(a*x**n*log(b*x)**p), x), x) + Simp(x**(m - n + S(1))*sin(a*x**n*log(b*x)**p)/(a*n), x)
    rule4959 = ReplacementRule(pattern4959, replacement4959)
    pattern4960 = Pattern(Integral(sin(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons7, cons27, cons148)
    def replacement4960(a, n, x, d, c):
        rubi.append(4960)
        return -Dist(S(1)/d, Subst(Int(sin(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
    rule4960 = ReplacementRule(pattern4960, replacement4960)
    pattern4961 = Pattern(Integral(cos(WC('a', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons7, cons27, cons148)
    def replacement4961(a, n, x, d, c):
        rubi.append(4961)
        return -Dist(S(1)/d, Subst(Int(cos(a*x)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
    rule4961 = ReplacementRule(pattern4961, replacement4961)
    pattern4962 = Pattern(Integral(sin((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons148, cons71)
    def replacement4962(e, a, n, b, x, d, c):
        rubi.append(4962)
        return -Dist(S(1)/d, Subst(Int(sin(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
    rule4962 = ReplacementRule(pattern4962, replacement4962)
    pattern4963 = Pattern(Integral(cos((x_*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons148, cons71)
    def replacement4963(e, a, n, b, x, d, c):
        rubi.append(4963)
        return -Dist(S(1)/d, Subst(Int(cos(b*e/d - e*x*(-a*d + b*c)/d)**n/x**S(2), x), x, S(1)/(c + d*x)), x)
    rule4963 = ReplacementRule(pattern4963, replacement4963)
    def With4964(u, n, x):
        lst = QuotientOfLinearsParts(u, x)
        rubi.append(4964)
        return Int(sin((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
    pattern4964 = Pattern(Integral(sin(u_)**WC('n', S(1)), x_), cons148, cons1717)
    rule4964 = ReplacementRule(pattern4964, With4964)
    def With4965(u, n, x):
        lst = QuotientOfLinearsParts(u, x)
        rubi.append(4965)
        return Int(cos((x*Part(lst, S(2)) + Part(lst, S(1)))/(x*Part(lst, S(4)) + Part(lst, S(3))))**n, x)
    pattern4965 = Pattern(Integral(cos(u_)**WC('n', S(1)), x_), cons148, cons1717)
    rule4965 = ReplacementRule(pattern4965, With4965)
    pattern4966 = Pattern(Integral(WC('u', S(1))*sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1680)
    def replacement4966(u, v, p, w, x, q):
        rubi.append(4966)
        return Int(u*sin(v)**(p + q), x)
    rule4966 = ReplacementRule(pattern4966, replacement4966)
    pattern4967 = Pattern(Integral(WC('u', S(1))*cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1680)
    def replacement4967(u, v, p, w, x, q):
        rubi.append(4967)
        return Int(u*cos(v)**(p + q), x)
    rule4967 = ReplacementRule(pattern4967, replacement4967)
    pattern4968 = Pattern(Integral(sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1718, cons555)
    def replacement4968(v, p, w, x, q):
        rubi.append(4968)
        return Int(ExpandTrigReduce(sin(v)**p*sin(w)**q, x), x)
    rule4968 = ReplacementRule(pattern4968, replacement4968)
    pattern4969 = Pattern(Integral(cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1718, cons555)
    def replacement4969(v, p, w, x, q):
        rubi.append(4969)
        return Int(ExpandTrigReduce(cos(v)**p*cos(w)**q, x), x)
    rule4969 = ReplacementRule(pattern4969, replacement4969)
    pattern4970 = Pattern(Integral(x_**WC('m', S(1))*sin(v_)**WC('p', S(1))*sin(w_)**WC('q', S(1)), x_), cons1719, cons1718)
    def replacement4970(m, v, p, w, x, q):
        rubi.append(4970)
        return Int(ExpandTrigReduce(x**m, sin(v)**p*sin(w)**q, x), x)
    rule4970 = ReplacementRule(pattern4970, replacement4970)
    pattern4971 = Pattern(Integral(x_**WC('m', S(1))*cos(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1719, cons1718)
    def replacement4971(m, v, p, w, x, q):
        rubi.append(4971)
        return Int(ExpandTrigReduce(x**m, cos(v)**p*cos(w)**q, x), x)
    rule4971 = ReplacementRule(pattern4971, replacement4971)
    pattern4972 = Pattern(Integral(WC('u', S(1))*sin(v_)**WC('p', S(1))*cos(w_)**WC('p', S(1)), x_), cons1680, cons38)
    def replacement4972(u, v, p, w, x):
        rubi.append(4972)
        return Dist(S(2)**(-p), Int(u*sin(S(2)*v)**p, x), x)
    rule4972 = ReplacementRule(pattern4972, replacement4972)
    pattern4973 = Pattern(Integral(sin(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons555, cons1718)
    def replacement4973(v, p, w, x, q):
        rubi.append(4973)
        return Int(ExpandTrigReduce(sin(v)**p*cos(w)**q, x), x)
    rule4973 = ReplacementRule(pattern4973, replacement4973)
    pattern4974 = Pattern(Integral(x_**WC('m', S(1))*sin(v_)**WC('p', S(1))*cos(w_)**WC('q', S(1)), x_), cons1719, cons1718)
    def replacement4974(m, v, p, w, x, q):
        rubi.append(4974)
        return Int(ExpandTrigReduce(x**m, sin(v)**p*cos(w)**q, x), x)
    rule4974 = ReplacementRule(pattern4974, replacement4974)
    pattern4975 = Pattern(Integral(sin(v_)*tan(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4975(w, n, x, v):
        rubi.append(4975)
        return Dist(cos(v - w), Int(tan(w)**(n + S(-1))*sec(w), x), x) - Int(cos(v)*tan(w)**(n + S(-1)), x)
    rule4975 = ReplacementRule(pattern4975, replacement4975)
    pattern4976 = Pattern(Integral(cos(v_)*cot(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4976(w, n, x, v):
        rubi.append(4976)
        return Dist(cos(v - w), Int(cot(w)**(n + S(-1))*csc(w), x), x) - Int(sin(v)*cot(w)**(n + S(-1)), x)
    rule4976 = ReplacementRule(pattern4976, replacement4976)
    pattern4977 = Pattern(Integral(sin(v_)*cot(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4977(w, n, x, v):
        rubi.append(4977)
        return Dist(sin(v - w), Int(cot(w)**(n + S(-1))*csc(w), x), x) + Int(cos(v)*cot(w)**(n + S(-1)), x)
    rule4977 = ReplacementRule(pattern4977, replacement4977)
    pattern4978 = Pattern(Integral(cos(v_)*tan(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4978(w, n, x, v):
        rubi.append(4978)
        return -Dist(sin(v - w), Int(tan(w)**(n + S(-1))*sec(w), x), x) + Int(sin(v)*tan(w)**(n + S(-1)), x)
    rule4978 = ReplacementRule(pattern4978, replacement4978)
    pattern4979 = Pattern(Integral(sin(v_)*sec(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4979(w, n, x, v):
        rubi.append(4979)
        return Dist(sin(v - w), Int(sec(w)**(n + S(-1)), x), x) + Dist(cos(v - w), Int(tan(w)*sec(w)**(n + S(-1)), x), x)
    rule4979 = ReplacementRule(pattern4979, replacement4979)
    pattern4980 = Pattern(Integral(cos(v_)*csc(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4980(w, n, x, v):
        rubi.append(4980)
        return -Dist(sin(v - w), Int(csc(w)**(n + S(-1)), x), x) + Dist(cos(v - w), Int(cot(w)*csc(w)**(n + S(-1)), x), x)
    rule4980 = ReplacementRule(pattern4980, replacement4980)
    pattern4981 = Pattern(Integral(sin(v_)*csc(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4981(w, n, x, v):
        rubi.append(4981)
        return Dist(sin(v - w), Int(cot(w)*csc(w)**(n + S(-1)), x), x) + Dist(cos(v - w), Int(csc(w)**(n + S(-1)), x), x)
    rule4981 = ReplacementRule(pattern4981, replacement4981)
    pattern4982 = Pattern(Integral(cos(v_)*sec(w_)**WC('n', S(1)), x_), cons87, cons88, cons1720)
    def replacement4982(w, n, x, v):
        rubi.append(4982)
        return -Dist(sin(v - w), Int(tan(w)*sec(w)**(n + S(-1)), x), x) + Dist(cos(v - w), Int(sec(w)**(n + S(-1)), x), x)
    rule4982 = ReplacementRule(pattern4982, replacement4982)
    pattern4983 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons1352)
    def replacement4983(e, n, a, m, b, f, x, d, c):
        rubi.append(4983)
        return Int((a + b*sin(S(2)*c + S(2)*d*x)/S(2))**n*(e + f*x)**m, x)
    rule4983 = ReplacementRule(pattern4983, replacement4983)
    pattern4984 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons7, cons27, cons1470, cons150, cons168, cons463, cons1721)
    def replacement4984(a, n, m, b, x, d, c):
        rubi.append(4984)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a - b*cos(S(2)*c + S(2)*d*x) + b)**n, x), x)
    rule4984 = ReplacementRule(pattern4984, replacement4984)
    pattern4985 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons2, cons3, cons7, cons27, cons1470, cons150, cons168, cons463, cons1721)
    def replacement4985(a, n, m, b, x, d, c):
        rubi.append(4985)
        return Dist(S(2)**(-n), Int(x**m*(S(2)*a + b*cos(S(2)*c + S(2)*d*x) + b)**n, x), x)
    rule4985 = ReplacementRule(pattern4985, replacement4985)
    pattern4986 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons13)
    def replacement4986(e, a, n, m, p, b, f, x, d, c):
        rubi.append(4986)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*sin(a + b*x**n)**p, x), x, c + d*x), x)
    rule4986 = ReplacementRule(pattern4986, replacement4986)
    pattern4987 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos((c_ + x_*WC('d', S(1)))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons62, cons13)
    def replacement4987(e, a, n, m, p, b, f, x, d, c):
        rubi.append(4987)
        return Dist(d**(-m + S(-1)), Subst(Int((-c*f + d*e + f*x)**m*cos(a + b*x**n)**p, x), x, c + d*x), x)
    rule4987 = ReplacementRule(pattern4987, replacement4987)
    pattern4988 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons62, cons1470, cons1722)
    def replacement4988(e, a, m, g, b, f, x, d, c):
        rubi.append(4988)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
    rule4988 = ReplacementRule(pattern4988, replacement4988)
    pattern4989 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(b_ + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons7, cons27, cons48, cons125, cons208, cons62)
    def replacement4989(e, m, g, b, f, x, d, c):
        rubi.append(4989)
        return Dist(S(2), Int((f + g*x)**m/(b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
    rule4989 = ReplacementRule(pattern4989, replacement4989)
    pattern4990 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(WC('a', S(1))*sec(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(0)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons62, cons1470, cons1722)
    def replacement4990(e, a, m, g, b, f, x, d, c):
        rubi.append(4990)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
    rule4990 = ReplacementRule(pattern4990, replacement4990)
    pattern4991 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(c_ + WC('b', S(1))*cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons3, cons7, cons27, cons48, cons125, cons208, cons62)
    def replacement4991(e, m, g, b, f, x, d, c):
        rubi.append(4991)
        return Dist(S(2), Int((f + g*x)**m/(b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
    rule4991 = ReplacementRule(pattern4991, replacement4991)
    pattern4992 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)/(WC('a', S(1))*csc(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('b', S(1))*cot(x_*WC('e', S(1)) + WC('d', S(0)))**S(2) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons62, cons1470, cons1722)
    def replacement4992(e, a, m, g, b, f, x, d, c):
        rubi.append(4992)
        return Dist(S(2), Int((f + g*x)**m/(S(2)*a + b + c + (b - c)*cos(S(2)*d + S(2)*e*x)), x), x)
    rule4992 = ReplacementRule(pattern4992, replacement4992)
    pattern4993 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons1723)
    def replacement4993(e, a, m, b, f, x, d, c):
        rubi.append(4993)
        return Int((e + f*x)**m*exp(I*(c + d*x))/(a - I*b*exp(I*(c + d*x)) - Rt(a**S(2) - b**S(2), S(2))), x) + Int((e + f*x)**m*exp(I*(c + d*x))/(a - I*b*exp(I*(c + d*x)) + Rt(a**S(2) - b**S(2), S(2))), x) - Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
    rule4993 = ReplacementRule(pattern4993, replacement4993)
    pattern4994 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons1723)
    def replacement4994(e, a, m, b, f, x, d, c):
        rubi.append(4994)
        return -Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(a + b*exp(I*(c + d*x)) - Rt(a**S(2) - b**S(2), S(2))), x), x) - Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(a + b*exp(I*(c + d*x)) + Rt(a**S(2) - b**S(2), S(2))), x), x) + Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
    rule4994 = ReplacementRule(pattern4994, replacement4994)
    pattern4995 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons1724)
    def replacement4995(e, a, m, b, f, x, d, c):
        rubi.append(4995)
        return Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + b*exp(I*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x), x) + Dist(I, Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + b*exp(I*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x), x) - Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
    rule4995 = ReplacementRule(pattern4995, replacement4995)
    pattern4996 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons1724)
    def replacement4996(e, a, m, b, f, x, d, c):
        rubi.append(4996)
        return Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + I*b*exp(I*(c + d*x)) - Rt(-a**S(2) + b**S(2), S(2))), x) + Int((e + f*x)**m*exp(I*(c + d*x))/(I*a + I*b*exp(I*(c + d*x)) + Rt(-a**S(2) + b**S(2), S(2))), x) + Simp(I*(e + f*x)**(m + S(1))/(b*f*(m + S(1))), x)
    rule4996 = ReplacementRule(pattern4996, replacement4996)
    pattern4997 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons85, cons165, cons1257)
    def replacement4997(e, n, a, m, b, f, x, d, c):
        rubi.append(4997)
        return Dist(S(1)/a, Int((e + f*x)**m*cos(c + d*x)**(n + S(-2)), x), x) - Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)*cos(c + d*x)**(n + S(-2)), x), x)
    rule4997 = ReplacementRule(pattern4997, replacement4997)
    pattern4998 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons85, cons165, cons1257)
    def replacement4998(e, a, n, m, b, f, x, d, c):
        rubi.append(4998)
        return Dist(S(1)/a, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2)), x), x) - Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))*cos(c + d*x), x), x)
    rule4998 = ReplacementRule(pattern4998, replacement4998)
    pattern4999 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons85, cons165, cons1259)
    def replacement4999(e, n, a, m, b, f, x, d, c):
        rubi.append(4999)
        return -Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)*cos(c + d*x)**(n + S(-2)), x), x) + Dist(a/b**S(2), Int((e + f*x)**m*cos(c + d*x)**(n + S(-2)), x), x) - Dist((a**S(2) - b**S(2))/b**S(2), Int((e + f*x)**m*cos(c + d*x)**(n + S(-2))/(a + b*sin(c + d*x)), x), x)
    rule4999 = ReplacementRule(pattern4999, replacement4999)
    pattern5000 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons62, cons85, cons165, cons1259)
    def replacement5000(e, a, n, m, b, f, x, d, c):
        rubi.append(5000)
        return -Dist(S(1)/b, Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))*cos(c + d*x), x), x) + Dist(a/b**S(2), Int((e + f*x)**m*sin(c + d*x)**(n + S(-2)), x), x) - Dist((a**S(2) - b**S(2))/b**S(2), Int((e + f*x)**m*sin(c + d*x)**(n + S(-2))/(a + b*cos(c + d*x)), x), x)
    rule5000 = ReplacementRule(pattern5000, replacement5000)
    pattern5001 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons1466)
    def replacement5001(e, a, B, b, A, f, x, d, c):
        rubi.append(5001)
        return Dist(B*f/(a*d), Int(cos(c + d*x)/(a + b*sin(c + d*x)), x), x) - Simp(B*(e + f*x)*cos(c + d*x)/(a*d*(a + b*sin(c + d*x))), x)
    rule5001 = ReplacementRule(pattern5001, replacement5001)
    pattern5002 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons34, cons35, cons1466)
    def replacement5002(e, a, B, b, A, f, x, d, c):
        rubi.append(5002)
        return -Dist(B*f/(a*d), Int(sin(c + d*x)/(a + b*cos(c + d*x)), x), x) + Simp(B*(e + f*x)*sin(c + d*x)/(a*d*(a + b*cos(c + d*x))), x)
    rule5002 = ReplacementRule(pattern5002, replacement5002)
    pattern5003 = Pattern(Integral((a_ + WC('b', S(1))*tan(v_))**WC('n', S(1))*sec(v_)**WC('m', S(1)), x_), cons2, cons3, cons150, cons1543, cons1473)
    def replacement5003(n, a, m, v, b, x):
        rubi.append(5003)
        return Int((a*cos(v) + b*sin(v))**n, x)
    rule5003 = ReplacementRule(pattern5003, replacement5003)
    pattern5004 = Pattern(Integral((a_ + WC('b', S(1))*cot(v_))**WC('n', S(1))*csc(v_)**WC('m', S(1)), x_), cons2, cons3, cons150, cons1543, cons1473)
    def replacement5004(n, a, m, v, b, x):
        rubi.append(5004)
        return Int((a*sin(v) + b*cos(v))**n, x)
    rule5004 = ReplacementRule(pattern5004, replacement5004)
    pattern5005 = Pattern(Integral(WC('u', S(1))*sin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons528)
    def replacement5005(a, n, m, b, u, x, d, c):
        rubi.append(5005)
        return Int(ExpandTrigReduce(u, sin(a + b*x)**m*sin(c + d*x)**n, x), x)
    rule5005 = ReplacementRule(pattern5005, replacement5005)
    pattern5006 = Pattern(Integral(WC('u', S(1))*cos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons528)
    def replacement5006(a, n, m, b, u, x, d, c):
        rubi.append(5006)
        return Int(ExpandTrigReduce(u, cos(a + b*x)**m*cos(c + d*x)**n, x), x)
    rule5006 = ReplacementRule(pattern5006, replacement5006)
    pattern5007 = Pattern(Integral(sec(c_ + x_*WC('d', S(1)))*sec(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1725, cons71)
    def replacement5007(a, b, x, d, c):
        rubi.append(5007)
        return Dist(csc((-a*d + b*c)/b), Int(tan(c + d*x), x), x) - Dist(csc((-a*d + b*c)/d), Int(tan(a + b*x), x), x)
    rule5007 = ReplacementRule(pattern5007, replacement5007)
    pattern5008 = Pattern(Integral(csc(c_ + x_*WC('d', S(1)))*csc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1725, cons71)
    def replacement5008(a, b, x, d, c):
        rubi.append(5008)
        return Dist(csc((-a*d + b*c)/b), Int(cot(a + b*x), x), x) - Dist(csc((-a*d + b*c)/d), Int(cot(c + d*x), x), x)
    rule5008 = ReplacementRule(pattern5008, replacement5008)
    pattern5009 = Pattern(Integral(tan(c_ + x_*WC('d', S(1)))*tan(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1725, cons71)
    def replacement5009(a, b, x, d, c):
        rubi.append(5009)
        return Dist(b*cos((-a*d + b*c)/d)/d, Int(sec(a + b*x)*sec(c + d*x), x), x) - Simp(b*x/d, x)
    rule5009 = ReplacementRule(pattern5009, replacement5009)
    pattern5010 = Pattern(Integral(cot(c_ + x_*WC('d', S(1)))*cot(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1725, cons71)
    def replacement5010(a, b, x, d, c):
        rubi.append(5010)
        return Dist(cos((-a*d + b*c)/d), Int(csc(a + b*x)*csc(c + d*x), x), x) - Simp(b*x/d, x)
    rule5010 = ReplacementRule(pattern5010, replacement5010)
    pattern5011 = Pattern(Integral((WC('a', S(1))*cos(v_) + WC('b', S(1))*sin(v_))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1431)
    def replacement5011(a, n, v, b, u, x):
        rubi.append(5011)
        return Int(u*(a*exp(-a*v/b))**n, x)
    rule5011 = ReplacementRule(pattern5011, replacement5011)
    return [rule4666, rule4667, rule4668, rule4669, rule4670, rule4671, rule4672, rule4673, rule4674, rule4675, rule4676, rule4677, rule4678, rule4679, rule4680, rule4681, rule4682, rule4683, rule4684, rule4685, rule4686, rule4687, rule4688, rule4689, rule4690, rule4691, rule4692, rule4693, rule4694, rule4695, rule4696, rule4697, rule4698, rule4699, rule4700, rule4701, rule4702, rule4703, rule4704, rule4705, rule4706, rule4707, rule4708, rule4709, rule4710, rule4711, rule4712, rule4713, rule4714, rule4715, rule4716, rule4717, rule4718, rule4719, rule4720, rule4721, rule4722, rule4723, rule4724, rule4725, rule4726, rule4727, rule4728, rule4729, rule4730, rule4731, rule4732, rule4733, rule4734, rule4735, rule4736, rule4737, rule4738, rule4739, rule4740, rule4741, rule4742, rule4743, rule4744, rule4745, rule4746, rule4747, rule4748, rule4749, rule4750, rule4751, rule4752, rule4753, rule4754, rule4755, rule4756, rule4757, rule4758, rule4759, rule4760, rule4761, rule4762, rule4763, rule4764, rule4765, rule4766, rule4767, rule4768, rule4769, rule4770, rule4771, rule4772, rule4773, rule4774, rule4775, rule4776, rule4777, rule4778, rule4779, rule4780, rule4781, rule4782, rule4783, rule4784, rule4785, rule4786, rule4787, rule4788, rule4789, rule4790, rule4791, rule4792, rule4793, rule4794, rule4795, rule4796, rule4797, rule4798, rule4799, rule4800, rule4801, rule4802, rule4803, rule4804, rule4805, rule4806, rule4807, rule4808, rule4809, rule4810, rule4811, rule4812, rule4813, rule4814, rule4815, rule4816, rule4817, rule4818, rule4819, rule4820, rule4821, rule4822, rule4823, rule4824, rule4825, rule4826, rule4827, rule4828, rule4829, rule4830, rule4831, rule4832, rule4833, rule4834, rule4835, rule4836, rule4837, rule4838, rule4839, rule4840, rule4841, rule4842, rule4843, rule4844, rule4845, rule4846, rule4847, rule4848, rule4849, rule4850, rule4851, rule4852, rule4853, rule4854, rule4855, rule4856, rule4857, rule4858, rule4859, rule4860, rule4861, rule4862, rule4863, rule4864, rule4865, rule4866, rule4867, rule4868, rule4869, rule4870, rule4871, rule4872, rule4873, rule4874, rule4875, rule4876, rule4877, rule4878, rule4879, rule4880, rule4881, rule4882, rule4883, rule4884, rule4885, rule4886, rule4887, rule4888, rule4889, rule4890, rule4891, rule4892, rule4893, rule4894, rule4895, rule4896, rule4897, rule4898, rule4899, rule4900, rule4901, rule4902, rule4903, rule4904, rule4905, rule4906, rule4907, rule4908, rule4909, rule4910, rule4911, rule4912, rule4913, rule4914, rule4915, rule4916, rule4917, rule4918, rule4919, rule4920, rule4921, rule4922, rule4923, rule4924, rule4925, rule4926, rule4927, rule4928, rule4929, rule4930, rule4931, rule4932, rule4933, rule4934, rule4935, rule4936, rule4937, rule4938, rule4939, rule4940, rule4941, rule4942, rule4943, rule4944, rule4945, rule4946, rule4947, rule4948, rule4949, rule4950, rule4951, rule4952, rule4953, rule4954, rule4955, rule4956, rule4957, rule4958, rule4959, rule4960, rule4961, rule4962, rule4963, rule4964, rule4965, rule4966, rule4967, rule4968, rule4969, rule4970, rule4971, rule4972, rule4973, rule4974, rule4975, rule4976, rule4977, rule4978, rule4979, rule4980, rule4981, rule4982, rule4983, rule4984, rule4985, rule4986, rule4987, rule4988, rule4989, rule4990, rule4991, rule4992, rule4993, rule4994, rule4995, rule4996, rule4997, rule4998, rule4999, rule5000, rule5001, rule5002, rule5003, rule5004, rule5005, rule5006, rule5007, rule5008, rule5009, rule5010, rule5011, ]
