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


def inverse_trig():
    from sympy.integrals.rubi.constraints import cons89, cons90, cons2, cons3, cons8, cons91, cons1581, cons4, cons150, cons68, cons29, cons19, cons64, cons1736, cons1737, cons1738, cons1739, cons270, cons50, cons586, cons1740, cons130, cons340, cons165, cons139, cons232, cons5, cons1741, cons1742, cons1743, cons1744, cons1745, cons40, cons1746, cons1572, cons338, cons1747, cons149, cons127, cons210, cons56, cons244, cons1748, cons349, cons1749, cons488, cons963, cons95, cons96, cons164, cons274, cons1750, cons20, cons168, cons276, cons1751, cons21, cons1752, cons240, cons239, cons1753, cons248, cons1754, cons1755, cons1756, cons1757, cons1758, cons211, cons927, cons466, cons86, cons1759, cons1760, cons721, cons170, cons1761, cons669, cons1762, cons269, cons719, cons1763, cons1610, cons14, cons152, cons1200, cons1275, cons1362, cons1764, cons1765, cons36, cons37, cons38, cons1766, cons1767, cons167, cons1444, cons1768, cons1769, cons1770, cons1232, cons1771, cons1772, cons1773, cons1774, cons342, cons1775, cons1776, cons1777, cons1778, cons1045, cons87, cons33, cons1779, cons1499, cons1780, cons13, cons1781, cons1782, cons1783, cons1784, cons242, cons243, cons148, cons1785, cons1512, cons1786, cons1154, cons321, cons1787, cons1788, cons1789, cons1790, cons1791, cons1792, cons1793, cons1794, cons1795, cons1796, cons1797, cons1798, cons603, cons1799, cons263, cons1800, cons1801, cons1802, cons1803, cons1804, cons1805, cons1806, cons1807, cons179, cons119, cons1808, cons1809, cons1810, cons1811, cons1812, cons1813, cons1814, cons1815, cons1816, cons1817, cons1818, cons1819, cons1582, cons1820, cons1821, cons1822, cons1823, cons1824, cons1825, cons1826, cons1827, cons1828, cons1829, cons1830, cons1831, cons1832, cons1096, cons1833, cons1834, cons1835, cons1836, cons1837, cons385, cons1838, cons1839, cons820, cons465, cons1840, cons1841, cons1842, cons1843, cons1844, cons1845, cons1846, cons69, cons1847, cons1848, cons1849, cons1850, cons1851, cons1852, cons1853, cons554, cons1148, cons1854, cons1855, cons1244, cons1245, cons1856, cons180, cons1857, cons1858, cons1301, cons1859, cons1860


    pattern5034 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons89, cons90)
    rule5034 = ReplacementRule(pattern5034, replacement5034)

    pattern5035 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons89, cons90)
    rule5035 = ReplacementRule(pattern5035, replacement5035)

    pattern5036 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons89, cons91)
    rule5036 = ReplacementRule(pattern5036, replacement5036)

    pattern5037 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons89, cons91)
    rule5037 = ReplacementRule(pattern5037, replacement5037)

    pattern5038 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule5038 = ReplacementRule(pattern5038, replacement5038)

    pattern5039 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule5039 = ReplacementRule(pattern5039, replacement5039)

    pattern5040 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons150)
    rule5040 = ReplacementRule(pattern5040, replacement5040)

    pattern5041 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons150)
    rule5041 = ReplacementRule(pattern5041, replacement5041)

    pattern5042 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons150, cons68)
    rule5042 = ReplacementRule(pattern5042, replacement5042)

    pattern5043 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons150, cons68)
    rule5043 = ReplacementRule(pattern5043, replacement5043)

    pattern5044 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons90)
    rule5044 = ReplacementRule(pattern5044, replacement5044)

    pattern5045 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons90)
    rule5045 = ReplacementRule(pattern5045, replacement5045)

    pattern5046 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1736)
    rule5046 = ReplacementRule(pattern5046, replacement5046)

    pattern5047 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1736)
    rule5047 = ReplacementRule(pattern5047, replacement5047)

    pattern5048 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1737)
    rule5048 = ReplacementRule(pattern5048, replacement5048)

    pattern5049 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1737)
    rule5049 = ReplacementRule(pattern5049, replacement5049)

    pattern5050 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons64)
    rule5050 = ReplacementRule(pattern5050, replacement5050)

    pattern5051 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons64)
    rule5051 = ReplacementRule(pattern5051, replacement5051)

    pattern5052 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule5052 = ReplacementRule(pattern5052, replacement5052)

    pattern5053 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule5053 = ReplacementRule(pattern5053, replacement5053)

    pattern5054 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule5054 = ReplacementRule(pattern5054, replacement5054)

    pattern5055 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule5055 = ReplacementRule(pattern5055, replacement5055)

    pattern5056 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons270, cons586)
    rule5056 = ReplacementRule(pattern5056, replacement5056)

    pattern5057 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons270, cons586)
    rule5057 = ReplacementRule(pattern5057, replacement5057)

    pattern5058 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1740)
    rule5058 = ReplacementRule(pattern5058, replacement5058)

    pattern5059 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1740)
    rule5059 = ReplacementRule(pattern5059, replacement5059)

    pattern5060 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule5060 = ReplacementRule(pattern5060, With5060)

    pattern5061 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule5061 = ReplacementRule(pattern5061, With5061)

    pattern5062 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule5062 = ReplacementRule(pattern5062, replacement5062)

    pattern5063 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule5063 = ReplacementRule(pattern5063, replacement5063)

    pattern5064 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons165)
    rule5064 = ReplacementRule(pattern5064, replacement5064)

    pattern5065 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons165)
    rule5065 = ReplacementRule(pattern5065, replacement5065)

    pattern5066 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons270)
    rule5066 = ReplacementRule(pattern5066, replacement5066)

    pattern5067 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons270)
    rule5067 = ReplacementRule(pattern5067, replacement5067)

    pattern5068 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule5068 = ReplacementRule(pattern5068, replacement5068)

    pattern5069 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule5069 = ReplacementRule(pattern5069, replacement5069)

    pattern5070 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons139, cons232)
    rule5070 = ReplacementRule(pattern5070, replacement5070)

    pattern5071 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons139, cons232)
    rule5071 = ReplacementRule(pattern5071, replacement5071)

    pattern5072 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5072 = ReplacementRule(pattern5072, replacement5072)

    pattern5073 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5073 = ReplacementRule(pattern5073, replacement5073)

    pattern5074 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons91)
    rule5074 = ReplacementRule(pattern5074, replacement5074)

    pattern5075 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons91)
    rule5075 = ReplacementRule(pattern5075, replacement5075)

    pattern5076 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1741, cons1742)
    rule5076 = ReplacementRule(pattern5076, replacement5076)

    pattern5077 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1741, cons1742)
    rule5077 = ReplacementRule(pattern5077, replacement5077)

    pattern5078 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1741, cons1743)
    rule5078 = ReplacementRule(pattern5078, replacement5078)

    pattern5079 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1741, cons1743)
    rule5079 = ReplacementRule(pattern5079, replacement5079)

    pattern5080 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1744, cons1745)
    rule5080 = ReplacementRule(pattern5080, With5080)

    pattern5081 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1744, cons1745)
    rule5081 = ReplacementRule(pattern5081, With5081)

    pattern5082 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1744, cons40, cons1746)
    rule5082 = ReplacementRule(pattern5082, replacement5082)

    pattern5083 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1744, cons40, cons1746)
    rule5083 = ReplacementRule(pattern5083, replacement5083)

    pattern5084 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5084 = ReplacementRule(pattern5084, replacement5084)

    pattern5085 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5085 = ReplacementRule(pattern5085, replacement5085)

    pattern5086 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons338, cons1747, cons149)
    rule5086 = ReplacementRule(pattern5086, replacement5086)

    pattern5087 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons338, cons1747, cons149)
    rule5087 = ReplacementRule(pattern5087, replacement5087)

    pattern5088 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5088 = ReplacementRule(pattern5088, replacement5088)

    pattern5089 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5089 = ReplacementRule(pattern5089, replacement5089)

    pattern5090 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons56)
    rule5090 = ReplacementRule(pattern5090, replacement5090)

    pattern5091 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons56)
    rule5091 = ReplacementRule(pattern5091, replacement5091)

    pattern5092 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5092 = ReplacementRule(pattern5092, replacement5092)

    pattern5093 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule5093 = ReplacementRule(pattern5093, replacement5093)

    pattern5094 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons90, cons244, cons68)
    rule5094 = ReplacementRule(pattern5094, replacement5094)

    pattern5095 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons90, cons244, cons68)
    rule5095 = ReplacementRule(pattern5095, replacement5095)

    pattern5096 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule5096 = ReplacementRule(pattern5096, replacement5096)

    pattern5097 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule5097 = ReplacementRule(pattern5097, replacement5097)

    pattern5098 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons130, cons1748)
    rule5098 = ReplacementRule(pattern5098, replacement5098)

    pattern5099 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons130, cons1748)
    rule5099 = ReplacementRule(pattern5099, replacement5099)

    pattern5100 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons130)
    rule5100 = ReplacementRule(pattern5100, With5100)

    pattern5101 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons130)
    rule5101 = ReplacementRule(pattern5101, With5101)

    pattern5102 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons349, cons1749, cons488, cons270)
    rule5102 = ReplacementRule(pattern5102, With5102)

    pattern5103 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons349, cons1749, cons488, cons270)
    rule5103 = ReplacementRule(pattern5103, With5103)

    pattern5104 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons963, cons1749)
    rule5104 = ReplacementRule(pattern5104, With5104)

    pattern5105 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons963, cons1749)
    rule5105 = ReplacementRule(pattern5105, With5105)

    pattern5106 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons95, cons90, cons96)
    rule5106 = ReplacementRule(pattern5106, replacement5106)

    pattern5107 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons95, cons90, cons96)
    rule5107 = ReplacementRule(pattern5107, replacement5107)

    pattern5108 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons164, cons90, cons165, cons96)
    rule5108 = ReplacementRule(pattern5108, replacement5108)

    pattern5109 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons164, cons90, cons165, cons96)
    rule5109 = ReplacementRule(pattern5109, replacement5109)

    pattern5110 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons90, cons274, cons1750)
    rule5110 = ReplacementRule(pattern5110, replacement5110)

    pattern5111 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons90, cons274, cons1750)
    rule5111 = ReplacementRule(pattern5111, replacement5111)

    pattern5112 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons165, cons274, cons1750)
    rule5112 = ReplacementRule(pattern5112, replacement5112)

    pattern5113 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons165, cons274, cons1750)
    rule5113 = ReplacementRule(pattern5113, replacement5113)

    pattern5114 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons96, cons20)
    rule5114 = ReplacementRule(pattern5114, replacement5114)

    pattern5115 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons96, cons20)
    rule5115 = ReplacementRule(pattern5115, replacement5115)

    pattern5116 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons164, cons90, cons139, cons168)
    rule5116 = ReplacementRule(pattern5116, replacement5116)

    pattern5117 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons164, cons90, cons139, cons168)
    rule5117 = ReplacementRule(pattern5117, replacement5117)

    pattern5118 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons139, cons276, cons1751)
    rule5118 = ReplacementRule(pattern5118, replacement5118)

    pattern5119 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons139, cons276, cons1751)
    rule5119 = ReplacementRule(pattern5119, replacement5119)

    pattern5120 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons95, cons90, cons168, cons20)
    rule5120 = ReplacementRule(pattern5120, replacement5120)

    pattern5121 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons95, cons90, cons168, cons20)
    rule5121 = ReplacementRule(pattern5121, replacement5121)

    pattern5122 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270, cons150, cons20)
    rule5122 = ReplacementRule(pattern5122, replacement5122)

    pattern5123 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270, cons150, cons20)
    rule5123 = ReplacementRule(pattern5123, replacement5123)

    pattern5124 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons270, cons21)
    rule5124 = ReplacementRule(pattern5124, replacement5124)

    pattern5125 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons270, cons21)
    rule5125 = ReplacementRule(pattern5125, replacement5125)

    pattern5126 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons90, cons1740, cons1752)
    rule5126 = ReplacementRule(pattern5126, replacement5126)

    pattern5127 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons90, cons1740, cons1752)
    rule5127 = ReplacementRule(pattern5127, replacement5127)

    pattern5128 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons168, cons240, cons20)
    rule5128 = ReplacementRule(pattern5128, replacement5128)

    pattern5129 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons168, cons240, cons20)
    rule5129 = ReplacementRule(pattern5129, replacement5129)

    pattern5130 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons91, cons239)
    rule5130 = ReplacementRule(pattern5130, replacement5130)

    pattern5131 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons91, cons239)
    rule5131 = ReplacementRule(pattern5131, replacement5131)

    pattern5132 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons91, cons270)
    rule5132 = ReplacementRule(pattern5132, replacement5132)

    pattern5133 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons89, cons91, cons270)
    rule5133 = ReplacementRule(pattern5133, replacement5133)

    pattern5134 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons89, cons91, cons20, cons1753, cons1741)
    rule5134 = ReplacementRule(pattern5134, replacement5134)

    pattern5135 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons89, cons91, cons20, cons1753, cons1741)
    rule5135 = ReplacementRule(pattern5135, replacement5135)

    pattern5136 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons248, cons1754, cons64, cons1742)
    rule5136 = ReplacementRule(pattern5136, replacement5136)

    pattern5137 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons248, cons1754, cons64, cons1742)
    rule5137 = ReplacementRule(pattern5137, replacement5137)

    pattern5138 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons248, cons1754, cons64, cons1743)
    rule5138 = ReplacementRule(pattern5138, replacement5138)

    pattern5139 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons248, cons1754, cons64, cons1743)
    rule5139 = ReplacementRule(pattern5139, replacement5139)

    pattern5140 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1739, cons270, cons963, cons1755, cons20, cons1756)
    rule5140 = ReplacementRule(pattern5140, replacement5140)

    pattern5141 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1739, cons270, cons963, cons1755, cons20, cons1756)
    rule5141 = ReplacementRule(pattern5141, replacement5141)

    pattern5142 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1744, cons56)
    rule5142 = ReplacementRule(pattern5142, replacement5142)

    pattern5143 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1744, cons56)
    rule5143 = ReplacementRule(pattern5143, replacement5143)

    pattern5144 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1744, cons40, cons1757)
    rule5144 = ReplacementRule(pattern5144, With5144)

    pattern5145 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1744, cons40, cons1757)
    rule5145 = ReplacementRule(pattern5145, With5145)

    pattern5146 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1744, cons150, cons40, cons20)
    rule5146 = ReplacementRule(pattern5146, replacement5146)

    pattern5147 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1744, cons150, cons40, cons20)
    rule5147 = ReplacementRule(pattern5147, replacement5147)

    pattern5148 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1758)
    rule5148 = ReplacementRule(pattern5148, replacement5148)

    pattern5149 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1758)
    rule5149 = ReplacementRule(pattern5149, replacement5149)

    pattern5150 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons338, cons1747, cons149)
    rule5150 = ReplacementRule(pattern5150, replacement5150)

    pattern5151 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons338, cons1747, cons149)
    rule5151 = ReplacementRule(pattern5151, replacement5151)

    pattern5152 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons150)
    rule5152 = ReplacementRule(pattern5152, replacement5152)

    pattern5153 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons150)
    rule5153 = ReplacementRule(pattern5153, replacement5153)

    pattern5154 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150, cons68)
    rule5154 = ReplacementRule(pattern5154, replacement5154)

    pattern5155 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150, cons68)
    rule5155 = ReplacementRule(pattern5155, replacement5155)

    pattern5156 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons64, cons89, cons91)
    rule5156 = ReplacementRule(pattern5156, replacement5156)

    pattern5157 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons64, cons89, cons91)
    rule5157 = ReplacementRule(pattern5157, replacement5157)

    pattern5158 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons64)
    rule5158 = ReplacementRule(pattern5158, replacement5158)

    pattern5159 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons64)
    rule5159 = ReplacementRule(pattern5159, replacement5159)

    pattern5160 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons927)
    rule5160 = ReplacementRule(pattern5160, With5160)

    pattern5161 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons927)
    rule5161 = ReplacementRule(pattern5161, With5161)

    pattern5162 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons927)
    rule5162 = ReplacementRule(pattern5162, replacement5162)

    pattern5163 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons927)
    rule5163 = ReplacementRule(pattern5163, replacement5163)

    pattern5164 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons927)
    rule5164 = ReplacementRule(pattern5164, With5164)

    pattern5165 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons927)
    rule5165 = ReplacementRule(pattern5165, With5165)

    pattern5166 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons466, cons86, cons1759)
    rule5166 = ReplacementRule(pattern5166, With5166)

    pattern5167 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons466, cons86, cons1759)
    rule5167 = ReplacementRule(pattern5167, With5167)

    pattern5168 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons466, cons1760)
    rule5168 = ReplacementRule(pattern5168, With5168)

    pattern5169 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons466, cons1760)
    rule5169 = ReplacementRule(pattern5169, With5169)

    pattern5170 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons927, cons150, cons20)
    rule5170 = ReplacementRule(pattern5170, replacement5170)

    pattern5171 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons927, cons150, cons20)
    rule5171 = ReplacementRule(pattern5171, replacement5171)

    pattern5172 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons721, cons270, cons170, cons1761)
    rule5172 = ReplacementRule(pattern5172, With5172)

    pattern5173 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons721, cons270, cons170, cons1761)
    rule5173 = ReplacementRule(pattern5173, With5173)

    pattern5174 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons669, cons270, cons150, cons170, cons1762)
    rule5174 = ReplacementRule(pattern5174, replacement5174)

    pattern5175 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons669, cons270, cons150, cons170, cons1762)
    rule5175 = ReplacementRule(pattern5175, replacement5175)

    pattern5176 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons270, cons150, cons269)
    rule5176 = ReplacementRule(pattern5176, replacement5176)

    pattern5177 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons270, cons150, cons269)
    rule5177 = ReplacementRule(pattern5177, replacement5177)

    pattern5178 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons963, cons270, cons150)
    rule5178 = ReplacementRule(pattern5178, replacement5178)

    pattern5179 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons963, cons270, cons150)
    rule5179 = ReplacementRule(pattern5179, replacement5179)

    pattern5180 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons719, cons270, cons150, cons269)
    rule5180 = ReplacementRule(pattern5180, replacement5180)

    pattern5181 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons719, cons270, cons150, cons269)
    rule5181 = ReplacementRule(pattern5181, replacement5181)

    pattern5182 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons270, cons170, cons89, cons91)
    rule5182 = ReplacementRule(pattern5182, replacement5182)

    pattern5183 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons270, cons170, cons89, cons91)
    rule5183 = ReplacementRule(pattern5183, replacement5183)

    pattern5184 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1739, cons20, cons270, cons1763)
    rule5184 = ReplacementRule(pattern5184, replacement5184)

    pattern5185 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1739, cons20, cons270, cons1763)
    rule5185 = ReplacementRule(pattern5185, replacement5185)

    pattern5186 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons721, cons270, cons150)
    rule5186 = ReplacementRule(pattern5186, replacement5186)

    pattern5187 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1739, cons20, cons721, cons270, cons150)
    rule5187 = ReplacementRule(pattern5187, replacement5187)

    pattern5188 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1739, cons20, cons349, cons1740)
    rule5188 = ReplacementRule(pattern5188, replacement5188)

    pattern5189 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1739, cons20, cons349, cons1740)
    rule5189 = ReplacementRule(pattern5189, replacement5189)

    pattern5190 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons1739, cons270, cons150)
    rule5190 = ReplacementRule(pattern5190, replacement5190)

    pattern5191 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons1739, cons270, cons150)
    rule5191 = ReplacementRule(pattern5191, replacement5191)

    pattern5192 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons1739, cons349, cons1740)
    rule5192 = ReplacementRule(pattern5192, replacement5192)

    pattern5193 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons1739, cons349, cons1740)
    rule5193 = ReplacementRule(pattern5193, replacement5193)

    pattern5194 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1610)
    rule5194 = ReplacementRule(pattern5194, With5194)

    pattern5195 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1610)
    rule5195 = ReplacementRule(pattern5195, With5195)

    pattern5196 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons20)
    rule5196 = ReplacementRule(pattern5196, replacement5196)

    pattern5197 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons20)
    rule5197 = ReplacementRule(pattern5197, replacement5197)

    pattern5198 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons14, CustomConstraint(With5198))
    rule5198 = ReplacementRule(pattern5198, replacement5198)

    pattern5199 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons14, CustomConstraint(With5199))
    rule5199 = ReplacementRule(pattern5199, replacement5199)

    pattern5200 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons927, cons1739, cons349, CustomConstraint(With5200))
    rule5200 = ReplacementRule(pattern5200, replacement5200)

    pattern5201 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons927, cons1739, cons349, CustomConstraint(With5201))
    rule5201 = ReplacementRule(pattern5201, replacement5201)

    pattern5202 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons927, cons1739, cons963, cons152, CustomConstraint(With5202))
    rule5202 = ReplacementRule(pattern5202, replacement5202)

    pattern5203 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons927, cons1739, cons963, cons152, CustomConstraint(With5203))
    rule5203 = ReplacementRule(pattern5203, replacement5203)

    pattern5204 = Pattern(Integral(RFx_*asin(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons1200, cons150, CustomConstraint(With5204))
    rule5204 = ReplacementRule(pattern5204, replacement5204)

    pattern5205 = Pattern(Integral(RFx_*acos(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons1200, cons150, CustomConstraint(With5205))
    rule5205 = ReplacementRule(pattern5205, replacement5205)

    pattern5206 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons1200, cons150)
    rule5206 = ReplacementRule(pattern5206, replacement5206)

    pattern5207 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons1200, cons150)
    rule5207 = ReplacementRule(pattern5207, replacement5207)

    pattern5208 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*asin(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons29, cons50, cons1200, cons150, cons1739, cons349, CustomConstraint(With5208))
    rule5208 = ReplacementRule(pattern5208, replacement5208)

    pattern5209 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*acos(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons29, cons50, cons1200, cons150, cons1739, cons349, CustomConstraint(With5209))
    rule5209 = ReplacementRule(pattern5209, replacement5209)

    pattern5210 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons1200, cons150, cons1739, cons349)
    rule5210 = ReplacementRule(pattern5210, replacement5210)

    pattern5211 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons1200, cons150, cons1739, cons349)
    rule5211 = ReplacementRule(pattern5211, replacement5211)

    pattern5212 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5212 = ReplacementRule(pattern5212, replacement5212)

    pattern5213 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5213 = ReplacementRule(pattern5213, replacement5213)

    pattern5214 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1275)
    rule5214 = ReplacementRule(pattern5214, replacement5214)

    pattern5215 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1275)
    rule5215 = ReplacementRule(pattern5215, replacement5215)

    pattern5216 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5216 = ReplacementRule(pattern5216, replacement5216)

    pattern5217 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule5217 = ReplacementRule(pattern5217, replacement5217)

    pattern5218 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1764, cons1765)
    rule5218 = ReplacementRule(pattern5218, replacement5218)

    pattern5219 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1764, cons1765)
    rule5219 = ReplacementRule(pattern5219, replacement5219)

    pattern5220 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1764, cons1765)
    rule5220 = ReplacementRule(pattern5220, replacement5220)

    pattern5221 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1764, cons1765)
    rule5221 = ReplacementRule(pattern5221, replacement5221)

    pattern5222 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1766)
    rule5222 = ReplacementRule(pattern5222, replacement5222)

    pattern5223 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule5223 = ReplacementRule(pattern5223, replacement5223)

    pattern5224 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule5224 = ReplacementRule(pattern5224, replacement5224)

    pattern5225 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons167)
    rule5225 = ReplacementRule(pattern5225, replacement5225)

    pattern5226 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons167)
    rule5226 = ReplacementRule(pattern5226, replacement5226)

    pattern5227 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1766)
    rule5227 = ReplacementRule(pattern5227, replacement5227)

    pattern5228 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule5228 = ReplacementRule(pattern5228, replacement5228)

    pattern5229 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule5229 = ReplacementRule(pattern5229, replacement5229)

    pattern5230 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1766)
    rule5230 = ReplacementRule(pattern5230, replacement5230)

    pattern5231 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule5231 = ReplacementRule(pattern5231, replacement5231)

    pattern5232 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule5232 = ReplacementRule(pattern5232, replacement5232)

    pattern5233 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**(S(-3)/2), x_), cons2, cons3, cons8, cons29, cons1766)
    rule5233 = ReplacementRule(pattern5233, replacement5233)

    pattern5234 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-3)/2), x_), cons2, cons3, cons29, cons1767)
    rule5234 = ReplacementRule(pattern5234, replacement5234)

    pattern5235 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-3)/2), x_), cons2, cons3, cons29, cons1767)
    rule5235 = ReplacementRule(pattern5235, replacement5235)

    pattern5236 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**(S(-2)), x_), cons2, cons3, cons8, cons29, cons1766)
    rule5236 = ReplacementRule(pattern5236, replacement5236)

    pattern5237 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-2)), x_), cons2, cons3, cons29, cons1767)
    rule5237 = ReplacementRule(pattern5237, replacement5237)

    pattern5238 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-2)), x_), cons2, cons3, cons29, cons1767)
    rule5238 = ReplacementRule(pattern5238, replacement5238)

    pattern5239 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asin(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons91, cons1444)
    rule5239 = ReplacementRule(pattern5239, replacement5239)

    pattern5240 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acos(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons91, cons1444)
    rule5240 = ReplacementRule(pattern5240, replacement5240)

    pattern5241 = Pattern(Integral(asin(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), cons2, cons5, cons150)
    rule5241 = ReplacementRule(pattern5241, replacement5241)

    pattern5242 = Pattern(Integral(acos(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), cons2, cons5, cons150)
    rule5242 = ReplacementRule(pattern5242, replacement5242)

    pattern5243 = Pattern(Integral(WC('u', S(1))*asin(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5243 = ReplacementRule(pattern5243, replacement5243)

    pattern5244 = Pattern(Integral(WC('u', S(1))*acos(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5244 = ReplacementRule(pattern5244, replacement5244)

    pattern5245 = Pattern(Integral(asin(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), cons3, cons4, cons1769)
    rule5245 = ReplacementRule(pattern5245, replacement5245)

    pattern5246 = Pattern(Integral(acos(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), cons3, cons4, cons1769)
    rule5246 = ReplacementRule(pattern5246, replacement5246)

    pattern5247 = Pattern(Integral(f_**(WC('c', S(1))*asin(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons8, cons127, cons150)
    rule5247 = ReplacementRule(pattern5247, replacement5247)

    pattern5248 = Pattern(Integral(f_**(WC('c', S(1))*acos(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons8, cons127, cons150)
    rule5248 = ReplacementRule(pattern5248, replacement5248)

    pattern5249 = Pattern(Integral(asin(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons1770)
    rule5249 = ReplacementRule(pattern5249, replacement5249)

    pattern5250 = Pattern(Integral(acos(x_**S(2)*WC('a', S(1)) + sqrt(c_ + x_**S(2)*WC('d', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons1770)
    rule5250 = ReplacementRule(pattern5250, replacement5250)

    pattern5251 = Pattern(Integral(asin(u_), x_), cons1232, cons1771)
    rule5251 = ReplacementRule(pattern5251, replacement5251)

    pattern5252 = Pattern(Integral(acos(u_), x_), cons1232, cons1771)
    rule5252 = ReplacementRule(pattern5252, replacement5252)

    pattern5253 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asin(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule5253 = ReplacementRule(pattern5253, replacement5253)

    pattern5254 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acos(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule5254 = ReplacementRule(pattern5254, replacement5254)

    pattern5255 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asin(u_)), x_), cons2, cons3, cons1232, cons1773, CustomConstraint(With5255))
    rule5255 = ReplacementRule(pattern5255, replacement5255)

    pattern5256 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acos(u_)), x_), cons2, cons3, cons1232, cons1774, CustomConstraint(With5256))
    rule5256 = ReplacementRule(pattern5256, replacement5256)

    pattern5257 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons150)
    rule5257 = ReplacementRule(pattern5257, replacement5257)

    pattern5258 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons150)
    rule5258 = ReplacementRule(pattern5258, replacement5258)

    pattern5259 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons4, cons342)
    rule5259 = ReplacementRule(pattern5259, replacement5259)

    pattern5260 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons342)
    rule5260 = ReplacementRule(pattern5260, replacement5260)

    pattern5261 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150)
    rule5261 = ReplacementRule(pattern5261, replacement5261)

    pattern5262 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150)
    rule5262 = ReplacementRule(pattern5262, replacement5262)

    pattern5263 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(d_ + x_*WC('e', S(1))), x_), cons8, cons29, cons50, cons1776, cons1777)
    rule5263 = ReplacementRule(pattern5263, replacement5263)

    pattern5264 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule5264 = ReplacementRule(pattern5264, replacement5264)

    pattern5265 = Pattern(Integral(acot(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule5265 = ReplacementRule(pattern5265, replacement5265)

    pattern5266 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5266 = ReplacementRule(pattern5266, replacement5266)

    pattern5267 = Pattern(Integral((a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5267 = ReplacementRule(pattern5267, replacement5267)

    pattern5268 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5268 = ReplacementRule(pattern5268, replacement5268)

    pattern5269 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5269 = ReplacementRule(pattern5269, replacement5269)

    pattern5270 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/x_, x_), cons2, cons3, cons8, cons87, cons167)
    rule5270 = ReplacementRule(pattern5270, replacement5270)

    pattern5271 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/x_, x_), cons2, cons3, cons8, cons87, cons167)
    rule5271 = ReplacementRule(pattern5271, replacement5271)

    pattern5272 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons19, cons87, cons167, cons68)
    rule5272 = ReplacementRule(pattern5272, replacement5272)

    pattern5273 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons19, cons87, cons167, cons68)
    rule5273 = ReplacementRule(pattern5273, replacement5273)

    pattern5274 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons466)
    rule5274 = ReplacementRule(pattern5274, replacement5274)

    pattern5275 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons466)
    rule5275 = ReplacementRule(pattern5275, replacement5275)

    pattern5276 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5276 = ReplacementRule(pattern5276, replacement5276)

    pattern5277 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5277 = ReplacementRule(pattern5277, replacement5277)

    pattern5278 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150, cons33, cons170)
    rule5278 = ReplacementRule(pattern5278, replacement5278)

    pattern5279 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150, cons33, cons170)
    rule5279 = ReplacementRule(pattern5279, replacement5279)

    pattern5280 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150)
    rule5280 = ReplacementRule(pattern5280, replacement5280)

    pattern5281 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150)
    rule5281 = ReplacementRule(pattern5281, replacement5281)

    pattern5282 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150, cons33, cons96)
    rule5282 = ReplacementRule(pattern5282, replacement5282)

    pattern5283 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1775, cons150, cons33, cons96)
    rule5283 = ReplacementRule(pattern5283, replacement5283)

    pattern5284 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1779)
    rule5284 = ReplacementRule(pattern5284, replacement5284)

    pattern5285 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1779)
    rule5285 = ReplacementRule(pattern5285, replacement5285)

    pattern5286 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5286 = ReplacementRule(pattern5286, replacement5286)

    pattern5287 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5287 = ReplacementRule(pattern5287, replacement5287)

    pattern5288 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons165)
    rule5288 = ReplacementRule(pattern5288, replacement5288)

    pattern5289 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons165)
    rule5289 = ReplacementRule(pattern5289, replacement5289)

    pattern5290 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons165, cons167)
    rule5290 = ReplacementRule(pattern5290, replacement5290)

    pattern5291 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons165, cons167)
    rule5291 = ReplacementRule(pattern5291, replacement5291)

    pattern5292 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780)
    rule5292 = ReplacementRule(pattern5292, replacement5292)

    pattern5293 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1780)
    rule5293 = ReplacementRule(pattern5293, replacement5293)

    pattern5294 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons586)
    rule5294 = ReplacementRule(pattern5294, replacement5294)

    pattern5295 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons586)
    rule5295 = ReplacementRule(pattern5295, replacement5295)

    pattern5296 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270)
    rule5296 = ReplacementRule(pattern5296, replacement5296)

    pattern5297 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270)
    rule5297 = ReplacementRule(pattern5297, replacement5297)

    pattern5298 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons270)
    rule5298 = ReplacementRule(pattern5298, replacement5298)

    pattern5299 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons270)
    rule5299 = ReplacementRule(pattern5299, replacement5299)

    pattern5300 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons1740)
    rule5300 = ReplacementRule(pattern5300, replacement5300)

    pattern5301 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons1740)
    rule5301 = ReplacementRule(pattern5301, replacement5301)

    pattern5302 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5302 = ReplacementRule(pattern5302, replacement5302)

    pattern5303 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5303 = ReplacementRule(pattern5303, replacement5303)

    pattern5304 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1780)
    rule5304 = ReplacementRule(pattern5304, replacement5304)

    pattern5305 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1780)
    rule5305 = ReplacementRule(pattern5305, replacement5305)

    pattern5306 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons139, cons232)
    rule5306 = ReplacementRule(pattern5306, replacement5306)

    pattern5307 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons139, cons232)
    rule5307 = ReplacementRule(pattern5307, replacement5307)

    pattern5308 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons167)
    rule5308 = ReplacementRule(pattern5308, replacement5308)

    pattern5309 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons167)
    rule5309 = ReplacementRule(pattern5309, replacement5309)

    pattern5310 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons139, cons167, cons232)
    rule5310 = ReplacementRule(pattern5310, replacement5310)

    pattern5311 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons139, cons167, cons232)
    rule5311 = ReplacementRule(pattern5311, replacement5311)

    pattern5312 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons139, cons91)
    rule5312 = ReplacementRule(pattern5312, replacement5312)

    pattern5313 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons139, cons91)
    rule5313 = ReplacementRule(pattern5313, replacement5313)

    pattern5314 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1781, cons1742)
    rule5314 = ReplacementRule(pattern5314, replacement5314)

    pattern5315 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1781, cons1743)
    rule5315 = ReplacementRule(pattern5315, replacement5315)

    pattern5316 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1781, cons40)
    rule5316 = ReplacementRule(pattern5316, replacement5316)

    pattern5317 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1781, cons149)
    rule5317 = ReplacementRule(pattern5317, replacement5317)

    pattern5318 = Pattern(Integral(ArcTan(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule5318 = ReplacementRule(pattern5318, replacement5318)

    pattern5319 = Pattern(Integral(acot(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule5319 = ReplacementRule(pattern5319, replacement5319)

    pattern5320 = Pattern(Integral((a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5320 = ReplacementRule(pattern5320, replacement5320)

    pattern5321 = Pattern(Integral((a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule5321 = ReplacementRule(pattern5321, replacement5321)

    pattern5322 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1782)
    rule5322 = ReplacementRule(pattern5322, With5322)

    pattern5323 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1782)
    rule5323 = ReplacementRule(pattern5323, With5323)

    pattern5324 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150)
    rule5324 = ReplacementRule(pattern5324, replacement5324)

    pattern5325 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150)
    rule5325 = ReplacementRule(pattern5325, replacement5325)

    pattern5326 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5326 = ReplacementRule(pattern5326, replacement5326)

    pattern5327 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5327 = ReplacementRule(pattern5327, replacement5327)

    pattern5328 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons168)
    rule5328 = ReplacementRule(pattern5328, replacement5328)

    pattern5329 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons168)
    rule5329 = ReplacementRule(pattern5329, replacement5329)

    pattern5330 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons96)
    rule5330 = ReplacementRule(pattern5330, replacement5330)

    pattern5331 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons96)
    rule5331 = ReplacementRule(pattern5331, replacement5331)

    pattern5332 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150)
    rule5332 = ReplacementRule(pattern5332, replacement5332)

    pattern5333 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150)
    rule5333 = ReplacementRule(pattern5333, replacement5333)

    pattern5334 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons342, cons586)
    rule5334 = ReplacementRule(pattern5334, replacement5334)

    pattern5335 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons342, cons586)
    rule5335 = ReplacementRule(pattern5335, replacement5335)

    pattern5336 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons168)
    rule5336 = ReplacementRule(pattern5336, replacement5336)

    pattern5337 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons168)
    rule5337 = ReplacementRule(pattern5337, replacement5337)

    pattern5338 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5338 = ReplacementRule(pattern5338, replacement5338)

    pattern5339 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5339 = ReplacementRule(pattern5339, replacement5339)

    pattern5340 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons96)
    rule5340 = ReplacementRule(pattern5340, replacement5340)

    pattern5341 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons96)
    rule5341 = ReplacementRule(pattern5341, replacement5341)

    pattern5342 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons89, cons91)
    rule5342 = ReplacementRule(pattern5342, replacement5342)

    pattern5343 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons89, cons91)
    rule5343 = ReplacementRule(pattern5343, replacement5343)

    pattern5344 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons20, cons1783)
    rule5344 = ReplacementRule(pattern5344, replacement5344)

    pattern5345 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons20, cons1783)
    rule5345 = ReplacementRule(pattern5345, replacement5345)

    pattern5346 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons56)
    rule5346 = ReplacementRule(pattern5346, replacement5346)

    pattern5347 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons56)
    rule5347 = ReplacementRule(pattern5347, replacement5347)

    pattern5348 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons91, cons1444)
    rule5348 = ReplacementRule(pattern5348, replacement5348)

    pattern5349 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons91, cons1444)
    rule5349 = ReplacementRule(pattern5349, replacement5349)

    pattern5350 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons139, cons1784)
    rule5350 = ReplacementRule(pattern5350, replacement5350)

    pattern5351 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons13, cons139, cons1784)
    rule5351 = ReplacementRule(pattern5351, replacement5351)

    pattern5352 = Pattern(Integral(x_**S(2)*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5352 = ReplacementRule(pattern5352, replacement5352)

    pattern5353 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5353 = ReplacementRule(pattern5353, replacement5353)

    pattern5354 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons242, cons13, cons139)
    rule5354 = ReplacementRule(pattern5354, replacement5354)

    pattern5355 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons242, cons13, cons139)
    rule5355 = ReplacementRule(pattern5355, replacement5355)

    pattern5356 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons242, cons340, cons139, cons167)
    rule5356 = ReplacementRule(pattern5356, replacement5356)

    pattern5357 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons242, cons340, cons139, cons167)
    rule5357 = ReplacementRule(pattern5357, replacement5357)

    pattern5358 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1780, cons242, cons89, cons91)
    rule5358 = ReplacementRule(pattern5358, replacement5358)

    pattern5359 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1780, cons242, cons89, cons91)
    rule5359 = ReplacementRule(pattern5359, replacement5359)

    pattern5360 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1780, cons244, cons89, cons90, cons68)
    rule5360 = ReplacementRule(pattern5360, replacement5360)

    pattern5361 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1780, cons244, cons89, cons90, cons68)
    rule5361 = ReplacementRule(pattern5361, replacement5361)

    pattern5362 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons243)
    rule5362 = ReplacementRule(pattern5362, replacement5362)

    pattern5363 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons243)
    rule5363 = ReplacementRule(pattern5363, replacement5363)

    pattern5364 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons150, cons40, cons148)
    rule5364 = ReplacementRule(pattern5364, replacement5364)

    pattern5365 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons150, cons40, cons148)
    rule5365 = ReplacementRule(pattern5365, replacement5365)

    pattern5366 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons13, cons165, cons150, cons1785)
    rule5366 = ReplacementRule(pattern5366, replacement5366)

    pattern5367 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1780, cons13, cons165, cons150, cons1785)
    rule5367 = ReplacementRule(pattern5367, replacement5367)

    pattern5368 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons168)
    rule5368 = ReplacementRule(pattern5368, replacement5368)

    pattern5369 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons168)
    rule5369 = ReplacementRule(pattern5369, replacement5369)

    pattern5370 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270)
    rule5370 = ReplacementRule(pattern5370, replacement5370)

    pattern5371 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270)
    rule5371 = ReplacementRule(pattern5371, replacement5371)

    pattern5372 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons270)
    rule5372 = ReplacementRule(pattern5372, replacement5372)

    pattern5373 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons270)
    rule5373 = ReplacementRule(pattern5373, replacement5373)

    pattern5374 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons1740)
    rule5374 = ReplacementRule(pattern5374, replacement5374)

    pattern5375 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150, cons1740)
    rule5375 = ReplacementRule(pattern5375, replacement5375)

    pattern5376 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5376 = ReplacementRule(pattern5376, replacement5376)

    pattern5377 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule5377 = ReplacementRule(pattern5377, replacement5377)

    pattern5378 = Pattern(Integral(x_**m_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons96, cons1512)
    rule5378 = ReplacementRule(pattern5378, replacement5378)

    pattern5379 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons95, cons90, cons96, cons1512)
    rule5379 = ReplacementRule(pattern5379, replacement5379)

    pattern5380 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons1786, cons139, cons168, cons1154)
    rule5380 = ReplacementRule(pattern5380, replacement5380)

    pattern5381 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons1786, cons139, cons168, cons1154)
    rule5381 = ReplacementRule(pattern5381, replacement5381)

    pattern5382 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons1786, cons139, cons269, cons1154)
    rule5382 = ReplacementRule(pattern5382, replacement5382)

    pattern5383 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons1786, cons139, cons269, cons1154)
    rule5383 = ReplacementRule(pattern5383, replacement5383)

    pattern5384 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons164, cons139, cons91, cons321)
    rule5384 = ReplacementRule(pattern5384, replacement5384)

    pattern5385 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons164, cons139, cons91, cons321)
    rule5385 = ReplacementRule(pattern5385, replacement5385)

    pattern5386 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons64, cons1787, cons1742)
    rule5386 = ReplacementRule(pattern5386, replacement5386)

    pattern5387 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons64, cons1787, cons1743)
    rule5387 = ReplacementRule(pattern5387, replacement5387)

    pattern5388 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons64, cons1787, cons40)
    rule5388 = ReplacementRule(pattern5388, replacement5388)

    pattern5389 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons64, cons1787, cons149)
    rule5389 = ReplacementRule(pattern5389, replacement5389)

    pattern5390 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5390 = ReplacementRule(pattern5390, replacement5390)

    pattern5391 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5391 = ReplacementRule(pattern5391, replacement5391)

    pattern5392 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule5392 = ReplacementRule(pattern5392, With5392)

    pattern5393 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule5393 = ReplacementRule(pattern5393, With5393)

    pattern5394 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1789)
    rule5394 = ReplacementRule(pattern5394, replacement5394)

    pattern5395 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1789)
    rule5395 = ReplacementRule(pattern5395, replacement5395)

    pattern5396 = Pattern(Integral(x_**WC('m', S(1))*(a_ + ArcTan(x_*WC('c', S(1)))*WC('b', S(1)))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1790)
    rule5396 = ReplacementRule(pattern5396, replacement5396)

    pattern5397 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*acot(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1790)
    rule5397 = ReplacementRule(pattern5397, replacement5397)

    pattern5398 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5398 = ReplacementRule(pattern5398, replacement5398)

    pattern5399 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5399 = ReplacementRule(pattern5399, replacement5399)

    pattern5400 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1791)
    rule5400 = ReplacementRule(pattern5400, replacement5400)

    pattern5401 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1791)
    rule5401 = ReplacementRule(pattern5401, replacement5401)

    pattern5402 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1792)
    rule5402 = ReplacementRule(pattern5402, replacement5402)

    pattern5403 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1792)
    rule5403 = ReplacementRule(pattern5403, replacement5403)

    pattern5404 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1793)
    rule5404 = ReplacementRule(pattern5404, replacement5404)

    pattern5405 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1793)
    rule5405 = ReplacementRule(pattern5405, replacement5405)

    pattern5406 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1794)
    rule5406 = ReplacementRule(pattern5406, replacement5406)

    pattern5407 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90, cons1794)
    rule5407 = ReplacementRule(pattern5407, replacement5407)

    pattern5408 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons1791)
    rule5408 = ReplacementRule(pattern5408, replacement5408)

    pattern5409 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons1791)
    rule5409 = ReplacementRule(pattern5409, replacement5409)

    pattern5410 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons1792)
    rule5410 = ReplacementRule(pattern5410, replacement5410)

    pattern5411 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons1792)
    rule5411 = ReplacementRule(pattern5411, replacement5411)

    pattern5412 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1780)
    rule5412 = ReplacementRule(pattern5412, replacement5412)

    pattern5413 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('m', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons152, cons1795)
    rule5413 = ReplacementRule(pattern5413, replacement5413)

    pattern5414 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons152, cons1796)
    rule5414 = ReplacementRule(pattern5414, replacement5414)

    pattern5415 = Pattern(Integral(ArcTan(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons8, cons29, cons87, cons1797)
    rule5415 = ReplacementRule(pattern5415, replacement5415)

    pattern5416 = Pattern(Integral(acot(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons8, cons29, cons87, cons1797)
    rule5416 = ReplacementRule(pattern5416, replacement5416)

    pattern5417 = Pattern(Integral((ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1798)
    rule5417 = ReplacementRule(pattern5417, replacement5417)

    pattern5418 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1798)
    rule5418 = ReplacementRule(pattern5418, replacement5418)

    pattern5419 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons603)
    rule5419 = ReplacementRule(pattern5419, replacement5419)

    pattern5420 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons603)
    rule5420 = ReplacementRule(pattern5420, replacement5420)

    pattern5421 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1799)
    rule5421 = ReplacementRule(pattern5421, With5421)

    pattern5422 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1799)
    rule5422 = ReplacementRule(pattern5422, With5422)

    pattern5423 = Pattern(Integral(x_**WC('m', S(1))*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons20, cons263)
    rule5423 = ReplacementRule(pattern5423, With5423)

    pattern5424 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons20, cons263)
    rule5424 = ReplacementRule(pattern5424, With5424)

    pattern5425 = Pattern(Integral(x_*(ArcTan(x_*WC('c', S(1)))*WC('b', S(1)) + WC('a', S(0)))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1800)
    rule5425 = ReplacementRule(pattern5425, replacement5425)

    pattern5426 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acot(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1800)
    rule5426 = ReplacementRule(pattern5426, replacement5426)

    pattern5427 = Pattern(Integral(exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons1801)
    rule5427 = ReplacementRule(pattern5427, replacement5427)

    pattern5428 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons19, cons1801)
    rule5428 = ReplacementRule(pattern5428, replacement5428)

    pattern5429 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons4, cons1802)
    rule5429 = ReplacementRule(pattern5429, replacement5429)

    pattern5430 = Pattern(Integral(x_**WC('m', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons19, cons4, cons1802)
    rule5430 = ReplacementRule(pattern5430, replacement5430)

    pattern5431 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1803, cons1804)
    rule5431 = ReplacementRule(pattern5431, replacement5431)

    pattern5432 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1803, cons1805)
    rule5432 = ReplacementRule(pattern5432, replacement5432)

    pattern5433 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1806, cons40)
    rule5433 = ReplacementRule(pattern5433, replacement5433)

    pattern5434 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1806, cons149, cons1807, cons179)
    rule5434 = ReplacementRule(pattern5434, replacement5434)

    pattern5435 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1806, cons149, cons1807, cons119)
    rule5435 = ReplacementRule(pattern5435, replacement5435)

    pattern5436 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1806, cons149)
    rule5436 = ReplacementRule(pattern5436, replacement5436)

    pattern5437 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1808, cons1809)
    rule5437 = ReplacementRule(pattern5437, replacement5437)

    pattern5438 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons139, cons1809, cons1810, cons248)
    rule5438 = ReplacementRule(pattern5438, replacement5438)

    pattern5439 = Pattern(Integral(exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons8, cons29, cons4, cons1808)
    rule5439 = ReplacementRule(pattern5439, replacement5439)

    pattern5440 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1808, cons40, cons1811, cons1812)
    rule5440 = ReplacementRule(pattern5440, replacement5440)

    pattern5441 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1804)
    rule5441 = ReplacementRule(pattern5441, replacement5441)

    pattern5442 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1808, cons1805, cons1813)
    rule5442 = ReplacementRule(pattern5442, replacement5442)

    pattern5443 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1808, cons1805, cons1814)
    rule5443 = ReplacementRule(pattern5443, replacement5443)

    pattern5444 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1805)
    rule5444 = ReplacementRule(pattern5444, replacement5444)

    pattern5445 = Pattern(Integral(x_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1808, cons1809)
    rule5445 = ReplacementRule(pattern5445, replacement5445)

    pattern5446 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons139, cons1809, cons248)
    rule5446 = ReplacementRule(pattern5446, replacement5446)

    pattern5447 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1808, cons1815, cons1809)
    rule5447 = ReplacementRule(pattern5447, replacement5447)

    pattern5448 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons139, cons1809, cons1810, cons248)
    rule5448 = ReplacementRule(pattern5448, replacement5448)

    pattern5449 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1808, cons1804, cons1811, cons1812)
    rule5449 = ReplacementRule(pattern5449, replacement5449)

    pattern5450 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1808, cons1804)
    rule5450 = ReplacementRule(pattern5450, replacement5450)

    pattern5451 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1808, cons1805, cons1813)
    rule5451 = ReplacementRule(pattern5451, replacement5451)

    pattern5452 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1808, cons1805, cons1814)
    rule5452 = ReplacementRule(pattern5452, replacement5452)

    pattern5453 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1808, cons1805)
    rule5453 = ReplacementRule(pattern5453, replacement5453)

    pattern5454 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1804)
    rule5454 = ReplacementRule(pattern5454, replacement5454)

    pattern5455 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1804, cons1807)
    rule5455 = ReplacementRule(pattern5455, replacement5455)

    pattern5456 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1805, cons1816)
    rule5456 = ReplacementRule(pattern5456, replacement5456)

    pattern5457 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons1817, cons40)
    rule5457 = ReplacementRule(pattern5457, replacement5457)

    pattern5458 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1817, cons149, cons1807, cons179)
    rule5458 = ReplacementRule(pattern5458, replacement5458)

    pattern5459 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*ArcTan(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons149, cons1807, cons119)
    rule5459 = ReplacementRule(pattern5459, replacement5459)

    pattern5460 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(ArcTan(x_*WC('a', S(1)))*WC('n', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons149, cons1816)
    rule5460 = ReplacementRule(pattern5460, replacement5460)

    pattern5461 = Pattern(Integral(exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5461 = ReplacementRule(pattern5461, replacement5461)

    pattern5462 = Pattern(Integral(x_**m_*exp(n_*ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons86, cons1818, cons1819)
    rule5462 = ReplacementRule(pattern5462, replacement5462)

    pattern5463 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(ArcTan((a_ + x_*WC('b', S(1)))*WC('c', S(1)))*WC('n', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule5463 = ReplacementRule(pattern5463, replacement5463)

    pattern5464 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1820, cons1821, cons1822)
    rule5464 = ReplacementRule(pattern5464, replacement5464)

    pattern5465 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(ArcTan(a_ + x_*WC('b', S(1)))*WC('n', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1820, cons1821, cons1823)
    rule5465 = ReplacementRule(pattern5465, replacement5465)

    pattern5466 = Pattern(Integral(WC('u', S(1))*exp(ArcTan(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))*WC('n', S(1))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5466 = ReplacementRule(pattern5466, replacement5466)

    pattern5467 = Pattern(Integral(WC('u', S(1))*exp(n_*acot(x_*WC('a', S(1)))), x_), cons2, cons1807)
    rule5467 = ReplacementRule(pattern5467, replacement5467)

    pattern5468 = Pattern(Integral(exp(n_*acot(x_*WC('a', S(1)))), x_), cons2, cons1801)
    rule5468 = ReplacementRule(pattern5468, replacement5468)

    pattern5469 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*acot(x_*WC('a', S(1)))), x_), cons2, cons1801, cons20)
    rule5469 = ReplacementRule(pattern5469, replacement5469)

    pattern5470 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons4, cons1809)
    rule5470 = ReplacementRule(pattern5470, replacement5470)

    pattern5471 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons4, cons1809, cons20)
    rule5471 = ReplacementRule(pattern5471, replacement5471)

    pattern5472 = Pattern(Integral(x_**m_*exp(n_*acot(x_*WC('a', S(1)))), x_), cons2, cons19, cons1801, cons21)
    rule5472 = ReplacementRule(pattern5472, replacement5472)

    pattern5473 = Pattern(Integral(x_**m_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons19, cons4, cons1816, cons21)
    rule5473 = ReplacementRule(pattern5473, replacement5473)

    pattern5474 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1803, cons1816, cons40)
    rule5474 = ReplacementRule(pattern5474, replacement5474)

    pattern5475 = Pattern(Integral((c_ + x_*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1803, cons1816, cons149)
    rule5475 = ReplacementRule(pattern5475, replacement5475)

    pattern5476 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1806, cons1816, cons1804)
    rule5476 = ReplacementRule(pattern5476, replacement5476)

    pattern5477 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1806, cons1816, cons1804, cons20)
    rule5477 = ReplacementRule(pattern5477, replacement5477)

    pattern5478 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1806, cons1816, cons1805)
    rule5478 = ReplacementRule(pattern5478, replacement5478)

    pattern5479 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1806, cons1816, cons1804, cons21)
    rule5479 = ReplacementRule(pattern5479, replacement5479)

    pattern5480 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1806, cons1816, cons1805)
    rule5480 = ReplacementRule(pattern5480, replacement5480)

    pattern5481 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons8, cons29, cons4, cons1808)
    rule5481 = ReplacementRule(pattern5481, replacement5481)

    pattern5482 = Pattern(Integral(exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1808, cons1802)
    rule5482 = ReplacementRule(pattern5482, replacement5482)

    pattern5483 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons139, cons232, cons1810, cons1824, cons1825)
    rule5483 = ReplacementRule(pattern5483, replacement5483)

    pattern5484 = Pattern(Integral(x_*exp(WC('n', S(1))*acot(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1808, cons1802)
    rule5484 = ReplacementRule(pattern5484, replacement5484)

    pattern5485 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons1826, cons232, cons1810, cons1824, cons1825)
    rule5485 = ReplacementRule(pattern5485, replacement5485)

    pattern5486 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons1815, cons1827)
    rule5486 = ReplacementRule(pattern5486, replacement5486)

    pattern5487 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons13, cons1826, cons1828, cons1810, cons1824, cons1825)
    rule5487 = ReplacementRule(pattern5487, replacement5487)

    pattern5488 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons20, cons1829, cons40)
    rule5488 = ReplacementRule(pattern5488, replacement5488)

    pattern5489 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1808, cons1816, cons40)
    rule5489 = ReplacementRule(pattern5489, replacement5489)

    pattern5490 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1808, cons1816, cons149)
    rule5490 = ReplacementRule(pattern5490, replacement5490)

    pattern5491 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons1816, cons1804, cons1830)
    rule5491 = ReplacementRule(pattern5491, replacement5491)

    pattern5492 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons1816, cons1804, cons1831)
    rule5492 = ReplacementRule(pattern5492, replacement5492)

    pattern5493 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons1816, cons1804, cons1831, cons20)
    rule5493 = ReplacementRule(pattern5493, replacement5493)

    pattern5494 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1817, cons1816, cons1804, cons1831, cons21)
    rule5494 = ReplacementRule(pattern5494, replacement5494)

    pattern5495 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*acot(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1817, cons1816, cons1805)
    rule5495 = ReplacementRule(pattern5495, replacement5495)

    pattern5496 = Pattern(Integral(WC('u', S(1))*exp(n_*acot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons1807)
    rule5496 = ReplacementRule(pattern5496, replacement5496)

    pattern5497 = Pattern(Integral(exp(WC('n', S(1))*acot((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1816)
    rule5497 = ReplacementRule(pattern5497, replacement5497)

    pattern5498 = Pattern(Integral(x_**m_*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons86, cons1818, cons1819)
    rule5498 = ReplacementRule(pattern5498, replacement5498)

    pattern5499 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1816)
    rule5499 = ReplacementRule(pattern5499, replacement5499)

    pattern5500 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1816, cons1820, cons1821, cons1822)
    rule5500 = ReplacementRule(pattern5500, replacement5500)

    pattern5501 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acot(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1816, cons1820, cons1821, cons1823)
    rule5501 = ReplacementRule(pattern5501, replacement5501)

    pattern5502 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*acot(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule5502 = ReplacementRule(pattern5502, replacement5502)

    pattern5503 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150)
    rule5503 = ReplacementRule(pattern5503, replacement5503)

    pattern5504 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150)
    rule5504 = ReplacementRule(pattern5504, replacement5504)

    pattern5505 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons342)
    rule5505 = ReplacementRule(pattern5505, replacement5505)

    pattern5506 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons342)
    rule5506 = ReplacementRule(pattern5506, replacement5506)

    pattern5507 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule5507 = ReplacementRule(pattern5507, replacement5507)

    pattern5508 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule5508 = ReplacementRule(pattern5508, replacement5508)

    pattern5509 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons342)
    rule5509 = ReplacementRule(pattern5509, replacement5509)

    pattern5510 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons342)
    rule5510 = ReplacementRule(pattern5510, replacement5510)

    pattern5511 = Pattern(Integral((ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1832, cons1765)
    rule5511 = ReplacementRule(pattern5511, replacement5511)

    pattern5512 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1832, cons1765)
    rule5512 = ReplacementRule(pattern5512, replacement5512)

    pattern5513 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(ArcTan(c_ + x_*WC('d', S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1832, cons1765)
    rule5513 = ReplacementRule(pattern5513, replacement5513)

    pattern5514 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1832, cons1765)
    rule5514 = ReplacementRule(pattern5514, replacement5514)

    pattern5515 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons89)
    rule5515 = ReplacementRule(pattern5515, replacement5515)

    pattern5516 = Pattern(Integral(acot(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons89)
    rule5516 = ReplacementRule(pattern5516, replacement5516)

    pattern5517 = Pattern(Integral(ArcTan(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons1096)
    rule5517 = ReplacementRule(pattern5517, replacement5517)

    pattern5518 = Pattern(Integral(acot(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons1096)
    rule5518 = ReplacementRule(pattern5518, replacement5518)

    pattern5519 = Pattern(Integral(ArcTan(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons4, cons1833)
    rule5519 = ReplacementRule(pattern5519, replacement5519)

    pattern5520 = Pattern(Integral(acot(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons4, cons1833)
    rule5520 = ReplacementRule(pattern5520, replacement5520)

    pattern5521 = Pattern(Integral(ArcTan(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), cons2, cons3, cons4, cons1833)
    rule5521 = ReplacementRule(pattern5521, replacement5521)

    pattern5522 = Pattern(Integral(acot(x_**n_*WC('b', S(1)) + WC('a', S(0)))/x_, x_), cons2, cons3, cons4, cons1833)
    rule5522 = ReplacementRule(pattern5522, replacement5522)

    pattern5523 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons95, cons1834, cons1835)
    rule5523 = ReplacementRule(pattern5523, replacement5523)

    pattern5524 = Pattern(Integral(x_**WC('m', S(1))*acot(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons95, cons1834, cons1835)
    rule5524 = ReplacementRule(pattern5524, replacement5524)

    pattern5525 = Pattern(Integral(ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons1836)
    rule5525 = ReplacementRule(pattern5525, replacement5525)

    pattern5526 = Pattern(Integral(acot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons1836)
    rule5526 = ReplacementRule(pattern5526, replacement5526)

    pattern5527 = Pattern(Integral(x_**WC('m', S(1))*ArcTan(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons20, cons170)
    rule5527 = ReplacementRule(pattern5527, replacement5527)

    pattern5528 = Pattern(Integral(x_**WC('m', S(1))*acot(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons20, cons170)
    rule5528 = ReplacementRule(pattern5528, replacement5528)

    pattern5529 = Pattern(Integral(ArcTan(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5529 = ReplacementRule(pattern5529, replacement5529)

    pattern5530 = Pattern(Integral(WC('u', S(1))*acot(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5530 = ReplacementRule(pattern5530, replacement5530)

    pattern5531 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons1837)
    rule5531 = ReplacementRule(pattern5531, replacement5531)

    pattern5532 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons1837)
    rule5532 = ReplacementRule(pattern5532, replacement5532)

    pattern5533 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons19, cons1837, cons68)
    rule5533 = ReplacementRule(pattern5533, replacement5533)

    pattern5534 = Pattern(Integral(acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons19, cons1837, cons68)
    rule5534 = ReplacementRule(pattern5534, replacement5534)

    pattern5535 = Pattern(Integral(ArcTan(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1837, cons385)
    rule5535 = ReplacementRule(pattern5535, replacement5535)

    pattern5536 = Pattern(Integral(acot(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1837, cons385)
    rule5536 = ReplacementRule(pattern5536, replacement5536)

    pattern5537 = Pattern(Integral(ArcTan(v_ + sqrt(w_)*WC('s', S(1)))*WC('u', S(1)), x_), cons1838, cons1839)
    rule5537 = ReplacementRule(pattern5537, replacement5537)

    pattern5538 = Pattern(Integral(WC('u', S(1))*acot(v_ + sqrt(w_)*WC('s', S(1))), x_), cons1838, cons1839)
    rule5538 = ReplacementRule(pattern5538, replacement5538)

    pattern5539 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), cons820, cons87, cons465, cons1840, cons1841, CustomConstraint(With5539))
    rule5539 = ReplacementRule(pattern5539, replacement5539)

    pattern5540 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), cons820, cons87, cons465, cons1840, cons1842, CustomConstraint(With5540))
    rule5540 = ReplacementRule(pattern5540, replacement5540)

    pattern5541 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1843)
    rule5541 = ReplacementRule(pattern5541, replacement5541)

    pattern5542 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1843)
    rule5542 = ReplacementRule(pattern5542, replacement5542)

    pattern5543 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1844)
    rule5543 = ReplacementRule(pattern5543, replacement5543)

    pattern5544 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1844)
    rule5544 = ReplacementRule(pattern5544, replacement5544)

    pattern5545 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1845)
    rule5545 = ReplacementRule(pattern5545, replacement5545)

    pattern5546 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1845)
    rule5546 = ReplacementRule(pattern5546, replacement5546)

    pattern5547 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1845)
    rule5547 = ReplacementRule(pattern5547, replacement5547)

    pattern5548 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1846)
    rule5548 = ReplacementRule(pattern5548, replacement5548)

    pattern5549 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1843)
    rule5549 = ReplacementRule(pattern5549, replacement5549)

    pattern5550 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1843)
    rule5550 = ReplacementRule(pattern5550, replacement5550)

    pattern5551 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1844)
    rule5551 = ReplacementRule(pattern5551, replacement5551)

    pattern5552 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1844)
    rule5552 = ReplacementRule(pattern5552, replacement5552)

    pattern5553 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1845)
    rule5553 = ReplacementRule(pattern5553, replacement5553)

    pattern5554 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1845)
    rule5554 = ReplacementRule(pattern5554, replacement5554)

    pattern5555 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1846)
    rule5555 = ReplacementRule(pattern5555, replacement5555)

    pattern5556 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1846)
    rule5556 = ReplacementRule(pattern5556, replacement5556)

    pattern5557 = Pattern(Integral(ArcTan(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule5557 = ReplacementRule(pattern5557, replacement5557)

    pattern5558 = Pattern(Integral(acot(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule5558 = ReplacementRule(pattern5558, replacement5558)

    pattern5559 = Pattern(Integral(ArcTan(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule5559 = ReplacementRule(pattern5559, replacement5559)

    pattern5560 = Pattern(Integral(acot(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule5560 = ReplacementRule(pattern5560, replacement5560)

    pattern5561 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5561 = ReplacementRule(pattern5561, replacement5561)

    pattern5562 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5562 = ReplacementRule(pattern5562, replacement5562)

    pattern5563 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5563 = ReplacementRule(pattern5563, replacement5563)

    pattern5564 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(S(1)/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule5564 = ReplacementRule(pattern5564, replacement5564)

    pattern5565 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1847)
    rule5565 = ReplacementRule(pattern5565, replacement5565)

    pattern5566 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1847)
    rule5566 = ReplacementRule(pattern5566, replacement5566)

    pattern5567 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1847)
    rule5567 = ReplacementRule(pattern5567, replacement5567)

    pattern5568 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1847)
    rule5568 = ReplacementRule(pattern5568, replacement5568)

    pattern5569 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1848)
    rule5569 = ReplacementRule(pattern5569, replacement5569)

    pattern5570 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1848)
    rule5570 = ReplacementRule(pattern5570, replacement5570)

    pattern5571 = Pattern(Integral(ArcTan(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1848)
    rule5571 = ReplacementRule(pattern5571, replacement5571)

    pattern5572 = Pattern(Integral(acot(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1848)
    rule5572 = ReplacementRule(pattern5572, replacement5572)

    pattern5573 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1847)
    rule5573 = ReplacementRule(pattern5573, replacement5573)

    pattern5574 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1847)
    rule5574 = ReplacementRule(pattern5574, replacement5574)

    pattern5575 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1847)
    rule5575 = ReplacementRule(pattern5575, replacement5575)

    pattern5576 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1847)
    rule5576 = ReplacementRule(pattern5576, replacement5576)

    pattern5577 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1848)
    rule5577 = ReplacementRule(pattern5577, replacement5577)

    pattern5578 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1848)
    rule5578 = ReplacementRule(pattern5578, replacement5578)

    pattern5579 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*ArcTan(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1848)
    rule5579 = ReplacementRule(pattern5579, replacement5579)

    pattern5580 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acot(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1848)
    rule5580 = ReplacementRule(pattern5580, replacement5580)

    pattern5581 = Pattern(Integral(ArcTan(u_), x_), cons1232)
    rule5581 = ReplacementRule(pattern5581, replacement5581)

    pattern5582 = Pattern(Integral(acot(u_), x_), cons1232)
    rule5582 = ReplacementRule(pattern5582, replacement5582)

    pattern5583 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(ArcTan(u_)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1849)
    rule5583 = ReplacementRule(pattern5583, replacement5583)

    pattern5584 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acot(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1849)
    rule5584 = ReplacementRule(pattern5584, replacement5584)

    pattern5585 = Pattern(Integral(v_*(ArcTan(u_)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons1232, cons1850, cons1851, CustomConstraint(With5585))
    rule5585 = ReplacementRule(pattern5585, replacement5585)

    pattern5586 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acot(u_)), x_), cons2, cons3, cons1232, cons1852, cons1853, CustomConstraint(With5586))
    rule5586 = ReplacementRule(pattern5586, replacement5586)

    pattern5587 = Pattern(Integral(ArcTan(v_)*log(w_)/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons554, cons1148, cons1854, cons1855)
    rule5587 = ReplacementRule(pattern5587, replacement5587)

    pattern5588 = Pattern(Integral(ArcTan(v_)*log(w_), x_), cons1244, cons1245)
    rule5588 = ReplacementRule(pattern5588, replacement5588)

    pattern5589 = Pattern(Integral(log(w_)*acot(v_), x_), cons1244, cons1245)
    rule5589 = ReplacementRule(pattern5589, replacement5589)

    pattern5590 = Pattern(Integral(u_*ArcTan(v_)*log(w_), x_), cons1244, cons1245, CustomConstraint(With5590))
    rule5590 = ReplacementRule(pattern5590, replacement5590)

    pattern5591 = Pattern(Integral(u_*log(w_)*acot(v_), x_), cons1244, cons1245, CustomConstraint(With5591))
    rule5591 = ReplacementRule(pattern5591, replacement5591)

    pattern5592 = Pattern(Integral(asec(x_*WC('c', S(1))), x_), cons8, cons8)
    rule5592 = ReplacementRule(pattern5592, replacement5592)

    pattern5593 = Pattern(Integral(acsc(x_*WC('c', S(1))), x_), cons8, cons8)
    rule5593 = ReplacementRule(pattern5593, replacement5593)

    pattern5594 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule5594 = ReplacementRule(pattern5594, replacement5594)

    pattern5595 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule5595 = ReplacementRule(pattern5595, replacement5595)

    pattern5596 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons14)
    rule5596 = ReplacementRule(pattern5596, replacement5596)

    pattern5597 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons14)
    rule5597 = ReplacementRule(pattern5597, replacement5597)

    pattern5598 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons68)
    rule5598 = ReplacementRule(pattern5598, replacement5598)

    pattern5599 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons68)
    rule5599 = ReplacementRule(pattern5599, replacement5599)

    pattern5600 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons20)
    rule5600 = ReplacementRule(pattern5600, replacement5600)

    pattern5601 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons20)
    rule5601 = ReplacementRule(pattern5601, replacement5601)

    pattern5602 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons1856)
    rule5602 = ReplacementRule(pattern5602, replacement5602)

    pattern5603 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons1856)
    rule5603 = ReplacementRule(pattern5603, replacement5603)

    pattern5604 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1745)
    rule5604 = ReplacementRule(pattern5604, With5604)

    pattern5605 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1745)
    rule5605 = ReplacementRule(pattern5605, With5605)

    pattern5606 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons40)
    rule5606 = ReplacementRule(pattern5606, replacement5606)

    pattern5607 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons40)
    rule5607 = ReplacementRule(pattern5607, replacement5607)

    pattern5608 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons180, cons1857)
    rule5608 = ReplacementRule(pattern5608, replacement5608)

    pattern5609 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons180, cons1857)
    rule5609 = ReplacementRule(pattern5609, replacement5609)

    pattern5610 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons1858)
    rule5610 = ReplacementRule(pattern5610, replacement5610)

    pattern5611 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons1858)
    rule5611 = ReplacementRule(pattern5611, replacement5611)

    pattern5612 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5612 = ReplacementRule(pattern5612, replacement5612)

    pattern5613 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule5613 = ReplacementRule(pattern5613, replacement5613)

    pattern5614 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5614 = ReplacementRule(pattern5614, replacement5614)

    pattern5615 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule5615 = ReplacementRule(pattern5615, replacement5615)

    pattern5616 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule5616 = ReplacementRule(pattern5616, With5616)

    pattern5617 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule5617 = ReplacementRule(pattern5617, With5617)

    pattern5618 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1301)
    rule5618 = ReplacementRule(pattern5618, replacement5618)

    pattern5619 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1301)
    rule5619 = ReplacementRule(pattern5619, replacement5619)

    pattern5620 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons180, cons1857)
    rule5620 = ReplacementRule(pattern5620, replacement5620)

    pattern5621 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons180, cons1857)
    rule5621 = ReplacementRule(pattern5621, replacement5621)

    pattern5622 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons1858)
    rule5622 = ReplacementRule(pattern5622, replacement5622)

    pattern5623 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons1858)
    rule5623 = ReplacementRule(pattern5623, replacement5623)

    pattern5624 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5624 = ReplacementRule(pattern5624, replacement5624)

    pattern5625 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule5625 = ReplacementRule(pattern5625, replacement5625)

    pattern5626 = Pattern(Integral(asec(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule5626 = ReplacementRule(pattern5626, replacement5626)

    pattern5627 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule5627 = ReplacementRule(pattern5627, replacement5627)

    pattern5628 = Pattern(Integral(asec(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons1833)
    rule5628 = ReplacementRule(pattern5628, replacement5628)

    pattern5629 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons1833)
    rule5629 = ReplacementRule(pattern5629, replacement5629)

    pattern5630 = Pattern(Integral(asec(a_ + x_*WC('b', S(1)))/x_, x_), cons2, cons3, cons69)
    rule5630 = ReplacementRule(pattern5630, replacement5630)

    pattern5631 = Pattern(Integral(acsc(a_ + x_*WC('b', S(1)))/x_, x_), cons2, cons3, cons69)
    rule5631 = ReplacementRule(pattern5631, replacement5631)

    pattern5632 = Pattern(Integral(x_**WC('m', S(1))*asec(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons20, cons68)
    rule5632 = ReplacementRule(pattern5632, replacement5632)

    pattern5633 = Pattern(Integral(x_**WC('m', S(1))*acsc(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons20, cons68)
    rule5633 = ReplacementRule(pattern5633, replacement5633)

    pattern5634 = Pattern(Integral(x_**WC('m', S(1))*asec(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons64)
    rule5634 = ReplacementRule(pattern5634, replacement5634)

    pattern5635 = Pattern(Integral(x_**WC('m', S(1))*acsc(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons64)
    rule5635 = ReplacementRule(pattern5635, replacement5635)

    pattern5636 = Pattern(Integral(WC('u', S(1))*asec(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5636 = ReplacementRule(pattern5636, replacement5636)

    pattern5637 = Pattern(Integral(WC('u', S(1))*acsc(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule5637 = ReplacementRule(pattern5637, replacement5637)

    pattern5638 = Pattern(Integral(f_**(WC('c', S(1))*asec(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons8, cons127, cons150)
    rule5638 = ReplacementRule(pattern5638, replacement5638)

    pattern5639 = Pattern(Integral(f_**(WC('c', S(1))*acsc(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons8, cons127, cons150)
    rule5639 = ReplacementRule(pattern5639, replacement5639)

    pattern5640 = Pattern(Integral(asec(u_), x_), cons1232, cons1771)
    rule5640 = ReplacementRule(pattern5640, replacement5640)

    pattern5641 = Pattern(Integral(acsc(u_), x_), cons1232, cons1771)
    rule5641 = ReplacementRule(pattern5641, replacement5641)

    pattern5642 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asec(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule5642 = ReplacementRule(pattern5642, replacement5642)

    pattern5643 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsc(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule5643 = ReplacementRule(pattern5643, replacement5643)

    pattern5644 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asec(u_)), x_), cons2, cons3, cons1232, cons1859, CustomConstraint(With5644))
    rule5644 = ReplacementRule(pattern5644, replacement5644)

    pattern5645 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acsc(u_)), x_), cons2, cons3, cons1232, cons1860, CustomConstraint(With5645))
    rule5645 = ReplacementRule(pattern5645, replacement5645)
    return [rule5034, rule5035, rule5036, rule5037, rule5038, rule5039, rule5040, rule5041, rule5042, rule5043, rule5044, rule5045, rule5046, rule5047, rule5048, rule5049, rule5050, rule5051, rule5052, rule5053, rule5054, rule5055, rule5056, rule5057, rule5058, rule5059, rule5060, rule5061, rule5062, rule5063, rule5064, rule5065, rule5066, rule5067, rule5068, rule5069, rule5070, rule5071, rule5072, rule5073, rule5074, rule5075, rule5076, rule5077, rule5078, rule5079, rule5080, rule5081, rule5082, rule5083, rule5084, rule5085, rule5086, rule5087, rule5088, rule5089, rule5090, rule5091, rule5092, rule5093, rule5094, rule5095, rule5096, rule5097, rule5098, rule5099, rule5100, rule5101, rule5102, rule5103, rule5104, rule5105, rule5106, rule5107, rule5108, rule5109, rule5110, rule5111, rule5112, rule5113, rule5114, rule5115, rule5116, rule5117, rule5118, rule5119, rule5120, rule5121, rule5122, rule5123, rule5124, rule5125, rule5126, rule5127, rule5128, rule5129, rule5130, rule5131, rule5132, rule5133, rule5134, rule5135, rule5136, rule5137, rule5138, rule5139, rule5140, rule5141, rule5142, rule5143, rule5144, rule5145, rule5146, rule5147, rule5148, rule5149, rule5150, rule5151, rule5152, rule5153, rule5154, rule5155, rule5156, rule5157, rule5158, rule5159, rule5160, rule5161, rule5162, rule5163, rule5164, rule5165, rule5166, rule5167, rule5168, rule5169, rule5170, rule5171, rule5172, rule5173, rule5174, rule5175, rule5176, rule5177, rule5178, rule5179, rule5180, rule5181, rule5182, rule5183, rule5184, rule5185, rule5186, rule5187, rule5188, rule5189, rule5190, rule5191, rule5192, rule5193, rule5194, rule5195, rule5196, rule5197, rule5198, rule5199, rule5200, rule5201, rule5202, rule5203, rule5204, rule5205, rule5206, rule5207, rule5208, rule5209, rule5210, rule5211, rule5212, rule5213, rule5214, rule5215, rule5216, rule5217, rule5218, rule5219, rule5220, rule5221, rule5222, rule5223, rule5224, rule5225, rule5226, rule5227, rule5228, rule5229, rule5230, rule5231, rule5232, rule5233, rule5234, rule5235, rule5236, rule5237, rule5238, rule5239, rule5240, rule5241, rule5242, rule5243, rule5244, rule5245, rule5246, rule5247, rule5248, rule5249, rule5250, rule5251, rule5252, rule5253, rule5254, rule5255, rule5256, rule5257, rule5258, rule5259, rule5260, rule5261, rule5262, rule5263, rule5264, rule5265, rule5266, rule5267, rule5268, rule5269, rule5270, rule5271, rule5272, rule5273, rule5274, rule5275, rule5276, rule5277, rule5278, rule5279, rule5280, rule5281, rule5282, rule5283, rule5284, rule5285, rule5286, rule5287, rule5288, rule5289, rule5290, rule5291, rule5292, rule5293, rule5294, rule5295, rule5296, rule5297, rule5298, rule5299, rule5300, rule5301, rule5302, rule5303, rule5304, rule5305, rule5306, rule5307, rule5308, rule5309, rule5310, rule5311, rule5312, rule5313, rule5314, rule5315, rule5316, rule5317, rule5318, rule5319, rule5320, rule5321, rule5322, rule5323, rule5324, rule5325, rule5326, rule5327, rule5328, rule5329, rule5330, rule5331, rule5332, rule5333, rule5334, rule5335, rule5336, rule5337, rule5338, rule5339, rule5340, rule5341, rule5342, rule5343, rule5344, rule5345, rule5346, rule5347, rule5348, rule5349, rule5350, rule5351, rule5352, rule5353, rule5354, rule5355, rule5356, rule5357, rule5358, rule5359, rule5360, rule5361, rule5362, rule5363, rule5364, rule5365, rule5366, rule5367, rule5368, rule5369, rule5370, rule5371, rule5372, rule5373, rule5374, rule5375, rule5376, rule5377, rule5378, rule5379, rule5380, rule5381, rule5382, rule5383, rule5384, rule5385, rule5386, rule5387, rule5388, rule5389, rule5390, rule5391, rule5392, rule5393, rule5394, rule5395, rule5396, rule5397, rule5398, rule5399, rule5400, rule5401, rule5402, rule5403, rule5404, rule5405, rule5406, rule5407, rule5408, rule5409, rule5410, rule5411, rule5412, rule5413, rule5414, rule5415, rule5416, rule5417, rule5418, rule5419, rule5420, rule5421, rule5422, rule5423, rule5424, rule5425, rule5426, rule5427, rule5428, rule5429, rule5430, rule5431, rule5432, rule5433, rule5434, rule5435, rule5436, rule5437, rule5438, rule5439, rule5440, rule5441, rule5442, rule5443, rule5444, rule5445, rule5446, rule5447, rule5448, rule5449, rule5450, rule5451, rule5452, rule5453, rule5454, rule5455, rule5456, rule5457, rule5458, rule5459, rule5460, rule5461, rule5462, rule5463, rule5464, rule5465, rule5466, rule5467, rule5468, rule5469, rule5470, rule5471, rule5472, rule5473, rule5474, rule5475, rule5476, rule5477, rule5478, rule5479, rule5480, rule5481, rule5482, rule5483, rule5484, rule5485, rule5486, rule5487, rule5488, rule5489, rule5490, rule5491, rule5492, rule5493, rule5494, rule5495, rule5496, rule5497, rule5498, rule5499, rule5500, rule5501, rule5502, rule5503, rule5504, rule5505, rule5506, rule5507, rule5508, rule5509, rule5510, rule5511, rule5512, rule5513, rule5514, rule5515, rule5516, rule5517, rule5518, rule5519, rule5520, rule5521, rule5522, rule5523, rule5524, rule5525, rule5526, rule5527, rule5528, rule5529, rule5530, rule5531, rule5532, rule5533, rule5534, rule5535, rule5536, rule5537, rule5538, rule5539, rule5540, rule5541, rule5542, rule5543, rule5544, rule5545, rule5546, rule5547, rule5548, rule5549, rule5550, rule5551, rule5552, rule5553, rule5554, rule5555, rule5556, rule5557, rule5558, rule5559, rule5560, rule5561, rule5562, rule5563, rule5564, rule5565, rule5566, rule5567, rule5568, rule5569, rule5570, rule5571, rule5572, rule5573, rule5574, rule5575, rule5576, rule5577, rule5578, rule5579, rule5580, rule5581, rule5582, rule5583, rule5584, rule5585, rule5586, rule5587, rule5588, rule5589, rule5590, rule5591, rule5592, rule5593, rule5594, rule5595, rule5596, rule5597, rule5598, rule5599, rule5600, rule5601, rule5602, rule5603, rule5604, rule5605, rule5606, rule5607, rule5608, rule5609, rule5610, rule5611, rule5612, rule5613, rule5614, rule5615, rule5616, rule5617, rule5618, rule5619, rule5620, rule5621, rule5622, rule5623, rule5624, rule5625, rule5626, rule5627, rule5628, rule5629, rule5630, rule5631, rule5632, rule5633, rule5634, rule5635, rule5636, rule5637, rule5638, rule5639, rule5640, rule5641, rule5642, rule5643, rule5644, rule5645, ]





def replacement5034(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*asin(c*x))**n, x)


def replacement5035(a, b, c, n, x):
    return Dist(b*c*n, Int(x*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*acos(c*x))**n, x)


def replacement5036(a, b, c, n, x):
    return Dist(c/(b*(n + S(1))), Int(x*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5037(a, b, c, n, x):
    return -Dist(c/(b*(n + S(1))), Int(x*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5038(a, b, c, n, x):
    return Dist(S(1)/(b*c), Subst(Int(x**n*cos(a/b - x/b), x), x, a + b*asin(c*x)), x)


def replacement5039(a, b, c, n, x):
    return Dist(S(1)/(b*c), Subst(Int(x**n*sin(a/b - x/b), x), x, a + b*acos(c*x)), x)


def replacement5040(a, b, c, n, x):
    return Subst(Int((a + b*x)**n/tan(x), x), x, asin(c*x))


def replacement5041(a, b, c, n, x):
    return -Subst(Int((a + b*x)**n*tan(x), x), x, acos(c*x))


def replacement5042(a, b, c, d, m, n, x):
    return -Dist(b*c*n/(d*(m + S(1))), Int((d*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*asin(c*x))**n/(d*(m + S(1))), x)


def replacement5043(a, b, c, d, m, n, x):
    return Dist(b*c*n/(d*(m + S(1))), Int((d*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*acos(c*x))**n/(d*(m + S(1))), x)


def replacement5044(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*asin(c*x))**n/(m + S(1)), x)


def replacement5045(a, b, c, m, n, x):
    return Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*acos(c*x))**n/(m + S(1)), x)


def replacement5046(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1))/(b*(n + S(1))), Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*sin(x)**S(2))*sin(x)**(m + S(-1)), x), x), x, asin(c*x)), x) + Simp(x**m*(a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5047(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1))/(b*(n + S(1))), Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m - (m + S(1))*cos(x)**S(2))*cos(x)**(m + S(-1)), x), x), x, acos(c*x)), x) - Simp(x**m*(a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5048(a, b, c, m, n, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(c*(m + S(1))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*asin(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**m*(a + b*asin(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5049(a, b, c, m, n, x):
    return Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(c*(m + S(1))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*acos(c*x))**(n + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Simp(x**m*(a + b*acos(c*x))**(n + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5050(a, b, c, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*sin(x)**m*cos(x), x), x, asin(c*x)), x)


def replacement5051(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*sin(x)*cos(x)**m, x), x, acos(c*x)), x)


def replacement5052(a, b, c, d, m, n, x):
    return Int((d*x)**m*(a + b*asin(c*x))**n, x)


def replacement5053(a, b, c, d, m, n, x):
    return Int((d*x)**m*(a + b*acos(c*x))**n, x)


def replacement5054(a, b, c, d, e, x):
    return Simp(log(a + b*asin(c*x))/(b*c*sqrt(d)), x)


def replacement5055(a, b, c, d, e, x):
    return -Simp(log(a + b*acos(c*x))/(b*c*sqrt(d)), x)


def replacement5056(a, b, c, d, e, n, x):
    return Simp((a + b*asin(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5057(a, b, c, d, e, n, x):
    return -Simp((a + b*acos(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5058(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def replacement5059(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def With5060(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5061(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5062(a, b, c, d, e, n, x):
    return Dist(sqrt(d + e*x**S(2))/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))), Int((a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))), Int(x*(a + b*asin(c*x))**(n + S(-1)), x), x) + Simp(x*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/S(2), x)


def replacement5063(a, b, c, d, e, n, x):
    return Dist(sqrt(d + e*x**S(2))/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))), Int((a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(b*c*n*sqrt(d + e*x**S(2))/(S(2)*sqrt(-c**S(2)*x**S(2) + S(1))), Int(x*(a + b*acos(c*x))**(n + S(-1)), x), x) + Simp(x*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/S(2), x)


def replacement5064(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*p + S(1)), Int(x*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x)


def replacement5065(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*p + S(1)), Int(x*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x)


def replacement5066(a, b, c, d, e, n, x):
    return -Dist(b*c*n/sqrt(d), Int(x*(a + b*asin(c*x))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x*(a + b*asin(c*x))**n/(d*sqrt(d + e*x**S(2))), x)


def replacement5067(a, b, c, d, e, n, x):
    return Dist(b*c*n/sqrt(d), Int(x*(a + b*acos(c*x))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp(x*(a + b*acos(c*x))**n/(d*sqrt(d + e*x**S(2))), x)


def replacement5068(a, b, c, d, e, n, x):
    return -Dist(b*c*n*sqrt(-c**S(2)*x**S(2) + S(1))/(d*sqrt(d + e*x**S(2))), Int(x*(a + b*asin(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*asin(c*x))**n/(d*sqrt(d + e*x**S(2))), x)


def replacement5069(a, b, c, d, e, n, x):
    return Dist(b*c*n*sqrt(-c**S(2)*x**S(2) + S(1))/(d*sqrt(d + e*x**S(2))), Int(x*(a + b*acos(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*acos(c*x))**n/(d*sqrt(d + e*x**S(2))), x)


def replacement5070(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*(p + S(1))), Int(x*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp(x*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement5071(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*(p + S(1))), Int(x*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp(x*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement5072(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*d), Subst(Int((a + b*x)**n/cos(x), x), x, asin(c*x)), x)


def replacement5073(a, b, c, d, e, n, x):
    return -Dist(S(1)/(c*d), Subst(Int((a + b*x)**n/sin(x), x), x, acos(c*x)), x)


def replacement5074(a, b, c, d, e, n, p, x):
    return Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*(n + S(1))), Int(x*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5075(a, b, c, d, e, n, p, x):
    return -Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*(n + S(1))), Int(x*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5076(a, b, c, d, e, n, p, x):
    return Dist(d**p/c, Subst(Int((a + b*x)**n*cos(x)**(S(2)*p + S(1)), x), x, asin(c*x)), x)


def replacement5077(a, b, c, d, e, n, p, x):
    return -Dist(d**p/c, Subst(Int((a + b*x)**n*sin(x)**(S(2)*p + S(1)), x), x, acos(c*x)), x)


def replacement5078(a, b, c, d, e, n, p, x):
    return Dist(d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(-c**S(2)*x**S(2) + S(1)), Int((a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5079(a, b, c, d, e, n, p, x):
    return Dist(d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(-c**S(2)*x**S(2) + S(1)), Int((a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def With5080(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5081(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5082(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n, (d + e*x**S(2))**p, x), x)


def replacement5083(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n, (d + e*x**S(2))**p, x), x)


def replacement5084(a, b, c, d, e, n, p, x):
    return Int((a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5085(a, b, c, d, e, n, p, x):
    return Int((a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5086(a, b, c, d, e, f, g, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((a + b*asin(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement5087(a, b, c, d, e, f, g, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((a + b*acos(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement5088(a, b, c, d, e, n, x):
    return -Dist(S(1)/e, Subst(Int((a + b*x)**n*tan(x), x), x, asin(c*x)), x)


def replacement5089(a, b, c, d, e, n, x):
    return Dist(S(1)/e, Subst(Int((a + b*x)**n/tan(x), x), x, acos(c*x)), x)


def replacement5090(a, b, c, d, e, n, p, x):
    return Dist(b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5091(a, b, c, d, e, n, p, x):
    return -Dist(b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5092(a, b, c, d, e, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*x)**n/(sin(x)*cos(x)), x), x, asin(c*x)), x)


def replacement5093(a, b, c, d, e, n, x):
    return -Dist(S(1)/d, Subst(Int((a + b*x)**n/(sin(x)*cos(x)), x), x, acos(c*x)), x)


def replacement5094(a, b, c, d, e, f, m, n, p, x):
    return -Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement5095(a, b, c, d, e, f, m, n, p, x):
    return Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement5096(a, b, c, d, e, p, x):
    return Dist(d, Int((a + b*asin(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x), x) - Dist(b*c*d**p/(S(2)*p), Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*asin(c*x))*(d + e*x**S(2))**p/(S(2)*p), x)


def replacement5097(a, b, c, d, e, p, x):
    return Dist(d, Int((a + b*acos(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x), x) + Dist(b*c*d**p/(S(2)*p), Int((-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*acos(c*x))*(d + e*x**S(2))**p/(S(2)*p), x)


def replacement5098(a, b, c, d, e, f, m, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asin(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement5099(a, b, c, d, e, f, m, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acos(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b*c*d**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def With5100(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5101(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def With5102(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
    return Dist(d**p*(a + b*asin(c*x)), u, x) - Dist(b*c*d**p, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x)


def With5103(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
    return Dist(d**p*(a + b*acos(c*x)), u, x) + Dist(b*c*d**p, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x)


def With5104(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
    return -Dist(b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(-c**S(2)*x**S(2) + S(1)), Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), Int(x**m*(d + e*x**S(2))**p, x), x)


def With5105(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(-c**S(2)*x**S(2) + S(1))**p, x)
    return Dist(b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(-c**S(2)*x**S(2) + S(1)), Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), Int(x**m*(d + e*x**S(2))**p, x), x)


def replacement5106(a, b, c, d, e, f, m, n, x):
    return Dist(c**S(2)*sqrt(d + e*x**S(2))/(f**S(2)*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))), x)


def replacement5107(a, b, c, d, e, f, m, n, x):
    return Dist(c**S(2)*sqrt(d + e*x**S(2))/(f**S(2)*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(1))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))), x)


def replacement5108(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement5109(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement5110(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(d + e*x**S(2))/((m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**m*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))), x)


def replacement5111(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(d + e*x**S(2))/((m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**m*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(2))*sqrt(-c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))), x)


def replacement5112(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(2)*d*p/(m + S(2)*p + S(1)), Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(2)*p + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))), x)


def replacement5113(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(2)*d*p/(m + S(2)*p + S(1)), Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(2)*p + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))), x)


def replacement5114(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement5115(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement5116(a, b, c, d, e, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5117(a, b, c, d, e, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5118(a, b, c, d, e, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))), x)


def replacement5119(a, b, c, d, e, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))), x)


def replacement5120(a, b, c, d, e, f, m, n, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*m), Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), x), x) + Dist(b*f*n*sqrt(-c**S(2)*x**S(2) + S(1))/(c*m*sqrt(d + e*x**S(2))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1)), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*sqrt(d + e*x**S(2))/(e*m), x)


def replacement5121(a, b, c, d, e, f, m, n, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*m), Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*f*n*sqrt(-c**S(2)*x**S(2) + S(1))/(c*m*sqrt(d + e*x**S(2))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1)), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*sqrt(d + e*x**S(2))/(e*m), x)


def replacement5122(a, b, c, d, e, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*sin(x)**m, x), x, asin(c*x)), x)


def replacement5123(a, b, c, d, e, m, n, x):
    return -Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*cos(x)**m, x), x, acos(c*x)), x)


def replacement5124(a, b, c, d, e, f, m, x):
    return Simp((f*x)**(m + S(1))*(a + b*asin(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))), x) - Simp(b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))), x)


def replacement5125(a, b, c, d, e, f, m, x):
    return Simp((f*x)**(m + S(1))*(a + b*acos(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))), x) + Simp(b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))), x)


def replacement5126(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((f*x)**m*(a + b*asin(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def replacement5127(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((f*x)**m*(a + b*acos(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def replacement5128(a, b, c, d, e, f, m, n, p, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x), x) + Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asin(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement5129(a, b, c, d, e, f, m, n, p, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x), x) - Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(-1))*(-c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acos(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement5130(a, b, c, d, e, f, m, n, p, x):
    return -Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5131(a, b, c, d, e, f, m, n, p, x):
    return Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) - Simp((f*x)**m*(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5132(a, b, c, d, e, f, m, n, x):
    return -Dist(f*m/(b*c*sqrt(d)*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1)), x), x) + Simp((f*x)**m*(a + b*asin(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5133(a, b, c, d, e, f, m, n, x):
    return Dist(f*m/(b*c*sqrt(d)*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1)), x), x) - Simp((f*x)**m*(a + b*acos(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5134(a, b, c, d, e, f, m, n, p, x):
    return -Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))/(b*f*(n + S(1))), Int((f*x)**(m + S(1))*(a + b*asin(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5135(a, b, c, d, e, f, m, n, p, x):
    return Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) - Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))/(b*f*(n + S(1))), Int((f*x)**(m + S(1))*(a + b*acos(c*x))**(n + S(1))*(-c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) - Simp((f*x)**m*(a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(-c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement5136(a, b, c, d, e, m, n, p, x):
    return Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sin(x)**m*cos(x)**(S(2)*p + S(1)), x), x, asin(c*x)), x)


def replacement5137(a, b, c, d, e, m, n, p, x):
    return -Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sin(x)**(S(2)*p + S(1))*cos(x)**m, x), x, acos(c*x)), x)


def replacement5138(a, b, c, d, e, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(x**m*(a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5139(a, b, c, d, e, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(x**m*(a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5140(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x)


def replacement5141(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x)


def replacement5142(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asin(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5143(a, b, c, d, e, p, x):
    return Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acos(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With5144(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5145(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5146(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x)


def replacement5147(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x)


def replacement5148(a, b, c, d, e, f, m, n, p, x):
    return Int((f*x)**m*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5149(a, b, c, d, e, f, m, n, p, x):
    return Int((f*x)**m*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5150(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((h*x)**m*(a + b*asin(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement5151(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((h*x)**m*(a + b*acos(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement5152(a, b, c, d, e, n, x):
    return Subst(Int((a + b*x)**n*cos(x)/(c*d + e*sin(x)), x), x, asin(c*x))


def replacement5153(a, b, c, d, e, n, x):
    return -Subst(Int((a + b*x)**n*sin(x)/(c*d + e*cos(x)), x), x, acos(c*x))


def replacement5154(a, b, c, d, e, m, n, x):
    return -Dist(b*c*n/(e*(m + S(1))), Int((a + b*asin(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asin(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement5155(a, b, c, d, e, m, n, x):
    return Dist(b*c*n/(e*(m + S(1))), Int((a + b*acos(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acos(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement5156(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x)**m, x), x)


def replacement5157(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x)**m, x), x)


def replacement5158(a, b, c, d, e, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(c*d + e*sin(x))**m*cos(x), x), x, asin(c*x)), x)


def replacement5159(a, b, c, d, e, m, n, x):
    return -Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(c*d + e*cos(x))**m*sin(x), x), x, acos(c*x)), x)


def With5160(Px, a, b, c, x):
    u = IntHide(Px, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5161(Px, a, b, c, x):
    u = IntHide(Px, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5162(Px, a, b, c, n, x):
    return Int(ExpandIntegrand(Px*(a + b*asin(c*x))**n, x), x)


def replacement5163(Px, a, b, c, n, x):
    return Int(ExpandIntegrand(Px*(a + b*acos(c*x))**n, x), x)


def With5164(Px, a, b, c, d, e, m, x):
    u = IntHide(Px*(d + e*x)**m, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5165(Px, a, b, c, d, e, m, x):
    u = IntHide(Px*(d + e*x)**m, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), u, x)


def With5166(a, b, c, d, e, f, g, m, n, p, x):
    u = IntHide((d + e*x)**m*(f + g*x)**p, x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*asin(c*x))**n, u, x)


def With5167(a, b, c, d, e, f, g, m, n, p, x):
    u = IntHide((d + e*x)**m*(f + g*x)**p, x)
    return Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*acos(c*x))**n, u, x)


def With5168(a, b, c, d, e, f, g, h, n, p, x):
    u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*asin(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*asin(c*x))**n, u, x)


def With5169(a, b, c, d, e, f, g, h, n, p, x):
    u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
    return Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*acos(c*x))**(n + S(-1))/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*acos(c*x))**n, u, x)


def replacement5170(Px, a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand(Px*(a + b*asin(c*x))**n*(d + e*x)**m, x), x)


def replacement5171(Px, a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand(Px*(a + b*acos(c*x))**n*(d + e*x)**m, x), x)


def With5172(a, b, c, d, e, f, g, m, p, x):
    u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5173(a, b, c, d, e, f, g, m, p, x):
    u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
    return Dist(b*c, Int(Dist(S(1)/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5174(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x)


def replacement5175(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x)


def replacement5176(a, b, c, d, e, f, g, m, n, x):
    return -Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5177(a, b, c, d, e, f, g, m, n, x):
    return Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5178(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x)


def replacement5179(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x)


def replacement5180(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int(ExpandIntegrand((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5181(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int(ExpandIntegrand((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5182(a, b, c, d, e, f, g, m, n, x):
    return -Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asin(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5183(a, b, c, d, e, f, g, m, n, x):
    return Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*acos(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5184(a, b, c, d, e, f, g, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*(c*f + g*sin(x))**m, x), x, asin(c*x)), x)


def replacement5185(a, b, c, d, e, f, g, m, n, x):
    return -Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*(c*f + g*cos(x))**m, x), x, acos(c*x)), x)


def replacement5186(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x)


def replacement5187(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x)


def replacement5188(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*asin(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5189(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*acos(c*x))**n*(f + g*x)**m*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5190(a, b, c, d, e, f, g, h, m, n, x):
    return -Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asin(c*x))**(n + S(1))/(f + g*x), x), x) + Simp((a + b*asin(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5191(a, b, c, d, e, f, g, h, m, n, x):
    return Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*acos(c*x))**(n + S(1))/(f + g*x), x), x) - Simp((a + b*acos(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))), x)


def replacement5192(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*asin(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x), x)


def replacement5193(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*acos(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x), x)


def With5194(a, b, c, d, e, f, g, m, x):
    u = IntHide((d + e*x)**m*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*asin(c*x), u, x)


def With5195(a, b, c, d, e, f, g, m, x):
    u = IntHide((d + e*x)**m*(f + g*x)**m, x)
    return Dist(b*c, Int(Dist(S(1)/sqrt(-c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*acos(c*x), u, x)


def replacement5196(a, b, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((a + b*asin(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x)


def replacement5197(a, b, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((a + b*acos(c*x))**n*(d + e*x)**m*(f + g*x)**m, x), x)


def With5198(a, b, c, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = IntHide(u, x)
    if InverseFunctionFreeQ(v, x):
        return True
    return False


def replacement5198(a, b, c, u, x):

    v = IntHide(u, x)
    return -Dist(b*c, Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asin(c*x), v, x)


def With5199(a, b, c, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = IntHide(u, x)
    if InverseFunctionFreeQ(v, x):
        return True
    return False


def replacement5199(a, b, c, u, x):

    v = IntHide(u, x)
    return Dist(b*c, Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acos(c*x), v, x)


def With5200(Px, a, b, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)
    if SumQ(u):
        return True
    return False


def replacement5200(Px, a, b, c, d, e, n, p, x):

    u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(d + e*x**S(2))**p, x)
    return Int(u, x)


def With5201(Px, a, b, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)
    if SumQ(u):
        return True
    return False


def replacement5201(Px, a, b, c, d, e, n, p, x):

    u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(d + e*x**S(2))**p, x)
    return Int(u, x)


def With5202(Px, a, b, c, d, e, f, g, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    if SumQ(u):
        return True
    return False


def replacement5202(Px, a, b, c, d, e, f, g, m, n, p, x):

    u = ExpandIntegrand(Px*(a + b*asin(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    return Int(u, x)


def With5203(Px, a, b, c, d, e, f, g, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    if SumQ(u):
        return True
    return False


def replacement5203(Px, a, b, c, d, e, f, g, m, n, p, x):

    u = ExpandIntegrand(Px*(a + b*acos(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    return Int(u, x)


def With5204(RFx, c, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(asin(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement5204(RFx, c, n, x):

    u = ExpandIntegrand(asin(c*x)**n, RFx, x)
    return Int(u, x)


def With5205(RFx, c, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(acos(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement5205(RFx, c, n, x):

    u = ExpandIntegrand(acos(c*x)**n, RFx, x)
    return Int(u, x)


def replacement5206(RFx, a, b, c, n, x):
    return Int(ExpandIntegrand(RFx*(a + b*asin(c*x))**n, x), x)


def replacement5207(RFx, a, b, c, n, x):
    return Int(ExpandIntegrand(RFx*(a + b*acos(c*x))**n, x), x)


def With5208(RFx, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((d + e*x**S(2))**p*asin(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement5208(RFx, c, d, e, n, p, x):

    u = ExpandIntegrand((d + e*x**S(2))**p*asin(c*x)**n, RFx, x)
    return Int(u, x)


def With5209(RFx, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((d + e*x**S(2))**p*acos(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement5209(RFx, c, d, e, n, p, x):

    u = ExpandIntegrand((d + e*x**S(2))**p*acos(c*x)**n, RFx, x)
    return Int(u, x)


def replacement5210(RFx, a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*asin(c*x))**n, x), x)


def replacement5211(RFx, a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*acos(c*x))**n, x), x)


def replacement5212(a, b, c, n, u, x):
    return Int(u*(a + b*asin(c*x))**n, x)


def replacement5213(a, b, c, n, u, x):
    return Int(u*(a + b*acos(c*x))**n, x)


def replacement5214(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*asin(x))**n, x), x, c + d*x), x)


def replacement5215(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acos(x))**n, x), x, c + d*x), x)


def replacement5216(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*asin(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5217(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acos(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5218(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*asin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x), x)


def replacement5219(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x), x)


def replacement5220(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*asin(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5221(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acos(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5222(a, b, c, d, x):
    return Simp(x*sqrt(a + b*asin(c + d*x**S(2))), x) + Simp(sqrt(Pi)*x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(sqrt(c/b)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x) - Simp(sqrt(Pi)*x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(sqrt(c/b)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x)


def replacement5223(a, b, d, x):
    return Simp(-S(2)*sqrt(a + b*acos(d*x**S(2) + S(1)))*sin(acos(d*x**S(2) + S(1))/S(2))**S(2)/(d*x), x) - Simp(S(2)*sqrt(Pi)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x*sqrt(S(1)/b)), x) + Simp(S(2)*sqrt(Pi)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x*sqrt(S(1)/b)), x)


def replacement5224(a, b, d, x):
    return Simp(S(2)*sqrt(a + b*acos(d*x**S(2) + S(-1)))*cos(acos(d*x**S(2) + S(-1))/S(2))**S(2)/(d*x), x) - Simp(S(2)*sqrt(Pi)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x*sqrt(S(1)/b)), x) - Simp(S(2)*sqrt(Pi)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x*sqrt(S(1)/b)), x)


def replacement5225(a, b, c, d, n, x):
    return -Dist(S(4)*b**S(2)*n*(n + S(-1)), Int((a + b*asin(c + d*x**S(2)))**(n + S(-2)), x), x) + Simp(x*(a + b*asin(c + d*x**S(2)))**n, x) + Simp(S(2)*b*n*(a + b*asin(c + d*x**S(2)))**(n + S(-1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x), x)


def replacement5226(a, b, c, d, n, x):
    return -Dist(S(4)*b**S(2)*n*(n + S(-1)), Int((a + b*acos(c + d*x**S(2)))**(n + S(-2)), x), x) + Simp(x*(a + b*acos(c + d*x**S(2)))**n, x) - Simp(S(2)*b*n*(a + b*acos(c + d*x**S(2)))**(n + S(-1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(d*x), x)


def replacement5227(a, b, c, d, x):
    return -Simp(x*(c*cos(a/(S(2)*b)) - sin(a/(S(2)*b)))*CosIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x) - Simp(x*(c*cos(a/(S(2)*b)) + sin(a/(S(2)*b)))*SinIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x)


def replacement5228(a, b, d, x):
    return Simp(sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(2)*b*sqrt(-d*x**S(2))), x) + Simp(sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(2)*b*sqrt(-d*x**S(2))), x)


def replacement5229(a, b, d, x):
    return Simp(sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x) - Simp(sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x)


def replacement5230(a, b, c, d, x):
    return -Simp(sqrt(Pi)*x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(a + b*asin(c + d*x**S(2)))/(sqrt(Pi)*sqrt(b*c)))/(sqrt(b*c)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x) - Simp(sqrt(Pi)*x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(a + b*asin(c + d*x**S(2)))/(sqrt(Pi)*sqrt(b*c)))/(sqrt(b*c)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x)


def replacement5231(a, b, d, x):
    return Simp(-S(2)*sqrt(Pi/b)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x), x) - Simp(S(2)*sqrt(Pi/b)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x), x)


def replacement5232(a, b, d, x):
    return Simp(S(2)*sqrt(Pi/b)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x), x) - Simp(S(2)*sqrt(Pi/b)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x), x)


def replacement5233(a, b, c, d, x):
    return -Simp(sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(b*d*x*sqrt(a + b*asin(c + d*x**S(2)))), x) + Simp(sqrt(Pi)*x*(c/b)**(S(3)/2)*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelS(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2))), x) - Simp(sqrt(Pi)*x*(c/b)**(S(3)/2)*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*FresnelC(sqrt(c/(Pi*b))*sqrt(a + b*asin(c + d*x**S(2))))/(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2))), x)


def replacement5234(a, b, d, x):
    return Simp(sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(b*d*x*sqrt(a + b*acos(d*x**S(2) + S(1)))), x) - Simp(S(2)*sqrt(Pi)*(S(1)/b)**(S(3)/2)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(a/(S(2)*b))*sin(acos(d*x**S(2) + S(1))/S(2))/(d*x), x) + Simp(S(2)*sqrt(Pi)*(S(1)/b)**(S(3)/2)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(1))))*sin(acos(d*x**S(2) + S(1))/S(2))*cos(a/(S(2)*b))/(d*x), x)


def replacement5235(a, b, d, x):
    return Simp(sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(b*d*x*sqrt(a + b*acos(d*x**S(2) + S(-1)))), x) - Simp(S(2)*sqrt(Pi)*(S(1)/b)**(S(3)/2)*FresnelC(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*cos(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x), x) - Simp(S(2)*sqrt(Pi)*(S(1)/b)**(S(3)/2)*FresnelS(sqrt(S(1)/(Pi*b))*sqrt(a + b*acos(d*x**S(2) + S(-1))))*sin(a/(S(2)*b))*cos(acos(d*x**S(2) + S(-1))/S(2))/(d*x), x)


def replacement5236(a, b, c, d, x):
    return Simp(x*(-c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*SinIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x) - Simp(x*(c*sin(a/(S(2)*b)) + cos(a/(S(2)*b)))*CosIntegral(c*(a + b*asin(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(-c*sin(asin(c + d*x**S(2))/S(2)) + cos(asin(c + d*x**S(2))/S(2)))), x) - Simp(sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(a + b*asin(c + d*x**S(2)))), x)


def replacement5237(a, b, d, x):
    return Simp(sqrt(-d**S(2)*x**S(4) - S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*acos(d*x**S(2) + S(1)))), x) + Simp(sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(-d*x**S(2))), x) - Simp(sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(-d*x**S(2))), x)


def replacement5238(a, b, d, x):
    return Simp(sqrt(-d**S(2)*x**S(4) + S(2)*d*x**S(2))/(S(2)*b*d*x*(a + b*acos(d*x**S(2) + S(-1)))), x) - Simp(sqrt(S(2))*x*CosIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*cos(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x) - Simp(sqrt(S(2))*x*SinIntegral((a + b*acos(d*x**S(2) + S(-1)))/(S(2)*b))*sin(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x)


def replacement5239(a, b, c, d, n, x):
    return -Dist(S(1)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), Int((a + b*asin(c + d*x**S(2)))**(n + S(2)), x), x) + Simp(x*(a + b*asin(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), x) + Simp((a + b*asin(c + d*x**S(2)))**(n + S(1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))), x)


def replacement5240(a, b, c, d, n, x):
    return -Dist(S(1)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), Int((a + b*acos(c + d*x**S(2)))**(n + S(2)), x), x) + Simp(x*(a + b*acos(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), x) - Simp((a + b*acos(c + d*x**S(2)))**(n + S(1))*sqrt(-S(2)*c*d*x**S(2) - d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))), x)


def replacement5241(a, n, p, x):
    return Dist(S(1)/p, Subst(Int(x**n/tan(x), x), x, asin(a*x**p)), x)


def replacement5242(a, n, p, x):
    return -Dist(S(1)/p, Subst(Int(x**n*tan(x), x), x, acos(a*x**p)), x)


def replacement5243(a, b, c, m, n, u, x):
    return Int(u*acsc(a/c + b*x**n/c)**m, x)


def replacement5244(a, b, c, m, n, u, x):
    return Int(u*asec(a/c + b*x**n/c)**m, x)


def replacement5245(b, n, x):
    return Dist(sqrt(-b*x**S(2))/(b*x), Subst(Int(asin(x)**n/sqrt(S(1) - x**S(2)), x), x, sqrt(b*x**S(2) + S(1))), x)


def replacement5246(b, n, x):
    return Dist(sqrt(-b*x**S(2))/(b*x), Subst(Int(acos(x)**n/sqrt(S(1) - x**S(2)), x), x, sqrt(b*x**S(2) + S(1))), x)


def replacement5247(a, b, c, f, n, u, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + sin(x)/b))*cos(x), x), x, asin(a + b*x)), x)


def replacement5248(a, b, c, f, n, u, x):
    return -Dist(S(1)/b, Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + cos(x)/b))*sin(x), x), x, acos(a + b*x)), x)


def replacement5249(a, b, c, d, x):
    return -Dist(x*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)/sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), Int(x*(S(2)*a*sqrt(c + d*x**S(2)) + b*d)/(sqrt(c + d*x**S(2))*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), x), x) + Simp(x*asin(a*x**S(2) + b*sqrt(c + d*x**S(2))), x)


def replacement5250(a, b, c, d, x):
    return Dist(x*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)/sqrt(-x**S(2)*(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), Int(x*(S(2)*a*sqrt(c + d*x**S(2)) + b*d)/(sqrt(c + d*x**S(2))*sqrt(a**S(2)*x**S(2) + S(2)*a*b*sqrt(c + d*x**S(2)) + b**S(2)*d)), x), x) + Simp(x*acos(a*x**S(2) + b*sqrt(c + d*x**S(2))), x)


def replacement5251(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/sqrt(S(1) - u**S(2)), x), x) + Simp(x*asin(u), x)


def replacement5252(u, x):
    return Int(SimplifyIntegrand(x*D(u, x)/sqrt(S(1) - u**S(2)), x), x) + Simp(x*acos(u), x)


def replacement5253(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(S(1) - u**S(2)), x), x), x) + Simp((a + b*asin(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement5254(a, b, c, d, m, u, x):
    return Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(S(1) - u**S(2)), x), x), x) + Simp((a + b*acos(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With5255(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5255(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/sqrt(S(1) - u**S(2)), x), x), x) + Dist(a + b*asin(u), w, x)


def With5256(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5256(a, b, u, v, x):

    w = IntHide(v, x)
    return Dist(b, Int(SimplifyIntegrand(w*D(u, x)/sqrt(S(1) - u**S(2)), x), x), x) + Dist(a + b*acos(u), w, x)


def replacement5257(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*ArcTan(c*x))**n, x)


def replacement5258(a, b, c, n, x):
    return Dist(b*c*n, Int(x*(a + b*acot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*acot(c*x))**n, x)


def replacement5259(a, b, c, n, x):
    return Int((a + b*ArcTan(c*x))**n, x)


def replacement5260(a, b, c, n, x):
    return Int((a + b*acot(c*x))**n, x)


def replacement5261(a, b, c, d, e, n, x):
    return Dist(b*c*n/e, Int((a + b*ArcTan(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x), x) - Simp((a + b*ArcTan(c*x))**n*log(S(2)*d/(d + e*x))/e, x)


def replacement5262(a, b, c, d, e, n, x):
    return -Dist(b*c*n/e, Int((a + b*acot(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x), x) - Simp((a + b*acot(c*x))**n*log(S(2)*d/(d + e*x))/e, x)


def replacement5263(c, d, e, x):
    return Simp(I*PolyLog(S(2), Simp(I*c*(d + e*x)/(I*c*d - e), x))/(S(2)*e), x) - Simp(I*PolyLog(S(2), Simp(I*c*(d + e*x)/(I*c*d + e), x))/(S(2)*e), x) - Simp(ArcTan(c*d/e)*log(d + e*x)/e, x)


def replacement5264(c, d, e, x):
    return Dist(I/S(2), Int(log(-I*c*x + S(1))/(d + e*x), x), x) - Dist(I/S(2), Int(log(I*c*x + S(1))/(d + e*x), x), x)


def replacement5265(c, d, e, x):
    return Dist(I/S(2), Int(log(S(1) - I/(c*x))/(d + e*x), x), x) - Dist(I/S(2), Int(log(S(1) + I/(c*x))/(d + e*x), x), x)


def replacement5266(a, b, c, d, e, x):
    return Dist(b, Int(ArcTan(c*x)/(d + e*x), x), x) + Simp(a*log(RemoveContent(d + e*x, x))/e, x)


def replacement5267(a, b, c, d, e, x):
    return Dist(b, Int(acot(c*x)/(d + e*x), x), x) + Simp(a*log(RemoveContent(d + e*x, x))/e, x)


def replacement5268(a, b, c, d, e, p, x):
    return -Dist(b*c/(e*(p + S(1))), Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*ArcTan(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))), x)


def replacement5269(a, b, c, d, e, p, x):
    return Dist(b*c/(e*(p + S(1))), Int((d + e*x)**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acot(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))), x)


def replacement5270(a, b, c, n, x):
    return -Dist(S(2)*b*c*n, Int((a + b*ArcTan(c*x))**(n + S(-1))*atanh(S(1) - S(2)*I/(-c*x + I))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(S(2)*(a + b*ArcTan(c*x))**n*atanh(S(1) - S(2)*I/(-c*x + I)), x)


def replacement5271(a, b, c, n, x):
    return Dist(S(2)*b*c*n, Int((a + b*acot(c*x))**(n + S(-1))*acoth(S(1) - S(2)*I/(-c*x + I))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(S(2)*(a + b*acot(c*x))**n*acoth(S(1) - S(2)*I/(-c*x + I)), x)


def replacement5272(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*ArcTan(c*x))**n/(m + S(1)), x)


def replacement5273(a, b, c, m, n, x):
    return Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*acot(c*x))**n/(m + S(1)), x)


def replacement5274(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x)


def replacement5275(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acot(c*x))**n*(d + e*x)**p, x), x)


def replacement5276(a, b, c, d, e, n, p, x):
    return Int((a + b*ArcTan(c*x))**n*(d + e*x)**p, x)


def replacement5277(a, b, c, d, e, n, p, x):
    return Int((a + b*acot(c*x))**n*(d + e*x)**p, x)


def replacement5278(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**n/(d + e*x), x), x)


def replacement5279(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-1))*(a + b*acot(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-1))*(a + b*acot(c*x))**n/(d + e*x), x), x)


def replacement5280(a, b, c, d, e, n, x):
    return -Dist(b*c*n/d, Int((a + b*ArcTan(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*ArcTan(c*x))**n*log(S(2)*e*x/(d + e*x))/d, x)


def replacement5281(a, b, c, d, e, n, x):
    return Dist(b*c*n/d, Int((a + b*acot(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acot(c*x))**n*log(S(2)*e*x/(d + e*x))/d, x)


def replacement5282(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*ArcTan(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(1))*(a + b*ArcTan(c*x))**n/(d + e*x), x), x)


def replacement5283(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acot(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(1))*(a + b*acot(c*x))**n/(d + e*x), x), x)


def replacement5284(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x), x)


def replacement5285(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*acot(c*x))**n*(d + e*x)**p, x), x)


def replacement5286(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x)**p, x)


def replacement5287(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acot(c*x))**n*(d + e*x)**p, x)


def replacement5288(a, b, c, d, e, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) - Simp(b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement5289(a, b, c, d, e, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*acot(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement5290(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b**S(2)*d*n*(n + S(-1))/(S(2)*p*(S(2)*p + S(1))), Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) - Simp(b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement5291(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(b**S(2)*d*n*(n + S(-1))/(S(2)*p*(S(2)*p + S(1))), Int((a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*acot(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*n*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement5292(a, b, c, d, e, x):
    return Simp(log(RemoveContent(a + b*ArcTan(c*x), x))/(b*c*d), x)


def replacement5293(a, b, c, d, e, x):
    return -Simp(log(RemoveContent(a + b*acot(c*x), x))/(b*c*d), x)


def replacement5294(a, b, c, d, e, n, x):
    return Simp((a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5295(a, b, c, d, e, n, x):
    return -Simp((a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5296(a, b, c, d, e, x):
    return Simp(I*b*PolyLog(S(2), -I*sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x) - Simp(I*b*PolyLog(S(2), I*sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x) + Simp(-S(2)*I*(a + b*ArcTan(c*x))*ArcTan(sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x)


def replacement5297(a, b, c, d, e, x):
    return -Simp(I*b*PolyLog(S(2), -I*sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x) + Simp(I*b*PolyLog(S(2), I*sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x) + Simp(-S(2)*I*(a + b*acot(c*x))*ArcTan(sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/(c*sqrt(d)), x)


def replacement5298(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*sqrt(d)), Subst(Int((a + b*x)**n/cos(x), x), x, ArcTan(c*x)), x)


def replacement5299(a, b, c, d, e, n, x):
    return -Dist(x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n/sin(x), x), x, acot(c*x)), x)


def replacement5300(a, b, c, d, e, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*ArcTan(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x)


def replacement5301(a, b, c, d, e, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*acot(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x)


def replacement5302(a, b, c, d, e, n, x):
    return -Dist(b*c*n/S(2), Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) + Simp(x*(a + b*ArcTan(c*x))**n/(S(2)*d*(d + e*x**S(2))), x) + Simp((a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))), x)


def replacement5303(a, b, c, d, e, n, x):
    return Dist(b*c*n/S(2), Int(x*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) + Simp(x*(a + b*acot(c*x))**n/(S(2)*d*(d + e*x**S(2))), x) - Simp((a + b*acot(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))), x)


def replacement5304(a, b, c, d, e, x):
    return Simp(b/(c*d*sqrt(d + e*x**S(2))), x) + Simp(x*(a + b*ArcTan(c*x))/(d*sqrt(d + e*x**S(2))), x)


def replacement5305(a, b, c, d, e, x):
    return -Simp(b/(c*d*sqrt(d + e*x**S(2))), x) + Simp(x*(a + b*acot(c*x))/(d*sqrt(d + e*x**S(2))), x)


def replacement5306(a, b, c, d, e, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) + Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement5307(a, b, c, d, e, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement5308(a, b, c, d, e, n, x):
    return -Dist(b**S(2)*n*(n + S(-1)), Int((a + b*ArcTan(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x), x) + Simp(x*(a + b*ArcTan(c*x))**n/(d*sqrt(d + e*x**S(2))), x) + Simp(b*n*(a + b*ArcTan(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))), x)


def replacement5309(a, b, c, d, e, n, x):
    return -Dist(b**S(2)*n*(n + S(-1)), Int((a + b*acot(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x), x) + Simp(x*(a + b*acot(c*x))**n/(d*sqrt(d + e*x**S(2))), x) - Simp(b*n*(a + b*acot(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))), x)


def replacement5310(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b**S(2)*n*(n + S(-1))/(S(4)*(p + S(1))**S(2)), Int((a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Simp(x*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x) + Simp(b*n*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x)


def replacement5311(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b**S(2)*n*(n + S(-1))/(S(4)*(p + S(1))**S(2)), Int((a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Simp(x*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x) - Simp(b*n*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x)


def replacement5312(a, b, c, d, e, n, p, x):
    return -Dist(S(2)*c*(p + S(1))/(b*(n + S(1))), Int(x*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5313(a, b, c, d, e, n, p, x):
    return Dist(S(2)*c*(p + S(1))/(b*(n + S(1))), Int(x*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) - Simp((a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5314(a, b, c, d, e, n, p, x):
    return Dist(d**p/c, Subst(Int((a + b*x)**n*cos(x)**(-S(2)*p + S(-2)), x), x, ArcTan(c*x)), x)


def replacement5315(a, b, c, d, e, n, p, x):
    return Dist(d**(p + S(1)/2)*sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5316(a, b, c, d, e, n, p, x):
    return -Dist(d**p/c, Subst(Int((a + b*x)**n*sin(x)**(-S(2)*p + S(-2)), x), x, acot(c*x)), x)


def replacement5317(a, b, c, d, e, n, p, x):
    return -Dist(d**(p + S(1)/2)*x*sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n*sin(x)**(-S(2)*p + S(-2)), x), x, acot(c*x)), x)


def replacement5318(c, d, e, x):
    return Dist(I/S(2), Int(log(-I*c*x + S(1))/(d + e*x**S(2)), x), x) - Dist(I/S(2), Int(log(I*c*x + S(1))/(d + e*x**S(2)), x), x)


def replacement5319(c, d, e, x):
    return Dist(I/S(2), Int(log(S(1) - I/(c*x))/(d + e*x**S(2)), x), x) - Dist(I/S(2), Int(log(S(1) + I/(c*x))/(d + e*x**S(2)), x), x)


def replacement5320(a, b, c, d, e, x):
    return Dist(a, Int(S(1)/(d + e*x**S(2)), x), x) + Dist(b, Int(ArcTan(c*x)/(d + e*x**S(2)), x), x)


def replacement5321(a, b, c, d, e, x):
    return Dist(a, Int(S(1)/(d + e*x**S(2)), x), x) + Dist(b, Int(acot(c*x)/(d + e*x**S(2)), x), x)


def With5322(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*ArcTan(c*x), u, x)


def With5323(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return Dist(b*c, Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acot(c*x), u, x)


def replacement5324(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5325(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5326(a, b, c, d, e, n, p, x):
    return Int((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5327(a, b, c, d, e, n, p, x):
    return Int((a + b*acot(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5328(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5329(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5330(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*ArcTan(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5331(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acot(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5332(a, b, c, d, e, n, x):
    return -Dist(S(1)/(c*d), Int((a + b*ArcTan(c*x))**n/(-c*x + I), x), x) - Simp(I*(a + b*ArcTan(c*x))**(n + S(1))/(b*e*(n + S(1))), x)


def replacement5333(a, b, c, d, e, n, x):
    return -Dist(S(1)/(c*d), Int((a + b*acot(c*x))**n/(-c*x + I), x), x) + Simp(I*(a + b*acot(c*x))**(n + S(1))/(b*e*(n + S(1))), x)


def replacement5334(a, b, c, d, e, n, x):
    return -Dist(S(1)/(b*c*d*(n + S(1))), Int((a + b*ArcTan(c*x))**(n + S(1)), x), x) + Simp(x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5335(a, b, c, d, e, n, x):
    return Dist(S(1)/(b*c*d*(n + S(1))), Int((a + b*acot(c*x))**(n + S(1)), x), x) - Simp(x*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5336(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5337(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5338(a, b, c, d, e, n, x):
    return Dist(I/d, Int((a + b*ArcTan(c*x))**n/(x*(c*x + I)), x), x) - Simp(I*(a + b*ArcTan(c*x))**(n + S(1))/(b*d*(n + S(1))), x)


def replacement5339(a, b, c, d, e, n, x):
    return Dist(I/d, Int((a + b*acot(c*x))**n/(x*(c*x + I)), x), x) + Simp(I*(a + b*acot(c*x))**(n + S(1))/(b*d*(n + S(1))), x)


def replacement5340(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*ArcTan(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5341(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acot(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acot(c*x))**n/(d + e*x**S(2)), x), x)


def replacement5342(a, b, c, d, e, m, n, x):
    return -Dist(m/(b*c*d*(n + S(1))), Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1)), x), x) + Simp(x**m*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5343(a, b, c, d, e, m, n, x):
    return Dist(m/(b*c*d*(n + S(1))), Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1)), x), x) - Simp(x**m*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement5344(a, b, c, d, e, m, x):
    return Int(ExpandIntegrand(a + b*ArcTan(c*x), x**m/(d + e*x**S(2)), x), x)


def replacement5345(a, b, c, d, e, m, x):
    return Int(ExpandIntegrand(a + b*acot(c*x), x**m/(d + e*x**S(2)), x), x)


def replacement5346(a, b, c, d, e, n, p, x):
    return -Dist(b*n/(S(2)*c*(p + S(1))), Int((a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5347(a, b, c, d, e, n, p, x):
    return Dist(b*n/(S(2)*c*(p + S(1))), Int((a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5348(a, b, c, d, e, n, x):
    return -Dist(S(4)/(b**S(2)*(n + S(1))*(n + S(2))), Int(x*(a + b*ArcTan(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x), x) - Simp((a + b*ArcTan(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))), x) + Simp(x*(a + b*ArcTan(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))), x)


def replacement5349(a, b, c, d, e, n, x):
    return -Dist(S(4)/(b**S(2)*(n + S(1))*(n + S(2))), Int(x*(a + b*acot(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x), x) - Simp((a + b*acot(c*x))**(n + S(2))*(-c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))), x) - Simp(x*(a + b*acot(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))), x)


def replacement5350(a, b, c, d, e, p, x):
    return -Dist(S(1)/(S(2)*c**S(2)*d*(p + S(1))), Int((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)), x) + Simp(x*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))), x)


def replacement5351(a, b, c, d, e, p, x):
    return -Dist(S(1)/(S(2)*c**S(2)*d*(p + S(1))), Int((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) + Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)), x) + Simp(x*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))), x)


def replacement5352(a, b, c, d, e, n, x):
    return Dist(b*n/(S(2)*c), Int(x*(a + b*ArcTan(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) + Simp((a + b*ArcTan(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))), x) - Simp(x*(a + b*ArcTan(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))), x)


def replacement5353(a, b, c, d, e, n, x):
    return -Dist(b*n/(S(2)*c), Int(x*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) - Simp((a + b*acot(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))), x) - Simp(x*(a + b*acot(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))), x)


def replacement5354(a, b, c, d, e, m, p, x):
    return Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) + Simp(b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x) - Simp(x**(m + S(-1))*(a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x)


def replacement5355(a, b, c, d, e, m, p, x):
    return Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x) - Simp(x**(m + S(-1))*(a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x)


def replacement5356(a, b, c, d, e, m, n, p, x):
    return -Dist(b**S(2)*n*(n + S(-1))/m**S(2), Int(x**m*(a + b*ArcTan(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) + Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(x**(m + S(-1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x) + Simp(b*n*x**m*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x)


def replacement5357(a, b, c, d, e, m, n, p, x):
    return -Dist(b**S(2)*n*(n + S(-1))/m**S(2), Int(x**m*(a + b*acot(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) + Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(x**(m + S(-1))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x) - Simp(b*n*x**m*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x)


def replacement5358(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5359(a, b, c, d, e, m, n, p, x):
    return Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) - Simp(x**m*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5360(a, b, c, d, e, m, n, p, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp(x**(m + S(1))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))), x)


def replacement5361(a, b, c, d, e, m, n, p, x):
    return Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp(x**(m + S(1))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))), x)


def replacement5362(a, b, c, d, e, m, x):
    return Dist(d/(m + S(2)), Int(x**m*(a + b*ArcTan(c*x))/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*d/(m + S(2)), Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*ArcTan(c*x))*sqrt(d + e*x**S(2))/(m + S(2)), x)


def replacement5363(a, b, c, d, e, m, x):
    return Dist(d/(m + S(2)), Int(x**m*(a + b*acot(c*x))/sqrt(d + e*x**S(2)), x), x) + Dist(b*c*d/(m + S(2)), Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acot(c*x))*sqrt(d + e*x**S(2))/(m + S(2)), x)


def replacement5364(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5365(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5366(a, b, c, d, e, m, n, p, x):
    return Dist(d, Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(c**S(2)*d, Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x)


def replacement5367(a, b, c, d, e, m, n, p, x):
    return Dist(d, Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) + Dist(c**S(2)*d, Int(x**(m + S(2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x)


def replacement5368(a, b, c, d, e, m, n, x):
    return -Dist((m + S(-1))/(c**S(2)*m), Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*n/(c*m), Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(-1))*(a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m), x)


def replacement5369(a, b, c, d, e, m, n, x):
    return -Dist((m + S(-1))/(c**S(2)*m), Int(x**(m + S(-2))*(a + b*acot(c*x))**n/sqrt(d + e*x**S(2)), x), x) + Dist(b*n/(c*m), Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(-1))*(a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m), x)


def replacement5370(a, b, c, d, e, x):
    return Simp(-S(2)*(a + b*ArcTan(c*x))*atanh(sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x) + Simp(I*b*PolyLog(S(2), -sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x) - Simp(I*b*PolyLog(S(2), sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x)


def replacement5371(a, b, c, d, e, x):
    return Simp(-S(2)*(a + b*acot(c*x))*atanh(sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x) - Simp(I*b*PolyLog(S(2), -sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x) + Simp(I*b*PolyLog(S(2), sqrt(I*c*x + S(1))/sqrt(-I*c*x + S(1)))/sqrt(d), x)


def replacement5372(a, b, c, d, e, n, x):
    return Dist(S(1)/sqrt(d), Subst(Int((a + b*x)**n/sin(x), x), x, ArcTan(c*x)), x)


def replacement5373(a, b, c, d, e, n, x):
    return -Dist(c*x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n/cos(x), x), x, acot(c*x)), x)


def replacement5374(a, b, c, d, e, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*ArcTan(c*x))**n/(x*sqrt(c**S(2)*x**S(2) + S(1))), x), x)


def replacement5375(a, b, c, d, e, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*acot(c*x))**n/(x*sqrt(c**S(2)*x**S(2) + S(1))), x), x)


def replacement5376(a, b, c, d, e, n, x):
    return Dist(b*c*n, Int((a + b*ArcTan(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x), x) - Simp((a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(d*x), x)


def replacement5377(a, b, c, d, e, n, x):
    return -Dist(b*c*n, Int((a + b*acot(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x), x) - Simp((a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(d*x), x)


def replacement5378(a, b, c, d, e, m, n, x):
    return -Dist(c**S(2)*(m + S(2))/(m + S(1)), Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*ArcTan(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))), x)


def replacement5379(a, b, c, d, e, m, n, x):
    return -Dist(c**S(2)*(m + S(2))/(m + S(1)), Int(x**(m + S(2))*(a + b*acot(c*x))**n/sqrt(d + e*x**S(2)), x), x) + Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acot(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))), x)


def replacement5380(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5381(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5382(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/d, Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5383(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement5384(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) - Dist(c*(m + S(2)*p + S(2))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*ArcTan(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5385(a, b, c, d, e, m, n, p, x):
    return Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Dist(c*(m + S(2)*p + S(2))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) - Simp(x**m*(a + b*acot(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement5386(a, b, c, d, e, m, n, p, x):
    return Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sin(x)**m*cos(x)**(-m - S(2)*p + S(-2)), x), x, ArcTan(c*x)), x)


def replacement5387(a, b, c, d, e, m, n, p, x):
    return Dist(d**(p + S(1)/2)*sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int(x**m*(a + b*ArcTan(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement5388(a, b, c, d, e, m, n, p, x):
    return -Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sin(x)**(-m - S(2)*p + S(-2))*cos(x)**m, x), x, acot(c*x)), x)


def replacement5389(a, b, c, d, e, m, n, p, x):
    return -Dist(c**(-m)*d**(p + S(1)/2)*x*sqrt((c**S(2)*x**S(2) + S(1))/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n*sin(x)**(-m - S(2)*p + S(-2))*cos(x)**m, x), x, acot(c*x)), x)


def replacement5390(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*ArcTan(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5391(a, b, c, d, e, p, x):
    return Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acot(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With5392(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*ArcTan(c*x), u, x)


def With5393(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return Dist(b*c, Int(SimplifyIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acot(c*x), u, x)


def replacement5394(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand((a + b*ArcTan(c*x))**n, x**m*(d + e*x**S(2))**p, x), x)


def replacement5395(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acot(c*x))**n, x**m*(d + e*x**S(2))**p, x), x)


def replacement5396(a, b, c, d, e, m, p, x):
    return Dist(a, Int(x**m*(d + e*x**S(2))**p, x), x) + Dist(b, Int(x**m*(d + e*x**S(2))**p*ArcTan(c*x), x), x)


def replacement5397(a, b, c, d, e, m, p, x):
    return Dist(a, Int(x**m*(d + e*x**S(2))**p, x), x) + Dist(b, Int(x**m*(d + e*x**S(2))**p*acot(c*x), x), x)


def replacement5398(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*ArcTan(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5399(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acot(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5400(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*ArcTan(c*x))**n*log(S(1) - u)/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*ArcTan(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x), x)


def replacement5401(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) + S(1)/u, x))/(d + e*x**S(2)), x), x)


def replacement5402(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*ArcTan(c*x))**n*log(S(1) - u)/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*ArcTan(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x), x)


def replacement5403(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*acot(c*x))**n*log(SimplifyIntegrand(S(1) + S(1)/u, x))/(d + e*x**S(2)), x), x)


def replacement5404(a, b, c, d, e, n, u, x):
    return -Dist(I*b*n/S(2), Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) + Simp(I*(a + b*ArcTan(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement5405(a, b, c, d, e, n, u, x):
    return Dist(I*b*n/S(2), Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) + Simp(I*(a + b*acot(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement5406(a, b, c, d, e, n, u, x):
    return Dist(I*b*n/S(2), Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) - Simp(I*(a + b*ArcTan(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement5407(a, b, c, d, e, n, u, x):
    return -Dist(I*b*n/S(2), Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) - Simp(I*(a + b*acot(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement5408(a, b, c, d, e, n, p, u, x):
    return Dist(I*b*n/S(2), Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) - Simp(I*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement5409(a, b, c, d, e, n, p, u, x):
    return -Dist(I*b*n/S(2), Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) - Simp(I*(a + b*acot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement5410(a, b, c, d, e, n, p, u, x):
    return -Dist(I*b*n/S(2), Int((a + b*ArcTan(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) + Simp(I*(a + b*ArcTan(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement5411(a, b, c, d, e, n, p, u, x):
    return Dist(I*b*n/S(2), Int((a + b*acot(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) + Simp(I*(a + b*acot(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement5412(a, b, c, d, e, x):
    return Simp((log(a + b*ArcTan(c*x)) - log(a + b*acot(c*x)))/(b*c*d*(S(2)*a + b*ArcTan(c*x) + b*acot(c*x))), x)


def replacement5413(a, b, c, d, e, m, n, x):
    return Dist(n/(m + S(1)), Int((a + b*ArcTan(c*x))**(n + S(-1))*(a + b*acot(c*x))**(m + S(1))/(d + e*x**S(2)), x), x) - Simp((a + b*ArcTan(c*x))**n*(a + b*acot(c*x))**(m + S(1))/(b*c*d*(m + S(1))), x)


def replacement5414(a, b, c, d, e, m, n, x):
    return Dist(n/(m + S(1)), Int((a + b*ArcTan(c*x))**(m + S(1))*(a + b*acot(c*x))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp((a + b*ArcTan(c*x))**(m + S(1))*(a + b*acot(c*x))**n/(b*c*d*(m + S(1))), x)


def replacement5415(a, c, d, n, x):
    return Dist(I/S(2), Int(log(-I*a*x + S(1))/(c + d*x**n), x), x) - Dist(I/S(2), Int(log(I*a*x + S(1))/(c + d*x**n), x), x)


def replacement5416(a, c, d, n, x):
    return Dist(I/S(2), Int(log(S(1) - I/(a*x))/(c + d*x**n), x), x) - Dist(I/S(2), Int(log(S(1) + I/(a*x))/(c + d*x**n), x), x)


def replacement5417(a, b, c, d, e, f, g, x):
    return -Dist(b*c, Int(x*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g, Int(x**S(2)*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x), x) + Simp(x*(a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2))), x)


def replacement5418(a, b, c, d, e, f, g, x):
    return Dist(b*c, Int(x*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g, Int(x**S(2)*(a + b*acot(c*x))/(f + g*x**S(2)), x), x) + Simp(x*(a + b*acot(c*x))*(d + e*log(f + g*x**S(2))), x)


def replacement5419(a, b, c, d, e, f, g, m, x):
    return -Dist(b*c/(m + S(1)), Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g/(m + S(1)), Int(x**(m + S(2))*(a + b*ArcTan(c*x))/(f + g*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)), x)


def replacement5420(a, b, c, d, e, f, g, m, x):
    return Dist(b*c/(m + S(1)), Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g/(m + S(1)), Int(x**(m + S(2))*(a + b*acot(c*x))/(f + g*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acot(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)), x)


def With5421(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*ArcTan(c*x), u, x)


def With5422(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
    return Dist(b*c, Int(ExpandIntegrand(u/(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acot(c*x), u, x)


def With5423(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(a + b*ArcTan(c*x)), x)
    return -Dist(S(2)*e*g, Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)


def With5424(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(a + b*acot(c*x)), x)
    return -Dist(S(2)*e*g, Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)


def replacement5425(a, b, c, d, e, f, g, x):
    return -Dist(b/c, Int((a + b*ArcTan(c*x))*(d + e*log(f + g*x**S(2))), x), x) + Dist(b*c*e, Int(x**S(2)*(a + b*ArcTan(c*x))/(c**S(2)*x**S(2) + S(1)), x), x) - Simp(e*x**S(2)*(a + b*ArcTan(c*x))**S(2)/S(2), x) + Simp((a + b*ArcTan(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g), x)


def replacement5426(a, b, c, d, e, f, g, x):
    return Dist(b/c, Int((a + b*acot(c*x))*(d + e*log(f + g*x**S(2))), x), x) - Dist(b*c*e, Int(x**S(2)*(a + b*acot(c*x))/(c**S(2)*x**S(2) + S(1)), x), x) - Simp(e*x**S(2)*(a + b*acot(c*x))**S(2)/S(2), x) + Simp((a + b*acot(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g), x)


def replacement5427(a, n, x):
    return Int((-I*a*x + S(1))**(I*n/S(2) + S(1)/2)*(I*a*x + S(1))**(-I*n/S(2) + S(1)/2)/sqrt(a**S(2)*x**S(2) + S(1)), x)


def replacement5428(a, m, n, x):
    return Int(x**m*(-I*a*x + S(1))**(I*n/S(2) + S(1)/2)*(I*a*x + S(1))**(-I*n/S(2) + S(1)/2)/sqrt(a**S(2)*x**S(2) + S(1)), x)


def replacement5429(a, n, x):
    return Int((-I*a*x + S(1))**(I*n/S(2))*(I*a*x + S(1))**(-I*n/S(2)), x)


def replacement5430(a, m, n, x):
    return Int(x**m*(-I*a*x + S(1))**(I*n/S(2))*(I*a*x + S(1))**(-I*n/S(2)), x)


def replacement5431(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(S(1) + d*x/c)**p*(-I*a*x + S(1))**(I*n/S(2))*(I*a*x + S(1))**(-I*n/S(2)), x), x)


def replacement5432(a, c, d, n, p, u, x):
    return Int(u*(c + d*x)**p*(-I*a*x + S(1))**(I*n/S(2))*(I*a*x + S(1))**(-I*n/S(2)), x)


def replacement5433(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5434(a, c, d, n, p, u, x):
    return Dist((S(-1))**(n/S(2))*c**p, Int(u*(S(1) - I/(a*x))**(-I*n/S(2))*(S(1) + I/(a*x))**(I*n/S(2))*(S(1) + d/(c*x))**p, x), x)


def replacement5435(a, c, d, n, p, u, x):
    return Int(u*(c + d/x)**p*(-I*a*x + S(1))**(I*n/S(2))*(I*a*x + S(1))**(-I*n/S(2)), x)


def replacement5436(a, c, d, n, p, u, x):
    return Dist(x**p*(c + d/x)**p*(c*x/d + S(1))**(-p), Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5437(a, c, d, n, x):
    return Simp((a*x + n)*exp(n*ArcTan(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))), x)


def replacement5438(a, c, d, n, p, x):
    return Dist(S(2)*(p + S(1))*(S(2)*p + S(3))/(c*(n**S(2) + S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) + n)*exp(n*ArcTan(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))), x)


def replacement5439(a, c, d, n, x):
    return Simp(exp(n*ArcTan(a*x))/(a*c*n), x)


def replacement5440(a, c, d, n, p, x):
    return Dist(c**p, Int((a**S(2)*x**S(2) + S(1))**(-I*n/S(2) + p)*(-I*a*x + S(1))**(I*n), x), x)


def replacement5441(a, c, d, n, p, x):
    return Dist(c**p, Int((-I*a*x + S(1))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n/S(2) + p), x), x)


def replacement5442(a, c, d, n, p, x):
    return Dist(c**(I*n/S(2)), Int((c + d*x**S(2))**(-I*n/S(2) + p)*(-I*a*x + S(1))**(I*n), x), x)


def replacement5443(a, c, d, n, p, x):
    return Dist(c**(-I*n/S(2)), Int((c + d*x**S(2))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n), x), x)


def replacement5444(a, c, d, n, p, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5445(a, c, d, n, x):
    return -Simp((-a*n*x + S(1))*exp(n*ArcTan(a*x))/(d*sqrt(c + d*x**S(2))*(n**S(2) + S(1))), x)


def replacement5446(a, c, d, n, p, x):
    return -Dist(a*c*n/(S(2)*d*(p + S(1))), Int((c + d*x**S(2))**p*exp(n*ArcTan(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x))/(S(2)*d*(p + S(1))), x)


def replacement5447(a, c, d, n, p, x):
    return -Simp((c + d*x**S(2))**(p + S(1))*(-a*n*x + S(1))*exp(n*ArcTan(a*x))/(a*d*n*(n**S(2) + S(1))), x)


def replacement5448(a, c, d, n, p, x):
    return Dist((n**S(2) - S(2)*p + S(-2))/(d*(n**S(2) + S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*ArcTan(a*x)), x), x) - Simp((c + d*x**S(2))**(p + S(1))*(-S(2)*a*x*(p + S(1)) + n)*exp(n*ArcTan(a*x))/(a*d*(n**S(2) + S(4)*(p + S(1))**S(2))), x)


def replacement5449(a, c, d, m, n, p, x):
    return Dist(c**p, Int(x**m*(a**S(2)*x**S(2) + S(1))**(-I*n/S(2) + p)*(-I*a*x + S(1))**(I*n), x), x)


def replacement5450(a, c, d, m, n, p, x):
    return Dist(c**p, Int(x**m*(-I*a*x + S(1))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n/S(2) + p), x), x)


def replacement5451(a, c, d, m, n, p, x):
    return Dist(c**(I*n/S(2)), Int(x**m*(c + d*x**S(2))**(-I*n/S(2) + p)*(-I*a*x + S(1))**(I*n), x), x)


def replacement5452(a, c, d, m, n, p, x):
    return Dist(c**(-I*n/S(2)), Int(x**m*(c + d*x**S(2))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n), x), x)


def replacement5453(a, c, d, m, n, p, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(x**m*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5454(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(-I*a*x + S(1))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n/S(2) + p), x), x)


def replacement5455(a, c, d, n, p, u, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-I*a*x + S(1))**(-FracPart(p))*(I*a*x + S(1))**(-FracPart(p)), Int(u*(-I*a*x + S(1))**(I*n/S(2) + p)*(I*a*x + S(1))**(-I*n/S(2) + p), x), x)


def replacement5456(a, c, d, n, p, u, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(u*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5457(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5458(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(S(1) - I/(a*x))**p*(S(1) + I/(a*x))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5459(a, c, d, n, p, u, x):
    return Dist(x**(S(2)*p)*(c + d/x**S(2))**p*(-I*a*x + S(1))**(-p)*(I*a*x + S(1))**(-p), Int(u*x**(-S(2)*p)*(-I*a*x + S(1))**p*(I*a*x + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5460(a, c, d, n, p, u, x):
    return Dist(x**(S(2)*p)*(c + d/x**S(2))**p*(a**S(2)*x**S(2) + S(1))**(-p), Int(u*x**(-S(2)*p)*(a**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5461(a, b, c, n, x):
    return Int((-I*a*c - I*b*c*x + S(1))**(I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(-I*n/S(2)), x)


def replacement5462(a, b, c, m, n, x):
    return Dist(S(4)*I**(-m)*b**(-m + S(-1))*c**(-m + S(-1))/n, Subst(Int(x**(-S(2)*I/n)*(S(1) + x**(-S(2)*I/n))**(-m + S(-2))*(-I*a*c + S(1) - x**(-S(2)*I/n)*(I*a*c + S(1)))**m, x), x, (-I*c*(a + b*x) + S(1))**(I*n/S(2))*(I*c*(a + b*x) + S(1))**(-I*n/S(2))), x)


def replacement5463(a, b, c, d, e, m, n, x):
    return Int((d + e*x)**m*(-I*a*c - I*b*c*x + S(1))**(I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(-I*n/S(2)), x)


def replacement5464(a, b, c, d, e, n, p, u, x):
    return Dist((c/(a**S(2) + S(1)))**p, Int(u*(-I*a - I*b*x + S(1))**(I*n/S(2) + p)*(I*a + I*b*x + S(1))**(-I*n/S(2) + p), x), x)


def replacement5465(a, b, c, d, e, n, p, u, x):
    return Dist((c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p), Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*ArcTan(a*x)), x), x)


def replacement5466(a, b, c, n, u, x):
    return Int(u*exp(n*acot(a/c + b*x/c)), x)


def replacement5467(a, n, u, x):
    return Dist((S(-1))**(I*n/S(2)), Int(u*exp(-n*ArcTan(a*x)), x), x)


def replacement5468(a, n, x):
    return -Subst(Int((S(1) - I*x/a)**(I*n/S(2) + S(1)/2)*(S(1) + I*x/a)**(-I*n/S(2) + S(1)/2)/(x**S(2)*sqrt(S(1) + x**S(2)/a**S(2))), x), x, S(1)/x)


def replacement5469(a, m, n, x):
    return -Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2) + S(1)/2)*(S(1) + I*x/a)**(-I*n/S(2) + S(1)/2)/sqrt(S(1) + x**S(2)/a**S(2)), x), x, S(1)/x)


def replacement5470(a, n, x):
    return -Subst(Int((S(1) - I*x/a)**(I*n/S(2))*(S(1) + I*x/a)**(-I*n/S(2))/x**S(2), x), x, S(1)/x)


def replacement5471(a, m, n, x):
    return -Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(n/S(2))*(S(1) + I*x/a)**(-n/S(2)), x), x, S(1)/x)


def replacement5472(a, m, n, x):
    return -Dist(x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2) + S(1)/2)*(S(1) + I*x/a)**(-I*n/S(2) + S(1)/2)/sqrt(S(1) + x**S(2)/a**S(2)), x), x, S(1)/x), x)


def replacement5473(a, m, n, x):
    return -Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(n/S(2))*(S(1) + I*x/a)**(-n/S(2)), x), x, S(1)/x)


def replacement5474(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acot(a*x)), x), x)


def replacement5475(a, c, d, n, p, u, x):
    return Dist(x**(-p)*(c + d*x)**p*(c/(d*x) + S(1))**(-p), Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acot(a*x)), x), x)


def replacement5476(a, c, d, n, p, x):
    return -Dist(c**p, Subst(Int((S(1) - I*x/a)**(I*n/S(2))*(S(1) + I*x/a)**(-I*n/S(2))*(S(1) + d*x/c)**p/x**S(2), x), x, S(1)/x), x)


def replacement5477(a, c, d, m, n, p, x):
    return -Dist(c**p, Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2))*(S(1) + I*x/a)**(-I*n/S(2))*(S(1) + d*x/c)**p, x), x, S(1)/x), x)


def replacement5478(a, c, d, n, p, x):
    return Dist((S(1) + d/(c*x))**(-p)*(c + d/x)**p, Int((S(1) + d/(c*x))**p*exp(n*acot(a*x)), x), x)


def replacement5479(a, c, d, m, n, p, x):
    return -Dist(c**p*x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2))*(S(1) + I*x/a)**(-I*n/S(2))*(S(1) + d*x/c)**p, x), x, S(1)/x), x)


def replacement5480(a, c, d, n, p, u, x):
    return Dist((S(1) + d/(c*x))**(-p)*(c + d/x)**p, Int(u*(S(1) + d/(c*x))**p*exp(n*acot(a*x)), x), x)


def replacement5481(a, c, d, n, x):
    return -Simp(exp(n*acot(a*x))/(a*c*n), x)


def replacement5482(a, c, d, n, x):
    return -Simp((-a*x + n)*exp(n*acot(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))), x)


def replacement5483(a, c, d, n, p, x):
    return Dist(S(2)*(p + S(1))*(S(2)*p + S(3))/(c*(n**S(2) + S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x), x) - Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acot(a*x))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))), x)


def replacement5484(a, c, d, n, x):
    return -Simp((a*n*x + S(1))*exp(n*acot(a*x))/(a**S(2)*c*sqrt(c + d*x**S(2))*(n**S(2) + S(1))), x)


def replacement5485(a, c, d, n, p, x):
    return Dist(n*(S(2)*p + S(3))/(a*c*(n**S(2) + S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(-a*n*x + S(2)*p + S(2))*exp(n*acot(a*x))/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))), x)


def replacement5486(a, c, d, n, p, x):
    return Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acot(a*x))/(a**S(3)*c*n**S(2)*(n**S(2) + S(1))), x)


def replacement5487(a, c, d, n, p, x):
    return Dist((n**S(2) - S(2)*p + S(-2))/(a**S(2)*c*(n**S(2) + S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acot(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acot(a*x))/(a**S(3)*c*(n**S(2) + S(4)*(p + S(1))**S(2))), x)


def replacement5488(a, c, d, m, n, p, x):
    return -Dist(a**(-m + S(-1))*c**p, Subst(Int((S(1)/tan(x))**(m + S(2)*p + S(2))*exp(n*x)*cos(x)**(-S(2)*p + S(-2)), x), x, acot(a*x)), x)


def replacement5489(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x), x)


def replacement5490(a, c, d, n, p, u, x):
    return Dist(x**(-S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d*x**S(2))**p, Int(u*x**(S(2)*p)*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x), x)


def replacement5491(a, c, d, n, p, u, x):
    return Dist(c**p*(I*a)**(-S(2)*p), Int(u*x**(-S(2)*p)*(I*a*x + S(-1))**(-I*n/S(2) + p)*(I*a*x + S(1))**(I*n/S(2) + p), x), x)


def replacement5492(a, c, d, n, p, x):
    return -Dist(c**p, Subst(Int((S(1) - I*x/a)**(I*n/S(2) + p)*(S(1) + I*x/a)**(-I*n/S(2) + p)/x**S(2), x), x, S(1)/x), x)


def replacement5493(a, c, d, m, n, p, x):
    return -Dist(c**p, Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2) + p)*(S(1) + I*x/a)**(-I*n/S(2) + p), x), x, S(1)/x), x)


def replacement5494(a, c, d, m, n, p, x):
    return -Dist(c**p*x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - I*x/a)**(I*n/S(2) + p)*(S(1) + I*x/a)**(-I*n/S(2) + p), x), x, S(1)/x), x)


def replacement5495(a, c, d, n, p, u, x):
    return Dist((S(1) + S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d/x**S(2))**p, Int(u*(S(1) + S(1)/(a**S(2)*x**S(2)))**p*exp(n*acot(a*x)), x), x)


def replacement5496(a, b, c, n, u, x):
    return Dist((S(-1))**(I*n/S(2)), Int(u*exp(-n*ArcTan(c*(a + b*x))), x), x)


def replacement5497(a, b, c, n, x):
    return Dist((I*c*(a + b*x))**(I*n/S(2))*(S(1) - I/(c*(a + b*x)))**(I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(-I*n/S(2)), Int((I*a*c + I*b*c*x + S(-1))**(-I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(I*n/S(2)), x), x)


def replacement5498(a, b, c, m, n, x):
    return Dist(S(4)*I**(-m)*b**(-m + S(-1))*c**(-m + S(-1))/n, Subst(Int(x**(-S(2)*I/n)*(S(-1) + x**(-S(2)*I/n))**(-m + S(-2))*(I*a*c + S(1) + x**(-S(2)*I/n)*(-I*a*c + S(1)))**m, x), x, (S(1) - I/(c*(a + b*x)))**(I*n/S(2))*(S(1) + I/(c*(a + b*x)))**(-I*n/S(2))), x)


def replacement5499(a, b, c, d, e, m, n, x):
    return Dist((I*c*(a + b*x))**(I*n/S(2))*(S(1) - I/(c*(a + b*x)))**(I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(-I*n/S(2)), Int((d + e*x)**m*(I*a*c + I*b*c*x + S(-1))**(-I*n/S(2))*(I*a*c + I*b*c*x + S(1))**(I*n/S(2)), x), x)


def replacement5500(a, b, c, d, e, n, p, u, x):
    return Dist((c/(a**S(2) + S(1)))**p*((I*a + I*b*x + S(1))/(I*a + I*b*x))**(I*n/S(2))*((I*a + I*b*x)/(I*a + I*b*x + S(1)))**(I*n/S(2))*(-I*a - I*b*x + S(1))**(I*n/S(2))*(I*a + I*b*x + S(-1))**(-I*n/S(2)), Int(u*(-I*a - I*b*x + S(1))**(-I*n/S(2) + p)*(I*a + I*b*x + S(1))**(I*n/S(2) + p), x), x)


def replacement5501(a, b, c, d, e, n, p, u, x):
    return Dist((c + d*x + e*x**S(2))**p*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**(-p), Int(u*(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2) + S(1))**p*exp(n*acot(a*x)), x), x)


def replacement5502(a, b, c, n, u, x):
    return Int(u*exp(n*ArcTan(a/c + b*x/c)), x)


def replacement5503(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*ArcTan(x))**n, x), x, c + d*x), x)


def replacement5504(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acot(x))**n, x), x, c + d*x), x)


def replacement5505(a, b, c, d, n, x):
    return Int((a + b*ArcTan(c + d*x))**n, x)


def replacement5506(a, b, c, d, n, x):
    return Int((a + b*acot(c + d*x))**n, x)


def replacement5507(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*ArcTan(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5508(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acot(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5509(a, b, c, d, e, f, m, n, x):
    return Int((a + b*ArcTan(c + d*x))**n*(e + f*x)**m, x)


def replacement5510(a, b, c, d, e, f, m, n, x):
    return Int((a + b*acot(c + d*x))**n*(e + f*x)**m, x)


def replacement5511(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x), x)


def replacement5512(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x), x)


def replacement5513(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*ArcTan(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5514(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acot(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement5515(a, b, c, d, n, x):
    return Dist(I/S(2), Int(log(-I*a - I*b*x + S(1))/(c + d*x**n), x), x) - Dist(I/S(2), Int(log(I*a + I*b*x + S(1))/(c + d*x**n), x), x)


def replacement5516(a, b, c, d, n, x):
    return Dist(I/S(2), Int(log((a + b*x - I)/(a + b*x))/(c + d*x**n), x), x) - Dist(I/S(2), Int(log((a + b*x + I)/(a + b*x))/(c + d*x**n), x), x)


def replacement5517(a, b, c, d, n, x):
    return Int(ArcTan(a + b*x)/(c + d*x**n), x)


def replacement5518(a, b, c, d, n, x):
    return Int(acot(a + b*x)/(c + d*x**n), x)


def replacement5519(a, b, n, x):
    return -Dist(b*n, Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x*ArcTan(a + b*x**n), x)


def replacement5520(a, b, n, x):
    return Dist(b*n, Int(x**n/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x*acot(a + b*x**n), x)


def replacement5521(a, b, n, x):
    return Dist(I/S(2), Int(log(-I*a - I*b*x**n + S(1))/x, x), x) - Dist(I/S(2), Int(log(I*a + I*b*x**n + S(1))/x, x), x)


def replacement5522(a, b, n, x):
    return Dist(I/S(2), Int(log(S(1) - I/(a + b*x**n))/x, x), x) - Dist(I/S(2), Int(log(S(1) + I/(a + b*x**n))/x, x), x)


def replacement5523(a, b, m, n, x):
    return -Dist(b*n/(m + S(1)), Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x**(m + S(1))*ArcTan(a + b*x**n)/(m + S(1)), x)


def replacement5524(a, b, m, n, x):
    return Dist(b*n/(m + S(1)), Int(x**(m + n)/(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x**(m + S(1))*acot(a + b*x**n)/(m + S(1)), x)


def replacement5525(a, b, c, d, f, x):
    return Dist(I/S(2), Int(log(-I*a - I*b*f**(c + d*x) + S(1)), x), x) - Dist(I/S(2), Int(log(I*a + I*b*f**(c + d*x) + S(1)), x), x)


def replacement5526(a, b, c, d, f, x):
    return Dist(I/S(2), Int(log(S(1) - I/(a + b*f**(c + d*x))), x), x) - Dist(I/S(2), Int(log(S(1) + I/(a + b*f**(c + d*x))), x), x)


def replacement5527(a, b, c, d, f, m, x):
    return Dist(I/S(2), Int(x**m*log(-I*a - I*b*f**(c + d*x) + S(1)), x), x) - Dist(I/S(2), Int(x**m*log(I*a + I*b*f**(c + d*x) + S(1)), x), x)


def replacement5528(a, b, c, d, f, m, x):
    return Dist(I/S(2), Int(x**m*log(S(1) - I/(a + b*f**(c + d*x))), x), x) - Dist(I/S(2), Int(x**m*log(S(1) + I/(a + b*f**(c + d*x))), x), x)


def replacement5529(a, b, c, m, n, u, x):
    return Int(u*acot(a/c + b*x**n/c)**m, x)


def replacement5530(a, b, c, m, n, u, x):
    return Int(u*ArcTan(a/c + b*x**n/c)**m, x)


def replacement5531(a, b, c, x):
    return Simp(log(ArcTan(c*x/sqrt(a + b*x**S(2))))/c, x)


def replacement5532(a, b, c, x):
    return -Simp(log(acot(c*x/sqrt(a + b*x**S(2))))/c, x)


def replacement5533(a, b, c, m, x):
    return Simp(ArcTan(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))), x)


def replacement5534(a, b, c, m, x):
    return -Simp(acot(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))), x)


def replacement5535(a, b, c, d, e, m, x):
    return Dist(sqrt(a + b*x**S(2))/sqrt(d + e*x**S(2)), Int(ArcTan(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x), x)


def replacement5536(a, b, c, d, e, m, x):
    return Dist(sqrt(a + b*x**S(2))/sqrt(d + e*x**S(2)), Int(acot(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x), x)


def replacement5537(s, u, v, w, x):
    return Dist(S(1)/2, Int(u*ArcTan(v), x), x) + Dist(Pi*s/S(4), Int(u, x), x)


def replacement5538(s, u, v, w, x):
    return -Dist(S(1)/2, Int(u*ArcTan(v), x), x) + Dist(Pi*s/S(4), Int(u, x), x)


def With5539(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        tmp = InverseFunctionOfLinear(u, x)
        res = And(Not(FalseQ(tmp)), SameQ(Head(tmp), ArcTan), ZeroQ(D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2)))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement5539(n, u, v, x):

    tmp = InverseFunctionOfLinear(u, x)
    return Dist((-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n/Coefficient(Part(tmp, S(1)), x, S(1)), Subst(Int(SimplifyIntegrand((S(1)/cos(x))**(S(2)*n + S(2))*SubstForInverseFunction(u, tmp, x), x), x), x, tmp), x)


def With5540(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        tmp = InverseFunctionOfLinear(u, x)
        res = And(Not(FalseQ(tmp)), SameQ(Head(tmp), ArcCot), ZeroQ(D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2)))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement5540(n, u, v, x):

    tmp = InverseFunctionOfLinear(u, x)
    return -Dist((-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n/Coefficient(Part(tmp, S(1)), x, S(1)), Subst(Int(SimplifyIntegrand((S(1)/sin(x))**(S(2)*n + S(2))*SubstForInverseFunction(u, tmp, x), x), x), x, tmp), x)


def replacement5541(a, b, c, d, x):
    return -Dist(I*b, Int(x/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp(x*ArcTan(c + d*tan(a + b*x)), x)


def replacement5542(a, b, c, d, x):
    return Dist(I*b, Int(x/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp(x*acot(c + d*tan(a + b*x)), x)


def replacement5543(a, b, c, d, x):
    return -Dist(I*b, Int(x/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp(x*ArcTan(c + d/tan(a + b*x)), x)


def replacement5544(a, b, c, d, x):
    return Dist(I*b, Int(x/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp(x*acot(c + d/tan(a + b*x)), x)


def replacement5545(a, b, c, d, x):
    return Dist(b*(-I*c - d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c + d + (-I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(b*(I*c + d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(I*c - d + (I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*ArcTan(c + d*tan(a + b*x)), x)


def replacement5546(a, b, c, d, x):
    return -Dist(b*(-I*c - d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c + d + (-I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(b*(I*c + d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(I*c - d + (I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*acot(c + d*tan(a + b*x)), x)


def replacement5547(a, b, c, d, x):
    return -Dist(b*(-I*c + d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c - d - (-I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(b*(I*c - d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(I*c + d - (I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*ArcTan(c + d/tan(a + b*x)), x)


def replacement5548(a, b, c, d, x):
    return Dist(b*(-I*c + d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c - d - (-I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(b*(I*c - d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(I*c + d - (I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*acot(c + d/tan(a + b*x)), x)


def replacement5549(a, b, c, d, e, f, m, x):
    return -Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement5550(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement5551(a, b, c, d, e, f, m, x):
    return -Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement5552(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement5553(a, b, c, d, e, f, m, x):
    return Dist(b*(-I*c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c + d + (-I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(b*(I*c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(I*c - d + (I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement5554(a, b, c, d, e, f, m, x):
    return -Dist(b*(-I*c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c + d + (-I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(b*(I*c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(I*c - d + (I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement5555(a, b, c, d, e, f, m, x):
    return -Dist(b*(-I*c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c - d - (-I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(b*(I*c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(I*c + d - (I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement5556(a, b, c, d, e, f, m, x):
    return Dist(b*(-I*c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-I*c - d - (-I*c + d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(b*(I*c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(I*c + d - (I*c - d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement5557(a, b, x):
    return -Dist(b, Int(x/cosh(S(2)*a + S(2)*b*x), x), x) + Simp(x*ArcTan(tanh(a + b*x)), x)


def replacement5558(a, b, x):
    return Dist(b, Int(x/cosh(S(2)*a + S(2)*b*x), x), x) + Simp(x*acot(tanh(a + b*x)), x)


def replacement5559(a, b, x):
    return Dist(b, Int(x/cosh(S(2)*a + S(2)*b*x), x), x) + Simp(x*ArcTan(S(1)/tanh(a + b*x)), x)


def replacement5560(a, b, x):
    return -Dist(b, Int(x/cosh(S(2)*a + S(2)*b*x), x), x) + Simp(x*acot(S(1)/tanh(a + b*x)), x)


def replacement5561(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cosh(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5562(a, b, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cosh(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*acot(tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5563(a, b, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cosh(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(S(1)/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5564(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cosh(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*acot(S(1)/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5565(a, b, c, d, x):
    return -Dist(b, Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*ArcTan(c + d*tanh(a + b*x)), x)


def replacement5566(a, b, c, d, x):
    return Dist(b, Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*acot(c + d*tanh(a + b*x)), x)


def replacement5567(a, b, c, d, x):
    return -Dist(b, Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*ArcTan(c + d/tanh(a + b*x)), x)


def replacement5568(a, b, c, d, x):
    return Dist(b, Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*acot(c + d/tanh(a + b*x)), x)


def replacement5569(a, b, c, d, x):
    return Dist(I*b*(-c - d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) - Dist(I*b*(c + d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp(x*ArcTan(c + d*tanh(a + b*x)), x)


def replacement5570(a, b, c, d, x):
    return -Dist(I*b*(-c - d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Dist(I*b*(c + d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp(x*acot(c + d*tanh(a + b*x)), x)


def replacement5571(a, b, c, d, x):
    return -Dist(I*b*(-c - d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Dist(I*b*(c + d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp(x*ArcTan(c + d/tanh(a + b*x)), x)


def replacement5572(a, b, c, d, x):
    return Dist(I*b*(-c - d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) - Dist(I*b*(c + d + I), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp(x*acot(c + d/tanh(a + b*x)), x)


def replacement5573(a, b, c, d, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5574(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5575(a, b, c, d, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5576(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5577(a, b, c, d, e, f, m, x):
    return Dist(I*b*(-c - d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) - Dist(I*b*(c + d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5578(a, b, c, d, e, f, m, x):
    return -Dist(I*b*(-c - d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Dist(I*b*(c + d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5579(a, b, c, d, e, f, m, x):
    return -Dist(I*b*(-c - d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Dist(I*b*(c + d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp((e + f*x)**(m + S(1))*ArcTan(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5580(a, b, c, d, e, f, m, x):
    return Dist(I*b*(-c - d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) - Dist(I*b*(c + d + I)/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + I)*exp(S(2)*a + S(2)*b*x) + I), x), x) + Simp((e + f*x)**(m + S(1))*acot(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement5581(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x) + Simp(x*ArcTan(u), x)


def replacement5582(u, x):
    return Int(SimplifyIntegrand(x*D(u, x)/(u**S(2) + S(1)), x), x) + Simp(x*acot(u), x)


def replacement5583(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x), x) + Simp((a + b*ArcTan(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement5584(a, b, c, d, m, u, x):
    return Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u**S(2) + S(1)), x), x), x) + Simp((a + b*acot(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With5585(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5585(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/(u**S(2) + S(1)), x), x), x) + Dist(a + b*ArcTan(u), w, x)


def With5586(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5586(a, b, u, v, x):

    w = IntHide(v, x)
    return Dist(b, Int(SimplifyIntegrand(w*D(u, x)/(u**S(2) + S(1)), x), x), x) + Dist(a + b*acot(u), w, x)


def replacement5587(a, b, v, w, x):
    return Dist(I/S(2), Int(log(w)*log(-I*v + S(1))/(a + b*x), x), x) - Dist(I/S(2), Int(log(w)*log(I*v + S(1))/(a + b*x), x), x)


def replacement5588(v, w, x):
    return -Int(SimplifyIntegrand(x*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(x*D(v, x)*log(w)/(v**S(2) + S(1)), x), x) + Simp(x*ArcTan(v)*log(w), x)


def replacement5589(v, w, x):
    return -Int(SimplifyIntegrand(x*D(w, x)*acot(v)/w, x), x) + Int(SimplifyIntegrand(x*D(v, x)*log(w)/(v**S(2) + S(1)), x), x) + Simp(x*log(w)*acot(v), x)


def With5590(u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    z = IntHide(u, x)
    if InverseFunctionFreeQ(z, x):
        return True
    return False


def replacement5590(u, v, w, x):

    z = IntHide(u, x)
    return Dist(ArcTan(v)*log(w), z, x) - Int(SimplifyIntegrand(z*ArcTan(v)*D(w, x)/w, x), x) - Int(SimplifyIntegrand(z*D(v, x)*log(w)/(v**S(2) + S(1)), x), x)


def With5591(u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    z = IntHide(u, x)
    if InverseFunctionFreeQ(z, x):
        return True
    return False


def replacement5591(u, v, w, x):

    z = IntHide(u, x)
    return Dist(log(w)*acot(v), z, x) - Int(SimplifyIntegrand(z*D(w, x)*acot(v)/w, x), x) + Int(SimplifyIntegrand(z*D(v, x)*log(w)/(v**S(2) + S(1)), x), x)


def replacement5592(c, x):
    return -Dist(S(1)/c, Int(S(1)/(x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x), x) + Simp(x*asec(c*x), x)


def replacement5593(c, x):
    return Dist(S(1)/c, Int(S(1)/(x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))), x), x) + Simp(x*acsc(c*x), x)


def replacement5594(a, b, c, n, x):
    return Dist(S(1)/c, Subst(Int((a + b*x)**n*tan(x)/cos(x), x), x, asec(c*x)), x)


def replacement5595(a, b, c, n, x):
    return -Dist(S(1)/c, Subst(Int((a + b*x)**n/(sin(x)*tan(x)), x), x, acsc(c*x)), x)


def replacement5596(a, b, c, x):
    return -Subst(Int((a + b*acos(x/c))/x, x), x, S(1)/x)


def replacement5597(a, b, c, x):
    return -Subst(Int((a + b*asin(x/c))/x, x), x, S(1)/x)


def replacement5598(a, b, c, m, x):
    return -Dist(b/(c*(m + S(1))), Int(x**(m + S(-1))/sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x), x) + Simp(x**(m + S(1))*(a + b*asec(c*x))/(m + S(1)), x)


def replacement5599(a, b, c, m, x):
    return Dist(b/(c*(m + S(1))), Int(x**(m + S(-1))/sqrt(S(1) - S(1)/(c**S(2)*x**S(2))), x), x) + Simp(x**(m + S(1))*(a + b*acsc(c*x))/(m + S(1)), x)


def replacement5600(a, b, c, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(S(1)/cos(x))**(m + S(1))*tan(x), x), x, asec(c*x)), x)


def replacement5601(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(S(1)/sin(x))**(m + S(1))/tan(x), x), x, acsc(c*x)), x)


def replacement5602(a, b, c, m, n, x):
    return Int(x**m*(a + b*asec(c*x))**n, x)


def replacement5603(a, b, c, m, n, x):
    return Int(x**m*(a + b*acsc(c*x))**n, x)


def With5604(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c*x/sqrt(c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*asec(c*x), u, x)


def With5605(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return Dist(b*c*x/sqrt(c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*acsc(c*x), u, x)


def replacement5606(a, b, c, d, e, n, p, x):
    return -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement5607(a, b, c, d, e, n, p, x):
    return -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement5608(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5609(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5610(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5611(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5612(a, b, c, d, e, n, p, x):
    return Int((a + b*asec(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5613(a, b, c, d, e, n, p, x):
    return Int((a + b*acsc(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5614(a, b, c, d, e, p, x):
    return -Dist(b*c*x/(S(2)*e*sqrt(c**S(2)*x**S(2))*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x) + Simp((a + b*asec(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement5615(a, b, c, d, e, p, x):
    return Dist(b*c*x/(S(2)*e*sqrt(c**S(2)*x**S(2))*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x) + Simp((a + b*acsc(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With5616(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c*x/sqrt(c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*asec(c*x), u, x)


def With5617(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return Dist(b*c*x/sqrt(c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*acsc(c*x), u, x)


def replacement5618(a, b, c, d, e, m, n, p, x):
    return -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement5619(a, b, c, d, e, m, n, p, x):
    return -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement5620(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5621(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5622(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acos(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5623(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asin(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement5624(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*asec(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5625(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acsc(c*x))**n*(d + e*x**S(2))**p, x)


def replacement5626(a, b, x):
    return -Int(S(1)/(sqrt(S(1) - S(1)/(a + b*x)**S(2))*(a + b*x)), x) + Simp((a + b*x)*asec(a + b*x)/b, x)


def replacement5627(a, b, x):
    return Int(S(1)/(sqrt(S(1) - S(1)/(a + b*x)**S(2))*(a + b*x)), x) + Simp((a + b*x)*acsc(a + b*x)/b, x)


def replacement5628(a, b, n, x):
    return Dist(S(1)/b, Subst(Int(x**n*tan(x)/cos(x), x), x, asec(a + b*x)), x)


def replacement5629(a, b, n, x):
    return -Dist(S(1)/b, Subst(Int(x**n/(sin(x)*tan(x)), x), x, acsc(a + b*x)), x)


def replacement5630(a, b, x):
    return -Simp(I*PolyLog(S(2), (S(1) - sqrt(S(1) - a**S(2)))*exp(I*asec(a + b*x))/a), x) - Simp(I*PolyLog(S(2), (sqrt(S(1) - a**S(2)) + S(1))*exp(I*asec(a + b*x))/a), x) + Simp(I*PolyLog(S(2), -exp(S(2)*I*asec(a + b*x)))/S(2), x) + Simp(log(S(1) - (S(1) - sqrt(S(1) - a**S(2)))*exp(I*asec(a + b*x))/a)*asec(a + b*x), x) + Simp(log(S(1) - (sqrt(S(1) - a**S(2)) + S(1))*exp(I*asec(a + b*x))/a)*asec(a + b*x), x) - Simp(log(exp(S(2)*I*asec(a + b*x)) + S(1))*asec(a + b*x), x)


def replacement5631(a, b, x):
    return Simp(I*PolyLog(S(2), I*(S(1) - sqrt(S(1) - a**S(2)))*exp(-I*acsc(a + b*x))/a), x) + Simp(I*PolyLog(S(2), I*(sqrt(S(1) - a**S(2)) + S(1))*exp(-I*acsc(a + b*x))/a), x) + Simp(I*PolyLog(S(2), exp(S(2)*I*acsc(a + b*x)))/S(2), x) + Simp(I*acsc(a + b*x)**S(2), x) + Simp(log(S(1) - I*(S(1) - sqrt(S(1) - a**S(2)))*exp(-I*acsc(a + b*x))/a)*acsc(a + b*x), x) + Simp(log(S(1) - I*(sqrt(S(1) - a**S(2)) + S(1))*exp(-I*acsc(a + b*x))/a)*acsc(a + b*x), x) - Simp(log(S(1) - exp(S(2)*I*acsc(a + b*x)))*acsc(a + b*x), x)


def replacement5632(a, b, m, x):
    return -Dist(b**(-m + S(-1))/(m + S(1)), Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(S(1) - x**S(2)), x), x, S(1)/(a + b*x)), x) - Simp(b**(-m + S(-1))*(-b**(m + S(1))*x**(m + S(1)) + (-a)**(m + S(1)))*asec(a + b*x)/(m + S(1)), x)


def replacement5633(a, b, m, x):
    return Dist(b**(-m + S(-1))/(m + S(1)), Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(S(1) - x**S(2)), x), x, S(1)/(a + b*x)), x) - Simp(b**(-m + S(-1))*(-b**(m + S(1))*x**(m + S(1)) + (-a)**(m + S(1)))*acsc(a + b*x)/(m + S(1)), x)


def replacement5634(a, b, m, n, x):
    return Dist(b**(-m + S(-1)), Subst(Int(x**n*(-a + S(1)/cos(x))**m*tan(x)/cos(x), x), x, asec(a + b*x)), x)


def replacement5635(a, b, m, n, x):
    return -Dist(b**(-m + S(-1)), Subst(Int(x**n*(-a + S(1)/sin(x))**m/(sin(x)*tan(x)), x), x, acsc(a + b*x)), x)


def replacement5636(a, b, c, m, n, u, x):
    return Int(u*acos(a/c + b*x**n/c)**m, x)


def replacement5637(a, b, c, m, n, u, x):
    return Int(u*asin(a/c + b*x**n/c)**m, x)


def replacement5638(a, b, c, f, n, u, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + S(1)/(b*cos(x))))*tan(x)/cos(x), x), x, asec(a + b*x)), x)


def replacement5639(a, b, c, f, n, u, x):
    return -Dist(S(1)/b, Subst(Int(f**(c*x**n)*ReplaceAll(u, Rule(x, -a/b + S(1)/(b*sin(x))))/(sin(x)*tan(x)), x), x, acsc(a + b*x)), x)


def replacement5640(u, x):
    return -Dist(u/sqrt(u**S(2)), Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Simp(x*asec(u), x)


def replacement5641(u, x):
    return Dist(u/sqrt(u**S(2)), Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Simp(x*acsc(u), x)


def replacement5642(a, b, c, d, m, u, x):
    return -Dist(b*u/(d*(m + S(1))*sqrt(u**S(2))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Simp((a + b*asec(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement5643(a, b, c, d, m, u, x):
    return Dist(b*u/(d*(m + S(1))*sqrt(u**S(2))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Simp((a + b*acsc(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With5644(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5644(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b*u/sqrt(u**S(2)), Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Dist(a + b*asec(u), w, x)


def With5645(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement5645(a, b, u, v, x):

    w = IntHide(v, x)
    return Dist(b*u/sqrt(u**S(2)), Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(u**S(2) + S(-1))), x), x), x) + Dist(a + b*acsc(u), w, x)
