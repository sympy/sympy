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


def inverse_hyperbolic():
    from sympy.integrals.rubi.constraints import cons89, cons90, cons2, cons3, cons8, cons91, cons1581, cons4, cons150, cons68, cons29, cons19, cons64, cons1736, cons1737, cons1738, cons1780, cons270, cons50, cons1894, cons1895, cons1896, cons1897, cons733, cons654, cons734, cons656, cons586, cons1740, cons1898, cons130, cons1739, cons340, cons165, cons40, cons349, cons139, cons232, cons669, cons5, cons1741, cons1742, cons963, cons1899, cons1743, cons1900, cons1745, cons1744, cons1746, cons1572, cons1901, cons338, cons1902, cons149, cons127, cons210, cons56, cons244, cons1748, cons1749, cons488, cons164, cons96, cons95, cons274, cons1750, cons20, cons168, cons276, cons1751, cons1752, cons21, cons240, cons239, cons1753, cons248, cons1754, cons1903, cons1755, cons1756, cons1904, cons1757, cons1758, cons1905, cons211, cons927, cons466, cons86, cons1759, cons1760, cons721, cons170, cons1761, cons1762, cons269, cons719, cons1763, cons1610, cons14, cons152, cons1200, cons1275, cons1362, cons1832, cons1765, cons36, cons37, cons38, cons1764, cons1906, cons167, cons1444, cons1767, cons1766, cons1768, cons1769, cons530, cons1232, cons1771, cons1772, cons1907, cons1908, cons87, cons806, cons33, cons342, cons1909, cons1910, cons1911, cons1778, cons1045, cons1779, cons1499, cons13, cons1781, cons1782, cons1783, cons1784, cons242, cons243, cons148, cons1785, cons1512, cons1786, cons1154, cons321, cons1787, cons1788, cons1789, cons1790, cons1912, cons1913, cons1914, cons1915, cons1795, cons1796, cons1916, cons1798, cons603, cons1799, cons263, cons1917, cons1484, cons1443, cons1918, cons1252, cons1919, cons1920, cons1804, cons1805, cons1921, cons745, cons179, cons119, cons1922, cons25, cons1923, cons1924, cons1925, cons1926, cons1927, cons676, cons1928, cons1929, cons1930, cons996, cons1582, cons1820, cons1931, cons1932, cons1933, cons1934, cons1935, cons1936, cons1937, cons1826, cons975, cons1938, cons1829, cons1939, cons1940, cons1096, cons1833, cons1834, cons1835, cons1836, cons1941, cons385, cons810, cons1588, cons820, cons465, cons1942, cons1943, cons1944, cons1945, cons1946, cons69, cons1947, cons1948, cons1949, cons1950, cons1849, cons1951, cons1952, cons1953, cons1954, cons1856, cons180, cons1857, cons1858, cons1301, cons1955, cons1956, cons1957, cons1958


    pattern6087 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons89, cons90)
    rule6087 = ReplacementRule(pattern6087, replacement6087)

    pattern6088 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons89, cons90)
    rule6088 = ReplacementRule(pattern6088, replacement6088)

    pattern6089 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons89, cons91)
    rule6089 = ReplacementRule(pattern6089, replacement6089)

    pattern6090 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons89, cons91)
    rule6090 = ReplacementRule(pattern6090, replacement6090)

    pattern6091 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule6091 = ReplacementRule(pattern6091, replacement6091)

    pattern6092 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule6092 = ReplacementRule(pattern6092, replacement6092)

    pattern6093 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons150)
    rule6093 = ReplacementRule(pattern6093, replacement6093)

    pattern6094 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/x_, x_), cons2, cons3, cons8, cons150)
    rule6094 = ReplacementRule(pattern6094, replacement6094)

    pattern6095 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons150, cons68)
    rule6095 = ReplacementRule(pattern6095, replacement6095)

    pattern6096 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons150, cons68)
    rule6096 = ReplacementRule(pattern6096, replacement6096)

    pattern6097 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons90)
    rule6097 = ReplacementRule(pattern6097, replacement6097)

    pattern6098 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons90)
    rule6098 = ReplacementRule(pattern6098, replacement6098)

    pattern6099 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1736)
    rule6099 = ReplacementRule(pattern6099, replacement6099)

    pattern6100 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1736)
    rule6100 = ReplacementRule(pattern6100, replacement6100)

    pattern6101 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1737)
    rule6101 = ReplacementRule(pattern6101, replacement6101)

    pattern6102 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons64, cons89, cons1737)
    rule6102 = ReplacementRule(pattern6102, replacement6102)

    pattern6103 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons64)
    rule6103 = ReplacementRule(pattern6103, replacement6103)

    pattern6104 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons64)
    rule6104 = ReplacementRule(pattern6104, replacement6104)

    pattern6105 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule6105 = ReplacementRule(pattern6105, replacement6105)

    pattern6106 = Pattern(Integral((x_*WC('d', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1738)
    rule6106 = ReplacementRule(pattern6106, replacement6106)

    pattern6107 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270)
    rule6107 = ReplacementRule(pattern6107, replacement6107)

    pattern6108 = Pattern(Integral(S(1)/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons1896, cons1897)
    rule6108 = ReplacementRule(pattern6108, replacement6108)

    pattern6109 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons270, cons586)
    rule6109 = ReplacementRule(pattern6109, replacement6109)

    pattern6110 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons1896, cons1897, cons586)
    rule6110 = ReplacementRule(pattern6110, replacement6110)

    pattern6111 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1740)
    rule6111 = ReplacementRule(pattern6111, replacement6111)

    pattern6112 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons1898)
    rule6112 = ReplacementRule(pattern6112, replacement6112)

    pattern6113 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons130)
    rule6113 = ReplacementRule(pattern6113, With6113)

    pattern6114 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule6114 = ReplacementRule(pattern6114, With6114)

    pattern6115 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons165, cons40)
    rule6115 = ReplacementRule(pattern6115, replacement6115)

    pattern6116 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule6116 = ReplacementRule(pattern6116, replacement6116)

    pattern6117 = Pattern(Integral(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons89, cons90)
    rule6117 = ReplacementRule(pattern6117, replacement6117)

    pattern6118 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons90, cons165)
    rule6118 = ReplacementRule(pattern6118, replacement6118)

    pattern6119 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons340, cons90, cons165, cons349)
    rule6119 = ReplacementRule(pattern6119, replacement6119)

    pattern6120 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons340, cons90, cons165)
    rule6120 = ReplacementRule(pattern6120, replacement6120)

    pattern6121 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons89, cons90)
    rule6121 = ReplacementRule(pattern6121, replacement6121)

    pattern6122 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/((d1_ + x_*WC('e1', S(1)))**(S(3)/2)*(d2_ + x_*WC('e2', S(1)))**(S(3)/2)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons89, cons90)
    rule6122 = ReplacementRule(pattern6122, replacement6122)

    pattern6123 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons90, cons139, cons40)
    rule6123 = ReplacementRule(pattern6123, replacement6123)

    pattern6124 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons340, cons90, cons139, cons232)
    rule6124 = ReplacementRule(pattern6124, replacement6124)

    pattern6125 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons340, cons90, cons139, cons232, cons669)
    rule6125 = ReplacementRule(pattern6125, replacement6125)

    pattern6126 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons340, cons90, cons139, cons232)
    rule6126 = ReplacementRule(pattern6126, replacement6126)

    pattern6127 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150)
    rule6127 = ReplacementRule(pattern6127, replacement6127)

    pattern6128 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule6128 = ReplacementRule(pattern6128, replacement6128)

    pattern6129 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons91, cons40)
    rule6129 = ReplacementRule(pattern6129, replacement6129)

    pattern6130 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons91)
    rule6130 = ReplacementRule(pattern6130, replacement6130)

    pattern6131 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons5, cons1894, cons1895, cons89, cons91, cons349)
    rule6131 = ReplacementRule(pattern6131, replacement6131)

    pattern6132 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons5, cons1894, cons1895, cons89, cons91)
    rule6132 = ReplacementRule(pattern6132, replacement6132)

    pattern6133 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1741, cons1742)
    rule6133 = ReplacementRule(pattern6133, replacement6133)

    pattern6134 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons130)
    rule6134 = ReplacementRule(pattern6134, replacement6134)

    pattern6135 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons963, cons1899)
    rule6135 = ReplacementRule(pattern6135, replacement6135)

    pattern6136 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons1741, cons1743)
    rule6136 = ReplacementRule(pattern6136, replacement6136)

    pattern6137 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons1741, cons1898)
    rule6137 = ReplacementRule(pattern6137, replacement6137)

    pattern6138 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1900, cons1745)
    rule6138 = ReplacementRule(pattern6138, With6138)

    pattern6139 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1744, cons1745)
    rule6139 = ReplacementRule(pattern6139, With6139)

    pattern6140 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1900, cons40, cons1746)
    rule6140 = ReplacementRule(pattern6140, replacement6140)

    pattern6141 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1744, cons40, cons1746)
    rule6141 = ReplacementRule(pattern6141, replacement6141)

    pattern6142 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule6142 = ReplacementRule(pattern6142, replacement6142)

    pattern6143 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons40)
    rule6143 = ReplacementRule(pattern6143, replacement6143)

    pattern6144 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons5, cons1901)
    rule6144 = ReplacementRule(pattern6144, replacement6144)

    pattern6145 = Pattern(Integral((d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons338, cons1902, cons149)
    rule6145 = ReplacementRule(pattern6145, replacement6145)

    pattern6146 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1739, cons149)
    rule6146 = ReplacementRule(pattern6146, replacement6146)

    pattern6147 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150)
    rule6147 = ReplacementRule(pattern6147, replacement6147)

    pattern6148 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule6148 = ReplacementRule(pattern6148, replacement6148)

    pattern6149 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons56, cons40)
    rule6149 = ReplacementRule(pattern6149, replacement6149)

    pattern6150 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1780, cons89, cons90, cons56)
    rule6150 = ReplacementRule(pattern6150, replacement6150)

    pattern6151 = Pattern(Integral(x_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons5, cons1894, cons1895, cons89, cons90, cons56, cons669)
    rule6151 = ReplacementRule(pattern6151, replacement6151)

    pattern6152 = Pattern(Integral(x_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons5, cons1894, cons1895, cons89, cons90, cons56)
    rule6152 = ReplacementRule(pattern6152, replacement6152)

    pattern6153 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons150)
    rule6153 = ReplacementRule(pattern6153, replacement6153)

    pattern6154 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule6154 = ReplacementRule(pattern6154, replacement6154)

    pattern6155 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons90, cons244, cons68, cons40)
    rule6155 = ReplacementRule(pattern6155, replacement6155)

    pattern6156 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1780, cons89, cons90, cons244, cons68)
    rule6156 = ReplacementRule(pattern6156, replacement6156)

    pattern6157 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons5, cons1894, cons1895, cons89, cons90, cons244, cons68, cons669)
    rule6157 = ReplacementRule(pattern6157, replacement6157)

    pattern6158 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons5, cons1894, cons1895, cons89, cons90, cons244, cons68)
    rule6158 = ReplacementRule(pattern6158, replacement6158)

    pattern6159 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons130)
    rule6159 = ReplacementRule(pattern6159, replacement6159)

    pattern6160 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons130)
    rule6160 = ReplacementRule(pattern6160, replacement6160)

    pattern6161 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons130, cons1748)
    rule6161 = ReplacementRule(pattern6161, replacement6161)

    pattern6162 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons130, cons1748)
    rule6162 = ReplacementRule(pattern6162, replacement6162)

    pattern6163 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons130)
    rule6163 = ReplacementRule(pattern6163, With6163)

    pattern6164 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons130)
    rule6164 = ReplacementRule(pattern6164, With6164)

    pattern6165 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons349, cons1749, cons488, cons270)
    rule6165 = ReplacementRule(pattern6165, With6165)

    pattern6166 = Pattern(Integral(x_**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons349, cons1749, cons488, cons1896, cons1897)
    rule6166 = ReplacementRule(pattern6166, With6166)

    pattern6167 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons963, cons1749)
    rule6167 = ReplacementRule(pattern6167, With6167)

    pattern6168 = Pattern(Integral(x_**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons963, cons1749)
    rule6168 = ReplacementRule(pattern6168, With6168)

    pattern6169 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons164, cons90, cons165, cons96, cons40)
    rule6169 = ReplacementRule(pattern6169, replacement6169)

    pattern6170 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons95, cons90, cons96)
    rule6170 = ReplacementRule(pattern6170, replacement6170)

    pattern6171 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons95, cons90, cons96)
    rule6171 = ReplacementRule(pattern6171, replacement6171)

    pattern6172 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons164, cons90, cons165, cons96)
    rule6172 = ReplacementRule(pattern6172, replacement6172)

    pattern6173 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons164, cons90, cons165, cons96, cons349)
    rule6173 = ReplacementRule(pattern6173, replacement6173)

    pattern6174 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons165, cons274, cons40, cons1750)
    rule6174 = ReplacementRule(pattern6174, replacement6174)

    pattern6175 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons89, cons90, cons274, cons1750)
    rule6175 = ReplacementRule(pattern6175, replacement6175)

    pattern6176 = Pattern(Integral((x_*WC('f', S(1)))**m_*sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons89, cons90, cons274, cons1750)
    rule6176 = ReplacementRule(pattern6176, replacement6176)

    pattern6177 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons340, cons90, cons165, cons274, cons1750)
    rule6177 = ReplacementRule(pattern6177, replacement6177)

    pattern6178 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons340, cons90, cons165, cons274, cons349, cons1750)
    rule6178 = ReplacementRule(pattern6178, replacement6178)

    pattern6179 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons96, cons20, cons40)
    rule6179 = ReplacementRule(pattern6179, replacement6179)

    pattern6180 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1780, cons95, cons90, cons96, cons20)
    rule6180 = ReplacementRule(pattern6180, replacement6180)

    pattern6181 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons5, cons1894, cons1895, cons95, cons90, cons96, cons20, cons669)
    rule6181 = ReplacementRule(pattern6181, replacement6181)

    pattern6182 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons5, cons1894, cons1895, cons95, cons90, cons96, cons20)
    rule6182 = ReplacementRule(pattern6182, replacement6182)

    pattern6183 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons95, cons90, cons139, cons168, cons40)
    rule6183 = ReplacementRule(pattern6183, replacement6183)

    pattern6184 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons164, cons90, cons139, cons168)
    rule6184 = ReplacementRule(pattern6184, replacement6184)

    pattern6185 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons164, cons90, cons139, cons168, cons669)
    rule6185 = ReplacementRule(pattern6185, replacement6185)

    pattern6186 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons164, cons90, cons139, cons149, cons168)
    rule6186 = ReplacementRule(pattern6186, replacement6186)

    pattern6187 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1739, cons340, cons90, cons139, cons276, cons40)
    rule6187 = ReplacementRule(pattern6187, replacement6187)

    pattern6188 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons340, cons90, cons139, cons276, cons1751)
    rule6188 = ReplacementRule(pattern6188, replacement6188)

    pattern6189 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons340, cons90, cons139, cons276, cons1752, cons669)
    rule6189 = ReplacementRule(pattern6189, replacement6189)

    pattern6190 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons340, cons90, cons139, cons276, cons1751)
    rule6190 = ReplacementRule(pattern6190, replacement6190)

    pattern6191 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons95, cons90, cons168, cons20)
    rule6191 = ReplacementRule(pattern6191, replacement6191)

    pattern6192 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons95, cons90, cons168, cons20)
    rule6192 = ReplacementRule(pattern6192, replacement6192)

    pattern6193 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1780, cons270, cons150, cons20)
    rule6193 = ReplacementRule(pattern6193, replacement6193)

    pattern6194 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1894, cons1895, cons150, cons1896, cons1897, cons20)
    rule6194 = ReplacementRule(pattern6194, replacement6194)

    pattern6195 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons270, cons21)
    rule6195 = ReplacementRule(pattern6195, replacement6195)

    pattern6196 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons1896, cons1897, cons21)
    rule6196 = ReplacementRule(pattern6196, replacement6196)

    pattern6197 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons89, cons90, cons1740, cons1752)
    rule6197 = ReplacementRule(pattern6197, replacement6197)

    pattern6198 = Pattern(Integral((x_*WC('f', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons89, cons90, cons1898, cons1752)
    rule6198 = ReplacementRule(pattern6198, replacement6198)

    pattern6199 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1739, cons95, cons90, cons168, cons240, cons40, cons20)
    rule6199 = ReplacementRule(pattern6199, replacement6199)

    pattern6200 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons1780, cons95, cons90, cons168, cons240, cons20)
    rule6200 = ReplacementRule(pattern6200, replacement6200)

    pattern6201 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons5, cons1894, cons1895, cons95, cons90, cons168, cons240, cons20, cons669)
    rule6201 = ReplacementRule(pattern6201, replacement6201)

    pattern6202 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons5, cons1894, cons1895, cons95, cons90, cons168, cons240, cons20)
    rule6202 = ReplacementRule(pattern6202, replacement6202)

    pattern6203 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1739, cons89, cons91, cons239, cons40)
    rule6203 = ReplacementRule(pattern6203, replacement6203)

    pattern6204 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1780, cons89, cons91, cons239)
    rule6204 = ReplacementRule(pattern6204, replacement6204)

    pattern6205 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons5, cons1894, cons1895, cons89, cons91, cons239, cons349)
    rule6205 = ReplacementRule(pattern6205, replacement6205)

    pattern6206 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons5, cons1894, cons1895, cons89, cons91, cons239)
    rule6206 = ReplacementRule(pattern6206, replacement6206)

    pattern6207 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons89, cons91, cons270)
    rule6207 = ReplacementRule(pattern6207, replacement6207)

    pattern6208 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons89, cons91, cons1896, cons1897)
    rule6208 = ReplacementRule(pattern6208, replacement6208)

    pattern6209 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1780, cons89, cons91, cons1740)
    rule6209 = ReplacementRule(pattern6209, replacement6209)

    pattern6210 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons1894, cons1895, cons89, cons91, cons1898)
    rule6210 = ReplacementRule(pattern6210, replacement6210)

    pattern6211 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1739, cons89, cons91, cons20, cons1753, cons130)
    rule6211 = ReplacementRule(pattern6211, replacement6211)

    pattern6212 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1780, cons89, cons91, cons20, cons1753, cons1741)
    rule6212 = ReplacementRule(pattern6212, replacement6212)

    pattern6213 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons1894, cons1895, cons89, cons91, cons20, cons1753, cons963)
    rule6213 = ReplacementRule(pattern6213, replacement6213)

    pattern6214 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons248, cons1754, cons64, cons1742)
    rule6214 = ReplacementRule(pattern6214, replacement6214)

    pattern6215 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons130, cons64)
    rule6215 = ReplacementRule(pattern6215, replacement6215)

    pattern6216 = Pattern(Integral(x_**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons669, cons1754, cons64, cons1899)
    rule6216 = ReplacementRule(pattern6216, replacement6216)

    pattern6217 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons248, cons1754, cons64, cons1743)
    rule6217 = ReplacementRule(pattern6217, replacement6217)

    pattern6218 = Pattern(Integral(x_**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons1894, cons1895, cons248, cons1754, cons64, cons1903)
    rule6218 = ReplacementRule(pattern6218, replacement6218)

    pattern6219 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1780, cons270, cons963, cons1755, cons20, cons1756)
    rule6219 = ReplacementRule(pattern6219, replacement6219)

    pattern6220 = Pattern(Integral((x_*WC('f', S(1)))**m_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons4, cons1894, cons1895, cons1896, cons1897, cons963, cons1755, cons20, cons1756)
    rule6220 = ReplacementRule(pattern6220, replacement6220)

    pattern6221 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1744, cons68, cons1904)
    rule6221 = ReplacementRule(pattern6221, replacement6221)

    pattern6222 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1900, cons56)
    rule6222 = ReplacementRule(pattern6222, replacement6222)

    pattern6223 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1744, cons56)
    rule6223 = ReplacementRule(pattern6223, replacement6223)

    pattern6224 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1900, cons40, cons1757)
    rule6224 = ReplacementRule(pattern6224, With6224)

    pattern6225 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1744, cons40, cons1757)
    rule6225 = ReplacementRule(pattern6225, With6225)

    pattern6226 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1900, cons150, cons40, cons20)
    rule6226 = ReplacementRule(pattern6226, replacement6226)

    pattern6227 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1744, cons150, cons40, cons20)
    rule6227 = ReplacementRule(pattern6227, replacement6227)

    pattern6228 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1758)
    rule6228 = ReplacementRule(pattern6228, replacement6228)

    pattern6229 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons40)
    rule6229 = ReplacementRule(pattern6229, replacement6229)

    pattern6230 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d1_ + x_*WC('e1', S(1)))**WC('p', S(1))*(d2_ + x_*WC('e2', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons19, cons4, cons5, cons1905)
    rule6230 = ReplacementRule(pattern6230, replacement6230)

    pattern6231 = Pattern(Integral((x_*WC('h', S(1)))**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons338, cons1902, cons149)
    rule6231 = ReplacementRule(pattern6231, replacement6231)

    pattern6232 = Pattern(Integral((x_*WC('f', S(1)))**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1739, cons149)
    rule6232 = ReplacementRule(pattern6232, replacement6232)

    pattern6233 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons150)
    rule6233 = ReplacementRule(pattern6233, replacement6233)

    pattern6234 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons150)
    rule6234 = ReplacementRule(pattern6234, replacement6234)

    pattern6235 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150, cons68)
    rule6235 = ReplacementRule(pattern6235, replacement6235)

    pattern6236 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150, cons68)
    rule6236 = ReplacementRule(pattern6236, replacement6236)

    pattern6237 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons64, cons89, cons91)
    rule6237 = ReplacementRule(pattern6237, replacement6237)

    pattern6238 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons64, cons89, cons91)
    rule6238 = ReplacementRule(pattern6238, replacement6238)

    pattern6239 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons64)
    rule6239 = ReplacementRule(pattern6239, replacement6239)

    pattern6240 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons64)
    rule6240 = ReplacementRule(pattern6240, replacement6240)

    pattern6241 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons927)
    rule6241 = ReplacementRule(pattern6241, With6241)

    pattern6242 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons927)
    rule6242 = ReplacementRule(pattern6242, With6242)

    pattern6243 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons927)
    rule6243 = ReplacementRule(pattern6243, replacement6243)

    pattern6244 = Pattern(Integral(Px_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons927)
    rule6244 = ReplacementRule(pattern6244, replacement6244)

    pattern6245 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons927)
    rule6245 = ReplacementRule(pattern6245, With6245)

    pattern6246 = Pattern(Integral(Px_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons927)
    rule6246 = ReplacementRule(pattern6246, With6246)

    pattern6247 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons466, cons86, cons1759)
    rule6247 = ReplacementRule(pattern6247, With6247)

    pattern6248 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons466, cons86, cons1759)
    rule6248 = ReplacementRule(pattern6248, With6248)

    pattern6249 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons466, cons1760)
    rule6249 = ReplacementRule(pattern6249, With6249)

    pattern6250 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_*(x_**S(2)*WC('h', S(1)) + x_*WC('g', S(1)) + WC('f', S(0)))**WC('p', S(1))/(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons466, cons1760)
    rule6250 = ReplacementRule(pattern6250, With6250)

    pattern6251 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons927, cons150, cons20)
    rule6251 = ReplacementRule(pattern6251, replacement6251)

    pattern6252 = Pattern(Integral(Px_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons927, cons150, cons20)
    rule6252 = ReplacementRule(pattern6252, replacement6252)

    pattern6253 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons721, cons270, cons170, cons1761)
    rule6253 = ReplacementRule(pattern6253, With6253)

    pattern6254 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons721, cons1896, cons1897, cons170, cons1761)
    rule6254 = ReplacementRule(pattern6254, With6254)

    pattern6255 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons669, cons270, cons150, cons170, cons1762)
    rule6255 = ReplacementRule(pattern6255, replacement6255)

    pattern6256 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons669, cons1896, cons1897, cons150, cons170, cons1762)
    rule6256 = ReplacementRule(pattern6256, replacement6256)

    pattern6257 = Pattern(Integral(sqrt(d_ + x_**S(2)*WC('e', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons270, cons150, cons269)
    rule6257 = ReplacementRule(pattern6257, replacement6257)

    pattern6258 = Pattern(Integral(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons1896, cons1897, cons150, cons269)
    rule6258 = ReplacementRule(pattern6258, replacement6258)

    pattern6259 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons963, cons270, cons150)
    rule6259 = ReplacementRule(pattern6259, replacement6259)

    pattern6260 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons963, cons1896, cons1897, cons150)
    rule6260 = ReplacementRule(pattern6260, replacement6260)

    pattern6261 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons719, cons270, cons150, cons269)
    rule6261 = ReplacementRule(pattern6261, replacement6261)

    pattern6262 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons719, cons1896, cons1897, cons150, cons269)
    rule6262 = ReplacementRule(pattern6262, replacement6262)

    pattern6263 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons270, cons170, cons89, cons91)
    rule6263 = ReplacementRule(pattern6263, replacement6263)

    pattern6264 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons1896, cons1897, cons170, cons89, cons91)
    rule6264 = ReplacementRule(pattern6264, replacement6264)

    pattern6265 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1780, cons20, cons270, cons1763)
    rule6265 = ReplacementRule(pattern6265, replacement6265)

    pattern6266 = Pattern(Integral((f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons4, cons1894, cons1895, cons20, cons1896, cons1897, cons1763)
    rule6266 = ReplacementRule(pattern6266, replacement6266)

    pattern6267 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1780, cons20, cons721, cons270, cons150)
    rule6267 = ReplacementRule(pattern6267, replacement6267)

    pattern6268 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons1894, cons1895, cons20, cons721, cons1896, cons1897, cons150)
    rule6268 = ReplacementRule(pattern6268, replacement6268)

    pattern6269 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1780, cons20, cons349, cons1740)
    rule6269 = ReplacementRule(pattern6269, replacement6269)

    pattern6270 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1739, cons20, cons349)
    rule6270 = ReplacementRule(pattern6270, replacement6270)

    pattern6271 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons4, cons1894, cons1895, cons20, cons349, cons1898)
    rule6271 = ReplacementRule(pattern6271, replacement6271)

    pattern6272 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons1780, cons270, cons150)
    rule6272 = ReplacementRule(pattern6272, replacement6272)

    pattern6273 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1)))/(sqrt(d1_ + x_*WC('e1', S(1)))*sqrt(d2_ + x_*WC('e2', S(1)))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons211, cons19, cons1894, cons1895, cons1896, cons1897, cons150)
    rule6273 = ReplacementRule(pattern6273, replacement6273)

    pattern6274 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons1780, cons349, cons1740)
    rule6274 = ReplacementRule(pattern6274, replacement6274)

    pattern6275 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons1739, cons349)
    rule6275 = ReplacementRule(pattern6275, replacement6275)

    pattern6276 = Pattern(Integral((d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*log((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))*WC('h', S(1))), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons211, cons19, cons4, cons1894, cons1895, cons349, cons1898)
    rule6276 = ReplacementRule(pattern6276, replacement6276)

    pattern6277 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1610)
    rule6277 = ReplacementRule(pattern6277, With6277)

    pattern6278 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1610)
    rule6278 = ReplacementRule(pattern6278, With6278)

    pattern6279 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons20)
    rule6279 = ReplacementRule(pattern6279, replacement6279)

    pattern6280 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons20)
    rule6280 = ReplacementRule(pattern6280, replacement6280)

    pattern6281 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons14, CustomConstraint(With6281))
    rule6281 = ReplacementRule(pattern6281, replacement6281)

    pattern6282 = Pattern(Integral(u_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons14, CustomConstraint(With6282))
    rule6282 = ReplacementRule(pattern6282, replacement6282)

    pattern6283 = Pattern(Integral(Px_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons927, cons1780, cons349, CustomConstraint(With6283))
    rule6283 = ReplacementRule(pattern6283, replacement6283)

    pattern6284 = Pattern(Integral(Px_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons4, cons927, cons1894, cons1895, cons349, CustomConstraint(With6284))
    rule6284 = ReplacementRule(pattern6284, replacement6284)

    pattern6285 = Pattern(Integral((f_ + (d_ + x_**S(2)*WC('e', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons927, cons1780, cons963, cons152, CustomConstraint(With6285))
    rule6285 = ReplacementRule(pattern6285, replacement6285)

    pattern6286 = Pattern(Integral((f_ + (d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*WC('g', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*WC('Px', S(1)), x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons127, cons210, cons927, cons1894, cons1895, cons963, cons152, CustomConstraint(With6286))
    rule6286 = ReplacementRule(pattern6286, replacement6286)

    pattern6287 = Pattern(Integral(RFx_*asinh(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons1200, cons150, CustomConstraint(With6287))
    rule6287 = ReplacementRule(pattern6287, replacement6287)

    pattern6288 = Pattern(Integral(RFx_*acosh(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons1200, cons150, CustomConstraint(With6288))
    rule6288 = ReplacementRule(pattern6288, replacement6288)

    pattern6289 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons1200, cons150)
    rule6289 = ReplacementRule(pattern6289, replacement6289)

    pattern6290 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons1200, cons150)
    rule6290 = ReplacementRule(pattern6290, replacement6290)

    pattern6291 = Pattern(Integral(RFx_*(d_ + x_**S(2)*WC('e', S(1)))**p_*asinh(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons29, cons50, cons1200, cons150, cons1780, cons349, CustomConstraint(With6291))
    rule6291 = ReplacementRule(pattern6291, replacement6291)

    pattern6292 = Pattern(Integral(RFx_*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_*acosh(x_*WC('c', S(1)))**WC('n', S(1)), x_), cons8, cons733, cons654, cons734, cons656, cons1200, cons150, cons1894, cons1895, cons349, CustomConstraint(With6292))
    rule6292 = ReplacementRule(pattern6292, replacement6292)

    pattern6293 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons1200, cons150, cons1780, cons349)
    rule6293 = ReplacementRule(pattern6293, replacement6293)

    pattern6294 = Pattern(Integral(RFx_*(a_ + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*(d1_ + x_*WC('e1', S(1)))**p_*(d2_ + x_*WC('e2', S(1)))**p_, x_), cons2, cons3, cons8, cons733, cons654, cons734, cons656, cons1200, cons150, cons1894, cons1895, cons349)
    rule6294 = ReplacementRule(pattern6294, replacement6294)

    pattern6295 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule6295 = ReplacementRule(pattern6295, replacement6295)

    pattern6296 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_*WC('c', S(1))))**WC('n', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule6296 = ReplacementRule(pattern6296, replacement6296)

    pattern6297 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1275)
    rule6297 = ReplacementRule(pattern6297, replacement6297)

    pattern6298 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1275)
    rule6298 = ReplacementRule(pattern6298, replacement6298)

    pattern6299 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule6299 = ReplacementRule(pattern6299, replacement6299)

    pattern6300 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule6300 = ReplacementRule(pattern6300, replacement6300)

    pattern6301 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1832, cons1765)
    rule6301 = ReplacementRule(pattern6301, replacement6301)

    pattern6302 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1764, cons1765)
    rule6302 = ReplacementRule(pattern6302, replacement6302)

    pattern6303 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1832, cons1765)
    rule6303 = ReplacementRule(pattern6303, replacement6303)

    pattern6304 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1764, cons1765)
    rule6304 = ReplacementRule(pattern6304, replacement6304)

    pattern6305 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1906)
    rule6305 = ReplacementRule(pattern6305, replacement6305)

    pattern6306 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1906, cons89, cons167)
    rule6306 = ReplacementRule(pattern6306, replacement6306)

    pattern6307 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1906)
    rule6307 = ReplacementRule(pattern6307, replacement6307)

    pattern6308 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons1906)
    rule6308 = ReplacementRule(pattern6308, replacement6308)

    pattern6309 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**(S(-3)/2), x_), cons2, cons3, cons8, cons29, cons1906)
    rule6309 = ReplacementRule(pattern6309, replacement6309)

    pattern6310 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**(S(-2)), x_), cons2, cons3, cons8, cons29, cons1906)
    rule6310 = ReplacementRule(pattern6310, replacement6310)

    pattern6311 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asinh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1906, cons89, cons91, cons1444)
    rule6311 = ReplacementRule(pattern6311, replacement6311)

    pattern6312 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule6312 = ReplacementRule(pattern6312, replacement6312)

    pattern6313 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule6313 = ReplacementRule(pattern6313, replacement6313)

    pattern6314 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons167)
    rule6314 = ReplacementRule(pattern6314, replacement6314)

    pattern6315 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule6315 = ReplacementRule(pattern6315, replacement6315)

    pattern6316 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule6316 = ReplacementRule(pattern6316, replacement6316)

    pattern6317 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1))), x_), cons2, cons3, cons29, cons1767)
    rule6317 = ReplacementRule(pattern6317, replacement6317)

    pattern6318 = Pattern(Integral(S(1)/sqrt(WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1))), x_), cons2, cons3, cons29, cons1767)
    rule6318 = ReplacementRule(pattern6318, replacement6318)

    pattern6319 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-3)/2), x_), cons2, cons3, cons29, cons1767)
    rule6319 = ReplacementRule(pattern6319, replacement6319)

    pattern6320 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-3)/2), x_), cons2, cons3, cons29, cons1767)
    rule6320 = ReplacementRule(pattern6320, replacement6320)

    pattern6321 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(1)))**(S(-2)), x_), cons2, cons3, cons29, cons1767)
    rule6321 = ReplacementRule(pattern6321, replacement6321)

    pattern6322 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(x_**S(2)*WC('d', S(1)) + S(-1)))**(S(-2)), x_), cons2, cons3, cons29, cons1767)
    rule6322 = ReplacementRule(pattern6322, replacement6322)

    pattern6323 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acosh(c_ + x_**S(2)*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons1766, cons89, cons91, cons1444)
    rule6323 = ReplacementRule(pattern6323, replacement6323)

    pattern6324 = Pattern(Integral(asinh(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), cons2, cons5, cons150)
    rule6324 = ReplacementRule(pattern6324, replacement6324)

    pattern6325 = Pattern(Integral(acosh(x_**p_*WC('a', S(1)))**WC('n', S(1))/x_, x_), cons2, cons5, cons150)
    rule6325 = ReplacementRule(pattern6325, replacement6325)

    pattern6326 = Pattern(Integral(WC('u', S(1))*asinh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6326 = ReplacementRule(pattern6326, replacement6326)

    pattern6327 = Pattern(Integral(WC('u', S(1))*acosh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6327 = ReplacementRule(pattern6327, replacement6327)

    pattern6328 = Pattern(Integral(asinh(sqrt(x_**S(2)*WC('b', S(1)) + S(-1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(-1)), x_), cons3, cons4, cons1769)
    rule6328 = ReplacementRule(pattern6328, replacement6328)

    pattern6329 = Pattern(Integral(acosh(sqrt(x_**S(2)*WC('b', S(1)) + S(1)))**WC('n', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + S(1)), x_), cons3, cons4, cons1769)
    rule6329 = ReplacementRule(pattern6329, replacement6329)

    pattern6330 = Pattern(Integral(f_**(WC('c', S(1))*asinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))), x_), cons2, cons3, cons8, cons127, cons150)
    rule6330 = ReplacementRule(pattern6330, replacement6330)

    pattern6331 = Pattern(Integral(f_**(WC('c', S(1))*acosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))), x_), cons2, cons3, cons8, cons127, cons150)
    rule6331 = ReplacementRule(pattern6331, replacement6331)

    pattern6332 = Pattern(Integral(f_**(WC('c', S(1))*asinh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*x_**WC('m', S(1)), x_), cons2, cons3, cons8, cons127, cons530)
    rule6332 = ReplacementRule(pattern6332, replacement6332)

    pattern6333 = Pattern(Integral(f_**(WC('c', S(1))*acosh(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)))*x_**WC('m', S(1)), x_), cons2, cons3, cons8, cons127, cons530)
    rule6333 = ReplacementRule(pattern6333, replacement6333)

    pattern6334 = Pattern(Integral(asinh(u_), x_), cons1232, cons1771)
    rule6334 = ReplacementRule(pattern6334, replacement6334)

    pattern6335 = Pattern(Integral(acosh(u_), x_), cons1232, cons1771)
    rule6335 = ReplacementRule(pattern6335, replacement6335)

    pattern6336 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asinh(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule6336 = ReplacementRule(pattern6336, replacement6336)

    pattern6337 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acosh(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule6337 = ReplacementRule(pattern6337, replacement6337)

    pattern6338 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asinh(u_)), x_), cons2, cons3, cons1232, cons1907, CustomConstraint(With6338))
    rule6338 = ReplacementRule(pattern6338, replacement6338)

    pattern6339 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acosh(u_)), x_), cons2, cons3, cons1232, cons1908, CustomConstraint(With6339))
    rule6339 = ReplacementRule(pattern6339, replacement6339)

    pattern6340 = Pattern(Integral(exp(WC('n', S(1))*asinh(u_)), x_), cons87, cons806)
    rule6340 = ReplacementRule(pattern6340, replacement6340)

    pattern6341 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*asinh(u_)), x_), cons33, cons87, cons806)
    rule6341 = ReplacementRule(pattern6341, replacement6341)

    pattern6342 = Pattern(Integral(exp(WC('n', S(1))*acosh(u_)), x_), cons87, cons806)
    rule6342 = ReplacementRule(pattern6342, replacement6342)

    pattern6343 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acosh(u_)), x_), cons33, cons87, cons806)
    rule6343 = ReplacementRule(pattern6343, replacement6343)

    pattern6344 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons150)
    rule6344 = ReplacementRule(pattern6344, replacement6344)

    pattern6345 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons150)
    rule6345 = ReplacementRule(pattern6345, replacement6345)

    pattern6346 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons342)
    rule6346 = ReplacementRule(pattern6346, replacement6346)

    pattern6347 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons342)
    rule6347 = ReplacementRule(pattern6347, replacement6347)

    pattern6348 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150)
    rule6348 = ReplacementRule(pattern6348, replacement6348)

    pattern6349 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150)
    rule6349 = ReplacementRule(pattern6349, replacement6349)

    pattern6350 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(d_ + x_*WC('e', S(1))), x_), cons8, cons29, cons50, cons1910, cons1911)
    rule6350 = ReplacementRule(pattern6350, replacement6350)

    pattern6351 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule6351 = ReplacementRule(pattern6351, replacement6351)

    pattern6352 = Pattern(Integral(acoth(x_*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule6352 = ReplacementRule(pattern6352, replacement6352)

    pattern6353 = Pattern(Integral((a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule6353 = ReplacementRule(pattern6353, replacement6353)

    pattern6354 = Pattern(Integral((a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule6354 = ReplacementRule(pattern6354, replacement6354)

    pattern6355 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6355 = ReplacementRule(pattern6355, replacement6355)

    pattern6356 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6356 = ReplacementRule(pattern6356, replacement6356)

    pattern6357 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/x_, x_), cons2, cons3, cons8, cons87, cons167)
    rule6357 = ReplacementRule(pattern6357, replacement6357)

    pattern6358 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/x_, x_), cons2, cons3, cons8, cons87, cons167)
    rule6358 = ReplacementRule(pattern6358, replacement6358)

    pattern6359 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons19, cons87, cons167, cons68)
    rule6359 = ReplacementRule(pattern6359, replacement6359)

    pattern6360 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons19, cons87, cons167, cons68)
    rule6360 = ReplacementRule(pattern6360, replacement6360)

    pattern6361 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons466)
    rule6361 = ReplacementRule(pattern6361, replacement6361)

    pattern6362 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons466)
    rule6362 = ReplacementRule(pattern6362, replacement6362)

    pattern6363 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons342)
    rule6363 = ReplacementRule(pattern6363, replacement6363)

    pattern6364 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons342)
    rule6364 = ReplacementRule(pattern6364, replacement6364)

    pattern6365 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150, cons33, cons170)
    rule6365 = ReplacementRule(pattern6365, replacement6365)

    pattern6366 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150, cons33, cons170)
    rule6366 = ReplacementRule(pattern6366, replacement6366)

    pattern6367 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150)
    rule6367 = ReplacementRule(pattern6367, replacement6367)

    pattern6368 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150)
    rule6368 = ReplacementRule(pattern6368, replacement6368)

    pattern6369 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150, cons33, cons96)
    rule6369 = ReplacementRule(pattern6369, replacement6369)

    pattern6370 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1909, cons150, cons33, cons96)
    rule6370 = ReplacementRule(pattern6370, replacement6370)

    pattern6371 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1779)
    rule6371 = ReplacementRule(pattern6371, replacement6371)

    pattern6372 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1779)
    rule6372 = ReplacementRule(pattern6372, replacement6372)

    pattern6373 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6373 = ReplacementRule(pattern6373, replacement6373)

    pattern6374 = Pattern(Integral(x_**WC('m', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6374 = ReplacementRule(pattern6374, replacement6374)

    pattern6375 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons165)
    rule6375 = ReplacementRule(pattern6375, replacement6375)

    pattern6376 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons165)
    rule6376 = ReplacementRule(pattern6376, replacement6376)

    pattern6377 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons165, cons167)
    rule6377 = ReplacementRule(pattern6377, replacement6377)

    pattern6378 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons165, cons167)
    rule6378 = ReplacementRule(pattern6378, replacement6378)

    pattern6379 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1739)
    rule6379 = ReplacementRule(pattern6379, replacement6379)

    pattern6380 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1739)
    rule6380 = ReplacementRule(pattern6380, replacement6380)

    pattern6381 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons586)
    rule6381 = ReplacementRule(pattern6381, replacement6381)

    pattern6382 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons586)
    rule6382 = ReplacementRule(pattern6382, replacement6382)

    pattern6383 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule6383 = ReplacementRule(pattern6383, replacement6383)

    pattern6384 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule6384 = ReplacementRule(pattern6384, replacement6384)

    pattern6385 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons270)
    rule6385 = ReplacementRule(pattern6385, replacement6385)

    pattern6386 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons270)
    rule6386 = ReplacementRule(pattern6386, replacement6386)

    pattern6387 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons1740)
    rule6387 = ReplacementRule(pattern6387, replacement6387)

    pattern6388 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons1740)
    rule6388 = ReplacementRule(pattern6388, replacement6388)

    pattern6389 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6389 = ReplacementRule(pattern6389, replacement6389)

    pattern6390 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6390 = ReplacementRule(pattern6390, replacement6390)

    pattern6391 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739)
    rule6391 = ReplacementRule(pattern6391, replacement6391)

    pattern6392 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739)
    rule6392 = ReplacementRule(pattern6392, replacement6392)

    pattern6393 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons139, cons232)
    rule6393 = ReplacementRule(pattern6393, replacement6393)

    pattern6394 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons139, cons232)
    rule6394 = ReplacementRule(pattern6394, replacement6394)

    pattern6395 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons167)
    rule6395 = ReplacementRule(pattern6395, replacement6395)

    pattern6396 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons167)
    rule6396 = ReplacementRule(pattern6396, replacement6396)

    pattern6397 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons139, cons167, cons232)
    rule6397 = ReplacementRule(pattern6397, replacement6397)

    pattern6398 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons139, cons167, cons232)
    rule6398 = ReplacementRule(pattern6398, replacement6398)

    pattern6399 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons139, cons91)
    rule6399 = ReplacementRule(pattern6399, replacement6399)

    pattern6400 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons340, cons139, cons91)
    rule6400 = ReplacementRule(pattern6400, replacement6400)

    pattern6401 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1781, cons1742)
    rule6401 = ReplacementRule(pattern6401, replacement6401)

    pattern6402 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1781, cons1743)
    rule6402 = ReplacementRule(pattern6402, replacement6402)

    pattern6403 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1781, cons40)
    rule6403 = ReplacementRule(pattern6403, replacement6403)

    pattern6404 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons1781, cons149)
    rule6404 = ReplacementRule(pattern6404, replacement6404)

    pattern6405 = Pattern(Integral(atanh(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule6405 = ReplacementRule(pattern6405, replacement6405)

    pattern6406 = Pattern(Integral(acoth(x_*WC('c', S(1)))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons8, cons29, cons50, cons1778)
    rule6406 = ReplacementRule(pattern6406, replacement6406)

    pattern6407 = Pattern(Integral((a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule6407 = ReplacementRule(pattern6407, replacement6407)

    pattern6408 = Pattern(Integral((a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule6408 = ReplacementRule(pattern6408, replacement6408)

    pattern6409 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1782)
    rule6409 = ReplacementRule(pattern6409, With6409)

    pattern6410 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1782)
    rule6410 = ReplacementRule(pattern6410, With6410)

    pattern6411 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150)
    rule6411 = ReplacementRule(pattern6411, replacement6411)

    pattern6412 = Pattern(Integral((d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150)
    rule6412 = ReplacementRule(pattern6412, replacement6412)

    pattern6413 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule6413 = ReplacementRule(pattern6413, replacement6413)

    pattern6414 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule6414 = ReplacementRule(pattern6414, replacement6414)

    pattern6415 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons168)
    rule6415 = ReplacementRule(pattern6415, replacement6415)

    pattern6416 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons168)
    rule6416 = ReplacementRule(pattern6416, replacement6416)

    pattern6417 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons96)
    rule6417 = ReplacementRule(pattern6417, replacement6417)

    pattern6418 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons95, cons90, cons96)
    rule6418 = ReplacementRule(pattern6418, replacement6418)

    pattern6419 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule6419 = ReplacementRule(pattern6419, replacement6419)

    pattern6420 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150)
    rule6420 = ReplacementRule(pattern6420, replacement6420)

    pattern6421 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons342, cons586)
    rule6421 = ReplacementRule(pattern6421, replacement6421)

    pattern6422 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons342, cons586)
    rule6422 = ReplacementRule(pattern6422, replacement6422)

    pattern6423 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons168)
    rule6423 = ReplacementRule(pattern6423, replacement6423)

    pattern6424 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons168)
    rule6424 = ReplacementRule(pattern6424, replacement6424)

    pattern6425 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6425 = ReplacementRule(pattern6425, replacement6425)

    pattern6426 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6426 = ReplacementRule(pattern6426, replacement6426)

    pattern6427 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons96)
    rule6427 = ReplacementRule(pattern6427, replacement6427)

    pattern6428 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons96)
    rule6428 = ReplacementRule(pattern6428, replacement6428)

    pattern6429 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons89, cons91)
    rule6429 = ReplacementRule(pattern6429, replacement6429)

    pattern6430 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons89, cons91)
    rule6430 = ReplacementRule(pattern6430, replacement6430)

    pattern6431 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons20, cons1783)
    rule6431 = ReplacementRule(pattern6431, replacement6431)

    pattern6432 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons20, cons1783)
    rule6432 = ReplacementRule(pattern6432, replacement6432)

    pattern6433 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons56)
    rule6433 = ReplacementRule(pattern6433, replacement6433)

    pattern6434 = Pattern(Integral(x_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons56)
    rule6434 = ReplacementRule(pattern6434, replacement6434)

    pattern6435 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons91, cons1444)
    rule6435 = ReplacementRule(pattern6435, replacement6435)

    pattern6436 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons91, cons1444)
    rule6436 = ReplacementRule(pattern6436, replacement6436)

    pattern6437 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons139, cons1784)
    rule6437 = ReplacementRule(pattern6437, replacement6437)

    pattern6438 = Pattern(Integral(x_**S(2)*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons13, cons139, cons1784)
    rule6438 = ReplacementRule(pattern6438, replacement6438)

    pattern6439 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6439 = ReplacementRule(pattern6439, replacement6439)

    pattern6440 = Pattern(Integral(x_**S(2)*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6440 = ReplacementRule(pattern6440, replacement6440)

    pattern6441 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons242, cons13, cons139)
    rule6441 = ReplacementRule(pattern6441, replacement6441)

    pattern6442 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons242, cons13, cons139)
    rule6442 = ReplacementRule(pattern6442, replacement6442)

    pattern6443 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons242, cons340, cons139, cons167)
    rule6443 = ReplacementRule(pattern6443, replacement6443)

    pattern6444 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons242, cons340, cons139, cons167)
    rule6444 = ReplacementRule(pattern6444, replacement6444)

    pattern6445 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1739, cons242, cons89, cons91)
    rule6445 = ReplacementRule(pattern6445, replacement6445)

    pattern6446 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1739, cons242, cons89, cons91)
    rule6446 = ReplacementRule(pattern6446, replacement6446)

    pattern6447 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1739, cons244, cons89, cons90, cons68)
    rule6447 = ReplacementRule(pattern6447, replacement6447)

    pattern6448 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1739, cons244, cons89, cons90, cons68)
    rule6448 = ReplacementRule(pattern6448, replacement6448)

    pattern6449 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons243)
    rule6449 = ReplacementRule(pattern6449, replacement6449)

    pattern6450 = Pattern(Integral(x_**m_*sqrt(d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons243)
    rule6450 = ReplacementRule(pattern6450, replacement6450)

    pattern6451 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons150, cons40, cons148)
    rule6451 = ReplacementRule(pattern6451, replacement6451)

    pattern6452 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons150, cons40, cons148)
    rule6452 = ReplacementRule(pattern6452, replacement6452)

    pattern6453 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons13, cons165, cons150, cons1785)
    rule6453 = ReplacementRule(pattern6453, replacement6453)

    pattern6454 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1739, cons13, cons165, cons150, cons1785)
    rule6454 = ReplacementRule(pattern6454, replacement6454)

    pattern6455 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons168)
    rule6455 = ReplacementRule(pattern6455, replacement6455)

    pattern6456 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons168)
    rule6456 = ReplacementRule(pattern6456, replacement6456)

    pattern6457 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule6457 = ReplacementRule(pattern6457, replacement6457)

    pattern6458 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons270)
    rule6458 = ReplacementRule(pattern6458, replacement6458)

    pattern6459 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons270)
    rule6459 = ReplacementRule(pattern6459, replacement6459)

    pattern6460 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**n_/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons270)
    rule6460 = ReplacementRule(pattern6460, replacement6460)

    pattern6461 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons1740)
    rule6461 = ReplacementRule(pattern6461, replacement6461)

    pattern6462 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons150, cons1740)
    rule6462 = ReplacementRule(pattern6462, replacement6462)

    pattern6463 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6463 = ReplacementRule(pattern6463, replacement6463)

    pattern6464 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/(x_**S(2)*sqrt(d_ + x_**S(2)*WC('e', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90)
    rule6464 = ReplacementRule(pattern6464, replacement6464)

    pattern6465 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons96, cons1512)
    rule6465 = ReplacementRule(pattern6465, replacement6465)

    pattern6466 = Pattern(Integral(x_**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))/sqrt(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons95, cons90, cons96, cons1512)
    rule6466 = ReplacementRule(pattern6466, replacement6466)

    pattern6467 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons1786, cons139, cons168, cons1154)
    rule6467 = ReplacementRule(pattern6467, replacement6467)

    pattern6468 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons1786, cons139, cons168, cons1154)
    rule6468 = ReplacementRule(pattern6468, replacement6468)

    pattern6469 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons1786, cons139, cons269, cons1154)
    rule6469 = ReplacementRule(pattern6469, replacement6469)

    pattern6470 = Pattern(Integral(x_**m_*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons1786, cons139, cons269, cons1154)
    rule6470 = ReplacementRule(pattern6470, replacement6470)

    pattern6471 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons164, cons139, cons91, cons321)
    rule6471 = ReplacementRule(pattern6471, replacement6471)

    pattern6472 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons164, cons139, cons91, cons321)
    rule6472 = ReplacementRule(pattern6472, replacement6472)

    pattern6473 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons64, cons1787, cons1742)
    rule6473 = ReplacementRule(pattern6473, replacement6473)

    pattern6474 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons64, cons1787, cons1743)
    rule6474 = ReplacementRule(pattern6474, replacement6474)

    pattern6475 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons64, cons1787, cons40)
    rule6475 = ReplacementRule(pattern6475, replacement6475)

    pattern6476 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**p_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons64, cons1787, cons149)
    rule6476 = ReplacementRule(pattern6476, replacement6476)

    pattern6477 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6477 = ReplacementRule(pattern6477, replacement6477)

    pattern6478 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6478 = ReplacementRule(pattern6478, replacement6478)

    pattern6479 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule6479 = ReplacementRule(pattern6479, With6479)

    pattern6480 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule6480 = ReplacementRule(pattern6480, With6480)

    pattern6481 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1789)
    rule6481 = ReplacementRule(pattern6481, replacement6481)

    pattern6482 = Pattern(Integral(x_**WC('m', S(1))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons1789)
    rule6482 = ReplacementRule(pattern6482, replacement6482)

    pattern6483 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1790)
    rule6483 = ReplacementRule(pattern6483, replacement6483)

    pattern6484 = Pattern(Integral(x_**WC('m', S(1))*(a_ + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(d_ + x_**S(2)*WC('e', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1790)
    rule6484 = ReplacementRule(pattern6484, replacement6484)

    pattern6485 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6485 = ReplacementRule(pattern6485, replacement6485)

    pattern6486 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6486 = ReplacementRule(pattern6486, replacement6486)

    pattern6487 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1912)
    rule6487 = ReplacementRule(pattern6487, replacement6487)

    pattern6488 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1912)
    rule6488 = ReplacementRule(pattern6488, replacement6488)

    pattern6489 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*atanh(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1913)
    rule6489 = ReplacementRule(pattern6489, replacement6489)

    pattern6490 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*acoth(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1913)
    rule6490 = ReplacementRule(pattern6490, replacement6490)

    pattern6491 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1914)
    rule6491 = ReplacementRule(pattern6491, replacement6491)

    pattern6492 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1914)
    rule6492 = ReplacementRule(pattern6492, replacement6492)

    pattern6493 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1915)
    rule6493 = ReplacementRule(pattern6493, replacement6493)

    pattern6494 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*log(u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons89, cons90, cons1915)
    rule6494 = ReplacementRule(pattern6494, replacement6494)

    pattern6495 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons1912)
    rule6495 = ReplacementRule(pattern6495, replacement6495)

    pattern6496 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons1912)
    rule6496 = ReplacementRule(pattern6496, replacement6496)

    pattern6497 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons1913)
    rule6497 = ReplacementRule(pattern6497, replacement6497)

    pattern6498 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*PolyLog(p_, u_)/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons1739, cons89, cons90, cons1913)
    rule6498 = ReplacementRule(pattern6498, replacement6498)

    pattern6499 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('e', S(1)))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))), x_), cons2, cons3, cons8, cons29, cons50, cons1739)
    rule6499 = ReplacementRule(pattern6499, replacement6499)

    pattern6500 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('n', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons152, cons1795)
    rule6500 = ReplacementRule(pattern6500, replacement6500)

    pattern6501 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**WC('m', S(1))/(d_ + x_**S(2)*WC('e', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons1739, cons152, cons1796)
    rule6501 = ReplacementRule(pattern6501, replacement6501)

    pattern6502 = Pattern(Integral(atanh(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons8, cons29, cons87, cons1916)
    rule6502 = ReplacementRule(pattern6502, replacement6502)

    pattern6503 = Pattern(Integral(acoth(x_*WC('a', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons8, cons29, cons87, cons1916)
    rule6503 = ReplacementRule(pattern6503, replacement6503)

    pattern6504 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1798)
    rule6504 = ReplacementRule(pattern6504, replacement6504)

    pattern6505 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1798)
    rule6505 = ReplacementRule(pattern6505, replacement6505)

    pattern6506 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons603)
    rule6506 = ReplacementRule(pattern6506, replacement6506)

    pattern6507 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons603)
    rule6507 = ReplacementRule(pattern6507, replacement6507)

    pattern6508 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1799)
    rule6508 = ReplacementRule(pattern6508, With6508)

    pattern6509 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1799)
    rule6509 = ReplacementRule(pattern6509, With6509)

    pattern6510 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons20, cons263)
    rule6510 = ReplacementRule(pattern6510, With6510)

    pattern6511 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))*(WC('d', S(0)) + WC('e', S(1))*log(x_**S(2)*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons20, cons263)
    rule6511 = ReplacementRule(pattern6511, With6511)

    pattern6512 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*atanh(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1917)
    rule6512 = ReplacementRule(pattern6512, replacement6512)

    pattern6513 = Pattern(Integral(x_*(WC('a', S(0)) + WC('b', S(1))*acoth(x_*WC('c', S(1))))**S(2)*(WC('d', S(0)) + WC('e', S(1))*log(f_ + x_**S(2)*WC('g', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1917)
    rule6513 = ReplacementRule(pattern6513, replacement6513)

    pattern6514 = Pattern(Integral(exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons1484)
    rule6514 = ReplacementRule(pattern6514, replacement6514)

    pattern6515 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons19, cons1484)
    rule6515 = ReplacementRule(pattern6515, replacement6515)

    pattern6516 = Pattern(Integral(exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons4, cons1443)
    rule6516 = ReplacementRule(pattern6516, replacement6516)

    pattern6517 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons19, cons4, cons1443)
    rule6517 = ReplacementRule(pattern6517, replacement6517)

    pattern6518 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1918, cons1252, cons248)
    rule6518 = ReplacementRule(pattern6518, replacement6518)

    pattern6519 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons50, cons127, cons19, cons5, cons1918, cons1252, cons1919, cons248)
    rule6519 = ReplacementRule(pattern6519, replacement6519)

    pattern6520 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1920, cons1804)
    rule6520 = ReplacementRule(pattern6520, replacement6520)

    pattern6521 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1920, cons1805)
    rule6521 = ReplacementRule(pattern6521, replacement6521)

    pattern6522 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1921, cons40)
    rule6522 = ReplacementRule(pattern6522, replacement6522)

    pattern6523 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1921, cons149, cons745, cons179)
    rule6523 = ReplacementRule(pattern6523, replacement6523)

    pattern6524 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1921, cons149, cons745, cons119)
    rule6524 = ReplacementRule(pattern6524, replacement6524)

    pattern6525 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1921, cons149)
    rule6525 = ReplacementRule(pattern6525, replacement6525)

    pattern6526 = Pattern(Integral(exp(n_*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1922, cons25)
    rule6526 = ReplacementRule(pattern6526, replacement6526)

    pattern6527 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons13, cons139, cons25, cons1923, cons248)
    rule6527 = ReplacementRule(pattern6527, replacement6527)

    pattern6528 = Pattern(Integral(exp(WC('n', S(1))*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924)
    rule6528 = ReplacementRule(pattern6528, replacement6528)

    pattern6529 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1922, cons40, cons1925, cons1926)
    rule6529 = ReplacementRule(pattern6529, replacement6529)

    pattern6530 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1922, cons40, cons1927, cons1926)
    rule6530 = ReplacementRule(pattern6530, replacement6530)

    pattern6531 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1804)
    rule6531 = ReplacementRule(pattern6531, replacement6531)

    pattern6532 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1922, cons1805, cons676)
    rule6532 = ReplacementRule(pattern6532, replacement6532)

    pattern6533 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1922, cons1805, cons1928)
    rule6533 = ReplacementRule(pattern6533, replacement6533)

    pattern6534 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1805)
    rule6534 = ReplacementRule(pattern6534, replacement6534)

    pattern6535 = Pattern(Integral(x_*exp(n_*atanh(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1922, cons25)
    rule6535 = ReplacementRule(pattern6535, replacement6535)

    pattern6536 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons13, cons139, cons25, cons248)
    rule6536 = ReplacementRule(pattern6536, replacement6536)

    pattern6537 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1929, cons25)
    rule6537 = ReplacementRule(pattern6537, replacement6537)

    pattern6538 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons13, cons139, cons25, cons1923, cons248)
    rule6538 = ReplacementRule(pattern6538, replacement6538)

    pattern6539 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1922, cons1804, cons1925, cons1926)
    rule6539 = ReplacementRule(pattern6539, replacement6539)

    pattern6540 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1922, cons1804, cons1927, cons1926)
    rule6540 = ReplacementRule(pattern6540, replacement6540)

    pattern6541 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1922, cons1804)
    rule6541 = ReplacementRule(pattern6541, replacement6541)

    pattern6542 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1922, cons1805, cons676)
    rule6542 = ReplacementRule(pattern6542, replacement6542)

    pattern6543 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons5, cons1922, cons1805, cons1928)
    rule6543 = ReplacementRule(pattern6543, replacement6543)

    pattern6544 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1922, cons1805, cons1924)
    rule6544 = ReplacementRule(pattern6544, replacement6544)

    pattern6545 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1804)
    rule6545 = ReplacementRule(pattern6545, replacement6545)

    pattern6546 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1805, cons745)
    rule6546 = ReplacementRule(pattern6546, replacement6546)

    pattern6547 = Pattern(Integral(u_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1805, cons1924)
    rule6547 = ReplacementRule(pattern6547, replacement6547)

    pattern6548 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1930, cons40)
    rule6548 = ReplacementRule(pattern6548, replacement6548)

    pattern6549 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1930, cons149, cons745, cons179)
    rule6549 = ReplacementRule(pattern6549, replacement6549)

    pattern6550 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(n_*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons149, cons745, cons119)
    rule6550 = ReplacementRule(pattern6550, replacement6550)

    pattern6551 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*atanh(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons149, cons1924)
    rule6551 = ReplacementRule(pattern6551, replacement6551)

    pattern6552 = Pattern(Integral(exp(WC('n', S(1))*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule6552 = ReplacementRule(pattern6552, replacement6552)

    pattern6553 = Pattern(Integral(x_**m_*exp(n_*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons86, cons89, cons996)
    rule6553 = ReplacementRule(pattern6553, replacement6553)

    pattern6554 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*atanh((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule6554 = ReplacementRule(pattern6554, replacement6554)

    pattern6555 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1820, cons1931, cons1932)
    rule6555 = ReplacementRule(pattern6555, replacement6555)

    pattern6556 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*atanh(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1820, cons1931, cons1933)
    rule6556 = ReplacementRule(pattern6556, replacement6556)

    pattern6557 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*atanh(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule6557 = ReplacementRule(pattern6557, replacement6557)

    pattern6558 = Pattern(Integral(WC('u', S(1))*exp(n_*acoth(x_*WC('a', S(1)))), x_), cons2, cons745)
    rule6558 = ReplacementRule(pattern6558, replacement6558)

    pattern6559 = Pattern(Integral(exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons1484)
    rule6559 = ReplacementRule(pattern6559, replacement6559)

    pattern6560 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons1484, cons20)
    rule6560 = ReplacementRule(pattern6560, replacement6560)

    pattern6561 = Pattern(Integral(exp(n_*acoth(x_*WC('a', S(1)))), x_), cons2, cons4, cons25)
    rule6561 = ReplacementRule(pattern6561, replacement6561)

    pattern6562 = Pattern(Integral(x_**WC('m', S(1))*exp(n_*acoth(x_*WC('a', S(1)))), x_), cons2, cons4, cons25, cons20)
    rule6562 = ReplacementRule(pattern6562, replacement6562)

    pattern6563 = Pattern(Integral(x_**m_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons19, cons1484, cons21)
    rule6563 = ReplacementRule(pattern6563, replacement6563)

    pattern6564 = Pattern(Integral(x_**m_*exp(n_*acoth(x_*WC('a', S(1)))), x_), cons2, cons19, cons4, cons25, cons21)
    rule6564 = ReplacementRule(pattern6564, replacement6564)

    pattern6565 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1918, cons1934, cons1924)
    rule6565 = ReplacementRule(pattern6565, replacement6565)

    pattern6566 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1920, cons1924, cons40)
    rule6566 = ReplacementRule(pattern6566, replacement6566)

    pattern6567 = Pattern(Integral((c_ + x_*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1920, cons1924, cons149)
    rule6567 = ReplacementRule(pattern6567, replacement6567)

    pattern6568 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1935, cons1252, cons1919, cons248)
    rule6568 = ReplacementRule(pattern6568, replacement6568)

    pattern6569 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons5, cons1935, cons1252, cons20, cons1936, cons248)
    rule6569 = ReplacementRule(pattern6569, replacement6569)

    pattern6570 = Pattern(Integral((c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1921, cons1924, cons1804)
    rule6570 = ReplacementRule(pattern6570, replacement6570)

    pattern6571 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1921, cons1924, cons1804, cons20)
    rule6571 = ReplacementRule(pattern6571, replacement6571)

    pattern6572 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_)**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1921, cons1924, cons1804, cons21)
    rule6572 = ReplacementRule(pattern6572, replacement6572)

    pattern6573 = Pattern(Integral((c_ + WC('d', S(1))/x_)**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1921, cons1924, cons1805)
    rule6573 = ReplacementRule(pattern6573, replacement6573)

    pattern6574 = Pattern(Integral(exp(WC('n', S(1))*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924)
    rule6574 = ReplacementRule(pattern6574, replacement6574)

    pattern6575 = Pattern(Integral(exp(n_*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1922, cons25)
    rule6575 = ReplacementRule(pattern6575, replacement6575)

    pattern6576 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons13, cons139, cons232, cons1923, cons1937)
    rule6576 = ReplacementRule(pattern6576, replacement6576)

    pattern6577 = Pattern(Integral(x_*exp(n_*acoth(x_*WC('a', S(1))))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons2, cons8, cons29, cons4, cons1922, cons25)
    rule6577 = ReplacementRule(pattern6577, replacement6577)

    pattern6578 = Pattern(Integral(x_*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons13, cons1826, cons232, cons1923, cons1937)
    rule6578 = ReplacementRule(pattern6578, replacement6578)

    pattern6579 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons1929, cons975)
    rule6579 = ReplacementRule(pattern6579, replacement6579)

    pattern6580 = Pattern(Integral(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons13, cons1826, cons1938, cons1923, cons1937)
    rule6580 = ReplacementRule(pattern6580, replacement6580)

    pattern6581 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**S(2)*WC('d', S(1)))**p_*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons20, cons13, cons1829, cons40)
    rule6581 = ReplacementRule(pattern6581, replacement6581)

    pattern6582 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons1922, cons1924, cons40)
    rule6582 = ReplacementRule(pattern6582, replacement6582)

    pattern6583 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1922, cons1924, cons149)
    rule6583 = ReplacementRule(pattern6583, replacement6583)

    pattern6584 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons1924, cons1804, cons1939)
    rule6584 = ReplacementRule(pattern6584, replacement6584)

    pattern6585 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons1924, cons1804, cons1940)
    rule6585 = ReplacementRule(pattern6585, replacement6585)

    pattern6586 = Pattern(Integral(x_**WC('m', S(1))*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons1924, cons1804, cons1940, cons20)
    rule6586 = ReplacementRule(pattern6586, replacement6586)

    pattern6587 = Pattern(Integral(x_**m_*(c_ + WC('d', S(1))/x_**S(2))**WC('p', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons1930, cons1924, cons1804, cons1940, cons21)
    rule6587 = ReplacementRule(pattern6587, replacement6587)

    pattern6588 = Pattern(Integral((c_ + WC('d', S(1))/x_**S(2))**p_*WC('u', S(1))*exp(WC('n', S(1))*acoth(x_*WC('a', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1930, cons1924, cons1805)
    rule6588 = ReplacementRule(pattern6588, replacement6588)

    pattern6589 = Pattern(Integral(WC('u', S(1))*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons745)
    rule6589 = ReplacementRule(pattern6589, replacement6589)

    pattern6590 = Pattern(Integral(exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons4, cons1924)
    rule6590 = ReplacementRule(pattern6590, replacement6590)

    pattern6591 = Pattern(Integral(x_**m_*exp(n_*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons86, cons89, cons996)
    rule6591 = ReplacementRule(pattern6591, replacement6591)

    pattern6592 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*exp(WC('n', S(1))*acoth((a_ + x_*WC('b', S(1)))*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1924)
    rule6592 = ReplacementRule(pattern6592, replacement6592)

    pattern6593 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1924, cons1820, cons1931, cons1932)
    rule6593 = ReplacementRule(pattern6593, replacement6593)

    pattern6594 = Pattern(Integral((c_ + x_**S(2)*WC('e', S(1)) + x_*WC('d', S(1)))**WC('p', S(1))*WC('u', S(1))*exp(WC('n', S(1))*acoth(a_ + x_*WC('b', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1924, cons1820, cons1931, cons1933)
    rule6594 = ReplacementRule(pattern6594, replacement6594)

    pattern6595 = Pattern(Integral(WC('u', S(1))*exp(WC('n', S(1))*acoth(WC('c', S(1))/(x_*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons4, cons1581)
    rule6595 = ReplacementRule(pattern6595, replacement6595)

    pattern6596 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150)
    rule6596 = ReplacementRule(pattern6596, replacement6596)

    pattern6597 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons150)
    rule6597 = ReplacementRule(pattern6597, replacement6597)

    pattern6598 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons342)
    rule6598 = ReplacementRule(pattern6598, replacement6598)

    pattern6599 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons342)
    rule6599 = ReplacementRule(pattern6599, replacement6599)

    pattern6600 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule6600 = ReplacementRule(pattern6600, replacement6600)

    pattern6601 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule6601 = ReplacementRule(pattern6601, replacement6601)

    pattern6602 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons342)
    rule6602 = ReplacementRule(pattern6602, replacement6602)

    pattern6603 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**m_*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons342)
    rule6603 = ReplacementRule(pattern6603, replacement6603)

    pattern6604 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1764, cons1765)
    rule6604 = ReplacementRule(pattern6604, replacement6604)

    pattern6605 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons36, cons37, cons38, cons4, cons5, cons1764, cons1765)
    rule6605 = ReplacementRule(pattern6605, replacement6605)

    pattern6606 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1764, cons1765)
    rule6606 = ReplacementRule(pattern6606, replacement6606)

    pattern6607 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(c_ + x_*WC('d', S(1))))**WC('n', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons1764, cons1765)
    rule6607 = ReplacementRule(pattern6607, replacement6607)

    pattern6608 = Pattern(Integral(atanh(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons89)
    rule6608 = ReplacementRule(pattern6608, replacement6608)

    pattern6609 = Pattern(Integral(acoth(a_ + x_*WC('b', S(1)))/(c_ + x_**WC('n', S(1))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons89)
    rule6609 = ReplacementRule(pattern6609, replacement6609)

    pattern6610 = Pattern(Integral(atanh(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons1096)
    rule6610 = ReplacementRule(pattern6610, replacement6610)

    pattern6611 = Pattern(Integral(acoth(a_ + x_*WC('b', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons1096)
    rule6611 = ReplacementRule(pattern6611, replacement6611)

    pattern6612 = Pattern(Integral(atanh(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons4, cons1833)
    rule6612 = ReplacementRule(pattern6612, replacement6612)

    pattern6613 = Pattern(Integral(acoth(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons4, cons1833)
    rule6613 = ReplacementRule(pattern6613, replacement6613)

    pattern6614 = Pattern(Integral(atanh(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), cons2, cons3, cons4, cons1833)
    rule6614 = ReplacementRule(pattern6614, replacement6614)

    pattern6615 = Pattern(Integral(acoth(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))/x_, x_), cons2, cons3, cons4, cons1833)
    rule6615 = ReplacementRule(pattern6615, replacement6615)

    pattern6616 = Pattern(Integral(x_**WC('m', S(1))*atanh(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons95, cons1834, cons1835)
    rule6616 = ReplacementRule(pattern6616, replacement6616)

    pattern6617 = Pattern(Integral(x_**WC('m', S(1))*acoth(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons95, cons1834, cons1835)
    rule6617 = ReplacementRule(pattern6617, replacement6617)

    pattern6618 = Pattern(Integral(atanh(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons1836)
    rule6618 = ReplacementRule(pattern6618, replacement6618)

    pattern6619 = Pattern(Integral(acoth(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons1836)
    rule6619 = ReplacementRule(pattern6619, replacement6619)

    pattern6620 = Pattern(Integral(x_**WC('m', S(1))*atanh(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons20, cons170)
    rule6620 = ReplacementRule(pattern6620, replacement6620)

    pattern6621 = Pattern(Integral(x_**WC('m', S(1))*acoth(f_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons127, cons20, cons170)
    rule6621 = ReplacementRule(pattern6621, replacement6621)

    pattern6622 = Pattern(Integral(WC('u', S(1))*atanh(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6622 = ReplacementRule(pattern6622, replacement6622)

    pattern6623 = Pattern(Integral(WC('u', S(1))*acoth(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6623 = ReplacementRule(pattern6623, replacement6623)

    pattern6624 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons1941)
    rule6624 = ReplacementRule(pattern6624, replacement6624)

    pattern6625 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0)))*acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))), x_), cons2, cons3, cons8, cons1941)
    rule6625 = ReplacementRule(pattern6625, replacement6625)

    pattern6626 = Pattern(Integral(atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons19, cons1941, cons68)
    rule6626 = ReplacementRule(pattern6626, replacement6626)

    pattern6627 = Pattern(Integral(acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons19, cons1941, cons68)
    rule6627 = ReplacementRule(pattern6627, replacement6627)

    pattern6628 = Pattern(Integral(atanh(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1941, cons385)
    rule6628 = ReplacementRule(pattern6628, replacement6628)

    pattern6629 = Pattern(Integral(acoth(x_*WC('c', S(1))/sqrt(x_**S(2)*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1))/sqrt(x_**S(2)*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1941, cons385)
    rule6629 = ReplacementRule(pattern6629, replacement6629)

    pattern6630 = Pattern(Integral((x_**S(2)*WC('d', S(1)) + WC('c', S(0)))**n_*atanh(x_*WC('a', S(1))), x_), cons2, cons8, cons29, cons810, cons1588)
    rule6630 = ReplacementRule(pattern6630, With6630)

    pattern6631 = Pattern(Integral((x_**S(2)*WC('d', S(1)) + WC('c', S(0)))**n_*acoth(x_*WC('a', S(1))), x_), cons2, cons8, cons29, cons810, cons1588)
    rule6631 = ReplacementRule(pattern6631, With6631)

    pattern6632 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), cons820, cons87, cons465, cons1942, cons1943, CustomConstraint(With6632))
    rule6632 = ReplacementRule(pattern6632, replacement6632)

    pattern6633 = Pattern(Integral(u_*v_**WC('n', S(1)), x_), cons820, cons87, cons465, cons1942, cons1944, CustomConstraint(With6633))
    rule6633 = ReplacementRule(pattern6633, replacement6633)

    pattern6634 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1945)
    rule6634 = ReplacementRule(pattern6634, replacement6634)

    pattern6635 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1945)
    rule6635 = ReplacementRule(pattern6635, replacement6635)

    pattern6636 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1945)
    rule6636 = ReplacementRule(pattern6636, replacement6636)

    pattern6637 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1945)
    rule6637 = ReplacementRule(pattern6637, replacement6637)

    pattern6638 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1946)
    rule6638 = ReplacementRule(pattern6638, replacement6638)

    pattern6639 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1946)
    rule6639 = ReplacementRule(pattern6639, replacement6639)

    pattern6640 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1946)
    rule6640 = ReplacementRule(pattern6640, replacement6640)

    pattern6641 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1946)
    rule6641 = ReplacementRule(pattern6641, replacement6641)

    pattern6642 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1945)
    rule6642 = ReplacementRule(pattern6642, replacement6642)

    pattern6643 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1945)
    rule6643 = ReplacementRule(pattern6643, replacement6643)

    pattern6644 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1945)
    rule6644 = ReplacementRule(pattern6644, replacement6644)

    pattern6645 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1945)
    rule6645 = ReplacementRule(pattern6645, replacement6645)

    pattern6646 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1946)
    rule6646 = ReplacementRule(pattern6646, replacement6646)

    pattern6647 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1946)
    rule6647 = ReplacementRule(pattern6647, replacement6647)

    pattern6648 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1946)
    rule6648 = ReplacementRule(pattern6648, replacement6648)

    pattern6649 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))/tanh(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1946)
    rule6649 = ReplacementRule(pattern6649, replacement6649)

    pattern6650 = Pattern(Integral(atanh(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule6650 = ReplacementRule(pattern6650, replacement6650)

    pattern6651 = Pattern(Integral(acoth(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule6651 = ReplacementRule(pattern6651, replacement6651)

    pattern6652 = Pattern(Integral(atanh(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule6652 = ReplacementRule(pattern6652, replacement6652)

    pattern6653 = Pattern(Integral(acoth(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons69)
    rule6653 = ReplacementRule(pattern6653, replacement6653)

    pattern6654 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule6654 = ReplacementRule(pattern6654, replacement6654)

    pattern6655 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule6655 = ReplacementRule(pattern6655, replacement6655)

    pattern6656 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule6656 = ReplacementRule(pattern6656, replacement6656)

    pattern6657 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(S(1)/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons50, cons127, cons64)
    rule6657 = ReplacementRule(pattern6657, replacement6657)

    pattern6658 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1947)
    rule6658 = ReplacementRule(pattern6658, replacement6658)

    pattern6659 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1947)
    rule6659 = ReplacementRule(pattern6659, replacement6659)

    pattern6660 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1948)
    rule6660 = ReplacementRule(pattern6660, replacement6660)

    pattern6661 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1948)
    rule6661 = ReplacementRule(pattern6661, replacement6661)

    pattern6662 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1949)
    rule6662 = ReplacementRule(pattern6662, replacement6662)

    pattern6663 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1949)
    rule6663 = ReplacementRule(pattern6663, replacement6663)

    pattern6664 = Pattern(Integral(atanh(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1950)
    rule6664 = ReplacementRule(pattern6664, replacement6664)

    pattern6665 = Pattern(Integral(acoth(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons1950)
    rule6665 = ReplacementRule(pattern6665, replacement6665)

    pattern6666 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1947)
    rule6666 = ReplacementRule(pattern6666, replacement6666)

    pattern6667 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1947)
    rule6667 = ReplacementRule(pattern6667, replacement6667)

    pattern6668 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1948)
    rule6668 = ReplacementRule(pattern6668, replacement6668)

    pattern6669 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1948)
    rule6669 = ReplacementRule(pattern6669, replacement6669)

    pattern6670 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1949)
    rule6670 = ReplacementRule(pattern6670, replacement6670)

    pattern6671 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1949)
    rule6671 = ReplacementRule(pattern6671, replacement6671)

    pattern6672 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*atanh(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1950)
    rule6672 = ReplacementRule(pattern6672, replacement6672)

    pattern6673 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*acoth(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons64, cons1950)
    rule6673 = ReplacementRule(pattern6673, replacement6673)

    pattern6674 = Pattern(Integral(atanh(u_), x_), cons1232)
    rule6674 = ReplacementRule(pattern6674, replacement6674)

    pattern6675 = Pattern(Integral(acoth(u_), x_), cons1232)
    rule6675 = ReplacementRule(pattern6675, replacement6675)

    pattern6676 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*atanh(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1849)
    rule6676 = ReplacementRule(pattern6676, replacement6676)

    pattern6677 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acoth(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1849)
    rule6677 = ReplacementRule(pattern6677, replacement6677)

    pattern6678 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*atanh(u_)), x_), cons2, cons3, cons1232, cons1951, cons1952, CustomConstraint(With6678))
    rule6678 = ReplacementRule(pattern6678, replacement6678)

    pattern6679 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acoth(u_)), x_), cons2, cons3, cons1232, cons1953, cons1954, CustomConstraint(With6679))
    rule6679 = ReplacementRule(pattern6679, replacement6679)

    pattern6680 = Pattern(Integral(asech(x_*WC('c', S(1))), x_), cons8, cons8)
    rule6680 = ReplacementRule(pattern6680, replacement6680)

    pattern6681 = Pattern(Integral(acsch(x_*WC('c', S(1))), x_), cons8, cons8)
    rule6681 = ReplacementRule(pattern6681, replacement6681)

    pattern6682 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule6682 = ReplacementRule(pattern6682, replacement6682)

    pattern6683 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons1581)
    rule6683 = ReplacementRule(pattern6683, replacement6683)

    pattern6684 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons14)
    rule6684 = ReplacementRule(pattern6684, replacement6684)

    pattern6685 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))/x_, x_), cons2, cons3, cons8, cons14)
    rule6685 = ReplacementRule(pattern6685, replacement6685)

    pattern6686 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons68)
    rule6686 = ReplacementRule(pattern6686, replacement6686)

    pattern6687 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons19, cons68)
    rule6687 = ReplacementRule(pattern6687, replacement6687)

    pattern6688 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons20)
    rule6688 = ReplacementRule(pattern6688, replacement6688)

    pattern6689 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**n_, x_), cons2, cons3, cons8, cons4, cons20)
    rule6689 = ReplacementRule(pattern6689, replacement6689)

    pattern6690 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons1856)
    rule6690 = ReplacementRule(pattern6690, replacement6690)

    pattern6691 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons1856)
    rule6691 = ReplacementRule(pattern6691, replacement6691)

    pattern6692 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1745)
    rule6692 = ReplacementRule(pattern6692, With6692)

    pattern6693 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons1745)
    rule6693 = ReplacementRule(pattern6693, With6693)

    pattern6694 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons40)
    rule6694 = ReplacementRule(pattern6694, replacement6694)

    pattern6695 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons40)
    rule6695 = ReplacementRule(pattern6695, replacement6695)

    pattern6696 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons180, cons1857)
    rule6696 = ReplacementRule(pattern6696, replacement6696)

    pattern6697 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons669, cons180, cons1857)
    rule6697 = ReplacementRule(pattern6697, replacement6697)

    pattern6698 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons669, cons1858)
    rule6698 = ReplacementRule(pattern6698, replacement6698)

    pattern6699 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons669, cons1858)
    rule6699 = ReplacementRule(pattern6699, replacement6699)

    pattern6700 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule6700 = ReplacementRule(pattern6700, replacement6700)

    pattern6701 = Pattern(Integral((x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule6701 = ReplacementRule(pattern6701, replacement6701)

    pattern6702 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6702 = ReplacementRule(pattern6702, replacement6702)

    pattern6703 = Pattern(Integral(x_*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons5, cons56)
    rule6703 = ReplacementRule(pattern6703, replacement6703)

    pattern6704 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule6704 = ReplacementRule(pattern6704, With6704)

    pattern6705 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons1788)
    rule6705 = ReplacementRule(pattern6705, With6705)

    pattern6706 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1301)
    rule6706 = ReplacementRule(pattern6706, replacement6706)

    pattern6707 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1301)
    rule6707 = ReplacementRule(pattern6707, replacement6707)

    pattern6708 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons180, cons1857)
    rule6708 = ReplacementRule(pattern6708, replacement6708)

    pattern6709 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons20, cons669, cons180, cons1857)
    rule6709 = ReplacementRule(pattern6709, replacement6709)

    pattern6710 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1739, cons20, cons669, cons1858)
    rule6710 = ReplacementRule(pattern6710, replacement6710)

    pattern6711 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**p_*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1780, cons20, cons669, cons1858)
    rule6711 = ReplacementRule(pattern6711, replacement6711)

    pattern6712 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6712 = ReplacementRule(pattern6712, replacement6712)

    pattern6713 = Pattern(Integral(x_**WC('m', S(1))*(x_**S(2)*WC('e', S(1)) + WC('d', S(0)))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(x_*WC('c', S(1))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule6713 = ReplacementRule(pattern6713, replacement6713)

    pattern6714 = Pattern(Integral(asech(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule6714 = ReplacementRule(pattern6714, replacement6714)

    pattern6715 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule6715 = ReplacementRule(pattern6715, replacement6715)

    pattern6716 = Pattern(Integral(asech(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons1833)
    rule6716 = ReplacementRule(pattern6716, replacement6716)

    pattern6717 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons1833)
    rule6717 = ReplacementRule(pattern6717, replacement6717)

    pattern6718 = Pattern(Integral(asech(a_ + x_*WC('b', S(1)))/x_, x_), cons2, cons3, cons69)
    rule6718 = ReplacementRule(pattern6718, replacement6718)

    pattern6719 = Pattern(Integral(acsch(a_ + x_*WC('b', S(1)))/x_, x_), cons2, cons3, cons69)
    rule6719 = ReplacementRule(pattern6719, replacement6719)

    pattern6720 = Pattern(Integral(x_**WC('m', S(1))*asech(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons20, cons68)
    rule6720 = ReplacementRule(pattern6720, replacement6720)

    pattern6721 = Pattern(Integral(x_**WC('m', S(1))*acsch(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons20, cons68)
    rule6721 = ReplacementRule(pattern6721, replacement6721)

    pattern6722 = Pattern(Integral(x_**WC('m', S(1))*asech(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons64)
    rule6722 = ReplacementRule(pattern6722, replacement6722)

    pattern6723 = Pattern(Integral(x_**WC('m', S(1))*acsch(a_ + x_*WC('b', S(1)))**n_, x_), cons2, cons3, cons4, cons64)
    rule6723 = ReplacementRule(pattern6723, replacement6723)

    pattern6724 = Pattern(Integral(WC('u', S(1))*asech(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6724 = ReplacementRule(pattern6724, replacement6724)

    pattern6725 = Pattern(Integral(WC('u', S(1))*acsch(WC('c', S(1))/(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons4, cons19, cons1768)
    rule6725 = ReplacementRule(pattern6725, replacement6725)

    pattern6726 = Pattern(Integral(exp(asech(x_*WC('a', S(1)))), x_), cons2, cons2)
    rule6726 = ReplacementRule(pattern6726, replacement6726)

    pattern6727 = Pattern(Integral(exp(asech(x_**p_*WC('a', S(1)))), x_), cons2, cons5, cons1955)
    rule6727 = ReplacementRule(pattern6727, replacement6727)

    pattern6728 = Pattern(Integral(exp(acsch(x_**WC('p', S(1))*WC('a', S(1)))), x_), cons2, cons5, cons1955)
    rule6728 = ReplacementRule(pattern6728, replacement6728)

    pattern6729 = Pattern(Integral(exp(WC('n', S(1))*asech(u_)), x_), cons87)
    rule6729 = ReplacementRule(pattern6729, replacement6729)

    pattern6730 = Pattern(Integral(exp(WC('n', S(1))*acsch(u_)), x_), cons87)
    rule6730 = ReplacementRule(pattern6730, replacement6730)

    pattern6731 = Pattern(Integral(exp(asech(x_**WC('p', S(1))*WC('a', S(1))))/x_, x_), cons2, cons5, cons1955)
    rule6731 = ReplacementRule(pattern6731, replacement6731)

    pattern6732 = Pattern(Integral(x_**WC('m', S(1))*exp(asech(x_**WC('p', S(1))*WC('a', S(1)))), x_), cons2, cons19, cons5, cons68)
    rule6732 = ReplacementRule(pattern6732, replacement6732)

    pattern6733 = Pattern(Integral(x_**WC('m', S(1))*exp(acsch(x_**WC('p', S(1))*WC('a', S(1)))), x_), cons2, cons19, cons5, cons1956)
    rule6733 = ReplacementRule(pattern6733, replacement6733)

    pattern6734 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*asech(u_)), x_), cons19, cons87)
    rule6734 = ReplacementRule(pattern6734, replacement6734)

    pattern6735 = Pattern(Integral(x_**WC('m', S(1))*exp(WC('n', S(1))*acsch(u_)), x_), cons19, cons87)
    rule6735 = ReplacementRule(pattern6735, replacement6735)

    pattern6736 = Pattern(Integral(asech(u_), x_), cons1232, cons1771)
    rule6736 = ReplacementRule(pattern6736, replacement6736)

    pattern6737 = Pattern(Integral(acsch(u_), x_), cons1232, cons1771)
    rule6737 = ReplacementRule(pattern6737, replacement6737)

    pattern6738 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*asech(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule6738 = ReplacementRule(pattern6738, replacement6738)

    pattern6739 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*acsch(u_)), x_), cons2, cons3, cons8, cons29, cons19, cons68, cons1232, cons1772, cons1771)
    rule6739 = ReplacementRule(pattern6739, replacement6739)

    pattern6740 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*asech(u_)), x_), cons2, cons3, cons1232, cons1957, CustomConstraint(With6740))
    rule6740 = ReplacementRule(pattern6740, replacement6740)

    pattern6741 = Pattern(Integral(v_*(WC('a', S(0)) + WC('b', S(1))*acsch(u_)), x_), cons2, cons3, cons1232, cons1958, CustomConstraint(With6741))
    rule6741 = ReplacementRule(pattern6741, replacement6741)
    return [rule6087, rule6088, rule6089, rule6090, rule6091, rule6092, rule6093, rule6094, rule6095, rule6096, rule6097, rule6098, rule6099, rule6100, rule6101, rule6102, rule6103, rule6104, rule6105, rule6106, rule6107, rule6108, rule6109, rule6110, rule6111, rule6112, rule6113, rule6114, rule6115, rule6116, rule6117, rule6118, rule6119, rule6120, rule6121, rule6122, rule6123, rule6124, rule6125, rule6126, rule6127, rule6128, rule6129, rule6130, rule6131, rule6132, rule6133, rule6134, rule6135, rule6136, rule6137, rule6138, rule6139, rule6140, rule6141, rule6142, rule6143, rule6144, rule6145, rule6146, rule6147, rule6148, rule6149, rule6150, rule6151, rule6152, rule6153, rule6154, rule6155, rule6156, rule6157, rule6158, rule6159, rule6160, rule6161, rule6162, rule6163, rule6164, rule6165, rule6166, rule6167, rule6168, rule6169, rule6170, rule6171, rule6172, rule6173, rule6174, rule6175, rule6176, rule6177, rule6178, rule6179, rule6180, rule6181, rule6182, rule6183, rule6184, rule6185, rule6186, rule6187, rule6188, rule6189, rule6190, rule6191, rule6192, rule6193, rule6194, rule6195, rule6196, rule6197, rule6198, rule6199, rule6200, rule6201, rule6202, rule6203, rule6204, rule6205, rule6206, rule6207, rule6208, rule6209, rule6210, rule6211, rule6212, rule6213, rule6214, rule6215, rule6216, rule6217, rule6218, rule6219, rule6220, rule6221, rule6222, rule6223, rule6224, rule6225, rule6226, rule6227, rule6228, rule6229, rule6230, rule6231, rule6232, rule6233, rule6234, rule6235, rule6236, rule6237, rule6238, rule6239, rule6240, rule6241, rule6242, rule6243, rule6244, rule6245, rule6246, rule6247, rule6248, rule6249, rule6250, rule6251, rule6252, rule6253, rule6254, rule6255, rule6256, rule6257, rule6258, rule6259, rule6260, rule6261, rule6262, rule6263, rule6264, rule6265, rule6266, rule6267, rule6268, rule6269, rule6270, rule6271, rule6272, rule6273, rule6274, rule6275, rule6276, rule6277, rule6278, rule6279, rule6280, rule6281, rule6282, rule6283, rule6284, rule6285, rule6286, rule6287, rule6288, rule6289, rule6290, rule6291, rule6292, rule6293, rule6294, rule6295, rule6296, rule6297, rule6298, rule6299, rule6300, rule6301, rule6302, rule6303, rule6304, rule6305, rule6306, rule6307, rule6308, rule6309, rule6310, rule6311, rule6312, rule6313, rule6314, rule6315, rule6316, rule6317, rule6318, rule6319, rule6320, rule6321, rule6322, rule6323, rule6324, rule6325, rule6326, rule6327, rule6328, rule6329, rule6330, rule6331, rule6332, rule6333, rule6334, rule6335, rule6336, rule6337, rule6338, rule6339, rule6340, rule6341, rule6342, rule6343, rule6344, rule6345, rule6346, rule6347, rule6348, rule6349, rule6350, rule6351, rule6352, rule6353, rule6354, rule6355, rule6356, rule6357, rule6358, rule6359, rule6360, rule6361, rule6362, rule6363, rule6364, rule6365, rule6366, rule6367, rule6368, rule6369, rule6370, rule6371, rule6372, rule6373, rule6374, rule6375, rule6376, rule6377, rule6378, rule6379, rule6380, rule6381, rule6382, rule6383, rule6384, rule6385, rule6386, rule6387, rule6388, rule6389, rule6390, rule6391, rule6392, rule6393, rule6394, rule6395, rule6396, rule6397, rule6398, rule6399, rule6400, rule6401, rule6402, rule6403, rule6404, rule6405, rule6406, rule6407, rule6408, rule6409, rule6410, rule6411, rule6412, rule6413, rule6414, rule6415, rule6416, rule6417, rule6418, rule6419, rule6420, rule6421, rule6422, rule6423, rule6424, rule6425, rule6426, rule6427, rule6428, rule6429, rule6430, rule6431, rule6432, rule6433, rule6434, rule6435, rule6436, rule6437, rule6438, rule6439, rule6440, rule6441, rule6442, rule6443, rule6444, rule6445, rule6446, rule6447, rule6448, rule6449, rule6450, rule6451, rule6452, rule6453, rule6454, rule6455, rule6456, rule6457, rule6458, rule6459, rule6460, rule6461, rule6462, rule6463, rule6464, rule6465, rule6466, rule6467, rule6468, rule6469, rule6470, rule6471, rule6472, rule6473, rule6474, rule6475, rule6476, rule6477, rule6478, rule6479, rule6480, rule6481, rule6482, rule6483, rule6484, rule6485, rule6486, rule6487, rule6488, rule6489, rule6490, rule6491, rule6492, rule6493, rule6494, rule6495, rule6496, rule6497, rule6498, rule6499, rule6500, rule6501, rule6502, rule6503, rule6504, rule6505, rule6506, rule6507, rule6508, rule6509, rule6510, rule6511, rule6512, rule6513, rule6514, rule6515, rule6516, rule6517, rule6518, rule6519, rule6520, rule6521, rule6522, rule6523, rule6524, rule6525, rule6526, rule6527, rule6528, rule6529, rule6530, rule6531, rule6532, rule6533, rule6534, rule6535, rule6536, rule6537, rule6538, rule6539, rule6540, rule6541, rule6542, rule6543, rule6544, rule6545, rule6546, rule6547, rule6548, rule6549, rule6550, rule6551, rule6552, rule6553, rule6554, rule6555, rule6556, rule6557, rule6558, rule6559, rule6560, rule6561, rule6562, rule6563, rule6564, rule6565, rule6566, rule6567, rule6568, rule6569, rule6570, rule6571, rule6572, rule6573, rule6574, rule6575, rule6576, rule6577, rule6578, rule6579, rule6580, rule6581, rule6582, rule6583, rule6584, rule6585, rule6586, rule6587, rule6588, rule6589, rule6590, rule6591, rule6592, rule6593, rule6594, rule6595, rule6596, rule6597, rule6598, rule6599, rule6600, rule6601, rule6602, rule6603, rule6604, rule6605, rule6606, rule6607, rule6608, rule6609, rule6610, rule6611, rule6612, rule6613, rule6614, rule6615, rule6616, rule6617, rule6618, rule6619, rule6620, rule6621, rule6622, rule6623, rule6624, rule6625, rule6626, rule6627, rule6628, rule6629, rule6630, rule6631, rule6632, rule6633, rule6634, rule6635, rule6636, rule6637, rule6638, rule6639, rule6640, rule6641, rule6642, rule6643, rule6644, rule6645, rule6646, rule6647, rule6648, rule6649, rule6650, rule6651, rule6652, rule6653, rule6654, rule6655, rule6656, rule6657, rule6658, rule6659, rule6660, rule6661, rule6662, rule6663, rule6664, rule6665, rule6666, rule6667, rule6668, rule6669, rule6670, rule6671, rule6672, rule6673, rule6674, rule6675, rule6676, rule6677, rule6678, rule6679, rule6680, rule6681, rule6682, rule6683, rule6684, rule6685, rule6686, rule6687, rule6688, rule6689, rule6690, rule6691, rule6692, rule6693, rule6694, rule6695, rule6696, rule6697, rule6698, rule6699, rule6700, rule6701, rule6702, rule6703, rule6704, rule6705, rule6706, rule6707, rule6708, rule6709, rule6710, rule6711, rule6712, rule6713, rule6714, rule6715, rule6716, rule6717, rule6718, rule6719, rule6720, rule6721, rule6722, rule6723, rule6724, rule6725, rule6726, rule6727, rule6728, rule6729, rule6730, rule6731, rule6732, rule6733, rule6734, rule6735, rule6736, rule6737, rule6738, rule6739, rule6740, rule6741, ]





def replacement6087(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*asinh(c*x))**n, x)


def replacement6088(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp(x*(a + b*acosh(c*x))**n, x)


def replacement6089(a, b, c, n, x):
    return -Dist(c/(b*(n + S(1))), Int(x*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6090(a, b, c, n, x):
    return -Dist(c/(b*(n + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6091(a, b, c, n, x):
    return Dist(S(1)/(b*c), Subst(Int(x**n*cosh(a/b - x/b), x), x, a + b*asinh(c*x)), x)


def replacement6092(a, b, c, n, x):
    return -Dist(S(1)/(b*c), Subst(Int(x**n*sinh(a/b - x/b), x), x, a + b*acosh(c*x)), x)


def replacement6093(a, b, c, n, x):
    return Subst(Int((a + b*x)**n/tanh(x), x), x, asinh(c*x))


def replacement6094(a, b, c, n, x):
    return Subst(Int((a + b*x)**n*tanh(x), x), x, acosh(c*x))


def replacement6095(a, b, c, d, m, n, x):
    return -Dist(b*c*n/(d*(m + S(1))), Int((d*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp((d*x)**(m + S(1))*(a + b*asinh(c*x))**n/(d*(m + S(1))), x)


def replacement6096(a, b, c, d, m, n, x):
    return -Dist(b*c*n/(d*(m + S(1))), Int((d*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp((d*x)**(m + S(1))*(a + b*acosh(c*x))**n/(d*(m + S(1))), x)


def replacement6097(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*asinh(c*x))**n/(m + S(1)), x)


def replacement6098(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp(x**(m + S(1))*(a + b*acosh(c*x))**n/(m + S(1)), x)


def replacement6099(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1))/(b*(n + S(1))), Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1)), (m + (m + S(1))*sinh(x)**S(2))*sinh(x)**(m + S(-1)), x), x), x, asinh(c*x)), x) + Simp(x**m*(a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6100(a, b, c, m, n, x):
    return Dist(c**(-m + S(-1))/(b*(n + S(1))), Subst(Int(ExpandTrigReduce((a + b*x)**(n + S(1))*(m - (m + S(1))*cosh(x)**S(2))*cosh(x)**(m + S(-1)), x), x), x, acosh(c*x)), x) + Simp(x**m*(a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6101(a, b, c, m, n, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) - Dist(c*(m + S(1))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*asinh(c*x))**(n + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**m*(a + b*asinh(c*x))**(n + S(1))*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6102(a, b, c, m, n, x):
    return Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) - Dist(c*(m + S(1))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp(x**m*(a + b*acosh(c*x))**(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6103(a, b, c, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*sinh(x)**m*cosh(x), x), x, asinh(c*x)), x)


def replacement6104(a, b, c, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*sinh(x)*cosh(x)**m, x), x, acosh(c*x)), x)


def replacement6105(a, b, c, d, m, n, x):
    return Int((d*x)**m*(a + b*asinh(c*x))**n, x)


def replacement6106(a, b, c, d, m, n, x):
    return Int((d*x)**m*(a + b*acosh(c*x))**n, x)


def replacement6107(a, b, c, d, e, x):
    return Simp(log(a + b*asinh(c*x))/(b*c*sqrt(d)), x)


def replacement6108(a, b, c, d1, d2, e1, e2, x):
    return Simp(log(a + b*acosh(c*x))/(b*c*sqrt(-d1*d2)), x)


def replacement6109(a, b, c, d, e, n, x):
    return Simp((a + b*asinh(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6110(a, b, c, d1, d2, e1, e2, n, x):
    return Simp((a + b*acosh(c*x))**(n + S(1))/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6111(a, b, c, d, e, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x)


def replacement6112(a, b, c, d1, d2, e1, e2, n, x):
    return Dist(sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), Int((a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x)


def With6113(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6114(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6115(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*n*(-d)**p/(S(2)*p + S(1)), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x)


def replacement6116(a, b, c, d, e, n, x):
    return Dist(sqrt(d + e*x**S(2))/(S(2)*sqrt(c**S(2)*x**S(2) + S(1))), Int((a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(S(2)*sqrt(c**S(2)*x**S(2) + S(1))), Int(x*(a + b*asinh(c*x))**(n + S(-1)), x), x) + Simp(x*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/S(2), x)


def replacement6117(a, b, c, d1, d2, e1, e2, n, x):
    return -Dist(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(S(2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) - Dist(b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(S(2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(-1)), x), x) + Simp(x*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/S(2), x)


def replacement6118(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*p + S(1)), Int(x*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x)


def replacement6119(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist(S(2)*d1*d2*p/(S(2)*p + S(1)), Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x), x) - Dist(b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/((S(2)*p + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(S(2)*p + S(1)), x)


def replacement6120(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist(S(2)*d1*d2*p/(S(2)*p + S(1)), Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x), x) - Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*p + S(1)), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp(x*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(S(2)*p + S(1)), x)


def replacement6121(a, b, c, d, e, n, x):
    return -Dist(b*c*n*sqrt(c**S(2)*x**S(2) + S(1))/(d*sqrt(d + e*x**S(2))), Int(x*(a + b*asinh(c*x))**(n + S(-1))/(c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*asinh(c*x))**n/(d*sqrt(d + e*x**S(2))), x)


def replacement6122(a, b, c, d1, d2, e1, e2, n, x):
    return Dist(b*c*n*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(d1*d2*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), Int(x*(a + b*acosh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*acosh(c*x))**n/(d1*d2*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), x)


def replacement6123(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*c*n*(-d)**p/(S(2)*p + S(2)), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) - Simp(x*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement6124(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*(p + S(1))), Int(x*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp(x*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement6125(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d1*d2*(p + S(1))), Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*c*n*(-d1*d2)**(p + S(1)/2)*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(S(2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(p + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) - Simp(x*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*(p + S(1))), x)


def replacement6126(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d1*d2*(p + S(1))), Int((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*(p + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) - Simp(x*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*(p + S(1))), x)


def replacement6127(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*d), Subst(Int((a + b*x)**n/cosh(x), x), x, asinh(c*x)), x)


def replacement6128(a, b, c, d, e, n, x):
    return -Dist(S(1)/(c*d), Subst(Int((a + b*x)**n/sinh(x), x), x, acosh(c*x)), x)


def replacement6129(a, b, c, d, e, n, p, x):
    return -Dist(c*(-d)**p*(S(2)*p + S(1))/(b*(n + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((-d)**p*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2)/(b*c*(n + S(1))), x)


def replacement6130(a, b, c, d, e, n, p, x):
    return -Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(S(2)*p + S(1))*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*(n + S(1))), Int(x*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6131(a, b, c, d1, d2, e1, e2, n, p, x):
    return -Dist(c*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(S(2)*p + S(1))/(b*(n + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6132(a, b, c, d1, d2, e1, e2, n, p, x):
    return -Dist(c*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(S(2)*p + S(1))*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(b*(n + S(1))), Int(x*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6133(a, b, c, d, e, n, p, x):
    return Dist(d**p/c, Subst(Int((a + b*x)**n*cosh(x)**(S(2)*p + S(1)), x), x, asinh(c*x)), x)


def replacement6134(a, b, c, d, e, n, p, x):
    return Dist((-d)**p/c, Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x)), x)


def replacement6135(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist((-d1*d2)**p/c, Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1)), x), x, acosh(c*x)), x)


def replacement6136(a, b, c, d, e, n, p, x):
    return Dist(d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(c**S(2)*x**S(2) + S(1)), Int((a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement6137(a, b, c, d1, d2, e1, e2, n, p, x):
    return Dist((-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def With6138(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6139(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6140(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n, (d + e*x**S(2))**p, x), x)


def replacement6141(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n, (d + e*x**S(2))**p, x), x)


def replacement6142(a, b, c, d, e, n, p, x):
    return Int((a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6143(a, b, c, d, e, n, p, x):
    return Int((a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6144(a, b, c, d1, d2, e1, e2, n, p, x):
    return Int((a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)


def replacement6145(a, b, c, d, e, f, g, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((a + b*asinh(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement6146(a, b, c, d, e, n, p, x):
    return Dist((-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def replacement6147(a, b, c, d, e, n, x):
    return Dist(S(1)/e, Subst(Int((a + b*x)**n*tanh(x), x), x, asinh(c*x)), x)


def replacement6148(a, b, c, d, e, n, x):
    return Dist(S(1)/e, Subst(Int((a + b*x)**n/tanh(x), x), x, acosh(c*x)), x)


def replacement6149(a, b, c, d, e, n, p, x):
    return -Dist(b*n*(-d)**p/(S(2)*c*(p + S(1))), Int((a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6150(a, b, c, d, e, n, p, x):
    return -Dist(b*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6151(a, b, c, d1, d2, e1, e2, n, p, x):
    return -Dist(b*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) + Simp((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))), x)


def replacement6152(a, b, c, d1, d2, e1, e2, n, p, x):
    return -Dist(b*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))), x)


def replacement6153(a, b, c, d, e, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*x)**n/(sinh(x)*cosh(x)), x), x, asinh(c*x)), x)


def replacement6154(a, b, c, d, e, n, x):
    return -Dist(S(1)/d, Subst(Int((a + b*x)**n/(sinh(x)*cosh(x)), x), x, acosh(c*x)), x)


def replacement6155(a, b, c, d, e, f, m, n, p, x):
    return Dist(b*c*n*(-d)**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement6156(a, b, c, d, e, f, m, n, p, x):
    return -Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement6157(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))), x)


def replacement6158(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))), x)


def replacement6159(a, b, c, d, e, p, x):
    return Dist(d, Int((a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x), x) - Dist(b*c*d**p/(S(2)*p), Int((c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*asinh(c*x))*(d + e*x**S(2))**p/(S(2)*p), x)


def replacement6160(a, b, c, d, e, p, x):
    return Dist(d, Int((a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(-1))/x, x), x) - Dist(b*c*(-d)**p/(S(2)*p), Int((c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((a + b*acosh(c*x))*(d + e*x**S(2))**p/(S(2)*p), x)


def replacement6161(a, b, c, d, e, f, m, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement6162(a, b, c, d, e, f, m, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*(-d)**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def With6163(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6164(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def With6165(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(c**S(2)*x**S(2) + S(1))**p, x)
    return Dist(d**p*(a + b*asinh(c*x)), u, x) - Dist(b*c*d**p, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x)


def With6166(a, b, c, d1, d2, e1, e2, m, p, x):
    u = IntHide(x**m*(c*x + S(-1))**p*(c*x + S(1))**p, x)
    return Dist((-d1*d2)**p*(a + b*acosh(c*x)), u, x) - Dist(b*c*(-d1*d2)**p, Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x)


def With6167(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(c**S(2)*x**S(2) + S(1))**p, x)
    return -Dist(b*c*d**(p + S(-1)/2)*sqrt(d + e*x**S(2))/sqrt(c**S(2)*x**S(2) + S(1)), Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), Int(x**m*(d + e*x**S(2))**p, x), x)


def With6168(a, b, c, d1, d2, e1, e2, m, p, x):
    u = IntHide(x**m*(c*x + S(-1))**p*(c*x + S(1))**p, x)
    return -Dist(b*c*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*acosh(c*x), Int(x**m*(d1 + e1*x)**p*(d2 + e2*x)**p, x), x)


def replacement6169(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*n*(-d)**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement6170(a, b, c, d, e, f, m, n, x):
    return -Dist(c**S(2)*sqrt(d + e*x**S(2))/(f**S(2)*(m + S(1))*sqrt(c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(1))*sqrt(c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(1))), x)


def replacement6171(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return -Dist(c**S(2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f**S(2)*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) - Dist(b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(1))), x)


def replacement6172(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(2)*e*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(1))), x)


def replacement6173(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return -Dist(S(2)*e1*e2*p/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x), x) - Dist(b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(1))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(f*(m + S(1))), x)


def replacement6174(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(2)*d*p/(m + S(2)*p + S(1)), Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*n*(-d)**p/(f*(m + S(2)*p + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))), x)


def replacement6175(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(d + e*x**S(2))/((m + S(2))*sqrt(c**S(2)*x**S(2) + S(1))), Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x) - Dist(b*c*n*sqrt(d + e*x**S(2))/(f*(m + S(2))*sqrt(c**S(2)*x**S(2) + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(f*(m + S(2))), x)


def replacement6176(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return -Dist(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/((m + S(2))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) - Dist(b*c*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(2))*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1)), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*(m + S(2))), x)


def replacement6177(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(2)*d*p/(m + S(2)*p + S(1)), Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(2)*p + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p/(f*(m + S(2)*p + S(1))), x)


def replacement6178(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(S(2)*d1*d2*p/(m + S(2)*p + S(1)), Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(-1))*(d2 + e2*x)**(p + S(-1)), x), x) - Dist(b*c*n*(-d1*d2)**(p + S(-1)/2)*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(f*sqrt(c*x + S(-1))*sqrt(c*x + S(1))*(m + S(2)*p + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p/(f*(m + S(2)*p + S(1))), x)


def replacement6179(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x), x) + Dist(b*c*n*(-d)**p/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement6180(a, b, c, d, e, f, m, n, p, x):
    return -Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x), x) - Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*f*(m + S(1))), x)


def replacement6181(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x), x) + Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))), x)


def replacement6182(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(c**S(2)*(m + S(2)*p + S(3))/(f**S(2)*(m + S(1))), Int((f*x)**(m + S(2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x), x) + Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(f*(m + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(d1*d2*f*(m + S(1))), x)


def replacement6183(a, b, c, d, e, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*f*n*(-d)**p/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6184(a, b, c, d, e, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6185(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e1*e2*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))), x)


def replacement6186(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(S(2)*e1*e2*(p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*c*(p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*e1*e2*(p + S(1))), x)


def replacement6187(a, b, c, d, e, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(b*c*n*(-d)**p/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))), x)


def replacement6188(a, b, c, d, e, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b*c*d**IntPart(p)*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*f*(p + S(1))), x)


def replacement6189(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d1*d2*(p + S(1))), Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*f*(p + S(1))), x)


def replacement6190(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist((m + S(2)*p + S(3))/(S(2)*d1*d2*(p + S(1))), Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1)), x), x) - Dist(b*c*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(S(2)*f*(p + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) - Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(S(2)*d1*d2*f*(p + S(1))), x)


def replacement6191(a, b, c, d, e, f, m, n, x):
    return -Dist(f**S(2)*(m + S(-1))/(c**S(2)*m), Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*f*n*sqrt(c**S(2)*x**S(2) + S(1))/(c*m*sqrt(d + e*x**S(2))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1)), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*sqrt(d + e*x**S(2))/(e*m), x)


def replacement6192(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*m), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), x), x) + Dist(b*f*n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(c*d1*d2*m*sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1)), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)/(e1*e2*m), x)


def replacement6193(a, b, c, d, e, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*sinh(x)**m, x), x, asinh(c*x)), x)


def replacement6194(a, b, c, d1, d2, e1, e2, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(-d1*d2), Subst(Int((a + b*x)**n*cosh(x)**m, x), x, acosh(c*x)), x)


def replacement6195(a, b, c, d, e, f, m, x):
    return Simp((f*x)**(m + S(1))*(a + b*asinh(c*x))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, -c**S(2)*x**S(2))/(sqrt(d)*f*(m + S(1))), x) - Simp(b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), -c**S(2)*x**S(2))/(sqrt(d)*f**S(2)*(m + S(1))*(m + S(2))), x)


def replacement6196(a, b, c, d1, d2, e1, e2, f, m, x):
    return Simp(b*c*(f*x)**(m + S(2))*HypergeometricPFQ(List(S(1), m/S(2) + S(1), m/S(2) + S(1)), List(m/S(2) + S(3)/2, m/S(2) + S(2)), c**S(2)*x**S(2))/(f**S(2)*sqrt(-d1*d2)*(m + S(1))*(m + S(2))), x) + Simp((f*x)**(m + S(1))*(a + b*acosh(c*x))*sqrt(-c**S(2)*x**S(2) + S(1))*Hypergeometric2F1(S(1)/2, m/S(2) + S(1)/2, m/S(2) + S(3)/2, c**S(2)*x**S(2))/(f*sqrt(d1 + e1*x)*sqrt(d2 + e2*x)*(m + S(1))), x)


def replacement6197(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x)


def replacement6198(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return Dist(sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x)


def replacement6199(a, b, c, d, e, f, m, n, p, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x), x) - Dist(b*f*n*(-d)**p/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement6200(a, b, c, d, e, f, m, n, p, x):
    return -Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x), x) - Dist(b*d**IntPart(p)*f*n*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*asinh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(e*(m + S(2)*p + S(1))), x)


def replacement6201(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x), x) - Dist(b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c**S(2)*x**S(2) + S(-1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(e1*e2*(m + S(2)*p + S(1))), x)


def replacement6202(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(f**S(2)*(m + S(-1))/(c**S(2)*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-2))*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x), x) - Dist(b*f*n*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(c*(m + S(2)*p + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(-1))*(c*x + S(-1))**(p + S(1)/2)*(c*x + S(1))**(p + S(1)/2), x), x) + Simp(f*(f*x)**(m + S(-1))*(a + b*acosh(c*x))**n*(d1 + e1*x)**(p + S(1))*(d2 + e2*x)**(p + S(1))/(e1*e2*(m + S(2)*p + S(1))), x)


def replacement6203(a, b, c, d, e, f, m, n, p, x):
    return Dist(f*m*(-d)**p/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6204(a, b, c, d, e, f, m, n, p, x):
    return -Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6205(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6206(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6207(a, b, c, d, e, f, m, n, x):
    return -Dist(f*m/(b*c*sqrt(d)*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1)), x), x) + Simp((f*x)**m*(a + b*asinh(c*x))**(n + S(1))/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6208(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return -Dist(f*m/(b*c*sqrt(-d1*d2)*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1)), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6209(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((f*x)**m*(a + b*asinh(c*x))**n/sqrt(c**S(2)*x**S(2) + S(1)), x), x)


def replacement6210(a, b, c, d1, d2, e1, e2, f, m, n, x):
    return Dist(sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), Int((f*x)**m*(a + b*acosh(c*x))**n/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x)


def replacement6211(a, b, c, d, e, f, m, n, p, x):
    return Dist(f*m*(-d)**p/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) - Dist(c*(-d)**p*(m + S(2)*p + S(1))/(b*f*(n + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))*(c*x + S(-1))**(p + S(-1)/2)*(c*x + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6212(a, b, c, d, e, f, m, n, p, x):
    return -Dist(d**IntPart(p)*f*m*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) - Dist(c*d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))/(b*f*(n + S(1))), Int((f*x)**(m + S(1))*(a + b*asinh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**p*sqrt(c**S(2)*x**S(2) + S(1))/(b*c*(n + S(1))), x)


def replacement6213(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Dist(f*m*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))/(b*c*(n + S(1))), Int((f*x)**(m + S(-1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) - Dist(c*(-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p))*(m + S(2)*p + S(1))/(b*f*(n + S(1))), Int((f*x)**(m + S(1))*(a + b*acosh(c*x))**(n + S(1))*(c**S(2)*x**S(2) + S(-1))**(p + S(-1)/2), x), x) + Simp((f*x)**m*(a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**p*(d2 + e2*x)**p*sqrt(c*x + S(-1))*sqrt(c*x + S(1))/(b*c*(n + S(1))), x)


def replacement6214(a, b, c, d, e, m, n, p, x):
    return Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sinh(x)**m*cosh(x)**(S(2)*p + S(1)), x), x, asinh(c*x)), x)


def replacement6215(a, b, c, d, e, m, n, p, x):
    return Dist(c**(-m + S(-1))*(-d)**p, Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1))*cosh(x)**m, x), x, acosh(c*x)), x)


def replacement6216(a, b, c, d1, d2, e1, e2, m, n, p, x):
    return Dist(c**(-m + S(-1))*(-d1*d2)**p, Subst(Int((a + b*x)**n*sinh(x)**(S(2)*p + S(1))*cosh(x)**m, x), x, acosh(c*x)), x)


def replacement6217(a, b, c, d, e, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(x**m*(a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement6218(a, b, c, d1, d2, e1, e2, m, n, p, x):
    return Dist((-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int(x**m*(a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def replacement6219(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), (f*x)**m*(d + e*x**S(2))**(p + S(1)/2), x), x)


def replacement6220(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), (f*x)**m*(d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2), x), x)


def replacement6221(a, b, c, d, e, f, m, x):
    return -Dist(b*c/(f*(m + S(1))*(m + S(3))), Int((f*x)**(m + S(1))*(d*(m + S(3)) + e*x**S(2)*(m + S(1)))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp(d*(f*x)**(m + S(1))*(a + b*acosh(c*x))/(f*(m + S(1))), x) + Simp(e*(f*x)**(m + S(3))*(a + b*acosh(c*x))/(f**S(3)*(m + S(3))), x)


def replacement6222(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asinh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6223(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp((a + b*acosh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With6224(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6225(a, b, c, d, e, f, m, p, x):
    u = IntHide((f*x)**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6226(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x)


def replacement6227(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n, (f*x)**m*(d + e*x**S(2))**p, x), x)


def replacement6228(a, b, c, d, e, f, m, n, p, x):
    return Int((f*x)**m*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6229(a, b, c, d, e, f, m, n, p, x):
    return Int((f*x)**m*(a + b*acosh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6230(a, b, c, d1, d2, e1, e2, f, m, n, p, x):
    return Int((f*x)**m*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)


def replacement6231(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist((d + e*x)**FracPart(p)*(f + g*x)**FracPart(p)*(d*f + e*g*x**S(2))**(-FracPart(p)), Int((h*x)**m*(a + b*asinh(c*x))**n*(d*f + e*g*x**S(2))**p, x), x)


def replacement6232(a, b, c, d, e, f, m, n, p, x):
    return Dist((-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int((f*x)**m*(a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def replacement6233(a, b, c, d, e, n, x):
    return Subst(Int((a + b*x)**n*cosh(x)/(c*d + e*sinh(x)), x), x, asinh(c*x))


def replacement6234(a, b, c, d, e, n, x):
    return Subst(Int((a + b*x)**n*sinh(x)/(c*d + e*cosh(x)), x), x, acosh(c*x))


def replacement6235(a, b, c, d, e, m, n, x):
    return -Dist(b*c*n/(e*(m + S(1))), Int((a + b*asinh(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*asinh(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement6236(a, b, c, d, e, m, n, x):
    return -Dist(b*c*n/(e*(m + S(1))), Int((a + b*acosh(c*x))**(n + S(-1))*(d + e*x)**(m + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x) + Simp((a + b*acosh(c*x))**n*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement6237(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n*(d + e*x)**m, x), x)


def replacement6238(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n*(d + e*x)**m, x), x)


def replacement6239(a, b, c, d, e, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(c*d + e*sinh(x))**m*cosh(x), x), x, asinh(c*x)), x)


def replacement6240(a, b, c, d, e, m, n, x):
    return Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(c*d + e*cosh(x))**m*sinh(x), x), x, acosh(c*x)), x)


def With6241(Px, a, b, c, x):
    u = IntHide(Px, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6242(Px, a, b, c, x):
    u = IntHide(Px, x)
    return -Dist(b*c*sqrt(-c**S(2)*x**S(2) + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6243(Px, a, b, c, n, x):
    return Int(ExpandIntegrand(Px*(a + b*asinh(c*x))**n, x), x)


def replacement6244(Px, a, b, c, n, x):
    return Int(ExpandIntegrand(Px*(a + b*acosh(c*x))**n, x), x)


def With6245(Px, a, b, c, d, e, m, x):
    u = IntHide(Px*(d + e*x)**m, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6246(Px, a, b, c, d, e, m, x):
    u = IntHide(Px*(d + e*x)**m, x)
    return -Dist(b*c*sqrt(-c**S(2)*x**S(2) + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(SimplifyIntegrand(u/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acosh(c*x), u, x)


def With6247(a, b, c, d, e, f, g, m, n, p, x):
    u = IntHide((d + e*x)**m*(f + g*x)**p, x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*asinh(c*x))**n, u, x)


def With6248(a, b, c, d, e, f, g, m, n, p, x):
    u = IntHide((d + e*x)**m*(f + g*x)**p, x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist((a + b*acosh(c*x))**n, u, x)


def With6249(a, b, c, d, e, f, g, h, n, p, x):
    u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*asinh(c*x))**(n + S(-1))/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist((a + b*asinh(c*x))**n, u, x)


def With6250(a, b, c, d, e, f, g, h, n, p, x):
    u = IntHide((f + g*x + h*x**S(2))**p/(d + e*x)**S(2), x)
    return -Dist(b*c*n, Int(SimplifyIntegrand(u*(a + b*acosh(c*x))**(n + S(-1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), x), x), x) + Dist((a + b*acosh(c*x))**n, u, x)


def replacement6251(Px, a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(d + e*x)**m, x), x)


def replacement6252(Px, a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(d + e*x)**m, x), x)


def With6253(a, b, c, d, e, f, g, m, p, x):
    u = IntHide((d + e*x**S(2))**p*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/sqrt(c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6254(a, b, c, d1, d2, e1, e2, f, g, m, p, x):
    u = IntHide((d1 + e1*x)**p*(d2 + e2*x)**p*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), u, x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6255(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n*(d + e*x**S(2))**p, (f + g*x)**m, x), x)


def replacement6256(a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, (f + g*x)**m, x), x)


def replacement6257(a, b, c, d, e, f, g, m, n, x):
    return -Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d*g*m + S(2)*e*f*x + e*g*x**S(2)*(m + S(2))), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6258(a, b, c, d1, d2, e1, e2, f, g, m, n, x):
    return -Dist(S(1)/(b*c*sqrt(-d1*d2)*(n + S(1))), Int((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1))*(d1*d2*g*m + S(2)*e1*e2*f*x + e1*e2*g*x**S(2)*(m + S(2))), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**m*(d1*d2 + e1*e2*x**S(2))/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6259(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n*sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(-1)/2)*(f + g*x)**m, x), x)


def replacement6260(a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n*sqrt(d1 + e1*x)*sqrt(d2 + e2*x), (d1 + e1*x)**(p + S(-1)/2)*(d2 + e2*x)**(p + S(-1)/2)*(f + g*x)**m, x), x)


def replacement6261(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(S(1)/(b*c*sqrt(d)*(n + S(1))), Int(ExpandIntegrand((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d + e*x**S(2))**(p + S(-1)/2)*(d*g*m + e*f*x*(S(2)*p + S(1)) + e*g*x**S(2)*(m + S(2)*p + S(1))), x), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6262(a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    return -Dist(S(1)/(b*c*sqrt(-d1*d2)*(n + S(1))), Int(ExpandIntegrand((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), (d1 + e1*x)**(p + S(-1)/2)*(d2 + e2*x)**(p + S(-1)/2)*(d1*d2*g*m + e1*e2*f*x*(S(2)*p + S(1)) + e1*e2*g*x**S(2)*(m + S(2)*p + S(1))), x), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*(d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2)*(f + g*x)**m/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6263(a, b, c, d, e, f, g, m, n, x):
    return -Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6264(a, b, c, d1, d2, e1, e2, f, g, m, n, x):
    return -Dist(g*m/(b*c*sqrt(-d1*d2)*(n + S(1))), Int((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**(m + S(-1)), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*(f + g*x)**m/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6265(a, b, c, d, e, f, g, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(d), Subst(Int((a + b*x)**n*(c*f + g*sinh(x))**m, x), x, asinh(c*x)), x)


def replacement6266(a, b, c, d1, d2, e1, e2, f, g, m, n, x):
    return Dist(c**(-m + S(-1))/sqrt(-d1*d2), Subst(Int((a + b*x)**n*(c*f + g*cosh(x))**m, x), x, acosh(c*x)), x)


def replacement6267(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n/sqrt(d + e*x**S(2)), (d + e*x**S(2))**(p + S(1)/2)*(f + g*x)**m, x), x)


def replacement6268(a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n/(sqrt(d1 + e1*x)*sqrt(d2 + e2*x)), (d1 + e1*x)**(p + S(1)/2)*(d2 + e2*x)**(p + S(1)/2)*(f + g*x)**m, x), x)


def replacement6269(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*asinh(c*x))**n*(f + g*x)**m*(c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement6270(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int((a + b*acosh(c*x))**n*(f + g*x)**m*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def replacement6271(a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    return Dist((-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(-c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*acosh(c*x))**n*(f + g*x)**m*(c*x + S(-1))**p*(c*x + S(1))**p, x), x)


def replacement6272(a, b, c, d, e, f, g, h, m, n, x):
    return -Dist(g*m/(b*c*sqrt(d)*(n + S(1))), Int((a + b*asinh(c*x))**(n + S(1))/(f + g*x), x), x) + Simp((a + b*asinh(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(d)*(n + S(1))), x)


def replacement6273(a, b, c, d1, d2, e1, e2, f, g, h, m, n, x):
    return -Dist(g*m/(b*c*sqrt(-d1*d2)*(n + S(1))), Int((a + b*acosh(c*x))**(n + S(1))/(f + g*x), x), x) + Simp((a + b*acosh(c*x))**(n + S(1))*log(h*(f + g*x)**m)/(b*c*sqrt(-d1*d2)*(n + S(1))), x)


def replacement6274(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(d**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((a + b*asinh(c*x))**n*(c**S(2)*x**S(2) + S(1))**p*log(h*(f + g*x)**m), x), x)


def replacement6275(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist((-d)**IntPart(p)*(d + e*x**S(2))**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p*log(h*(f + g*x)**m), x), x)


def replacement6276(a, b, c, d1, d2, e1, e2, f, g, h, m, n, p, x):
    return Dist((-d1*d2)**IntPart(p)*(d1 + e1*x)**FracPart(p)*(d2 + e2*x)**FracPart(p)*(c*x + S(-1))**(-FracPart(p))*(c*x + S(1))**(-FracPart(p)), Int((a + b*acosh(c*x))**n*(c*x + S(-1))**p*(c*x + S(1))**p*log(h*(f + g*x)**m), x), x)


def With6277(a, b, c, d, e, f, g, m, x):
    u = IntHide((d + e*x)**m*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/sqrt(c**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(a + b*asinh(c*x), u, x)


def With6278(a, b, c, d, e, f, g, m, x):
    u = IntHide((d + e*x)**m*(f + g*x)**m, x)
    return -Dist(b*c, Int(Dist(S(1)/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), u, x), x), x) + Dist(a + b*acosh(c*x), u, x)


def replacement6279(a, b, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((a + b*asinh(c*x))**n, (d + e*x)**m*(f + g*x)**m, x), x)


def replacement6280(a, b, c, d, e, f, g, m, n, x):
    return Int(ExpandIntegrand((a + b*acosh(c*x))**n, (d + e*x)**m*(f + g*x)**m, x), x)


def With6281(a, b, c, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = IntHide(u, x)
    if InverseFunctionFreeQ(v, x):
        return True
    return False


def replacement6281(a, b, c, u, x):

    v = IntHide(u, x)
    return -Dist(b*c, Int(SimplifyIntegrand(v/sqrt(c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(c*x), v, x)


def With6282(a, b, c, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = IntHide(u, x)
    if InverseFunctionFreeQ(v, x):
        return True
    return False


def replacement6282(a, b, c, u, x):

    v = IntHide(u, x)
    return -Dist(b*c*sqrt(-c**S(2)*x**S(2) + S(1))/(sqrt(c*x + S(-1))*sqrt(c*x + S(1))), Int(SimplifyIntegrand(v/sqrt(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acosh(c*x), v, x)


def With6283(Px, a, b, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)
    if SumQ(u):
        return True
    return False


def replacement6283(Px, a, b, c, d, e, n, p, x):

    u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(d + e*x**S(2))**p, x)
    return Int(u, x)


def With6284(Px, a, b, c, d1, d2, e1, e2, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)
    if SumQ(u):
        return True
    return False


def replacement6284(Px, a, b, c, d1, d2, e1, e2, n, p, x):

    u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(d1 + e1*x)**p*(d2 + e2*x)**p, x)
    return Int(u, x)


def With6285(Px, a, b, c, d, e, f, g, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    if SumQ(u):
        return True
    return False


def replacement6285(Px, a, b, c, d, e, f, g, m, n, p, x):

    u = ExpandIntegrand(Px*(a + b*asinh(c*x))**n*(f + g*(d + e*x**S(2))**p)**m, x)
    return Int(u, x)


def With6286(Px, a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(f + g*(d1 + e1*x)**p*(d2 + e2*x)**p)**m, x)
    if SumQ(u):
        return True
    return False


def replacement6286(Px, a, b, c, d1, d2, e1, e2, f, g, m, n, p, x):

    u = ExpandIntegrand(Px*(a + b*acosh(c*x))**n*(f + g*(d1 + e1*x)**p*(d2 + e2*x)**p)**m, x)
    return Int(u, x)


def With6287(RFx, c, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(asinh(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement6287(RFx, c, n, x):

    u = ExpandIntegrand(asinh(c*x)**n, RFx, x)
    return Int(u, x)


def With6288(RFx, c, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand(acosh(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement6288(RFx, c, n, x):

    u = ExpandIntegrand(acosh(c*x)**n, RFx, x)
    return Int(u, x)


def replacement6289(RFx, a, b, c, n, x):
    return Int(ExpandIntegrand(RFx*(a + b*asinh(c*x))**n, x), x)


def replacement6290(RFx, a, b, c, n, x):
    return Int(ExpandIntegrand(RFx*(a + b*acosh(c*x))**n, x), x)


def With6291(RFx, c, d, e, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((d + e*x**S(2))**p*asinh(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement6291(RFx, c, d, e, n, p, x):

    u = ExpandIntegrand((d + e*x**S(2))**p*asinh(c*x)**n, RFx, x)
    return Int(u, x)


def With6292(RFx, c, d1, d2, e1, e2, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    u = ExpandIntegrand((d1 + e1*x)**p*(d2 + e2*x)**p*acosh(c*x)**n, RFx, x)
    if SumQ(u):
        return True
    return False


def replacement6292(RFx, c, d1, d2, e1, e2, n, p, x):

    u = ExpandIntegrand((d1 + e1*x)**p*(d2 + e2*x)**p*acosh(c*x)**n, RFx, x)
    return Int(u, x)


def replacement6293(RFx, a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((d + e*x**S(2))**p, RFx*(a + b*asinh(c*x))**n, x), x)


def replacement6294(RFx, a, b, c, d1, d2, e1, e2, n, p, x):
    return Int(ExpandIntegrand((d1 + e1*x)**p*(d2 + e2*x)**p, RFx*(a + b*acosh(c*x))**n, x), x)


def replacement6295(a, b, c, n, u, x):
    return Int(u*(a + b*asinh(c*x))**n, x)


def replacement6296(a, b, c, n, u, x):
    return Int(u*(a + b*acosh(c*x))**n, x)


def replacement6297(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*asinh(x))**n, x), x, c + d*x), x)


def replacement6298(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acosh(x))**n, x), x, c + d*x), x)


def replacement6299(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*asinh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6300(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acosh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6301(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*asinh(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x), x)


def replacement6302(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acosh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x), x)


def replacement6303(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*asinh(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6304(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acosh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6305(a, b, c, d, x):
    return Simp(x*sqrt(a + b*asinh(c + d*x**S(2))), x) - Simp(sqrt(Pi)*x*(-c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*FresnelC(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(sqrt(-c/b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x) + Simp(sqrt(Pi)*x*(c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*FresnelS(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(sqrt(-c/b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x)


def replacement6306(a, b, c, d, n, x):
    return Dist(S(4)*b**S(2)*n*(n + S(-1)), Int((a + b*asinh(c + d*x**S(2)))**(n + S(-2)), x), x) + Simp(x*(a + b*asinh(c + d*x**S(2)))**n, x) - Simp(S(2)*b*n*(a + b*asinh(c + d*x**S(2)))**(n + S(-1))*sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(d*x), x)


def replacement6307(a, b, c, d, x):
    return Simp(x*(-c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*SinhIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x) + Simp(x*(c*cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*CoshIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(2)*b*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x)


def replacement6308(a, b, c, d, x):
    return Simp(sqrt(S(2))*sqrt(Pi)*x*(c + S(-1))*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*asinh(c + d*x**S(2)))/(S(2)*sqrt(b)))/(S(4)*sqrt(b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x) + Simp(sqrt(S(2))*sqrt(Pi)*x*(c + S(1))*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*asinh(c + d*x**S(2)))/(S(2)*sqrt(b)))/(S(4)*sqrt(b)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x)


def replacement6309(a, b, c, d, x):
    return -Simp(sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(b*d*x*sqrt(a + b*asinh(c + d*x**S(2)))), x) - Simp(sqrt(Pi)*x*(-c/b)**(S(3)/2)*(-c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*FresnelC(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2))), x) + Simp(sqrt(Pi)*x*(-c/b)**(S(3)/2)*(c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*FresnelS(sqrt(-c/(Pi*b))*sqrt(a + b*asinh(c + d*x**S(2))))/(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2))), x)


def replacement6310(a, b, c, d, x):
    return Simp(x*(-c*sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*CoshIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x) + Simp(x*(c*cosh(a/(S(2)*b)) - sinh(a/(S(2)*b)))*SinhIntegral((a + b*asinh(c + d*x**S(2)))/(S(2)*b))/(S(4)*b**S(2)*(c*sinh(asinh(c + d*x**S(2))/S(2)) + cosh(asinh(c + d*x**S(2))/S(2)))), x) - Simp(sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(S(2)*b*d*x*(a + b*asinh(c + d*x**S(2)))), x)


def replacement6311(a, b, c, d, n, x):
    return Dist(S(1)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), Int((a + b*asinh(c + d*x**S(2)))**(n + S(2)), x), x) - Simp(x*(a + b*asinh(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), x) + Simp((a + b*asinh(c + d*x**S(2)))**(n + S(1))*sqrt(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(S(2)*b*d*x*(n + S(1))), x)


def replacement6312(a, b, d, x):
    return Simp(S(2)*sqrt(a + b*acosh(d*x**S(2) + S(1)))*sinh(acosh(d*x**S(2) + S(1))/S(2))**S(2)/(d*x), x) - Simp(sqrt(S(2))*sqrt(Pi)*sqrt(b)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*d*x), x) + Simp(sqrt(S(2))*sqrt(Pi)*sqrt(b)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*d*x), x)


def replacement6313(a, b, d, x):
    return Simp(S(2)*sqrt(a + b*acosh(d*x**S(2) + S(-1)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))**S(2)/(d*x), x) - Simp(sqrt(S(2))*sqrt(Pi)*sqrt(b)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*d*x), x) - Simp(sqrt(S(2))*sqrt(Pi)*sqrt(b)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*d*x), x)


def replacement6314(a, b, c, d, n, x):
    return Dist(S(4)*b**S(2)*n*(n + S(-1)), Int((a + b*acosh(c + d*x**S(2)))**(n + S(-2)), x), x) + Simp(x*(a + b*acosh(c + d*x**S(2)))**n, x) - Simp(S(2)*b*n*(a + b*acosh(c + d*x**S(2)))**(n + S(-1))*(S(2)*c*d*x**S(2) + d**S(2)*x**S(4))/(d*x*sqrt(c + d*x**S(2) + S(-1))*sqrt(c + d*x**S(2) + S(1))), x)


def replacement6315(a, b, d, x):
    return Simp(sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*cosh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x) - Simp(sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x)


def replacement6316(a, b, d, x):
    return -Simp(sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x) + Simp(sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*cosh(a/(S(2)*b))/(S(2)*b*sqrt(d*x**S(2))), x)


def replacement6317(a, b, d, x):
    return Simp(sqrt(S(2))*sqrt(Pi)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*sqrt(b)*d*x), x) + Simp(sqrt(S(2))*sqrt(Pi)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*sqrt(b)*d*x), x)


def replacement6318(a, b, d, x):
    return Simp(sqrt(S(2))*sqrt(Pi)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*sqrt(b)*d*x), x) - Simp(sqrt(S(2))*sqrt(Pi)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*sqrt(b)*d*x), x)


def replacement6319(a, b, d, x):
    return -Simp(sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(2))/(b*d*x*sqrt(a + b*acosh(d*x**S(2) + S(1)))), x) + Simp(sqrt(S(2))*sqrt(Pi)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*b**(S(3)/2)*d*x), x) - Simp(sqrt(S(2))*sqrt(Pi)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(1)))/(S(2)*sqrt(b)))*sinh(acosh(d*x**S(2) + S(1))/S(2))/(S(2)*b**(S(3)/2)*d*x), x)


def replacement6320(a, b, d, x):
    return -Simp(sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(-2))/(b*d*x*sqrt(a + b*acosh(d*x**S(2) + S(-1)))), x) + Simp(sqrt(S(2))*sqrt(Pi)*(-sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erfi(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*b**(S(3)/2)*d*x), x) + Simp(sqrt(S(2))*sqrt(Pi)*(sinh(a/(S(2)*b)) + cosh(a/(S(2)*b)))*Erf(sqrt(S(2))*sqrt(a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*sqrt(b)))*cosh(acosh(d*x**S(2) + S(-1))/S(2))/(S(2)*b**(S(3)/2)*d*x), x)


def replacement6321(a, b, d, x):
    return -Simp(sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x) + Simp(sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(1)))/(S(2)*b))*cosh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x) - Simp(sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(2))/(S(2)*b*d*x*(a + b*acosh(d*x**S(2) + S(1)))), x)


def replacement6322(a, b, d, x):
    return Simp(sqrt(S(2))*x*CoshIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*cosh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x) - Simp(sqrt(S(2))*x*SinhIntegral((a + b*acosh(d*x**S(2) + S(-1)))/(S(2)*b))*sinh(a/(S(2)*b))/(S(4)*b**S(2)*sqrt(d*x**S(2))), x) - Simp(sqrt(d*x**S(2))*sqrt(d*x**S(2) + S(-2))/(S(2)*b*d*x*(a + b*acosh(d*x**S(2) + S(-1)))), x)


def replacement6323(a, b, c, d, n, x):
    return Dist(S(1)/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), Int((a + b*acosh(c + d*x**S(2)))**(n + S(2)), x), x) - Simp(x*(a + b*acosh(c + d*x**S(2)))**(n + S(2))/(S(4)*b**S(2)*(n + S(1))*(n + S(2))), x) + Simp((a + b*acosh(c + d*x**S(2)))**(n + S(1))*(S(2)*c*x**S(2) + d*x**S(4))/(S(2)*b*x*(n + S(1))*sqrt(c + d*x**S(2) + S(-1))*sqrt(c + d*x**S(2) + S(1))), x)


def replacement6324(a, n, p, x):
    return Dist(S(1)/p, Subst(Int(x**n/tanh(x), x), x, asinh(a*x**p)), x)


def replacement6325(a, n, p, x):
    return Dist(S(1)/p, Subst(Int(x**n*tanh(x), x), x, acosh(a*x**p)), x)


def replacement6326(a, b, c, m, n, u, x):
    return Int(u*acsch(a/c + b*x**n/c)**m, x)


def replacement6327(a, b, c, m, n, u, x):
    return Int(u*asech(a/c + b*x**n/c)**m, x)


def replacement6328(b, n, x):
    return Dist(sqrt(b*x**S(2))/(b*x), Subst(Int(asinh(x)**n/sqrt(x**S(2) + S(1)), x), x, sqrt(b*x**S(2) + S(-1))), x)


def replacement6329(b, n, x):
    return Dist(sqrt(sqrt(b*x**S(2) + S(1)) + S(-1))*sqrt(sqrt(b*x**S(2) + S(1)) + S(1))/(b*x), Subst(Int(acosh(x)**n/(sqrt(x + S(-1))*sqrt(x + S(1))), x), x, sqrt(b*x**S(2) + S(1))), x)


def replacement6330(a, b, c, f, n, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*cosh(x), x), x, asinh(a + b*x)), x)


def replacement6331(a, b, c, f, n, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*sinh(x), x), x, acosh(a + b*x)), x)


def replacement6332(a, b, c, f, m, n, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*(-a/b + sinh(x)/b)**m*cosh(x), x), x, asinh(a + b*x)), x)


def replacement6333(a, b, c, f, m, n, x):
    return Dist(S(1)/b, Subst(Int(f**(c*x**n)*(-a/b + cosh(x)/b)**m*sinh(x), x), x, acosh(a + b*x)), x)


def replacement6334(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/sqrt(u**S(2) + S(1)), x), x) + Simp(x*asinh(u), x)


def replacement6335(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x) + Simp(x*acosh(u), x)


def replacement6336(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/sqrt(u**S(2) + S(1)), x), x), x) + Simp((a + b*asinh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement6337(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x), x) + Simp((a + b*acosh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With6338(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6338(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/sqrt(u**S(2) + S(1)), x), x), x) + Dist(a + b*asinh(u), w, x)


def With6339(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6339(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/(sqrt(u + S(-1))*sqrt(u + S(1))), x), x), x) + Dist(a + b*acosh(u), w, x)


def replacement6340(n, u, x):
    return Int((u + sqrt(u**S(2) + S(1)))**n, x)


def replacement6341(m, n, u, x):
    return Int(x**m*(u + sqrt(u**S(2) + S(1)))**n, x)


def replacement6342(n, u, x):
    return Int((u + sqrt(u + S(-1))*sqrt(u + S(1)))**n, x)


def replacement6343(m, n, u, x):
    return Int(x**m*(u + sqrt(u + S(-1))*sqrt(u + S(1)))**n, x)


def replacement6344(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*atanh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*atanh(c*x))**n, x)


def replacement6345(a, b, c, n, x):
    return -Dist(b*c*n, Int(x*(a + b*acoth(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x*(a + b*acoth(c*x))**n, x)


def replacement6346(a, b, c, n, x):
    return Int((a + b*atanh(c*x))**n, x)


def replacement6347(a, b, c, n, x):
    return Int((a + b*acoth(c*x))**n, x)


def replacement6348(a, b, c, d, e, n, x):
    return Dist(b*c*n/e, Int((a + b*atanh(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x), x) - Simp((a + b*atanh(c*x))**n*log(S(2)*d/(d + e*x))/e, x)


def replacement6349(a, b, c, d, e, n, x):
    return Dist(b*c*n/e, Int((a + b*acoth(c*x))**(n + S(-1))*log(S(2)*d/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x), x) - Simp((a + b*acoth(c*x))**n*log(S(2)*d/(d + e*x))/e, x)


def replacement6350(c, d, e, x):
    return -Simp(PolyLog(S(2), Simp(c*(d + e*x)/(c*d - e), x))/(S(2)*e), x) + Simp(PolyLog(S(2), Simp(c*(d + e*x)/(c*d + e), x))/(S(2)*e), x) - Simp(log(d + e*x)*atanh(c*d/e)/e, x)


def replacement6351(c, d, e, x):
    return -Dist(S(1)/2, Int(log(-c*x + S(1))/(d + e*x), x), x) + Dist(S(1)/2, Int(log(c*x + S(1))/(d + e*x), x), x)


def replacement6352(c, d, e, x):
    return -Dist(S(1)/2, Int(log(S(1) - S(1)/(c*x))/(d + e*x), x), x) + Dist(S(1)/2, Int(log(S(1) + S(1)/(c*x))/(d + e*x), x), x)


def replacement6353(a, b, c, d, e, x):
    return Dist(b, Int(atanh(c*x)/(d + e*x), x), x) + Simp(a*log(RemoveContent(d + e*x, x))/e, x)


def replacement6354(a, b, c, d, e, x):
    return Dist(b, Int(acoth(c*x)/(d + e*x), x), x) + Simp(a*log(RemoveContent(d + e*x, x))/e, x)


def replacement6355(a, b, c, d, e, p, x):
    return -Dist(b*c/(e*(p + S(1))), Int((d + e*x)**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*atanh(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))), x)


def replacement6356(a, b, c, d, e, p, x):
    return -Dist(b*c/(e*(p + S(1))), Int((d + e*x)**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acoth(c*x))*(d + e*x)**(p + S(1))/(e*(p + S(1))), x)


def replacement6357(a, b, c, n, x):
    return -Dist(S(2)*b*c*n, Int((a + b*atanh(c*x))**(n + S(-1))*atanh(S(1) - S(2)/(-c*x + S(1)))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(S(2)*(a + b*atanh(c*x))**n*atanh(S(1) - S(2)/(-c*x + S(1))), x)


def replacement6358(a, b, c, n, x):
    return -Dist(S(2)*b*c*n, Int((a + b*acoth(c*x))**(n + S(-1))*acoth(S(1) - S(2)/(-c*x + S(1)))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(S(2)*(a + b*acoth(c*x))**n*acoth(S(1) - S(2)/(-c*x + S(1))), x)


def replacement6359(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*atanh(c*x))**n/(m + S(1)), x)


def replacement6360(a, b, c, m, n, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp(x**(m + S(1))*(a + b*acoth(c*x))**n/(m + S(1)), x)


def replacement6361(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*atanh(c*x))**n*(d + e*x)**p, x), x)


def replacement6362(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acoth(c*x))**n*(d + e*x)**p, x), x)


def replacement6363(a, b, c, d, e, n, p, x):
    return Int((a + b*atanh(c*x))**n*(d + e*x)**p, x)


def replacement6364(a, b, c, d, e, n, p, x):
    return Int((a + b*acoth(c*x))**n*(d + e*x)**p, x)


def replacement6365(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-1))*(a + b*atanh(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-1))*(a + b*atanh(c*x))**n/(d + e*x), x), x)


def replacement6366(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-1))*(a + b*acoth(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-1))*(a + b*acoth(c*x))**n/(d + e*x), x), x)


def replacement6367(a, b, c, d, e, n, x):
    return -Dist(b*c*n/d, Int((a + b*atanh(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*atanh(c*x))**n*log(S(2)*e*x/(d + e*x))/d, x)


def replacement6368(a, b, c, d, e, n, x):
    return -Dist(b*c*n/d, Int((a + b*acoth(c*x))**(n + S(-1))*log(S(2)*e*x/(d + e*x))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acoth(c*x))**n*log(S(2)*e*x/(d + e*x))/d, x)


def replacement6369(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*atanh(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(1))*(a + b*atanh(c*x))**n/(d + e*x), x), x)


def replacement6370(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acoth(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(1))*(a + b*acoth(c*x))**n/(d + e*x), x), x)


def replacement6371(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*atanh(c*x))**n*(d + e*x)**p, x), x)


def replacement6372(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*acoth(c*x))**n*(d + e*x)**p, x), x)


def replacement6373(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*atanh(c*x))**n*(d + e*x)**p, x)


def replacement6374(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acoth(c*x))**n*(d + e*x)**p, x)


def replacement6375(a, b, c, d, e, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*atanh(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement6376(a, b, c, d, e, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*acoth(c*x))*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement6377(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b**S(2)*d*n*(n + S(-1))/(S(2)*p*(S(2)*p + S(1))), Int((a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*n*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement6378(a, b, c, d, e, n, p, x):
    return Dist(S(2)*d*p/(S(2)*p + S(1)), Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(b**S(2)*d*n*(n + S(-1))/(S(2)*p*(S(2)*p + S(1))), Int((a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**(p + S(-1)), x), x) + Simp(x*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p/(S(2)*p + S(1)), x) + Simp(b*n*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p/(S(2)*c*p*(S(2)*p + S(1))), x)


def replacement6379(a, b, c, d, e, x):
    return Simp(log(RemoveContent(a + b*atanh(c*x), x))/(b*c*d), x)


def replacement6380(a, b, c, d, e, x):
    return Simp(log(RemoveContent(a + b*acoth(c*x), x))/(b*c*d), x)


def replacement6381(a, b, c, d, e, n, x):
    return Simp((a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6382(a, b, c, d, e, n, x):
    return Simp((a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6383(a, b, c, d, e, x):
    return Simp(-S(2)*(a + b*atanh(c*x))*ArcTan(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x) - Simp(I*b*PolyLog(S(2), -I*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x) + Simp(I*b*PolyLog(S(2), I*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x)


def replacement6384(a, b, c, d, e, x):
    return Simp(-S(2)*(a + b*acoth(c*x))*ArcTan(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x) - Simp(I*b*PolyLog(S(2), -I*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x) + Simp(I*b*PolyLog(S(2), I*sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/(c*sqrt(d)), x)


def replacement6385(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*sqrt(d)), Subst(Int((a + b*x)**n/cosh(x), x), x, atanh(c*x)), x)


def replacement6386(a, b, c, d, e, n, x):
    return -Dist(x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n/sinh(x), x), x, acoth(c*x)), x)


def replacement6387(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*atanh(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def replacement6388(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*acoth(c*x))**n/sqrt(-c**S(2)*x**S(2) + S(1)), x), x)


def replacement6389(a, b, c, d, e, n, x):
    return -Dist(b*c*n/S(2), Int(x*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) + Simp(x*(a + b*atanh(c*x))**n/(S(2)*d*(d + e*x**S(2))), x) + Simp((a + b*atanh(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))), x)


def replacement6390(a, b, c, d, e, n, x):
    return -Dist(b*c*n/S(2), Int(x*(a + b*acoth(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) + Simp(x*(a + b*acoth(c*x))**n/(S(2)*d*(d + e*x**S(2))), x) + Simp((a + b*acoth(c*x))**(n + S(1))/(S(2)*b*c*d**S(2)*(n + S(1))), x)


def replacement6391(a, b, c, d, e, x):
    return -Simp(b/(c*d*sqrt(d + e*x**S(2))), x) + Simp(x*(a + b*atanh(c*x))/(d*sqrt(d + e*x**S(2))), x)


def replacement6392(a, b, c, d, e, x):
    return -Simp(b/(c*d*sqrt(d + e*x**S(2))), x) + Simp(x*(a + b*acoth(c*x))/(d*sqrt(d + e*x**S(2))), x)


def replacement6393(a, b, c, d, e, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement6394(a, b, c, d, e, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x)


def replacement6395(a, b, c, d, e, n, x):
    return Dist(b**S(2)*n*(n + S(-1)), Int((a + b*atanh(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x), x) + Simp(x*(a + b*atanh(c*x))**n/(d*sqrt(d + e*x**S(2))), x) - Simp(b*n*(a + b*atanh(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))), x)


def replacement6396(a, b, c, d, e, n, x):
    return Dist(b**S(2)*n*(n + S(-1)), Int((a + b*acoth(c*x))**(n + S(-2))/(d + e*x**S(2))**(S(3)/2), x), x) + Simp(x*(a + b*acoth(c*x))**n/(d*sqrt(d + e*x**S(2))), x) - Simp(b*n*(a + b*acoth(c*x))**(n + S(-1))/(c*d*sqrt(d + e*x**S(2))), x)


def replacement6397(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b**S(2)*n*(n + S(-1))/(S(4)*(p + S(1))**S(2)), Int((a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Simp(x*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x) - Simp(b*n*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x)


def replacement6398(a, b, c, d, e, n, p, x):
    return Dist((S(2)*p + S(3))/(S(2)*d*(p + S(1))), Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Dist(b**S(2)*n*(n + S(-1))/(S(4)*(p + S(1))**S(2)), Int((a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Simp(x*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*d*(p + S(1))), x) - Simp(b*n*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(S(4)*c*d*(p + S(1))**S(2)), x)


def replacement6399(a, b, c, d, e, n, p, x):
    return Dist(S(2)*c*(p + S(1))/(b*(n + S(1))), Int(x*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6400(a, b, c, d, e, n, p, x):
    return Dist(S(2)*c*(p + S(1))/(b*(n + S(1))), Int(x*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6401(a, b, c, d, e, n, p, x):
    return Dist(d**p/c, Subst(Int((a + b*x)**n*cosh(x)**(-S(2)*p + S(-2)), x), x, atanh(c*x)), x)


def replacement6402(a, b, c, d, e, n, p, x):
    return Dist(d**(p + S(1)/2)*sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*atanh(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement6403(a, b, c, d, e, n, p, x):
    return -Dist((-d)**p/c, Subst(Int((a + b*x)**n*sinh(x)**(-S(2)*p + S(-2)), x), x, acoth(c*x)), x)


def replacement6404(a, b, c, d, e, n, p, x):
    return -Dist(x*(-d)**(p + S(1)/2)*sqrt((c**S(2)*x**S(2) + S(-1))/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n*sinh(x)**(-S(2)*p + S(-2)), x), x, acoth(c*x)), x)


def replacement6405(c, d, e, x):
    return -Dist(S(1)/2, Int(log(-c*x + S(1))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int(log(c*x + S(1))/(d + e*x**S(2)), x), x)


def replacement6406(c, d, e, x):
    return -Dist(S(1)/2, Int(log(S(1) - S(1)/(c*x))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int(log(S(1) + S(1)/(c*x))/(d + e*x**S(2)), x), x)


def replacement6407(a, b, c, d, e, x):
    return Dist(a, Int(S(1)/(d + e*x**S(2)), x), x) + Dist(b, Int(atanh(c*x)/(d + e*x**S(2)), x), x)


def replacement6408(a, b, c, d, e, x):
    return Dist(a, Int(S(1)/(d + e*x**S(2)), x), x) + Dist(b, Int(acoth(c*x)/(d + e*x**S(2)), x), x)


def With6409(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*atanh(c*x), u, x)


def With6410(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acoth(c*x), u, x)


def replacement6411(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6412(a, b, c, d, e, n, p, x):
    return Int(ExpandIntegrand((a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6413(a, b, c, d, e, n, p, x):
    return Int((a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6414(a, b, c, d, e, n, p, x):
    return Int((a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6415(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6416(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6417(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*atanh(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6418(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acoth(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6419(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*d), Int((a + b*atanh(c*x))**n/(-c*x + S(1)), x), x) + Simp((a + b*atanh(c*x))**(n + S(1))/(b*e*(n + S(1))), x)


def replacement6420(a, b, c, d, e, n, x):
    return Dist(S(1)/(c*d), Int((a + b*acoth(c*x))**n/(-c*x + S(1)), x), x) + Simp((a + b*acoth(c*x))**(n + S(1))/(b*e*(n + S(1))), x)


def replacement6421(a, b, c, d, e, n, x):
    return -Dist(S(1)/(b*c*d*(n + S(1))), Int((a + b*atanh(c*x))**(n + S(1)), x), x) + Simp(x*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6422(a, b, c, d, e, n, x):
    return -Dist(S(1)/(b*c*d*(n + S(1))), Int((a + b*acoth(c*x))**(n + S(1)), x), x) - Simp(x*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6423(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6424(a, b, c, d, e, m, n, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n, x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6425(a, b, c, d, e, n, x):
    return Dist(S(1)/d, Int((a + b*atanh(c*x))**n/(x*(c*x + S(1))), x), x) + Simp((a + b*atanh(c*x))**(n + S(1))/(b*d*(n + S(1))), x)


def replacement6426(a, b, c, d, e, n, x):
    return Dist(S(1)/d, Int((a + b*acoth(c*x))**n/(x*(c*x + S(1))), x), x) + Simp((a + b*acoth(c*x))**(n + S(1))/(b*d*(n + S(1))), x)


def replacement6427(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*atanh(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*atanh(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6428(a, b, c, d, e, m, n, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acoth(c*x))**n, x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acoth(c*x))**n/(d + e*x**S(2)), x), x)


def replacement6429(a, b, c, d, e, m, n, x):
    return -Dist(m/(b*c*d*(n + S(1))), Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1)), x), x) + Simp(x**m*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6430(a, b, c, d, e, m, n, x):
    return -Dist(m/(b*c*d*(n + S(1))), Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1)), x), x) + Simp(x**m*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(n + S(1))), x)


def replacement6431(a, b, c, d, e, m, x):
    return Int(ExpandIntegrand(a + b*atanh(c*x), x**m/(d + e*x**S(2)), x), x)


def replacement6432(a, b, c, d, e, m, x):
    return Int(ExpandIntegrand(a + b*acoth(c*x), x**m/(d + e*x**S(2)), x), x)


def replacement6433(a, b, c, d, e, n, p, x):
    return Dist(b*n/(S(2)*c*(p + S(1))), Int((a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6434(a, b, c, d, e, n, p, x):
    return Dist(b*n/(S(2)*c*(p + S(1))), Int((a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp((a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6435(a, b, c, d, e, n, x):
    return Dist(S(4)/(b**S(2)*(n + S(1))*(n + S(2))), Int(x*(a + b*atanh(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x), x) + Simp((a + b*atanh(c*x))**(n + S(2))*(c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))), x) + Simp(x*(a + b*atanh(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))), x)


def replacement6436(a, b, c, d, e, n, x):
    return Dist(S(4)/(b**S(2)*(n + S(1))*(n + S(2))), Int(x*(a + b*acoth(c*x))**(n + S(2))/(d + e*x**S(2))**S(2), x), x) + Simp((a + b*acoth(c*x))**(n + S(2))*(c**S(2)*x**S(2) + S(1))/(b**S(2)*e*(d + e*x**S(2))*(n + S(1))*(n + S(2))), x) + Simp(x*(a + b*acoth(c*x))**(n + S(1))/(b*c*d*(d + e*x**S(2))*(n + S(1))), x)


def replacement6437(a, b, c, d, e, p, x):
    return Dist(S(1)/(S(2)*c**S(2)*d*(p + S(1))), Int((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))), x)


def replacement6438(a, b, c, d, e, p, x):
    return Dist(S(1)/(S(2)*c**S(2)*d*(p + S(1))), Int((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*(d + e*x**S(2))**(p + S(1))/(S(4)*c**S(3)*d*(p + S(1))**S(2)), x) - Simp(x*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*c**S(2)*d*(p + S(1))), x)


def replacement6439(a, b, c, d, e, n, x):
    return -Dist(b*n/(S(2)*c), Int(x*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) - Simp((a + b*atanh(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))), x) + Simp(x*(a + b*atanh(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))), x)


def replacement6440(a, b, c, d, e, n, x):
    return -Dist(b*n/(S(2)*c), Int(x*(a + b*acoth(c*x))**(n + S(-1))/(d + e*x**S(2))**S(2), x), x) - Simp((a + b*acoth(c*x))**(n + S(1))/(S(2)*b*c**S(3)*d**S(2)*(n + S(1))), x) + Simp(x*(a + b*acoth(c*x))**n/(S(2)*c**S(2)*d*(d + e*x**S(2))), x)


def replacement6441(a, b, c, d, e, m, p, x):
    return -Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x) + Simp(x**(m + S(-1))*(a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x)


def replacement6442(a, b, c, d, e, m, p, x):
    return -Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1)), x), x) - Simp(b*x**m*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x) + Simp(x**(m + S(-1))*(a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x)


def replacement6443(a, b, c, d, e, m, n, p, x):
    return Dist(b**S(2)*n*(n + S(-1))/m**S(2), Int(x**m*(a + b*atanh(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Simp(x**(m + S(-1))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x) - Simp(b*n*x**m*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x)


def replacement6444(a, b, c, d, e, m, n, p, x):
    return Dist(b**S(2)*n*(n + S(-1))/m**S(2), Int(x**m*(a + b*acoth(c*x))**(n + S(-2))*(d + e*x**S(2))**p, x), x) - Dist((m + S(-1))/(c**S(2)*d*m), Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) + Simp(x**(m + S(-1))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(c**S(2)*d*m), x) - Simp(b*n*x**m*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**(p + S(1))/(c*d*m**S(2)), x)


def replacement6445(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6446(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6447(a, b, c, d, e, m, n, p, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp(x**(m + S(1))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))), x)


def replacement6448(a, b, c, d, e, m, n, p, x):
    return -Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))*(d + e*x**S(2))**p, x), x) + Simp(x**(m + S(1))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1))/(d*(m + S(1))), x)


def replacement6449(a, b, c, d, e, m, x):
    return Dist(d/(m + S(2)), Int(x**m*(a + b*atanh(c*x))/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*d/(m + S(2)), Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*atanh(c*x))*sqrt(d + e*x**S(2))/(m + S(2)), x)


def replacement6450(a, b, c, d, e, m, x):
    return Dist(d/(m + S(2)), Int(x**m*(a + b*acoth(c*x))/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*d/(m + S(2)), Int(x**(m + S(1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acoth(c*x))*sqrt(d + e*x**S(2))/(m + S(2)), x)


def replacement6451(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6452(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6453(a, b, c, d, e, m, n, p, x):
    return Dist(d, Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(c**S(2)*d, Int(x**(m + S(2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x)


def replacement6454(a, b, c, d, e, m, n, p, x):
    return Dist(d, Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x) - Dist(c**S(2)*d, Int(x**(m + S(2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(-1)), x), x)


def replacement6455(a, b, c, d, e, m, n, x):
    return Dist((m + S(-1))/(c**S(2)*m), Int(x**(m + S(-2))*(a + b*atanh(c*x))**n/sqrt(d + e*x**S(2)), x), x) + Dist(b*n/(c*m), Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) - Simp(x**(m + S(-1))*(a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m), x)


def replacement6456(a, b, c, d, e, m, n, x):
    return Dist((m + S(-1))/(c**S(2)*m), Int(x**(m + S(-2))*(a + b*acoth(c*x))**n/sqrt(d + e*x**S(2)), x), x) + Dist(b*n/(c*m), Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) - Simp(x**(m + S(-1))*(a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(c**S(2)*d*m), x)


def replacement6457(a, b, c, d, e, x):
    return Simp(b*PolyLog(S(2), -sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x) - Simp(b*PolyLog(S(2), sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x) + Simp(-S(2)*(a + b*atanh(c*x))*atanh(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x)


def replacement6458(a, b, c, d, e, x):
    return Simp(b*PolyLog(S(2), -sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x) - Simp(b*PolyLog(S(2), sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x) + Simp(-S(2)*(a + b*acoth(c*x))*atanh(sqrt(-c*x + S(1))/sqrt(c*x + S(1)))/sqrt(d), x)


def replacement6459(a, b, c, d, e, n, x):
    return Dist(S(1)/sqrt(d), Subst(Int((a + b*x)**n/sinh(x), x), x, atanh(c*x)), x)


def replacement6460(a, b, c, d, e, n, x):
    return -Dist(c*x*sqrt(S(1) - S(1)/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n/cosh(x), x), x, acoth(c*x)), x)


def replacement6461(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*atanh(c*x))**n/(x*sqrt(-c**S(2)*x**S(2) + S(1))), x), x)


def replacement6462(a, b, c, d, e, n, x):
    return Dist(sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int((a + b*acoth(c*x))**n/(x*sqrt(-c**S(2)*x**S(2) + S(1))), x), x)


def replacement6463(a, b, c, d, e, n, x):
    return Dist(b*c*n, Int((a + b*atanh(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x), x) - Simp((a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(d*x), x)


def replacement6464(a, b, c, d, e, n, x):
    return Dist(b*c*n, Int((a + b*acoth(c*x))**(n + S(-1))/(x*sqrt(d + e*x**S(2))), x), x) - Simp((a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(d*x), x)


def replacement6465(a, b, c, d, e, m, n, x):
    return Dist(c**S(2)*(m + S(2))/(m + S(1)), Int(x**(m + S(2))*(a + b*atanh(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*atanh(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))), x)


def replacement6466(a, b, c, d, e, m, n, x):
    return Dist(c**S(2)*(m + S(2))/(m + S(1)), Int(x**(m + S(2))*(a + b*acoth(c*x))**n/sqrt(d + e*x**S(2)), x), x) - Dist(b*c*n/(m + S(1)), Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(-1))/sqrt(d + e*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acoth(c*x))**n*sqrt(d + e*x**S(2))/(d*(m + S(1))), x)


def replacement6467(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6468(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(d/e, Int(x**(m + S(-2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6469(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/d, Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6470(a, b, c, d, e, m, n, p, x):
    return Dist(S(1)/d, Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**(p + S(1)), x), x) - Dist(e/d, Int(x**(m + S(2))*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x), x)


def replacement6471(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Dist(c*(m + S(2)*p + S(2))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*atanh(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6472(a, b, c, d, e, m, n, p, x):
    return -Dist(m/(b*c*(n + S(1))), Int(x**(m + S(-1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Dist(c*(m + S(2)*p + S(2))/(b*(n + S(1))), Int(x**(m + S(1))*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**p, x), x) + Simp(x**m*(a + b*acoth(c*x))**(n + S(1))*(d + e*x**S(2))**(p + S(1))/(b*c*d*(n + S(1))), x)


def replacement6473(a, b, c, d, e, m, n, p, x):
    return Dist(c**(-m + S(-1))*d**p, Subst(Int((a + b*x)**n*sinh(x)**m*cosh(x)**(-m - S(2)*p + S(-2)), x), x, atanh(c*x)), x)


def replacement6474(a, b, c, d, e, m, n, p, x):
    return Dist(d**(p + S(1)/2)*sqrt(-c**S(2)*x**S(2) + S(1))/sqrt(d + e*x**S(2)), Int(x**m*(a + b*atanh(c*x))**n*(-c**S(2)*x**S(2) + S(1))**p, x), x)


def replacement6475(a, b, c, d, e, m, n, p, x):
    return -Dist(c**(-m + S(-1))*(-d)**p, Subst(Int((a + b*x)**n*sinh(x)**(-m - S(2)*p + S(-2))*cosh(x)**m, x), x, acoth(c*x)), x)


def replacement6476(a, b, c, d, e, m, n, p, x):
    return -Dist(c**(-m)*x*(-d)**(p + S(1)/2)*sqrt((c**S(2)*x**S(2) + S(-1))/(c**S(2)*x**S(2)))/sqrt(d + e*x**S(2)), Subst(Int((a + b*x)**n*sinh(x)**(-m - S(2)*p + S(-2))*cosh(x)**m, x), x, acoth(c*x)), x)


def replacement6477(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*atanh(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6478(a, b, c, d, e, p, x):
    return -Dist(b*c/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(-c**S(2)*x**S(2) + S(1)), x), x) + Simp((a + b*acoth(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With6479(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*atanh(c*x), u, x)


def With6480(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c, Int(SimplifyIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acoth(c*x), u, x)


def replacement6481(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand((a + b*atanh(c*x))**n, x**m*(d + e*x**S(2))**p, x), x)


def replacement6482(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand((a + b*acoth(c*x))**n, x**m*(d + e*x**S(2))**p, x), x)


def replacement6483(a, b, c, d, e, m, p, x):
    return Dist(a, Int(x**m*(d + e*x**S(2))**p, x), x) + Dist(b, Int(x**m*(d + e*x**S(2))**p*atanh(c*x), x), x)


def replacement6484(a, b, c, d, e, m, p, x):
    return Dist(a, Int(x**m*(d + e*x**S(2))**p, x), x) + Dist(b, Int(x**m*(d + e*x**S(2))**p*acoth(c*x), x), x)


def replacement6485(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*atanh(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6486(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acoth(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6487(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*atanh(c*x))**n*log(S(1) - u)/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*atanh(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x), x)


def replacement6488(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) + S(1)/u, x))/(d + e*x**S(2)), x), x)


def replacement6489(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*atanh(c*x))**n*log(S(1) - u)/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*atanh(c*x))**n*log(u + S(1))/(d + e*x**S(2)), x), x)


def replacement6490(a, b, c, d, e, n, u, x):
    return -Dist(S(1)/2, Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) - S(1)/u, x))/(d + e*x**S(2)), x), x) + Dist(S(1)/2, Int((a + b*acoth(c*x))**n*log(SimplifyIntegrand(S(1) + S(1)/u, x))/(d + e*x**S(2)), x), x)


def replacement6491(a, b, c, d, e, n, u, x):
    return -Dist(b*n/S(2), Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) + Simp((a + b*atanh(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement6492(a, b, c, d, e, n, u, x):
    return -Dist(b*n/S(2), Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) + Simp((a + b*acoth(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement6493(a, b, c, d, e, n, u, x):
    return Dist(b*n/S(2), Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) - Simp((a + b*atanh(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement6494(a, b, c, d, e, n, u, x):
    return Dist(b*n/S(2), Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(S(2), Together(S(1) - u))/(d + e*x**S(2)), x), x) - Simp((a + b*acoth(c*x))**n*PolyLog(S(2), Together(S(1) - u))/(S(2)*c*d), x)


def replacement6495(a, b, c, d, e, n, p, u, x):
    return Dist(b*n/S(2), Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) - Simp((a + b*atanh(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement6496(a, b, c, d, e, n, p, u, x):
    return Dist(b*n/S(2), Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) - Simp((a + b*acoth(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement6497(a, b, c, d, e, n, p, u, x):
    return -Dist(b*n/S(2), Int((a + b*atanh(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) + Simp((a + b*atanh(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement6498(a, b, c, d, e, n, p, u, x):
    return -Dist(b*n/S(2), Int((a + b*acoth(c*x))**(n + S(-1))*PolyLog(p + S(1), u)/(d + e*x**S(2)), x), x) + Simp((a + b*acoth(c*x))**n*PolyLog(p + S(1), u)/(S(2)*c*d), x)


def replacement6499(a, b, c, d, e, x):
    return Simp((-log(a + b*acoth(c*x)) + log(a + b*atanh(c*x)))/(b**S(2)*c*d*(acoth(c*x) - atanh(c*x))), x)


def replacement6500(a, b, c, d, e, m, n, x):
    return -Dist(n/(m + S(1)), Int((a + b*acoth(c*x))**(m + S(1))*(a + b*atanh(c*x))**(n + S(-1))/(d + e*x**S(2)), x), x) + Simp((a + b*acoth(c*x))**(m + S(1))*(a + b*atanh(c*x))**n/(b*c*d*(m + S(1))), x)


def replacement6501(a, b, c, d, e, m, n, x):
    return -Dist(n/(m + S(1)), Int((a + b*acoth(c*x))**(n + S(-1))*(a + b*atanh(c*x))**(m + S(1))/(d + e*x**S(2)), x), x) + Simp((a + b*acoth(c*x))**n*(a + b*atanh(c*x))**(m + S(1))/(b*c*d*(m + S(1))), x)


def replacement6502(a, c, d, n, x):
    return -Dist(S(1)/2, Int(log(-a*x + S(1))/(c + d*x**n), x), x) + Dist(S(1)/2, Int(log(a*x + S(1))/(c + d*x**n), x), x)


def replacement6503(a, c, d, n, x):
    return -Dist(S(1)/2, Int(log(S(1) - S(1)/(a*x))/(c + d*x**n), x), x) + Dist(S(1)/2, Int(log(S(1) + S(1)/(a*x))/(c + d*x**n), x), x)


def replacement6504(a, b, c, d, e, f, g, x):
    return -Dist(b*c, Int(x*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g, Int(x**S(2)*(a + b*atanh(c*x))/(f + g*x**S(2)), x), x) + Simp(x*(a + b*atanh(c*x))*(d + e*log(f + g*x**S(2))), x)


def replacement6505(a, b, c, d, e, f, g, x):
    return -Dist(b*c, Int(x*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g, Int(x**S(2)*(a + b*acoth(c*x))/(f + g*x**S(2)), x), x) + Simp(x*(a + b*acoth(c*x))*(d + e*log(f + g*x**S(2))), x)


def replacement6506(a, b, c, d, e, f, g, m, x):
    return -Dist(b*c/(m + S(1)), Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g/(m + S(1)), Int(x**(m + S(2))*(a + b*atanh(c*x))/(f + g*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*atanh(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)), x)


def replacement6507(a, b, c, d, e, f, g, m, x):
    return -Dist(b*c/(m + S(1)), Int(x**(m + S(1))*(d + e*log(f + g*x**S(2)))/(-c**S(2)*x**S(2) + S(1)), x), x) - Dist(S(2)*e*g/(m + S(1)), Int(x**(m + S(2))*(a + b*acoth(c*x))/(f + g*x**S(2)), x), x) + Simp(x**(m + S(1))*(a + b*acoth(c*x))*(d + e*log(f + g*x**S(2)))/(m + S(1)), x)


def With6508(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*atanh(c*x), u, x)


def With6509(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(d + e*log(f + g*x**S(2))), x)
    return -Dist(b*c, Int(ExpandIntegrand(u/(-c**S(2)*x**S(2) + S(1)), x), x), x) + Dist(a + b*acoth(c*x), u, x)


def With6510(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(a + b*atanh(c*x)), x)
    return -Dist(S(2)*e*g, Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)


def With6511(a, b, c, d, e, f, g, m, x):
    u = IntHide(x**m*(a + b*acoth(c*x)), x)
    return -Dist(S(2)*e*g, Int(ExpandIntegrand(u*x/(f + g*x**S(2)), x), x), x) + Dist(d + e*log(f + g*x**S(2)), u, x)


def replacement6512(a, b, c, d, e, f, g, x):
    return Dist(b/c, Int((a + b*atanh(c*x))*(d + e*log(f + g*x**S(2))), x), x) + Dist(b*c*e, Int(x**S(2)*(a + b*atanh(c*x))/(-c**S(2)*x**S(2) + S(1)), x), x) - Simp(e*x**S(2)*(a + b*atanh(c*x))**S(2)/S(2), x) + Simp((a + b*atanh(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g), x)


def replacement6513(a, b, c, d, e, f, g, x):
    return Dist(b/c, Int((a + b*acoth(c*x))*(d + e*log(f + g*x**S(2))), x), x) + Dist(b*c*e, Int(x**S(2)*(a + b*acoth(c*x))/(-c**S(2)*x**S(2) + S(1)), x), x) - Simp(e*x**S(2)*(a + b*acoth(c*x))**S(2)/S(2), x) + Simp((a + b*acoth(c*x))**S(2)*(d + e*log(f + g*x**S(2)))*(f + g*x**S(2))/(S(2)*g), x)


def replacement6514(a, n, x):
    return Int((-a*x + S(1))**(S(1)/2 - n/S(2))*(a*x + S(1))**(n/S(2) + S(1)/2)/sqrt(-a**S(2)*x**S(2) + S(1)), x)


def replacement6515(a, m, n, x):
    return Int(x**m*(-a*x + S(1))**(S(1)/2 - n/S(2))*(a*x + S(1))**(n/S(2) + S(1)/2)/sqrt(-a**S(2)*x**S(2) + S(1)), x)


def replacement6516(a, n, x):
    return Int((-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x)


def replacement6517(a, m, n, x):
    return Int(x**m*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x)


def replacement6518(a, c, d, n, p, x):
    return Dist(c**n, Int((c + d*x)**(-n + p)*(-a**S(2)*x**S(2) + S(1))**(n/S(2)), x), x)


def replacement6519(a, c, d, e, f, m, n, p, x):
    return Dist(c**n, Int((c + d*x)**(-n + p)*(e + f*x)**m*(-a**S(2)*x**S(2) + S(1))**(n/S(2)), x), x)


def replacement6520(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(S(1) + d*x/c)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x), x)


def replacement6521(a, c, d, n, p, u, x):
    return Int(u*(c + d*x)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x)


def replacement6522(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6523(a, c, d, n, p, u, x):
    return Dist((S(-1))**(n/S(2))*c**p, Int(u*(S(1) - S(1)/(a*x))**(-n/S(2))*(S(1) + S(1)/(a*x))**(n/S(2))*(S(1) + d/(c*x))**p, x), x)


def replacement6524(a, c, d, n, p, u, x):
    return Int(u*(c + d/x)**p*(-a*x + S(1))**(-n/S(2))*(a*x + S(1))**(n/S(2)), x)


def replacement6525(a, c, d, n, p, u, x):
    return Dist(x**p*(c + d/x)**p*(c*x/d + S(1))**(-p), Int(u*x**(-p)*(c*x/d + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6526(a, c, d, n, x):
    return Simp((-a*x + n)*exp(n*atanh(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))), x)


def replacement6527(a, c, d, n, p, x):
    return -Dist(S(2)*(p + S(1))*(S(2)*p + S(3))/(c*(n**S(2) - S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*atanh(a*x))/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))), x)


def replacement6528(a, c, d, n, x):
    return Simp(exp(n*atanh(a*x))/(a*c*n), x)


def replacement6529(a, c, d, n, p, x):
    return Dist(c**p, Int((a*x + S(1))**n*(-a**S(2)*x**S(2) + S(1))**(-n/S(2) + p), x), x)


def replacement6530(a, c, d, n, p, x):
    return Dist(c**p, Int((-a*x + S(1))**(-n)*(-a**S(2)*x**S(2) + S(1))**(n/S(2) + p), x), x)


def replacement6531(a, c, d, n, p, x):
    return Dist(c**p, Int((-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x), x)


def replacement6532(a, c, d, n, p, x):
    return Dist(c**(n/S(2)), Int((c + d*x**S(2))**(-n/S(2) + p)*(a*x + S(1))**n, x), x)


def replacement6533(a, c, d, n, p, x):
    return Dist(c**(-n/S(2)), Int((c + d*x**S(2))**(n/S(2) + p)*(-a*x + S(1))**(-n), x), x)


def replacement6534(a, c, d, n, p, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int((-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6535(a, c, d, n, x):
    return Simp((-a*n*x + S(1))*exp(n*atanh(a*x))/(d*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))), x)


def replacement6536(a, c, d, n, p, x):
    return -Dist(a*c*n/(S(2)*d*(p + S(1))), Int((c + d*x**S(2))**p*exp(n*atanh(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x))/(S(2)*d*(p + S(1))), x)


def replacement6537(a, c, d, n, p, x):
    return Simp((c + d*x**S(2))**(p + S(1))*(-a*n*x + S(1))*exp(n*atanh(a*x))/(a*d*n*(n**S(2) + S(-1))), x)


def replacement6538(a, c, d, n, p, x):
    return Dist((n**S(2) + S(2)*p + S(2))/(d*(n**S(2) - S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*atanh(a*x)), x), x) - Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*atanh(a*x))/(a*d*(n**S(2) - S(4)*(p + S(1))**S(2))), x)


def replacement6539(a, c, d, m, n, p, x):
    return Dist(c**p, Int(x**m*(a*x + S(1))**n*(-a**S(2)*x**S(2) + S(1))**(-n/S(2) + p), x), x)


def replacement6540(a, c, d, m, n, p, x):
    return Dist(c**p, Int(x**m*(-a*x + S(1))**(-n)*(-a**S(2)*x**S(2) + S(1))**(n/S(2) + p), x), x)


def replacement6541(a, c, d, m, n, p, x):
    return Dist(c**p, Int(x**m*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x), x)


def replacement6542(a, c, d, m, n, p, x):
    return Dist(c**(n/S(2)), Int(x**m*(c + d*x**S(2))**(-n/S(2) + p)*(a*x + S(1))**n, x), x)


def replacement6543(a, c, d, m, n, p, x):
    return Dist(c**(-n/S(2)), Int(x**m*(c + d*x**S(2))**(n/S(2) + p)*(-a*x + S(1))**(-n), x), x)


def replacement6544(a, c, d, m, n, p, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(x**m*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6545(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x), x)


def replacement6546(a, c, d, n, p, u, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a*x + S(1))**(-FracPart(p))*(a*x + S(1))**(-FracPart(p)), Int(u*(-a*x + S(1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x), x)


def replacement6547(a, c, d, n, p, u, x):
    return Dist(c**IntPart(p)*(c + d*x**S(2))**FracPart(p)*(-a**S(2)*x**S(2) + S(1))**(-FracPart(p)), Int(u*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6548(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(-S(2)*p)*(-a**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6549(a, c, d, n, p, u, x):
    return Dist(c**p, Int(u*(S(1) - S(1)/(a*x))**p*(S(1) + S(1)/(a*x))**p*exp(n*atanh(a*x)), x), x)


def replacement6550(a, c, d, n, p, u, x):
    return Dist(x**(S(2)*p)*(c + d/x**S(2))**p*(-a*x + S(1))**(-p)*(a*x + S(1))**(-p), Int(u*x**(-S(2)*p)*(-a*x + S(1))**p*(a*x + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6551(a, c, d, n, p, u, x):
    return Dist(x**(S(2)*p)*(c + d/x**S(2))**p*(c*x**S(2)/d + S(1))**(-p), Int(u*x**(-S(2)*p)*(c*x**S(2)/d + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6552(a, b, c, n, x):
    return Int((-a*c - b*c*x + S(1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x)


def replacement6553(a, b, c, m, n, x):
    return Dist(S(4)*b**(-m + S(-1))*c**(-m + S(-1))/n, Subst(Int(x**(S(2)/n)*(x**(S(2)/n) + S(1))**(-m + S(-2))*(-a*c + x**(S(2)/n)*(-a*c + S(1)) + S(-1))**m, x), x, (-c*(a + b*x) + S(1))**(-n/S(2))*(c*(a + b*x) + S(1))**(n/S(2))), x)


def replacement6554(a, b, c, d, e, m, n, x):
    return Int((d + e*x)**m*(-a*c - b*c*x + S(1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x)


def replacement6555(a, b, c, d, e, n, p, u, x):
    return Dist((c/(S(1) - a**S(2)))**p, Int(u*(-a - b*x + S(1))**(-n/S(2) + p)*(a + b*x + S(1))**(n/S(2) + p), x), x)


def replacement6556(a, b, c, d, e, n, p, u, x):
    return Dist((c + d*x + e*x**S(2))**p*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**(-p), Int(u*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**p*exp(n*atanh(a*x)), x), x)


def replacement6557(a, b, c, n, u, x):
    return Int(u*exp(n*acoth(a/c + b*x/c)), x)


def replacement6558(a, n, u, x):
    return Dist((S(-1))**(n/S(2)), Int(u*exp(n*atanh(a*x)), x), x)


def replacement6559(a, n, x):
    return -Subst(Int((S(1) - x/a)**(S(1)/2 - n/S(2))*(S(1) + x/a)**(n/S(2) + S(1)/2)/(x**S(2)*sqrt(S(1) - x**S(2)/a**S(2))), x), x, S(1)/x)


def replacement6560(a, m, n, x):
    return -Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(S(1)/2 - n/S(2))*(S(1) + x/a)**(n/S(2) + S(1)/2)/sqrt(S(1) - x**S(2)/a**S(2)), x), x, S(1)/x)


def replacement6561(a, n, x):
    return -Subst(Int((S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))/x**S(2), x), x, S(1)/x)


def replacement6562(a, m, n, x):
    return -Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2)), x), x, S(1)/x)


def replacement6563(a, m, n, x):
    return -Dist(x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(S(1)/2 - n/S(2))*(S(1) + x/a)**(n/S(2) + S(1)/2)/sqrt(S(1) - x**S(2)/a**S(2)), x), x, S(1)/x), x)


def replacement6564(a, m, n, x):
    return -Dist(x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2)), x), x, S(1)/x), x)


def replacement6565(a, c, d, n, p, x):
    return Simp((c + d*x)**p*(a*x + S(1))*exp(n*acoth(a*x))/(a*(p + S(1))), x)


def replacement6566(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acoth(a*x)), x), x)


def replacement6567(a, c, d, n, p, u, x):
    return Dist(x**(-p)*(c + d*x)**p*(c/(d*x) + S(1))**(-p), Int(u*x**p*(c/(d*x) + S(1))**p*exp(n*acoth(a*x)), x), x)


def replacement6568(a, c, d, n, p, x):
    return -Dist(c**n, Subst(Int((S(1) - x**S(2)/a**S(2))**(n/S(2))*(c + d*x)**(-n + p)/x**S(2), x), x, S(1)/x), x)


def replacement6569(a, c, d, m, n, p, x):
    return -Dist(c**n, Subst(Int(x**(-m + S(-2))*(S(1) - x**S(2)/a**S(2))**(n/S(2))*(c + d*x)**(-n + p), x), x, S(1)/x), x)


def replacement6570(a, c, d, n, p, x):
    return -Dist(c**p, Subst(Int((S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p/x**S(2), x), x, S(1)/x), x)


def replacement6571(a, c, d, m, n, p, x):
    return -Dist(c**p, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p, x), x, S(1)/x), x)


def replacement6572(a, c, d, m, n, p, x):
    return -Dist(c**p*x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2))*(S(1) + x/a)**(n/S(2))*(S(1) + d*x/c)**p, x), x, S(1)/x), x)


def replacement6573(a, c, d, n, p, u, x):
    return Dist((S(1) + d/(c*x))**(-p)*(c + d/x)**p, Int(u*(S(1) + d/(c*x))**p*exp(n*acoth(a*x)), x), x)


def replacement6574(a, c, d, n, x):
    return Simp(exp(n*acoth(a*x))/(a*c*n), x)


def replacement6575(a, c, d, n, x):
    return Simp((-a*x + n)*exp(n*acoth(a*x))/(a*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))), x)


def replacement6576(a, c, d, n, p, x):
    return -Dist(S(2)*(p + S(1))*(S(2)*p + S(3))/(c*(n**S(2) - S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acoth(a*x))/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))), x)


def replacement6577(a, c, d, n, x):
    return -Simp((-a*n*x + S(1))*exp(n*acoth(a*x))/(a**S(2)*c*sqrt(c + d*x**S(2))*(n**S(2) + S(-1))), x)


def replacement6578(a, c, d, n, p, x):
    return -Dist(n*(S(2)*p + S(3))/(a*c*(n**S(2) - S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(a*n*x + S(2)*p + S(2))*exp(n*acoth(a*x))/(a**S(2)*c*(n**S(2) - S(4)*(p + S(1))**S(2))), x)


def replacement6579(a, c, d, n, p, x):
    return -Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acoth(a*x))/(a**S(3)*c*n**S(2)*(n**S(2) + S(-1))), x)


def replacement6580(a, c, d, n, p, x):
    return -Dist((n**S(2) + S(2)*p + S(2))/(a**S(2)*c*(n**S(2) - S(4)*(p + S(1))**S(2))), Int((c + d*x**S(2))**(p + S(1))*exp(n*acoth(a*x)), x), x) + Simp((c + d*x**S(2))**(p + S(1))*(S(2)*a*x*(p + S(1)) + n)*exp(n*acoth(a*x))/(a**S(3)*c*(n**S(2) - S(4)*(p + S(1))**S(2))), x)


def replacement6581(a, c, d, m, n, p, x):
    return -Dist(a**(-m + S(-1))*(-c)**p, Subst(Int((S(1)/tanh(x))**(m + S(2)*p + S(2))*exp(n*x)*cosh(x)**(-S(2)*p + S(-2)), x), x, acoth(a*x)), x)


def replacement6582(a, c, d, n, p, u, x):
    return Dist(d**p, Int(u*x**(S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x), x)


def replacement6583(a, c, d, n, p, u, x):
    return Dist(x**(-S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**(-p)*(c + d*x**S(2))**p, Int(u*x**(S(2)*p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x), x)


def replacement6584(a, c, d, n, p, u, x):
    return Dist(a**(-S(2)*p)*c**p, Int(u*x**(-S(2)*p)*(a*x + S(-1))**(-n/S(2) + p)*(a*x + S(1))**(n/S(2) + p), x), x)


def replacement6585(a, c, d, n, p, x):
    return -Dist(c**p, Subst(Int((S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p)/x**S(2), x), x, S(1)/x), x)


def replacement6586(a, c, d, m, n, p, x):
    return -Dist(c**p, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p), x), x, S(1)/x), x)


def replacement6587(a, c, d, m, n, p, x):
    return -Dist(c**p*x**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(S(1) - x/a)**(-n/S(2) + p)*(S(1) + x/a)**(n/S(2) + p), x), x, S(1)/x), x)


def replacement6588(a, c, d, n, p, u, x):
    return Dist(c**IntPart(p)*(S(1) - S(1)/(a**S(2)*x**S(2)))**(-FracPart(p))*(c + d/x**S(2))**FracPart(p), Int(u*(S(1) - S(1)/(a**S(2)*x**S(2)))**p*exp(n*acoth(a*x)), x), x)


def replacement6589(a, b, c, n, u, x):
    return Dist((S(-1))**(n/S(2)), Int(u*exp(n*atanh(c*(a + b*x))), x), x)


def replacement6590(a, b, c, n, x):
    return Dist((c*(a + b*x))**(n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2))*(a*c + b*c*x + S(1))**(-n/S(2)), Int((a*c + b*c*x + S(-1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x), x)


def replacement6591(a, b, c, m, n, x):
    return Dist(-S(4)*b**(-m + S(-1))*c**(-m + S(-1))/n, Subst(Int(x**(S(2)/n)*(x**(S(2)/n) + S(-1))**(-m + S(-2))*(a*c + x**(S(2)/n)*(-a*c + S(1)) + S(1))**m, x), x, (S(1) - S(1)/(c*(a + b*x)))**(-n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2))), x)


def replacement6592(a, b, c, d, e, m, n, x):
    return Dist((c*(a + b*x))**(n/S(2))*(S(1) + S(1)/(c*(a + b*x)))**(n/S(2))*(a*c + b*c*x + S(1))**(-n/S(2)), Int((d + e*x)**m*(a*c + b*c*x + S(-1))**(-n/S(2))*(a*c + b*c*x + S(1))**(n/S(2)), x), x)


def replacement6593(a, b, c, d, e, n, p, u, x):
    return Dist((c/(S(1) - a**S(2)))**p*((a + b*x + S(1))/(a + b*x))**(n/S(2))*((a + b*x)/(a + b*x + S(1)))**(n/S(2))*(-a - b*x + S(1))**(n/S(2))*(a + b*x + S(-1))**(-n/S(2)), Int(u*(-a - b*x + S(1))**(-n/S(2) + p)*(a + b*x + S(1))**(n/S(2) + p), x), x)


def replacement6594(a, b, c, d, e, n, p, u, x):
    return Dist((c + d*x + e*x**S(2))**p*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**(-p), Int(u*(-a**S(2) - S(2)*a*b*x - b**S(2)*x**S(2) + S(1))**p*exp(n*acoth(a*x)), x), x)


def replacement6595(a, b, c, n, u, x):
    return Int(u*exp(n*atanh(a/c + b*x/c)), x)


def replacement6596(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*atanh(x))**n, x), x, c + d*x), x)


def replacement6597(a, b, c, d, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acoth(x))**n, x), x, c + d*x), x)


def replacement6598(a, b, c, d, n, x):
    return Int((a + b*atanh(c + d*x))**n, x)


def replacement6599(a, b, c, d, n, x):
    return Int((a + b*acoth(c + d*x))**n, x)


def replacement6600(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*atanh(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6601(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/d, Subst(Int((a + b*acoth(x))**n*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6602(a, b, c, d, e, f, m, n, x):
    return Int((a + b*atanh(c + d*x))**n*(e + f*x)**m, x)


def replacement6603(a, b, c, d, e, f, m, n, x):
    return Int((a + b*acoth(c + d*x))**n*(e + f*x)**m, x)


def replacement6604(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*atanh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p, x), x, c + d*x), x)


def replacement6605(A, B, C, a, b, c, d, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acoth(x))**n*(C*x**S(2)/d**S(2) + C/d**S(2))**p, x), x, c + d*x), x)


def replacement6606(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*atanh(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6607(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*acoth(x))**n*(C*x**S(2)/d**S(2) - C/d**S(2))**p*(f*x/d + (-c*f + d*e)/d)**m, x), x, c + d*x), x)


def replacement6608(a, b, c, d, n, x):
    return -Dist(S(1)/2, Int(log(-a - b*x + S(1))/(c + d*x**n), x), x) + Dist(S(1)/2, Int(log(a + b*x + S(1))/(c + d*x**n), x), x)


def replacement6609(a, b, c, d, n, x):
    return -Dist(S(1)/2, Int(log((a + b*x + S(-1))/(a + b*x))/(c + d*x**n), x), x) + Dist(S(1)/2, Int(log((a + b*x + S(1))/(a + b*x))/(c + d*x**n), x), x)


def replacement6610(a, b, c, d, n, x):
    return Int(atanh(a + b*x)/(c + d*x**n), x)


def replacement6611(a, b, c, d, n, x):
    return Int(acoth(a + b*x)/(c + d*x**n), x)


def replacement6612(a, b, n, x):
    return -Dist(b*n, Int(x**n/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x*atanh(a + b*x**n), x)


def replacement6613(a, b, n, x):
    return -Dist(b*n, Int(x**n/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x*acoth(a + b*x**n), x)


def replacement6614(a, b, n, x):
    return -Dist(S(1)/2, Int(log(-a - b*x**n + S(1))/x, x), x) + Dist(S(1)/2, Int(log(a + b*x**n + S(1))/x, x), x)


def replacement6615(a, b, n, x):
    return -Dist(S(1)/2, Int(log(S(1) - S(1)/(a + b*x**n))/x, x), x) + Dist(S(1)/2, Int(log(S(1) + S(1)/(a + b*x**n))/x, x), x)


def replacement6616(a, b, m, n, x):
    return -Dist(b*n/(m + S(1)), Int(x**(m + n)/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x**(m + S(1))*atanh(a + b*x**n)/(m + S(1)), x)


def replacement6617(a, b, m, n, x):
    return -Dist(b*n/(m + S(1)), Int(x**(m + n)/(-a**S(2) - S(2)*a*b*x**n - b**S(2)*x**(S(2)*n) + S(1)), x), x) + Simp(x**(m + S(1))*acoth(a + b*x**n)/(m + S(1)), x)


def replacement6618(a, b, c, d, f, x):
    return -Dist(S(1)/2, Int(log(-a - b*f**(c + d*x) + S(1)), x), x) + Dist(S(1)/2, Int(log(a + b*f**(c + d*x) + S(1)), x), x)


def replacement6619(a, b, c, d, f, x):
    return -Dist(S(1)/2, Int(log(S(1) - S(1)/(a + b*f**(c + d*x))), x), x) + Dist(S(1)/2, Int(log(S(1) + S(1)/(a + b*f**(c + d*x))), x), x)


def replacement6620(a, b, c, d, f, m, x):
    return -Dist(S(1)/2, Int(x**m*log(-a - b*f**(c + d*x) + S(1)), x), x) + Dist(S(1)/2, Int(x**m*log(a + b*f**(c + d*x) + S(1)), x), x)


def replacement6621(a, b, c, d, f, m, x):
    return -Dist(S(1)/2, Int(x**m*log(S(1) - S(1)/(a + b*f**(c + d*x))), x), x) + Dist(S(1)/2, Int(x**m*log(S(1) + S(1)/(a + b*f**(c + d*x))), x), x)


def replacement6622(a, b, c, m, n, u, x):
    return Int(u*acoth(a/c + b*x**n/c)**m, x)


def replacement6623(a, b, c, m, n, u, x):
    return Int(u*atanh(a/c + b*x**n/c)**m, x)


def replacement6624(a, b, c, x):
    return Simp(log(atanh(c*x/sqrt(a + b*x**S(2))))/c, x)


def replacement6625(a, b, c, x):
    return -Simp(log(acoth(c*x/sqrt(a + b*x**S(2))))/c, x)


def replacement6626(a, b, c, m, x):
    return Simp(atanh(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))), x)


def replacement6627(a, b, c, m, x):
    return -Simp(acoth(c*x/sqrt(a + b*x**S(2)))**(m + S(1))/(c*(m + S(1))), x)


def replacement6628(a, b, c, d, e, m, x):
    return Dist(sqrt(a + b*x**S(2))/sqrt(d + e*x**S(2)), Int(atanh(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x), x)


def replacement6629(a, b, c, d, e, m, x):
    return Dist(sqrt(a + b*x**S(2))/sqrt(d + e*x**S(2)), Int(acoth(c*x/sqrt(a + b*x**S(2)))**m/sqrt(a + b*x**S(2)), x), x)


def With6630(a, c, d, n, x):
    u = IntHide((c + d*x**S(2))**n, x)
    return -Dist(a, Int(Dist(S(1)/(-a**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(atanh(a*x), u, x)


def With6631(a, c, d, n, x):
    u = IntHide((c + d*x**S(2))**n, x)
    return -Dist(a, Int(Dist(S(1)/(-a**S(2)*x**S(2) + S(1)), u, x), x), x) + Dist(acoth(a*x), u, x)


def With6632(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        tmp = InverseFunctionOfLinear(u, x)
        res = And(Not(FalseQ(tmp)), SameQ(Head(tmp), ArcTanh), ZeroQ(-D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2)))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6632(n, u, v, x):

    tmp = InverseFunctionOfLinear(u, x)
    return Dist((-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n/Coefficient(Part(tmp, S(1)), x, S(1)), Subst(Int(SimplifyIntegrand((S(1)/cosh(x))**(S(2)*n + S(2))*SubstForInverseFunction(u, tmp, x), x), x), x, tmp), x)


def With6633(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        tmp = InverseFunctionOfLinear(u, x)
        res = And(Not(FalseQ(tmp)), SameQ(Head(tmp), ArcCoth), ZeroQ(-D(v, x)**S(2) + Discriminant(v, x)*Part(tmp, S(1))**S(2)))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6633(n, u, v, x):

    tmp = InverseFunctionOfLinear(u, x)
    return Dist((-Discriminant(v, x)/(S(4)*Coefficient(v, x, S(2))))**n/Coefficient(Part(tmp, S(1)), x, S(1)), Subst(Int(SimplifyIntegrand((-S(1)/sinh(x)**S(2))**(n + S(1))*SubstForInverseFunction(u, tmp, x), x), x), x, tmp), x)


def replacement6634(a, b, c, d, x):
    return Dist(b, Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*atanh(c + d*tanh(a + b*x)), x)


def replacement6635(a, b, c, d, x):
    return Dist(b, Int(x/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*acoth(c + d*tanh(a + b*x)), x)


def replacement6636(a, b, c, d, x):
    return Dist(b, Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*atanh(c + d/tanh(a + b*x)), x)


def replacement6637(a, b, c, d, x):
    return Dist(b, Int(x/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp(x*acoth(c + d/tanh(a + b*x)), x)


def replacement6638(a, b, c, d, x):
    return Dist(b*(-c - d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) - Dist(b*(c + d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp(x*atanh(c + d*tanh(a + b*x)), x)


def replacement6639(a, b, c, d, x):
    return Dist(b*(-c - d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) - Dist(b*(c + d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp(x*acoth(c + d*tanh(a + b*x)), x)


def replacement6640(a, b, c, d, x):
    return -Dist(b*(-c - d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Dist(b*(c + d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp(x*atanh(c + d/tanh(a + b*x)), x)


def replacement6641(a, b, c, d, x):
    return -Dist(b*(-c - d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Dist(b*(c + d + S(1)), Int(x*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp(x*acoth(c + d/tanh(a + b*x)), x)


def replacement6642(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6643(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6644(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6645(a, b, c, d, e, f, m, x):
    return Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*a + S(2)*b*x) + c - d), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6646(a, b, c, d, e, f, m, x):
    return Dist(b*(-c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) - Dist(b*(c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6647(a, b, c, d, e, f, m, x):
    return Dist(b*(-c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d + (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) - Dist(b*(c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d + (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d*tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6648(a, b, c, d, e, f, m, x):
    return -Dist(b*(-c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Dist(b*(c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6649(a, b, c, d, e, f, m, x):
    return -Dist(b*(-c - d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(-c + d - (-c - d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Dist(b*(c + d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*a + S(2)*b*x)/(c - d - (c + d + S(1))*exp(S(2)*a + S(2)*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d/tanh(a + b*x))/(f*(m + S(1))), x)


def replacement6650(a, b, x):
    return -Dist(b, Int(x/cos(S(2)*a + S(2)*b*x), x), x) + Simp(x*atanh(tan(a + b*x)), x)


def replacement6651(a, b, x):
    return -Dist(b, Int(x/cos(S(2)*a + S(2)*b*x), x), x) + Simp(x*acoth(tan(a + b*x)), x)


def replacement6652(a, b, x):
    return -Dist(b, Int(x/cos(S(2)*a + S(2)*b*x), x), x) + Simp(x*atanh(S(1)/tan(a + b*x)), x)


def replacement6653(a, b, x):
    return -Dist(b, Int(x/cos(S(2)*a + S(2)*b*x), x), x) + Simp(x*acoth(S(1)/tan(a + b*x)), x)


def replacement6654(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cos(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*atanh(tan(a + b*x))/(f*(m + S(1))), x)


def replacement6655(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cos(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*acoth(tan(a + b*x))/(f*(m + S(1))), x)


def replacement6656(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cos(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*atanh(S(1)/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6657(a, b, e, f, m, x):
    return -Dist(b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/cos(S(2)*a + S(2)*b*x), x), x) + Simp((e + f*x)**(m + S(1))*acoth(S(1)/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6658(a, b, c, d, x):
    return Dist(I*b, Int(x/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp(x*atanh(c + d*tan(a + b*x)), x)


def replacement6659(a, b, c, d, x):
    return Dist(I*b, Int(x/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp(x*acoth(c + d*tan(a + b*x)), x)


def replacement6660(a, b, c, d, x):
    return Dist(I*b, Int(x/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp(x*atanh(c + d/tan(a + b*x)), x)


def replacement6661(a, b, c, d, x):
    return Dist(I*b, Int(x/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp(x*acoth(c + d/tan(a + b*x)), x)


def replacement6662(a, b, c, d, x):
    return Dist(I*b*(-c + I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-c - I*d + (-c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(I*b*(c - I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(c + I*d + (c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*atanh(c + d*tan(a + b*x)), x)


def replacement6663(a, b, c, d, x):
    return Dist(I*b*(-c + I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-c - I*d + (-c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(I*b*(c - I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(c + I*d + (c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*acoth(c + d*tan(a + b*x)), x)


def replacement6664(a, b, c, d, x):
    return -Dist(I*b*(-c - I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-c + I*d - (-c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(I*b*(c + I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(c - I*d - (c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*atanh(c + d/tan(a + b*x)), x)


def replacement6665(a, b, c, d, x):
    return -Dist(I*b*(-c - I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(-c + I*d - (-c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(I*b*(c + I*d + S(1)), Int(x*exp(S(2)*I*a + S(2)*I*b*x)/(c - I*d - (c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp(x*acoth(c + d/tan(a + b*x)), x)


def replacement6666(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement6667(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(c*exp(S(2)*I*a + S(2)*I*b*x) + c + I*d), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement6668(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6669(a, b, c, d, e, f, m, x):
    return Dist(I*b/(f*(m + S(1))), Int((e + f*x)**(m + S(1))/(-c*exp(S(2)*I*a + S(2)*I*b*x) + c - I*d), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6670(a, b, c, d, e, f, m, x):
    return Dist(I*b*(-c + I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-c - I*d + (-c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(I*b*(c - I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(c + I*d + (c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement6671(a, b, c, d, e, f, m, x):
    return Dist(I*b*(-c + I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-c - I*d + (-c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) - Dist(I*b*(c - I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(c + I*d + (c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d*tan(a + b*x))/(f*(m + S(1))), x)


def replacement6672(a, b, c, d, e, f, m, x):
    return -Dist(I*b*(-c - I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-c + I*d - (-c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(I*b*(c + I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(c - I*d - (c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*atanh(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6673(a, b, c, d, e, f, m, x):
    return -Dist(I*b*(-c - I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(-c + I*d - (-c - I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Dist(I*b*(c + I*d + S(1))/(f*(m + S(1))), Int((e + f*x)**(m + S(1))*exp(S(2)*I*a + S(2)*I*b*x)/(c - I*d - (c + I*d + S(1))*exp(S(2)*I*a + S(2)*I*b*x) + S(1)), x), x) + Simp((e + f*x)**(m + S(1))*acoth(c + d/tan(a + b*x))/(f*(m + S(1))), x)


def replacement6674(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/(S(1) - u**S(2)), x), x) + Simp(x*atanh(u), x)


def replacement6675(u, x):
    return -Int(SimplifyIntegrand(x*D(u, x)/(S(1) - u**S(2)), x), x) + Simp(x*acoth(u), x)


def replacement6676(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(S(1) - u**S(2)), x), x), x) + Simp((a + b*atanh(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement6677(a, b, c, d, m, u, x):
    return -Dist(b/(d*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(S(1) - u**S(2)), x), x), x) + Simp((a + b*acoth(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With6678(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6678(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/(S(1) - u**S(2)), x), x), x) + Dist(a + b*atanh(u), w, x)


def With6679(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6679(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b, Int(SimplifyIntegrand(w*D(u, x)/(S(1) - u**S(2)), x), x), x) + Dist(a + b*acoth(u), w, x)


def replacement6680(c, x):
    return Dist(sqrt(c*x + S(1))*sqrt(S(1)/(c*x + S(1))), Int(S(1)/(sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x) + Simp(x*asech(c*x), x)


def replacement6681(c, x):
    return Dist(S(1)/c, Int(S(1)/(x*sqrt(S(1) + S(1)/(c**S(2)*x**S(2)))), x), x) + Simp(x*acsch(c*x), x)


def replacement6682(a, b, c, n, x):
    return -Dist(S(1)/c, Subst(Int((a + b*x)**n*tanh(x)/cosh(x), x), x, asech(c*x)), x)


def replacement6683(a, b, c, n, x):
    return -Dist(S(1)/c, Subst(Int((a + b*x)**n/(sinh(x)*tanh(x)), x), x, acsch(c*x)), x)


def replacement6684(a, b, c, x):
    return -Subst(Int((a + b*acosh(x/c))/x, x), x, S(1)/x)


def replacement6685(a, b, c, x):
    return -Subst(Int((a + b*asinh(x/c))/x, x), x, S(1)/x)


def replacement6686(a, b, c, m, x):
    return Dist(b*sqrt(c*x + S(1))*sqrt(S(1)/(c*x + S(1)))/(m + S(1)), Int(x**m/(sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x) + Simp(x**(m + S(1))*(a + b*asech(c*x))/(m + S(1)), x)


def replacement6687(a, b, c, m, x):
    return Dist(b/(c*(m + S(1))), Int(x**(m + S(-1))/sqrt(S(1) + S(1)/(c**S(2)*x**S(2))), x), x) + Simp(x**(m + S(1))*(a + b*acsch(c*x))/(m + S(1)), x)


def replacement6688(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(S(1)/cosh(x))**(m + S(1))*tanh(x), x), x, asech(c*x)), x)


def replacement6689(a, b, c, m, n, x):
    return -Dist(c**(-m + S(-1)), Subst(Int((a + b*x)**n*(S(1)/sinh(x))**(m + S(1))/tanh(x), x), x, acsch(c*x)), x)


def replacement6690(a, b, c, m, n, x):
    return Int(x**m*(a + b*asech(c*x))**n, x)


def replacement6691(a, b, c, m, n, x):
    return Int(x**m*(a + b*acsch(c*x))**n, x)


def With6692(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return Dist(b*sqrt(c*x + S(1))*sqrt(S(1)/(c*x + S(1))), Int(SimplifyIntegrand(u/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*asech(c*x), u, x)


def With6693(a, b, c, d, e, p, x):
    u = IntHide((d + e*x**S(2))**p, x)
    return -Dist(b*c*x/sqrt(-c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*acsch(c*x), u, x)


def replacement6694(a, b, c, d, e, n, p, x):
    return -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement6695(a, b, c, d, e, n, p, x):
    return -Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement6696(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6697(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6698(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6699(a, b, c, d, e, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6700(a, b, c, d, e, n, p, x):
    return Int((a + b*asech(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6701(a, b, c, d, e, n, p, x):
    return Int((a + b*acsch(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6702(a, b, c, d, e, p, x):
    return Dist(b*sqrt(c*x + S(1))*sqrt(S(1)/(c*x + S(1)))/(S(2)*e*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x) + Simp((a + b*asech(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def replacement6703(a, b, c, d, e, p, x):
    return -Dist(b*c*x/(S(2)*e*sqrt(-c**S(2)*x**S(2))*(p + S(1))), Int((d + e*x**S(2))**(p + S(1))/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x), x) + Simp((a + b*acsch(c*x))*(d + e*x**S(2))**(p + S(1))/(S(2)*e*(p + S(1))), x)


def With6704(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return Dist(b*sqrt(c*x + S(1))*sqrt(S(1)/(c*x + S(1))), Int(SimplifyIntegrand(u/(x*sqrt(-c*x + S(1))*sqrt(c*x + S(1))), x), x), x) + Dist(a + b*asech(c*x), u, x)


def With6705(a, b, c, d, e, m, p, x):
    u = IntHide(x**m*(d + e*x**S(2))**p, x)
    return -Dist(b*c*x/sqrt(-c**S(2)*x**S(2)), Int(SimplifyIntegrand(u/(x*sqrt(-c**S(2)*x**S(2) + S(-1))), x), x), x) + Dist(a + b*acsch(c*x), u, x)


def replacement6706(a, b, c, d, e, m, n, p, x):
    return -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement6707(a, b, c, d, e, m, n, p, x):
    return -Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x)


def replacement6708(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6709(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(x**S(2))/x, Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6710(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*acosh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6711(a, b, c, d, e, m, n, p, x):
    return -Dist(sqrt(d + e*x**S(2))/(x*sqrt(d/x**S(2) + e)), Subst(Int(x**(-m - S(2)*p + S(-2))*(a + b*asinh(x/c))**n*(d*x**S(2) + e)**p, x), x, S(1)/x), x)


def replacement6712(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*asech(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6713(a, b, c, d, e, m, n, p, x):
    return Int(x**m*(a + b*acsch(c*x))**n*(d + e*x**S(2))**p, x)


def replacement6714(a, b, x):
    return Int(sqrt((-a - b*x + S(1))/(a + b*x + S(1)))/(-a - b*x + S(1)), x) + Simp((a + b*x)*asech(a + b*x)/b, x)


def replacement6715(a, b, x):
    return Int(S(1)/(sqrt(S(1) + (a + b*x)**(S(-2)))*(a + b*x)), x) + Simp((a + b*x)*acsch(a + b*x)/b, x)


def replacement6716(a, b, n, x):
    return -Dist(S(1)/b, Subst(Int(x**n*tanh(x)/cosh(x), x), x, asech(a + b*x)), x)


def replacement6717(a, b, n, x):
    return -Dist(S(1)/b, Subst(Int(x**n/(sinh(x)*tanh(x)), x), x, acsch(a + b*x)), x)


def replacement6718(a, b, x):
    return Simp(log(S(1) - (S(1) - sqrt(S(1) - a**S(2)))*exp(-asech(a + b*x))/a)*asech(a + b*x), x) + Simp(log(S(1) - (sqrt(S(1) - a**S(2)) + S(1))*exp(-asech(a + b*x))/a)*asech(a + b*x), x) - Simp(log(S(1) + exp(-S(2)*asech(a + b*x)))*asech(a + b*x), x) - Simp(PolyLog(S(2), (S(1) - sqrt(S(1) - a**S(2)))*exp(-asech(a + b*x))/a), x) - Simp(PolyLog(S(2), (sqrt(S(1) - a**S(2)) + S(1))*exp(-asech(a + b*x))/a), x) + Simp(PolyLog(S(2), -exp(-S(2)*asech(a + b*x)))/S(2), x)


def replacement6719(a, b, x):
    return Simp(log(S(1) + (S(1) - sqrt(a**S(2) + S(1)))*exp(acsch(a + b*x))/a)*acsch(a + b*x), x) + Simp(log(S(1) + (sqrt(a**S(2) + S(1)) + S(1))*exp(acsch(a + b*x))/a)*acsch(a + b*x), x) - Simp(log(S(1) - exp(-S(2)*acsch(a + b*x)))*acsch(a + b*x), x) + Simp(PolyLog(S(2), -(S(1) - sqrt(a**S(2) + S(1)))*exp(acsch(a + b*x))/a), x) + Simp(PolyLog(S(2), -(sqrt(a**S(2) + S(1)) + S(1))*exp(acsch(a + b*x))/a), x) + Simp(PolyLog(S(2), exp(-S(2)*acsch(a + b*x)))/S(2), x) - Simp(acsch(a + b*x)**S(2), x)


def replacement6720(a, b, m, x):
    return Dist(b**(-m + S(-1))/(m + S(1)), Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/(sqrt(x + S(-1))*sqrt(x + S(1))), x), x, S(1)/(a + b*x)), x) - Simp(b**(-m + S(-1))*(-b**(m + S(1))*x**(m + S(1)) + (-a)**(m + S(1)))*asech(a + b*x)/(m + S(1)), x)


def replacement6721(a, b, m, x):
    return Dist(b**(-m + S(-1))/(m + S(1)), Subst(Int(x**(-m + S(-1))*((-a*x)**(m + S(1)) - (-a*x + S(1))**(m + S(1)))/sqrt(x**S(2) + S(1)), x), x, S(1)/(a + b*x)), x) - Simp(b**(-m + S(-1))*(-b**(m + S(1))*x**(m + S(1)) + (-a)**(m + S(1)))*acsch(a + b*x)/(m + S(1)), x)


def replacement6722(a, b, m, n, x):
    return -Dist(b**(-m + S(-1)), Subst(Int(x**n*(-a + S(1)/cosh(x))**m*tanh(x)/cosh(x), x), x, asech(a + b*x)), x)


def replacement6723(a, b, m, n, x):
    return -Dist(b**(-m + S(-1)), Subst(Int(x**n*(-a + S(1)/sinh(x))**m/(sinh(x)*tanh(x)), x), x, acsch(a + b*x)), x)


def replacement6724(a, b, c, m, n, u, x):
    return Int(u*acosh(a/c + b*x**n/c)**m, x)


def replacement6725(a, b, c, m, n, u, x):
    return Int(u*asinh(a/c + b*x**n/c)**m, x)


def replacement6726(a, x):
    return Dist(S(1)/a, Int(sqrt((-a*x + S(1))/(a*x + S(1)))/(x*(-a*x + S(1))), x), x) + Simp(log(x)/a, x) + Simp(x*exp(asech(a*x)), x)


def replacement6727(a, p, x):
    return Dist(p/a, Int(x**(-p), x), x) + Dist(p*sqrt(a*x**p + S(1))*sqrt(S(1)/(a*x**p + S(1)))/a, Int(x**(-p)/(sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1))), x), x) + Simp(x*exp(asech(a*x**p)), x)


def replacement6728(a, p, x):
    return Dist(S(1)/a, Int(x**(-p), x), x) + Int(sqrt(S(1) + x**(-S(2)*p)/a**S(2)), x)


def replacement6729(n, u, x):
    return Int((sqrt((S(1) - u)/(u + S(1))) + sqrt((S(1) - u)/(u + S(1)))/u + S(1)/u)**n, x)


def replacement6730(n, u, x):
    return Int((sqrt(S(1) + u**(S(-2))) + S(1)/u)**n, x)


def replacement6731(a, p, x):
    return Dist(sqrt(a*x**p + S(1))*sqrt(S(1)/(a*x**p + S(1)))/a, Int(x**(-p + S(-1))*sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1)), x), x) - Simp(x**(-p)/(a*p), x)


def replacement6732(a, m, p, x):
    return Dist(p/(a*(m + S(1))), Int(x**(m - p), x), x) + Dist(p*sqrt(a*x**p + S(1))*sqrt(S(1)/(a*x**p + S(1)))/(a*(m + S(1))), Int(x**(m - p)/(sqrt(-a*x**p + S(1))*sqrt(a*x**p + S(1))), x), x) + Simp(x**(m + S(1))*exp(asech(a*x**p))/(m + S(1)), x)


def replacement6733(a, m, p, x):
    return Dist(S(1)/a, Int(x**(m - p), x), x) + Int(x**m*sqrt(S(1) + x**(-S(2)*p)/a**S(2)), x)


def replacement6734(m, n, u, x):
    return Int(x**m*(sqrt((S(1) - u)/(u + S(1))) + sqrt((S(1) - u)/(u + S(1)))/u + S(1)/u)**n, x)


def replacement6735(m, n, u, x):
    return Int(x**m*(sqrt(S(1) + u**(S(-2))) + S(1)/u)**n, x)


def replacement6736(u, x):
    return Dist(sqrt(S(1) - u**S(2))/(u*sqrt(S(-1) + S(1)/u)*sqrt(S(1) + S(1)/u)), Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(S(1) - u**S(2))), x), x), x) + Simp(x*asech(u), x)


def replacement6737(u, x):
    return -Dist(u/sqrt(-u**S(2)), Int(SimplifyIntegrand(x*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x), x) + Simp(x*acsch(u), x)


def replacement6738(a, b, c, d, m, u, x):
    return Dist(b*sqrt(S(1) - u**S(2))/(d*u*sqrt(S(-1) + S(1)/u)*sqrt(S(1) + S(1)/u)*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(S(1) - u**S(2))), x), x), x) + Simp((a + b*asech(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement6739(a, b, c, d, m, u, x):
    return -Dist(b*u/(d*sqrt(-u**S(2))*(m + S(1))), Int(SimplifyIntegrand((c + d*x)**(m + S(1))*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x), x) + Simp((a + b*acsch(u))*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With6740(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6740(a, b, u, v, x):

    w = IntHide(v, x)
    return Dist(b*sqrt(S(1) - u**S(2))/(u*sqrt(S(-1) + S(1)/u)*sqrt(S(1) + S(1)/u)), Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(S(1) - u**S(2))), x), x), x) + Dist(a + b*asech(u), w, x)


def With6741(a, b, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    w = IntHide(v, x)
    if InverseFunctionFreeQ(w, x):
        return True
    return False


def replacement6741(a, b, u, v, x):

    w = IntHide(v, x)
    return -Dist(b*u/sqrt(-u**S(2)), Int(SimplifyIntegrand(w*D(u, x)/(u*sqrt(-u**S(2) + S(-1))), x), x), x) + Dist(a + b*acsch(u), w, x)
