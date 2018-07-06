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

def special_functions(rubi):
    from sympy.integrals.rubi.constraints import cons67, cons2, cons3, cons66, cons21, cons1256, cons7, cons27, cons17, cons166, cons1726, cons1727, cons94, cons261, cons1728, cons1729, cons62, cons1730, cons1731, cons1732, cons247, cons1733, cons1734, cons1735, cons1736, cons4, cons1247, cons18, cons1351, cons1737, cons1738, cons168, cons1739, cons1740, cons31, cons1741, cons1742, cons1743, cons800, cons87, cons88, cons5, cons50, cons89, cons383, cons48, cons1744, cons1745, cons1746, cons52, cons1747, cons1091, cons125, cons1235, cons13, cons137, cons1371, cons1748, cons1749, cons196, cons1750, cons1751, cons1752, cons150, cons463, cons1753, cons163, cons948, cons949, cons1754, cons1755, cons803, cons1756, cons1757, cons1758, cons1759, cons338, cons1760, cons1761, cons1762, cons1763, cons1764, cons1765, cons38, cons1766, cons347, cons1767, cons1768, cons1769, cons1770, cons1771, cons1772, cons1773

    pattern5012 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5012(a, x, b):
        rubi.append(5012)
        return Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erf(a + b*x)/b, x)
    rule5012 = ReplacementRule(pattern5012, replacement5012)
    pattern5013 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5013(a, x, b):
        rubi.append(5013)
        return -Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfc(a + b*x)/b, x)
    rule5013 = ReplacementRule(pattern5013, replacement5013)
    pattern5014 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5014(a, x, b):
        rubi.append(5014)
        return -Simp(exp((a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfi(a + b*x)/b, x)
    rule5014 = ReplacementRule(pattern5014, replacement5014)
    pattern5015 = Pattern(Integral(Erf(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5015(x, b):
        rubi.append(5015)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -b**S(2)*x**S(2))/sqrt(Pi), x)
    rule5015 = ReplacementRule(pattern5015, replacement5015)
    pattern5016 = Pattern(Integral(Erfc(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5016(x, b):
        rubi.append(5016)
        return -Int(Erf(b*x)/x, x) + Simp(log(x), x)
    rule5016 = ReplacementRule(pattern5016, replacement5016)
    pattern5017 = Pattern(Integral(Erfi(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5017(x, b):
        rubi.append(5017)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), b**S(2)*x**S(2))/sqrt(Pi), x)
    rule5017 = ReplacementRule(pattern5017, replacement5017)
    pattern5018 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5018(x, a, m, b):
        rubi.append(5018)
        return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)/(m + S(1)), x)
    rule5018 = ReplacementRule(pattern5018, replacement5018)
    pattern5019 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5019(x, a, m, b):
        rubi.append(5019)
        return Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)/(m + S(1)), x)
    rule5019 = ReplacementRule(pattern5019, replacement5019)
    pattern5020 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5020(x, a, m, b):
        rubi.append(5020)
        return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp((a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)/(m + S(1)), x)
    rule5020 = ReplacementRule(pattern5020, replacement5020)
    pattern5021 = Pattern(Integral(x_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5021(a, b, x, d, c):
        rubi.append(5021)
        return -Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5021 = ReplacementRule(pattern5021, replacement5021)
    pattern5022 = Pattern(Integral(x_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5022(a, b, x, d, c):
        rubi.append(5022)
        return Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5022 = ReplacementRule(pattern5022, replacement5022)
    pattern5023 = Pattern(Integral(x_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5023(a, b, x, d, c):
        rubi.append(5023)
        return -Dist(b/(sqrt(Pi)*d), Int(exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5023 = ReplacementRule(pattern5023, replacement5023)
    pattern5024 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement5024(a, m, b, x, d, c):
        rubi.append(5024)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5024 = ReplacementRule(pattern5024, replacement5024)
    pattern5025 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement5025(a, m, b, x, d, c):
        rubi.append(5025)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5025 = ReplacementRule(pattern5025, replacement5025)
    pattern5026 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement5026(a, m, b, x, d, c):
        rubi.append(5026)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(-1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule5026 = ReplacementRule(pattern5026, replacement5026)
    pattern5027 = Pattern(Integral(Erf(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1726)
    def replacement5027(x, d, b, c):
        rubi.append(5027)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)
    rule5027 = ReplacementRule(pattern5027, replacement5027)
    pattern5028 = Pattern(Integral(Erfc(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1726)
    def replacement5028(x, d, b, c):
        rubi.append(5028)
        return Int(exp(c + d*x**S(2))/x, x) - Int(Erf(b*x)*exp(c + d*x**S(2))/x, x)
    rule5028 = ReplacementRule(pattern5028, replacement5028)
    pattern5029 = Pattern(Integral(Erfi(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1727)
    def replacement5029(x, d, b, c):
        rubi.append(5029)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)
    rule5029 = ReplacementRule(pattern5029, replacement5029)
    pattern5030 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5030(a, m, b, x, d, c):
        rubi.append(5030)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule5030 = ReplacementRule(pattern5030, replacement5030)
    pattern5031 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5031(a, m, b, x, d, c):
        rubi.append(5031)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule5031 = ReplacementRule(pattern5031, replacement5031)
    pattern5032 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5032(a, m, b, x, d, c):
        rubi.append(5032)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule5032 = ReplacementRule(pattern5032, replacement5032)
    pattern5033 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5033(a, x, b):
        rubi.append(5033)
        return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erf(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erf(a + b*x)**S(2)/b, x)
    rule5033 = ReplacementRule(pattern5033, replacement5033)
    pattern5034 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5034(a, x, b):
        rubi.append(5034)
        return Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfc(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfc(a + b*x)**S(2)/b, x)
    rule5034 = ReplacementRule(pattern5034, replacement5034)
    pattern5035 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5035(a, x, b):
        rubi.append(5035)
        return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfi(a + b*x)*exp((a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfi(a + b*x)**S(2)/b, x)
    rule5035 = ReplacementRule(pattern5035, replacement5035)
    pattern5036 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons261, cons1728)
    def replacement5036(x, m, b):
        rubi.append(5036)
        return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erf(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erf(b*x)**S(2)/(m + S(1)), x)
    rule5036 = ReplacementRule(pattern5036, replacement5036)
    pattern5037 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1729, cons1728)
    def replacement5037(x, m, b):
        rubi.append(5037)
        return Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfc(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(b*x)**S(2)/(m + S(1)), x)
    rule5037 = ReplacementRule(pattern5037, replacement5037)
    pattern5038 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1729, cons1728)
    def replacement5038(x, m, b):
        rubi.append(5038)
        return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfi(b*x)*exp(b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(b*x)**S(2)/(m + S(1)), x)
    rule5038 = ReplacementRule(pattern5038, replacement5038)
    pattern5039 = Pattern(Integral(x_**WC('m', S(1))*Erf(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5039(x, a, m, b):
        rubi.append(5039)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erf(x)**S(2), x), x, a + b*x), x)
    rule5039 = ReplacementRule(pattern5039, replacement5039)
    pattern5040 = Pattern(Integral(x_**WC('m', S(1))*Erfc(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5040(x, a, m, b):
        rubi.append(5040)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfc(x)**S(2), x), x, a + b*x), x)
    rule5040 = ReplacementRule(pattern5040, replacement5040)
    pattern5041 = Pattern(Integral(x_**WC('m', S(1))*Erfi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5041(x, a, m, b):
        rubi.append(5041)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfi(x)**S(2), x), x, a + b*x), x)
    rule5041 = ReplacementRule(pattern5041, replacement5041)
    pattern5042 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5042(a, x, b):
        rubi.append(5042)
        return Simp(cos(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelS(a + b*x)/b, x)
    rule5042 = ReplacementRule(pattern5042, replacement5042)
    pattern5043 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5043(a, x, b):
        rubi.append(5043)
        return -Simp(sin(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelC(a + b*x)/b, x)
    rule5043 = ReplacementRule(pattern5043, replacement5043)
    pattern5044 = Pattern(Integral(FresnelS(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5044(x, b):
        rubi.append(5044)
        return Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)
    rule5044 = ReplacementRule(pattern5044, replacement5044)
    pattern5045 = Pattern(Integral(FresnelC(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5045(x, b):
        rubi.append(5045)
        return Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)
    rule5045 = ReplacementRule(pattern5045, replacement5045)
    pattern5046 = Pattern(Integral(x_**WC('m', S(1))*FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5046(x, a, m, b):
        rubi.append(5046)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(a + b*x)/(m + S(1)), x)
    rule5046 = ReplacementRule(pattern5046, replacement5046)
    pattern5047 = Pattern(Integral(x_**WC('m', S(1))*FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5047(x, a, m, b):
        rubi.append(5047)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(a + b*x)/(m + S(1)), x)
    rule5047 = ReplacementRule(pattern5047, replacement5047)
    pattern5048 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5048(a, x, b):
        rubi.append(5048)
        return -Dist(S(2), Int((a + b*x)*FresnelS(a + b*x)*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelS(a + b*x)**S(2)/b, x)
    rule5048 = ReplacementRule(pattern5048, replacement5048)
    pattern5049 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5049(a, x, b):
        rubi.append(5049)
        return -Dist(S(2), Int((a + b*x)*FresnelC(a + b*x)*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelC(a + b*x)**S(2)/b, x)
    rule5049 = ReplacementRule(pattern5049, replacement5049)
    pattern5050 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1729, cons1730)
    def replacement5050(x, m, b):
        rubi.append(5050)
        return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)**S(2)/(m + S(1)), x)
    rule5050 = ReplacementRule(pattern5050, replacement5050)
    pattern5051 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1729, cons1730)
    def replacement5051(x, m, b):
        rubi.append(5051)
        return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)**S(2)/(m + S(1)), x)
    rule5051 = ReplacementRule(pattern5051, replacement5051)
    pattern5052 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731)
    def replacement5052(x, b, c):
        rubi.append(5052)
        return Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) - Simp(FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5052 = ReplacementRule(pattern5052, replacement5052)
    pattern5053 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731)
    def replacement5053(x, b, c):
        rubi.append(5053)
        return -Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) + Simp(FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5053 = ReplacementRule(pattern5053, replacement5053)
    pattern5054 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons166, cons1732)
    def replacement5054(x, m, b, c):
        rubi.append(5054)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**(m + S(-1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5054 = ReplacementRule(pattern5054, replacement5054)
    pattern5055 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons166, cons1732)
    def replacement5055(x, m, b, c):
        rubi.append(5055)
        return -Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(-1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5055 = ReplacementRule(pattern5055, replacement5055)
    pattern5056 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons247, cons1733)
    def replacement5056(x, m, b, c):
        rubi.append(5056)
        return Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule5056 = ReplacementRule(pattern5056, replacement5056)
    pattern5057 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons247, cons1733)
    def replacement5057(x, m, b, c):
        rubi.append(5057)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule5057 = ReplacementRule(pattern5057, replacement5057)
    pattern5058 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731)
    def replacement5058(x, b, c):
        rubi.append(5058)
        return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) - Simp(x/(S(2)*Pi*b), x) + Simp(FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5058 = ReplacementRule(pattern5058, replacement5058)
    pattern5059 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731)
    def replacement5059(x, b, c):
        rubi.append(5059)
        return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) + Simp(x/(S(2)*Pi*b), x) - Simp(FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5059 = ReplacementRule(pattern5059, replacement5059)
    pattern5060 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons166, cons1734)
    def replacement5060(x, m, b, c):
        rubi.append(5060)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**m/(S(2)*Pi*b*m), x) + Simp(x**(m + S(-1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5060 = ReplacementRule(pattern5060, replacement5060)
    pattern5061 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons166, cons1734)
    def replacement5061(x, m, b, c):
        rubi.append(5061)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**m/(S(2)*Pi*b*m), x) - Simp(x**(m + S(-1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule5061 = ReplacementRule(pattern5061, replacement5061)
    pattern5062 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons94, cons1735)
    def replacement5062(x, m, b, c):
        rubi.append(5062)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule5062 = ReplacementRule(pattern5062, replacement5062)
    pattern5063 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1731, cons17, cons94, cons1735)
    def replacement5063(x, m, b, c):
        rubi.append(5063)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule5063 = ReplacementRule(pattern5063, replacement5063)
    pattern5064 = Pattern(Integral(ExpIntegralE(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1736)
    def replacement5064(x, a, n, b):
        rubi.append(5064)
        return -Simp(ExpIntegralE(n + S(1), a + b*x)/b, x)
    rule5064 = ReplacementRule(pattern5064, replacement5064)
    pattern5065 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1247, cons62)
    def replacement5065(x, n, m, b):
        rubi.append(5065)
        return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), b*x)/b, x)
    rule5065 = ReplacementRule(pattern5065, replacement5065)
    pattern5066 = Pattern(Integral(ExpIntegralE(S(1), x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5066(x, b):
        rubi.append(5066)
        return -Simp(EulerGamma*log(x), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x), x) - Simp(log(b*x)**S(2)/S(2), x)
    rule5066 = ReplacementRule(pattern5066, replacement5066)
    pattern5067 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1247, cons17, cons94)
    def replacement5067(n, x, m, b):
        rubi.append(5067)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + S(1)), x)
    rule5067 = ReplacementRule(pattern5067, replacement5067)
    pattern5068 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons21, cons4, cons1247, cons18)
    def replacement5068(n, x, m, b):
        rubi.append(5068)
        return -Simp(x**(m + S(1))*HypergeometricPFQ(List(m + S(1), m + S(1)), List(m + S(2), m + S(2)), -b*x)/(m + S(1))**S(2), x) + Simp(x**m*(b*x)**(-m)*Gamma(m + S(1))*log(x)/b, x)
    rule5068 = ReplacementRule(pattern5068, replacement5068)
    pattern5069 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons21, cons4, cons1351)
    def replacement5069(x, n, m, b):
        rubi.append(5069)
        return -Simp(x**(m + S(1))*ExpIntegralE(-m, b*x)/(m + n), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + n), x)
    rule5069 = ReplacementRule(pattern5069, replacement5069)
    pattern5070 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons1737)
    def replacement5070(n, a, m, b, x):
        rubi.append(5070)
        return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), a + b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), a + b*x)/b, x)
    rule5070 = ReplacementRule(pattern5070, replacement5070)
    pattern5071 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons1738, cons66)
    def replacement5071(n, a, m, b, x):
        rubi.append(5071)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, a + b*x)/(m + S(1)), x)
    rule5071 = ReplacementRule(pattern5071, replacement5071)
    pattern5072 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5072(a, x, b):
        rubi.append(5072)
        return -Simp(exp(a + b*x)/b, x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)/b, x)
    rule5072 = ReplacementRule(pattern5072, replacement5072)
    pattern5073 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5073(x, a, m, b):
        rubi.append(5073)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*exp(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)/(m + S(1)), x)
    rule5073 = ReplacementRule(pattern5073, replacement5073)
    pattern5074 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5074(a, x, b):
        rubi.append(5074)
        return -Dist(S(2), Int(ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)**S(2)/b, x)
    rule5074 = ReplacementRule(pattern5074, replacement5074)
    pattern5075 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement5075(x, m, b):
        rubi.append(5075)
        return -Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(b*x)*exp(b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(b*x)**S(2)/(m + S(1)), x)
    rule5075 = ReplacementRule(pattern5075, replacement5075)
    pattern5076 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5076(x, a, m, b):
        rubi.append(5076)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*ExpIntegralEi(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*ExpIntegralEi(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule5076 = ReplacementRule(pattern5076, replacement5076)
    pattern5077 = Pattern(Integral(ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5077(a, b, x, d, c):
        rubi.append(5077)
        return -Dist(d/b, Int(exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)
    rule5077 = ReplacementRule(pattern5077, replacement5077)
    pattern5078 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5078(a, m, b, x, d, c):
        rubi.append(5078)
        return -Dist(d/b, Int(x**m*exp(a + c + x*(b + d))/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) + Simp(x**m*ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)
    rule5078 = ReplacementRule(pattern5078, replacement5078)
    pattern5079 = Pattern(Integral(x_**m_*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5079(a, m, b, x, d, c):
        rubi.append(5079)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x)/(m + S(1)), x)
    rule5079 = ReplacementRule(pattern5079, replacement5079)
    pattern5080 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5080(a, x, b):
        rubi.append(5080)
        return -Simp(ExpIntegralEi(S(2)*log(a + b*x))/b, x) + Simp((a + b*x)*LogIntegral(a + b*x)/b, x)
    rule5080 = ReplacementRule(pattern5080, replacement5080)
    pattern5081 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5081(x, b):
        rubi.append(5081)
        return -Simp(b*x, x) + Simp(LogIntegral(b*x)*log(b*x), x)
    rule5081 = ReplacementRule(pattern5081, replacement5081)
    pattern5082 = Pattern(Integral(x_**WC('m', S(1))*LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5082(x, a, m, b):
        rubi.append(5082)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))/log(a + b*x), x), x) + Simp(x**(m + S(1))*LogIntegral(a + b*x)/(m + S(1)), x)
    rule5082 = ReplacementRule(pattern5082, replacement5082)
    pattern5083 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5083(a, x, b):
        rubi.append(5083)
        return Simp(cos(a + b*x)/b, x) + Simp((a + b*x)*SinIntegral(a + b*x)/b, x)
    rule5083 = ReplacementRule(pattern5083, replacement5083)
    pattern5084 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5084(a, x, b):
        rubi.append(5084)
        return -Simp(sin(a + b*x)/b, x) + Simp((a + b*x)*CosIntegral(a + b*x)/b, x)
    rule5084 = ReplacementRule(pattern5084, replacement5084)
    pattern5085 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5085(x, b):
        rubi.append(5085)
        return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x)
    rule5085 = ReplacementRule(pattern5085, replacement5085)
    pattern5086 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5086(x, b):
        rubi.append(5086)
        return Simp(EulerGamma*log(x), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)
    rule5086 = ReplacementRule(pattern5086, replacement5086)
    pattern5087 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5087(x, a, m, b):
        rubi.append(5087)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)/(m + S(1)), x)
    rule5087 = ReplacementRule(pattern5087, replacement5087)
    pattern5088 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5088(x, a, m, b):
        rubi.append(5088)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)/(m + S(1)), x)
    rule5088 = ReplacementRule(pattern5088, replacement5088)
    pattern5089 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5089(a, x, b):
        rubi.append(5089)
        return -Dist(S(2), Int(SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp((a + b*x)*SinIntegral(a + b*x)**S(2)/b, x)
    rule5089 = ReplacementRule(pattern5089, replacement5089)
    pattern5090 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5090(a, x, b):
        rubi.append(5090)
        return -Dist(S(2), Int(CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp((a + b*x)*CosIntegral(a + b*x)**S(2)/b, x)
    rule5090 = ReplacementRule(pattern5090, replacement5090)
    pattern5091 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement5091(x, m, b):
        rubi.append(5091)
        return -Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(b*x)*sin(b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(b*x)**S(2)/(m + S(1)), x)
    rule5091 = ReplacementRule(pattern5091, replacement5091)
    pattern5092 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement5092(x, m, b):
        rubi.append(5092)
        return -Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(b*x)*cos(b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(b*x)**S(2)/(m + S(1)), x)
    rule5092 = ReplacementRule(pattern5092, replacement5092)
    pattern5093 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5093(x, a, m, b):
        rubi.append(5093)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule5093 = ReplacementRule(pattern5093, replacement5093)
    pattern5094 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5094(x, a, m, b):
        rubi.append(5094)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CosIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CosIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule5094 = ReplacementRule(pattern5094, replacement5094)
    pattern5095 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5095(a, b, x, d, c):
        rubi.append(5095)
        return Dist(d/b, Int(sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) - Simp(SinIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule5095 = ReplacementRule(pattern5095, replacement5095)
    pattern5096 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5096(a, b, x, d, c):
        rubi.append(5096)
        return -Dist(d/b, Int(sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(CosIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule5096 = ReplacementRule(pattern5096, replacement5096)
    pattern5097 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5097(a, m, b, x, d, c):
        rubi.append(5097)
        return Dist(d/b, Int(x**m*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*SinIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule5097 = ReplacementRule(pattern5097, replacement5097)
    pattern5098 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5098(a, m, b, x, d, c):
        rubi.append(5098)
        return -Dist(d/b, Int(x**m*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*CosIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule5098 = ReplacementRule(pattern5098, replacement5098)
    pattern5099 = Pattern(Integral(x_**m_*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5099(a, m, b, x, d, c):
        rubi.append(5099)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)
    rule5099 = ReplacementRule(pattern5099, replacement5099)
    pattern5100 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5100(a, m, b, x, d, c):
        rubi.append(5100)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)
    rule5100 = ReplacementRule(pattern5100, replacement5100)
    pattern5101 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5101(a, b, x, d, c):
        rubi.append(5101)
        return -Dist(d/b, Int(sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(SinIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule5101 = ReplacementRule(pattern5101, replacement5101)
    pattern5102 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5102(a, b, x, d, c):
        rubi.append(5102)
        return Dist(d/b, Int(cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Simp(CosIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule5102 = ReplacementRule(pattern5102, replacement5102)
    pattern5103 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5103(a, m, b, x, d, c):
        rubi.append(5103)
        return -Dist(d/b, Int(x**m*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*SinIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule5103 = ReplacementRule(pattern5103, replacement5103)
    pattern5104 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5104(a, m, b, x, d, c):
        rubi.append(5104)
        return Dist(d/b, Int(x**m*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*CosIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule5104 = ReplacementRule(pattern5104, replacement5104)
    pattern5105 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5105(a, m, b, x, d, c):
        rubi.append(5105)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)
    rule5105 = ReplacementRule(pattern5105, replacement5105)
    pattern5106 = Pattern(Integral(x_**m_*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5106(a, m, b, x, d, c):
        rubi.append(5106)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)
    rule5106 = ReplacementRule(pattern5106, replacement5106)
    pattern5107 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5107(a, x, b):
        rubi.append(5107)
        return -Simp(Cosh(a + b*x)/b, x) + Simp((a + b*x)*SinhIntegral(a + b*x)/b, x)
    rule5107 = ReplacementRule(pattern5107, replacement5107)
    pattern5108 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5108(a, x, b):
        rubi.append(5108)
        return -Simp(sinh(a + b*x)/b, x) + Simp((a + b*x)*CoshIntegral(a + b*x)/b, x)
    rule5108 = ReplacementRule(pattern5108, replacement5108)
    pattern5109 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5109(x, b):
        rubi.append(5109)
        return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x)
    rule5109 = ReplacementRule(pattern5109, replacement5109)
    pattern5110 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement5110(x, b):
        rubi.append(5110)
        return Simp(EulerGamma*log(x), x) - Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)
    rule5110 = ReplacementRule(pattern5110, replacement5110)
    pattern5111 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5111(x, a, m, b):
        rubi.append(5111)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)/(m + S(1)), x)
    rule5111 = ReplacementRule(pattern5111, replacement5111)
    pattern5112 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement5112(x, a, m, b):
        rubi.append(5112)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)/(m + S(1)), x)
    rule5112 = ReplacementRule(pattern5112, replacement5112)
    pattern5113 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5113(a, x, b):
        rubi.append(5113)
        return -Dist(S(2), Int(SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp((a + b*x)*SinhIntegral(a + b*x)**S(2)/b, x)
    rule5113 = ReplacementRule(pattern5113, replacement5113)
    pattern5114 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement5114(a, x, b):
        rubi.append(5114)
        return -Dist(S(2), Int(Cosh(a + b*x)*CoshIntegral(a + b*x), x), x) + Simp((a + b*x)*CoshIntegral(a + b*x)**S(2)/b, x)
    rule5114 = ReplacementRule(pattern5114, replacement5114)
    pattern5115 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement5115(x, m, b):
        rubi.append(5115)
        return -Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(b*x)*sinh(b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(b*x)**S(2)/(m + S(1)), x)
    rule5115 = ReplacementRule(pattern5115, replacement5115)
    pattern5116 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement5116(x, m, b):
        rubi.append(5116)
        return -Dist(S(2)/(m + S(1)), Int(x**m*Cosh(b*x)*CoshIntegral(b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(b*x)**S(2)/(m + S(1)), x)
    rule5116 = ReplacementRule(pattern5116, replacement5116)
    pattern5117 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5117(x, a, m, b):
        rubi.append(5117)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinhIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinhIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule5117 = ReplacementRule(pattern5117, replacement5117)
    pattern5118 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement5118(x, a, m, b):
        rubi.append(5118)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CoshIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*Cosh(a + b*x)*CoshIntegral(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CoshIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule5118 = ReplacementRule(pattern5118, replacement5118)
    pattern5119 = Pattern(Integral(SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5119(a, b, x, d, c):
        rubi.append(5119)
        return -Dist(d/b, Int(Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(Cosh(a + b*x)*SinhIntegral(c + d*x)/b, x)
    rule5119 = ReplacementRule(pattern5119, replacement5119)
    pattern5120 = Pattern(Integral(Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5120(a, b, x, d, c):
        rubi.append(5120)
        return -Dist(d/b, Int(Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) + Simp(CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule5120 = ReplacementRule(pattern5120, replacement5120)
    pattern5121 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons168)
    def replacement5121(a, m, b, x, d, c):
        rubi.append(5121)
        return -Dist(d/b, Int(x**m*Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*Cosh(a + b*x)*SinhIntegral(c + d*x), x), x) + Simp(x**m*Cosh(a + b*x)*SinhIntegral(c + d*x)/b, x)
    rule5121 = ReplacementRule(pattern5121, replacement5121)
    pattern5122 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons168)
    def replacement5122(a, m, b, x, d, c):
        rubi.append(5122)
        return -Dist(d/b, Int(x**m*Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule5122 = ReplacementRule(pattern5122, replacement5122)
    pattern5123 = Pattern(Integral(x_**m_*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5123(a, m, b, x, d, c):
        rubi.append(5123)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*SinhIntegral(c + d*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)
    rule5123 = ReplacementRule(pattern5123, replacement5123)
    pattern5124 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5124(a, m, b, x, d, c):
        rubi.append(5124)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*Cosh(a + b*x)*CoshIntegral(c + d*x)/(m + S(1)), x)
    rule5124 = ReplacementRule(pattern5124, replacement5124)
    pattern5125 = Pattern(Integral(Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5125(a, b, x, d, c):
        rubi.append(5125)
        return -Dist(d/b, Int(sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule5125 = ReplacementRule(pattern5125, replacement5125)
    pattern5126 = Pattern(Integral(CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1256)
    def replacement5126(a, b, x, d, c):
        rubi.append(5126)
        return -Dist(d/b, Int(Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) + Simp(Cosh(a + b*x)*CoshIntegral(c + d*x)/b, x)
    rule5126 = ReplacementRule(pattern5126, replacement5126)
    pattern5127 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5127(a, m, b, x, d, c):
        rubi.append(5127)
        return -Dist(d/b, Int(x**m*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule5127 = ReplacementRule(pattern5127, replacement5127)
    pattern5128 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement5128(a, m, b, x, d, c):
        rubi.append(5128)
        return -Dist(d/b, Int(x**m*Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*Cosh(a + b*x)*CoshIntegral(c + d*x), x), x) + Simp(x**m*Cosh(a + b*x)*CoshIntegral(c + d*x)/b, x)
    rule5128 = ReplacementRule(pattern5128, replacement5128)
    pattern5129 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5129(a, m, b, x, d, c):
        rubi.append(5129)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*Cosh(a + b*x)*SinhIntegral(c + d*x)/(m + S(1)), x)
    rule5129 = ReplacementRule(pattern5129, replacement5129)
    pattern5130 = Pattern(Integral(x_**m_*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement5130(a, m, b, x, d, c):
        rubi.append(5130)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*CoshIntegral(c + d*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)
    rule5130 = ReplacementRule(pattern5130, replacement5130)
    pattern5131 = Pattern(Integral(Gamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5131(x, a, n, b):
        rubi.append(5131)
        return -Simp(Gamma(n + S(1), a + b*x)/b, x) + Simp((a + b*x)*Gamma(n, a + b*x)/b, x)
    rule5131 = ReplacementRule(pattern5131, replacement5131)
    pattern5132 = Pattern(Integral(Gamma(n_, b_*x_)/x_, x_), cons3, cons4, cons1739)
    def replacement5132(n, x, b):
        rubi.append(5132)
        return Simp(Gamma(n)*log(x), x) - Simp((b*x)**n*HypergeometricPFQ(List(n, n), List(n + S(1), n + S(1)), -b*x)/n**S(2), x)
    rule5132 = ReplacementRule(pattern5132, replacement5132)
    pattern5133 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, b_*x_), x_), cons3, cons21, cons4, cons66)
    def replacement5133(x, n, m, b):
        rubi.append(5133)
        return Simp(x**(m + S(1))*Gamma(n, b*x)/(m + S(1)), x) - Simp(x**m*(b*x)**(-m)*Gamma(m + n + S(1), b*x)/(b*(m + S(1))), x)
    rule5133 = ReplacementRule(pattern5133, replacement5133)
    def With5134(n, a, m, b, x):
        _UseGamma = True
        rubi.append(5134)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*(a + b*x)**(n + S(-1))*exp(-a - b*x), x), x) + Simp(x**(m + S(1))*Gamma(n, a + b*x)/(m + S(1)), x)
    pattern5134 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons1740, cons66)
    rule5134 = ReplacementRule(pattern5134, With5134)
    pattern5135 = Pattern(Integral(LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5135(a, x, b):
        rubi.append(5135)
        return Simp(PolyGamma(S(-2), a + b*x)/b, x)
    rule5135 = ReplacementRule(pattern5135, replacement5135)
    pattern5136 = Pattern(Integral(x_**WC('m', S(1))*LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons31, cons168)
    def replacement5136(x, a, m, b):
        rubi.append(5136)
        return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(S(-2), a + b*x), x), x) + Simp(x**m*PolyGamma(S(-2), a + b*x)/b, x)
    rule5136 = ReplacementRule(pattern5136, replacement5136)
    pattern5137 = Pattern(Integral(PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1736)
    def replacement5137(x, a, n, b):
        rubi.append(5137)
        return Simp(PolyGamma(n + S(-1), a + b*x)/b, x)
    rule5137 = ReplacementRule(pattern5137, replacement5137)
    pattern5138 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons31, cons168)
    def replacement5138(a, n, m, b, x):
        rubi.append(5138)
        return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(n + S(-1), a + b*x), x), x) + Simp(x**m*PolyGamma(n + S(-1), a + b*x)/b, x)
    rule5138 = ReplacementRule(pattern5138, replacement5138)
    pattern5139 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons31, cons94)
    def replacement5139(a, n, m, b, x):
        rubi.append(5139)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*PolyGamma(n + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*PolyGamma(n, a + b*x)/(m + S(1)), x)
    rule5139 = ReplacementRule(pattern5139, replacement5139)
    pattern5140 = Pattern(Integral(Gamma(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1736)
    def replacement5140(x, a, n, b):
        rubi.append(5140)
        return Simp(Gamma(a + b*x)**n/(b*n), x)
    rule5140 = ReplacementRule(pattern5140, replacement5140)
    pattern5141 = Pattern(Integral(Factorial(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons4, cons1741)
    def replacement5141(a, n, b, x, c):
        rubi.append(5141)
        return Simp(Factorial(a + b*x)**n/(b*n), x)
    rule5141 = ReplacementRule(pattern5141, replacement5141)
    pattern5142 = Pattern(Integral(Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement5142(a, x, b):
        rubi.append(5142)
        return Int(PolyGamma(S(1), a + b*x), x)
    rule5142 = ReplacementRule(pattern5142, replacement5142)
    pattern5143 = Pattern(Integral(Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1742, cons1743)
    def replacement5143(a, x, b, s):
        rubi.append(5143)
        return -Simp(Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)
    rule5143 = ReplacementRule(pattern5143, replacement5143)
    pattern5144 = Pattern(Integral(x_**WC('m', S(1))*Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons31)
    def replacement5144(x, a, m, b):
        rubi.append(5144)
        return Int(x**m*PolyGamma(S(1), a + b*x), x)
    rule5144 = ReplacementRule(pattern5144, replacement5144)
    pattern5145 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1742, cons1743, cons31, cons168)
    def replacement5145(a, m, b, x, s):
        rubi.append(5145)
        return Dist(m/(b*(s + S(-1))), Int(x**(m + S(-1))*Zeta(s + S(-1), a + b*x), x), x) - Simp(x**m*Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)
    rule5145 = ReplacementRule(pattern5145, replacement5145)
    pattern5146 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1742, cons1743, cons31, cons94)
    def replacement5146(a, m, b, x, s):
        rubi.append(5146)
        return Dist(b*s/(m + S(1)), Int(x**(m + S(1))*Zeta(s + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*Zeta(s, a + b*x)/(m + S(1)), x)
    rule5146 = ReplacementRule(pattern5146, replacement5146)
    pattern5147 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons50, cons87, cons88)
    def replacement5147(a, n, p, b, x, q):
        rubi.append(5147)
        return -Dist(p*q, Int(PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n, a*(b*x**p)**q), x)
    rule5147 = ReplacementRule(pattern5147, replacement5147)
    pattern5148 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons50, cons87, cons89)
    def replacement5148(a, n, p, b, x, q):
        rubi.append(5148)
        return -Dist(S(1)/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule5148 = ReplacementRule(pattern5148, replacement5148)
    pattern5149 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons383)
    def replacement5149(e, a, n, p, b, x, d, c):
        rubi.append(5149)
        return Simp(PolyLog(n + S(1), c*(a + b*x)**p)/(e*p), x)
    rule5149 = ReplacementRule(pattern5149, replacement5149)
    pattern5150 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))/x_, x_), cons2, cons3, cons4, cons5, cons50, cons1744)
    def replacement5150(a, n, p, b, x, q):
        rubi.append(5150)
        return Simp(PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule5150 = ReplacementRule(pattern5150, replacement5150)
    pattern5151 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons21, cons5, cons50, cons66, cons87, cons88)
    def replacement5151(a, n, m, p, b, x, q):
        rubi.append(5151)
        return -Dist(p*q/(m + S(1)), Int(x**m*PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n, a*(b*x**p)**q)/(m + S(1)), x)
    rule5151 = ReplacementRule(pattern5151, replacement5151)
    pattern5152 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons21, cons5, cons50, cons66, cons87, cons89)
    def replacement5152(a, n, m, p, b, x, q):
        rubi.append(5152)
        return -Dist((m + S(1))/(p*q), Int(x**m*PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule5152 = ReplacementRule(pattern5152, replacement5152)
    pattern5153 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))*log(x_**WC('m', S(1))*WC('c', S(1)))**WC('r', S(1))/x_, x_), cons2, cons3, cons7, cons21, cons4, cons50, cons52, cons1745, cons1746)
    def replacement5153(a, n, m, p, b, r, x, q, c):
        rubi.append(5153)
        return -Dist(m*r/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**(r + S(-1))/x, x), x) + Simp(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**r/(p*q), x)
    rule5153 = ReplacementRule(pattern5153, replacement5153)
    pattern5154 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons5, cons87, cons88)
    def replacement5154(a, n, p, b, x, c):
        rubi.append(5154)
        return -Dist(p, Int(PolyLog(n + S(-1), c*(a + b*x)**p), x), x) + Dist(a*p, Int(PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x*PolyLog(n, c*(a + b*x)**p), x)
    rule5154 = ReplacementRule(pattern5154, replacement5154)
    pattern5155 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons21, cons5, cons87, cons88, cons62)
    def replacement5155(a, n, m, p, b, x, c):
        rubi.append(5155)
        return -Dist(b*p/(m + S(1)), Int(x**(m + S(1))*PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x**(m + S(1))*PolyLog(n, c*(a + b*x)**p)/(m + S(1)), x)
    rule5155 = ReplacementRule(pattern5155, replacement5155)
    pattern5156 = Pattern(Integral(PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons4, cons5, cons1747)
    def replacement5156(a, n, p, b, F, x, d, c):
        rubi.append(5156)
        return Simp(PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)
    rule5156 = ReplacementRule(pattern5156, replacement5156)
    pattern5157 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons168)
    def replacement5157(e, a, n, m, p, b, F, f, x, d, c):
        rubi.append(5157)
        return -Dist(f*m/(b*c*p*log(F)), Int((e + f*x)**(m + S(-1))*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p), x), x) + Simp((e + f*x)**m*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)
    rule5157 = ReplacementRule(pattern5157, replacement5157)
    def With5158(u, n, x, v):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        w = DerivativeDivides(v, u*v, x)
        if Not(FalseQ(w)):
            return True
        return False
    pattern5158 = Pattern(Integral(u_*PolyLog(n_, v_), x_), cons4, cons4, CustomConstraint(With5158))
    def replacement5158(u, n, x, v):

        w = DerivativeDivides(v, u*v, x)
        rubi.append(5158)
        return Simp(w*PolyLog(n + S(1), v), x)
    rule5158 = ReplacementRule(pattern5158, replacement5158)
    def With5159(n, u, v, w, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = DerivativeDivides(v, u*v, x)
        if Not(FalseQ(z)):
            return True
        return False
    pattern5159 = Pattern(Integral(u_*PolyLog(n_, v_)*log(w_), x_), cons4, cons1235, CustomConstraint(With5159))
    def replacement5159(n, u, v, w, x):

        z = DerivativeDivides(v, u*v, x)
        rubi.append(5159)
        return -Int(SimplifyIntegrand(z*D(w, x)*PolyLog(n + S(1), v)/w, x), x) + Simp(z*PolyLog(n + S(1), v)*log(w), x)
    rule5159 = ReplacementRule(pattern5159, replacement5159)
    pattern5160 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons13, cons137)
    def replacement5160(a, p, b, x, c):
        rubi.append(5160)
        return Dist(p/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*(p + S(1))), x)
    rule5160 = ReplacementRule(pattern5160, replacement5160)
    pattern5161 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons1371)
    def replacement5161(a, p, b, x, c):
        rubi.append(5161)
        return -Dist(p, Int((c*ProductLog(a + b*x))**p/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/b, x)
    rule5161 = ReplacementRule(pattern5161, replacement5161)
    pattern5162 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons5, cons62)
    def replacement5162(a, m, p, b, x, c):
        rubi.append(5162)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p, (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule5162 = ReplacementRule(pattern5162, replacement5162)
    pattern5163 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons4, cons5, cons1748)
    def replacement5163(n, a, p, x, c):
        rubi.append(5163)
        return -Dist(n*p, Int((c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p, x)
    rule5163 = ReplacementRule(pattern5163, replacement5163)
    pattern5164 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons4, cons1749)
    def replacement5164(n, a, p, x, c):
        rubi.append(5164)
        return Dist(n*p/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(n*p + S(1)), x)
    rule5164 = ReplacementRule(pattern5164, replacement5164)
    pattern5165 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons5, cons196)
    def replacement5165(n, a, p, x, c):
        rubi.append(5165)
        return -Subst(Int((c*ProductLog(a*x**(-n)))**p/x**S(2), x), x, S(1)/x)
    rule5165 = ReplacementRule(pattern5165, replacement5165)
    pattern5166 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons4, cons5, cons66, cons1750)
    def replacement5166(n, a, m, p, x, c):
        rubi.append(5166)
        return -Dist(n*p/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + S(1)), x)
    rule5166 = ReplacementRule(pattern5166, replacement5166)
    pattern5167 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons4, cons5, cons1751)
    def replacement5167(n, a, m, p, x, c):
        rubi.append(5167)
        return Dist(n*p/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + n*p + S(1)), x)
    rule5167 = ReplacementRule(pattern5167, replacement5167)
    pattern5168 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons1752)
    def replacement5168(a, m, p, x, c):
        rubi.append(5168)
        return Dist(S(1)/c, Int(x**m*(c*ProductLog(a*x))**(p + S(1))/(ProductLog(a*x) + S(1)), x), x) + Int(x**m*(c*ProductLog(a*x))**p/(ProductLog(a*x) + S(1)), x)
    rule5168 = ReplacementRule(pattern5168, replacement5168)
    pattern5169 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons5, cons150, cons463, cons66)
    def replacement5169(a, n, m, p, x, c):
        rubi.append(5169)
        return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p, x), x, S(1)/x)
    rule5169 = ReplacementRule(pattern5169, replacement5169)
    pattern5170 = Pattern(Integral(S(1)/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons1753)
    def replacement5170(a, x, d, b):
        rubi.append(5170)
        return Simp((a + b*x)/(b*d*ProductLog(a + b*x)), x)
    rule5170 = ReplacementRule(pattern5170, replacement5170)
    pattern5171 = Pattern(Integral(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons1753)
    def replacement5171(a, x, d, b):
        rubi.append(5171)
        return -Int(S(1)/(d*ProductLog(a + b*x) + d), x) + Simp(d*x, x)
    rule5171 = ReplacementRule(pattern5171, replacement5171)
    pattern5172 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons13, cons163)
    def replacement5172(a, p, b, x, d, c):
        rubi.append(5172)
        return -Dist(c*p, Int((c*ProductLog(a + b*x))**(p + S(-1))/(d*ProductLog(a + b*x) + d), x), x) + Simp(c*(c*ProductLog(a + b*x))**(p + S(-1))*(a + b*x)/(b*d), x)
    rule5172 = ReplacementRule(pattern5172, replacement5172)
    pattern5173 = Pattern(Integral(S(1)/((d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))*ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons27, cons1753)
    def replacement5173(a, x, d, b):
        rubi.append(5173)
        return Simp(ExpIntegralEi(ProductLog(a + b*x))/(b*d), x)
    rule5173 = ReplacementRule(pattern5173, replacement5173)
    pattern5174 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons948)
    def replacement5174(a, b, x, d, c):
        rubi.append(5174)
        return Simp(Erfi(sqrt(c*ProductLog(a + b*x))/Rt(c, S(2)))*Rt(Pi*c, S(2))/(b*c*d), x)
    rule5174 = ReplacementRule(pattern5174, replacement5174)
    pattern5175 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons949)
    def replacement5175(a, b, x, d, c):
        rubi.append(5175)
        return Simp(Erf(sqrt(c*ProductLog(a + b*x))/Rt(-c, S(2)))*Rt(-Pi*c, S(2))/(b*c*d), x)
    rule5175 = ReplacementRule(pattern5175, replacement5175)
    pattern5176 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons13, cons137)
    def replacement5176(a, p, b, x, d, c):
        rubi.append(5176)
        return -Dist(S(1)/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(d*ProductLog(a + b*x) + d), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*d*(p + S(1))), x)
    rule5176 = ReplacementRule(pattern5176, replacement5176)
    pattern5177 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons5, cons1754)
    def replacement5177(a, p, b, x, d, c):
        rubi.append(5177)
        return Simp((-ProductLog(a + b*x))**(-p)*(c*ProductLog(a + b*x))**p*Gamma(p + S(1), -ProductLog(a + b*x))/(b*d), x)
    rule5177 = ReplacementRule(pattern5177, replacement5177)
    pattern5178 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons62)
    def replacement5178(a, m, b, x, d):
        rubi.append(5178)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand(S(1)/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule5178 = ReplacementRule(pattern5178, replacement5178)
    pattern5179 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons5, cons62)
    def replacement5179(a, m, p, b, x, d, c):
        rubi.append(5179)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule5179 = ReplacementRule(pattern5179, replacement5179)
    pattern5180 = Pattern(Integral(S(1)/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons196)
    def replacement5180(a, x, d, n):
        rubi.append(5180)
        return -Subst(Int(S(1)/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)
    rule5180 = ReplacementRule(pattern5180, replacement5180)
    pattern5181 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons4, cons5, cons1755)
    def replacement5181(a, n, p, x, d, c):
        rubi.append(5181)
        return Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)
    rule5181 = ReplacementRule(pattern5181, replacement5181)
    pattern5182 = Pattern(Integral(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons803, cons1756)
    def replacement5182(a, n, p, x, d):
        rubi.append(5182)
        return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)
    rule5182 = ReplacementRule(pattern5182, replacement5182)
    pattern5183 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons803, cons1757, cons1758)
    def replacement5183(a, n, p, x, d, c):
        rubi.append(5183)
        return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(c*n, S(2)))*Rt(Pi*c*n, S(2))/(d*n), x)
    rule5183 = ReplacementRule(pattern5183, replacement5183)
    pattern5184 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons803, cons1757, cons1759)
    def replacement5184(a, n, p, x, d, c):
        rubi.append(5184)
        return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(-c*n, S(2)))*Rt(-Pi*c*n, S(2))/(d*n), x)
    rule5184 = ReplacementRule(pattern5184, replacement5184)
    pattern5185 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons338, cons88, cons1760)
    def replacement5185(a, n, p, x, d, c):
        rubi.append(5185)
        return -Dist(c*(n*(p + S(-1)) + S(1)), Int((c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)
    rule5185 = ReplacementRule(pattern5185, replacement5185)
    pattern5186 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons338, cons88, cons1761)
    def replacement5186(a, n, p, x, d, c):
        rubi.append(5186)
        return -Dist(S(1)/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(d*(n*p + S(1))), x)
    rule5186 = ReplacementRule(pattern5186, replacement5186)
    pattern5187 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons5, cons196)
    def replacement5187(a, n, p, x, d, c):
        rubi.append(5187)
        return -Subst(Int((c*ProductLog(a*x**(-n)))**p/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)
    rule5187 = ReplacementRule(pattern5187, replacement5187)
    pattern5188 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons31, cons168)
    def replacement5188(x, a, m, d):
        rubi.append(5188)
        return -Dist(m/(m + S(1)), Int(x**m/((d*ProductLog(a*x) + d)*ProductLog(a*x)), x), x) + Simp(x**(m + S(1))/(d*(m + S(1))*ProductLog(a*x)), x)
    rule5188 = ReplacementRule(pattern5188, replacement5188)
    pattern5189 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons27, cons1762)
    def replacement5189(a, x, d):
        rubi.append(5189)
        return Simp(log(ProductLog(a*x))/d, x)
    rule5189 = ReplacementRule(pattern5189, replacement5189)
    pattern5190 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons31, cons94)
    def replacement5190(x, a, m, d):
        rubi.append(5190)
        return -Int(x**m*ProductLog(a*x)/(d*ProductLog(a*x) + d), x) + Simp(x**(m + S(1))/(d*(m + S(1))), x)
    rule5190 = ReplacementRule(pattern5190, replacement5190)
    pattern5191 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons21, cons18)
    def replacement5191(x, a, m, d):
        rubi.append(5191)
        return Simp(x**m*(-(m + S(1))*ProductLog(a*x))**(-m)*Gamma(m + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)
    rule5191 = ReplacementRule(pattern5191, replacement5191)
    pattern5192 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons27, cons4, cons1763)
    def replacement5192(x, a, n, d):
        rubi.append(5192)
        return Simp(log(ProductLog(a*x**n))/(d*n), x)
    rule5192 = ReplacementRule(pattern5192, replacement5192)
    pattern5193 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons150, cons463, cons66)
    def replacement5193(a, n, m, x, d):
        rubi.append(5193)
        return -Subst(Int(x**(-m + S(-2))/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)
    rule5193 = ReplacementRule(pattern5193, replacement5193)
    pattern5194 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons7, cons27, cons4, cons5, cons1764)
    def replacement5194(a, n, p, x, d, c):
        rubi.append(5194)
        return Simp((c*ProductLog(a*x**n))**p/(d*n*p), x)
    rule5194 = ReplacementRule(pattern5194, replacement5194)
    pattern5195 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1765)
    def replacement5195(n, a, m, p, x, d, c):
        rubi.append(5195)
        return Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)
    rule5195 = ReplacementRule(pattern5195, replacement5195)
    pattern5196 = Pattern(Integral(x_**WC('m', S(1))*ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons21, cons4, cons38, cons1766)
    def replacement5196(a, n, m, p, x, d):
        rubi.append(5196)
        return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)
    rule5196 = ReplacementRule(pattern5196, replacement5196)
    pattern5197 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons66, cons347, cons1767, cons1768)
    def replacement5197(n, a, m, p, x, d, c):
        rubi.append(5197)
        return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(c/(p + S(-1)/2), S(2)))*Rt(Pi*c/(p + S(-1)/2), S(2))/(d*n), x)
    rule5197 = ReplacementRule(pattern5197, replacement5197)
    pattern5198 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons66, cons347, cons1767, cons1769)
    def replacement5198(n, a, m, p, x, d, c):
        rubi.append(5198)
        return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(-c/(p + S(-1)/2), S(2)))*Rt(-Pi*c/(p + S(-1)/2), S(2))/(d*n), x)
    rule5198 = ReplacementRule(pattern5198, replacement5198)
    pattern5199 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1770, cons1771)
    def replacement5199(n, a, m, p, x, d, c):
        rubi.append(5199)
        return -Dist(c*(m + n*(p + S(-1)) + S(1))/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)
    rule5199 = ReplacementRule(pattern5199, replacement5199)
    pattern5200 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1770, cons1772)
    def replacement5200(n, a, m, p, x, d, c):
        rubi.append(5200)
        return -Dist((m + S(1))/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(d*(m + n*p + S(1))), x)
    rule5200 = ReplacementRule(pattern5200, replacement5200)
    pattern5201 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons5, cons66)
    def replacement5201(a, m, p, x, d, c):
        rubi.append(5201)
        return Simp(x**m*(c*ProductLog(a*x))**p*(-(m + S(1))*ProductLog(a*x))**(-m - p)*Gamma(m + p + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)
    rule5201 = ReplacementRule(pattern5201, replacement5201)
    pattern5202 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons5, cons66, cons150, cons463)
    def replacement5202(n, a, m, p, x, d, c):
        rubi.append(5202)
        return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)
    rule5202 = ReplacementRule(pattern5202, replacement5202)
    pattern5203 = Pattern(Integral(u_, x_), cons1773)
    def replacement5203(u, x):
        rubi.append(5203)
        return Subst(Int(SimplifyIntegrand((x + S(1))*SubstFor(ProductLog(x), u, x)*exp(x), x), x), x, ProductLog(x))
    rule5203 = ReplacementRule(pattern5203, replacement5203)
    return [rule5012, rule5013, rule5014, rule5015, rule5016, rule5017, rule5018, rule5019, rule5020, rule5021, rule5022, rule5023, rule5024, rule5025, rule5026, rule5027, rule5028, rule5029, rule5030, rule5031, rule5032, rule5033, rule5034, rule5035, rule5036, rule5037, rule5038, rule5039, rule5040, rule5041, rule5042, rule5043, rule5044, rule5045, rule5046, rule5047, rule5048, rule5049, rule5050, rule5051, rule5052, rule5053, rule5054, rule5055, rule5056, rule5057, rule5058, rule5059, rule5060, rule5061, rule5062, rule5063, rule5064, rule5065, rule5066, rule5067, rule5068, rule5069, rule5070, rule5071, rule5072, rule5073, rule5074, rule5075, rule5076, rule5077, rule5078, rule5079, rule5080, rule5081, rule5082, rule5083, rule5084, rule5085, rule5086, rule5087, rule5088, rule5089, rule5090, rule5091, rule5092, rule5093, rule5094, rule5095, rule5096, rule5097, rule5098, rule5099, rule5100, rule5101, rule5102, rule5103, rule5104, rule5105, rule5106, rule5107, rule5108, rule5109, rule5110, rule5111, rule5112, rule5113, rule5114, rule5115, rule5116, rule5117, rule5118, rule5119, rule5120, rule5121, rule5122, rule5123, rule5124, rule5125, rule5126, rule5127, rule5128, rule5129, rule5130, rule5131, rule5132, rule5133, rule5134, rule5135, rule5136, rule5137, rule5138, rule5139, rule5140, rule5141, rule5142, rule5143, rule5144, rule5145, rule5146, rule5147, rule5148, rule5149, rule5150, rule5151, rule5152, rule5153, rule5154, rule5155, rule5156, rule5157, rule5158, rule5159, rule5160, rule5161, rule5162, rule5163, rule5164, rule5165, rule5166, rule5167, rule5168, rule5169, rule5170, rule5171, rule5172, rule5173, rule5174, rule5175, rule5176, rule5177, rule5178, rule5179, rule5180, rule5181, rule5182, rule5183, rule5184, rule5185, rule5186, rule5187, rule5188, rule5189, rule5190, rule5191, rule5192, rule5193, rule5194, rule5195, rule5196, rule5197, rule5198, rule5199, rule5200, rule5201, rule5202, rule5203, ]
