from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
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
        Zeta, ProductLog, HypergeometricPFQ, exp, log, DerivativeDivides
    )
    from sympy import Integral, S, sqrt, And, Or, Integer, Float, Mod, I, Abs 
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
    from sympy.integrals.rubi.constraints import cons67, cons2, cons3, cons66, cons21, cons1245, cons7, cons27, cons17, cons166, cons1246, cons1247, cons94, cons261, cons1248, cons1249, cons62, cons1250, cons1251, cons1252, cons247, cons1253, cons1254, cons1255, cons1256, cons4, cons1257, cons18, cons1258, cons1259, cons1260, cons168, cons1261, cons1262, cons31, cons1263, cons1264, cons1265, cons800, cons87, cons88, cons5, cons50, cons89, cons383, cons48, cons1266, cons1267, cons1268, cons52, cons1269, cons1091, cons125, cons1239, cons13, cons137, cons1270, cons1271, cons1272, cons196, cons1273, cons1274, cons1275, cons150, cons463, cons1276, cons163, cons948, cons949, cons1277, cons1278, cons803, cons1279, cons1280, cons1281, cons1282, cons338, cons1283, cons1284, cons1285, cons1286, cons1287, cons1288, cons38, cons1289, cons347, cons1290, cons1291, cons1292, cons1293, cons1294, cons1295, cons1296

    pattern2145 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2145(x, a, b):
        rubi.append(2145)
        return Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erf(a + b*x)/b, x)
    rule2145 = ReplacementRule(pattern2145, replacement2145)
    pattern2146 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2146(x, a, b):
        rubi.append(2146)
        return -Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfc(a + b*x)/b, x)
    rule2146 = ReplacementRule(pattern2146, replacement2146)
    pattern2147 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2147(x, a, b):
        rubi.append(2147)
        return -Simp(exp((a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfi(a + b*x)/b, x)
    rule2147 = ReplacementRule(pattern2147, replacement2147)
    pattern2148 = Pattern(Integral(Erf(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2148(x, b):
        rubi.append(2148)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -b**S(2)*x**S(2))/sqrt(Pi), x)
    rule2148 = ReplacementRule(pattern2148, replacement2148)
    pattern2149 = Pattern(Integral(Erfc(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2149(x, b):
        rubi.append(2149)
        return -Int(Erf(b*x)/x, x) + Simp(log(x), x)
    rule2149 = ReplacementRule(pattern2149, replacement2149)
    pattern2150 = Pattern(Integral(Erfi(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2150(x, b):
        rubi.append(2150)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), b**S(2)*x**S(2))/sqrt(Pi), x)
    rule2150 = ReplacementRule(pattern2150, replacement2150)
    pattern2151 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2151(x, a, b, m):
        rubi.append(2151)
        return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)/(m + S(1)), x)
    rule2151 = ReplacementRule(pattern2151, replacement2151)
    pattern2152 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2152(x, a, b, m):
        rubi.append(2152)
        return Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)/(m + S(1)), x)
    rule2152 = ReplacementRule(pattern2152, replacement2152)
    pattern2153 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2153(x, a, b, m):
        rubi.append(2153)
        return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp((a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)/(m + S(1)), x)
    rule2153 = ReplacementRule(pattern2153, replacement2153)
    pattern2154 = Pattern(Integral(x_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2154(a, x, d, c, b):
        rubi.append(2154)
        return -Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2154 = ReplacementRule(pattern2154, replacement2154)
    pattern2155 = Pattern(Integral(x_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2155(a, x, d, c, b):
        rubi.append(2155)
        return Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2155 = ReplacementRule(pattern2155, replacement2155)
    pattern2156 = Pattern(Integral(x_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2156(a, x, d, c, b):
        rubi.append(2156)
        return -Dist(b/(sqrt(Pi)*d), Int(exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2156 = ReplacementRule(pattern2156, replacement2156)
    pattern2157 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement2157(a, x, d, m, c, b):
        rubi.append(2157)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2157 = ReplacementRule(pattern2157, replacement2157)
    pattern2158 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement2158(a, x, d, m, c, b):
        rubi.append(2158)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2158 = ReplacementRule(pattern2158, replacement2158)
    pattern2159 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons166)
    def replacement2159(a, x, d, m, c, b):
        rubi.append(2159)
        return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(-1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)
    rule2159 = ReplacementRule(pattern2159, replacement2159)
    pattern2160 = Pattern(Integral(Erf(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1246)
    def replacement2160(x, d, c, b):
        rubi.append(2160)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)
    rule2160 = ReplacementRule(pattern2160, replacement2160)
    pattern2161 = Pattern(Integral(Erfc(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1246)
    def replacement2161(x, d, c, b):
        rubi.append(2161)
        return Int(exp(c + d*x**S(2))/x, x) - Int(Erf(b*x)*exp(c + d*x**S(2))/x, x)
    rule2161 = ReplacementRule(pattern2161, replacement2161)
    pattern2162 = Pattern(Integral(Erfi(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1247)
    def replacement2162(x, d, c, b):
        rubi.append(2162)
        return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)
    rule2162 = ReplacementRule(pattern2162, replacement2162)
    pattern2163 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2163(a, x, d, m, c, b):
        rubi.append(2163)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule2163 = ReplacementRule(pattern2163, replacement2163)
    pattern2164 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2164(a, x, d, m, c, b):
        rubi.append(2164)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule2164 = ReplacementRule(pattern2164, replacement2164)
    pattern2165 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2165(a, x, d, m, c, b):
        rubi.append(2165)
        return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)
    rule2165 = ReplacementRule(pattern2165, replacement2165)
    pattern2166 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2166(x, a, b):
        rubi.append(2166)
        return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erf(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erf(a + b*x)**S(2)/b, x)
    rule2166 = ReplacementRule(pattern2166, replacement2166)
    pattern2167 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2167(x, a, b):
        rubi.append(2167)
        return Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfc(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfc(a + b*x)**S(2)/b, x)
    rule2167 = ReplacementRule(pattern2167, replacement2167)
    pattern2168 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2168(x, a, b):
        rubi.append(2168)
        return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfi(a + b*x)*exp((a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfi(a + b*x)**S(2)/b, x)
    rule2168 = ReplacementRule(pattern2168, replacement2168)
    pattern2169 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons261, cons1248)
    def replacement2169(x, b, m):
        rubi.append(2169)
        return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erf(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erf(b*x)**S(2)/(m + S(1)), x)
    rule2169 = ReplacementRule(pattern2169, replacement2169)
    pattern2170 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1249, cons1248)
    def replacement2170(x, b, m):
        rubi.append(2170)
        return Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfc(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(b*x)**S(2)/(m + S(1)), x)
    rule2170 = ReplacementRule(pattern2170, replacement2170)
    pattern2171 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1249, cons1248)
    def replacement2171(x, b, m):
        rubi.append(2171)
        return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfi(b*x)*exp(b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(b*x)**S(2)/(m + S(1)), x)
    rule2171 = ReplacementRule(pattern2171, replacement2171)
    pattern2172 = Pattern(Integral(x_**WC('m', S(1))*Erf(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2172(x, a, b, m):
        rubi.append(2172)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erf(x)**S(2), x), x, a + b*x), x)
    rule2172 = ReplacementRule(pattern2172, replacement2172)
    pattern2173 = Pattern(Integral(x_**WC('m', S(1))*Erfc(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2173(x, a, b, m):
        rubi.append(2173)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfc(x)**S(2), x), x, a + b*x), x)
    rule2173 = ReplacementRule(pattern2173, replacement2173)
    pattern2174 = Pattern(Integral(x_**WC('m', S(1))*Erfi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2174(x, a, b, m):
        rubi.append(2174)
        return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfi(x)**S(2), x), x, a + b*x), x)
    rule2174 = ReplacementRule(pattern2174, replacement2174)
    pattern2175 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2175(x, a, b):
        rubi.append(2175)
        return Simp(cos(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelS(a + b*x)/b, x)
    rule2175 = ReplacementRule(pattern2175, replacement2175)
    pattern2176 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2176(x, a, b):
        rubi.append(2176)
        return -Simp(sin(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelC(a + b*x)/b, x)
    rule2176 = ReplacementRule(pattern2176, replacement2176)
    pattern2177 = Pattern(Integral(FresnelS(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2177(x, b):
        rubi.append(2177)
        return Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)
    rule2177 = ReplacementRule(pattern2177, replacement2177)
    pattern2178 = Pattern(Integral(FresnelC(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2178(x, b):
        rubi.append(2178)
        return Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)
    rule2178 = ReplacementRule(pattern2178, replacement2178)
    pattern2179 = Pattern(Integral(x_**WC('m', S(1))*FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2179(x, a, b, m):
        rubi.append(2179)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(a + b*x)/(m + S(1)), x)
    rule2179 = ReplacementRule(pattern2179, replacement2179)
    pattern2180 = Pattern(Integral(x_**WC('m', S(1))*FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2180(x, a, b, m):
        rubi.append(2180)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(a + b*x)/(m + S(1)), x)
    rule2180 = ReplacementRule(pattern2180, replacement2180)
    pattern2181 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2181(x, a, b):
        rubi.append(2181)
        return -Dist(S(2), Int((a + b*x)*FresnelS(a + b*x)*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelS(a + b*x)**S(2)/b, x)
    rule2181 = ReplacementRule(pattern2181, replacement2181)
    pattern2182 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2182(x, a, b):
        rubi.append(2182)
        return -Dist(S(2), Int((a + b*x)*FresnelC(a + b*x)*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelC(a + b*x)**S(2)/b, x)
    rule2182 = ReplacementRule(pattern2182, replacement2182)
    pattern2183 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1249, cons1250)
    def replacement2183(x, b, m):
        rubi.append(2183)
        return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)**S(2)/(m + S(1)), x)
    rule2183 = ReplacementRule(pattern2183, replacement2183)
    pattern2184 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))**S(2), x_), cons3, cons17, cons1249, cons1250)
    def replacement2184(x, b, m):
        rubi.append(2184)
        return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)**S(2)/(m + S(1)), x)
    rule2184 = ReplacementRule(pattern2184, replacement2184)
    pattern2185 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251)
    def replacement2185(x, c, b):
        rubi.append(2185)
        return Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) - Simp(FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2185 = ReplacementRule(pattern2185, replacement2185)
    pattern2186 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251)
    def replacement2186(x, c, b):
        rubi.append(2186)
        return -Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) + Simp(FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2186 = ReplacementRule(pattern2186, replacement2186)
    pattern2187 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons166, cons1252)
    def replacement2187(x, c, b, m):
        rubi.append(2187)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**(m + S(-1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2187 = ReplacementRule(pattern2187, replacement2187)
    pattern2188 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons166, cons1252)
    def replacement2188(x, c, m, b):
        rubi.append(2188)
        return -Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(-1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2188 = ReplacementRule(pattern2188, replacement2188)
    pattern2189 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons247, cons1253)
    def replacement2189(x, c, b, m):
        rubi.append(2189)
        return Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule2189 = ReplacementRule(pattern2189, replacement2189)
    pattern2190 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons247, cons1253)
    def replacement2190(x, c, m, b):
        rubi.append(2190)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule2190 = ReplacementRule(pattern2190, replacement2190)
    pattern2191 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251)
    def replacement2191(x, c, b):
        rubi.append(2191)
        return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) - Simp(x/(S(2)*Pi*b), x) + Simp(FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2191 = ReplacementRule(pattern2191, replacement2191)
    pattern2192 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251)
    def replacement2192(x, c, b):
        rubi.append(2192)
        return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) + Simp(x/(S(2)*Pi*b), x) - Simp(FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2192 = ReplacementRule(pattern2192, replacement2192)
    pattern2193 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons166, cons1254)
    def replacement2193(x, c, m, b):
        rubi.append(2193)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**m/(S(2)*Pi*b*m), x) + Simp(x**(m + S(-1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2193 = ReplacementRule(pattern2193, replacement2193)
    pattern2194 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons166, cons1254)
    def replacement2194(x, c, b, m):
        rubi.append(2194)
        return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**m/(S(2)*Pi*b*m), x) - Simp(x**(m + S(-1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)
    rule2194 = ReplacementRule(pattern2194, replacement2194)
    pattern2195 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons94, cons1255)
    def replacement2195(x, c, m, b):
        rubi.append(2195)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule2195 = ReplacementRule(pattern2195, replacement2195)
    pattern2196 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons7, cons1251, cons17, cons94, cons1255)
    def replacement2196(x, c, b, m):
        rubi.append(2196)
        return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)
    rule2196 = ReplacementRule(pattern2196, replacement2196)
    pattern2197 = Pattern(Integral(ExpIntegralE(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1256)
    def replacement2197(x, a, b, n):
        rubi.append(2197)
        return -Simp(ExpIntegralE(n + S(1), a + b*x)/b, x)
    rule2197 = ReplacementRule(pattern2197, replacement2197)
    pattern2198 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1257, cons62)
    def replacement2198(x, b, m, n):
        rubi.append(2198)
        return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), b*x)/b, x)
    rule2198 = ReplacementRule(pattern2198, replacement2198)
    pattern2199 = Pattern(Integral(ExpIntegralE(S(1), x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2199(x, b):
        rubi.append(2199)
        return -Simp(EulerGamma*log(x), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x), x) - Simp(log(b*x)**S(2)/S(2), x)
    rule2199 = ReplacementRule(pattern2199, replacement2199)
    pattern2200 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1257, cons17, cons94)
    def replacement2200(x, b, m, n):
        rubi.append(2200)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + S(1)), x)
    rule2200 = ReplacementRule(pattern2200, replacement2200)
    pattern2201 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons21, cons4, cons1257, cons18)
    def replacement2201(x, b, m, n):
        rubi.append(2201)
        return -Simp(x**(m + S(1))*HypergeometricPFQ(List(m + S(1), m + S(1)), List(m + S(2), m + S(2)), -b*x)/(m + S(1))**S(2), x) + Simp(x**m*(b*x)**(-m)*Gamma(m + S(1))*log(x)/b, x)
    rule2201 = ReplacementRule(pattern2201, replacement2201)
    pattern2202 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons21, cons4, cons1258)
    def replacement2202(x, b, m, n):
        rubi.append(2202)
        return -Simp(x**(m + S(1))*ExpIntegralE(-m, b*x)/(m + n), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + n), x)
    rule2202 = ReplacementRule(pattern2202, replacement2202)
    pattern2203 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons1259)
    def replacement2203(a, n, x, m, b):
        rubi.append(2203)
        return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), a + b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), a + b*x)/b, x)
    rule2203 = ReplacementRule(pattern2203, replacement2203)
    pattern2204 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons1260, cons66)
    def replacement2204(a, n, x, m, b):
        rubi.append(2204)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, a + b*x)/(m + S(1)), x)
    rule2204 = ReplacementRule(pattern2204, replacement2204)
    pattern2205 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2205(x, a, b):
        rubi.append(2205)
        return -Simp(exp(a + b*x)/b, x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)/b, x)
    rule2205 = ReplacementRule(pattern2205, replacement2205)
    pattern2206 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2206(x, a, b, m):
        rubi.append(2206)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*exp(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)/(m + S(1)), x)
    rule2206 = ReplacementRule(pattern2206, replacement2206)
    pattern2207 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2207(x, a, b):
        rubi.append(2207)
        return -Dist(S(2), Int(ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)**S(2)/b, x)
    rule2207 = ReplacementRule(pattern2207, replacement2207)
    pattern2208 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement2208(x, b, m):
        rubi.append(2208)
        return -Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(b*x)*exp(b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(b*x)**S(2)/(m + S(1)), x)
    rule2208 = ReplacementRule(pattern2208, replacement2208)
    pattern2209 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2209(x, a, b, m):
        rubi.append(2209)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*ExpIntegralEi(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*ExpIntegralEi(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule2209 = ReplacementRule(pattern2209, replacement2209)
    pattern2210 = Pattern(Integral(ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2210(a, x, d, c, b):
        rubi.append(2210)
        return -Dist(d/b, Int(exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)
    rule2210 = ReplacementRule(pattern2210, replacement2210)
    pattern2211 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2211(a, x, d, m, c, b):
        rubi.append(2211)
        return -Dist(d/b, Int(x**m*exp(a + c + x*(b + d))/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) + Simp(x**m*ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)
    rule2211 = ReplacementRule(pattern2211, replacement2211)
    pattern2212 = Pattern(Integral(x_**m_*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2212(a, x, d, m, c, b):
        rubi.append(2212)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x)/(m + S(1)), x)
    rule2212 = ReplacementRule(pattern2212, replacement2212)
    pattern2213 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2213(x, a, b):
        rubi.append(2213)
        return -Simp(ExpIntegralEi(S(2)*log(a + b*x))/b, x) + Simp((a + b*x)*LogIntegral(a + b*x)/b, x)
    rule2213 = ReplacementRule(pattern2213, replacement2213)
    pattern2214 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2214(x, b):
        rubi.append(2214)
        return -Simp(b*x, x) + Simp(LogIntegral(b*x)*log(b*x), x)
    rule2214 = ReplacementRule(pattern2214, replacement2214)
    pattern2215 = Pattern(Integral(x_**WC('m', S(1))*LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2215(x, a, b, m):
        rubi.append(2215)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))/log(a + b*x), x), x) + Simp(x**(m + S(1))*LogIntegral(a + b*x)/(m + S(1)), x)
    rule2215 = ReplacementRule(pattern2215, replacement2215)
    pattern2216 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2216(x, a, b):
        rubi.append(2216)
        return Simp(cos(a + b*x)/b, x) + Simp((a + b*x)*SinIntegral(a + b*x)/b, x)
    rule2216 = ReplacementRule(pattern2216, replacement2216)
    pattern2217 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2217(x, a, b):
        rubi.append(2217)
        return -Simp(sin(a + b*x)/b, x) + Simp((a + b*x)*CosIntegral(a + b*x)/b, x)
    rule2217 = ReplacementRule(pattern2217, replacement2217)
    pattern2218 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2218(x, b):
        rubi.append(2218)
        return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x)
    rule2218 = ReplacementRule(pattern2218, replacement2218)
    pattern2219 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2219(x, b):
        rubi.append(2219)
        return Simp(EulerGamma*log(x), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)
    rule2219 = ReplacementRule(pattern2219, replacement2219)
    pattern2220 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2220(x, a, b, m):
        rubi.append(2220)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)/(m + S(1)), x)
    rule2220 = ReplacementRule(pattern2220, replacement2220)
    pattern2221 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2221(x, a, b, m):
        rubi.append(2221)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)/(m + S(1)), x)
    rule2221 = ReplacementRule(pattern2221, replacement2221)
    pattern2222 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2222(x, a, b):
        rubi.append(2222)
        return -Dist(S(2), Int(SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp((a + b*x)*SinIntegral(a + b*x)**S(2)/b, x)
    rule2222 = ReplacementRule(pattern2222, replacement2222)
    pattern2223 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2223(x, a, b):
        rubi.append(2223)
        return -Dist(S(2), Int(CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp((a + b*x)*CosIntegral(a + b*x)**S(2)/b, x)
    rule2223 = ReplacementRule(pattern2223, replacement2223)
    pattern2224 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement2224(x, b, m):
        rubi.append(2224)
        return -Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(b*x)*sin(b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(b*x)**S(2)/(m + S(1)), x)
    rule2224 = ReplacementRule(pattern2224, replacement2224)
    pattern2225 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement2225(x, b, m):
        rubi.append(2225)
        return -Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(b*x)*cos(b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(b*x)**S(2)/(m + S(1)), x)
    rule2225 = ReplacementRule(pattern2225, replacement2225)
    pattern2226 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2226(x, a, b, m):
        rubi.append(2226)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule2226 = ReplacementRule(pattern2226, replacement2226)
    pattern2227 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2227(x, a, b, m):
        rubi.append(2227)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CosIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CosIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule2227 = ReplacementRule(pattern2227, replacement2227)
    pattern2228 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2228(a, x, d, c, b):
        rubi.append(2228)
        return Dist(d/b, Int(sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) - Simp(SinIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule2228 = ReplacementRule(pattern2228, replacement2228)
    pattern2229 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2229(a, x, d, c, b):
        rubi.append(2229)
        return -Dist(d/b, Int(sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(CosIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule2229 = ReplacementRule(pattern2229, replacement2229)
    pattern2230 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2230(a, x, d, m, c, b):
        rubi.append(2230)
        return Dist(d/b, Int(x**m*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*SinIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule2230 = ReplacementRule(pattern2230, replacement2230)
    pattern2231 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2231(a, x, d, m, c, b):
        rubi.append(2231)
        return -Dist(d/b, Int(x**m*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*CosIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule2231 = ReplacementRule(pattern2231, replacement2231)
    pattern2232 = Pattern(Integral(x_**m_*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2232(a, x, d, m, c, b):
        rubi.append(2232)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)
    rule2232 = ReplacementRule(pattern2232, replacement2232)
    pattern2233 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2233(a, x, d, m, c, b):
        rubi.append(2233)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)
    rule2233 = ReplacementRule(pattern2233, replacement2233)
    pattern2234 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2234(a, x, d, c, b):
        rubi.append(2234)
        return -Dist(d/b, Int(sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(SinIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule2234 = ReplacementRule(pattern2234, replacement2234)
    pattern2235 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2235(a, x, d, c, b):
        rubi.append(2235)
        return Dist(d/b, Int(cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Simp(CosIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule2235 = ReplacementRule(pattern2235, replacement2235)
    pattern2236 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2236(a, x, d, m, c, b):
        rubi.append(2236)
        return -Dist(d/b, Int(x**m*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*SinIntegral(c + d*x)*sin(a + b*x)/b, x)
    rule2236 = ReplacementRule(pattern2236, replacement2236)
    pattern2237 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2237(a, x, d, m, c, b):
        rubi.append(2237)
        return Dist(d/b, Int(x**m*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*CosIntegral(c + d*x)*cos(a + b*x)/b, x)
    rule2237 = ReplacementRule(pattern2237, replacement2237)
    pattern2238 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2238(a, x, d, m, c, b):
        rubi.append(2238)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)
    rule2238 = ReplacementRule(pattern2238, replacement2238)
    pattern2239 = Pattern(Integral(x_**m_*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2239(a, x, d, m, c, b):
        rubi.append(2239)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)
    rule2239 = ReplacementRule(pattern2239, replacement2239)
    pattern2240 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2240(x, a, b):
        rubi.append(2240)
        return -Simp(Cosh(a + b*x)/b, x) + Simp((a + b*x)*SinhIntegral(a + b*x)/b, x)
    rule2240 = ReplacementRule(pattern2240, replacement2240)
    pattern2241 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2241(x, a, b):
        rubi.append(2241)
        return -Simp(sinh(a + b*x)/b, x) + Simp((a + b*x)*CoshIntegral(a + b*x)/b, x)
    rule2241 = ReplacementRule(pattern2241, replacement2241)
    pattern2242 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2242(x, b):
        rubi.append(2242)
        return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x)
    rule2242 = ReplacementRule(pattern2242, replacement2242)
    pattern2243 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    def replacement2243(x, b):
        rubi.append(2243)
        return Simp(EulerGamma*log(x), x) - Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)
    rule2243 = ReplacementRule(pattern2243, replacement2243)
    pattern2244 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2244(x, a, b, m):
        rubi.append(2244)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)/(m + S(1)), x)
    rule2244 = ReplacementRule(pattern2244, replacement2244)
    pattern2245 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons66)
    def replacement2245(x, a, b, m):
        rubi.append(2245)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)/(m + S(1)), x)
    rule2245 = ReplacementRule(pattern2245, replacement2245)
    pattern2246 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2246(x, a, b):
        rubi.append(2246)
        return -Dist(S(2), Int(SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp((a + b*x)*SinhIntegral(a + b*x)**S(2)/b, x)
    rule2246 = ReplacementRule(pattern2246, replacement2246)
    pattern2247 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons67)
    def replacement2247(x, a, b):
        rubi.append(2247)
        return -Dist(S(2), Int(Cosh(a + b*x)*CoshIntegral(a + b*x), x), x) + Simp((a + b*x)*CoshIntegral(a + b*x)**S(2)/b, x)
    rule2247 = ReplacementRule(pattern2247, replacement2247)
    pattern2248 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement2248(x, b, m):
        rubi.append(2248)
        return -Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(b*x)*sinh(b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(b*x)**S(2)/(m + S(1)), x)
    rule2248 = ReplacementRule(pattern2248, replacement2248)
    pattern2249 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons62)
    def replacement2249(x, b, m):
        rubi.append(2249)
        return -Dist(S(2)/(m + S(1)), Int(x**m*Cosh(b*x)*CoshIntegral(b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(b*x)**S(2)/(m + S(1)), x)
    rule2249 = ReplacementRule(pattern2249, replacement2249)
    pattern2250 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2250(x, a, b, m):
        rubi.append(2250)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinhIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinhIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule2250 = ReplacementRule(pattern2250, replacement2250)
    pattern2251 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons62)
    def replacement2251(x, a, b, m):
        rubi.append(2251)
        return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CoshIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*Cosh(a + b*x)*CoshIntegral(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CoshIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)
    rule2251 = ReplacementRule(pattern2251, replacement2251)
    pattern2252 = Pattern(Integral(SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2252(a, x, d, c, b):
        rubi.append(2252)
        return -Dist(d/b, Int(Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(Cosh(a + b*x)*SinhIntegral(c + d*x)/b, x)
    rule2252 = ReplacementRule(pattern2252, replacement2252)
    pattern2253 = Pattern(Integral(Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2253(a, x, d, c, b):
        rubi.append(2253)
        return -Dist(d/b, Int(Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) + Simp(CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule2253 = ReplacementRule(pattern2253, replacement2253)
    pattern2254 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons168)
    def replacement2254(a, x, d, m, c, b):
        rubi.append(2254)
        return -Dist(d/b, Int(x**m*Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*Cosh(a + b*x)*SinhIntegral(c + d*x), x), x) + Simp(x**m*Cosh(a + b*x)*SinhIntegral(c + d*x)/b, x)
    rule2254 = ReplacementRule(pattern2254, replacement2254)
    pattern2255 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons168)
    def replacement2255(a, x, d, m, c, b):
        rubi.append(2255)
        return -Dist(d/b, Int(x**m*Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule2255 = ReplacementRule(pattern2255, replacement2255)
    pattern2256 = Pattern(Integral(x_**m_*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2256(a, x, d, m, c, b):
        rubi.append(2256)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*SinhIntegral(c + d*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)
    rule2256 = ReplacementRule(pattern2256, replacement2256)
    pattern2257 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2257(a, x, d, m, c, b):
        rubi.append(2257)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*Cosh(a + b*x)*CoshIntegral(c + d*x)/(m + S(1)), x)
    rule2257 = ReplacementRule(pattern2257, replacement2257)
    pattern2258 = Pattern(Integral(Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2258(a, x, d, c, b):
        rubi.append(2258)
        return -Dist(d/b, Int(sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule2258 = ReplacementRule(pattern2258, replacement2258)
    pattern2259 = Pattern(Integral(CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons1245)
    def replacement2259(a, x, d, c, b):
        rubi.append(2259)
        return -Dist(d/b, Int(Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) + Simp(Cosh(a + b*x)*CoshIntegral(c + d*x)/b, x)
    rule2259 = ReplacementRule(pattern2259, replacement2259)
    pattern2260 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2260(a, x, d, m, c, b):
        rubi.append(2260)
        return -Dist(d/b, Int(x**m*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)
    rule2260 = ReplacementRule(pattern2260, replacement2260)
    pattern2261 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons62)
    def replacement2261(a, x, d, m, c, b):
        rubi.append(2261)
        return -Dist(d/b, Int(x**m*Cosh(a + b*x)*Cosh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*Cosh(a + b*x)*CoshIntegral(c + d*x), x), x) + Simp(x**m*Cosh(a + b*x)*CoshIntegral(c + d*x)/b, x)
    rule2261 = ReplacementRule(pattern2261, replacement2261)
    pattern2262 = Pattern(Integral(x_**WC('m', S(1))*Cosh(x_*WC('b', S(1)) + WC('a', S(0)))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2262(a, x, d, m, c, b):
        rubi.append(2262)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*Cosh(a + b*x)*SinhIntegral(c + d*x)/(m + S(1)), x)
    rule2262 = ReplacementRule(pattern2262, replacement2262)
    pattern2263 = Pattern(Integral(x_**m_*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons17, cons94)
    def replacement2263(a, x, d, m, c, b):
        rubi.append(2263)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*Cosh(a + b*x)*CoshIntegral(c + d*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*Cosh(c + d*x)*sinh(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)
    rule2263 = ReplacementRule(pattern2263, replacement2263)
    pattern2264 = Pattern(Integral(Gamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2264(x, a, b, n):
        rubi.append(2264)
        return -Simp(Gamma(n + S(1), a + b*x)/b, x) + Simp((a + b*x)*Gamma(n, a + b*x)/b, x)
    rule2264 = ReplacementRule(pattern2264, replacement2264)
    pattern2265 = Pattern(Integral(Gamma(n_, b_*x_)/x_, x_), cons3, cons4, cons1261)
    def replacement2265(x, n, b):
        rubi.append(2265)
        return Simp(Gamma(n)*log(x), x) - Simp((b*x)**n*HypergeometricPFQ(List(n, n), List(n + S(1), n + S(1)), -b*x)/n**S(2), x)
    rule2265 = ReplacementRule(pattern2265, replacement2265)
    pattern2266 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, b_*x_), x_), cons3, cons21, cons4, cons66)
    def replacement2266(x, n, m, b):
        rubi.append(2266)
        return Simp(x**(m + S(1))*Gamma(n, b*x)/(m + S(1)), x) - Simp(x**m*(b*x)**(-m)*Gamma(m + n + S(1), b*x)/(b*(m + S(1))), x)
    rule2266 = ReplacementRule(pattern2266, replacement2266)
    def With2267(a, n, x, m, b):
        _UseGamma = True
        rubi.append(2267)
        return Dist(b/(m + S(1)), Int(x**(m + S(1))*(a + b*x)**(n + S(-1))*exp(-a - b*x), x), x) + Simp(x**(m + S(1))*Gamma(n, a + b*x)/(m + S(1)), x)
    pattern2267 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons21, cons4, cons1262, cons66)
    rule2267 = ReplacementRule(pattern2267, With2267)
    pattern2268 = Pattern(Integral(LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2268(x, a, b):
        rubi.append(2268)
        return Simp(PolyGamma(S(-2), a + b*x)/b, x)
    rule2268 = ReplacementRule(pattern2268, replacement2268)
    pattern2269 = Pattern(Integral(x_**WC('m', S(1))*LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons31, cons168)
    def replacement2269(x, a, b, m):
        rubi.append(2269)
        return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(S(-2), a + b*x), x), x) + Simp(x**m*PolyGamma(S(-2), a + b*x)/b, x)
    rule2269 = ReplacementRule(pattern2269, replacement2269)
    pattern2270 = Pattern(Integral(PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1256)
    def replacement2270(x, a, b, n):
        rubi.append(2270)
        return Simp(PolyGamma(n + S(-1), a + b*x)/b, x)
    rule2270 = ReplacementRule(pattern2270, replacement2270)
    pattern2271 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons31, cons168)
    def replacement2271(a, n, x, m, b):
        rubi.append(2271)
        return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(n + S(-1), a + b*x), x), x) + Simp(x**m*PolyGamma(n + S(-1), a + b*x)/b, x)
    rule2271 = ReplacementRule(pattern2271, replacement2271)
    pattern2272 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons31, cons94)
    def replacement2272(a, n, x, m, b):
        rubi.append(2272)
        return -Dist(b/(m + S(1)), Int(x**(m + S(1))*PolyGamma(n + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*PolyGamma(n, a + b*x)/(m + S(1)), x)
    rule2272 = ReplacementRule(pattern2272, replacement2272)
    pattern2273 = Pattern(Integral(Gamma(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1256)
    def replacement2273(x, a, b, n):
        rubi.append(2273)
        return Simp(Gamma(a + b*x)**n/(b*n), x)
    rule2273 = ReplacementRule(pattern2273, replacement2273)
    pattern2274 = Pattern(Integral(Factorial(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons7, cons4, cons1263)
    def replacement2274(a, n, x, c, b):
        rubi.append(2274)
        return Simp(Factorial(a + b*x)**n/(b*n), x)
    rule2274 = ReplacementRule(pattern2274, replacement2274)
    pattern2275 = Pattern(Integral(Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons67)
    def replacement2275(x, a, b):
        rubi.append(2275)
        return Int(PolyGamma(S(1), a + b*x), x)
    rule2275 = ReplacementRule(pattern2275, replacement2275)
    pattern2276 = Pattern(Integral(Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1264, cons1265)
    def replacement2276(x, a, b, s):
        rubi.append(2276)
        return -Simp(Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)
    rule2276 = ReplacementRule(pattern2276, replacement2276)
    pattern2277 = Pattern(Integral(x_**WC('m', S(1))*Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons31)
    def replacement2277(x, a, b, m):
        rubi.append(2277)
        return Int(x**m*PolyGamma(S(1), a + b*x), x)
    rule2277 = ReplacementRule(pattern2277, replacement2277)
    pattern2278 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1264, cons1265, cons31, cons168)
    def replacement2278(a, x, m, b, s):
        rubi.append(2278)
        return Dist(m/(b*(s + S(-1))), Int(x**(m + S(-1))*Zeta(s + S(-1), a + b*x), x), x) - Simp(x**m*Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)
    rule2278 = ReplacementRule(pattern2278, replacement2278)
    pattern2279 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons800, cons1264, cons1265, cons31, cons94)
    def replacement2279(a, x, m, b, s):
        rubi.append(2279)
        return Dist(b*s/(m + S(1)), Int(x**(m + S(1))*Zeta(s + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*Zeta(s, a + b*x)/(m + S(1)), x)
    rule2279 = ReplacementRule(pattern2279, replacement2279)
    pattern2280 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons50, cons87, cons88)
    def replacement2280(a, n, p, x, q, b):
        rubi.append(2280)
        return -Dist(p*q, Int(PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n, a*(b*x**p)**q), x)
    rule2280 = ReplacementRule(pattern2280, replacement2280)
    pattern2281 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons50, cons87, cons89)
    def replacement2281(a, n, p, x, q, b):
        rubi.append(2281)
        return -Dist(S(1)/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule2281 = ReplacementRule(pattern2281, replacement2281)
    pattern2282 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons383)
    def replacement2282(a, n, p, x, d, c, b, e):
        rubi.append(2282)
        return Simp(PolyLog(n + S(1), c*(a + b*x)**p)/(e*p), x)
    rule2282 = ReplacementRule(pattern2282, replacement2282)
    pattern2283 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))/x_, x_), cons2, cons3, cons4, cons5, cons50, cons1266)
    def replacement2283(a, n, p, x, q, b):
        rubi.append(2283)
        return Simp(PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule2283 = ReplacementRule(pattern2283, replacement2283)
    pattern2284 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons21, cons5, cons50, cons66, cons87, cons88)
    def replacement2284(a, n, p, x, m, b, q):
        rubi.append(2284)
        return -Dist(p*q/(m + S(1)), Int(x**m*PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n, a*(b*x**p)**q)/(m + S(1)), x)
    rule2284 = ReplacementRule(pattern2284, replacement2284)
    pattern2285 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons21, cons5, cons50, cons66, cons87, cons89)
    def replacement2285(a, n, p, x, m, b, q):
        rubi.append(2285)
        return -Dist((m + S(1))/(p*q), Int(x**m*PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)
    rule2285 = ReplacementRule(pattern2285, replacement2285)
    pattern2286 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))*log(x_**WC('m', S(1))*WC('c', S(1)))**WC('r', S(1))/x_, x_), cons2, cons3, cons7, cons21, cons4, cons50, cons52, cons1267, cons1268)
    def replacement2286(a, r, n, p, x, m, q, c, b):
        rubi.append(2286)
        return -Dist(m*r/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**(r + S(-1))/x, x), x) + Simp(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**r/(p*q), x)
    rule2286 = ReplacementRule(pattern2286, replacement2286)
    pattern2287 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons5, cons87, cons88)
    def replacement2287(a, n, p, x, c, b):
        rubi.append(2287)
        return -Dist(p, Int(PolyLog(n + S(-1), c*(a + b*x)**p), x), x) + Dist(a*p, Int(PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x*PolyLog(n, c*(a + b*x)**p), x)
    rule2287 = ReplacementRule(pattern2287, replacement2287)
    pattern2288 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons7, cons21, cons5, cons87, cons88, cons62)
    def replacement2288(a, n, p, x, m, c, b):
        rubi.append(2288)
        return -Dist(b*p/(m + S(1)), Int(x**(m + S(1))*PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x**(m + S(1))*PolyLog(n, c*(a + b*x)**p)/(m + S(1)), x)
    rule2288 = ReplacementRule(pattern2288, replacement2288)
    pattern2289 = Pattern(Integral(PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons4, cons5, cons1269)
    def replacement2289(a, n, p, x, d, F, c, b):
        rubi.append(2289)
        return Simp(PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)
    rule2289 = ReplacementRule(pattern2289, replacement2289)
    pattern2290 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons168)
    def replacement2290(f, a, n, p, x, d, m, F, c, b, e):
        rubi.append(2290)
        return -Dist(f*m/(b*c*p*log(F)), Int((e + f*x)**(m + S(-1))*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p), x), x) + Simp((e + f*x)**m*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)
    rule2290 = ReplacementRule(pattern2290, replacement2290)
    def With2291(x, n, v, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        w = DerivativeDivides(v, u*v, x)
        if Not(FalseQ(w)):
            return True
        return False
    pattern2291 = Pattern(Integral(u_*PolyLog(n_, v_), x_), cons4, cons4, CustomConstraint(With2291))
    def replacement2291(x, n, v, u):
        
        w = DerivativeDivides(v, u*v, x)
        rubi.append(2291)
        return Simp(w*PolyLog(n + S(1), v), x)
    rule2291 = ReplacementRule(pattern2291, replacement2291)
    def With2292(n, x, w, v, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = DerivativeDivides(v, u*v, x)
        if Not(FalseQ(z)):
            return True
        return False
    pattern2292 = Pattern(Integral(u_*PolyLog(n_, v_)*log(w_), x_), cons4, cons1239, CustomConstraint(With2292))
    def replacement2292(n, x, w, v, u):
        
        z = DerivativeDivides(v, u*v, x)
        rubi.append(2292)
        return -Int(SimplifyIntegrand(z*D(w, x)*PolyLog(n + S(1), v)/w, x), x) + Simp(z*PolyLog(n + S(1), v)*log(w), x)
    rule2292 = ReplacementRule(pattern2292, replacement2292)
    pattern2293 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons13, cons137)
    def replacement2293(a, p, x, c, b):
        rubi.append(2293)
        return Dist(p/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*(p + S(1))), x)
    rule2293 = ReplacementRule(pattern2293, replacement2293)
    pattern2294 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons1270)
    def replacement2294(a, p, x, c, b):
        rubi.append(2294)
        return -Dist(p, Int((c*ProductLog(a + b*x))**p/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/b, x)
    rule2294 = ReplacementRule(pattern2294, replacement2294)
    pattern2295 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons5, cons62)
    def replacement2295(a, p, x, m, c, b):
        rubi.append(2295)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p, (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule2295 = ReplacementRule(pattern2295, replacement2295)
    pattern2296 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons4, cons5, cons1271)
    def replacement2296(a, n, p, x, c):
        rubi.append(2296)
        return -Dist(n*p, Int((c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p, x)
    rule2296 = ReplacementRule(pattern2296, replacement2296)
    pattern2297 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons4, cons1272)
    def replacement2297(a, n, p, x, c):
        rubi.append(2297)
        return Dist(n*p/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(n*p + S(1)), x)
    rule2297 = ReplacementRule(pattern2297, replacement2297)
    pattern2298 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons5, cons196)
    def replacement2298(a, n, p, x, c):
        rubi.append(2298)
        return -Subst(Int((c*ProductLog(a*x**(-n)))**p/x**S(2), x), x, S(1)/x)
    rule2298 = ReplacementRule(pattern2298, replacement2298)
    pattern2299 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons4, cons5, cons66, cons1273)
    def replacement2299(a, n, p, x, m, c):
        rubi.append(2299)
        return -Dist(n*p/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + S(1)), x)
    rule2299 = ReplacementRule(pattern2299, replacement2299)
    pattern2300 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons4, cons5, cons1274)
    def replacement2300(a, n, p, x, m, c):
        rubi.append(2300)
        return Dist(n*p/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + n*p + S(1)), x)
    rule2300 = ReplacementRule(pattern2300, replacement2300)
    pattern2301 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons21, cons1275)
    def replacement2301(a, p, x, m, c):
        rubi.append(2301)
        return Dist(S(1)/c, Int(x**m*(c*ProductLog(a*x))**(p + S(1))/(ProductLog(a*x) + S(1)), x), x) + Int(x**m*(c*ProductLog(a*x))**p/(ProductLog(a*x) + S(1)), x)
    rule2301 = ReplacementRule(pattern2301, replacement2301)
    pattern2302 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons7, cons5, cons150, cons463, cons66)
    def replacement2302(a, n, p, x, m, c):
        rubi.append(2302)
        return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p, x), x, S(1)/x)
    rule2302 = ReplacementRule(pattern2302, replacement2302)
    pattern2303 = Pattern(Integral(S(1)/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons1276)
    def replacement2303(x, a, d, b):
        rubi.append(2303)
        return Simp((a + b*x)/(b*d*ProductLog(a + b*x)), x)
    rule2303 = ReplacementRule(pattern2303, replacement2303)
    pattern2304 = Pattern(Integral(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons1276)
    def replacement2304(x, d, a, b):
        rubi.append(2304)
        return -Int(S(1)/(d*ProductLog(a + b*x) + d), x) + Simp(d*x, x)
    rule2304 = ReplacementRule(pattern2304, replacement2304)
    pattern2305 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons13, cons163)
    def replacement2305(a, p, x, d, c, b):
        rubi.append(2305)
        return -Dist(c*p, Int((c*ProductLog(a + b*x))**(p + S(-1))/(d*ProductLog(a + b*x) + d), x), x) + Simp(c*(c*ProductLog(a + b*x))**(p + S(-1))*(a + b*x)/(b*d), x)
    rule2305 = ReplacementRule(pattern2305, replacement2305)
    pattern2306 = Pattern(Integral(S(1)/((d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))*ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons27, cons1276)
    def replacement2306(x, d, a, b):
        rubi.append(2306)
        return Simp(ExpIntegralEi(ProductLog(a + b*x))/(b*d), x)
    rule2306 = ReplacementRule(pattern2306, replacement2306)
    pattern2307 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons948)
    def replacement2307(a, x, d, c, b):
        rubi.append(2307)
        return Simp(Erfi(sqrt(c*ProductLog(a + b*x))/Rt(c, S(2)))*Rt(Pi*c, S(2))/(b*c*d), x)
    rule2307 = ReplacementRule(pattern2307, replacement2307)
    pattern2308 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons949)
    def replacement2308(a, x, d, c, b):
        rubi.append(2308)
        return Simp(Erf(sqrt(c*ProductLog(a + b*x))/Rt(-c, S(2)))*Rt(-Pi*c, S(2))/(b*c*d), x)
    rule2308 = ReplacementRule(pattern2308, replacement2308)
    pattern2309 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons13, cons137)
    def replacement2309(a, p, x, d, c, b):
        rubi.append(2309)
        return -Dist(S(1)/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(d*ProductLog(a + b*x) + d), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*d*(p + S(1))), x)
    rule2309 = ReplacementRule(pattern2309, replacement2309)
    pattern2310 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons5, cons1277)
    def replacement2310(a, p, x, d, c, b):
        rubi.append(2310)
        return Simp((-ProductLog(a + b*x))**(-p)*(c*ProductLog(a + b*x))**p*Gamma(p + S(1), -ProductLog(a + b*x))/(b*d), x)
    rule2310 = ReplacementRule(pattern2310, replacement2310)
    pattern2311 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons27, cons62)
    def replacement2311(a, x, d, m, b):
        rubi.append(2311)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand(S(1)/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule2311 = ReplacementRule(pattern2311, replacement2311)
    pattern2312 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons7, cons27, cons5, cons62)
    def replacement2312(a, p, x, d, m, c, b):
        rubi.append(2312)
        return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)
    rule2312 = ReplacementRule(pattern2312, replacement2312)
    pattern2313 = Pattern(Integral(S(1)/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons196)
    def replacement2313(x, a, d, n):
        rubi.append(2313)
        return -Subst(Int(S(1)/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)
    rule2313 = ReplacementRule(pattern2313, replacement2313)
    pattern2314 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons4, cons5, cons1278)
    def replacement2314(a, n, p, x, d, c):
        rubi.append(2314)
        return Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)
    rule2314 = ReplacementRule(pattern2314, replacement2314)
    pattern2315 = Pattern(Integral(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons803, cons1279)
    def replacement2315(a, n, p, x, d):
        rubi.append(2315)
        return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)
    rule2315 = ReplacementRule(pattern2315, replacement2315)
    pattern2316 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons803, cons1280, cons1281)
    def replacement2316(a, n, p, x, d, c):
        rubi.append(2316)
        return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(c*n, S(2)))*Rt(Pi*c*n, S(2))/(d*n), x)
    rule2316 = ReplacementRule(pattern2316, replacement2316)
    pattern2317 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons803, cons1280, cons1282)
    def replacement2317(a, n, p, x, d, c):
        rubi.append(2317)
        return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(-c*n, S(2)))*Rt(-Pi*c*n, S(2))/(d*n), x)
    rule2317 = ReplacementRule(pattern2317, replacement2317)
    pattern2318 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons338, cons88, cons1283)
    def replacement2318(a, n, p, x, d, c):
        rubi.append(2318)
        return -Dist(c*(n*(p + S(-1)) + S(1)), Int((c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)
    rule2318 = ReplacementRule(pattern2318, replacement2318)
    pattern2319 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons338, cons88, cons1284)
    def replacement2319(a, n, p, x, d, c):
        rubi.append(2319)
        return -Dist(S(1)/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(d*(n*p + S(1))), x)
    rule2319 = ReplacementRule(pattern2319, replacement2319)
    pattern2320 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons5, cons196)
    def replacement2320(a, n, p, x, d, c):
        rubi.append(2320)
        return -Subst(Int((c*ProductLog(a*x**(-n)))**p/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)
    rule2320 = ReplacementRule(pattern2320, replacement2320)
    pattern2321 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons31, cons168)
    def replacement2321(x, a, d, m):
        rubi.append(2321)
        return -Dist(m/(m + S(1)), Int(x**m/((d*ProductLog(a*x) + d)*ProductLog(a*x)), x), x) + Simp(x**(m + S(1))/(d*(m + S(1))*ProductLog(a*x)), x)
    rule2321 = ReplacementRule(pattern2321, replacement2321)
    pattern2322 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons27, cons1285)
    def replacement2322(x, a, d):
        rubi.append(2322)
        return Simp(log(ProductLog(a*x))/d, x)
    rule2322 = ReplacementRule(pattern2322, replacement2322)
    pattern2323 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons31, cons94)
    def replacement2323(x, a, d, m):
        rubi.append(2323)
        return -Int(x**m*ProductLog(a*x)/(d*ProductLog(a*x) + d), x) + Simp(x**(m + S(1))/(d*(m + S(1))), x)
    rule2323 = ReplacementRule(pattern2323, replacement2323)
    pattern2324 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons21, cons18)
    def replacement2324(x, a, d, m):
        rubi.append(2324)
        return Simp(x**m*(-(m + S(1))*ProductLog(a*x))**(-m)*Gamma(m + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)
    rule2324 = ReplacementRule(pattern2324, replacement2324)
    pattern2325 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons27, cons4, cons1286)
    def replacement2325(x, a, d, n):
        rubi.append(2325)
        return Simp(log(ProductLog(a*x**n))/(d*n), x)
    rule2325 = ReplacementRule(pattern2325, replacement2325)
    pattern2326 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons150, cons463, cons66)
    def replacement2326(a, n, x, d, m):
        rubi.append(2326)
        return -Subst(Int(x**(-m + S(-2))/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)
    rule2326 = ReplacementRule(pattern2326, replacement2326)
    pattern2327 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons7, cons27, cons4, cons5, cons1287)
    def replacement2327(a, n, p, x, d, c):
        rubi.append(2327)
        return Simp((c*ProductLog(a*x**n))**p/(d*n*p), x)
    rule2327 = ReplacementRule(pattern2327, replacement2327)
    pattern2328 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1288)
    def replacement2328(a, n, p, x, d, m, c):
        rubi.append(2328)
        return Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)
    rule2328 = ReplacementRule(pattern2328, replacement2328)
    pattern2329 = Pattern(Integral(x_**WC('m', S(1))*ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons27, cons21, cons4, cons38, cons1289)
    def replacement2329(a, n, p, x, d, m):
        rubi.append(2329)
        return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)
    rule2329 = ReplacementRule(pattern2329, replacement2329)
    pattern2330 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons66, cons347, cons1290, cons1291)
    def replacement2330(a, n, p, x, d, m, c):
        rubi.append(2330)
        return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(c/(p + S(-1)/2), S(2)))*Rt(Pi*c/(p + S(-1)/2), S(2))/(d*n), x)
    rule2330 = ReplacementRule(pattern2330, replacement2330)
    pattern2331 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons66, cons347, cons1290, cons1292)
    def replacement2331(a, n, p, x, d, m, c):
        rubi.append(2331)
        return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(-c/(p + S(-1)/2), S(2)))*Rt(-Pi*c/(p + S(-1)/2), S(2))/(d*n), x)
    rule2331 = ReplacementRule(pattern2331, replacement2331)
    pattern2332 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1293, cons1294)
    def replacement2332(a, n, p, x, d, m, c):
        rubi.append(2332)
        return -Dist(c*(m + n*(p + S(-1)) + S(1))/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)
    rule2332 = ReplacementRule(pattern2332, replacement2332)
    pattern2333 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons4, cons5, cons66, cons1293, cons1295)
    def replacement2333(a, n, p, x, d, m, c):
        rubi.append(2333)
        return -Dist((m + S(1))/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(d*(m + n*p + S(1))), x)
    rule2333 = ReplacementRule(pattern2333, replacement2333)
    pattern2334 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons21, cons5, cons66)
    def replacement2334(a, p, x, d, m, c):
        rubi.append(2334)
        return Simp(x**m*(c*ProductLog(a*x))**p*(-(m + S(1))*ProductLog(a*x))**(-m - p)*Gamma(m + p + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)
    rule2334 = ReplacementRule(pattern2334, replacement2334)
    pattern2335 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons7, cons27, cons5, cons66, cons150, cons463)
    def replacement2335(a, n, p, x, d, m, c):
        rubi.append(2335)
        return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)
    rule2335 = ReplacementRule(pattern2335, replacement2335)
    pattern2336 = Pattern(Integral(u_, x_), cons1296)
    def replacement2336(x, u):
        rubi.append(2336)
        return Subst(Int(SimplifyIntegrand((x + S(1))*SubstFor(ProductLog(x), u, x)*exp(x), x), x), x, ProductLog(x))
    rule2336 = ReplacementRule(pattern2336, replacement2336)
    return [rule2145, rule2146, rule2147, rule2148, rule2149, rule2150, rule2151, rule2152, rule2153, rule2154, rule2155, rule2156, rule2157, rule2158, rule2159, rule2160, rule2161, rule2162, rule2163, rule2164, rule2165, rule2166, rule2167, rule2168, rule2169, rule2170, rule2171, rule2172, rule2173, rule2174, rule2175, rule2176, rule2177, rule2178, rule2179, rule2180, rule2181, rule2182, rule2183, rule2184, rule2185, rule2186, rule2187, rule2188, rule2189, rule2190, rule2191, rule2192, rule2193, rule2194, rule2195, rule2196, rule2197, rule2198, rule2199, rule2200, rule2201, rule2202, rule2203, rule2204, rule2205, rule2206, rule2207, rule2208, rule2209, rule2210, rule2211, rule2212, rule2213, rule2214, rule2215, rule2216, rule2217, rule2218, rule2219, rule2220, rule2221, rule2222, rule2223, rule2224, rule2225, rule2226, rule2227, rule2228, rule2229, rule2230, rule2231, rule2232, rule2233, rule2234, rule2235, rule2236, rule2237, rule2238, rule2239, rule2240, rule2241, rule2242, rule2243, rule2244, rule2245, rule2246, rule2247, rule2248, rule2249, rule2250, rule2251, rule2252, rule2253, rule2254, rule2255, rule2256, rule2257, rule2258, rule2259, rule2260, rule2261, rule2262, rule2263, rule2264, rule2265, rule2266, rule2267, rule2268, rule2269, rule2270, rule2271, rule2272, rule2273, rule2274, rule2275, rule2276, rule2277, rule2278, rule2279, rule2280, rule2281, rule2282, rule2283, rule2284, rule2285, rule2286, rule2287, rule2288, rule2289, rule2290, rule2291, rule2292, rule2293, rule2294, rule2295, rule2296, rule2297, rule2298, rule2299, rule2300, rule2301, rule2302, rule2303, rule2304, rule2305, rule2306, rule2307, rule2308, rule2309, rule2310, rule2311, rule2312, rule2313, rule2314, rule2315, rule2316, rule2317, rule2318, rule2319, rule2320, rule2321, rule2322, rule2323, rule2324, rule2325, rule2326, rule2327, rule2328, rule2329, rule2330, rule2331, rule2332, rule2333, rule2334, rule2335, rule2336, ]
