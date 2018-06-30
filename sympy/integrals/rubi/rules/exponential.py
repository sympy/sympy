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
        SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Erf, Gamma,
        FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ,
        _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify,
        _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum,
        _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux,
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist, Sum_doit, PolynomialQuotient, Floor,
        PolynomialRemainder, Factor, DerivativeDivides, Rule, ExpIntegralEi, IntHide, log, exp
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
    i, ii , Pqq, Q, R, r = symbols('i ii Pqq Q R r')
    _UseGamma = False

def exponential(rubi):
    from sympy.integrals.rubi.constraints import cons31, cons168, cons515, cons1090, cons1091, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons94, cons17, cons18, cons21, cons1092, cons128, cons2, cons244, cons137, cons552, cons1093, cons1094, cons5, cons380, cons54, cons1095, cons1096, cons1097, cons209, cons224, cons796, cons797, cons50, cons1098, cons804, cons1099, cons812, cons1100, cons1101, cons1102, cons1103, cons584, cons1104, cons1105, cons479, cons480, cons1106, cons196, cons23, cons1107, cons53, cons1108, cons1109, cons1110, cons1111, cons85, cons1112, cons356, cons531, cons1113, cons1114, cons535, cons93, cons1115, cons1116, cons176, cons367, cons166, cons744, cons68, cons840, cons1117, cons1118, cons1119, cons25, cons71, cons1120, cons1121, cons1122, cons818, cons1123, cons1124, cons1125, cons1126, cons819, cons1127, cons1128, cons1129, cons1130, cons148, cons810, cons811, cons1131, cons1132, cons52, cons800, cons1133, cons1134, cons1135, cons813, cons1136, cons226, cons62, cons1137, cons1138, cons1139, cons1140, cons1141, cons1142, cons1143, cons463, cons1144, cons43, cons448, cons1145, cons1146, cons1147, cons1017

    pattern1882 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons31, cons168, cons515, cons1090)
    def replacement1882(x, n, c, b, F, g, f, m, e, d):
        rubi.append(1882)
        return -Dist(d*m/(f*g*n*log(F)), Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(-1)), x), x) + Simp((F**(g*(e + f*x))*b)**n*(c + d*x)**m/(f*g*n*log(F)), x)
    rule1882 = ReplacementRule(pattern1882, replacement1882)
    pattern1883 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons1091, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons31, cons94, cons515, cons1090)
    def replacement1883(x, n, c, b, F, g, f, m, e, d):
        rubi.append(1883)
        return -Dist(f*g*n*log(F)/(d*(m + S(1))), Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1)), x), x) + Simp((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)
    rule1883 = ReplacementRule(pattern1883, replacement1883)
    pattern1884 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons1091, cons7, cons27, cons48, cons125, cons208, cons1090)
    def replacement1884(x, c, F, g, f, e, d):
        rubi.append(1884)
        return Simp(F**(g*(-c*f/d + e))*ExpIntegralEi(f*g*(c + d*x)*log(F)/d)/d, x)
    rule1884 = ReplacementRule(pattern1884, replacement1884)
    pattern1885 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons7, cons27, cons48, cons125, cons208, cons17)
    def replacement1885(x, c, F, g, f, m, e, d):
        rubi.append(1885)
        return Simp(F**(g*(-c*f/d + e))*f**(-m + S(-1))*g**(-m + S(-1))*(-d)**m*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)*log(F)**(-m + S(-1)), x)
    rule1885 = ReplacementRule(pattern1885, replacement1885)
    pattern1886 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons1091, cons7, cons27, cons48, cons125, cons208, cons1090)
    def replacement1886(x, c, F, g, f, e, d):
        rubi.append(1886)
        return Dist(S(2)/d, Subst(Int(F**(g*(-c*f/d + e) + f*g*x**S(2)/d), x), x, sqrt(c + d*x)), x)
    rule1886 = ReplacementRule(pattern1886, replacement1886)
    pattern1887 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons1091, cons7, cons27, cons48, cons125, cons208, cons21, cons18)
    def replacement1887(x, c, F, g, f, m, e, d):
        rubi.append(1887)
        return -Simp(F**(g*(-c*f/d + e))*(-f*g*log(F)/d)**(-IntPart(m) + S(-1))*(-f*g*(c + d*x)*log(F)/d)**(-FracPart(m))*(c + d*x)**FracPart(m)*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)/d, x)
    rule1887 = ReplacementRule(pattern1887, replacement1887)
    pattern1888 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons1092)
    def replacement1888(x, c, n, b, F, g, f, m, e, d):
        rubi.append(1888)
        return Dist(F**(-g*n*(e + f*x))*(F**(g*(e + f*x))*b)**n, Int(F**(g*n*(e + f*x))*(c + d*x)**m, x), x)
    rule1888 = ReplacementRule(pattern1888, replacement1888)
    pattern1889 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons128)
    def replacement1889(x, n, c, b, F, g, f, p, a, m, e, d):
        rubi.append(1889)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b*(F**(g*(e + f*x)))**n)**p, x), x)
    rule1889 = ReplacementRule(pattern1889, replacement1889)
    pattern1890 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons31, cons168)
    def replacement1890(x, n, c, b, F, g, f, a, m, e, d):
        rubi.append(1890)
        return Dist(d*m/(a*f*g*n*log(F)), Int((c + d*x)**(m + S(-1))*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1)), x), x) - Simp((c + d*x)**m*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1))/(a*f*g*n*log(F)), x)
    rule1890 = ReplacementRule(pattern1890, replacement1890)
    def With1891(x, n, c, b, F, g, f, p, a, m, e, d):
        u = IntHide((a + b*(F**(g*(e + f*x)))**n)**p, x)
        rubi.append(1891)
        return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)
    pattern1891 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**p_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons244, cons168, cons137)
    rule1891 = ReplacementRule(pattern1891, With1891)
    pattern1892 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1091, cons2, cons3, cons208, cons4, cons5, cons552, cons1093, cons1094, cons17)
    def replacement1892(u, x, n, b, g, p, a, m, v, F):
        rubi.append(1892)
        return Int((a + b*(F**(g*ExpandToSum(v, x)))**n)**p*NormalizePowerOfLinear(u, x)**m, x)
    rule1892 = ReplacementRule(pattern1892, replacement1892)
    def With1893(u, x, n, b, g, p, a, m, v, F):
        uu = NormalizePowerOfLinear(u, x)
        z = Symbol('z')
        z = If(And(PowerQ(uu), FreeQ(Part(uu, 2), x)), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
        return Simp(uu**m*Int(z*(a + b*(F**(g*ExpandToSum(v, x)))**n)**p, x)/z, x)
    pattern1893 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1091, cons2, cons3, cons208, cons21, cons4, cons5, cons552, cons1093, cons1094, cons18)
    rule1893 = ReplacementRule(pattern1893, With1893)
    pattern1894 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons380)
    def replacement1894(x, n, c, b, F, g, f, p, a, m, e, d):
        rubi.append(1894)
        return Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m, x)
    rule1894 = ReplacementRule(pattern1894, replacement1894)
    pattern1895 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons4, cons31, cons168)
    def replacement1895(x, n, c, b, F, g, f, a, m, e, d):
        rubi.append(1895)
        return -Dist(d*m/(b*f*g*n*log(F)), Int((c + d*x)**(m + S(-1))*log(S(1) + b*(F**(g*(e + f*x)))**n/a), x), x) + Simp((c + d*x)**m*log(S(1) + b*(F**(g*(e + f*x)))**n/a)/(b*f*g*n*log(F)), x)
    rule1895 = ReplacementRule(pattern1895, replacement1895)
    pattern1896 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons54)
    def replacement1896(x, n, c, b, F, g, f, p, a, m, e, d):
        rubi.append(1896)
        return -Dist(d*m/(b*f*g*n*(p + S(1))*log(F)), Int((a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**(m + S(-1)), x), x) + Simp((a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**m/(b*f*g*n*(p + S(1))*log(F)), x)
    rule1896 = ReplacementRule(pattern1896, replacement1896)
    pattern1897 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons21, cons4, cons5, cons1095)
    def replacement1897(x, n, c, b, F, g, f, p, a, m, e, d):
        rubi.append(1897)
        return Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x)
    rule1897 = ReplacementRule(pattern1897, replacement1897)
    pattern1898 = Pattern(Integral((G_**((x_*WC('i', S(1)) + WC('h', S(0)))*WC('j', S(1)))*WC('k', S(1)))**WC('q', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons796, cons797, cons21, cons4, cons5, cons50, cons1096, cons1097)
    def replacement1898(h, j, x, n, c, b, k, F, m, g, f, p, G, q, a, i, e, d):
        rubi.append(1898)
        return Dist((G**(j*(h + i*x))*k)**q*(F**(g*(e + f*x)))**(-n), Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x), x)
    rule1898 = ReplacementRule(pattern1898, replacement1898)
    pattern1899 = Pattern(Integral((F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons4, cons1098)
    def replacement1899(x, n, b, a, c, F):
        rubi.append(1899)
        return Simp((F**(c*(a + b*x)))**n/(b*c*n*log(F)), x)
    rule1899 = ReplacementRule(pattern1899, replacement1899)
    pattern1900 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), cons1091, cons7, cons804, cons552, cons1099)
    def replacement1900(u, x, c, v, F):
        rubi.append(1900)
        return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*u, x), x)
    rule1900 = ReplacementRule(pattern1900, replacement1900)
    pattern1901 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), cons1091, cons7, cons804, cons552, cons1090)
    def replacement1901(u, x, c, v, F):
        rubi.append(1901)
        return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), u, x), x)
    rule1901 = ReplacementRule(pattern1901, replacement1901)
    pattern1902 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1091, cons7, cons21, cons812, cons1100)
    def replacement1902(u, x, w, m, c, v, F):
        rubi.append(1902)
        return Simp(F**(c*v)*u**(m + S(1))*Coefficient(w, x, S(1))/(c*Coefficient(u, x, S(1))*Coefficient(v, x, S(1))*log(F)), x)
    rule1902 = ReplacementRule(pattern1902, replacement1902)
    pattern1903 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1091, cons7, cons1101, cons552, cons1093, cons17, cons1099)
    def replacement1903(u, x, w, m, c, v, F):
        rubi.append(1903)
        return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*NormalizePowerOfLinear(u, x)**m, x), x)
    rule1903 = ReplacementRule(pattern1903, replacement1903)
    pattern1904 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1091, cons7, cons1101, cons552, cons1093, cons17, cons1090)
    def replacement1904(u, x, w, m, c, v, F):
        rubi.append(1904)
        return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), w*NormalizePowerOfLinear(u, x)**m, x), x)
    rule1904 = ReplacementRule(pattern1904, replacement1904)
    def With1905(u, x, w, m, c, v, F):
        uu = NormalizePowerOfLinear(u, x)
        z = Symbol('z')
        z = If(And(PowerQ(uu), FreeQ(Part(uu, 2), x)), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
        return Simp(uu**m*Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*z, x), x)/z, x)
    pattern1905 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1091, cons7, cons21, cons1101, cons552, cons1093, cons18)
    rule1905 = ReplacementRule(pattern1905, With1905)
    pattern1906 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1102, cons1103, cons584)
    def replacement1906(h, x, n, b, F, g, f, a, c, d, e):
        rubi.append(1906)
        return Simp(F**(c*(a + b*x))*e*x*log(d*x)**(n + S(1))/(n + S(1)), x)
    rule1906 = ReplacementRule(pattern1906, replacement1906)
    pattern1907 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons1104, cons1103, cons584)
    def replacement1907(h, x, n, b, F, g, f, a, m, c, d, e):
        rubi.append(1907)
        return Simp(F**(c*(a + b*x))*e*x**(m + S(1))*log(d*x)**(n + S(1))/(n + S(1)), x)
    rule1907 = ReplacementRule(pattern1907, replacement1907)
    pattern1908 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons1105)
    def replacement1908(x, b, d, a, c, F):
        rubi.append(1908)
        return Simp(F**(a + b*(c + d*x))/(b*d*log(F)), x)
    rule1908 = ReplacementRule(pattern1908, replacement1908)
    pattern1909 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons479)
    def replacement1909(x, b, d, a, c, F):
        rubi.append(1909)
        return Simp(F**a*sqrt(Pi)*Erfi((c + d*x)*Rt(b*log(F), S(2)))/(S(2)*d*Rt(b*log(F), S(2))), x)
    rule1909 = ReplacementRule(pattern1909, replacement1909)
    pattern1910 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons480)
    def replacement1910(x, b, d, a, c, F):
        rubi.append(1910)
        return Simp(F**a*sqrt(Pi)*Erf((c + d*x)*Rt(-b*log(F), S(2)))/(S(2)*d*Rt(-b*log(F), S(2))), x)
    rule1910 = ReplacementRule(pattern1910, replacement1910)
    pattern1911 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons1106, cons196)
    def replacement1911(x, n, b, d, a, c, F):
        rubi.append(1911)
        return -Dist(b*n*log(F), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**n, x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)/d, x)
    rule1911 = ReplacementRule(pattern1911, replacement1911)
    def With1912(x, n, b, d, a, c, F):
        k = Denominator(n)
        rubi.append(1912)
        return Dist(k/d, Subst(Int(F**(a + b*x**(k*n))*x**(k + S(-1)), x), x, (c + d*x)**(S(1)/k)), x)
    pattern1912 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons1106, cons23)
    rule1912 = ReplacementRule(pattern1912, With1912)
    pattern1913 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons4, cons1107)
    def replacement1913(x, n, b, d, a, c, F):
        rubi.append(1913)
        return -Simp(F**a*(-b*(c + d*x)**n*log(F))**(-S(1)/n)*(c + d*x)*Gamma(S(1)/n, -b*(c + d*x)**n*log(F))/(d*n), x)
    rule1913 = ReplacementRule(pattern1913, replacement1913)
    pattern1914 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons53, cons1108)
    def replacement1914(x, n, b, F, f, a, m, c, d, e):
        rubi.append(1914)
        return Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(-n)*(e + f*x)**n/(b*f*n*log(F)), x)
    rule1914 = ReplacementRule(pattern1914, replacement1914)
    pattern1915 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons1108)
    def replacement1915(x, n, b, F, f, a, c, d, e):
        rubi.append(1915)
        return Simp(F**a*ExpIntegralEi(b*(c + d*x)**n*log(F))/(f*n), x)
    rule1915 = ReplacementRule(pattern1915, replacement1915)
    pattern1916 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons21, cons4, cons1109)
    def replacement1916(x, n, b, F, a, m, c, d):
        rubi.append(1916)
        return Dist(S(1)/(d*(m + S(1))), Subst(Int(F**(a + b*x**S(2)), x), x, (c + d*x)**(m + S(1))), x)
    rule1916 = ReplacementRule(pattern1916, replacement1916)
    pattern1917 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons31, cons1110, cons1111, cons85, cons1112)
    def replacement1917(x, n, b, F, a, m, c, d):
        rubi.append(1917)
        return -Dist((m - n + S(1))/(b*n*log(F)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)), x)
    rule1917 = ReplacementRule(pattern1917, replacement1917)
    pattern1918 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons21, cons4, cons1110, cons1111, cons356, cons531)
    def replacement1918(x, n, b, F, a, m, c, d):
        rubi.append(1918)
        return -Dist((m - n + S(1))/(b*n*log(F)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)), x)
    rule1918 = ReplacementRule(pattern1918, replacement1918)
    pattern1919 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons31, cons1110, cons1113, cons85, cons1114)
    def replacement1919(x, n, b, F, a, m, c, d):
        rubi.append(1919)
        return -Dist(b*n*log(F)/(m + S(1)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)
    rule1919 = ReplacementRule(pattern1919, replacement1919)
    pattern1920 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons21, cons4, cons1110, cons1113, cons356, cons535)
    def replacement1920(x, n, b, F, a, m, c, d):
        rubi.append(1920)
        return -Dist(b*n*log(F)/(m + S(1)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)
    rule1920 = ReplacementRule(pattern1920, replacement1920)
    def With1921(x, n, b, F, a, m, c, d):
        k = Denominator(n)
        rubi.append(1921)
        return Dist(k/d, Subst(Int(F**(a + b*x**(k*n))*x**(k*(m + S(1)) + S(-1)), x), x, (c + d*x)**(S(1)/k)), x)
    pattern1921 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons93, cons1110, cons1111, cons23)
    rule1921 = ReplacementRule(pattern1921, With1921)
    pattern1922 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons1108, cons1110, cons1115, cons18, cons1116)
    def replacement1922(x, n, b, F, f, a, m, c, d, e):
        rubi.append(1922)
        return Dist((c + d*x)**(-m)*(e + f*x)**m, Int(F**(a + b*(c + d*x)**n)*(c + d*x)**m, x), x)
    rule1922 = ReplacementRule(pattern1922, replacement1922)
    pattern1923 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons1108)
    def replacement1923(x, n, b, F, f, a, m, c, d, e):
        rubi.append(1923)
        return -Simp(F**a*(-b*(c + d*x)**n*log(F))**(-(m + S(1))/n)*(e + f*x)**(m + S(1))*Gamma((m + S(1))/n, -b*(c + d*x)**n*log(F))/(f*n), x)
    rule1923 = ReplacementRule(pattern1923, replacement1923)
    pattern1924 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons176, cons367, cons166)
    def replacement1924(x, b, F, f, a, m, c, d, e):
        rubi.append(1924)
        return Dist((-c*f + d*e)/d, Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-1)), x), x) - Dist(f**S(2)*(m + S(-1))/(S(2)*b*d**S(2)*log(F)), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-2)), x), x) + Simp(F**(a + b*(c + d*x)**S(2))*f*(e + f*x)**(m + S(-1))/(S(2)*b*d**S(2)*log(F)), x)
    rule1924 = ReplacementRule(pattern1924, replacement1924)
    pattern1925 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons176, cons31, cons94)
    def replacement1925(x, b, F, f, a, m, c, d, e):
        rubi.append(1925)
        return -Dist(S(2)*b*d**S(2)*log(F)/(f**S(2)*(m + S(1))), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(2)), x), x) + Dist(S(2)*b*d*(-c*f + d*e)*log(F)/(f**S(2)*(m + S(1))), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1)), x), x) + Simp(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)
    rule1925 = ReplacementRule(pattern1925, replacement1925)
    pattern1926 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons176, cons85, cons744, cons31, cons94)
    def replacement1926(x, n, b, F, f, a, m, c, d, e):
        rubi.append(1926)
        return -Dist(b*d*n*log(F)/(f*(m + S(1))), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(n + S(-1))*(e + f*x)**(m + S(1)), x), x) + Simp(F**(a + b*(c + d*x)**n)*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)
    rule1926 = ReplacementRule(pattern1926, replacement1926)
    pattern1927 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons176)
    def replacement1927(x, b, F, f, a, c, d, e):
        rubi.append(1927)
        return Dist(d/f, Int(F**(a + b/(c + d*x))/(c + d*x), x), x) - Dist((-c*f + d*e)/f, Int(F**(a + b/(c + d*x))/((c + d*x)*(e + f*x)), x), x)
    rule1927 = ReplacementRule(pattern1927, replacement1927)
    pattern1928 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons176, cons17, cons94)
    def replacement1928(x, b, F, f, a, m, c, d, e):
        rubi.append(1928)
        return Dist(b*d*log(F)/(f*(m + S(1))), Int(F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(c + d*x)**S(2), x), x) + Simp(F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)
    rule1928 = ReplacementRule(pattern1928, replacement1928)
    pattern1929 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons176)
    def replacement1929(x, n, b, F, f, a, c, d, e):
        rubi.append(1929)
        return Int(F**(a + b*(c + d*x)**n)/(e + f*x), x)
    rule1929 = ReplacementRule(pattern1929, replacement1929)
    pattern1930 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), cons1091, cons21, cons68, cons840, cons1117)
    def replacement1930(u, x, m, v, F):
        rubi.append(1930)
        return Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x)
    rule1930 = ReplacementRule(pattern1930, replacement1930)
    pattern1931 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*u_, x_), cons1091, cons2, cons3, cons7, cons27, cons4, cons804)
    def replacement1931(u, x, n, b, d, a, c, F):
        rubi.append(1931)
        return Int(ExpandLinearProduct(F**(a + b*(c + d*x)**n), u, c, d, x), x)
    rule1931 = ReplacementRule(pattern1931, replacement1931)
    pattern1932 = Pattern(Integral(F_**(v_*WC('b', S(1)) + WC('a', S(0)))*WC('u', S(1)), x_), cons1091, cons2, cons3, cons804, cons1118, cons1119)
    def replacement1932(u, x, b, a, v, F):
        rubi.append(1932)
        return Int(F**(a + b*NormalizePowerOfLinear(v, x))*u, x)
    rule1932 = ReplacementRule(pattern1932, replacement1932)
    pattern1933 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/((x_*WC('f', S(1)) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons1108)
    def replacement1933(h, x, b, F, g, f, a, c, d, e):
        rubi.append(1933)
        return -Dist(d/(f*(-c*h + d*g)), Subst(Int(F**(a + b*d*x/(-c*h + d*g) - b*h/(-c*h + d*g))/x, x), x, (g + h*x)/(c + d*x)), x)
    rule1933 = ReplacementRule(pattern1933, replacement1933)
    pattern1934 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons25)
    def replacement1934(h, x, c, b, F, g, f, a, m, e, d):
        rubi.append(1934)
        return Dist(F**(b*f/d + e), Int((g + h*x)**m, x), x)
    rule1934 = ReplacementRule(pattern1934, replacement1934)
    pattern1935 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons71, cons1120)
    def replacement1935(h, x, c, b, F, g, f, a, m, e, d):
        rubi.append(1935)
        return Int(F**(-f*(-a*d + b*c)/(d*(c + d*x)) + (b*f + d*e)/d)*(g + h*x)**m, x)
    rule1935 = ReplacementRule(pattern1935, replacement1935)
    pattern1936 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons71, cons1121)
    def replacement1936(h, x, c, b, F, g, f, a, e, d):
        rubi.append(1936)
        return Dist(d/h, Int(F**(e + f*(a + b*x)/(c + d*x))/(c + d*x), x), x) - Dist((-c*h + d*g)/h, Int(F**(e + f*(a + b*x)/(c + d*x))/((c + d*x)*(g + h*x)), x), x)
    rule1936 = ReplacementRule(pattern1936, replacement1936)
    pattern1937 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons71, cons1121, cons17, cons94)
    def replacement1937(h, x, c, b, F, g, f, a, m, e, d):
        rubi.append(1937)
        return -Dist(f*(-a*d + b*c)*log(F)/(h*(m + S(1))), Int(F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(c + d*x)**S(2), x), x) + Simp(F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)
    rule1937 = ReplacementRule(pattern1937, replacement1937)
    pattern1938 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons1120)
    def replacement1938(h, j, x, c, b, F, g, f, a, i, e, d):
        rubi.append(1938)
        return -Dist(d/(h*(-c*j + d*i)), Subst(Int(F**(e - f*x*(-a*d + b*c)/(-c*j + d*i) + f*(-a*j + b*i)/(-c*j + d*i))/x, x), x, (i + j*x)/(c + d*x)), x)
    rule1938 = ReplacementRule(pattern1938, replacement1938)
    pattern1939 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons1122)
    def replacement1939(x, b, a, c, F):
        rubi.append(1939)
        return Dist(F**(a - b**S(2)/(S(4)*c)), Int(F**((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x)
    rule1939 = ReplacementRule(pattern1939, replacement1939)
    pattern1940 = Pattern(Integral(F_**v_, x_), cons1091, cons818, cons1123)
    def replacement1940(x, v, F):
        rubi.append(1940)
        return Int(F**ExpandToSum(v, x), x)
    rule1940 = ReplacementRule(pattern1940, replacement1940)
    pattern1941 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1124)
    def replacement1941(x, b, F, a, c, d, e):
        rubi.append(1941)
        return Simp(F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)), x)
    rule1941 = ReplacementRule(pattern1941, replacement1941)
    pattern1942 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1124, cons31, cons166)
    def replacement1942(x, b, F, a, m, c, d, e):
        rubi.append(1942)
        return -Dist(e**S(2)*(m + S(-1))/(S(2)*c*log(F)), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)), x)
    rule1942 = ReplacementRule(pattern1942, replacement1942)
    pattern1943 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1124)
    def replacement1943(x, b, F, a, c, d, e):
        rubi.append(1943)
        return Simp(F**(a - b**S(2)/(S(4)*c))*ExpIntegralEi((b + S(2)*c*x)**S(2)*log(F)/(S(4)*c))/(S(2)*e), x)
    rule1943 = ReplacementRule(pattern1943, replacement1943)
    pattern1944 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1124, cons31, cons94)
    def replacement1944(x, b, F, a, m, c, d, e):
        rubi.append(1944)
        return -Dist(S(2)*c*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)
    rule1944 = ReplacementRule(pattern1944, replacement1944)
    pattern1945 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1125)
    def replacement1945(x, b, F, a, c, d, e):
        rubi.append(1945)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(F**(a + b*x + c*x**S(2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)), x)
    rule1945 = ReplacementRule(pattern1945, replacement1945)
    pattern1946 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1125, cons31, cons166)
    def replacement1946(x, b, F, a, m, c, d, e):
        rubi.append(1946)
        return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-1)), x), x) - Dist(e**S(2)*(m + S(-1))/(S(2)*c*log(F)), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)), x)
    rule1946 = ReplacementRule(pattern1946, replacement1946)
    pattern1947 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1125, cons31, cons94)
    def replacement1947(x, b, F, a, m, c, d, e):
        rubi.append(1947)
        return -Dist(S(2)*c*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x), x) - Dist((b*e - S(2)*c*d)*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1)), x), x) + Simp(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)
    rule1947 = ReplacementRule(pattern1947, replacement1947)
    pattern1948 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons21, cons1126)
    def replacement1948(x, b, F, a, m, c, d, e):
        rubi.append(1948)
        return Int(F**(a + b*x + c*x**S(2))*(d + e*x)**m, x)
    rule1948 = ReplacementRule(pattern1948, replacement1948)
    pattern1949 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), cons1091, cons21, cons68, cons818, cons819)
    def replacement1949(u, x, m, v, F):
        rubi.append(1949)
        return Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x)
    rule1949 = ReplacementRule(pattern1949, replacement1949)
    def With1950(x, c, n, b, F, a, m, e, v, d):
        u = IntHide(F**(e*(c + d*x))*(F**v*b + a)**n, x)
        rubi.append(1950)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)
    pattern1950 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**v_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons1127, cons31, cons168, cons196)
    rule1950 = ReplacementRule(pattern1950, With1950)
    def With1951(h, x, c, n, b, F, g, f, G, a, e, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        m = FullSimplify(g*h*log(G)/(d*e*log(F)))
        if And(RationalQ(m), GreaterEqual(Abs(m), S(1))):
            return True
        return False
    pattern1951 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1128, CustomConstraint(With1951))
    def replacement1951(h, x, c, n, b, F, g, f, G, a, e, d):
        
        m = FullSimplify(g*h*log(G)/(d*e*log(F)))
        rubi.append(1951)
        return Dist(G**(-c*g*h/d + f*h)*Denominator(m)/(d*e*log(F)), Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m))), x)
    rule1951 = ReplacementRule(pattern1951, replacement1951)
    def With1952(h, x, c, n, b, F, g, f, G, a, e, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        m = FullSimplify(d*e*log(F)/(g*h*log(G)))
        if And(RationalQ(m), Greater(Abs(m), S(1))):
            return True
        return False
    pattern1952 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1128, CustomConstraint(With1952))
    def replacement1952(h, x, c, n, b, F, g, f, G, a, e, d):
        
        m = FullSimplify(d*e*log(F)/(g*h*log(G)))
        rubi.append(1952)
        return Dist(Denominator(m)/(g*h*log(G)), Subst(Int(x**(Denominator(m) + S(-1))*(F**(c*e - d*e*f/g)*b*x**Numerator(m) + a)**n, x), x, G**(h*(f + g*x)/Denominator(m))), x)
    rule1952 = ReplacementRule(pattern1952, replacement1952)
    pattern1953 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons1130, cons148)
    def replacement1953(h, x, c, n, b, F, g, f, G, a, e, d):
        rubi.append(1953)
        return Int(G**(f*h)*G**(g*h*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x)
    rule1953 = ReplacementRule(pattern1953, replacement1953)
    pattern1954 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons1130, cons196)
    def replacement1954(h, x, c, n, b, F, g, f, G, a, e, d):
        rubi.append(1954)
        return Simp(G**(h*(f + g*x))*a**n*Hypergeometric2F1(-n, g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G)), x)
    rule1954 = ReplacementRule(pattern1954, replacement1954)
    pattern1955 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons1130, cons23)
    def replacement1955(h, x, c, n, b, F, g, f, G, a, e, d):
        rubi.append(1955)
        return Simp(G**(h*(f + g*x))*(F**(e*(c + d*x))*b + a)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1) + g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(a*g*h*log(G)), x)
    rule1955 = ReplacementRule(pattern1955, replacement1955)
    pattern1956 = Pattern(Integral(G_**(u_*WC('h', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons2, cons3, cons48, cons209, cons4, cons810, cons811)
    def replacement1956(u, h, x, n, b, G, a, e, v, F):
        rubi.append(1956)
        return Int(G**(h*ExpandToSum(u, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x)
    rule1956 = ReplacementRule(pattern1956, replacement1956)
    def With1957(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        m = FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))
        if RationalQ(m):
            return True
        return False
    pattern1957 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons1132, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons52, cons800, cons1133, cons4, cons1131, CustomConstraint(With1957))
    def replacement1957(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        
        m = FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))
        rubi.append(1957)
        return Dist(G**(-c*g*h/d + f*h)*H**(-c*s*t/d + r*t)*Denominator(m)/(d*e*log(F)), Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m))), x)
    rule1957 = ReplacementRule(pattern1957, replacement1957)
    pattern1958 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons1132, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons52, cons800, cons1133, cons1134, cons85)
    def replacement1958(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        rubi.append(1958)
        return Dist(G**(h*(-c*g/d + f)), Int(H**(t*(r + s*x))*(b + F**(-e*(c + d*x))*a)**n, x), x)
    rule1958 = ReplacementRule(pattern1958, replacement1958)
    pattern1959 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1091, cons1129, cons1132, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons52, cons800, cons1133, cons1135, cons148)
    def replacement1959(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        rubi.append(1959)
        return Int(G**(f*h)*G**(g*h*x)*H**(r*t)*H**(s*t*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x)
    rule1959 = ReplacementRule(pattern1959, replacement1959)
    pattern1960 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons1132, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons52, cons800, cons1133, cons1135, cons196)
    def replacement1960(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        rubi.append(1960)
        return Simp(G**(h*(f + g*x))*H**(t*(r + s*x))*a**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)), x)
    rule1960 = ReplacementRule(pattern1960, replacement1960)
    pattern1961 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons1132, cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons52, cons800, cons1133, cons4, cons1135, cons23)
    def replacement1961(h, t, x, c, n, b, r, F, g, f, s, G, a, e, H, d):
        rubi.append(1961)
        return Simp(G**(h*(f + g*x))*H**(t*(r + s*x))*((F**(e*(c + d*x))*b + a)/a)**(-n)*(F**(e*(c + d*x))*b + a)**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)), x)
    rule1961 = ReplacementRule(pattern1961, replacement1961)
    pattern1962 = Pattern(Integral(G_**(u_*WC('h', S(1)))*H_**(w_*WC('t', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1091, cons1129, cons1132, cons2, cons3, cons48, cons209, cons1133, cons4, cons812, cons813)
    def replacement1962(u, h, t, x, n, b, w, G, a, e, H, v, F):
        rubi.append(1962)
        return Int(G**(h*ExpandToSum(u, x))*H**(t*ExpandToSum(w, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x)
    rule1962 = ReplacementRule(pattern1962, replacement1962)
    pattern1963 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons5, cons54)
    def replacement1963(x, c, n, b, F, p, a, e, d):
        rubi.append(1963)
        return -Dist(a*n/(b*d*e*log(F)), Int(x**(n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x), x) + Simp((F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)), x)
    rule1963 = ReplacementRule(pattern1963, replacement1963)
    pattern1964 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons21, cons4, cons5, cons54)
    def replacement1964(x, c, n, b, F, p, a, m, e, d):
        rubi.append(1964)
        return -Dist(a*n/(b*d*e*log(F)), Int(x**(m + n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x), x) - Dist(m/(b*d*e*(p + S(1))*log(F)), Int(x**(m + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1)), x), x) + Simp(x**m*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)), x)
    rule1964 = ReplacementRule(pattern1964, replacement1964)
    def With1965(u, x, b, g, f, a, m, c, v, F):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1965)
        return Dist(S(2)*c/q, Int((f + g*x)**m/(S(2)*F**u*c + b - q), x), x) - Dist(S(2)*c/q, Int((f + g*x)**m/(S(2)*F**u*c + b + q), x), x)
    pattern1965 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons125, cons208, cons1136, cons68, cons226, cons62)
    rule1965 = ReplacementRule(pattern1965, With1965)
    def With1966(u, x, b, g, f, a, m, c, v, F):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1966)
        return Dist(S(2)*c/q, Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b - q), x), x) - Dist(S(2)*c/q, Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b + q), x), x)
    pattern1966 = Pattern(Integral(F_**u_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons125, cons208, cons1136, cons68, cons226, cons62)
    rule1966 = ReplacementRule(pattern1966, With1966)
    def With1967(u, h, x, b, m, g, f, a, i, c, v, F):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        rubi.append(1967)
        return -Dist(-i + (-b*i + S(2)*c*h)/q, Int((f + g*x)**m/(S(2)*F**u*c + b + q), x), x) + Dist(i + (-b*i + S(2)*c*h)/q, Int((f + g*x)**m/(S(2)*F**u*c + b - q), x), x)
    pattern1967 = Pattern(Integral((F_**u_*WC('i', S(1)) + h_)*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons125, cons208, cons209, cons224, cons1136, cons68, cons226, cons62)
    rule1967 = ReplacementRule(pattern1967, With1967)
    def With1968(x, b, F, a, m, c, v, d):
        u = IntHide(S(1)/(F**v*b + F**(c + d*x)*a), x)
        rubi.append(1968)
        return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Simp(u*x**m, x)
    pattern1968 = Pattern(Integral(x_**WC('m', S(1))/(F_**v_*WC('b', S(1)) + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons1137, cons31, cons168)
    rule1968 = ReplacementRule(pattern1968, With1968)
    pattern1969 = Pattern(Integral(u_/(F_**v_*WC('b', S(1)) + F_**w_*WC('c', S(1)) + a_), x_), cons1091, cons2, cons3, cons7, cons552, cons1138, cons1139, cons1140)
    def replacement1969(u, x, b, w, a, c, v, F):
        rubi.append(1969)
        return Int(F**v*u/(F**(S(2)*v)*b + F**v*a + c), x)
    rule1969 = ReplacementRule(pattern1969, replacement1969)
    pattern1970 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons208, cons4, cons1141)
    def replacement1970(x, n, c, b, F, g, a, e, d):
        rubi.append(1970)
        return Int(ExpandIntegrand(F**(g*(d + e*x)**n), S(1)/(a + b*x + c*x**S(2)), x), x)
    rule1970 = ReplacementRule(pattern1970, replacement1970)
    pattern1971 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons1091, cons2, cons7, cons27, cons48, cons208, cons4, cons1142)
    def replacement1971(x, n, F, g, a, c, d, e):
        rubi.append(1971)
        return Int(ExpandIntegrand(F**(g*(d + e*x)**n), S(1)/(a + c*x**S(2)), x), x)
    rule1971 = ReplacementRule(pattern1971, replacement1971)
    pattern1972 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(c_*x_**S(2) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons208, cons4, cons804, cons17)
    def replacement1972(u, x, n, c, b, F, g, a, m, e, d):
        rubi.append(1972)
        return Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + b*x + c*x**S(2)), x), x)
    rule1972 = ReplacementRule(pattern1972, replacement1972)
    pattern1973 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(a_ + c_*x_**S(2)), x_), cons1091, cons2, cons7, cons27, cons48, cons208, cons4, cons804, cons17)
    def replacement1973(u, x, n, c, F, g, a, m, e, d):
        rubi.append(1973)
        return Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + c*x**S(2)), x), x)
    rule1973 = ReplacementRule(pattern1973, replacement1973)
    pattern1974 = Pattern(Integral(F_**((x_**S(4)*WC('b', S(1)) + WC('a', S(0)))/x_**S(2)), x_), cons1091, cons2, cons3, cons1143)
    def replacement1974(a, b, x, F):
        rubi.append(1974)
        return -Simp(sqrt(Pi)*Erf((-x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*Exp(-S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))), x) + Simp(sqrt(Pi)*Erf((x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*Exp(S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))), x)
    rule1974 = ReplacementRule(pattern1974, replacement1974)
    pattern1975 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('m', S(1)) + exp(x_))**n_, x_), cons93, cons168, cons463, cons1144)
    def replacement1975(m, x, n):
        rubi.append(1975)
        return Dist(m, Int(x**(m + S(-1))*(x**m + exp(x))**n, x), x) + Int((x**m + exp(x))**(n + S(1)), x) - Simp((x**m + exp(x))**(n + S(1))/(n + S(1)), x)
    rule1975 = ReplacementRule(pattern1975, replacement1975)
    pattern1976 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons43)
    def replacement1976(x, c, n, b, F, a, e, d):
        rubi.append(1976)
        return Dist(S(1)/(d*e*n*log(F)), Subst(Int(log(a + b*x)/x, x), x, (F**(e*(c + d*x)))**n), x)
    rule1976 = ReplacementRule(pattern1976, replacement1976)
    pattern1977 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1091, cons2, cons3, cons7, cons27, cons48, cons4, cons448)
    def replacement1977(x, c, n, b, F, a, e, d):
        rubi.append(1977)
        return -Dist(b*d*e*n*log(F), Int(x*(F**(e*(c + d*x)))**n/(a + b*(F**(e*(c + d*x)))**n), x), x) + Simp(x*log(a + b*(F**(e*(c + d*x)))**n), x)
    rule1977 = ReplacementRule(pattern1977, replacement1977)
    pattern1978 = Pattern(Integral((F_**v_*WC('a', S(1)))**n_*WC('u', S(1)), x_), cons1091, cons2, cons4, cons23)
    def replacement1978(u, x, n, a, v, F):
        rubi.append(1978)
        return Dist(F**(-n*v)*(F**v*a)**n, Int(F**(n*v)*u, x), x)
    rule1978 = ReplacementRule(pattern1978, replacement1978)
    def With1979(u, x):
        v = FunctionOfExponential(u, x)
        rubi.append(1979)
        return Dist(v/D(v, x), Subst(Int(FunctionOfExponentialFunction(u, x)/x, x), x, v), x)
    pattern1979 = Pattern(Integral(u_, x_), cons1145)
    rule1979 = ReplacementRule(pattern1979, With1979)
    pattern1980 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1091, cons2, cons3, cons4, cons196, cons1146)
    def replacement1980(u, x, n, b, w, a, v, F):
        rubi.append(1980)
        return Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x)
    rule1980 = ReplacementRule(pattern1980, replacement1980)
    pattern1981 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1091, cons1129, cons2, cons3, cons4, cons196, cons1146)
    def replacement1981(u, x, n, b, w, G, a, v, F):
        rubi.append(1981)
        return Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x)
    rule1981 = ReplacementRule(pattern1981, replacement1981)
    pattern1982 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1091, cons2, cons3, cons4, cons23, cons1146)
    def replacement1982(u, x, n, b, w, a, v, F):
        rubi.append(1982)
        return Dist(F**(-n*v)*(F**v*a + F**w*b)**n*(F**ExpandToSum(-v + w, x)*b + a)**(-n), Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x), x)
    rule1982 = ReplacementRule(pattern1982, replacement1982)
    pattern1983 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1091, cons1129, cons2, cons3, cons4, cons23, cons1146)
    def replacement1983(u, x, n, b, w, G, a, v, F):
        rubi.append(1983)
        return Dist(F**(-n*v)*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**(-n)*(F**v*a + G**w*b)**n, Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x), x)
    rule1983 = ReplacementRule(pattern1983, replacement1983)
    pattern1984 = Pattern(Integral(F_**v_*G_**w_*WC('u', S(1)), x_), cons1091, cons1129, cons1147)
    def replacement1984(u, x, w, G, v, F):
        rubi.append(1984)
        return Int(u*NormalizeIntegrand(exp(v*log(F) + w*log(G)), x), x)
    rule1984 = ReplacementRule(pattern1984, replacement1984)
    def With1985(u, x, w, y, v, F):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = v*y/(D(u, x)*log(F))
        if ZeroQ(-w*y + D(z, x)):
            return True
        return False
    pattern1985 = Pattern(Integral(F_**u_*(v_ + w_)*WC('y', S(1)), x_), cons1091, cons1091, CustomConstraint(With1985))
    def replacement1985(u, x, w, y, v, F):
        
        z = v*y/(D(u, x)*log(F))
        rubi.append(1985)
        return Simp(F**u*z, x)
    rule1985 = ReplacementRule(pattern1985, replacement1985)
    def With1986(u, x, n, w, v, F):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        z = v*D(u, x)*log(F) + (n + S(1))*D(v, x)
        if And(Equal(Exponent(w, x), Exponent(z, x)), ZeroQ(w*Coefficient(z, x, Exponent(z, x)) - z*Coefficient(w, x, Exponent(w, x)))):
            return True
        return False
    pattern1986 = Pattern(Integral(F_**u_*v_**WC('n', S(1))*w_, x_), cons1091, cons4, cons804, cons1017, cons1101, CustomConstraint(With1986))
    def replacement1986(u, x, n, w, v, F):
        
        z = v*D(u, x)*log(F) + (n + S(1))*D(v, x)
        rubi.append(1986)
        return Simp(F**u*v**(n + S(1))*Coefficient(w, x, Exponent(w, x))/Coefficient(z, x, Exponent(z, x)), x)
    rule1986 = ReplacementRule(pattern1986, replacement1986)
    return [rule1882, rule1883, rule1884, rule1885, rule1886, rule1887, rule1888, rule1889, rule1890, rule1891, rule1892, rule1893, rule1894, rule1895, rule1896, rule1897, rule1898, rule1899, rule1900, rule1901, rule1902, rule1903, rule1904, rule1905, rule1906, rule1907, rule1908, rule1909, rule1910, rule1911, rule1912, rule1913, rule1914, rule1915, rule1916, rule1917, rule1918, rule1919, rule1920, rule1921, rule1922, rule1923, rule1924, rule1925, rule1926, rule1927, rule1928, rule1929, rule1930, rule1931, rule1932, rule1933, rule1934, rule1935, rule1936, rule1937, rule1938, rule1939, rule1940, rule1941, rule1942, rule1943, rule1944, rule1945, rule1946, rule1947, rule1948, rule1949, rule1950, rule1951, rule1952, rule1953, rule1954, rule1955, rule1956, rule1957, rule1958, rule1959, rule1960, rule1961, rule1962, rule1963, rule1964, rule1965, rule1966, rule1967, rule1968, rule1969, rule1970, rule1971, rule1972, rule1973, rule1974, rule1975, rule1976, rule1977, rule1978, rule1979, rule1980, rule1981, rule1982, rule1983, rule1984, rule1985, rule1986, ]
