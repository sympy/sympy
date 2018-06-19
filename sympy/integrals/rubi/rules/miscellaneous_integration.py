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
        PolynomialRemainder, Factor, DerivativeDivides, Rule
    )
    from sympy import Integral, S, sqrt, And, Or, Integer, Float, Mod
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii , Pqq, Q, R, r = symbols('i ii Pqq Q R r')
    _UseGamma = False

def miscellaneous_integration(rubi):
    from sympy.integrals.rubi.constraints import cons147, cons1090, cons2, cons3, cons7, cons4, cons5, cons386, cons27, cons50, cons1091, cons1092, cons1093, cons1094, cons48, cons125, cons208, cons34, cons35, cons36, cons1095, cons1096, cons66, cons21, cons84, cons1037, cons1036, cons38, cons1097, cons10, cons1098, cons1099, cons1100, cons209, cons1101, cons1102, cons1103, cons46, cons1104, cons1105, cons1106, cons1107, cons52, cons1108, cons800, cons1109, cons17, cons1110, cons586, cons1111, cons1112, cons1113, cons1114, cons1115, cons1116, cons1117, cons1118, cons1119, cons667, cons196, cons1120, cons840, cons1121, cons18, cons1122, cons148, cons45, cons1123, cons1124, cons1125, cons261, cons1126, cons367, cons1127, cons67, cons1128, cons744, cons1129, cons165, cons1130, cons1131, cons1132, cons1133, cons1134, cons347

    pattern1882 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons147, cons1090)
    def replacement1882(u, n, p, b, x, c, a):
        rubi.append(1882)
        return Dist(c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p)), Int(u*(a + b*x)**(n*p), x), x)
    rule1882 = ReplacementRule(pattern1882, replacement1882)
    pattern1883 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons147, cons386)
    def replacement1883(u, p, b, q, x, c, d, a):
        rubi.append(1883)
        return Dist((c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q), Int(u*(a + b*x)**(p*q), x), x)
    rule1883 = ReplacementRule(pattern1883, replacement1883)
    pattern1884 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons147, cons386)
    def replacement1884(u, n, p, b, q, x, c, d, a):
        rubi.append(1884)
        return Dist((c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q), Int(u*(a + b*x)**(n*p*q), x), x)
    rule1884 = ReplacementRule(pattern1884, replacement1884)
    pattern1885 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1095, cons1091, cons1092, cons1093, cons1094)
    def replacement1885(n, F, B, b, f, C, A, x, c, d, e, a, g):
        rubi.append(1885)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule1885 = ReplacementRule(pattern1885, replacement1885)
    pattern1886 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1095, cons1091, cons1096)
    def replacement1886(n, F, b, C, A, x, c, e, a, g):
        rubi.append(1886)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule1886 = ReplacementRule(pattern1886, replacement1886)
    pattern1887 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1095, cons1091, cons1092, cons1093, cons1094)
    def replacement1887(n, F, B, b, f, C, A, x, d, c, e, a, g):
        rubi.append(1887)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule1887 = ReplacementRule(pattern1887, replacement1887)
    pattern1888 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1095, cons1091, cons1096)
    def replacement1888(n, F, b, C, A, x, e, c, a, g):
        rubi.append(1888)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule1888 = ReplacementRule(pattern1888, replacement1888)
    def With1889(x, y, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1889 = Pattern(Integral(u_/y_, x_), CustomConstraint(With1889))
    def replacement1889(x, y, u):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1889)
        return Simp(q*log(RemoveContent(y, x)), x)
    rule1889 = ReplacementRule(pattern1889, replacement1889)
    def With1890(x, y, u, w):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(w*y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1890 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(With1890))
    def replacement1890(x, y, u, w):
        
        q = DerivativeDivides(w*y, u, x)
        rubi.append(1890)
        return Simp(q*log(RemoveContent(w*y, x)), x)
    rule1890 = ReplacementRule(pattern1890, replacement1890)
    def With1891(x, y, m, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1891 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons21, cons66, CustomConstraint(With1891))
    def replacement1891(x, y, m, u):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1891)
        return Simp(q*y**(m + S(1))/(m + S(1)), x)
    rule1891 = ReplacementRule(pattern1891, replacement1891)
    def With1892(y, u, n, x, m, z):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1892 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons21, cons4, cons66, CustomConstraint(With1892))
    def replacement1892(y, u, n, x, m, z):
        
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        rubi.append(1892)
        return Simp(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), x)
    rule1892 = ReplacementRule(pattern1892, replacement1892)
    def With1893(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = SimplifyIntegrand(u, x)
        if SimplerIntegrandQ(v, u, x):
            return True
        return False
    pattern1893 = Pattern(Integral(u_, x_), CustomConstraint(With1893))
    def replacement1893(x, u):
        
        v = SimplifyIntegrand(u, x)
        rubi.append(1893)
        return Int(v, x)
    rule1893 = ReplacementRule(pattern1893, replacement1893)
    pattern1894 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1037)
    def replacement1894(u, n, b, x, m, e, c, f, d, a):
        rubi.append(1894)
        return Dist((a*e**S(2) - c*f**S(2))**m, Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule1894 = ReplacementRule(pattern1894, replacement1894)
    pattern1895 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1036)
    def replacement1895(u, n, b, x, m, e, c, f, d, a):
        rubi.append(1895)
        return Dist((b*e**S(2) - d*f**S(2))**m, Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule1895 = ReplacementRule(pattern1895, replacement1895)
    pattern1896 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), cons2, cons21, cons4, cons38, cons1097, cons10)
    def replacement1896(u, n, p, w, v, x, m, a):
        rubi.append(1896)
        return Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x)
    rule1896 = ReplacementRule(pattern1896, replacement1896)
    def With1897(y, u, n, b, v, x, m, c, d, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1897 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1098, CustomConstraint(With1897))
    def replacement1897(y, u, n, b, v, x, m, c, d, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1897)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), x)
    rule1897 = ReplacementRule(pattern1897, replacement1897)
    def With1898(y, u, n, p, a, b, w, v, x, m, e, d, f, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1898 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons1098, cons1099, CustomConstraint(With1898))
    def replacement1898(y, u, n, p, a, b, w, v, x, m, e, d, f, c):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1898)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), x)
    rule1898 = ReplacementRule(pattern1898, replacement1898)
    def With1899(y, u, n, p, a, h, b, z, w, q, v, x, m, e, d, f, c, g):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern1899 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons1098, cons1099, cons1100, CustomConstraint(With1899))
    def replacement1899(y, u, n, p, a, h, b, z, w, q, v, x, m, e, d, f, c, g):
        
        r = DerivativeDivides(y, u, x)
        rubi.append(1899)
        return Dist(r, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), x)
    rule1899 = ReplacementRule(pattern1899, replacement1899)
    def With1900(y, u, n, b, x, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1900 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1101, CustomConstraint(With1900))
    def replacement1900(y, u, n, b, x, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1900)
        return Dist(a, Int(u, x), x) + Dist(b*q, Subst(Int(x**n, x), x, y), x)
    rule1900 = ReplacementRule(pattern1900, replacement1900)
    def With1901(y, u, n, p, b, x, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1901 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1102, CustomConstraint(With1901))
    def replacement1901(y, u, n, p, b, x, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1901)
        return Dist(q, Subst(Int((a + b*x**n)**p, x), x, y), x)
    rule1901 = ReplacementRule(pattern1901, replacement1901)
    def With1902(y, u, p, n, b, v, x, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern1902 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons1103, CustomConstraint(With1902))
    def replacement1902(y, u, p, n, b, v, x, m, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(1902)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n)**p, x), x, y), x)
    rule1902 = ReplacementRule(pattern1902, replacement1902)
    def With1903(y, u, n, p, b, v, x, c, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1903 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons1098, CustomConstraint(With1903))
    def replacement1903(y, u, n, p, b, v, x, c, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1903)
        return Dist(q, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule1903 = ReplacementRule(pattern1903, replacement1903)
    def With1904(y, u, p, n, B, b, w, v, A, x, c, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1904 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons46, cons1098, cons1099, CustomConstraint(With1904))
    def replacement1904(y, u, p, n, B, b, w, v, A, x, c, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1904)
        return Dist(q, Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule1904 = ReplacementRule(pattern1904, replacement1904)
    def With1905(y, u, p, n, B, w, A, x, c, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1905 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons4, cons5, cons46, cons1099, CustomConstraint(With1905))
    def replacement1905(y, u, p, n, B, w, A, x, c, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1905)
        return Dist(q, Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule1905 = ReplacementRule(pattern1905, replacement1905)
    def With1906(y, u, p, n, b, w, v, x, m, c, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern1906 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons1099, CustomConstraint(With1906))
    def replacement1906(y, u, p, n, b, w, v, x, m, c, n2, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(1906)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule1906 = ReplacementRule(pattern1906, replacement1906)
    def With1907(y, u, p, n, B, b, w, v, A, x, m, c, z, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern1907 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1098, cons1099, CustomConstraint(With1907))
    def replacement1907(y, u, p, n, B, b, w, v, A, x, m, c, z, n2, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(1907)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule1907 = ReplacementRule(pattern1907, replacement1907)
    def With1908(y, u, p, n, B, w, A, x, m, c, z, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern1908 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1099, CustomConstraint(With1908))
    def replacement1908(y, u, p, n, B, w, A, x, m, c, z, n2, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(1908)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule1908 = ReplacementRule(pattern1908, replacement1908)
    def With1909(y, u, p, n, b, v, x, m, c, d, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1909 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons1098, CustomConstraint(With1909))
    def replacement1909(y, u, p, n, b, v, x, m, c, d, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(1909)
        return Dist(q, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), x)
    rule1909 = ReplacementRule(pattern1909, replacement1909)
    def With1910(y, u, p, n, a, b, w, q, v, x, m, e, d, f, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern1910 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons1098, cons1099, CustomConstraint(With1910))
    def replacement1910(y, u, p, n, a, b, w, q, v, x, m, e, d, f, c):
        
        r = DerivativeDivides(y, u, x)
        rubi.append(1910)
        return Dist(r, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), x)
    rule1910 = ReplacementRule(pattern1910, replacement1910)
    def With1911(v, x, u, F):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1911 = Pattern(Integral(F_**v_*u_, x_), cons1095, cons1095, CustomConstraint(With1911))
    def replacement1911(v, x, u, F):
        
        q = DerivativeDivides(v, u, x)
        rubi.append(1911)
        return Simp(F**v*q/log(F), x)
    rule1911 = ReplacementRule(pattern1911, replacement1911)
    def With1912(u, F, w, v, x, m):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern1912 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), cons1095, cons21, cons1104, CustomConstraint(With1912))
    def replacement1912(u, F, w, v, x, m):
        
        q = DerivativeDivides(v, u, x)
        rubi.append(1912)
        return Dist(q, Subst(Int(F**x*x**m, x), x, v), x)
    rule1912 = ReplacementRule(pattern1912, replacement1912)
    def With1913(u, p, b, w, v, x, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(v*D(w, x) + w*D(v, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1913 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons38, CustomConstraint(With1913))
    def replacement1913(u, p, b, w, v, x, m, a):
        
        c = u/(v*D(w, x) + w*D(v, x))
        rubi.append(1913)
        return Dist(c, Subst(Int((a + b*x**p)**m, x), x, v*w), x)
    rule1913 = ReplacementRule(pattern1913, replacement1913)
    def With1914(u, p, b, w, q, v, x, m, r, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1914 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1105, cons1106, cons1107, CustomConstraint(With1914))
    def replacement1914(u, p, b, w, q, v, x, m, r, a):
        
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(1914)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w), x)
    rule1914 = ReplacementRule(pattern1914, replacement1914)
    def With1915(u, p, b, w, q, v, s, x, r, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1915 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1108, cons1106, cons1107, CustomConstraint(With1915))
    def replacement1915(u, p, b, w, q, v, s, x, r, m, a):
        
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(1915)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1))), x)
    rule1915 = ReplacementRule(pattern1915, replacement1915)
    def With1916(u, p, b, w, q, v, x, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1916 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons1109, cons38, cons17, CustomConstraint(With1916))
    def replacement1916(u, p, b, w, q, v, x, m, a):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(1916)
        return Dist(c*p, Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))), x)
    rule1916 = ReplacementRule(pattern1916, replacement1916)
    def With1917(u, p, b, w, q, v, x, r, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1917 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1110, cons586, cons17, CustomConstraint(With1917))
    def replacement1917(u, p, b, w, q, v, x, r, m, a):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(1917)
        return -Dist(c*q, Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w), x)
    rule1917 = ReplacementRule(pattern1917, replacement1917)
    def With1918(u, p, b, w, q, v, s, x, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1918 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons800, cons1111, cons1112, cons1113, cons17, CustomConstraint(With1918))
    def replacement1918(u, p, b, w, q, v, s, x, m, a):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(1918)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1))), x)
    rule1918 = ReplacementRule(pattern1918, replacement1918)
    def With1919(u, p, b, w, q, v, s, x, r, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern1919 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1114, cons1112, cons1113, cons17, CustomConstraint(With1919))
    def replacement1919(u, p, b, w, q, v, s, x, r, m, a):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(1919)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1))), x)
    rule1919 = ReplacementRule(pattern1919, replacement1919)
    pattern1920 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons21, cons66, cons1115)
    def replacement1920(x, m, u):
        rubi.append(1920)
        return Dist(S(1)/(m + S(1)), Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1))), x)
    rule1920 = ReplacementRule(pattern1920, replacement1920)
    def With1921(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfLinear(u, x)
        try:
            res = And(Not(FalseQ(lst)), SubstForFractionalPowerQ(u, Part(lst, S(3)), x))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern1921 = Pattern(Integral(u_, x_), CustomConstraint(With1921))
    def replacement1921(x, u):
        
        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(1921)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule1921 = ReplacementRule(pattern1921, replacement1921)
    def With1922(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern1922 = Pattern(Integral(u_, x_), CustomConstraint(With1922))
    def replacement1922(x, u):
        
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        rubi.append(1922)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule1922 = ReplacementRule(pattern1922, replacement1922)
    pattern1923 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons50, cons147, cons10, cons1116, cons1117)
    def replacement1923(u, n, p, w, q, v, x, m, z, a):
        rubi.append(1923)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p), Int(u*v**(m*p)*w**(n*p)*z**(p*q), x), x)
    rule1923 = ReplacementRule(pattern1923, replacement1923)
    pattern1924 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons147, cons10, cons1116)
    def replacement1924(u, n, p, w, v, x, m, a):
        rubi.append(1924)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p), Int(u*v**(m*p)*w**(n*p), x), x)
    rule1924 = ReplacementRule(pattern1924, replacement1924)
    pattern1925 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons5, cons147, cons10, cons1118, cons1119)
    def replacement1925(u, p, v, x, m, a):
        rubi.append(1925)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p), Int(u*v**(m*p), x), x)
    rule1925 = ReplacementRule(pattern1925, replacement1925)
    pattern1926 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons667, cons196, cons1120)
    def replacement1926(u, n, p, b, x, a):
        rubi.append(1926)
        return Dist(FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b)), Int(u*x**(n*p)*(a*x**(-n) + b)**p, x), x)
    rule1926 = ReplacementRule(pattern1926, replacement1926)
    pattern1927 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons147, cons196, cons840, cons1121)
    def replacement1927(u, n, p, b, v, x, a):
        rubi.append(1927)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b)**p, x), x)
    rule1927 = ReplacementRule(pattern1927, replacement1927)
    pattern1928 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons21, cons5, cons147, cons196, cons840)
    def replacement1928(u, n, p, b, v, x, m, a):
        rubi.append(1928)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x), x)
    rule1928 = ReplacementRule(pattern1928, replacement1928)
    def With1929(u, b, x, s, r, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        if Not(EqQ(v, S(1))):
            return True
        return False
    pattern1929 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons21, cons52, cons800, cons18, cons1122, CustomConstraint(With1929))
    def replacement1929(u, b, x, s, r, m, a):
        
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        rubi.append(1929)
        return Dist(v, Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), x)
    rule1929 = ReplacementRule(pattern1929, replacement1929)
    def With1930(u, n, b, x, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        if SumQ(v):
            return True
        return False
    pattern1930 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons148, CustomConstraint(With1930))
    def replacement1930(u, n, b, x, a):
        
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        rubi.append(1930)
        return Int(v, x)
    rule1930 = ReplacementRule(pattern1930, replacement1930)
    pattern1931 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38, cons1123)
    def replacement1931(u, n, p, a, b, x, n2, c):
        rubi.append(1931)
        return Dist(S(4)**(-p)*c**(-p), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1931 = ReplacementRule(pattern1931, replacement1931)
    pattern1932 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons1123)
    def replacement1932(u, n, p, a, b, x, n2, c):
        rubi.append(1932)
        return Dist((b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p, Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1932 = ReplacementRule(pattern1932, replacement1932)
    def With1933(u, n, a, b, x, n2, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        if SumQ(v):
            return True
        return False
    pattern1933 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons46, cons148, CustomConstraint(With1933))
    def replacement1933(u, n, a, b, x, n2, c):
        
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        rubi.append(1933)
        return Int(v, x)
    rule1933 = ReplacementRule(pattern1933, replacement1933)
    pattern1934 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons1124)
    def replacement1934(u, n, b, x, m, c, a):
        rubi.append(1934)
        return Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x)
    rule1934 = ReplacementRule(pattern1934, replacement1934)
    def With1935(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfLinear(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern1935 = Pattern(Integral(u_, x_), CustomConstraint(With1935))
    def replacement1935(x, u):
        
        lst = FunctionOfLinear(u, x)
        rubi.append(1935)
        return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)
    rule1935 = ReplacementRule(pattern1935, replacement1935)
    def With1936(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = PowerVariableExpn(u, S(0), x)
        try:
            res = And(Not(FalseQ(lst)), NonzeroQ(Part(lst, S(2))))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern1936 = Pattern(Integral(u_/x_, x_), cons1125, cons1120, CustomConstraint(With1936))
    def replacement1936(x, u):
        
        lst = PowerVariableExpn(u, S(0), x)
        rubi.append(1936)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule1936 = ReplacementRule(pattern1936, replacement1936)
    def With1937(x, m, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = PowerVariableExpn(u, m + S(1), x)
        try:
            res = And(Not(FalseQ(lst)), NonzeroQ(-m + Part(lst, S(2)) + S(-1)))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern1937 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons17, cons261, cons1125, cons1126, CustomConstraint(With1937))
    def replacement1937(x, m, u):
        
        lst = PowerVariableExpn(u, m + S(1), x)
        rubi.append(1937)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule1937 = ReplacementRule(pattern1937, replacement1937)
    def With1938(x, m, u):
        k = Denominator(m)
        rubi.append(1938)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(S(1)/k)), x)
    pattern1938 = Pattern(Integral(u_*x_**m_, x_), cons367)
    rule1938 = ReplacementRule(pattern1938, With1938)
    def With1939(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern1939 = Pattern(Integral(u_, x_), cons1127, CustomConstraint(With1939))
    def replacement1939(x, u):
        
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        rubi.append(1939)
        return Dist(S(2), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), x)
    rule1939 = ReplacementRule(pattern1939, replacement1939)
    pattern1940 = Pattern(Integral(S(1)/(a_ + v_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement1940(v, b, a, x):
        rubi.append(1940)
        return Dist(S(1)/(S(2)*a), Int(Together(S(1)/(-v/Rt(-a/b, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*a), Int(Together(S(1)/(v/Rt(-a/b, S(2)) + S(1))), x), x)
    rule1940 = ReplacementRule(pattern1940, replacement1940)
    pattern1941 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1128, cons744)
    def replacement1941(n, b, v, x, a):
        rubi.append(1941)
        return Dist(S(2)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x)
    rule1941 = ReplacementRule(pattern1941, replacement1941)
    pattern1942 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1129, cons165)
    def replacement1942(n, b, v, x, a):
        rubi.append(1942)
        return Dist(S(1)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x)
    rule1942 = ReplacementRule(pattern1942, replacement1942)
    pattern1943 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons148, cons1130)
    def replacement1943(u, n, b, v, x, a):
        rubi.append(1943)
        return Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x)
    rule1943 = ReplacementRule(pattern1943, replacement1943)
    def With1944(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = NormalizeIntegrand(u, x)
        if UnsameQ(v, u):
            return True
        return False
    pattern1944 = Pattern(Integral(u_, x_), CustomConstraint(With1944))
    def replacement1944(x, u):
        
        v = NormalizeIntegrand(u, x)
        rubi.append(1944)
        return Int(v, x)
    rule1944 = ReplacementRule(pattern1944, replacement1944)
    def With1945(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = ExpandIntegrand(u, x)
        if SumQ(v):
            return True
        return False
    pattern1945 = Pattern(Integral(u_, x_), CustomConstraint(With1945))
    def replacement1945(x, u):
        
        v = ExpandIntegrand(u, x)
        rubi.append(1945)
        return Int(v, x)
    rule1945 = ReplacementRule(pattern1945, replacement1945)
    pattern1946 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons1131, cons1132, cons1133, cons1134)
    def replacement1946(u, p, n, b, q, x, m, c, d, a):
        rubi.append(1946)
        return Dist(x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q, Int(u*x**(m*p), x), x)
    rule1946 = ReplacementRule(pattern1946, replacement1946)
    pattern1947 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons347)
    def replacement1947(u, n, p, a, b, x, n2, c):
        rubi.append(1947)
        return Dist((S(4)*c)**(-p + S(1)/2)*sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule1947 = ReplacementRule(pattern1947, replacement1947)
    def With1948(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfLinear(u, x)
        try:
            res = Not(FalseQ(lst))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern1948 = Pattern(Integral(u_, x_), CustomConstraint(With1948))
    def replacement1948(x, u):
        
        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(1948)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule1948 = ReplacementRule(pattern1948, replacement1948)
    # pattern1949 = Pattern(Integral(u_, x_))
    # def replacement1949(x, u):
    #     rubi.append(1949)
    #     return Int(u, x)
    # rule1949 = ReplacementRule(pattern1949, replacement1949)
    return [rule1882, rule1883, rule1884, rule1885, rule1886, rule1887, rule1888, rule1889, rule1890, rule1891, rule1892, rule1893, rule1894, rule1895, rule1896, rule1897, rule1898, rule1899, rule1900, rule1901, rule1902, rule1903, rule1904, rule1905, rule1906, rule1907, rule1908, rule1909, rule1910, rule1911, rule1912, rule1913, rule1914, rule1915, rule1916, rule1917, rule1918, rule1919, rule1920, rule1921, rule1922, rule1923, rule1924, rule1925, rule1926, rule1927, rule1928, rule1929, rule1930, rule1931, rule1932, rule1933, rule1934, rule1935, rule1936, rule1937, rule1938, rule1939, rule1940, rule1941, rule1942, rule1943, rule1944, rule1945, rule1946, rule1947, rule1948, ]
