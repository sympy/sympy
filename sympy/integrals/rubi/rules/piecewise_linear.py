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
        Zeta, ProductLog, DerivativeDivides, HypergeometricPFQ, IntHide, OneQ, Null, exp, log, Discriminant,
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

def piecewise_linear(rubi, matcher, load_rule):
    from sympy.integrals.rubi.constraints import cons1090, cons21, cons1091, cons_with_1883, cons87, cons88, cons1092, cons_with_1884, cons_with_1885, cons_with_1886, cons_with_1887, cons89, cons_with_1888, cons23, cons_with_1889, cons_with_1890, cons_with_1891, cons72, cons66, cons4, cons_with_1892, cons1093, cons_with_1893, cons214, cons683, cons100, cons101, cons_with_1894, cons1094, cons1095, cons_with_1895, cons31, cons94, cons_with_1896, cons356, cons1096, cons_with_1897, cons18, cons_with_1898, cons1097, cons2, cons3

    def With1882(m, u, x):
        c = D(u, x)
        rubi.append(1882)
        return Dist(S(1)/c, Subst(Int(x**m, x), x, u), x)
    pattern1882 = Pattern(Integral(u_**WC('m', S(1)), x_), cons21, cons1090)
    rule1882 = ReplacementRule(pattern1882, With1882)

    if load_rule:
        matcher.add(pattern1882, 1882)

    pattern1883 = Pattern(Integral(v_/u_, x_), cons1091, cons_with_1883)
    def replacement1883(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1883)
        return -Dist((-a*v + b*u)/a, Int(S(1)/u, x), x) + Simp(b*x/a, x)
    rule1883 = ReplacementRule(pattern1883, replacement1883)

    if load_rule:
        matcher.add(pattern1883, 1883)

    pattern1884 = Pattern(Integral(v_**n_/u_, x_), cons1091, cons87, cons88, cons1092, cons_with_1884)
    def replacement1884(x, n, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1884)
        return -Dist((-a*v + b*u)/a, Int(v**(n + S(-1))/u, x), x) + Simp(v**n/(a*n), x)
    rule1884 = ReplacementRule(pattern1884, replacement1884)

    if load_rule:
        matcher.add(pattern1884, 1884)

    pattern1885 = Pattern(Integral(S(1)/(u_*v_), x_), cons1091, cons_with_1885)
    def replacement1885(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1885)
        return -Dist(a/(-a*v + b*u), Int(S(1)/u, x), x) + Dist(b/(-a*v + b*u), Int(S(1)/v, x), x)
    rule1885 = ReplacementRule(pattern1885, replacement1885)

    if load_rule:
        matcher.add(pattern1885, 1885)

    pattern1886 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), cons1091, cons_with_1886)
    def replacement1886(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1886)
        return Simp(S(2)*ArcTan(sqrt(v)/Rt((-a*v + b*u)/a, S(2)))/(a*Rt((-a*v + b*u)/a, S(2))), x)
    rule1886 = ReplacementRule(pattern1886, replacement1886)

    if load_rule:
        matcher.add(pattern1886, 1886)

    pattern1887 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), cons1091, cons_with_1887)
    def replacement1887(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1887)
        return Simp(-S(2)*atanh(sqrt(v)/Rt(-(-a*v + b*u)/a, S(2)))/(a*Rt(-(-a*v + b*u)/a, S(2))), x)
    rule1887 = ReplacementRule(pattern1887, replacement1887)

    if load_rule:
        matcher.add(pattern1887, 1887)

    pattern1888 = Pattern(Integral(v_**n_/u_, x_), cons1091, cons87, cons89, cons_with_1888)
    def replacement1888(x, n, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1888)
        return -Dist(a/(-a*v + b*u), Int(v**(n + S(1))/u, x), x) + Simp(v**(n + S(1))/((n + S(1))*(-a*v + b*u)), x)
    rule1888 = ReplacementRule(pattern1888, replacement1888)

    if load_rule:
        matcher.add(pattern1888, 1888)

    pattern1889 = Pattern(Integral(v_**n_/u_, x_), cons1091, cons23, cons_with_1889)
    def replacement1889(x, n, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1889)
        return Simp(v**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), -a*v/(-a*v + b*u))/((n + S(1))*(-a*v + b*u)), x)
    rule1889 = ReplacementRule(pattern1889, replacement1889)

    if load_rule:
        matcher.add(pattern1889, 1889)

    pattern1890 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), cons1091, cons_with_1890)
    def replacement1890(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1890)
        return Simp(S(2)*atanh(sqrt(u)*Rt(a*b, S(2))/(a*sqrt(v)))/Rt(a*b, S(2)), x)
    rule1890 = ReplacementRule(pattern1890, replacement1890)

    if load_rule:
        matcher.add(pattern1890, 1890)

    pattern1891 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), cons1091, cons_with_1891)
    def replacement1891(x, v, u):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1891)
        return Simp(S(2)*ArcTan(sqrt(u)*Rt(-a*b, S(2))/(a*sqrt(v)))/Rt(-a*b, S(2)), x)
    rule1891 = ReplacementRule(pattern1891, replacement1891)

    if load_rule:
        matcher.add(pattern1891, 1891)

    pattern1892 = Pattern(Integral(u_**m_*v_**n_, x_), cons21, cons4, cons1091, cons72, cons66, cons_with_1892)
    def replacement1892(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1892)
        return -Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)
    rule1892 = ReplacementRule(pattern1892, replacement1892)

    if load_rule:
        matcher.add(pattern1892, 1892)

    pattern1893 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), cons21, cons4, cons1091, cons66, cons1093, cons_with_1893)
    def replacement1893(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1893)
        return -Dist(b*n/(a*(m + S(1))), Int(u**(m + S(1))*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + S(1))), x)
    rule1893 = ReplacementRule(pattern1893, replacement1893)

    if load_rule:
        matcher.add(pattern1893, 1893)

    pattern1894 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), cons1091, cons214, cons87, cons88, cons683, cons100, cons101, cons_with_1894)
    def replacement1894(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1894)
        return -Dist(n*(-a*v + b*u)/(a*(m + n + S(1))), Int(u**m*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + n + S(1))), x)
    rule1894 = ReplacementRule(pattern1894, replacement1894)

    if load_rule:
        matcher.add(pattern1894, 1894)

    pattern1895 = Pattern(Integral(u_**m_*v_**n_, x_), cons1091, cons683, cons1094, cons1095, cons_with_1895)
    def replacement1895(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1895)
        return -Dist(n*(-a*v + b*u)/(a*(m + n + S(1))), Int(u**m*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + n + S(1))), x)
    rule1895 = ReplacementRule(pattern1895, replacement1895)

    if load_rule:
        matcher.add(pattern1895, 1895)

    pattern1896 = Pattern(Integral(u_**m_*v_**n_, x_), cons1091, cons214, cons31, cons94, cons_with_1896)
    def replacement1896(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1896)
        return Dist(b*(m + n + S(2))/((m + S(1))*(-a*v + b*u)), Int(u**(m + S(1))*v**n, x), x) - Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)
    rule1896 = ReplacementRule(pattern1896, replacement1896)

    if load_rule:
        matcher.add(pattern1896, 1896)

    pattern1897 = Pattern(Integral(u_**m_*v_**n_, x_), cons1091, cons356, cons1096, cons_with_1897)
    def replacement1897(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1897)
        return Dist(b*(m + n + S(2))/((m + S(1))*(-a*v + b*u)), Int(u**(m + S(1))*v**n, x), x) - Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)
    rule1897 = ReplacementRule(pattern1897, replacement1897)

    if load_rule:
        matcher.add(pattern1897, 1897)

    pattern1898 = Pattern(Integral(u_**m_*v_**n_, x_), cons1091, cons18, cons23, cons_with_1898)
    def replacement1898(v, u, x, m, n):

        a = D(u, x)
        b = D(v, x)
        rubi.append(1898)
        return Simp(u**m*v**(n + S(1))*(b*u/(-a*v + b*u))**(-m)*Hypergeometric2F1(-m, n + S(1), n + S(2), -a*v/(-a*v + b*u))/(b*(n + S(1))), x)
    rule1898 = ReplacementRule(pattern1898, replacement1898)

    if load_rule:
        matcher.add(pattern1898, 1898)
    def With1899(u, x, b, a, n):
        c = D(u, x)
        rubi.append(1899)
        return -Dist(c*n/b, Int(u**(n + S(-1))*(a + b*x)*log(a + b*x), x), x) - Int(u**n, x) + Simp(u**n*(a + b*x)*log(a + b*x)/b, x)
    pattern1899 = Pattern(Integral(u_**WC('n', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons1090, cons1097, cons87, cons88)
    rule1899 = ReplacementRule(pattern1899, With1899)

    if load_rule:
        matcher.add(pattern1899, 1899)
    def With1900(u, m, x, b, a, n):
        c = D(u, x)
        rubi.append(1900)
        return -Dist(c*n/(b*(m + S(1))), Int(u**(n + S(-1))*(a + b*x)**(m + S(1))*log(a + b*x), x), x) - Dist(S(1)/(m + S(1)), Int(u**n*(a + b*x)**m, x), x) + Simp(u**n*(a + b*x)**(m + S(1))*log(a + b*x)/(b*(m + S(1))), x)
    pattern1900 = Pattern(Integral(u_**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons21, cons1090, cons1097, cons87, cons88, cons66)
    rule1900 = ReplacementRule(pattern1900, With1900)

    if load_rule:
        matcher.add(pattern1900, 1900)
    return matcher, [rule1882, rule1883, rule1884, rule1885, rule1886, rule1887, rule1888, rule1889, rule1890, rule1891, rule1892, rule1893, rule1894, rule1895, rule1896, rule1897, rule1898, rule1899, rule1900, ]
