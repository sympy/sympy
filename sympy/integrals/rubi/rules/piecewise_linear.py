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
    from sympy.core.add import Add
    from sympy.core.mod import Mod
    from sympy.core.mul import Mul
    from sympy.core import EulerGamma
    from sympy.core.numbers import (Float, I, Integer)
    from sympy.core.power import Pow
    from sympy.core.singleton import S
    from sympy.functions.elementary.complexes import (Abs, sign)
    from sympy.functions.elementary.miscellaneous import sqrt
    from sympy.integrals.integrals import Integral
    from sympy.logic.boolalg import (And, Or)
    from sympy.simplify.simplify import simplify
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec, atan2)
    from sympy.core.numbers import pi as Pi

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii, Pqq, Q, R, r, C, k, u = symbols('i ii Pqq Q R r C k u')
    _UseGamma = False
    ShowSteps = False
    StepCounter = None


def piecewise_linear():
    from sympy.integrals.rubi.constraints import cons1092, cons19, cons1093, cons89, cons90, cons1094, cons91, cons25, cons74, cons68, cons4, cons1095, cons216, cons685, cons102, cons103, cons1096, cons1097, cons33, cons96, cons358, cons1098, cons21, cons1099, cons2, cons3


    pattern1885 = Pattern(Integral(u_**WC('m', S(1)), x_), cons19, cons1092)
    rule1885 = ReplacementRule(pattern1885, With1885)

    pattern1886 = Pattern(Integral(v_/u_, x_), cons1093, CustomConstraint(With1886))
    rule1886 = ReplacementRule(pattern1886, replacement1886)

    pattern1887 = Pattern(Integral(v_**n_/u_, x_), cons1093, cons89, cons90, cons1094, CustomConstraint(With1887))
    rule1887 = ReplacementRule(pattern1887, replacement1887)

    pattern1888 = Pattern(Integral(S(1)/(u_*v_), x_), cons1093, CustomConstraint(With1888))
    rule1888 = ReplacementRule(pattern1888, replacement1888)

    pattern1889 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), cons1093, CustomConstraint(With1889))
    rule1889 = ReplacementRule(pattern1889, replacement1889)

    pattern1890 = Pattern(Integral(S(1)/(u_*sqrt(v_)), x_), cons1093, CustomConstraint(With1890))
    rule1890 = ReplacementRule(pattern1890, replacement1890)

    pattern1891 = Pattern(Integral(v_**n_/u_, x_), cons1093, cons89, cons91, CustomConstraint(With1891))
    rule1891 = ReplacementRule(pattern1891, replacement1891)

    pattern1892 = Pattern(Integral(v_**n_/u_, x_), cons1093, cons25, CustomConstraint(With1892))
    rule1892 = ReplacementRule(pattern1892, replacement1892)

    pattern1893 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), cons1093, CustomConstraint(With1893))
    rule1893 = ReplacementRule(pattern1893, replacement1893)

    pattern1894 = Pattern(Integral(S(1)/(sqrt(u_)*sqrt(v_)), x_), cons1093, CustomConstraint(With1894))
    rule1894 = ReplacementRule(pattern1894, replacement1894)

    pattern1895 = Pattern(Integral(u_**m_*v_**n_, x_), cons19, cons4, cons1093, cons74, cons68, CustomConstraint(With1895))
    rule1895 = ReplacementRule(pattern1895, replacement1895)

    pattern1896 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), cons19, cons4, cons1093, cons68, cons1095, CustomConstraint(With1896))
    rule1896 = ReplacementRule(pattern1896, replacement1896)

    pattern1897 = Pattern(Integral(u_**m_*v_**WC('n', S(1)), x_), cons1093, cons216, cons89, cons90, cons685, cons102, cons103, CustomConstraint(With1897))
    rule1897 = ReplacementRule(pattern1897, replacement1897)

    pattern1898 = Pattern(Integral(u_**m_*v_**n_, x_), cons1093, cons685, cons1096, cons1097, CustomConstraint(With1898))
    rule1898 = ReplacementRule(pattern1898, replacement1898)

    pattern1899 = Pattern(Integral(u_**m_*v_**n_, x_), cons1093, cons216, cons33, cons96, CustomConstraint(With1899))
    rule1899 = ReplacementRule(pattern1899, replacement1899)

    pattern1900 = Pattern(Integral(u_**m_*v_**n_, x_), cons1093, cons358, cons1098, CustomConstraint(With1900))
    rule1900 = ReplacementRule(pattern1900, replacement1900)

    pattern1901 = Pattern(Integral(u_**m_*v_**n_, x_), cons1093, cons21, cons25, CustomConstraint(With1901))
    rule1901 = ReplacementRule(pattern1901, replacement1901)

    pattern1902 = Pattern(Integral(u_**WC('n', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons1092, cons1099, cons89, cons90)
    rule1902 = ReplacementRule(pattern1902, With1902)

    pattern1903 = Pattern(Integral(u_**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*log(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons1092, cons1099, cons89, cons90, cons68)
    rule1903 = ReplacementRule(pattern1903, With1903)
    return [rule1885, rule1886, rule1887, rule1888, rule1889, rule1890, rule1891, rule1892, rule1893, rule1894, rule1895, rule1896, rule1897, rule1898, rule1899, rule1900, rule1901, rule1902, rule1903, ]





def With1885(m, u, x):
    c = D(u, x)
    return Dist(S(1)/c, Subst(Int(x**m, x), x, u), x)


def With1886(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1886(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist((-a*v + b*u)/a, Int(S(1)/u, x), x) + Simp(b*x/a, x)


def With1887(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1887(n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist((-a*v + b*u)/a, Int(v**(n + S(-1))/u, x), x) + Simp(v**n/(a*n), x)


def With1888(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1888(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist(a/(-a*v + b*u), Int(S(1)/u, x), x) + Dist(b/(-a*v + b*u), Int(S(1)/v, x), x)


def With1889(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if And(NonzeroQ(-a*v + b*u), PosQ((-a*v + b*u)/a)):
        return True
    return False


def replacement1889(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(S(2)*ArcTan(sqrt(v)/Rt((-a*v + b*u)/a, S(2)))/(a*Rt((-a*v + b*u)/a, S(2))), x)


def With1890(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if And(NonzeroQ(-a*v + b*u), NegQ((-a*v + b*u)/a)):
        return True
    return False


def replacement1890(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(-S(2)*atanh(sqrt(v)/Rt(-(-a*v + b*u)/a, S(2)))/(a*Rt(-(-a*v + b*u)/a, S(2))), x)


def With1891(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1891(n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist(a/(-a*v + b*u), Int(v**(n + S(1))/u, x), x) + Simp(v**(n + S(1))/((n + S(1))*(-a*v + b*u)), x)


def With1892(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1892(n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(v**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), -a*v/(-a*v + b*u))/((n + S(1))*(-a*v + b*u)), x)


def With1893(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if And(NonzeroQ(-a*v + b*u), PosQ(a*b)):
        return True
    return False


def replacement1893(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(S(2)*atanh(sqrt(u)*Rt(a*b, S(2))/(a*sqrt(v)))/Rt(a*b, S(2)), x)


def With1894(u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if And(NonzeroQ(-a*v + b*u), NegQ(a*b)):
        return True
    return False


def replacement1894(u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(S(2)*ArcTan(sqrt(u)*Rt(-a*b, S(2))/(a*sqrt(v)))/Rt(-a*b, S(2)), x)


def With1895(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1895(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)


def With1896(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1896(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist(b*n/(a*(m + S(1))), Int(u**(m + S(1))*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + S(1))), x)


def With1897(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1897(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist(n*(-a*v + b*u)/(a*(m + n + S(1))), Int(u**m*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + n + S(1))), x)


def With1898(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1898(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return -Dist(n*(-a*v + b*u)/(a*(m + n + S(1))), Int(u**m*v**(n + S(-1)), x), x) + Simp(u**(m + S(1))*v**n/(a*(m + n + S(1))), x)


def With1899(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1899(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Dist(b*(m + n + S(2))/((m + S(1))*(-a*v + b*u)), Int(u**(m + S(1))*v**n, x), x) - Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)


def With1900(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1900(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Dist(b*(m + n + S(2))/((m + S(1))*(-a*v + b*u)), Int(u**(m + S(1))*v**n, x), x) - Simp(u**(m + S(1))*v**(n + S(1))/((m + S(1))*(-a*v + b*u)), x)


def With1901(m, n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    a = D(u, x)
    b = D(v, x)
    if NonzeroQ(-a*v + b*u):
        return True
    return False


def replacement1901(m, n, u, v, x):

    a = D(u, x)
    b = D(v, x)
    return Simp(u**m*v**(n + S(1))*(b*u/(-a*v + b*u))**(-m)*Hypergeometric2F1(-m, n + S(1), n + S(2), -a*v/(-a*v + b*u))/(b*(n + S(1))), x)


def With1902(a, b, n, u, x):
    c = D(u, x)
    return -Dist(c*n/b, Int(u**(n + S(-1))*(a + b*x)*log(a + b*x), x), x) - Int(u**n, x) + Simp(u**n*(a + b*x)*log(a + b*x)/b, x)


def With1903(a, b, m, n, u, x):
    c = D(u, x)
    return -Dist(c*n/(b*(m + S(1))), Int(u**(n + S(-1))*(a + b*x)**(m + S(1))*log(a + b*x), x), x) - Dist(S(1)/(m + S(1)), Int(u**n*(a + b*x)**m, x), x) + Simp(u**n*(a + b*x)**(m + S(1))*log(a + b*x)/(b*(m + S(1))), x)
