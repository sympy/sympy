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


def miscellaneous_integration():
    from sympy.integrals.rubi.constraints import cons149, cons2004, cons2, cons3, cons8, cons4, cons5, cons388, cons29, cons52, cons2005, cons2006, cons2007, cons2008, cons50, cons127, cons210, cons36, cons37, cons38, cons1101, cons2009, cons68, cons19, cons86, cons1039, cons1038, cons40, cons2010, cons10, cons2011, cons2012, cons2013, cons211, cons1833, cons1246, cons2014, cons48, cons2015, cons2016, cons2017, cons2018, cons54, cons2019, cons802, cons2020, cons20, cons2021, cons588, cons2022, cons2023, cons2024, cons2025, cons2026, cons2027, cons2028, cons2029, cons2030, cons669, cons198, cons2031, cons842, cons2032, cons21, cons2033, cons150, cons47, cons2034, cons1856, cons1249, cons263, cons2035, cons369, cons2036, cons69, cons1481, cons746, cons1484, cons167, cons2037, cons2038, cons1678, cons1257, cons2039, cons349


    pattern6934 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons149, cons2004)
    rule6934 = ReplacementRule(pattern6934, replacement6934)

    pattern6935 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons52, cons149, cons388)
    rule6935 = ReplacementRule(pattern6935, replacement6935)

    pattern6936 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons52, cons149, cons388)
    rule6936 = ReplacementRule(pattern6936, replacement6936)

    pattern6937 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons36, cons37, cons38, cons1101, cons2005, cons2006, cons2007, cons2008)
    rule6937 = ReplacementRule(pattern6937, replacement6937)

    pattern6938 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons50, cons210, cons36, cons38, cons1101, cons2005, cons2009)
    rule6938 = ReplacementRule(pattern6938, replacement6938)

    pattern6939 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons36, cons37, cons38, cons1101, cons2005, cons2006, cons2007, cons2008)
    rule6939 = ReplacementRule(pattern6939, replacement6939)

    pattern6940 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons8, cons50, cons210, cons36, cons38, cons1101, cons2005, cons2009)
    rule6940 = ReplacementRule(pattern6940, replacement6940)

    pattern6941 = Pattern(Integral(u_/y_, x_), CustomConstraint(With6941))
    rule6941 = ReplacementRule(pattern6941, replacement6941)

    pattern6942 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(With6942))
    rule6942 = ReplacementRule(pattern6942, replacement6942)

    pattern6943 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons19, cons68, CustomConstraint(With6943))
    rule6943 = ReplacementRule(pattern6943, replacement6943)

    pattern6944 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons19, cons4, cons68, CustomConstraint(With6944))
    rule6944 = ReplacementRule(pattern6944, replacement6944)

    pattern6945 = Pattern(Integral(u_, x_), CustomConstraint(With6945))
    rule6945 = ReplacementRule(pattern6945, replacement6945)

    pattern6946 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons86, cons1039)
    rule6946 = ReplacementRule(pattern6946, replacement6946)

    pattern6947 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons86, cons1038)
    rule6947 = ReplacementRule(pattern6947, replacement6947)

    pattern6948 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), cons2, cons19, cons4, cons40, cons2010, cons10)
    rule6948 = ReplacementRule(pattern6948, replacement6948)

    pattern6949 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons2011, CustomConstraint(With6949))
    rule6949 = ReplacementRule(pattern6949, replacement6949)

    pattern6950 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons2011, cons2012, CustomConstraint(With6950))
    rule6950 = ReplacementRule(pattern6950, replacement6950)

    pattern6951 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons52, cons2011, cons2012, cons2013, CustomConstraint(With6951))
    rule6951 = ReplacementRule(pattern6951, replacement6951)

    pattern6952 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1833, CustomConstraint(With6952))
    rule6952 = ReplacementRule(pattern6952, replacement6952)

    pattern6953 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1246, CustomConstraint(With6953))
    rule6953 = ReplacementRule(pattern6953, replacement6953)

    pattern6954 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons2014, CustomConstraint(With6954))
    rule6954 = ReplacementRule(pattern6954, replacement6954)

    pattern6955 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons8, cons4, cons5, cons48, cons2011, CustomConstraint(With6955))
    rule6955 = ReplacementRule(pattern6955, replacement6955)

    pattern6956 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons4, cons5, cons48, cons2011, cons2012, CustomConstraint(With6956))
    rule6956 = ReplacementRule(pattern6956, replacement6956)

    pattern6957 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons8, cons36, cons37, cons4, cons5, cons48, cons2012, CustomConstraint(With6957))
    rule6957 = ReplacementRule(pattern6957, replacement6957)

    pattern6958 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons19, cons4, cons5, cons48, cons2012, CustomConstraint(With6958))
    rule6958 = ReplacementRule(pattern6958, replacement6958)

    pattern6959 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons36, cons37, cons19, cons4, cons5, cons48, cons2011, cons2012, CustomConstraint(With6959))
    rule6959 = ReplacementRule(pattern6959, replacement6959)

    pattern6960 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons8, cons36, cons37, cons19, cons4, cons5, cons48, cons2012, CustomConstraint(With6960))
    rule6960 = ReplacementRule(pattern6960, replacement6960)

    pattern6961 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons2011, CustomConstraint(With6961))
    rule6961 = ReplacementRule(pattern6961, replacement6961)

    pattern6962 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons52, cons2011, cons2012, CustomConstraint(With6962))
    rule6962 = ReplacementRule(pattern6962, replacement6962)

    pattern6963 = Pattern(Integral(F_**v_*u_, x_), cons1101, cons1101, CustomConstraint(With6963))
    rule6963 = ReplacementRule(pattern6963, replacement6963)

    pattern6964 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), cons1101, cons19, cons2015, CustomConstraint(With6964))
    rule6964 = ReplacementRule(pattern6964, replacement6964)

    pattern6965 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons40, CustomConstraint(With6965))
    rule6965 = ReplacementRule(pattern6965, replacement6965)

    pattern6966 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons54, cons2016, cons2017, cons2018, CustomConstraint(With6966))
    rule6966 = ReplacementRule(pattern6966, replacement6966)

    pattern6967 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons54, cons802, cons2019, cons2017, cons2018, CustomConstraint(With6967))
    rule6967 = ReplacementRule(pattern6967, replacement6967)

    pattern6968 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons2020, cons40, cons20, CustomConstraint(With6968))
    rule6968 = ReplacementRule(pattern6968, replacement6968)

    pattern6969 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons54, cons2021, cons588, cons20, CustomConstraint(With6969))
    rule6969 = ReplacementRule(pattern6969, replacement6969)

    pattern6970 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons802, cons2022, cons2023, cons2024, cons20, CustomConstraint(With6970))
    rule6970 = ReplacementRule(pattern6970, replacement6970)

    pattern6971 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons19, cons5, cons52, cons54, cons802, cons2025, cons2023, cons2024, cons20, CustomConstraint(With6971))
    rule6971 = ReplacementRule(pattern6971, replacement6971)

    pattern6972 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons19, cons68, cons2026)
    rule6972 = ReplacementRule(pattern6972, replacement6972)

    pattern6973 = Pattern(Integral(u_, x_), CustomConstraint(With6973))
    rule6973 = ReplacementRule(pattern6973, replacement6973)

    pattern6974 = Pattern(Integral(u_, x_), CustomConstraint(With6974))
    rule6974 = ReplacementRule(pattern6974, replacement6974)

    pattern6975 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons19, cons4, cons5, cons52, cons149, cons10, cons2027, cons2028)
    rule6975 = ReplacementRule(pattern6975, replacement6975)

    pattern6976 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons19, cons4, cons5, cons149, cons10, cons2027)
    rule6976 = ReplacementRule(pattern6976, replacement6976)

    pattern6977 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons19, cons5, cons149, cons10, cons2029, cons2030)
    rule6977 = ReplacementRule(pattern6977, replacement6977)

    pattern6978 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons669, cons198, cons2031)
    rule6978 = ReplacementRule(pattern6978, replacement6978)

    pattern6979 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons149, cons198, cons842, cons2032)
    rule6979 = ReplacementRule(pattern6979, replacement6979)

    pattern6980 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons19, cons5, cons149, cons198, cons842)
    rule6980 = ReplacementRule(pattern6980, replacement6980)

    pattern6981 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons19, cons54, cons802, cons21, cons2033, CustomConstraint(With6981))
    rule6981 = ReplacementRule(pattern6981, replacement6981)

    pattern6982 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons150, CustomConstraint(With6982))
    rule6982 = ReplacementRule(pattern6982, replacement6982)

    pattern6983 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons4, cons48, cons47, cons40, cons2034)
    rule6983 = ReplacementRule(pattern6983, replacement6983)

    pattern6984 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons47, cons149, cons2034)
    rule6984 = ReplacementRule(pattern6984, replacement6984)

    pattern6985 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons48, cons150, CustomConstraint(With6985))
    rule6985 = ReplacementRule(pattern6985, replacement6985)

    pattern6986 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons8, cons19, cons4, cons1856)
    rule6986 = ReplacementRule(pattern6986, replacement6986)

    pattern6987 = Pattern(Integral(u_, x_), CustomConstraint(With6987))
    rule6987 = ReplacementRule(pattern6987, replacement6987)

    pattern6988 = Pattern(Integral(u_/x_, x_), cons1249, cons2031, CustomConstraint(With6988))
    rule6988 = ReplacementRule(pattern6988, replacement6988)

    pattern6989 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons20, cons263, cons1249, cons2035, CustomConstraint(With6989))
    rule6989 = ReplacementRule(pattern6989, replacement6989)

    pattern6990 = Pattern(Integral(u_*x_**m_, x_), cons369)
    rule6990 = ReplacementRule(pattern6990, With6990)

    pattern6991 = Pattern(Integral(u_, x_), cons2036, CustomConstraint(With6991))
    rule6991 = ReplacementRule(pattern6991, replacement6991)

    pattern6992 = Pattern(Integral(S(1)/(a_ + v_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule6992 = ReplacementRule(pattern6992, replacement6992)

    pattern6993 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1481, cons746)
    rule6993 = ReplacementRule(pattern6993, replacement6993)

    pattern6994 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1484, cons167)
    rule6994 = ReplacementRule(pattern6994, replacement6994)

    pattern6995 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons150, cons2037)
    rule6995 = ReplacementRule(pattern6995, replacement6995)

    pattern6996 = Pattern(Integral(u_, x_), CustomConstraint(With6996))
    rule6996 = ReplacementRule(pattern6996, replacement6996)

    pattern6997 = Pattern(Integral(u_, x_), CustomConstraint(With6997))
    rule6997 = ReplacementRule(pattern6997, replacement6997)

    pattern6998 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons52, cons2038, cons1678, cons1257, cons2039)
    rule6998 = ReplacementRule(pattern6998, replacement6998)

    pattern6999 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons4, cons5, cons48, cons47, cons349)
    rule6999 = ReplacementRule(pattern6999, replacement6999)

    pattern7000 = Pattern(Integral(u_, x_), CustomConstraint(With7000))
    rule7000 = ReplacementRule(pattern7000, replacement7000)

    pattern7001 = Pattern(Integral(u_, x_))
    rule7001 = ReplacementRule(pattern7001, replacement7001)
    return [rule6934, rule6935, rule6936, rule6937, rule6938, rule6939, rule6940, rule6941, rule6942, rule6943, rule6944, rule6945, rule6946, rule6947, rule6948, rule6949, rule6950, rule6951, rule6952, rule6953, rule6954, rule6955, rule6956, rule6957, rule6958, rule6959, rule6960, rule6961, rule6962, rule6963, rule6964, rule6965, rule6966, rule6967, rule6968, rule6969, rule6970, rule6971, rule6972, rule6973, rule6974, rule6975, rule6976, rule6977, rule6978, rule6979, rule6980, rule6981, rule6982, rule6983, rule6984, rule6985, rule6986, rule6987, rule6988, rule6989, rule6990, rule6991, rule6992, rule6993, rule6994, rule6995, rule6996, rule6997, rule6998, rule6999, rule7000, rule7001, ]





def replacement6934(a, b, c, n, p, u, x):
    return Dist(c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p)), Int(u*(a + b*x)**(n*p), x), x)


def replacement6935(a, b, c, d, p, q, u, x):
    return Dist((c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q), Int(u*(a + b*x)**(p*q), x), x)


def replacement6936(a, b, c, d, n, p, q, u, x):
    return Dist((c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q), Int(u*(a + b*x)**(n*p*q), x), x)


def replacement6937(A, B, C, F, a, b, c, d, e, f, g, n, x):
    return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)


def replacement6938(A, C, F, a, b, c, e, g, n, x):
    return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)


def replacement6939(A, B, C, F, a, b, c, d, e, f, g, n, x):
    return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)


def replacement6940(A, C, F, a, b, c, e, g, n, x):
    return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)


def With6941(u, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6941(u, x, y):

    q = DerivativeDivides(y, u, x)
    return Simp(q*log(RemoveContent(y, x)), x)


def With6942(u, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(w*y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6942(u, w, x, y):

    q = DerivativeDivides(w*y, u, x)
    return Simp(q*log(RemoveContent(w*y, x)), x)


def With6943(m, u, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6943(m, u, x, y):

    q = DerivativeDivides(y, u, x)
    return Simp(q*y**(m + S(1))/(m + S(1)), x)


def With6944(m, n, u, x, y, z):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6944(m, n, u, x, y, z):

    q = DerivativeDivides(y*z, u*z**(-m + n), x)
    return Simp(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), x)


def With6945(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = SimplifyIntegrand(u, x)
    if SimplerIntegrandQ(v, u, x):
        return True
    return False


def replacement6945(u, x):

    v = SimplifyIntegrand(u, x)
    return Int(v, x)


def replacement6946(a, b, c, d, e, f, m, n, u, x):
    return Dist((a*e**S(2) - c*f**S(2))**m, Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)


def replacement6947(a, b, c, d, e, f, m, n, u, x):
    return Dist((b*e**S(2) - d*f**S(2))**m, Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)


def replacement6948(a, m, n, p, u, v, w, x):
    return Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x)


def With6949(a, b, c, d, m, n, u, v, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6949(a, b, c, d, m, n, u, v, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), x)


def With6950(a, b, c, d, e, f, m, n, p, u, v, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6950(a, b, c, d, e, f, m, n, p, u, v, w, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), x)


def With6951(a, b, c, d, e, f, g, h, m, n, p, q, u, v, w, x, y, z):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        r = DerivativeDivides(y, u, x)
        res = Not(FalseQ(r))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6951(a, b, c, d, e, f, g, h, m, n, p, q, u, v, w, x, y, z):

    r = DerivativeDivides(y, u, x)
    return Dist(r, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), x)


def With6952(a, b, n, u, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6952(a, b, n, u, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(a, Int(u, x), x) + Dist(b*q, Subst(Int(x**n, x), x, y), x)


def With6953(a, b, n, p, u, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6953(a, b, n, p, u, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((a + b*x**n)**p, x), x, y), x)


def With6954(a, b, m, n, p, u, v, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        res = And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6954(a, b, m, n, p, u, v, x, y):

    q = Symbol('q')
    r = Symbol('r')
    r = Divides(y**m, v**m, x)
    q = DerivativeDivides(y, u, x)
    return Dist(q*r, Subst(Int(x**m*(a + b*x**n)**p, x), x, y), x)


def With6955(a, b, c, n, n2, p, u, v, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6955(a, b, c, n, n2, p, u, v, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)


def With6956(A, B, a, b, c, n, n2, p, u, v, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6956(A, B, a, b, c, n, n2, p, u, v, w, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)


def With6957(A, B, a, c, n, n2, p, u, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6957(A, B, a, c, n, n2, p, u, w, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)


def With6958(a, b, c, m, n, n2, p, u, v, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        res = And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6958(a, b, c, m, n, n2, p, u, v, w, x, y):

    q = Symbol('q')
    r = Symbol('r')
    r = Divides(y**m, v**m, x)
    q = DerivativeDivides(y, u, x)
    return Dist(q*r, Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)


def With6959(A, B, a, b, c, m, n, n2, p, u, v, w, x, y, z):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        res = And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6959(A, B, a, b, c, m, n, n2, p, u, v, w, x, y, z):

    q = Symbol('q')
    r = Symbol('r')
    r = Divides(y**m, z**m, x)
    q = DerivativeDivides(y, u, x)
    return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)


def With6960(A, B, a, c, m, n, n2, p, u, w, x, y, z):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        res = And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x)))))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6960(A, B, a, c, m, n, n2, p, u, w, x, y, z):

    q = Symbol('q')
    r = Symbol('r')
    r = Divides(y**m, z**m, x)
    q = DerivativeDivides(y, u, x)
    return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)


def With6961(a, b, c, d, m, n, p, u, v, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(y, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6961(a, b, c, d, m, n, p, u, v, x, y):

    q = DerivativeDivides(y, u, x)
    return Dist(q, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), x)


def With6962(a, b, c, d, e, f, m, n, p, q, u, v, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        r = DerivativeDivides(y, u, x)
        res = Not(FalseQ(r))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6962(a, b, c, d, e, f, m, n, p, q, u, v, w, x, y):

    r = DerivativeDivides(y, u, x)
    return Dist(r, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), x)


def With6963(F, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(v, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6963(F, u, v, x):

    q = DerivativeDivides(v, u, x)
    return Simp(F**v*q/log(F), x)


def With6964(F, m, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        q = DerivativeDivides(v, u, x)
        res = Not(FalseQ(q))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6964(F, m, u, v, w, x):

    q = DerivativeDivides(v, u, x)
    return Dist(q, Subst(Int(F**x*x**m, x), x, v), x)


def With6965(a, b, m, p, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(v*D(w, x) + w*D(v, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6965(a, b, m, p, u, v, w, x):

    c = u/(v*D(w, x) + w*D(v, x))
    return Dist(c, Subst(Int((a + b*x**p)**m, x), x, v*w), x)


def With6966(a, b, m, p, q, r, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) + q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6966(a, b, m, p, q, r, u, v, w, x):

    c = u/(p*w*D(v, x) + q*v*D(w, x))
    return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w), x)


def With6967(a, b, m, p, q, r, s, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) + q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6967(a, b, m, p, q, r, s, u, v, w, x):

    c = u/(p*w*D(v, x) + q*v*D(w, x))
    return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1))), x)


def With6968(a, b, m, p, q, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) - q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6968(a, b, m, p, q, u, v, w, x):

    c = u/(p*w*D(v, x) - q*v*D(w, x))
    return Dist(c*p, Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))), x)


def With6969(a, b, m, p, q, r, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) - q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6969(a, b, m, p, q, r, u, v, w, x):

    c = u/(p*w*D(v, x) - q*v*D(w, x))
    return -Dist(c*q, Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w), x)


def With6970(a, b, m, p, q, s, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) - q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6970(a, b, m, p, q, s, u, v, w, x):

    c = u/(p*w*D(v, x) - q*v*D(w, x))
    return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1))), x)


def With6971(a, b, m, p, q, r, s, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    c = u/(p*w*D(v, x) - q*v*D(w, x))
    if FreeQ(c, x):
        return True
    return False


def replacement6971(a, b, m, p, q, r, s, u, v, w, x):

    c = u/(p*w*D(v, x) - q*v*D(w, x))
    return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1))), x)


def replacement6972(m, u, x):
    return Dist(S(1)/(m + S(1)), Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1))), x)


def With6973(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = SubstForFractionalPowerOfLinear(u, x)
        res = And(Not(FalseQ(lst)), SubstForFractionalPowerQ(u, Part(lst, S(3)), x))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6973(u, x):

    lst = SubstForFractionalPowerOfLinear(u, x)
    return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)


def With6974(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6974(u, x):

    lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
    return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)


def replacement6975(a, m, n, p, q, u, v, w, x, z):
    return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p), Int(u*v**(m*p)*w**(n*p)*z**(p*q), x), x)


def replacement6976(a, m, n, p, u, v, w, x):
    return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p), Int(u*v**(m*p)*w**(n*p), x), x)


def replacement6977(a, m, p, u, v, x):
    return Dist(a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p), Int(u*v**(m*p), x), x)


def replacement6978(a, b, n, p, u, x):
    return Dist(FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b)), Int(u*x**(n*p)*(a*x**(-n) + b)**p, x), x)


def replacement6979(a, b, n, p, u, v, x):
    return Dist(v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b)**p, x), x)


def replacement6980(a, b, m, n, p, u, v, x):
    return Dist(v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x), x)


def With6981(a, b, m, r, s, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
    if Not(EqQ(v, S(1))):
        return True
    return False


def replacement6981(a, b, m, r, s, u, x):

    v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
    return Dist(v, Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), x)


def With6982(a, b, n, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = RationalFunctionExpand(u/(a + b*x**n), x)
    if SumQ(v):
        return True
    return False


def replacement6982(a, b, n, u, x):

    v = RationalFunctionExpand(u/(a + b*x**n), x)
    return Int(v, x)


def replacement6983(a, b, c, n, n2, p, u, x):
    return Dist(S(4)**(-p)*c**(-p), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)


def replacement6984(a, b, c, n, n2, p, u, x):
    return Dist((b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p, Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)


def With6985(a, b, c, n, n2, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
    if SumQ(v):
        return True
    return False


def replacement6985(a, b, c, n, n2, u, x):

    v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
    return Int(v, x)


def replacement6986(a, b, c, m, n, u, x):
    return Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x)


def With6987(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = FunctionOfLinear(u, x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6987(u, x):

    lst = FunctionOfLinear(u, x)
    return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)


def With6988(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = PowerVariableExpn(u, S(0), x)
        res = And(Not(FalseQ(lst)), NonzeroQ(Part(lst, S(2))))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6988(u, x):

    lst = PowerVariableExpn(u, S(0), x)
    return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)


def With6989(m, u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = PowerVariableExpn(u, m + S(1), x)
        res = And(Not(FalseQ(lst)), NonzeroQ(-m + Part(lst, S(2)) + S(-1)))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6989(m, u, x):

    lst = PowerVariableExpn(u, m + S(1), x)
    return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)


def With6990(m, u, x):
    k = Denominator(m)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(S(1)/k)), x)


def With6991(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6991(u, x):

    lst = FunctionOfSquareRootOfQuadratic(u, x)
    return Dist(S(2), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), x)


def replacement6992(a, b, v, x):
    return Dist(S(1)/(S(2)*a), Int(Together(S(1)/(-v/Rt(-a/b, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*a), Int(Together(S(1)/(v/Rt(-a/b, S(2)) + S(1))), x), x)


def replacement6993(a, b, n, v, x):
    return Dist(S(2)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x)


def replacement6994(a, b, n, v, x):
    return Dist(S(1)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x)


def replacement6995(a, b, n, u, v, x):
    return Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x)


def With6996(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = NormalizeIntegrand(u, x)
    if UnsameQ(v, u):
        return True
    return False


def replacement6996(u, x):

    v = NormalizeIntegrand(u, x)
    return Int(v, x)


def With6997(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    v = ExpandIntegrand(u, x)
    if SumQ(v):
        return True
    return False


def replacement6997(u, x):

    v = ExpandIntegrand(u, x)
    return Int(v, x)


def replacement6998(a, b, c, d, m, n, p, q, u, x):
    return Dist(x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q, Int(u*x**(m*p), x), x)


def replacement6999(a, b, c, n, n2, p, u, x):
    return Dist((S(4)*c)**(S(1)/2 - p)*sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)


def With7000(u, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        lst = SubstForFractionalPowerOfLinear(u, x)
        res = Not(FalseQ(lst))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement7000(u, x):

    lst = SubstForFractionalPowerOfLinear(u, x)
    return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)


def replacement7001(u, x):
    return Int(u, x)
