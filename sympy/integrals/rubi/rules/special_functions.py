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


def special_functions():
    from sympy.integrals.rubi.constraints import cons69, cons2, cons3, cons68, cons19, cons1266, cons8, cons29, cons20, cons168, cons1959, cons1960, cons96, cons263, cons1961, cons1834, cons64, cons1962, cons1963, cons1964, cons249, cons1965, cons1966, cons1967, cons1833, cons4, cons1257, cons21, cons1361, cons1968, cons1969, cons170, cons1970, cons1971, cons33, cons1972, cons1973, cons1974, cons802, cons89, cons90, cons5, cons52, cons91, cons385, cons50, cons1975, cons1976, cons1977, cons54, cons1978, cons1101, cons127, cons1245, cons13, cons139, cons1381, cons1979, cons1980, cons198, cons1981, cons1982, cons1983, cons152, cons465, cons1767, cons165, cons950, cons951, cons1984, cons1985, cons805, cons1986, cons1987, cons1988, cons1989, cons340, cons1990, cons1991, cons1992, cons1993, cons1994, cons1995, cons40, cons1996, cons349, cons1997, cons1998, cons1999, cons2000, cons2001, cons2002, cons2003


    pattern6742 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6742 = ReplacementRule(pattern6742, replacement6742)

    pattern6743 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6743 = ReplacementRule(pattern6743, replacement6743)

    pattern6744 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6744 = ReplacementRule(pattern6744, replacement6744)

    pattern6745 = Pattern(Integral(Erf(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6745 = ReplacementRule(pattern6745, replacement6745)

    pattern6746 = Pattern(Integral(Erfc(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6746 = ReplacementRule(pattern6746, replacement6746)

    pattern6747 = Pattern(Integral(Erfi(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6747 = ReplacementRule(pattern6747, replacement6747)

    pattern6748 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6748 = ReplacementRule(pattern6748, replacement6748)

    pattern6749 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6749 = ReplacementRule(pattern6749, replacement6749)

    pattern6750 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6750 = ReplacementRule(pattern6750, replacement6750)

    pattern6751 = Pattern(Integral(x_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6751 = ReplacementRule(pattern6751, replacement6751)

    pattern6752 = Pattern(Integral(x_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6752 = ReplacementRule(pattern6752, replacement6752)

    pattern6753 = Pattern(Integral(x_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6753 = ReplacementRule(pattern6753, replacement6753)

    pattern6754 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons168)
    rule6754 = ReplacementRule(pattern6754, replacement6754)

    pattern6755 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons168)
    rule6755 = ReplacementRule(pattern6755, replacement6755)

    pattern6756 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons168)
    rule6756 = ReplacementRule(pattern6756, replacement6756)

    pattern6757 = Pattern(Integral(Erf(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1959)
    rule6757 = ReplacementRule(pattern6757, replacement6757)

    pattern6758 = Pattern(Integral(Erfc(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1959)
    rule6758 = ReplacementRule(pattern6758, replacement6758)

    pattern6759 = Pattern(Integral(Erfi(x_*WC('b', S(1)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0)))/x_, x_), cons3, cons1960)
    rule6759 = ReplacementRule(pattern6759, replacement6759)

    pattern6760 = Pattern(Integral(x_**m_*Erf(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6760 = ReplacementRule(pattern6760, replacement6760)

    pattern6761 = Pattern(Integral(x_**m_*Erfc(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6761 = ReplacementRule(pattern6761, replacement6761)

    pattern6762 = Pattern(Integral(x_**m_*Erfi(x_*WC('b', S(1)) + WC('a', S(0)))*exp(x_**S(2)*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6762 = ReplacementRule(pattern6762, replacement6762)

    pattern6763 = Pattern(Integral(Erf(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6763 = ReplacementRule(pattern6763, replacement6763)

    pattern6764 = Pattern(Integral(Erfc(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6764 = ReplacementRule(pattern6764, replacement6764)

    pattern6765 = Pattern(Integral(Erfi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6765 = ReplacementRule(pattern6765, replacement6765)

    pattern6766 = Pattern(Integral(x_**WC('m', S(1))*Erf(x_*WC('b', S(1)))**S(2), x_), cons3, cons20, cons263, cons1961)
    rule6766 = ReplacementRule(pattern6766, replacement6766)

    pattern6767 = Pattern(Integral(x_**WC('m', S(1))*Erfc(x_*WC('b', S(1)))**S(2), x_), cons3, cons20, cons1834, cons1961)
    rule6767 = ReplacementRule(pattern6767, replacement6767)

    pattern6768 = Pattern(Integral(x_**WC('m', S(1))*Erfi(x_*WC('b', S(1)))**S(2), x_), cons3, cons20, cons1834, cons1961)
    rule6768 = ReplacementRule(pattern6768, replacement6768)

    pattern6769 = Pattern(Integral(x_**WC('m', S(1))*Erf(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6769 = ReplacementRule(pattern6769, replacement6769)

    pattern6770 = Pattern(Integral(x_**WC('m', S(1))*Erfc(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6770 = ReplacementRule(pattern6770, replacement6770)

    pattern6771 = Pattern(Integral(x_**WC('m', S(1))*Erfi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6771 = ReplacementRule(pattern6771, replacement6771)

    pattern6772 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6772 = ReplacementRule(pattern6772, replacement6772)

    pattern6773 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6773 = ReplacementRule(pattern6773, replacement6773)

    pattern6774 = Pattern(Integral(FresnelS(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6774 = ReplacementRule(pattern6774, replacement6774)

    pattern6775 = Pattern(Integral(FresnelC(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6775 = ReplacementRule(pattern6775, replacement6775)

    pattern6776 = Pattern(Integral(x_**WC('m', S(1))*FresnelS(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6776 = ReplacementRule(pattern6776, replacement6776)

    pattern6777 = Pattern(Integral(x_**WC('m', S(1))*FresnelC(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6777 = ReplacementRule(pattern6777, replacement6777)

    pattern6778 = Pattern(Integral(FresnelS(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6778 = ReplacementRule(pattern6778, replacement6778)

    pattern6779 = Pattern(Integral(FresnelC(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6779 = ReplacementRule(pattern6779, replacement6779)

    pattern6780 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))**S(2), x_), cons3, cons20, cons1834, cons1962)
    rule6780 = ReplacementRule(pattern6780, replacement6780)

    pattern6781 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))**S(2), x_), cons3, cons20, cons1834, cons1962)
    rule6781 = ReplacementRule(pattern6781, replacement6781)

    pattern6782 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963)
    rule6782 = ReplacementRule(pattern6782, replacement6782)

    pattern6783 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963)
    rule6783 = ReplacementRule(pattern6783, replacement6783)

    pattern6784 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons168, cons1964)
    rule6784 = ReplacementRule(pattern6784, replacement6784)

    pattern6785 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons168, cons1964)
    rule6785 = ReplacementRule(pattern6785, replacement6785)

    pattern6786 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons249, cons1965)
    rule6786 = ReplacementRule(pattern6786, replacement6786)

    pattern6787 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons249, cons1965)
    rule6787 = ReplacementRule(pattern6787, replacement6787)

    pattern6788 = Pattern(Integral(x_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963)
    rule6788 = ReplacementRule(pattern6788, replacement6788)

    pattern6789 = Pattern(Integral(x_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963)
    rule6789 = ReplacementRule(pattern6789, replacement6789)

    pattern6790 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons168, cons1966)
    rule6790 = ReplacementRule(pattern6790, replacement6790)

    pattern6791 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons168, cons1966)
    rule6791 = ReplacementRule(pattern6791, replacement6791)

    pattern6792 = Pattern(Integral(x_**m_*FresnelS(x_*WC('b', S(1)))*cos(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons96, cons1967)
    rule6792 = ReplacementRule(pattern6792, replacement6792)

    pattern6793 = Pattern(Integral(x_**m_*FresnelC(x_*WC('b', S(1)))*sin(x_**S(2)*WC('c', S(1))), x_), cons3, cons8, cons1963, cons20, cons96, cons1967)
    rule6793 = ReplacementRule(pattern6793, replacement6793)

    pattern6794 = Pattern(Integral(ExpIntegralE(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1833)
    rule6794 = ReplacementRule(pattern6794, replacement6794)

    pattern6795 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1257, cons64)
    rule6795 = ReplacementRule(pattern6795, replacement6795)

    pattern6796 = Pattern(Integral(ExpIntegralE(S(1), x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6796 = ReplacementRule(pattern6796, replacement6796)

    pattern6797 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons1257, cons20, cons96)
    rule6797 = ReplacementRule(pattern6797, replacement6797)

    pattern6798 = Pattern(Integral(x_**m_*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons19, cons4, cons1257, cons21)
    rule6798 = ReplacementRule(pattern6798, replacement6798)

    pattern6799 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, x_*WC('b', S(1))), x_), cons3, cons19, cons4, cons1361)
    rule6799 = ReplacementRule(pattern6799, replacement6799)

    pattern6800 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons4, cons1968)
    rule6800 = ReplacementRule(pattern6800, replacement6800)

    pattern6801 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralE(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons1969, cons68)
    rule6801 = ReplacementRule(pattern6801, replacement6801)

    pattern6802 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6802 = ReplacementRule(pattern6802, replacement6802)

    pattern6803 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6803 = ReplacementRule(pattern6803, replacement6803)

    pattern6804 = Pattern(Integral(ExpIntegralEi(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6804 = ReplacementRule(pattern6804, replacement6804)

    pattern6805 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('b', S(1)))**S(2), x_), cons3, cons64)
    rule6805 = ReplacementRule(pattern6805, replacement6805)

    pattern6806 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6806 = ReplacementRule(pattern6806, replacement6806)

    pattern6807 = Pattern(Integral(ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6807 = ReplacementRule(pattern6807, replacement6807)

    pattern6808 = Pattern(Integral(x_**WC('m', S(1))*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6808 = ReplacementRule(pattern6808, replacement6808)

    pattern6809 = Pattern(Integral(x_**m_*ExpIntegralEi(x_*WC('d', S(1)) + WC('c', S(0)))*exp(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6809 = ReplacementRule(pattern6809, replacement6809)

    pattern6810 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6810 = ReplacementRule(pattern6810, replacement6810)

    pattern6811 = Pattern(Integral(LogIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6811 = ReplacementRule(pattern6811, replacement6811)

    pattern6812 = Pattern(Integral(x_**WC('m', S(1))*LogIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6812 = ReplacementRule(pattern6812, replacement6812)

    pattern6813 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6813 = ReplacementRule(pattern6813, replacement6813)

    pattern6814 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6814 = ReplacementRule(pattern6814, replacement6814)

    pattern6815 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6815 = ReplacementRule(pattern6815, replacement6815)

    pattern6816 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6816 = ReplacementRule(pattern6816, replacement6816)

    pattern6817 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6817 = ReplacementRule(pattern6817, replacement6817)

    pattern6818 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6818 = ReplacementRule(pattern6818, replacement6818)

    pattern6819 = Pattern(Integral(SinIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6819 = ReplacementRule(pattern6819, replacement6819)

    pattern6820 = Pattern(Integral(CosIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6820 = ReplacementRule(pattern6820, replacement6820)

    pattern6821 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons64)
    rule6821 = ReplacementRule(pattern6821, replacement6821)

    pattern6822 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons64)
    rule6822 = ReplacementRule(pattern6822, replacement6822)

    pattern6823 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6823 = ReplacementRule(pattern6823, replacement6823)

    pattern6824 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6824 = ReplacementRule(pattern6824, replacement6824)

    pattern6825 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6825 = ReplacementRule(pattern6825, replacement6825)

    pattern6826 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6826 = ReplacementRule(pattern6826, replacement6826)

    pattern6827 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6827 = ReplacementRule(pattern6827, replacement6827)

    pattern6828 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6828 = ReplacementRule(pattern6828, replacement6828)

    pattern6829 = Pattern(Integral(x_**m_*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6829 = ReplacementRule(pattern6829, replacement6829)

    pattern6830 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6830 = ReplacementRule(pattern6830, replacement6830)

    pattern6831 = Pattern(Integral(SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6831 = ReplacementRule(pattern6831, replacement6831)

    pattern6832 = Pattern(Integral(CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6832 = ReplacementRule(pattern6832, replacement6832)

    pattern6833 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6833 = ReplacementRule(pattern6833, replacement6833)

    pattern6834 = Pattern(Integral(x_**WC('m', S(1))*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6834 = ReplacementRule(pattern6834, replacement6834)

    pattern6835 = Pattern(Integral(x_**WC('m', S(1))*SinIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6835 = ReplacementRule(pattern6835, replacement6835)

    pattern6836 = Pattern(Integral(x_**m_*CosIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sin(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6836 = ReplacementRule(pattern6836, replacement6836)

    pattern6837 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6837 = ReplacementRule(pattern6837, replacement6837)

    pattern6838 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6838 = ReplacementRule(pattern6838, replacement6838)

    pattern6839 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6839 = ReplacementRule(pattern6839, replacement6839)

    pattern6840 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)))/x_, x_), cons3, cons3)
    rule6840 = ReplacementRule(pattern6840, replacement6840)

    pattern6841 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6841 = ReplacementRule(pattern6841, replacement6841)

    pattern6842 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons68)
    rule6842 = ReplacementRule(pattern6842, replacement6842)

    pattern6843 = Pattern(Integral(SinhIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6843 = ReplacementRule(pattern6843, replacement6843)

    pattern6844 = Pattern(Integral(CoshIntegral(x_*WC('b', S(1)) + WC('a', S(0)))**S(2), x_), cons2, cons3, cons69)
    rule6844 = ReplacementRule(pattern6844, replacement6844)

    pattern6845 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons64)
    rule6845 = ReplacementRule(pattern6845, replacement6845)

    pattern6846 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('b', S(1)))**S(2), x_), cons3, cons64)
    rule6846 = ReplacementRule(pattern6846, replacement6846)

    pattern6847 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6847 = ReplacementRule(pattern6847, replacement6847)

    pattern6848 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(a_ + x_*WC('b', S(1)))**S(2), x_), cons2, cons3, cons64)
    rule6848 = ReplacementRule(pattern6848, replacement6848)

    pattern6849 = Pattern(Integral(SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6849 = ReplacementRule(pattern6849, replacement6849)

    pattern6850 = Pattern(Integral(CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6850 = ReplacementRule(pattern6850, replacement6850)

    pattern6851 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons170)
    rule6851 = ReplacementRule(pattern6851, replacement6851)

    pattern6852 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons170)
    rule6852 = ReplacementRule(pattern6852, replacement6852)

    pattern6853 = Pattern(Integral(x_**m_*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6853 = ReplacementRule(pattern6853, replacement6853)

    pattern6854 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6854 = ReplacementRule(pattern6854, replacement6854)

    pattern6855 = Pattern(Integral(SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6855 = ReplacementRule(pattern6855, replacement6855)

    pattern6856 = Pattern(Integral(CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule6856 = ReplacementRule(pattern6856, replacement6856)

    pattern6857 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6857 = ReplacementRule(pattern6857, replacement6857)

    pattern6858 = Pattern(Integral(x_**WC('m', S(1))*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons64)
    rule6858 = ReplacementRule(pattern6858, replacement6858)

    pattern6859 = Pattern(Integral(x_**WC('m', S(1))*SinhIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*cosh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6859 = ReplacementRule(pattern6859, replacement6859)

    pattern6860 = Pattern(Integral(x_**m_*CoshIntegral(x_*WC('d', S(1)) + WC('c', S(0)))*sinh(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons20, cons96)
    rule6860 = ReplacementRule(pattern6860, replacement6860)

    pattern6861 = Pattern(Integral(Gamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6861 = ReplacementRule(pattern6861, replacement6861)

    pattern6862 = Pattern(Integral(Gamma(n_, b_*x_)/x_, x_), cons3, cons4, cons1970)
    rule6862 = ReplacementRule(pattern6862, replacement6862)

    pattern6863 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, b_*x_), x_), cons3, cons19, cons4, cons68)
    rule6863 = ReplacementRule(pattern6863, replacement6863)

    pattern6864 = Pattern(Integral(x_**WC('m', S(1))*Gamma(n_, a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons19, cons4, cons1971, cons68)
    rule6864 = ReplacementRule(pattern6864, With6864)

    pattern6865 = Pattern(Integral(LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6865 = ReplacementRule(pattern6865, replacement6865)

    pattern6866 = Pattern(Integral(x_**WC('m', S(1))*LogGamma(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons33, cons170)
    rule6866 = ReplacementRule(pattern6866, replacement6866)

    pattern6867 = Pattern(Integral(PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1833)
    rule6867 = ReplacementRule(pattern6867, replacement6867)

    pattern6868 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons33, cons170)
    rule6868 = ReplacementRule(pattern6868, replacement6868)

    pattern6869 = Pattern(Integral(x_**WC('m', S(1))*PolyGamma(n_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons33, cons96)
    rule6869 = ReplacementRule(pattern6869, replacement6869)

    pattern6870 = Pattern(Integral(Gamma(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons1833)
    rule6870 = ReplacementRule(pattern6870, replacement6870)

    pattern6871 = Pattern(Integral(Factorial(x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))*PolyGamma(S(0), x_*WC('b', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons4, cons1972)
    rule6871 = ReplacementRule(pattern6871, replacement6871)

    pattern6872 = Pattern(Integral(Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons69)
    rule6872 = ReplacementRule(pattern6872, replacement6872)

    pattern6873 = Pattern(Integral(Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons802, cons1973, cons1974)
    rule6873 = ReplacementRule(pattern6873, replacement6873)

    pattern6874 = Pattern(Integral(x_**WC('m', S(1))*Zeta(S(2), x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons33)
    rule6874 = ReplacementRule(pattern6874, replacement6874)

    pattern6875 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons802, cons1973, cons1974, cons33, cons170)
    rule6875 = ReplacementRule(pattern6875, replacement6875)

    pattern6876 = Pattern(Integral(x_**WC('m', S(1))*Zeta(s_, x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons802, cons1973, cons1974, cons33, cons96)
    rule6876 = ReplacementRule(pattern6876, replacement6876)

    pattern6877 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons52, cons89, cons90)
    rule6877 = ReplacementRule(pattern6877, replacement6877)

    pattern6878 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons5, cons52, cons89, cons91)
    rule6878 = ReplacementRule(pattern6878, replacement6878)

    pattern6879 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons385)
    rule6879 = ReplacementRule(pattern6879, replacement6879)

    pattern6880 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))/x_, x_), cons2, cons3, cons4, cons5, cons52, cons1975)
    rule6880 = ReplacementRule(pattern6880, replacement6880)

    pattern6881 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons19, cons5, cons52, cons68, cons89, cons90)
    rule6881 = ReplacementRule(pattern6881, replacement6881)

    pattern6882 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1))), x_), cons2, cons3, cons19, cons5, cons52, cons68, cons89, cons91)
    rule6882 = ReplacementRule(pattern6882, replacement6882)

    pattern6883 = Pattern(Integral(PolyLog(n_, (x_**WC('p', S(1))*WC('b', S(1)))**WC('q', S(1))*WC('a', S(1)))*log(x_**WC('m', S(1))*WC('c', S(1)))**WC('r', S(1))/x_, x_), cons2, cons3, cons8, cons19, cons4, cons52, cons54, cons1976, cons1977)
    rule6883 = ReplacementRule(pattern6883, replacement6883)

    pattern6884 = Pattern(Integral(PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons5, cons89, cons90)
    rule6884 = ReplacementRule(pattern6884, replacement6884)

    pattern6885 = Pattern(Integral(x_**WC('m', S(1))*PolyLog(n_, (x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('c', S(1))), x_), cons2, cons3, cons8, cons19, cons5, cons89, cons90, cons64)
    rule6885 = ReplacementRule(pattern6885, replacement6885)

    pattern6886 = Pattern(Integral(PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons4, cons5, cons1978)
    rule6886 = ReplacementRule(pattern6886, replacement6886)

    pattern6887 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*PolyLog(n_, (F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('p', S(1))*WC('d', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons33, cons170)
    rule6887 = ReplacementRule(pattern6887, replacement6887)

    pattern6888 = Pattern(Integral(u_*PolyLog(n_, v_), x_), cons4, cons4, CustomConstraint(With6888))
    rule6888 = ReplacementRule(pattern6888, replacement6888)

    pattern6889 = Pattern(Integral(u_*PolyLog(n_, v_)*log(w_), x_), cons4, cons1245, CustomConstraint(With6889))
    rule6889 = ReplacementRule(pattern6889, replacement6889)

    pattern6890 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons13, cons139)
    rule6890 = ReplacementRule(pattern6890, replacement6890)

    pattern6891 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons1381)
    rule6891 = ReplacementRule(pattern6891, replacement6891)

    pattern6892 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons5, cons64)
    rule6892 = ReplacementRule(pattern6892, replacement6892)

    pattern6893 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons4, cons5, cons1979)
    rule6893 = ReplacementRule(pattern6893, replacement6893)

    pattern6894 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons4, cons1980)
    rule6894 = ReplacementRule(pattern6894, replacement6894)

    pattern6895 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons5, cons198)
    rule6895 = ReplacementRule(pattern6895, replacement6895)

    pattern6896 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons19, cons4, cons5, cons68, cons1981)
    rule6896 = ReplacementRule(pattern6896, replacement6896)

    pattern6897 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons19, cons4, cons5, cons1982)
    rule6897 = ReplacementRule(pattern6897, replacement6897)

    pattern6898 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons19, cons1983)
    rule6898 = ReplacementRule(pattern6898, replacement6898)

    pattern6899 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons8, cons5, cons152, cons465, cons68)
    rule6899 = ReplacementRule(pattern6899, replacement6899)

    pattern6900 = Pattern(Integral(S(1)/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons29, cons1767)
    rule6900 = ReplacementRule(pattern6900, replacement6900)

    pattern6901 = Pattern(Integral(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons29, cons1767)
    rule6901 = ReplacementRule(pattern6901, replacement6901)

    pattern6902 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons13, cons165)
    rule6902 = ReplacementRule(pattern6902, replacement6902)

    pattern6903 = Pattern(Integral(S(1)/((d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))*ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons29, cons1767)
    rule6903 = ReplacementRule(pattern6903, replacement6903)

    pattern6904 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons950)
    rule6904 = ReplacementRule(pattern6904, replacement6904)

    pattern6905 = Pattern(Integral(S(1)/(sqrt(ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons951)
    rule6905 = ReplacementRule(pattern6905, replacement6905)

    pattern6906 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons13, cons139)
    rule6906 = ReplacementRule(pattern6906, replacement6906)

    pattern6907 = Pattern(Integral((ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('b', S(1)) + WC('a', S(0)))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons5, cons1984)
    rule6907 = ReplacementRule(pattern6907, replacement6907)

    pattern6908 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons29, cons64)
    rule6908 = ReplacementRule(pattern6908, replacement6908)

    pattern6909 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(a_ + x_*WC('b', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(a_ + x_*WC('b', S(1)))*WC('d', S(1))), x_), cons2, cons3, cons8, cons29, cons5, cons64)
    rule6909 = ReplacementRule(pattern6909, replacement6909)

    pattern6910 = Pattern(Integral(S(1)/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons198)
    rule6910 = ReplacementRule(pattern6910, replacement6910)

    pattern6911 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons4, cons5, cons1985)
    rule6911 = ReplacementRule(pattern6911, replacement6911)

    pattern6912 = Pattern(Integral(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons805, cons1986)
    rule6912 = ReplacementRule(pattern6912, replacement6912)

    pattern6913 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons805, cons1987, cons1988)
    rule6913 = ReplacementRule(pattern6913, replacement6913)

    pattern6914 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons805, cons1987, cons1989)
    rule6914 = ReplacementRule(pattern6914, replacement6914)

    pattern6915 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons340, cons90, cons1990)
    rule6915 = ReplacementRule(pattern6915, replacement6915)

    pattern6916 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons340, cons90, cons1991)
    rule6916 = ReplacementRule(pattern6916, replacement6916)

    pattern6917 = Pattern(Integral((ProductLog(x_**n_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons5, cons198)
    rule6917 = ReplacementRule(pattern6917, replacement6917)

    pattern6918 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons33, cons170)
    rule6918 = ReplacementRule(pattern6918, replacement6918)

    pattern6919 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons29, cons1992)
    rule6919 = ReplacementRule(pattern6919, replacement6919)

    pattern6920 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons33, cons96)
    rule6920 = ReplacementRule(pattern6920, replacement6920)

    pattern6921 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons19, cons21)
    rule6921 = ReplacementRule(pattern6921, replacement6921)

    pattern6922 = Pattern(Integral(S(1)/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons29, cons4, cons1993)
    rule6922 = ReplacementRule(pattern6922, replacement6922)

    pattern6923 = Pattern(Integral(x_**WC('m', S(1))/(d_ + ProductLog(x_**n_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons152, cons465, cons68)
    rule6923 = ReplacementRule(pattern6923, replacement6923)

    pattern6924 = Pattern(Integral((ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(x_*(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1)))), x_), cons2, cons8, cons29, cons4, cons5, cons1994)
    rule6924 = ReplacementRule(pattern6924, replacement6924)

    pattern6925 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons68, cons1995)
    rule6925 = ReplacementRule(pattern6925, replacement6925)

    pattern6926 = Pattern(Integral(x_**WC('m', S(1))*ProductLog(x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons29, cons19, cons4, cons40, cons1996)
    rule6926 = ReplacementRule(pattern6926, replacement6926)

    pattern6927 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons68, cons349, cons1997, cons1998)
    rule6927 = ReplacementRule(pattern6927, replacement6927)

    pattern6928 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**p_/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons68, cons349, cons1997, cons1999)
    rule6928 = ReplacementRule(pattern6928, replacement6928)

    pattern6929 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons68, cons2000, cons2001)
    rule6929 = ReplacementRule(pattern6929, replacement6929)

    pattern6930 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons4, cons5, cons68, cons2000, cons2002)
    rule6930 = ReplacementRule(pattern6930, replacement6930)

    pattern6931 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons19, cons5, cons68)
    rule6931 = ReplacementRule(pattern6931, replacement6931)

    pattern6932 = Pattern(Integral(x_**WC('m', S(1))*(ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('c', S(1)))**WC('p', S(1))/(d_ + ProductLog(x_**WC('n', S(1))*WC('a', S(1)))*WC('d', S(1))), x_), cons2, cons8, cons29, cons5, cons68, cons152, cons465)
    rule6932 = ReplacementRule(pattern6932, replacement6932)

    pattern6933 = Pattern(Integral(u_, x_), cons2003)
    rule6933 = ReplacementRule(pattern6933, replacement6933)
    return [rule6742, rule6743, rule6744, rule6745, rule6746, rule6747, rule6748, rule6749, rule6750, rule6751, rule6752, rule6753, rule6754, rule6755, rule6756, rule6757, rule6758, rule6759, rule6760, rule6761, rule6762, rule6763, rule6764, rule6765, rule6766, rule6767, rule6768, rule6769, rule6770, rule6771, rule6772, rule6773, rule6774, rule6775, rule6776, rule6777, rule6778, rule6779, rule6780, rule6781, rule6782, rule6783, rule6784, rule6785, rule6786, rule6787, rule6788, rule6789, rule6790, rule6791, rule6792, rule6793, rule6794, rule6795, rule6796, rule6797, rule6798, rule6799, rule6800, rule6801, rule6802, rule6803, rule6804, rule6805, rule6806, rule6807, rule6808, rule6809, rule6810, rule6811, rule6812, rule6813, rule6814, rule6815, rule6816, rule6817, rule6818, rule6819, rule6820, rule6821, rule6822, rule6823, rule6824, rule6825, rule6826, rule6827, rule6828, rule6829, rule6830, rule6831, rule6832, rule6833, rule6834, rule6835, rule6836, rule6837, rule6838, rule6839, rule6840, rule6841, rule6842, rule6843, rule6844, rule6845, rule6846, rule6847, rule6848, rule6849, rule6850, rule6851, rule6852, rule6853, rule6854, rule6855, rule6856, rule6857, rule6858, rule6859, rule6860, rule6861, rule6862, rule6863, rule6864, rule6865, rule6866, rule6867, rule6868, rule6869, rule6870, rule6871, rule6872, rule6873, rule6874, rule6875, rule6876, rule6877, rule6878, rule6879, rule6880, rule6881, rule6882, rule6883, rule6884, rule6885, rule6886, rule6887, rule6888, rule6889, rule6890, rule6891, rule6892, rule6893, rule6894, rule6895, rule6896, rule6897, rule6898, rule6899, rule6900, rule6901, rule6902, rule6903, rule6904, rule6905, rule6906, rule6907, rule6908, rule6909, rule6910, rule6911, rule6912, rule6913, rule6914, rule6915, rule6916, rule6917, rule6918, rule6919, rule6920, rule6921, rule6922, rule6923, rule6924, rule6925, rule6926, rule6927, rule6928, rule6929, rule6930, rule6931, rule6932, rule6933, ]





def replacement6742(a, b, x):
    return Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erf(a + b*x)/b, x)


def replacement6743(a, b, x):
    return -Simp(exp(-(a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfc(a + b*x)/b, x)


def replacement6744(a, b, x):
    return -Simp(exp((a + b*x)**S(2))/(sqrt(Pi)*b), x) + Simp((a + b*x)*Erfi(a + b*x)/b, x)


def replacement6745(b, x):
    return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -b**S(2)*x**S(2))/sqrt(Pi), x)


def replacement6746(b, x):
    return -Int(Erf(b*x)/x, x) + Simp(log(x), x)


def replacement6747(b, x):
    return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), b**S(2)*x**S(2))/sqrt(Pi), x)


def replacement6748(a, b, m, x):
    return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)/(m + S(1)), x)


def replacement6749(a, b, m, x):
    return Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-(a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)/(m + S(1)), x)


def replacement6750(a, b, m, x):
    return -Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp((a + b*x)**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)/(m + S(1)), x)


def replacement6751(a, b, c, d, x):
    return -Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6752(a, b, c, d, x):
    return Dist(b/(sqrt(Pi)*d), Int(exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6753(a, b, c, d, x):
    return -Dist(b/(sqrt(Pi)*d), Int(exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6754(a, b, c, d, m, x):
    return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erf(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6755(a, b, c, d, m, x):
    return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(-1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6756(a, b, c, d, m, x):
    return -Dist((m + S(-1))/(S(2)*d), Int(x**(m + S(-2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(b/(sqrt(Pi)*d), Int(x**(m + S(-1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(-1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(S(2)*d), x)


def replacement6757(b, c, d, x):
    return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)


def replacement6758(b, c, d, x):
    return Int(exp(c + d*x**S(2))/x, x) - Int(Erf(b*x)*exp(c + d*x**S(2))/x, x)


def replacement6759(b, c, d, x):
    return Simp(S(2)*b*x*HypergeometricPFQ(List(S(1)/2, S(1)), List(S(3)/2, S(3)/2), d*x**S(2))*exp(c)/sqrt(Pi), x)


def replacement6760(a, b, c, d, m, x):
    return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erf(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erf(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)


def replacement6761(a, b, c, d, m, x):
    return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfc(a + b*x)*exp(c + d*x**S(2)), x), x) + Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(-a**S(2) - S(2)*a*b*x + c - x**S(2)*(b**S(2) - d)), x), x) + Simp(x**(m + S(1))*Erfc(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)


def replacement6762(a, b, c, d, m, x):
    return -Dist(S(2)*d/(m + S(1)), Int(x**(m + S(2))*Erfi(a + b*x)*exp(c + d*x**S(2)), x), x) - Dist(S(2)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*exp(a**S(2) + S(2)*a*b*x + c + x**S(2)*(b**S(2) + d)), x), x) + Simp(x**(m + S(1))*Erfi(a + b*x)*exp(c + d*x**S(2))/(m + S(1)), x)


def replacement6763(a, b, x):
    return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erf(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erf(a + b*x)**S(2)/b, x)


def replacement6764(a, b, x):
    return Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfc(a + b*x)*exp(-(a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfc(a + b*x)**S(2)/b, x)


def replacement6765(a, b, x):
    return -Dist(S(4)/sqrt(Pi), Int((a + b*x)*Erfi(a + b*x)*exp((a + b*x)**S(2)), x), x) + Simp((a + b*x)*Erfi(a + b*x)**S(2)/b, x)


def replacement6766(b, m, x):
    return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erf(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erf(b*x)**S(2)/(m + S(1)), x)


def replacement6767(b, m, x):
    return Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfc(b*x)*exp(-b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfc(b*x)**S(2)/(m + S(1)), x)


def replacement6768(b, m, x):
    return -Dist(S(4)*b/(sqrt(Pi)*(m + S(1))), Int(x**(m + S(1))*Erfi(b*x)*exp(b**S(2)*x**S(2)), x), x) + Simp(x**(m + S(1))*Erfi(b*x)**S(2)/(m + S(1)), x)


def replacement6769(a, b, m, x):
    return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erf(x)**S(2), x), x, a + b*x), x)


def replacement6770(a, b, m, x):
    return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfc(x)**S(2), x), x, a + b*x), x)


def replacement6771(a, b, m, x):
    return Dist(S(1)/b, Subst(Int((-a/b + x/b)**m*Erfi(x)**S(2), x), x, a + b*x), x)


def replacement6772(a, b, x):
    return Simp(cos(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelS(a + b*x)/b, x)


def replacement6773(a, b, x):
    return -Simp(sin(Pi*(a + b*x)**S(2)/S(2))/(Pi*b), x) + Simp((a + b*x)*FresnelC(a + b*x)/b, x)


def replacement6774(b, x):
    return Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)


def replacement6775(b, x):
    return Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), -I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1)/2, S(1)/2), List(S(3)/2, S(3)/2), I*Pi*b**S(2)*x**S(2)/S(2))/S(2), x)


def replacement6776(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(a + b*x)/(m + S(1)), x)


def replacement6777(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(a + b*x)/(m + S(1)), x)


def replacement6778(a, b, x):
    return -Dist(S(2), Int((a + b*x)*FresnelS(a + b*x)*sin(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelS(a + b*x)**S(2)/b, x)


def replacement6779(a, b, x):
    return -Dist(S(2), Int((a + b*x)*FresnelC(a + b*x)*cos(Pi*(a + b*x)**S(2)/S(2)), x), x) + Simp((a + b*x)*FresnelC(a + b*x)**S(2)/b, x)


def replacement6780(b, m, x):
    return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)**S(2)/(m + S(1)), x)


def replacement6781(b, m, x):
    return -Dist(S(2)*b/(m + S(1)), Int(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)**S(2)/(m + S(1)), x)


def replacement6782(b, c, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) - Simp(FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6783(b, c, x):
    return -Dist(S(1)/(S(2)*Pi*b), Int(sin(Pi*b**S(2)*x**S(2)), x), x) + Simp(FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6784(b, c, m, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**(m + S(-1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6785(b, c, m, x):
    return -Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(-1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6786(b, c, m, x):
    return Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)


def replacement6787(b, c, m, x):
    return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(b*x**(m + S(2))/(S(2)*(m + S(1))*(m + S(2))), x) + Simp(x**(m + S(1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)


def replacement6788(b, c, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) - Simp(x/(S(2)*Pi*b), x) + Simp(FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6789(b, c, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(cos(Pi*b**S(2)*x**S(2)), x), x) + Simp(x/(S(2)*Pi*b), x) - Simp(FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6790(b, c, m, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) - Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) - Simp(x**m/(S(2)*Pi*b*m), x) + Simp(x**(m + S(-1))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6791(b, c, m, x):
    return Dist(S(1)/(S(2)*Pi*b), Int(x**(m + S(-1))*cos(Pi*b**S(2)*x**S(2)), x), x) + Dist((m + S(-1))/(Pi*b**S(2)), Int(x**(m + S(-2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**m/(S(2)*Pi*b*m), x) - Simp(x**(m + S(-1))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(Pi*b**S(2)), x)


def replacement6792(b, c, m, x):
    return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) + Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelS(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelS(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)


def replacement6793(b, c, m, x):
    return -Dist(b/(S(2)*m + S(2)), Int(x**(m + S(1))*sin(Pi*b**S(2)*x**S(2)), x), x) - Dist(Pi*b**S(2)/(m + S(1)), Int(x**(m + S(2))*FresnelC(b*x)*cos(Pi*b**S(2)*x**S(2)/S(2)), x), x) + Simp(x**(m + S(1))*FresnelC(b*x)*sin(Pi*b**S(2)*x**S(2)/S(2))/(m + S(1)), x)


def replacement6794(a, b, n, x):
    return -Simp(ExpIntegralE(n + S(1), a + b*x)/b, x)


def replacement6795(b, m, n, x):
    return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), b*x)/b, x)


def replacement6796(b, x):
    return -Simp(EulerGamma*log(x), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x), x) - Simp(log(b*x)**S(2)/S(2), x)


def replacement6797(b, m, n, x):
    return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + S(1)), x)


def replacement6798(b, m, n, x):
    return -Simp(x**(m + S(1))*HypergeometricPFQ(List(m + S(1), m + S(1)), List(m + S(2), m + S(2)), -b*x)/(m + S(1))**S(2), x) + Simp(x**m*(b*x)**(-m)*Gamma(m + S(1))*log(x)/b, x)


def replacement6799(b, m, n, x):
    return -Simp(x**(m + S(1))*ExpIntegralE(-m, b*x)/(m + n), x) + Simp(x**(m + S(1))*ExpIntegralE(n, b*x)/(m + n), x)


def replacement6800(a, b, m, n, x):
    return Dist(m/b, Int(x**(m + S(-1))*ExpIntegralE(n + S(1), a + b*x), x), x) - Simp(x**m*ExpIntegralE(n + S(1), a + b*x)/b, x)


def replacement6801(a, b, m, n, x):
    return Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralE(n + S(-1), a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralE(n, a + b*x)/(m + S(1)), x)


def replacement6802(a, b, x):
    return -Simp(exp(a + b*x)/b, x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)/b, x)


def replacement6803(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*exp(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)/(m + S(1)), x)


def replacement6804(a, b, x):
    return -Dist(S(2), Int(ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp((a + b*x)*ExpIntegralEi(a + b*x)**S(2)/b, x)


def replacement6805(b, m, x):
    return -Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(b*x)*exp(b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(b*x)**S(2)/(m + S(1)), x)


def replacement6806(a, b, m, x):
    return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*ExpIntegralEi(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*ExpIntegralEi(a + b*x)*exp(a + b*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*ExpIntegralEi(a + b*x)**S(2)/(b*(m + S(1))), x)


def replacement6807(a, b, c, d, x):
    return -Dist(d/b, Int(exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)


def replacement6808(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*exp(a + c + x*(b + d))/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) + Simp(x**m*ExpIntegralEi(c + d*x)*exp(a + b*x)/b, x)


def replacement6809(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*exp(a + c + x*(b + d))/(c + d*x), x), x) + Simp(x**(m + S(1))*ExpIntegralEi(c + d*x)*exp(a + b*x)/(m + S(1)), x)


def replacement6810(a, b, x):
    return -Simp(ExpIntegralEi(S(2)*log(a + b*x))/b, x) + Simp((a + b*x)*LogIntegral(a + b*x)/b, x)


def replacement6811(b, x):
    return -Simp(b*x, x) + Simp(LogIntegral(b*x)*log(b*x), x)


def replacement6812(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))/log(a + b*x), x), x) + Simp(x**(m + S(1))*LogIntegral(a + b*x)/(m + S(1)), x)


def replacement6813(a, b, x):
    return Simp(cos(a + b*x)/b, x) + Simp((a + b*x)*SinIntegral(a + b*x)/b, x)


def replacement6814(a, b, x):
    return -Simp(sin(a + b*x)/b, x) + Simp((a + b*x)*CosIntegral(a + b*x)/b, x)


def replacement6815(b, x):
    return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x)


def replacement6816(b, x):
    return Simp(EulerGamma*log(x), x) - Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -I*b*x)/S(2), x) + Simp(I*b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), I*b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)


def replacement6817(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)/(m + S(1)), x)


def replacement6818(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)/(m + S(1)), x)


def replacement6819(a, b, x):
    return -Dist(S(2), Int(SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp((a + b*x)*SinIntegral(a + b*x)**S(2)/b, x)


def replacement6820(a, b, x):
    return -Dist(S(2), Int(CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp((a + b*x)*CosIntegral(a + b*x)**S(2)/b, x)


def replacement6821(b, m, x):
    return -Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(b*x)*sin(b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(b*x)**S(2)/(m + S(1)), x)


def replacement6822(b, m, x):
    return -Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(b*x)*cos(b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(b*x)**S(2)/(m + S(1)), x)


def replacement6823(a, b, m, x):
    return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinIntegral(a + b*x)*sin(a + b*x), x), x) + Simp(x**(m + S(1))*SinIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)


def replacement6824(a, b, m, x):
    return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CosIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*CosIntegral(a + b*x)*cos(a + b*x), x), x) + Simp(x**(m + S(1))*CosIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CosIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)


def replacement6825(a, b, c, d, x):
    return Dist(d/b, Int(sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) - Simp(SinIntegral(c + d*x)*cos(a + b*x)/b, x)


def replacement6826(a, b, c, d, x):
    return -Dist(d/b, Int(sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(CosIntegral(c + d*x)*sin(a + b*x)/b, x)


def replacement6827(a, b, c, d, m, x):
    return Dist(d/b, Int(x**m*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*SinIntegral(c + d*x)*cos(a + b*x)/b, x)


def replacement6828(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*CosIntegral(c + d*x)*sin(a + b*x)/b, x)


def replacement6829(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)


def replacement6830(a, b, c, d, m, x):
    return Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)


def replacement6831(a, b, c, d, x):
    return -Dist(d/b, Int(sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) + Simp(SinIntegral(c + d*x)*sin(a + b*x)/b, x)


def replacement6832(a, b, c, d, x):
    return Dist(d/b, Int(cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) - Simp(CosIntegral(c + d*x)*cos(a + b*x)/b, x)


def replacement6833(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*sin(a + b*x)*sin(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) + Simp(x**m*SinIntegral(c + d*x)*sin(a + b*x)/b, x)


def replacement6834(a, b, c, d, m, x):
    return Dist(d/b, Int(x**m*cos(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Dist(m/b, Int(x**(m + S(-1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Simp(x**m*CosIntegral(c + d*x)*cos(a + b*x)/b, x)


def replacement6835(a, b, c, d, m, x):
    return Dist(b/(m + S(1)), Int(x**(m + S(1))*SinIntegral(c + d*x)*sin(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(c + d*x)*cos(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinIntegral(c + d*x)*cos(a + b*x)/(m + S(1)), x)


def replacement6836(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CosIntegral(c + d*x)*cos(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sin(a + b*x)*cos(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CosIntegral(c + d*x)*sin(a + b*x)/(m + S(1)), x)


def replacement6837(a, b, x):
    return -Simp(cosh(a + b*x)/b, x) + Simp((a + b*x)*SinhIntegral(a + b*x)/b, x)


def replacement6838(a, b, x):
    return -Simp(sinh(a + b*x)/b, x) + Simp((a + b*x)*CoshIntegral(a + b*x)/b, x)


def replacement6839(b, x):
    return Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x)


def replacement6840(b, x):
    return Simp(EulerGamma*log(x), x) - Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), -b*x)/S(2), x) + Simp(b*x*HypergeometricPFQ(List(S(1), S(1), S(1)), List(S(2), S(2), S(2)), b*x)/S(2), x) + Simp(log(b*x)**S(2)/S(2), x)


def replacement6841(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)/(m + S(1)), x)


def replacement6842(a, b, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*cosh(a + b*x)/(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)/(m + S(1)), x)


def replacement6843(a, b, x):
    return -Dist(S(2), Int(SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp((a + b*x)*SinhIntegral(a + b*x)**S(2)/b, x)


def replacement6844(a, b, x):
    return -Dist(S(2), Int(CoshIntegral(a + b*x)*cosh(a + b*x), x), x) + Simp((a + b*x)*CoshIntegral(a + b*x)**S(2)/b, x)


def replacement6845(b, m, x):
    return -Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(b*x)*sinh(b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(b*x)**S(2)/(m + S(1)), x)


def replacement6846(b, m, x):
    return -Dist(S(2)/(m + S(1)), Int(x**m*CoshIntegral(b*x)*cosh(b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(b*x)**S(2)/(m + S(1)), x)


def replacement6847(a, b, m, x):
    return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*SinhIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*SinhIntegral(a + b*x)*sinh(a + b*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*SinhIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)


def replacement6848(a, b, m, x):
    return -Dist(a*m/(b*(m + S(1))), Int(x**(m + S(-1))*CoshIntegral(a + b*x)**S(2), x), x) - Dist(S(2)/(m + S(1)), Int(x**m*CoshIntegral(a + b*x)*cosh(a + b*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(a + b*x)**S(2)/(m + S(1)), x) + Simp(a*x**m*CoshIntegral(a + b*x)**S(2)/(b*(m + S(1))), x)


def replacement6849(a, b, c, d, x):
    return -Dist(d/b, Int(sinh(c + d*x)*cosh(a + b*x)/(c + d*x), x), x) + Simp(SinhIntegral(c + d*x)*cosh(a + b*x)/b, x)


def replacement6850(a, b, c, d, x):
    return -Dist(d/b, Int(sinh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) + Simp(CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)


def replacement6851(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*sinh(c + d*x)*cosh(a + b*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinhIntegral(c + d*x)*cosh(a + b*x), x), x) + Simp(x**m*SinhIntegral(c + d*x)*cosh(a + b*x)/b, x)


def replacement6852(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*sinh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*CoshIntegral(c + d*x)*sinh(a + b*x)/b, x)


def replacement6853(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinhIntegral(c + d*x)*cosh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)


def replacement6854(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*cosh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(c + d*x)*cosh(a + b*x)/(m + S(1)), x)


def replacement6855(a, b, c, d, x):
    return -Dist(d/b, Int(sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) + Simp(SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)


def replacement6856(a, b, c, d, x):
    return -Dist(d/b, Int(cosh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) + Simp(CoshIntegral(c + d*x)*cosh(a + b*x)/b, x)


def replacement6857(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*sinh(a + b*x)*sinh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) + Simp(x**m*SinhIntegral(c + d*x)*sinh(a + b*x)/b, x)


def replacement6858(a, b, c, d, m, x):
    return -Dist(d/b, Int(x**m*cosh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) - Dist(m/b, Int(x**(m + S(-1))*CoshIntegral(c + d*x)*cosh(a + b*x), x), x) + Simp(x**m*CoshIntegral(c + d*x)*cosh(a + b*x)/b, x)


def replacement6859(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*SinhIntegral(c + d*x)*sinh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sinh(c + d*x)*cosh(a + b*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*SinhIntegral(c + d*x)*cosh(a + b*x)/(m + S(1)), x)


def replacement6860(a, b, c, d, m, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*CoshIntegral(c + d*x)*cosh(a + b*x), x), x) - Dist(d/(m + S(1)), Int(x**(m + S(1))*sinh(a + b*x)*cosh(c + d*x)/(c + d*x), x), x) + Simp(x**(m + S(1))*CoshIntegral(c + d*x)*sinh(a + b*x)/(m + S(1)), x)


def replacement6861(a, b, n, x):
    return -Simp(Gamma(n + S(1), a + b*x)/b, x) + Simp((a + b*x)*Gamma(n, a + b*x)/b, x)


def replacement6862(b, n, x):
    return Simp(Gamma(n)*log(x), x) - Simp((b*x)**n*HypergeometricPFQ(List(n, n), List(n + S(1), n + S(1)), -b*x)/n**S(2), x)


def replacement6863(b, m, n, x):
    return Simp(x**(m + S(1))*Gamma(n, b*x)/(m + S(1)), x) - Simp(x**m*(b*x)**(-m)*Gamma(m + n + S(1), b*x)/(b*(m + S(1))), x)


def With6864(a, b, m, n, x):
    _UseGamma = True
    return Dist(b/(m + S(1)), Int(x**(m + S(1))*(a + b*x)**(n + S(-1))*exp(-a - b*x), x), x) + Simp(x**(m + S(1))*Gamma(n, a + b*x)/(m + S(1)), x)


def replacement6865(a, b, x):
    return Simp(PolyGamma(S(-2), a + b*x)/b, x)


def replacement6866(a, b, m, x):
    return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(S(-2), a + b*x), x), x) + Simp(x**m*PolyGamma(S(-2), a + b*x)/b, x)


def replacement6867(a, b, n, x):
    return Simp(PolyGamma(n + S(-1), a + b*x)/b, x)


def replacement6868(a, b, m, n, x):
    return -Dist(m/b, Int(x**(m + S(-1))*PolyGamma(n + S(-1), a + b*x), x), x) + Simp(x**m*PolyGamma(n + S(-1), a + b*x)/b, x)


def replacement6869(a, b, m, n, x):
    return -Dist(b/(m + S(1)), Int(x**(m + S(1))*PolyGamma(n + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*PolyGamma(n, a + b*x)/(m + S(1)), x)


def replacement6870(a, b, n, x):
    return Simp(Gamma(a + b*x)**n/(b*n), x)


def replacement6871(a, b, c, n, x):
    return Simp(Factorial(a + b*x)**n/(b*n), x)


def replacement6872(a, b, x):
    return Int(PolyGamma(S(1), a + b*x), x)


def replacement6873(a, b, s, x):
    return -Simp(Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)


def replacement6874(a, b, m, x):
    return Int(x**m*PolyGamma(S(1), a + b*x), x)


def replacement6875(a, b, m, s, x):
    return Dist(m/(b*(s + S(-1))), Int(x**(m + S(-1))*Zeta(s + S(-1), a + b*x), x), x) - Simp(x**m*Zeta(s + S(-1), a + b*x)/(b*(s + S(-1))), x)


def replacement6876(a, b, m, s, x):
    return Dist(b*s/(m + S(1)), Int(x**(m + S(1))*Zeta(s + S(1), a + b*x), x), x) + Simp(x**(m + S(1))*Zeta(s, a + b*x)/(m + S(1)), x)


def replacement6877(a, b, n, p, q, x):
    return -Dist(p*q, Int(PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n, a*(b*x**p)**q), x)


def replacement6878(a, b, n, p, q, x):
    return -Dist(S(1)/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)


def replacement6879(a, b, c, d, e, n, p, x):
    return Simp(PolyLog(n + S(1), c*(a + b*x)**p)/(e*p), x)


def replacement6880(a, b, n, p, q, x):
    return Simp(PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)


def replacement6881(a, b, m, n, p, q, x):
    return -Dist(p*q/(m + S(1)), Int(x**m*PolyLog(n + S(-1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n, a*(b*x**p)**q)/(m + S(1)), x)


def replacement6882(a, b, m, n, p, q, x):
    return -Dist((m + S(1))/(p*q), Int(x**m*PolyLog(n + S(1), a*(b*x**p)**q), x), x) + Simp(x**(m + S(1))*PolyLog(n + S(1), a*(b*x**p)**q)/(p*q), x)


def replacement6883(a, b, c, m, n, p, q, r, x):
    return -Dist(m*r/(p*q), Int(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**(r + S(-1))/x, x), x) + Simp(PolyLog(n + S(1), a*(b*x**p)**q)*log(c*x**m)**r/(p*q), x)


def replacement6884(a, b, c, n, p, x):
    return -Dist(p, Int(PolyLog(n + S(-1), c*(a + b*x)**p), x), x) + Dist(a*p, Int(PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x*PolyLog(n, c*(a + b*x)**p), x)


def replacement6885(a, b, c, m, n, p, x):
    return -Dist(b*p/(m + S(1)), Int(x**(m + S(1))*PolyLog(n + S(-1), c*(a + b*x)**p)/(a + b*x), x), x) + Simp(x**(m + S(1))*PolyLog(n, c*(a + b*x)**p)/(m + S(1)), x)


def replacement6886(F, a, b, c, d, n, p, x):
    return Simp(PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)


def replacement6887(F, a, b, c, d, e, f, m, n, p, x):
    return -Dist(f*m/(b*c*p*log(F)), Int((e + f*x)**(m + S(-1))*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p), x), x) + Simp((e + f*x)**m*PolyLog(n + S(1), d*(F**(c*(a + b*x)))**p)/(b*c*p*log(F)), x)


def With6888(n, u, v, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        w = DerivativeDivides(v, u*v, x)
        res = Not(FalseQ(w))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6888(n, u, v, x):

    w = DerivativeDivides(v, u*v, x)
    return Simp(w*PolyLog(n + S(1), v), x)


def With6889(n, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    try:
        z = DerivativeDivides(v, u*v, x)
        res = Not(FalseQ(z))
    except (TypeError, AttributeError):
        return False
    if res:
        return True
    return False


def replacement6889(n, u, v, w, x):

    z = DerivativeDivides(v, u*v, x)
    return -Int(SimplifyIntegrand(z*D(w, x)*PolyLog(n + S(1), v)/w, x), x) + Simp(z*PolyLog(n + S(1), v)*log(w), x)


def replacement6890(a, b, c, p, x):
    return Dist(p/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*(p + S(1))), x)


def replacement6891(a, b, c, p, x):
    return -Dist(p, Int((c*ProductLog(a + b*x))**p/(ProductLog(a + b*x) + S(1)), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/b, x)


def replacement6892(a, b, c, m, p, x):
    return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p, (-a/b + x/b)**m, x), x), x, a + b*x), x)


def replacement6893(a, c, n, p, x):
    return -Dist(n*p, Int((c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p, x)


def replacement6894(a, c, n, p, x):
    return Dist(n*p/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(n*p + S(1)), x)


def replacement6895(a, c, n, p, x):
    return -Subst(Int((c*ProductLog(a*x**(-n)))**p/x**S(2), x), x, S(1)/x)


def replacement6896(a, c, m, n, p, x):
    return -Dist(n*p/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**p/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + S(1)), x)


def replacement6897(a, c, m, n, p, x):
    return Dist(n*p/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(ProductLog(a*x**n) + S(1)), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(m + n*p + S(1)), x)


def replacement6898(a, c, m, p, x):
    return Dist(S(1)/c, Int(x**m*(c*ProductLog(a*x))**(p + S(1))/(ProductLog(a*x) + S(1)), x), x) + Int(x**m*(c*ProductLog(a*x))**p/(ProductLog(a*x) + S(1)), x)


def replacement6899(a, c, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p, x), x, S(1)/x)


def replacement6900(a, b, d, x):
    return Simp((a + b*x)/(b*d*ProductLog(a + b*x)), x)


def replacement6901(a, b, d, x):
    return -Int(S(1)/(d*ProductLog(a + b*x) + d), x) + Simp(d*x, x)


def replacement6902(a, b, c, d, p, x):
    return -Dist(c*p, Int((c*ProductLog(a + b*x))**(p + S(-1))/(d*ProductLog(a + b*x) + d), x), x) + Simp(c*(c*ProductLog(a + b*x))**(p + S(-1))*(a + b*x)/(b*d), x)


def replacement6903(a, b, d, x):
    return Simp(ExpIntegralEi(ProductLog(a + b*x))/(b*d), x)


def replacement6904(a, b, c, d, x):
    return Simp(Erfi(sqrt(c*ProductLog(a + b*x))/Rt(c, S(2)))*Rt(Pi*c, S(2))/(b*c*d), x)


def replacement6905(a, b, c, d, x):
    return Simp(Erf(sqrt(c*ProductLog(a + b*x))/Rt(-c, S(2)))*Rt(-Pi*c, S(2))/(b*c*d), x)


def replacement6906(a, b, c, d, p, x):
    return -Dist(S(1)/(c*(p + S(1))), Int((c*ProductLog(a + b*x))**(p + S(1))/(d*ProductLog(a + b*x) + d), x), x) + Simp((c*ProductLog(a + b*x))**p*(a + b*x)/(b*d*(p + S(1))), x)


def replacement6907(a, b, c, d, p, x):
    return Simp((-ProductLog(a + b*x))**(-p)*(c*ProductLog(a + b*x))**p*Gamma(p + S(1), -ProductLog(a + b*x))/(b*d), x)


def replacement6908(a, b, d, m, x):
    return Dist(S(1)/b, Subst(Int(ExpandIntegrand(S(1)/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)


def replacement6909(a, b, c, d, m, p, x):
    return Dist(S(1)/b, Subst(Int(ExpandIntegrand((c*ProductLog(x))**p/(d*ProductLog(x) + d), (-a/b + x/b)**m, x), x), x, a + b*x), x)


def replacement6910(a, d, n, x):
    return -Subst(Int(S(1)/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)


def replacement6911(a, c, d, n, p, x):
    return Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)


def replacement6912(a, d, n, p, x):
    return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)


def replacement6913(a, c, d, n, p, x):
    return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(c*n, S(2)))*Rt(Pi*c*n, S(2))/(d*n), x)


def replacement6914(a, c, d, n, p, x):
    return Simp(a**(-S(1)/n)*c**(-S(1)/n)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(-c*n, S(2)))*Rt(-Pi*c*n, S(2))/(d*n), x)


def replacement6915(a, c, d, n, p, x):
    return -Dist(c*(n*(p + S(-1)) + S(1)), Int((c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x*(c*ProductLog(a*x**n))**(p + S(-1))/d, x)


def replacement6916(a, c, d, n, p, x):
    return -Dist(S(1)/(c*(n*p + S(1))), Int((c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x*(c*ProductLog(a*x**n))**p/(d*(n*p + S(1))), x)


def replacement6917(a, c, d, n, p, x):
    return -Subst(Int((c*ProductLog(a*x**(-n)))**p/(x**S(2)*(d*ProductLog(a*x**(-n)) + d)), x), x, S(1)/x)


def replacement6918(a, d, m, x):
    return -Dist(m/(m + S(1)), Int(x**m/((d*ProductLog(a*x) + d)*ProductLog(a*x)), x), x) + Simp(x**(m + S(1))/(d*(m + S(1))*ProductLog(a*x)), x)


def replacement6919(a, d, x):
    return Simp(log(ProductLog(a*x))/d, x)


def replacement6920(a, d, m, x):
    return -Int(x**m*ProductLog(a*x)/(d*ProductLog(a*x) + d), x) + Simp(x**(m + S(1))/(d*(m + S(1))), x)


def replacement6921(a, d, m, x):
    return Simp(x**m*(-(m + S(1))*ProductLog(a*x))**(-m)*Gamma(m + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)


def replacement6922(a, d, n, x):
    return Simp(log(ProductLog(a*x**n))/(d*n), x)


def replacement6923(a, d, m, n, x):
    return -Subst(Int(x**(-m + S(-2))/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)


def replacement6924(a, c, d, n, p, x):
    return Simp((c*ProductLog(a*x**n))**p/(d*n*p), x)


def replacement6925(a, c, d, m, n, p, x):
    return Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)


def replacement6926(a, d, m, n, p, x):
    return Simp(a**p*ExpIntegralEi(-p*ProductLog(a*x**n))/(d*n), x)


def replacement6927(a, c, d, m, n, p, x):
    return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erf(sqrt(c*ProductLog(a*x**n))/Rt(c/(p + S(-1)/2), S(2)))*Rt(Pi*c/(p + S(-1)/2), S(2))/(d*n), x)


def replacement6928(a, c, d, m, n, p, x):
    return Simp(a**(p + S(-1)/2)*c**(p + S(-1)/2)*Erfi(sqrt(c*ProductLog(a*x**n))/Rt(-c/(p + S(-1)/2), S(2)))*Rt(-Pi*c/(p + S(-1)/2), S(2))/(d*n), x)


def replacement6929(a, c, d, m, n, p, x):
    return -Dist(c*(m + n*(p + S(-1)) + S(1))/(m + S(1)), Int(x**m*(c*ProductLog(a*x**n))**(p + S(-1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(c*x**(m + S(1))*(c*ProductLog(a*x**n))**(p + S(-1))/(d*(m + S(1))), x)


def replacement6930(a, c, d, m, n, p, x):
    return -Dist((m + S(1))/(c*(m + n*p + S(1))), Int(x**m*(c*ProductLog(a*x**n))**(p + S(1))/(d*ProductLog(a*x**n) + d), x), x) + Simp(x**(m + S(1))*(c*ProductLog(a*x**n))**p/(d*(m + n*p + S(1))), x)


def replacement6931(a, c, d, m, p, x):
    return Simp(x**m*(c*ProductLog(a*x))**p*(-(m + S(1))*ProductLog(a*x))**(-m - p)*Gamma(m + p + S(1), -(m + S(1))*ProductLog(a*x))*exp(-m*ProductLog(a*x))/(a*d*(m + S(1))), x)


def replacement6932(a, c, d, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(c*ProductLog(a*x**(-n)))**p/(d*ProductLog(a*x**(-n)) + d), x), x, S(1)/x)


def replacement6933(u, x):
    return Subst(Int(SimplifyIntegrand((x + S(1))*SubstFor(ProductLog(x), u, x)*exp(x), x), x), x, ProductLog(x))
