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


def exponential():
    from sympy.integrals.rubi.constraints import cons33, cons170, cons517, cons1100, cons1101, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons96, cons20, cons21, cons19, cons1102, cons130, cons2, cons246, cons139, cons554, cons1103, cons1104, cons5, cons382, cons56, cons1105, cons1106, cons1107, cons211, cons226, cons798, cons799, cons52, cons1108, cons806, cons1109, cons814, cons1110, cons1111, cons1112, cons1113, cons586, cons1114, cons1115, cons481, cons482, cons1116, cons198, cons25, cons1117, cons55, cons1118, cons1119, cons1120, cons1121, cons87, cons1122, cons358, cons533, cons1123, cons1124, cons537, cons95, cons1125, cons1126, cons178, cons369, cons168, cons746, cons70, cons842, cons1127, cons1128, cons1129, cons27, cons73, cons1130, cons1131, cons1132, cons820, cons1133, cons1134, cons1135, cons1136, cons821, cons1137, cons1138, cons1139, cons1140, cons150, cons812, cons813, cons1141, cons1142, cons54, cons802, cons1143, cons1144, cons1145, cons815, cons1146, cons228, cons64, cons1147, cons1148, cons1149, cons1150, cons1151, cons1152, cons1153, cons465, cons1154, cons45, cons450, cons1155, cons1156, cons1157, cons1019


    pattern1904 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons33, cons170, cons517, cons1100)
    rule1904 = ReplacementRule(pattern1904, replacement1904)

    pattern1905 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons1101, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons33, cons96, cons517, cons1100)
    rule1905 = ReplacementRule(pattern1905, replacement1905)

    pattern1906 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons1101, cons8, cons29, cons50, cons127, cons210, cons1100)
    rule1906 = ReplacementRule(pattern1906, replacement1906)

    pattern1907 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons8, cons29, cons50, cons127, cons210, cons20)
    rule1907 = ReplacementRule(pattern1907, replacement1907)

    pattern1908 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons1101, cons8, cons29, cons50, cons127, cons210, cons1100)
    rule1908 = ReplacementRule(pattern1908, replacement1908)

    pattern1909 = Pattern(Integral(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons1101, cons8, cons29, cons50, cons127, cons210, cons19, cons21)
    rule1909 = ReplacementRule(pattern1909, replacement1909)

    pattern1910 = Pattern(Integral((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1)))*WC('b', S(1)))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons1102)
    rule1910 = ReplacementRule(pattern1910, replacement1910)

    pattern1911 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons130)
    rule1911 = ReplacementRule(pattern1911, replacement1911)

    pattern1912 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons33, cons170)
    rule1912 = ReplacementRule(pattern1912, replacement1912)

    pattern1913 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**p_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons246, cons170, cons139)
    rule1913 = ReplacementRule(pattern1913, With1913)

    pattern1914 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1101, cons2, cons3, cons210, cons4, cons5, cons554, cons1103, cons1104, cons20)
    rule1914 = ReplacementRule(pattern1914, replacement1914)

    pattern1915 = Pattern(Integral(u_**WC('m', S(1))*((F_**(v_*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1101, cons2, cons3, cons210, cons19, cons4, cons5, cons554, cons1103, cons1104, cons21)
    rule1915 = ReplacementRule(pattern1915, With1915)

    pattern1916 = Pattern(Integral((a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule1916 = ReplacementRule(pattern1916, replacement1916)

    pattern1917 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))/(a_ + (F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons33, cons170)
    rule1917 = ReplacementRule(pattern1917, replacement1917)

    pattern1918 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons56)
    rule1918 = ReplacementRule(pattern1918, replacement1918)

    pattern1919 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons1105)
    rule1919 = ReplacementRule(pattern1919, replacement1919)

    pattern1920 = Pattern(Integral((G_**((x_*WC('i', S(1)) + WC('h', S(0)))*WC('j', S(1)))*WC('k', S(1)))**WC('q', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*((F_**((x_*WC('f', S(1)) + WC('e', S(0)))*WC('g', S(1))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons798, cons799, cons19, cons4, cons5, cons52, cons1106, cons1107)
    rule1920 = ReplacementRule(pattern1920, replacement1920)

    pattern1921 = Pattern(Integral((F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1))))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons4, cons1108)
    rule1921 = ReplacementRule(pattern1921, replacement1921)

    pattern1922 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), cons1101, cons8, cons806, cons554, cons1109)
    rule1922 = ReplacementRule(pattern1922, replacement1922)

    pattern1923 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_, x_), cons1101, cons8, cons806, cons554, cons1100)
    rule1923 = ReplacementRule(pattern1923, replacement1923)

    pattern1924 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1101, cons8, cons19, cons814, cons1110)
    rule1924 = ReplacementRule(pattern1924, replacement1924)

    pattern1925 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1101, cons8, cons1111, cons554, cons1103, cons20, cons1109)
    rule1925 = ReplacementRule(pattern1925, replacement1925)

    pattern1926 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1101, cons8, cons1111, cons554, cons1103, cons20, cons1100)
    rule1926 = ReplacementRule(pattern1926, replacement1926)

    pattern1927 = Pattern(Integral(F_**(v_*WC('c', S(1)))*u_**WC('m', S(1))*w_, x_), cons1101, cons8, cons19, cons1111, cons554, cons1103, cons21)
    rule1927 = ReplacementRule(pattern1927, With1927)

    pattern1928 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1112, cons1113, cons586)
    rule1928 = ReplacementRule(pattern1928, replacement1928)

    pattern1929 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('c', S(1)))*x_**WC('m', S(1))*(e_ + (x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1))*log(x_*WC('d', S(1))))*log(x_*WC('d', S(1)))**WC('n', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons1114, cons1113, cons586)
    rule1929 = ReplacementRule(pattern1929, replacement1929)

    pattern1930 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons1115)
    rule1930 = ReplacementRule(pattern1930, replacement1930)

    pattern1931 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons481)
    rule1931 = ReplacementRule(pattern1931, replacement1931)

    pattern1932 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons482)
    rule1932 = ReplacementRule(pattern1932, replacement1932)

    pattern1933 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons1116, cons198)
    rule1933 = ReplacementRule(pattern1933, replacement1933)

    pattern1934 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons1116, cons25)
    rule1934 = ReplacementRule(pattern1934, With1934)

    pattern1935 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons4, cons1117)
    rule1935 = ReplacementRule(pattern1935, replacement1935)

    pattern1936 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons55, cons1118)
    rule1936 = ReplacementRule(pattern1936, replacement1936)

    pattern1937 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1118)
    rule1937 = ReplacementRule(pattern1937, replacement1937)

    pattern1938 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons19, cons4, cons1119)
    rule1938 = ReplacementRule(pattern1938, replacement1938)

    pattern1939 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons33, cons1120, cons1121, cons87, cons1122)
    rule1939 = ReplacementRule(pattern1939, replacement1939)

    pattern1940 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons19, cons4, cons1120, cons1121, cons358, cons533)
    rule1940 = ReplacementRule(pattern1940, replacement1940)

    pattern1941 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons33, cons1120, cons1123, cons87, cons1124)
    rule1941 = ReplacementRule(pattern1941, replacement1941)

    pattern1942 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons19, cons4, cons1120, cons1123, cons358, cons537)
    rule1942 = ReplacementRule(pattern1942, replacement1942)

    pattern1943 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons95, cons1120, cons1121, cons25)
    rule1943 = ReplacementRule(pattern1943, With1943)

    pattern1944 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1118, cons1120, cons1125, cons21, cons1126)
    rule1944 = ReplacementRule(pattern1944, replacement1944)

    pattern1945 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1118)
    rule1945 = ReplacementRule(pattern1945, replacement1945)

    pattern1946 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons178, cons369, cons168)
    rule1946 = ReplacementRule(pattern1946, replacement1946)

    pattern1947 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**S(2)*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons178, cons33, cons96)
    rule1947 = ReplacementRule(pattern1947, replacement1947)

    pattern1948 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons178, cons87, cons746, cons33, cons96)
    rule1948 = ReplacementRule(pattern1948, replacement1948)

    pattern1949 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons178)
    rule1949 = ReplacementRule(pattern1949, replacement1949)

    pattern1950 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))*(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons178, cons20, cons96)
    rule1950 = ReplacementRule(pattern1950, replacement1950)

    pattern1951 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons178)
    rule1951 = ReplacementRule(pattern1951, replacement1951)

    pattern1952 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), cons1101, cons19, cons70, cons842, cons1127)
    rule1952 = ReplacementRule(pattern1952, replacement1952)

    pattern1953 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))**n_*WC('b', S(1)) + WC('a', S(0)))*u_, x_), cons1101, cons2, cons3, cons8, cons29, cons4, cons806)
    rule1953 = ReplacementRule(pattern1953, replacement1953)

    pattern1954 = Pattern(Integral(F_**(v_*WC('b', S(1)) + WC('a', S(0)))*WC('u', S(1)), x_), cons1101, cons2, cons3, cons806, cons1128, cons1129)
    rule1954 = ReplacementRule(pattern1954, replacement1954)

    pattern1955 = Pattern(Integral(F_**(WC('a', S(0)) + WC('b', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))))/((x_*WC('f', S(1)) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons1118)
    rule1955 = ReplacementRule(pattern1955, replacement1955)

    pattern1956 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons27)
    rule1956 = ReplacementRule(pattern1956, replacement1956)

    pattern1957 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons73, cons1130)
    rule1957 = ReplacementRule(pattern1957, replacement1957)

    pattern1958 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons73, cons1131)
    rule1958 = ReplacementRule(pattern1958, replacement1958)

    pattern1959 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons73, cons1131, cons20, cons96)
    rule1959 = ReplacementRule(pattern1959, replacement1959)

    pattern1960 = Pattern(Integral(F_**((x_*WC('b', S(1)) + WC('a', S(0)))*WC('f', S(1))/(x_*WC('d', S(1)) + WC('c', S(0))) + WC('e', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_*WC('j', S(1)) + WC('i', S(0)))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons1130)
    rule1960 = ReplacementRule(pattern1960, replacement1960)

    pattern1961 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons1132)
    rule1961 = ReplacementRule(pattern1961, replacement1961)

    pattern1962 = Pattern(Integral(F_**v_, x_), cons1101, cons820, cons1133)
    rule1962 = ReplacementRule(pattern1962, replacement1962)

    pattern1963 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1134)
    rule1963 = ReplacementRule(pattern1963, replacement1963)

    pattern1964 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1134, cons33, cons168)
    rule1964 = ReplacementRule(pattern1964, replacement1964)

    pattern1965 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1134)
    rule1965 = ReplacementRule(pattern1965, replacement1965)

    pattern1966 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1134, cons33, cons96)
    rule1966 = ReplacementRule(pattern1966, replacement1966)

    pattern1967 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1135)
    rule1967 = ReplacementRule(pattern1967, replacement1967)

    pattern1968 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1135, cons33, cons168)
    rule1968 = ReplacementRule(pattern1968, replacement1968)

    pattern1969 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1135, cons33, cons96)
    rule1969 = ReplacementRule(pattern1969, replacement1969)

    pattern1970 = Pattern(Integral(F_**(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons19, cons1136)
    rule1970 = ReplacementRule(pattern1970, replacement1970)

    pattern1971 = Pattern(Integral(F_**v_*u_**WC('m', S(1)), x_), cons1101, cons19, cons70, cons820, cons821)
    rule1971 = ReplacementRule(pattern1971, replacement1971)

    pattern1972 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**v_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons1137, cons33, cons170, cons198)
    rule1972 = ReplacementRule(pattern1972, With1972)

    pattern1973 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1138, CustomConstraint(With1973))
    rule1973 = ReplacementRule(pattern1973, replacement1973)

    pattern1974 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1138, CustomConstraint(With1974))
    rule1974 = ReplacementRule(pattern1974, replacement1974)

    pattern1975 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons1140, cons150)
    rule1975 = ReplacementRule(pattern1975, replacement1975)

    pattern1976 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons1140, cons198)
    rule1976 = ReplacementRule(pattern1976, replacement1976)

    pattern1977 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons1140, cons25)
    rule1977 = ReplacementRule(pattern1977, replacement1977)

    pattern1978 = Pattern(Integral(G_**(u_*WC('h', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons2, cons3, cons50, cons211, cons4, cons812, cons813)
    rule1978 = ReplacementRule(pattern1978, replacement1978)

    pattern1979 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons1142, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons54, cons802, cons1143, cons4, cons1141, CustomConstraint(With1979))
    rule1979 = ReplacementRule(pattern1979, replacement1979)

    pattern1980 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons1142, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons54, cons802, cons1143, cons1144, cons87)
    rule1980 = ReplacementRule(pattern1980, replacement1980)

    pattern1981 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**WC('n', S(1)), x_), cons1101, cons1139, cons1142, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons54, cons802, cons1143, cons1145, cons150)
    rule1981 = ReplacementRule(pattern1981, replacement1981)

    pattern1982 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons1142, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons54, cons802, cons1143, cons1145, cons198)
    rule1982 = ReplacementRule(pattern1982, replacement1982)

    pattern1983 = Pattern(Integral(G_**((x_*WC('g', S(1)) + WC('f', S(0)))*WC('h', S(1)))*H_**((x_*WC('s', S(1)) + WC('r', S(0)))*WC('t', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons1142, cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons54, cons802, cons1143, cons4, cons1145, cons25)
    rule1983 = ReplacementRule(pattern1983, replacement1983)

    pattern1984 = Pattern(Integral(G_**(u_*WC('h', S(1)))*H_**(w_*WC('t', S(1)))*(F_**(v_*WC('e', S(1)))*WC('b', S(1)) + a_)**n_, x_), cons1101, cons1139, cons1142, cons2, cons3, cons50, cons211, cons1143, cons4, cons814, cons815)
    rule1984 = ReplacementRule(pattern1984, replacement1984)

    pattern1985 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons56)
    rule1985 = ReplacementRule(pattern1985, replacement1985)

    pattern1986 = Pattern(Integral(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*x_**WC('m', S(1))*(F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1)))*WC('b', S(1)) + x_**WC('n', S(1))*WC('a', S(1)))**WC('p', S(1)), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons56)
    rule1986 = ReplacementRule(pattern1986, replacement1986)

    pattern1987 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons127, cons210, cons1146, cons70, cons228, cons64)
    rule1987 = ReplacementRule(pattern1987, With1987)

    pattern1988 = Pattern(Integral(F_**u_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons127, cons210, cons1146, cons70, cons228, cons64)
    rule1988 = ReplacementRule(pattern1988, With1988)

    pattern1989 = Pattern(Integral((F_**u_*WC('i', S(1)) + h_)*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('m', S(1))/(F_**u_*WC('b', S(1)) + F_**v_*WC('c', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons127, cons210, cons211, cons226, cons1146, cons70, cons228, cons64)
    rule1989 = ReplacementRule(pattern1989, With1989)

    pattern1990 = Pattern(Integral(x_**WC('m', S(1))/(F_**v_*WC('b', S(1)) + F_**(x_*WC('d', S(1)) + WC('c', S(0)))*WC('a', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons1147, cons33, cons170)
    rule1990 = ReplacementRule(pattern1990, With1990)

    pattern1991 = Pattern(Integral(u_/(F_**v_*WC('b', S(1)) + F_**w_*WC('c', S(1)) + a_), x_), cons1101, cons2, cons3, cons8, cons554, cons1148, cons1149, cons1150)
    rule1991 = ReplacementRule(pattern1991, replacement1991)

    pattern1992 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons210, cons4, cons1151)
    rule1992 = ReplacementRule(pattern1992, replacement1992)

    pattern1993 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons1101, cons2, cons8, cons29, cons50, cons210, cons4, cons1152)
    rule1993 = ReplacementRule(pattern1993, replacement1993)

    pattern1994 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(c_*x_**S(2) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons210, cons4, cons806, cons20)
    rule1994 = ReplacementRule(pattern1994, replacement1994)

    pattern1995 = Pattern(Integral(F_**((x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1))*WC('g', S(1)))*u_**WC('m', S(1))/(a_ + c_*x_**S(2)), x_), cons1101, cons2, cons8, cons29, cons50, cons210, cons4, cons806, cons20)
    rule1995 = ReplacementRule(pattern1995, replacement1995)

    pattern1996 = Pattern(Integral(F_**((x_**S(4)*WC('b', S(1)) + WC('a', S(0)))/x_**S(2)), x_), cons1101, cons2, cons3, cons1153)
    rule1996 = ReplacementRule(pattern1996, replacement1996)

    pattern1997 = Pattern(Integral(x_**WC('m', S(1))*(x_**WC('m', S(1)) + exp(x_))**n_, x_), cons95, cons170, cons465, cons1154)
    rule1997 = ReplacementRule(pattern1997, replacement1997)

    pattern1998 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons45)
    rule1998 = ReplacementRule(pattern1998, replacement1998)

    pattern1999 = Pattern(Integral(log(a_ + (F_**((x_*WC('d', S(1)) + WC('c', S(0)))*WC('e', S(1))))**WC('n', S(1))*WC('b', S(1))), x_), cons1101, cons2, cons3, cons8, cons29, cons50, cons4, cons450)
    rule1999 = ReplacementRule(pattern1999, replacement1999)

    pattern2000 = Pattern(Integral((F_**v_*WC('a', S(1)))**n_*WC('u', S(1)), x_), cons1101, cons2, cons4, cons25)
    rule2000 = ReplacementRule(pattern2000, replacement2000)

    pattern2001 = Pattern(Integral(u_, x_), cons1155)
    rule2001 = ReplacementRule(pattern2001, With2001)

    pattern2002 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1101, cons2, cons3, cons4, cons198, cons1156)
    rule2002 = ReplacementRule(pattern2002, replacement2002)

    pattern2003 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1101, cons1139, cons2, cons3, cons4, cons198, cons1156)
    rule2003 = ReplacementRule(pattern2003, replacement2003)

    pattern2004 = Pattern(Integral((F_**v_*WC('a', S(1)) + F_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1101, cons2, cons3, cons4, cons25, cons1156)
    rule2004 = ReplacementRule(pattern2004, replacement2004)

    pattern2005 = Pattern(Integral((F_**v_*WC('a', S(1)) + G_**w_*WC('b', S(1)))**n_*WC('u', S(1)), x_), cons1101, cons1139, cons2, cons3, cons4, cons25, cons1156)
    rule2005 = ReplacementRule(pattern2005, replacement2005)

    pattern2006 = Pattern(Integral(F_**v_*G_**w_*WC('u', S(1)), x_), cons1101, cons1139, cons1157)
    rule2006 = ReplacementRule(pattern2006, replacement2006)

    pattern2007 = Pattern(Integral(F_**u_*(v_ + w_)*WC('y', S(1)), x_), cons1101, cons1101, CustomConstraint(With2007))
    rule2007 = ReplacementRule(pattern2007, replacement2007)

    pattern2008 = Pattern(Integral(F_**u_*v_**WC('n', S(1))*w_, x_), cons1101, cons4, cons806, cons1019, cons1111, CustomConstraint(With2008))
    rule2008 = ReplacementRule(pattern2008, replacement2008)
    return [rule1904, rule1905, rule1906, rule1907, rule1908, rule1909, rule1910, rule1911, rule1912, rule1913, rule1914, rule1915, rule1916, rule1917, rule1918, rule1919, rule1920, rule1921, rule1922, rule1923, rule1924, rule1925, rule1926, rule1927, rule1928, rule1929, rule1930, rule1931, rule1932, rule1933, rule1934, rule1935, rule1936, rule1937, rule1938, rule1939, rule1940, rule1941, rule1942, rule1943, rule1944, rule1945, rule1946, rule1947, rule1948, rule1949, rule1950, rule1951, rule1952, rule1953, rule1954, rule1955, rule1956, rule1957, rule1958, rule1959, rule1960, rule1961, rule1962, rule1963, rule1964, rule1965, rule1966, rule1967, rule1968, rule1969, rule1970, rule1971, rule1972, rule1973, rule1974, rule1975, rule1976, rule1977, rule1978, rule1979, rule1980, rule1981, rule1982, rule1983, rule1984, rule1985, rule1986, rule1987, rule1988, rule1989, rule1990, rule1991, rule1992, rule1993, rule1994, rule1995, rule1996, rule1997, rule1998, rule1999, rule2000, rule2001, rule2002, rule2003, rule2004, rule2005, rule2006, rule2007, rule2008, ]





def replacement1904(F, b, c, d, e, f, g, m, n, x):
    return -Dist(d*m/(f*g*n*log(F)), Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(-1)), x), x) + Simp((F**(g*(e + f*x))*b)**n*(c + d*x)**m/(f*g*n*log(F)), x)


def replacement1905(F, b, c, d, e, f, g, m, n, x):
    return -Dist(f*g*n*log(F)/(d*(m + S(1))), Int((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1)), x), x) + Simp((F**(g*(e + f*x))*b)**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement1906(F, c, d, e, f, g, x):
    return Simp(F**(g*(-c*f/d + e))*ExpIntegralEi(f*g*(c + d*x)*log(F)/d)/d, x)


def replacement1907(F, c, d, e, f, g, m, x):
    return Simp(F**(g*(-c*f/d + e))*f**(-m + S(-1))*g**(-m + S(-1))*(-d)**m*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)*log(F)**(-m + S(-1)), x)


def replacement1908(F, c, d, e, f, g, x):
    return Dist(S(2)/d, Subst(Int(F**(g*(-c*f/d + e) + f*g*x**S(2)/d), x), x, sqrt(c + d*x)), x)


def replacement1909(F, c, d, e, f, g, m, x):
    return -Simp(F**(g*(-c*f/d + e))*(-f*g*log(F)/d)**(-IntPart(m) + S(-1))*(-f*g*(c + d*x)*log(F)/d)**(-FracPart(m))*(c + d*x)**FracPart(m)*Gamma(m + S(1), -f*g*(c + d*x)*log(F)/d)/d, x)


def replacement1910(F, b, c, d, e, f, g, m, n, x):
    return Dist(F**(-g*n*(e + f*x))*(F**(g*(e + f*x))*b)**n, Int(F**(g*n*(e + f*x))*(c + d*x)**m, x), x)


def replacement1911(F, a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a + b*(F**(g*(e + f*x)))**n)**p, x), x)


def replacement1912(F, a, b, c, d, e, f, g, m, n, x):
    return Dist(d*m/(a*f*g*n*log(F)), Int((c + d*x)**(m + S(-1))*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1)), x), x) - Simp((c + d*x)**m*log(a*(F**(g*(e + f*x)))**(-n)/b + S(1))/(a*f*g*n*log(F)), x)


def With1913(F, a, b, c, d, e, f, g, m, n, p, x):
    u = IntHide((a + b*(F**(g*(e + f*x)))**n)**p, x)
    return -Dist(d*m, Int(u*(c + d*x)**(m + S(-1)), x), x) + Dist((c + d*x)**m, u, x)


def replacement1914(F, a, b, g, m, n, p, u, v, x):
    return Int((a + b*(F**(g*ExpandToSum(v, x)))**n)**p*NormalizePowerOfLinear(u, x)**m, x)


def With1915(F, a, b, g, m, n, p, u, v, x):
    uu = NormalizePowerOfLinear(u, x)
    z = Symbol('z')
    z = If(And(PowerQ(uu), FreeQ(Part(uu, S(2)), x)), Part(uu, S(1))**(m*Part(uu, S(2))), uu**m)
    z = If(And(PowerQ(uu), FreeQ(Part(uu, 2), x)), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
    return Simp(uu**m*Int(z*(a + b*(F**(g*ExpandToSum(v, x)))**n)**p, x)/z, x)


def replacement1916(F, a, b, c, d, e, f, g, m, n, p, x):
    return Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m, x)


def replacement1917(F, a, b, c, d, e, f, g, m, n, x):
    return -Dist(d*m/(b*f*g*n*log(F)), Int((c + d*x)**(m + S(-1))*log(S(1) + b*(F**(g*(e + f*x)))**n/a), x), x) + Simp((c + d*x)**m*log(S(1) + b*(F**(g*(e + f*x)))**n/a)/(b*f*g*n*log(F)), x)


def replacement1918(F, a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(d*m/(b*f*g*n*(p + S(1))*log(F)), Int((a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**(m + S(-1)), x), x) + Simp((a + b*(F**(g*(e + f*x)))**n)**(p + S(1))*(c + d*x)**m/(b*f*g*n*(p + S(1))*log(F)), x)


def replacement1919(F, a, b, c, d, e, f, g, m, n, p, x):
    return Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x)


def replacement1920(F, G, a, b, c, d, e, f, g, h, i, j, k, m, n, p, q, x):
    return Dist((G**(j*(h + i*x))*k)**q*(F**(g*(e + f*x)))**(-n), Int((a + b*(F**(g*(e + f*x)))**n)**p*(c + d*x)**m*(F**(g*(e + f*x)))**n, x), x)


def replacement1921(F, a, b, c, n, x):
    return Simp((F**(c*(a + b*x)))**n/(b*c*n*log(F)), x)


def replacement1922(F, c, u, v, x):
    return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*u, x), x)


def replacement1923(F, c, u, v, x):
    return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), u, x), x)


def replacement1924(F, c, m, u, v, w, x):
    return Simp(F**(c*v)*u**(m + S(1))*Coefficient(w, x, S(1))/(c*Coefficient(u, x, S(1))*Coefficient(v, x, S(1))*log(F)), x)


def replacement1925(F, c, m, u, v, w, x):
    return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*NormalizePowerOfLinear(u, x)**m, x), x)


def replacement1926(F, c, m, u, v, w, x):
    return Int(ExpandIntegrand(F**(c*ExpandToSum(v, x)), w*NormalizePowerOfLinear(u, x)**m, x), x)


def With1927(F, c, m, u, v, w, x):
    uu = NormalizePowerOfLinear(u, x)
    z = Symbol('z')
    z = If(And(PowerQ(uu), FreeQ(Part(uu, S(2)), x)), Part(uu, S(1))**(m*Part(uu, S(2))), uu**m)
    z = If(And(PowerQ(uu), FreeQ(Part(uu, 2), x)), Part(uu, 1)**(m*Part(uu, 2)), uu**m)
    return Simp(uu**m*Int(ExpandIntegrand(F**(c*ExpandToSum(v, x))*w*z, x), x)/z, x)


def replacement1928(F, a, b, c, d, e, f, g, h, n, x):
    return Simp(F**(c*(a + b*x))*e*x*log(d*x)**(n + S(1))/(n + S(1)), x)


def replacement1929(F, a, b, c, d, e, f, g, h, m, n, x):
    return Simp(F**(c*(a + b*x))*e*x**(m + S(1))*log(d*x)**(n + S(1))/(n + S(1)), x)


def replacement1930(F, a, b, c, d, x):
    return Simp(F**(a + b*(c + d*x))/(b*d*log(F)), x)


def replacement1931(F, a, b, c, d, x):
    return Simp(F**a*sqrt(Pi)*Erfi((c + d*x)*Rt(b*log(F), S(2)))/(S(2)*d*Rt(b*log(F), S(2))), x)


def replacement1932(F, a, b, c, d, x):
    return Simp(F**a*sqrt(Pi)*Erf((c + d*x)*Rt(-b*log(F), S(2)))/(S(2)*d*Rt(-b*log(F), S(2))), x)


def replacement1933(F, a, b, c, d, n, x):
    return -Dist(b*n*log(F), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**n, x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)/d, x)


def With1934(F, a, b, c, d, n, x):
    k = Denominator(n)
    return Dist(k/d, Subst(Int(F**(a + b*x**(k*n))*x**(k + S(-1)), x), x, (c + d*x)**(S(1)/k)), x)


def replacement1935(F, a, b, c, d, n, x):
    return -Simp(F**a*(-b*(c + d*x)**n*log(F))**(-S(1)/n)*(c + d*x)*Gamma(S(1)/n, -b*(c + d*x)**n*log(F))/(d*n), x)


def replacement1936(F, a, b, c, d, e, f, m, n, x):
    return Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(-n)*(e + f*x)**n/(b*f*n*log(F)), x)


def replacement1937(F, a, b, c, d, e, f, n, x):
    return Simp(F**a*ExpIntegralEi(b*(c + d*x)**n*log(F))/(f*n), x)


def replacement1938(F, a, b, c, d, m, n, x):
    return Dist(S(1)/(d*(m + S(1))), Subst(Int(F**(a + b*x**S(2)), x), x, (c + d*x)**(m + S(1))), x)


def replacement1939(F, a, b, c, d, m, n, x):
    return -Dist((m - n + S(1))/(b*n*log(F)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)), x)


def replacement1940(F, a, b, c, d, m, n, x):
    return -Dist((m - n + S(1))/(b*n*log(F)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m - n + S(1))/(b*d*n*log(F)), x)


def replacement1941(F, a, b, c, d, m, n, x):
    return -Dist(b*n*log(F)/(m + S(1)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def replacement1942(F, a, b, c, d, m, n, x):
    return -Dist(b*n*log(F)/(m + S(1)), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + n), x), x) + Simp(F**(a + b*(c + d*x)**n)*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)


def With1943(F, a, b, c, d, m, n, x):
    k = Denominator(n)
    return Dist(k/d, Subst(Int(F**(a + b*x**(k*n))*x**(k*(m + S(1)) + S(-1)), x), x, (c + d*x)**(S(1)/k)), x)


def replacement1944(F, a, b, c, d, e, f, m, n, x):
    return Dist((c + d*x)**(-m)*(e + f*x)**m, Int(F**(a + b*(c + d*x)**n)*(c + d*x)**m, x), x)


def replacement1945(F, a, b, c, d, e, f, m, n, x):
    return -Simp(F**a*(-b*(c + d*x)**n*log(F))**(-(m + S(1))/n)*(e + f*x)**(m + S(1))*Gamma((m + S(1))/n, -b*(c + d*x)**n*log(F))/(f*n), x)


def replacement1946(F, a, b, c, d, e, f, m, x):
    return Dist((-c*f + d*e)/d, Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-1)), x), x) - Dist(f**S(2)*(m + S(-1))/(S(2)*b*d**S(2)*log(F)), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(-2)), x), x) + Simp(F**(a + b*(c + d*x)**S(2))*f*(e + f*x)**(m + S(-1))/(S(2)*b*d**S(2)*log(F)), x)


def replacement1947(F, a, b, c, d, e, f, m, x):
    return -Dist(S(2)*b*d**S(2)*log(F)/(f**S(2)*(m + S(1))), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(2)), x), x) + Dist(S(2)*b*d*(-c*f + d*e)*log(F)/(f**S(2)*(m + S(1))), Int(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1)), x), x) + Simp(F**(a + b*(c + d*x)**S(2))*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)


def replacement1948(F, a, b, c, d, e, f, m, n, x):
    return -Dist(b*d*n*log(F)/(f*(m + S(1))), Int(F**(a + b*(c + d*x)**n)*(c + d*x)**(n + S(-1))*(e + f*x)**(m + S(1)), x), x) + Simp(F**(a + b*(c + d*x)**n)*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)


def replacement1949(F, a, b, c, d, e, f, x):
    return Dist(d/f, Int(F**(a + b/(c + d*x))/(c + d*x), x), x) - Dist((-c*f + d*e)/f, Int(F**(a + b/(c + d*x))/((c + d*x)*(e + f*x)), x), x)


def replacement1950(F, a, b, c, d, e, f, m, x):
    return Dist(b*d*log(F)/(f*(m + S(1))), Int(F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(c + d*x)**S(2), x), x) + Simp(F**(a + b/(c + d*x))*(e + f*x)**(m + S(1))/(f*(m + S(1))), x)


def replacement1951(F, a, b, c, d, e, f, n, x):
    return Int(F**(a + b*(c + d*x)**n)/(e + f*x), x)


def replacement1952(F, m, u, v, x):
    return Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x)


def replacement1953(F, a, b, c, d, n, u, x):
    return Int(ExpandLinearProduct(F**(a + b*(c + d*x)**n), u, c, d, x), x)


def replacement1954(F, a, b, u, v, x):
    return Int(F**(a + b*NormalizePowerOfLinear(v, x))*u, x)


def replacement1955(F, a, b, c, d, e, f, g, h, x):
    return -Dist(d/(f*(-c*h + d*g)), Subst(Int(F**(a + b*d*x/(-c*h + d*g) - b*h/(-c*h + d*g))/x, x), x, (g + h*x)/(c + d*x)), x)


def replacement1956(F, a, b, c, d, e, f, g, h, m, x):
    return Dist(F**(b*f/d + e), Int((g + h*x)**m, x), x)


def replacement1957(F, a, b, c, d, e, f, g, h, m, x):
    return Int(F**(-f*(-a*d + b*c)/(d*(c + d*x)) + (b*f + d*e)/d)*(g + h*x)**m, x)


def replacement1958(F, a, b, c, d, e, f, g, h, x):
    return Dist(d/h, Int(F**(e + f*(a + b*x)/(c + d*x))/(c + d*x), x), x) - Dist((-c*h + d*g)/h, Int(F**(e + f*(a + b*x)/(c + d*x))/((c + d*x)*(g + h*x)), x), x)


def replacement1959(F, a, b, c, d, e, f, g, h, m, x):
    return -Dist(f*(-a*d + b*c)*log(F)/(h*(m + S(1))), Int(F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(c + d*x)**S(2), x), x) + Simp(F**(e + f*(a + b*x)/(c + d*x))*(g + h*x)**(m + S(1))/(h*(m + S(1))), x)


def replacement1960(F, a, b, c, d, e, f, g, h, i, j, x):
    return -Dist(d/(h*(-c*j + d*i)), Subst(Int(F**(e - f*x*(-a*d + b*c)/(-c*j + d*i) + f*(-a*j + b*i)/(-c*j + d*i))/x, x), x, (i + j*x)/(c + d*x)), x)


def replacement1961(F, a, b, c, x):
    return Dist(F**(a - b**S(2)/(S(4)*c)), Int(F**((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x)


def replacement1962(F, v, x):
    return Int(F**ExpandToSum(v, x), x)


def replacement1963(F, a, b, c, d, e, x):
    return Simp(F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)), x)


def replacement1964(F, a, b, c, d, e, m, x):
    return -Dist(e**S(2)*(m + S(-1))/(S(2)*c*log(F)), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)), x)


def replacement1965(F, a, b, c, d, e, x):
    return Simp(F**(a - b**S(2)/(S(4)*c))*ExpIntegralEi((b + S(2)*c*x)**S(2)*log(F)/(S(4)*c))/(S(2)*e), x)


def replacement1966(F, a, b, c, d, e, m, x):
    return -Dist(S(2)*c*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement1967(F, a, b, c, d, e, x):
    return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(F**(a + b*x + c*x**S(2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e/(S(2)*c*log(F)), x)


def replacement1968(F, a, b, c, d, e, m, x):
    return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-1)), x), x) - Dist(e**S(2)*(m + S(-1))/(S(2)*c*log(F)), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(-2)), x), x) + Simp(F**(a + b*x + c*x**S(2))*e*(d + e*x)**(m + S(-1))/(S(2)*c*log(F)), x)


def replacement1969(F, a, b, c, d, e, m, x):
    return -Dist(S(2)*c*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(2)), x), x) - Dist((b*e - S(2)*c*d)*log(F)/(e**S(2)*(m + S(1))), Int(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1)), x), x) + Simp(F**(a + b*x + c*x**S(2))*(d + e*x)**(m + S(1))/(e*(m + S(1))), x)


def replacement1970(F, a, b, c, d, e, m, x):
    return Int(F**(a + b*x + c*x**S(2))*(d + e*x)**m, x)


def replacement1971(F, m, u, v, x):
    return Int(F**ExpandToSum(v, x)*ExpandToSum(u, x)**m, x)


def With1972(F, a, b, c, d, e, m, n, v, x):
    u = IntHide(F**(e*(c + d*x))*(F**v*b + a)**n, x)
    return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Dist(x**m, u, x)


def With1973(F, G, a, b, c, d, e, f, g, h, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    m = FullSimplify(g*h*log(G)/(d*e*log(F)))
    if And(RationalQ(m), GreaterEqual(Abs(m), S(1))):
        return True
    return False


def replacement1973(F, G, a, b, c, d, e, f, g, h, n, x):

    m = FullSimplify(g*h*log(G)/(d*e*log(F)))
    return Dist(G**(-c*g*h/d + f*h)*Denominator(m)/(d*e*log(F)), Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m))), x)


def With1974(F, G, a, b, c, d, e, f, g, h, n, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    m = FullSimplify(d*e*log(F)/(g*h*log(G)))
    if And(RationalQ(m), Greater(Abs(m), S(1))):
        return True
    return False


def replacement1974(F, G, a, b, c, d, e, f, g, h, n, x):

    m = FullSimplify(d*e*log(F)/(g*h*log(G)))
    return Dist(Denominator(m)/(g*h*log(G)), Subst(Int(x**(Denominator(m) + S(-1))*(F**(c*e - d*e*f/g)*b*x**Numerator(m) + a)**n, x), x, G**(h*(f + g*x)/Denominator(m))), x)


def replacement1975(F, G, a, b, c, d, e, f, g, h, n, x):
    return Int(G**(f*h)*G**(g*h*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x)


def replacement1976(F, G, a, b, c, d, e, f, g, h, n, x):
    return Simp(G**(h*(f + g*x))*a**n*Hypergeometric2F1(-n, g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G)), x)


def replacement1977(F, G, a, b, c, d, e, f, g, h, n, x):
    return Simp(G**(h*(f + g*x))*(F**(e*(c + d*x))*b + a)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1) + g*h*log(G)/(d*e*log(F)), S(1) + g*h*log(G)/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(a*g*h*log(G)), x)


def replacement1978(F, G, a, b, e, h, n, u, v, x):
    return Int(G**(h*ExpandToSum(u, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x)


def With1979(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    m = FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))
    if RationalQ(m):
        return True
    return False


def replacement1979(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):

    m = FullSimplify((g*h*log(G) + s*t*log(H))/(d*e*log(F)))
    return Dist(G**(-c*g*h/d + f*h)*H**(-c*s*t/d + r*t)*Denominator(m)/(d*e*log(F)), Subst(Int(x**(Numerator(m) + S(-1))*(a + b*x**Denominator(m))**n, x), x, F**(e*(c + d*x)/Denominator(m))), x)


def replacement1980(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):
    return Dist(G**(h*(-c*g/d + f)), Int(H**(t*(r + s*x))*(b + F**(-e*(c + d*x))*a)**n, x), x)


def replacement1981(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):
    return Int(G**(f*h)*G**(g*h*x)*H**(r*t)*H**(s*t*x)*(F**(c*e)*F**(d*e*x)*b + a)**n, x)


def replacement1982(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):
    return Simp(G**(h*(f + g*x))*H**(t*(r + s*x))*a**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)), x)


def replacement1983(F, G, H, a, b, c, d, e, f, g, h, n, r, s, t, x):
    return Simp(G**(h*(f + g*x))*H**(t*(r + s*x))*((F**(e*(c + d*x))*b + a)/a)**(-n)*(F**(e*(c + d*x))*b + a)**n*Hypergeometric2F1(-n, (g*h*log(G) + s*t*log(H))/(d*e*log(F)), S(1) + (g*h*log(G) + s*t*log(H))/(d*e*log(F)), -F**(e*(c + d*x))*b/a)/(g*h*log(G) + s*t*log(H)), x)


def replacement1984(F, G, H, a, b, e, h, n, t, u, v, w, x):
    return Int(G**(h*ExpandToSum(u, x))*H**(t*ExpandToSum(w, x))*(F**(e*ExpandToSum(v, x))*b + a)**n, x)


def replacement1985(F, a, b, c, d, e, n, p, x):
    return -Dist(a*n/(b*d*e*log(F)), Int(x**(n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x), x) + Simp((F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)), x)


def replacement1986(F, a, b, c, d, e, m, n, p, x):
    return -Dist(a*n/(b*d*e*log(F)), Int(x**(m + n + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**p, x), x) - Dist(m/(b*d*e*(p + S(1))*log(F)), Int(x**(m + S(-1))*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1)), x), x) + Simp(x**m*(F**(e*(c + d*x))*b + a*x**n)**(p + S(1))/(b*d*e*(p + S(1))*log(F)), x)


def With1987(F, a, b, c, f, g, m, u, v, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int((f + g*x)**m/(S(2)*F**u*c + b - q), x), x) - Dist(S(2)*c/q, Int((f + g*x)**m/(S(2)*F**u*c + b + q), x), x)


def With1988(F, a, b, c, f, g, m, u, v, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b - q), x), x) - Dist(S(2)*c/q, Int(F**u*(f + g*x)**m/(S(2)*F**u*c + b + q), x), x)


def With1989(F, a, b, c, f, g, h, i, m, u, v, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return -Dist(-i + (-b*i + S(2)*c*h)/q, Int((f + g*x)**m/(S(2)*F**u*c + b + q), x), x) + Dist(i + (-b*i + S(2)*c*h)/q, Int((f + g*x)**m/(S(2)*F**u*c + b - q), x), x)


def With1990(F, a, b, c, d, m, v, x):
    u = IntHide(S(1)/(F**v*b + F**(c + d*x)*a), x)
    return -Dist(m, Int(u*x**(m + S(-1)), x), x) + Simp(u*x**m, x)


def replacement1991(F, a, b, c, u, v, w, x):
    return Int(F**v*u/(F**(S(2)*v)*b + F**v*a + c), x)


def replacement1992(F, a, b, c, d, e, g, n, x):
    return Int(ExpandIntegrand(F**(g*(d + e*x)**n), S(1)/(a + b*x + c*x**S(2)), x), x)


def replacement1993(F, a, c, d, e, g, n, x):
    return Int(ExpandIntegrand(F**(g*(d + e*x)**n), S(1)/(a + c*x**S(2)), x), x)


def replacement1994(F, a, b, c, d, e, g, m, n, u, x):
    return Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + b*x + c*x**S(2)), x), x)


def replacement1995(F, a, c, d, e, g, m, n, u, x):
    return Int(ExpandIntegrand(F**(g*(d + e*x)**n), u**m/(a + c*x**S(2)), x), x)


def replacement1996(F, a, b, x):
    return -Simp(sqrt(Pi)*Erf((-x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*exp(-S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))), x) + Simp(sqrt(Pi)*Erf((x**S(2)*sqrt(-b*log(F)) + sqrt(-a*log(F)))/x)*exp(S(2)*sqrt(-a*log(F))*sqrt(-b*log(F)))/(S(4)*sqrt(-b*log(F))), x)


def replacement1997(m, n, x):
    return Dist(m, Int(x**(m + S(-1))*(x**m + exp(x))**n, x), x) + Int((x**m + exp(x))**(n + S(1)), x) - Simp((x**m + exp(x))**(n + S(1))/(n + S(1)), x)


def replacement1998(F, a, b, c, d, e, n, x):
    return Dist(S(1)/(d*e*n*log(F)), Subst(Int(log(a + b*x)/x, x), x, (F**(e*(c + d*x)))**n), x)


def replacement1999(F, a, b, c, d, e, n, x):
    return -Dist(b*d*e*n*log(F), Int(x*(F**(e*(c + d*x)))**n/(a + b*(F**(e*(c + d*x)))**n), x), x) + Simp(x*log(a + b*(F**(e*(c + d*x)))**n), x)


def replacement2000(F, a, n, u, v, x):
    return Dist(F**(-n*v)*(F**v*a)**n, Int(F**(n*v)*u, x), x)


def With2001(u, x):
    v = FunctionOfExponential(u, x)
    return Dist(v/D(v, x), Subst(Int(FunctionOfExponentialFunction(u, x)/x, x), x, v), x)


def replacement2002(F, a, b, n, u, v, w, x):
    return Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x)


def replacement2003(F, G, a, b, n, u, v, w, x):
    return Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x)


def replacement2004(F, a, b, n, u, v, w, x):
    return Dist(F**(-n*v)*(F**v*a + F**w*b)**n*(F**ExpandToSum(-v + w, x)*b + a)**(-n), Int(F**(n*v)*u*(F**ExpandToSum(-v + w, x)*b + a)**n, x), x)


def replacement2005(F, G, a, b, n, u, v, w, x):
    return Dist(F**(-n*v)*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**(-n)*(F**v*a + G**w*b)**n, Int(F**(n*v)*u*(a + b*exp(ExpandToSum(-v*log(F) + w*log(G), x)))**n, x), x)


def replacement2006(F, G, u, v, w, x):
    return Int(u*NormalizeIntegrand(exp(v*log(F) + w*log(G)), x), x)


def With2007(F, u, v, w, x, y):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    z = v*y/(D(u, x)*log(F))
    if ZeroQ(-w*y + D(z, x)):
        return True
    return False


def replacement2007(F, u, v, w, x, y):

    z = v*y/(D(u, x)*log(F))
    return Simp(F**u*z, x)


def With2008(F, n, u, v, w, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    z = v*D(u, x)*log(F) + (n + S(1))*D(v, x)
    if And(Equal(Exponent(w, x), Exponent(z, x)), ZeroQ(w*Coefficient(z, x, Exponent(z, x)) - z*Coefficient(w, x, Exponent(w, x)))):
        return True
    return False


def replacement2008(F, n, u, v, w, x):

    z = v*D(u, x)*log(F) + (n + S(1))*D(v, x)
    return Simp(F**u*v**(n + S(1))*Coefficient(w, x, Exponent(w, x))/Coefficient(z, x, Exponent(z, x)), x)
