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


def secant():
    from sympy.integrals.rubi.constraints import cons1583, cons1504, cons2, cons3, cons50, cons127, cons19, cons4, cons1584, cons1514, cons1585, cons21, cons95, cons168, cons91, cons1172, cons96, cons167, cons33, cons1586, cons89, cons1361, cons25, cons1257, cons1260, cons676, cons8, cons29, cons810, cons1263, cons1587, cons1266, cons1267, cons1588, cons545, cons45, cons450, cons1269, cons746, cons1589, cons1256, cons64, cons1425, cons517, cons1322, cons1323, cons1590, cons1591, cons1592, cons1332, cons113, cons157, cons1593, cons1521, cons1338, cons1594, cons465, cons87, cons1595, cons79, cons170, cons274, cons1335, cons1596, cons1336, cons1597, cons1327, cons1598, cons1599, cons1600, cons1601, cons1555, cons1602, cons1359, cons1603, cons1604, cons1605, cons1606, cons20, cons210, cons5, cons1276, cons1607, cons1608, cons1310, cons149, cons1230, cons1509, cons150, cons1517, cons1609, cons198, cons1313, cons1610, cons1611, cons1582, cons72, cons1612, cons1613, cons81, cons1614, cons1615, cons1306, cons73, cons1414, cons1411, cons1325, cons1324, cons1616, cons82, cons1362, cons1423, cons1317, cons1233, cons1617, cons1316, cons1268, cons1618, cons152, cons1619, cons1620, cons1621, cons1622, cons1623, cons40, cons1624, cons1417, cons382, cons1430, cons36, cons37, cons1247, cons1571, cons1625, cons1626, cons1627, cons1628, cons1629, cons34, cons1551, cons1630, cons1631, cons348, cons90, cons1329, cons1632, cons1633, cons1258, cons1634, cons377, cons35, cons38, cons1435, cons1635, cons1636, cons1433, cons1637, cons1638, cons1639, cons1640, cons1641, cons1642, cons685, cons1643, cons1644, cons1645, cons1456, cons1480, cons56, cons1482, cons1481, cons1483, cons378, cons48, cons47, cons228, cons1646, cons530, cons812, cons813, cons1575, cons1497, cons70, cons71, cons825, cons826, cons1576, cons1578, cons1499, cons1579, cons1647


    pattern3920 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1583, cons1504)
    rule3920 = ReplacementRule(pattern3920, replacement3920)

    pattern3921 = Pattern(Integral((S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons50, cons127, cons1584)
    rule3921 = ReplacementRule(pattern3921, replacement3921)

    pattern3922 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1514, cons1585, cons21)
    rule3922 = ReplacementRule(pattern3922, replacement3922)

    pattern3923 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1514, cons1585, cons21)
    rule3923 = ReplacementRule(pattern3923, replacement3923)

    pattern3924 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons168, cons91, cons1172)
    rule3924 = ReplacementRule(pattern3924, replacement3924)

    pattern3925 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons96, cons167, cons1172)
    rule3925 = ReplacementRule(pattern3925, replacement3925)

    pattern3926 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons168, cons1172, cons1586)
    rule3926 = ReplacementRule(pattern3926, replacement3926)

    pattern3927 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons167, cons1172)
    rule3927 = ReplacementRule(pattern3927, replacement3927)

    pattern3928 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons96, cons1361, cons1172)
    rule3928 = ReplacementRule(pattern3928, replacement3928)

    pattern3929 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons91, cons1361, cons1172)
    rule3929 = ReplacementRule(pattern3929, replacement3929)

    pattern3930 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons25, cons1257)
    rule3930 = ReplacementRule(pattern3930, replacement3930)

    pattern3931 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1260)
    rule3931 = ReplacementRule(pattern3931, replacement3931)

    pattern3932 = Pattern(Integral((S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons8, cons29, cons676)
    rule3932 = ReplacementRule(pattern3932, replacement3932)

    pattern3933 = Pattern(Integral((S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons8, cons29, cons676)
    rule3933 = ReplacementRule(pattern3933, replacement3933)

    pattern3934 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons167, cons810)
    rule3934 = ReplacementRule(pattern3934, replacement3934)

    pattern3935 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons167, cons810)
    rule3935 = ReplacementRule(pattern3935, replacement3935)

    pattern3936 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons91, cons810)
    rule3936 = ReplacementRule(pattern3936, replacement3936)

    pattern3937 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons91, cons810)
    rule3937 = ReplacementRule(pattern3937, replacement3937)

    pattern3938 = Pattern(Integral(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons1263)
    rule3938 = ReplacementRule(pattern3938, replacement3938)

    pattern3939 = Pattern(Integral(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons1263)
    rule3939 = ReplacementRule(pattern3939, replacement3939)

    pattern3940 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons1587)
    rule3940 = ReplacementRule(pattern3940, replacement3940)

    pattern3941 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons1587)
    rule3941 = ReplacementRule(pattern3941, replacement3941)

    pattern3942 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons25)
    rule3942 = ReplacementRule(pattern3942, replacement3942)

    pattern3943 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons25)
    rule3943 = ReplacementRule(pattern3943, replacement3943)

    pattern3944 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule3944 = ReplacementRule(pattern3944, replacement3944)

    pattern3945 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule3945 = ReplacementRule(pattern3945, replacement3945)

    pattern3946 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule3946 = ReplacementRule(pattern3946, replacement3946)

    pattern3947 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule3947 = ReplacementRule(pattern3947, replacement3947)

    pattern3948 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons167, cons810)
    rule3948 = ReplacementRule(pattern3948, replacement3948)

    pattern3949 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons167, cons810)
    rule3949 = ReplacementRule(pattern3949, replacement3949)

    pattern3950 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule3950 = ReplacementRule(pattern3950, replacement3950)

    pattern3951 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule3951 = ReplacementRule(pattern3951, replacement3951)

    pattern3952 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons1588, cons810)
    rule3952 = ReplacementRule(pattern3952, replacement3952)

    pattern3953 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons1588, cons810)
    rule3953 = ReplacementRule(pattern3953, replacement3953)

    pattern3954 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons45)
    rule3954 = ReplacementRule(pattern3954, replacement3954)

    pattern3955 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons45)
    rule3955 = ReplacementRule(pattern3955, replacement3955)

    pattern3956 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons450)
    rule3956 = ReplacementRule(pattern3956, replacement3956)

    pattern3957 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons450)
    rule3957 = ReplacementRule(pattern3957, replacement3957)

    pattern3958 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3958 = ReplacementRule(pattern3958, replacement3958)

    pattern3959 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3959 = ReplacementRule(pattern3959, replacement3959)

    pattern3960 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3960 = ReplacementRule(pattern3960, replacement3960)

    pattern3961 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3961 = ReplacementRule(pattern3961, replacement3961)

    pattern3962 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons746, cons810)
    rule3962 = ReplacementRule(pattern3962, replacement3962)

    pattern3963 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons746, cons810)
    rule3963 = ReplacementRule(pattern3963, replacement3963)

    pattern3964 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3964 = ReplacementRule(pattern3964, replacement3964)

    pattern3965 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3965 = ReplacementRule(pattern3965, replacement3965)

    pattern3966 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3966 = ReplacementRule(pattern3966, replacement3966)

    pattern3967 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule3967 = ReplacementRule(pattern3967, replacement3967)

    pattern3968 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons91, cons810)
    rule3968 = ReplacementRule(pattern3968, replacement3968)

    pattern3969 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons91, cons810)
    rule3969 = ReplacementRule(pattern3969, replacement3969)

    pattern3970 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1269, cons545)
    rule3970 = ReplacementRule(pattern3970, replacement3970)

    pattern3971 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1269, cons545)
    rule3971 = ReplacementRule(pattern3971, replacement3971)

    pattern3972 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1589)
    rule3972 = ReplacementRule(pattern3972, replacement3972)

    pattern3973 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1589)
    rule3973 = ReplacementRule(pattern3973, replacement3973)

    pattern3974 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1589)
    rule3974 = ReplacementRule(pattern3974, replacement3974)

    pattern3975 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1589)
    rule3975 = ReplacementRule(pattern3975, replacement3975)

    pattern3976 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons1256)
    rule3976 = ReplacementRule(pattern3976, replacement3976)

    pattern3977 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons1256)
    rule3977 = ReplacementRule(pattern3977, replacement3977)

    pattern3978 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(3)), x_), cons2, cons3, cons50, cons127, cons1256)
    rule3978 = ReplacementRule(pattern3978, replacement3978)

    pattern3979 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(3)), x_), cons2, cons3, cons50, cons127, cons1256)
    rule3979 = ReplacementRule(pattern3979, replacement3979)

    pattern3980 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons64, cons89)
    rule3980 = ReplacementRule(pattern3980, replacement3980)

    pattern3981 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons64, cons89)
    rule3981 = ReplacementRule(pattern3981, replacement3981)

    pattern3982 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3982 = ReplacementRule(pattern3982, replacement3982)

    pattern3983 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3983 = ReplacementRule(pattern3983, replacement3983)

    pattern3984 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1425, cons517)
    rule3984 = ReplacementRule(pattern3984, replacement3984)

    pattern3985 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1425, cons517)
    rule3985 = ReplacementRule(pattern3985, replacement3985)

    pattern3986 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3986 = ReplacementRule(pattern3986, replacement3986)

    pattern3987 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3987 = ReplacementRule(pattern3987, replacement3987)

    pattern3988 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3988 = ReplacementRule(pattern3988, replacement3988)

    pattern3989 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1267)
    rule3989 = ReplacementRule(pattern3989, replacement3989)

    pattern3990 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322, cons517)
    rule3990 = ReplacementRule(pattern3990, replacement3990)

    pattern3991 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322, cons517)
    rule3991 = ReplacementRule(pattern3991, replacement3991)

    pattern3992 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule3992 = ReplacementRule(pattern3992, replacement3992)

    pattern3993 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule3993 = ReplacementRule(pattern3993, replacement3993)

    pattern3994 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule3994 = ReplacementRule(pattern3994, replacement3994)

    pattern3995 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule3995 = ReplacementRule(pattern3995, replacement3995)

    pattern3996 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule3996 = ReplacementRule(pattern3996, replacement3996)

    pattern3997 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule3997 = ReplacementRule(pattern3997, replacement3997)

    pattern3998 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule3998 = ReplacementRule(pattern3998, replacement3998)

    pattern3999 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule3999 = ReplacementRule(pattern3999, replacement3999)

    pattern4000 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1590)
    rule4000 = ReplacementRule(pattern4000, replacement4000)

    pattern4001 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1590)
    rule4001 = ReplacementRule(pattern4001, replacement4001)

    pattern4002 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1591)
    rule4002 = ReplacementRule(pattern4002, replacement4002)

    pattern4003 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1591)
    rule4003 = ReplacementRule(pattern4003, replacement4003)

    pattern4004 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons167, cons810)
    rule4004 = ReplacementRule(pattern4004, replacement4004)

    pattern4005 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons167, cons810)
    rule4005 = ReplacementRule(pattern4005, replacement4005)

    pattern4006 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4006 = ReplacementRule(pattern4006, replacement4006)

    pattern4007 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4007 = ReplacementRule(pattern4007, replacement4007)

    pattern4008 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons1592, cons810)
    rule4008 = ReplacementRule(pattern4008, replacement4008)

    pattern4009 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons1592, cons810)
    rule4009 = ReplacementRule(pattern4009, replacement4009)

    pattern4010 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267)
    rule4010 = ReplacementRule(pattern4010, replacement4010)

    pattern4011 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267)
    rule4011 = ReplacementRule(pattern4011, replacement4011)

    pattern4012 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332, cons45)
    rule4012 = ReplacementRule(pattern4012, replacement4012)

    pattern4013 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332, cons45)
    rule4013 = ReplacementRule(pattern4013, replacement4013)

    pattern4014 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4014 = ReplacementRule(pattern4014, replacement4014)

    pattern4015 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4015 = ReplacementRule(pattern4015, replacement4015)

    pattern4016 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1257, cons33, cons1425, cons517)
    rule4016 = ReplacementRule(pattern4016, replacement4016)

    pattern4017 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1257, cons33, cons1425, cons517)
    rule4017 = ReplacementRule(pattern4017, replacement4017)

    pattern4018 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1257, cons33, cons1322, cons517)
    rule4018 = ReplacementRule(pattern4018, replacement4018)

    pattern4019 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1257, cons33, cons1322, cons517)
    rule4019 = ReplacementRule(pattern4019, replacement4019)

    pattern4020 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons113, cons1322)
    rule4020 = ReplacementRule(pattern4020, replacement4020)

    pattern4021 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons113, cons1322)
    rule4021 = ReplacementRule(pattern4021, replacement4021)

    pattern4022 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons157, cons1323)
    rule4022 = ReplacementRule(pattern4022, replacement4022)

    pattern4023 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons157, cons1323)
    rule4023 = ReplacementRule(pattern4023, replacement4023)

    pattern4024 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons168, cons1593, cons517)
    rule4024 = ReplacementRule(pattern4024, replacement4024)

    pattern4025 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons168, cons1593, cons517)
    rule4025 = ReplacementRule(pattern4025, replacement4025)

    pattern4026 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons168, cons1521, cons517)
    rule4026 = ReplacementRule(pattern4026, replacement4026)

    pattern4027 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons168, cons1521, cons517)
    rule4027 = ReplacementRule(pattern4027, replacement4027)

    pattern4028 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons96, cons1338, cons1594)
    rule4028 = ReplacementRule(pattern4028, replacement4028)

    pattern4029 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons96, cons1338, cons1594)
    rule4029 = ReplacementRule(pattern4029, replacement4029)

    pattern4030 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons96, cons746, cons1594)
    rule4030 = ReplacementRule(pattern4030, replacement4030)

    pattern4031 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons95, cons96, cons746, cons1594)
    rule4031 = ReplacementRule(pattern4031, replacement4031)

    pattern4032 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons96, cons1594)
    rule4032 = ReplacementRule(pattern4032, replacement4032)

    pattern4033 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons96, cons1594)
    rule4033 = ReplacementRule(pattern4033, replacement4033)

    pattern4034 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons167)
    rule4034 = ReplacementRule(pattern4034, replacement4034)

    pattern4035 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons167)
    rule4035 = ReplacementRule(pattern4035, replacement4035)

    pattern4036 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons465)
    rule4036 = ReplacementRule(pattern4036, replacement4036)

    pattern4037 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons465)
    rule4037 = ReplacementRule(pattern4037, replacement4037)

    pattern4038 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267)
    rule4038 = ReplacementRule(pattern4038, replacement4038)

    pattern4039 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267)
    rule4039 = ReplacementRule(pattern4039, replacement4039)

    pattern4040 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4040 = ReplacementRule(pattern4040, replacement4040)

    pattern4041 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267)
    rule4041 = ReplacementRule(pattern4041, replacement4041)

    pattern4042 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons746, cons810)
    rule4042 = ReplacementRule(pattern4042, replacement4042)

    pattern4043 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons746, cons810)
    rule4043 = ReplacementRule(pattern4043, replacement4043)

    pattern4044 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons465, cons810)
    rule4044 = ReplacementRule(pattern4044, replacement4044)

    pattern4045 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons89, cons465, cons810)
    rule4045 = ReplacementRule(pattern4045, replacement4045)

    pattern4046 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1267, cons89, cons746, cons1521, cons87)
    rule4046 = ReplacementRule(pattern4046, replacement4046)

    pattern4047 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1267, cons89, cons746, cons1521, cons87)
    rule4047 = ReplacementRule(pattern4047, replacement4047)

    pattern4048 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons25, cons1590)
    rule4048 = ReplacementRule(pattern4048, replacement4048)

    pattern4049 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons25, cons1590)
    rule4049 = ReplacementRule(pattern4049, replacement4049)

    pattern4050 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons25, cons1595)
    rule4050 = ReplacementRule(pattern4050, replacement4050)

    pattern4051 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons25, cons1595)
    rule4051 = ReplacementRule(pattern4051, replacement4051)

    pattern4052 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45)
    rule4052 = ReplacementRule(pattern4052, replacement4052)

    pattern4053 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45)
    rule4053 = ReplacementRule(pattern4053, replacement4053)

    pattern4054 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons450)
    rule4054 = ReplacementRule(pattern4054, replacement4054)

    pattern4055 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons450)
    rule4055 = ReplacementRule(pattern4055, replacement4055)

    pattern4056 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4056 = ReplacementRule(pattern4056, replacement4056)

    pattern4057 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4057 = ReplacementRule(pattern4057, replacement4057)

    pattern4058 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons168, cons517)
    rule4058 = ReplacementRule(pattern4058, replacement4058)

    pattern4059 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons168, cons517)
    rule4059 = ReplacementRule(pattern4059, replacement4059)

    pattern4060 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4060 = ReplacementRule(pattern4060, replacement4060)

    pattern4061 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4061 = ReplacementRule(pattern4061, replacement4061)

    pattern4062 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4062 = ReplacementRule(pattern4062, replacement4062)

    pattern4063 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4063 = ReplacementRule(pattern4063, replacement4063)

    pattern4064 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96, cons517)
    rule4064 = ReplacementRule(pattern4064, replacement4064)

    pattern4065 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96, cons517)
    rule4065 = ReplacementRule(pattern4065, replacement4065)

    pattern4066 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons79)
    rule4066 = ReplacementRule(pattern4066, replacement4066)

    pattern4067 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons79)
    rule4067 = ReplacementRule(pattern4067, replacement4067)

    pattern4068 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons170)
    rule4068 = ReplacementRule(pattern4068, replacement4068)

    pattern4069 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons170)
    rule4069 = ReplacementRule(pattern4069, replacement4069)

    pattern4070 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96)
    rule4070 = ReplacementRule(pattern4070, replacement4070)

    pattern4071 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96)
    rule4071 = ReplacementRule(pattern4071, replacement4071)

    pattern4072 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4072 = ReplacementRule(pattern4072, replacement4072)

    pattern4073 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4073 = ReplacementRule(pattern4073, replacement4073)

    pattern4074 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1269)
    rule4074 = ReplacementRule(pattern4074, replacement4074)

    pattern4075 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1269)
    rule4075 = ReplacementRule(pattern4075, replacement4075)

    pattern4076 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96)
    rule4076 = ReplacementRule(pattern4076, replacement4076)

    pattern4077 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96)
    rule4077 = ReplacementRule(pattern4077, replacement4077)

    pattern4078 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons274)
    rule4078 = ReplacementRule(pattern4078, replacement4078)

    pattern4079 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(3), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons274)
    rule4079 = ReplacementRule(pattern4079, replacement4079)

    pattern4080 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1335, cons1596)
    rule4080 = ReplacementRule(pattern4080, replacement4080)

    pattern4081 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1335, cons1596)
    rule4081 = ReplacementRule(pattern4081, replacement4081)

    pattern4082 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons33, cons1335, cons1336, cons1597)
    rule4082 = ReplacementRule(pattern4082, replacement4082)

    pattern4083 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons33, cons1335, cons1336, cons1597)
    rule4083 = ReplacementRule(pattern4083, replacement4083)

    pattern4084 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1327, cons1172)
    rule4084 = ReplacementRule(pattern4084, replacement4084)

    pattern4085 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1327, cons1172)
    rule4085 = ReplacementRule(pattern4085, replacement4085)

    pattern4086 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1338, cons1172)
    rule4086 = ReplacementRule(pattern4086, replacement4086)

    pattern4087 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1338, cons1172)
    rule4087 = ReplacementRule(pattern4087, replacement4087)

    pattern4088 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1598)
    rule4088 = ReplacementRule(pattern4088, replacement4088)

    pattern4089 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons96, cons1598)
    rule4089 = ReplacementRule(pattern4089, replacement4089)

    pattern4090 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1599)
    rule4090 = ReplacementRule(pattern4090, replacement4090)

    pattern4091 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1599)
    rule4091 = ReplacementRule(pattern4091, replacement4091)

    pattern4092 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons33, cons96, cons1172)
    rule4092 = ReplacementRule(pattern4092, replacement4092)

    pattern4093 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons33, cons96, cons1172)
    rule4093 = ReplacementRule(pattern4093, replacement4093)

    pattern4094 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4094 = ReplacementRule(pattern4094, replacement4094)

    pattern4095 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4095 = ReplacementRule(pattern4095, replacement4095)

    pattern4096 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4096 = ReplacementRule(pattern4096, replacement4096)

    pattern4097 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4097 = ReplacementRule(pattern4097, replacement4097)

    pattern4098 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4098 = ReplacementRule(pattern4098, replacement4098)

    pattern4099 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4099 = ReplacementRule(pattern4099, replacement4099)

    pattern4100 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1600)
    rule4100 = ReplacementRule(pattern4100, replacement4100)

    pattern4101 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1600)
    rule4101 = ReplacementRule(pattern4101, replacement4101)

    pattern4102 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4102 = ReplacementRule(pattern4102, replacement4102)

    pattern4103 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4103 = ReplacementRule(pattern4103, replacement4103)

    pattern4104 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons810)
    rule4104 = ReplacementRule(pattern4104, replacement4104)

    pattern4105 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons810)
    rule4105 = ReplacementRule(pattern4105, replacement4105)

    pattern4106 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4106 = ReplacementRule(pattern4106, replacement4106)

    pattern4107 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4107 = ReplacementRule(pattern4107, replacement4107)

    pattern4108 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons167, cons810)
    rule4108 = ReplacementRule(pattern4108, replacement4108)

    pattern4109 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons167, cons810)
    rule4109 = ReplacementRule(pattern4109, replacement4109)

    pattern4110 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4110 = ReplacementRule(pattern4110, replacement4110)

    pattern4111 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4111 = ReplacementRule(pattern4111, replacement4111)

    pattern4112 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons810)
    rule4112 = ReplacementRule(pattern4112, replacement4112)

    pattern4113 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons810)
    rule4113 = ReplacementRule(pattern4113, replacement4113)

    pattern4114 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4114 = ReplacementRule(pattern4114, replacement4114)

    pattern4115 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4115 = ReplacementRule(pattern4115, replacement4115)

    pattern4116 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4116 = ReplacementRule(pattern4116, replacement4116)

    pattern4117 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4117 = ReplacementRule(pattern4117, replacement4117)

    pattern4118 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons746, cons810)
    rule4118 = ReplacementRule(pattern4118, replacement4118)

    pattern4119 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons746, cons810)
    rule4119 = ReplacementRule(pattern4119, replacement4119)

    pattern4120 = Pattern(Integral(cos(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4120 = ReplacementRule(pattern4120, replacement4120)

    pattern4121 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1269)
    rule4121 = ReplacementRule(pattern4121, replacement4121)

    pattern4122 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4122 = ReplacementRule(pattern4122, replacement4122)

    pattern4123 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4123 = ReplacementRule(pattern4123, replacement4123)

    pattern4124 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons91, cons810)
    rule4124 = ReplacementRule(pattern4124, replacement4124)

    pattern4125 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons91, cons810)
    rule4125 = ReplacementRule(pattern4125, replacement4125)

    pattern4126 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons1601)
    rule4126 = ReplacementRule(pattern4126, replacement4126)

    pattern4127 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons89, cons1588, cons1601)
    rule4127 = ReplacementRule(pattern4127, replacement4127)

    pattern4128 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons89, cons1600, cons1555, cons1602)
    rule4128 = ReplacementRule(pattern4128, replacement4128)

    pattern4129 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons89, cons1600, cons1555, cons1602)
    rule4129 = ReplacementRule(pattern4129, replacement4129)

    pattern4130 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1359, cons1603, cons1521, cons1336)
    rule4130 = ReplacementRule(pattern4130, replacement4130)

    pattern4131 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1359, cons1603, cons1521, cons1336)
    rule4131 = ReplacementRule(pattern4131, replacement4131)

    pattern4132 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1604, cons1605, cons1521, cons1555)
    rule4132 = ReplacementRule(pattern4132, replacement4132)

    pattern4133 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons95, cons1604, cons1605, cons1521, cons1555)
    rule4133 = ReplacementRule(pattern4133, replacement4133)

    pattern4134 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4134 = ReplacementRule(pattern4134, replacement4134)

    pattern4135 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule4135 = ReplacementRule(pattern4135, replacement4135)

    pattern4136 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons64)
    rule4136 = ReplacementRule(pattern4136, replacement4136)

    pattern4137 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons64)
    rule4137 = ReplacementRule(pattern4137, replacement4137)

    pattern4138 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1606)
    rule4138 = ReplacementRule(pattern4138, replacement4138)

    pattern4139 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1606)
    rule4139 = ReplacementRule(pattern4139, replacement4139)

    pattern4140 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons20)
    rule4140 = ReplacementRule(pattern4140, replacement4140)

    pattern4141 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons20)
    rule4141 = ReplacementRule(pattern4141, replacement4141)

    pattern4142 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1267)
    rule4142 = ReplacementRule(pattern4142, replacement4142)

    pattern4143 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1267)
    rule4143 = ReplacementRule(pattern4143, replacement4143)

    pattern4144 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1269)
    rule4144 = ReplacementRule(pattern4144, replacement4144)

    pattern4145 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1269)
    rule4145 = ReplacementRule(pattern4145, replacement4145)

    pattern4146 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1607)
    rule4146 = ReplacementRule(pattern4146, replacement4146)

    pattern4147 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1607)
    rule4147 = ReplacementRule(pattern4147, replacement4147)

    pattern4148 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1608)
    rule4148 = ReplacementRule(pattern4148, replacement4148)

    pattern4149 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1608)
    rule4149 = ReplacementRule(pattern4149, replacement4149)

    pattern4150 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1310)
    rule4150 = ReplacementRule(pattern4150, replacement4150)

    pattern4151 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1310)
    rule4151 = ReplacementRule(pattern4151, replacement4151)

    pattern4152 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule4152 = ReplacementRule(pattern4152, replacement4152)

    pattern4153 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule4153 = ReplacementRule(pattern4153, replacement4153)

    pattern4154 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1230, cons1267, cons87)
    rule4154 = ReplacementRule(pattern4154, replacement4154)

    pattern4155 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1230, cons1267, cons87)
    rule4155 = ReplacementRule(pattern4155, replacement4155)

    pattern4156 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1230, cons1267, cons25)
    rule4156 = ReplacementRule(pattern4156, replacement4156)

    pattern4157 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1230, cons1267, cons25)
    rule4157 = ReplacementRule(pattern4157, replacement4157)

    pattern4158 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168)
    rule4158 = ReplacementRule(pattern4158, replacement4158)

    pattern4159 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168)
    rule4159 = ReplacementRule(pattern4159, replacement4159)

    pattern4160 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96)
    rule4160 = ReplacementRule(pattern4160, replacement4160)

    pattern4161 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96)
    rule4161 = ReplacementRule(pattern4161, replacement4161)

    pattern4162 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))/tan(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule4162 = ReplacementRule(pattern4162, replacement4162)

    pattern4163 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))*tan(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1266)
    rule4163 = ReplacementRule(pattern4163, replacement4163)

    pattern4164 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule4164 = ReplacementRule(pattern4164, replacement4164)

    pattern4165 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule4165 = ReplacementRule(pattern4165, replacement4165)

    pattern4166 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1230, cons1269)
    rule4166 = ReplacementRule(pattern4166, replacement4166)

    pattern4167 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1230, cons1269)
    rule4167 = ReplacementRule(pattern4167, replacement4167)

    pattern4168 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150)
    rule4168 = ReplacementRule(pattern4168, replacement4168)

    pattern4169 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons150)
    rule4169 = ReplacementRule(pattern4169, replacement4169)

    pattern4170 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1267, cons1517, cons1609)
    rule4170 = ReplacementRule(pattern4170, replacement4170)

    pattern4171 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1267, cons1517, cons1609)
    rule4171 = ReplacementRule(pattern4171, replacement4171)

    pattern4172 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1267, cons198)
    rule4172 = ReplacementRule(pattern4172, replacement4172)

    pattern4173 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1267, cons198)
    rule4173 = ReplacementRule(pattern4173, replacement4173)

    pattern4174 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1267, cons25)
    rule4174 = ReplacementRule(pattern4174, replacement4174)

    pattern4175 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1267, cons25)
    rule4175 = ReplacementRule(pattern4175, replacement4175)

    pattern4176 = Pattern(Integral(sqrt(WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))/(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269)
    rule4176 = ReplacementRule(pattern4176, replacement4176)

    pattern4177 = Pattern(Integral(sqrt(WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))/(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269)
    rule4177 = ReplacementRule(pattern4177, replacement4177)

    pattern4178 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_/(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269, cons1313)
    rule4178 = ReplacementRule(pattern4178, replacement4178)

    pattern4179 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_/(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269, cons1313)
    rule4179 = ReplacementRule(pattern4179, replacement4179)

    pattern4180 = Pattern(Integral(S(1)/(sqrt(WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons1269)
    rule4180 = ReplacementRule(pattern4180, replacement4180)

    pattern4181 = Pattern(Integral(S(1)/(sqrt(WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons1269)
    rule4181 = ReplacementRule(pattern4181, replacement4181)

    pattern4182 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_/(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269, cons1610)
    rule4182 = ReplacementRule(pattern4182, replacement4182)

    pattern4183 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_/(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1269, cons1610)
    rule4183 = ReplacementRule(pattern4183, replacement4183)

    pattern4184 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons1269)
    rule4184 = ReplacementRule(pattern4184, replacement4184)

    pattern4185 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_/tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons1269)
    rule4185 = ReplacementRule(pattern4185, replacement4185)

    pattern4186 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1269, cons150)
    rule4186 = ReplacementRule(pattern4186, replacement4186)

    pattern4187 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1269, cons150)
    rule4187 = ReplacementRule(pattern4187, replacement4187)

    pattern4188 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1269, cons87, cons20, cons1611)
    rule4188 = ReplacementRule(pattern4188, replacement4188)

    pattern4189 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1269, cons87, cons20, cons1611)
    rule4189 = ReplacementRule(pattern4189, replacement4189)

    pattern4190 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule4190 = ReplacementRule(pattern4190, replacement4190)

    pattern4191 = Pattern(Integral((WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule4191 = ReplacementRule(pattern4191, replacement4191)

    pattern4192 = Pattern(Integral((WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**p_)**m_*(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons21)
    rule4192 = ReplacementRule(pattern4192, replacement4192)

    pattern4193 = Pattern(Integral(((S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**p_*WC('e', S(1)))**m_*(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons21)
    rule4193 = ReplacementRule(pattern4193, replacement4193)

    pattern4194 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons64, cons198, cons1612)
    rule4194 = ReplacementRule(pattern4194, replacement4194)

    pattern4195 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons64, cons198, cons1612)
    rule4195 = ReplacementRule(pattern4195, replacement4195)

    pattern4196 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons20, cons89, cons1613)
    rule4196 = ReplacementRule(pattern4196, replacement4196)

    pattern4197 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons20, cons89, cons1613)
    rule4197 = ReplacementRule(pattern4197, replacement4197)

    pattern4198 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons81)
    rule4198 = ReplacementRule(pattern4198, replacement4198)

    pattern4199 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons81)
    rule4199 = ReplacementRule(pattern4199, replacement4199)

    pattern4200 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1614)
    rule4200 = ReplacementRule(pattern4200, replacement4200)

    pattern4201 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1614)
    rule4201 = ReplacementRule(pattern4201, replacement4201)

    pattern4202 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4202 = ReplacementRule(pattern4202, replacement4202)

    pattern4203 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4203 = ReplacementRule(pattern4203, replacement4203)

    pattern4204 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4204 = ReplacementRule(pattern4204, replacement4204)

    pattern4205 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4205 = ReplacementRule(pattern4205, replacement4205)

    pattern4206 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1615)
    rule4206 = ReplacementRule(pattern4206, replacement4206)

    pattern4207 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1615)
    rule4207 = ReplacementRule(pattern4207, replacement4207)

    pattern4208 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4208 = ReplacementRule(pattern4208, replacement4208)

    pattern4209 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons89, cons1592)
    rule4209 = ReplacementRule(pattern4209, replacement4209)

    pattern4210 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1306, cons1257)
    rule4210 = ReplacementRule(pattern4210, replacement4210)

    pattern4211 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1306, cons1257)
    rule4211 = ReplacementRule(pattern4211, replacement4211)

    pattern4212 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267)
    rule4212 = ReplacementRule(pattern4212, replacement4212)

    pattern4213 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267)
    rule4213 = ReplacementRule(pattern4213, replacement4213)

    pattern4214 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72)
    rule4214 = ReplacementRule(pattern4214, replacement4214)

    pattern4215 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72)
    rule4215 = ReplacementRule(pattern4215, replacement4215)

    pattern4216 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1414)
    rule4216 = ReplacementRule(pattern4216, replacement4216)

    pattern4217 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1414)
    rule4217 = ReplacementRule(pattern4217, replacement4217)

    pattern4218 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4218 = ReplacementRule(pattern4218, replacement4218)

    pattern4219 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4219 = ReplacementRule(pattern4219, replacement4219)

    pattern4220 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4220 = ReplacementRule(pattern4220, replacement4220)

    pattern4221 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4221 = ReplacementRule(pattern4221, replacement4221)

    pattern4222 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons168, cons1267, cons517)
    rule4222 = ReplacementRule(pattern4222, replacement4222)

    pattern4223 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons168, cons1267, cons517)
    rule4223 = ReplacementRule(pattern4223, replacement4223)

    pattern4224 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons168, cons1269, cons517)
    rule4224 = ReplacementRule(pattern4224, replacement4224)

    pattern4225 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons168, cons1269, cons517)
    rule4225 = ReplacementRule(pattern4225, replacement4225)

    pattern4226 = Pattern(Integral((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4226 = ReplacementRule(pattern4226, replacement4226)

    pattern4227 = Pattern(Integral((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4227 = ReplacementRule(pattern4227, replacement4227)

    pattern4228 = Pattern(Integral((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4228 = ReplacementRule(pattern4228, replacement4228)

    pattern4229 = Pattern(Integral((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4229 = ReplacementRule(pattern4229, replacement4229)

    pattern4230 = Pattern(Integral((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4230 = ReplacementRule(pattern4230, replacement4230)

    pattern4231 = Pattern(Integral((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4231 = ReplacementRule(pattern4231, replacement4231)

    pattern4232 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1267, cons517)
    rule4232 = ReplacementRule(pattern4232, replacement4232)

    pattern4233 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1267, cons517)
    rule4233 = ReplacementRule(pattern4233, replacement4233)

    pattern4234 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1269, cons517)
    rule4234 = ReplacementRule(pattern4234, replacement4234)

    pattern4235 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1269, cons517)
    rule4235 = ReplacementRule(pattern4235, replacement4235)

    pattern4236 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons79)
    rule4236 = ReplacementRule(pattern4236, replacement4236)

    pattern4237 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons79)
    rule4237 = ReplacementRule(pattern4237, replacement4237)

    pattern4238 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4238 = ReplacementRule(pattern4238, replacement4238)

    pattern4239 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4239 = ReplacementRule(pattern4239, replacement4239)

    pattern4240 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4240 = ReplacementRule(pattern4240, replacement4240)

    pattern4241 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4241 = ReplacementRule(pattern4241, replacement4241)

    pattern4242 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4242 = ReplacementRule(pattern4242, replacement4242)

    pattern4243 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4243 = ReplacementRule(pattern4243, replacement4243)

    pattern4244 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4244 = ReplacementRule(pattern4244, replacement4244)

    pattern4245 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4245 = ReplacementRule(pattern4245, replacement4245)

    pattern4246 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4246 = ReplacementRule(pattern4246, replacement4246)

    pattern4247 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4247 = ReplacementRule(pattern4247, replacement4247)

    pattern4248 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4248 = ReplacementRule(pattern4248, replacement4248)

    pattern4249 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule4249 = ReplacementRule(pattern4249, replacement4249)

    pattern4250 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4250 = ReplacementRule(pattern4250, replacement4250)

    pattern4251 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4251 = ReplacementRule(pattern4251, replacement4251)

    pattern4252 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4252 = ReplacementRule(pattern4252, replacement4252)

    pattern4253 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4253 = ReplacementRule(pattern4253, replacement4253)

    pattern4254 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4254 = ReplacementRule(pattern4254, replacement4254)

    pattern4255 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4255 = ReplacementRule(pattern4255, replacement4255)

    pattern4256 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4256 = ReplacementRule(pattern4256, replacement4256)

    pattern4257 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4257 = ReplacementRule(pattern4257, replacement4257)

    pattern4258 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4258 = ReplacementRule(pattern4258, replacement4258)

    pattern4259 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4259 = ReplacementRule(pattern4259, replacement4259)

    pattern4260 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4260 = ReplacementRule(pattern4260, replacement4260)

    pattern4261 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4261 = ReplacementRule(pattern4261, replacement4261)

    pattern4262 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4262 = ReplacementRule(pattern4262, replacement4262)

    pattern4263 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule4263 = ReplacementRule(pattern4263, replacement4263)

    pattern4264 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4264 = ReplacementRule(pattern4264, replacement4264)

    pattern4265 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4265 = ReplacementRule(pattern4265, replacement4265)

    pattern4266 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1325)
    rule4266 = ReplacementRule(pattern4266, replacement4266)

    pattern4267 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1325)
    rule4267 = ReplacementRule(pattern4267, replacement4267)

    pattern4268 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1267, cons1325, cons1306)
    rule4268 = ReplacementRule(pattern4268, replacement4268)

    pattern4269 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1267, cons1325, cons1306)
    rule4269 = ReplacementRule(pattern4269, replacement4269)

    pattern4270 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons20, cons87, cons1616)
    rule4270 = ReplacementRule(pattern4270, replacement4270)

    pattern4271 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons20, cons87, cons1616)
    rule4271 = ReplacementRule(pattern4271, replacement4271)

    pattern4272 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons81, cons82, cons1616)
    rule4272 = ReplacementRule(pattern4272, replacement4272)

    pattern4273 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons81, cons82, cons1616)
    rule4273 = ReplacementRule(pattern4273, replacement4273)

    pattern4274 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1257, cons79)
    rule4274 = ReplacementRule(pattern4274, replacement4274)

    pattern4275 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1257, cons79)
    rule4275 = ReplacementRule(pattern4275, replacement4275)

    pattern4276 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule4276 = ReplacementRule(pattern4276, replacement4276)

    pattern4277 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons150)
    rule4277 = ReplacementRule(pattern4277, replacement4277)

    pattern4278 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule4278 = ReplacementRule(pattern4278, replacement4278)

    pattern4279 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule4279 = ReplacementRule(pattern4279, replacement4279)

    pattern4280 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons25, cons20)
    rule4280 = ReplacementRule(pattern4280, replacement4280)

    pattern4281 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons25, cons20)
    rule4281 = ReplacementRule(pattern4281, replacement4281)

    pattern4282 = Pattern(Integral(((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25)
    rule4282 = ReplacementRule(pattern4282, replacement4282)

    pattern4283 = Pattern(Integral(((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25)
    rule4283 = ReplacementRule(pattern4283, replacement4283)

    pattern4284 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons157, cons1423)
    rule4284 = ReplacementRule(pattern4284, replacement4284)

    pattern4285 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons157, cons1423)
    rule4285 = ReplacementRule(pattern4285, replacement4285)

    pattern4286 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1317, cons1423, cons1233, cons1617)
    rule4286 = ReplacementRule(pattern4286, replacement4286)

    pattern4287 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1317, cons1423, cons1233, cons1617)
    rule4287 = ReplacementRule(pattern4287, replacement4287)

    pattern4288 = Pattern(Integral(sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule4288 = ReplacementRule(pattern4288, replacement4288)

    pattern4289 = Pattern(Integral(sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule4289 = ReplacementRule(pattern4289, replacement4289)

    pattern4290 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1316)
    rule4290 = ReplacementRule(pattern4290, replacement4290)

    pattern4291 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1316)
    rule4291 = ReplacementRule(pattern4291, replacement4291)

    pattern4292 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1268, cons33, cons1322)
    rule4292 = ReplacementRule(pattern4292, replacement4292)

    pattern4293 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1268, cons33, cons1322)
    rule4293 = ReplacementRule(pattern4293, replacement4293)

    pattern4294 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1268, cons1323, cons1618)
    rule4294 = ReplacementRule(pattern4294, replacement4294)

    pattern4295 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1268, cons1323, cons1618)
    rule4295 = ReplacementRule(pattern4295, replacement4295)

    pattern4296 = Pattern(Integral((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons150)
    rule4296 = ReplacementRule(pattern4296, replacement4296)

    pattern4297 = Pattern(Integral((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons150)
    rule4297 = ReplacementRule(pattern4297, replacement4297)

    pattern4298 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons150, cons1322, cons517)
    rule4298 = ReplacementRule(pattern4298, replacement4298)

    pattern4299 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons150, cons1322, cons517)
    rule4299 = ReplacementRule(pattern4299, replacement4299)

    pattern4300 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons152, cons1619, cons1620)
    rule4300 = ReplacementRule(pattern4300, replacement4300)

    pattern4301 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons152, cons1619, cons1620)
    rule4301 = ReplacementRule(pattern4301, replacement4301)

    pattern4302 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons81)
    rule4302 = ReplacementRule(pattern4302, replacement4302)

    pattern4303 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons81)
    rule4303 = ReplacementRule(pattern4303, replacement4303)

    pattern4304 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1621)
    rule4304 = ReplacementRule(pattern4304, replacement4304)

    pattern4305 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1621)
    rule4305 = ReplacementRule(pattern4305, replacement4305)

    pattern4306 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267)
    rule4306 = ReplacementRule(pattern4306, replacement4306)

    pattern4307 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267)
    rule4307 = ReplacementRule(pattern4307, replacement4307)

    pattern4308 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons152, cons1619, cons1620)
    rule4308 = ReplacementRule(pattern4308, replacement4308)

    pattern4309 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons152, cons1619, cons1620)
    rule4309 = ReplacementRule(pattern4309, replacement4309)

    pattern4310 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons81)
    rule4310 = ReplacementRule(pattern4310, replacement4310)

    pattern4311 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons81)
    rule4311 = ReplacementRule(pattern4311, replacement4311)

    pattern4312 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267)
    rule4312 = ReplacementRule(pattern4312, replacement4312)

    pattern4313 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267)
    rule4313 = ReplacementRule(pattern4313, replacement4313)

    pattern4314 = Pattern(Integral(sqrt(WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4314 = ReplacementRule(pattern4314, replacement4314)

    pattern4315 = Pattern(Integral(sqrt(WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4315 = ReplacementRule(pattern4315, replacement4315)

    pattern4316 = Pattern(Integral(sqrt(WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4316 = ReplacementRule(pattern4316, replacement4316)

    pattern4317 = Pattern(Integral(sqrt(WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4317 = ReplacementRule(pattern4317, replacement4317)

    pattern4318 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4318 = ReplacementRule(pattern4318, replacement4318)

    pattern4319 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule4319 = ReplacementRule(pattern4319, replacement4319)

    pattern4320 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4320 = ReplacementRule(pattern4320, replacement4320)

    pattern4321 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4321 = ReplacementRule(pattern4321, replacement4321)

    pattern4322 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4322 = ReplacementRule(pattern4322, replacement4322)

    pattern4323 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4323 = ReplacementRule(pattern4323, replacement4323)

    pattern4324 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4324 = ReplacementRule(pattern4324, replacement4324)

    pattern4325 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4325 = ReplacementRule(pattern4325, replacement4325)

    pattern4326 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4326 = ReplacementRule(pattern4326, replacement4326)

    pattern4327 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4327 = ReplacementRule(pattern4327, replacement4327)

    pattern4328 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4328 = ReplacementRule(pattern4328, replacement4328)

    pattern4329 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4329 = ReplacementRule(pattern4329, replacement4329)

    pattern4330 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4330 = ReplacementRule(pattern4330, replacement4330)

    pattern4331 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4331 = ReplacementRule(pattern4331, replacement4331)

    pattern4332 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4332 = ReplacementRule(pattern4332, replacement4332)

    pattern4333 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4333 = ReplacementRule(pattern4333, replacement4333)

    pattern4334 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4334 = ReplacementRule(pattern4334, replacement4334)

    pattern4335 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4335 = ReplacementRule(pattern4335, replacement4335)

    pattern4336 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4336 = ReplacementRule(pattern4336, replacement4336)

    pattern4337 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1411)
    rule4337 = ReplacementRule(pattern4337, replacement4337)

    pattern4338 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4338 = ReplacementRule(pattern4338, replacement4338)

    pattern4339 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4339 = ReplacementRule(pattern4339, replacement4339)

    pattern4340 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4340 = ReplacementRule(pattern4340, replacement4340)

    pattern4341 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule4341 = ReplacementRule(pattern4341, replacement4341)

    pattern4342 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4342 = ReplacementRule(pattern4342, replacement4342)

    pattern4343 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(5)/2)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269)
    rule4343 = ReplacementRule(pattern4343, replacement4343)

    pattern4344 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4344 = ReplacementRule(pattern4344, replacement4344)

    pattern4345 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4345 = ReplacementRule(pattern4345, replacement4345)

    pattern4346 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4346 = ReplacementRule(pattern4346, replacement4346)

    pattern4347 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule4347 = ReplacementRule(pattern4347, replacement4347)

    pattern4348 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4348 = ReplacementRule(pattern4348, replacement4348)

    pattern4349 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4349 = ReplacementRule(pattern4349, replacement4349)

    pattern4350 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4350 = ReplacementRule(pattern4350, replacement4350)

    pattern4351 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule4351 = ReplacementRule(pattern4351, replacement4351)

    pattern4352 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4352 = ReplacementRule(pattern4352, replacement4352)

    pattern4353 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4353 = ReplacementRule(pattern4353, replacement4353)

    pattern4354 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4354 = ReplacementRule(pattern4354, replacement4354)

    pattern4355 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule4355 = ReplacementRule(pattern4355, replacement4355)

    pattern4356 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4356 = ReplacementRule(pattern4356, replacement4356)

    pattern4357 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule4357 = ReplacementRule(pattern4357, replacement4357)

    pattern4358 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1267, cons1325, cons1622)
    rule4358 = ReplacementRule(pattern4358, replacement4358)

    pattern4359 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1267, cons1325, cons1622)
    rule4359 = ReplacementRule(pattern4359, replacement4359)

    pattern4360 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons20, cons87)
    rule4360 = ReplacementRule(pattern4360, replacement4360)

    pattern4361 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons20, cons87)
    rule4361 = ReplacementRule(pattern4361, replacement4361)

    pattern4362 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons73, cons1623, cons20)
    rule4362 = ReplacementRule(pattern4362, replacement4362)

    pattern4363 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons73, cons1623, cons20)
    rule4363 = ReplacementRule(pattern4363, replacement4363)

    pattern4364 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1623, cons21)
    rule4364 = ReplacementRule(pattern4364, replacement4364)

    pattern4365 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1623, cons21)
    rule4365 = ReplacementRule(pattern4365, replacement4365)

    pattern4366 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1306, cons1609, cons40, cons1624)
    rule4366 = ReplacementRule(pattern4366, replacement4366)

    pattern4367 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1306, cons1609, cons40, cons1624)
    rule4367 = ReplacementRule(pattern4367, replacement4367)

    pattern4368 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1417)
    rule4368 = ReplacementRule(pattern4368, replacement4368)

    pattern4369 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1417)
    rule4369 = ReplacementRule(pattern4369, replacement4369)

    pattern4370 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule4370 = ReplacementRule(pattern4370, replacement4370)

    pattern4371 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule4371 = ReplacementRule(pattern4371, replacement4371)

    pattern4372 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430)
    rule4372 = ReplacementRule(pattern4372, replacement4372)

    pattern4373 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430)
    rule4373 = ReplacementRule(pattern4373, replacement4373)

    pattern4374 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons89, cons1588)
    rule4374 = ReplacementRule(pattern4374, replacement4374)

    pattern4375 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons89, cons1588)
    rule4375 = ReplacementRule(pattern4375, replacement4375)

    pattern4376 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1571)
    rule4376 = ReplacementRule(pattern4376, replacement4376)

    pattern4377 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1571)
    rule4377 = ReplacementRule(pattern4377, replacement4377)

    pattern4378 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247)
    rule4378 = ReplacementRule(pattern4378, replacement4378)

    pattern4379 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247)
    rule4379 = ReplacementRule(pattern4379, replacement4379)

    pattern4380 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1267, cons1625)
    rule4380 = ReplacementRule(pattern4380, replacement4380)

    pattern4381 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1267, cons1625)
    rule4381 = ReplacementRule(pattern4381, replacement4381)

    pattern4382 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1267, cons1626, cons33, cons1322)
    rule4382 = ReplacementRule(pattern4382, replacement4382)

    pattern4383 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1267, cons1626, cons33, cons1322)
    rule4383 = ReplacementRule(pattern4383, replacement4383)

    pattern4384 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1267, cons1626, cons1323)
    rule4384 = ReplacementRule(pattern4384, replacement4384)

    pattern4385 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1267, cons1626, cons1323)
    rule4385 = ReplacementRule(pattern4385, replacement4385)

    pattern4386 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons33, cons170)
    rule4386 = ReplacementRule(pattern4386, replacement4386)

    pattern4387 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons33, cons170)
    rule4387 = ReplacementRule(pattern4387, replacement4387)

    pattern4388 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons33, cons96)
    rule4388 = ReplacementRule(pattern4388, replacement4388)

    pattern4389 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons33, cons96)
    rule4389 = ReplacementRule(pattern4389, replacement4389)

    pattern4390 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1269, cons1627)
    rule4390 = ReplacementRule(pattern4390, replacement4390)

    pattern4391 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1269, cons1627)
    rule4391 = ReplacementRule(pattern4391, replacement4391)

    pattern4392 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1269, cons1628)
    rule4392 = ReplacementRule(pattern4392, replacement4392)

    pattern4393 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1269, cons1628)
    rule4393 = ReplacementRule(pattern4393, replacement4393)

    pattern4394 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons1627, cons79)
    rule4394 = ReplacementRule(pattern4394, replacement4394)

    pattern4395 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons1247, cons1269, cons1627, cons79)
    rule4395 = ReplacementRule(pattern4395, replacement4395)

    pattern4396 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1269)
    rule4396 = ReplacementRule(pattern4396, replacement4396)

    pattern4397 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons36, cons37, cons50, cons127, cons19, cons1247, cons1269)
    rule4397 = ReplacementRule(pattern4397, replacement4397)

    pattern4398 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247, cons1267, cons33, cons1322)
    rule4398 = ReplacementRule(pattern4398, replacement4398)

    pattern4399 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247, cons1267, cons33, cons1322)
    rule4399 = ReplacementRule(pattern4399, replacement4399)

    pattern4400 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247, cons1269, cons33, cons96)
    rule4400 = ReplacementRule(pattern4400, replacement4400)

    pattern4401 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1247, cons1269, cons33, cons96)
    rule4401 = ReplacementRule(pattern4401, replacement4401)

    pattern4402 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons19, cons1247, cons274)
    rule4402 = ReplacementRule(pattern4402, replacement4402)

    pattern4403 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons19, cons1247, cons274)
    rule4403 = ReplacementRule(pattern4403, replacement4403)

    pattern4404 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1267, cons157, cons1629)
    rule4404 = ReplacementRule(pattern4404, replacement4404)

    pattern4405 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1267, cons157, cons1629)
    rule4405 = ReplacementRule(pattern4405, replacement4405)

    pattern4406 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons157, cons33, cons34)
    rule4406 = ReplacementRule(pattern4406, replacement4406)

    pattern4407 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons157, cons33, cons34)
    rule4407 = ReplacementRule(pattern4407, replacement4407)

    pattern4408 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1267, cons157, cons1551)
    rule4408 = ReplacementRule(pattern4408, replacement4408)

    pattern4409 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1267, cons157, cons1551)
    rule4409 = ReplacementRule(pattern4409, replacement4409)

    pattern4410 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons1630)
    rule4410 = ReplacementRule(pattern4410, replacement4410)

    pattern4411 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons1630)
    rule4411 = ReplacementRule(pattern4411, replacement4411)

    pattern4412 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons1631, cons89, cons465)
    rule4412 = ReplacementRule(pattern4412, replacement4412)

    pattern4413 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons1631, cons89, cons465)
    rule4413 = ReplacementRule(pattern4413, replacement4413)

    pattern4414 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons1631, cons1233)
    rule4414 = ReplacementRule(pattern4414, replacement4414)

    pattern4415 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons1631, cons1233)
    rule4415 = ReplacementRule(pattern4415, replacement4415)

    pattern4416 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons95, cons1425, cons91)
    rule4416 = ReplacementRule(pattern4416, replacement4416)

    pattern4417 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons95, cons1425, cons91)
    rule4417 = ReplacementRule(pattern4417, replacement4417)

    pattern4418 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons33, cons1425, cons348)
    rule4418 = ReplacementRule(pattern4418, replacement4418)

    pattern4419 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons33, cons1425, cons348)
    rule4419 = ReplacementRule(pattern4419, replacement4419)

    pattern4420 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons95, cons1322, cons90)
    rule4420 = ReplacementRule(pattern4420, replacement4420)

    pattern4421 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1267, cons95, cons1322, cons90)
    rule4421 = ReplacementRule(pattern4421, replacement4421)

    pattern4422 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons33, cons1322, cons1329)
    rule4422 = ReplacementRule(pattern4422, replacement4422)

    pattern4423 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1267, cons33, cons1322, cons1329)
    rule4423 = ReplacementRule(pattern4423, replacement4423)

    pattern4424 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267, cons89, cons167)
    rule4424 = ReplacementRule(pattern4424, replacement4424)

    pattern4425 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267, cons89, cons167)
    rule4425 = ReplacementRule(pattern4425, replacement4425)

    pattern4426 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267, cons89, cons465)
    rule4426 = ReplacementRule(pattern4426, replacement4426)

    pattern4427 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267, cons89, cons465)
    rule4427 = ReplacementRule(pattern4427, replacement4427)

    pattern4428 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267)
    rule4428 = ReplacementRule(pattern4428, replacement4428)

    pattern4429 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1267)
    rule4429 = ReplacementRule(pattern4429, replacement4429)

    pattern4430 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons168, cons1588)
    rule4430 = ReplacementRule(pattern4430, replacement4430)

    pattern4431 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons168, cons1588)
    rule4431 = ReplacementRule(pattern4431, replacement4431)

    pattern4432 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269, cons33, cons168, cons1632)
    rule4432 = ReplacementRule(pattern4432, replacement4432)

    pattern4433 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269, cons33, cons168, cons1632)
    rule4433 = ReplacementRule(pattern4433, replacement4433)

    pattern4434 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons96, cons1327)
    rule4434 = ReplacementRule(pattern4434, replacement4434)

    pattern4435 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons96, cons1327)
    rule4435 = ReplacementRule(pattern4435, replacement4435)

    pattern4436 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons96, cons167)
    rule4436 = ReplacementRule(pattern4436, replacement4436)

    pattern4437 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons96, cons167)
    rule4437 = ReplacementRule(pattern4437, replacement4437)

    pattern4438 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269, cons33, cons96, cons1633)
    rule4438 = ReplacementRule(pattern4438, replacement4438)

    pattern4439 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269, cons33, cons96, cons1633)
    rule4439 = ReplacementRule(pattern4439, replacement4439)

    pattern4440 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons1258, cons90)
    rule4440 = ReplacementRule(pattern4440, replacement4440)

    pattern4441 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons1258, cons90)
    rule4441 = ReplacementRule(pattern4441, replacement4441)

    pattern4442 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons1258, cons1588)
    rule4442 = ReplacementRule(pattern4442, replacement4442)

    pattern4443 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269, cons95, cons1258, cons1588)
    rule4443 = ReplacementRule(pattern4443, replacement4443)

    pattern4444 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1269, cons89, cons167, cons1361, cons1634)
    rule4444 = ReplacementRule(pattern4444, replacement4444)

    pattern4445 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1269, cons89, cons167, cons1361, cons1634)
    rule4445 = ReplacementRule(pattern4445, replacement4445)

    pattern4446 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1269, cons89, cons1588)
    rule4446 = ReplacementRule(pattern4446, replacement4446)

    pattern4447 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons1247, cons1269, cons89, cons1588)
    rule4447 = ReplacementRule(pattern4447, replacement4447)

    pattern4448 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4448 = ReplacementRule(pattern4448, replacement4448)

    pattern4449 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4449 = ReplacementRule(pattern4449, replacement4449)

    pattern4450 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4450 = ReplacementRule(pattern4450, replacement4450)

    pattern4451 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4451 = ReplacementRule(pattern4451, replacement4451)

    pattern4452 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4452 = ReplacementRule(pattern4452, replacement4452)

    pattern4453 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1247, cons1269)
    rule4453 = ReplacementRule(pattern4453, replacement4453)

    pattern4454 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269)
    rule4454 = ReplacementRule(pattern4454, replacement4454)

    pattern4455 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons4, cons1247, cons1269)
    rule4455 = ReplacementRule(pattern4455, replacement4455)

    pattern4456 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1269)
    rule4456 = ReplacementRule(pattern4456, replacement4456)

    pattern4457 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons1247, cons1269)
    rule4457 = ReplacementRule(pattern4457, replacement4457)

    pattern4458 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons5, cons72, cons1267, cons377)
    rule4458 = ReplacementRule(pattern4458, replacement4458)

    pattern4459 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons5, cons72, cons1267, cons377)
    rule4459 = ReplacementRule(pattern4459, replacement4459)

    pattern4460 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule4460 = ReplacementRule(pattern4460, replacement4460)

    pattern4461 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule4461 = ReplacementRule(pattern4461, replacement4461)

    pattern4462 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule4462 = ReplacementRule(pattern4462, replacement4462)

    pattern4463 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule4463 = ReplacementRule(pattern4463, replacement4463)

    pattern4464 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1635)
    rule4464 = ReplacementRule(pattern4464, replacement4464)

    pattern4465 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1635)
    rule4465 = ReplacementRule(pattern4465, replacement4465)

    pattern4466 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons1636, cons33, cons34)
    rule4466 = ReplacementRule(pattern4466, replacement4466)

    pattern4467 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons1636, cons33, cons34)
    rule4467 = ReplacementRule(pattern4467, replacement4467)

    pattern4468 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1636, cons1551)
    rule4468 = ReplacementRule(pattern4468, replacement4468)

    pattern4469 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1636, cons1551)
    rule4469 = ReplacementRule(pattern4469, replacement4469)

    pattern4470 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons37, cons38, cons19, cons1433)
    rule4470 = ReplacementRule(pattern4470, replacement4470)

    pattern4471 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons37, cons38, cons19, cons1433)
    rule4471 = ReplacementRule(pattern4471, replacement4471)

    pattern4472 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1637)
    rule4472 = ReplacementRule(pattern4472, replacement4472)

    pattern4473 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1637)
    rule4473 = ReplacementRule(pattern4473, replacement4473)

    pattern4474 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1638)
    rule4474 = ReplacementRule(pattern4474, replacement4474)

    pattern4475 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1638)
    rule4475 = ReplacementRule(pattern4475, replacement4475)

    pattern4476 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1639)
    rule4476 = ReplacementRule(pattern4476, replacement4476)

    pattern4477 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1639)
    rule4477 = ReplacementRule(pattern4477, replacement4477)

    pattern4478 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1638)
    rule4478 = ReplacementRule(pattern4478, replacement4478)

    pattern4479 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1638)
    rule4479 = ReplacementRule(pattern4479, replacement4479)

    pattern4480 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1639)
    rule4480 = ReplacementRule(pattern4480, replacement4480)

    pattern4481 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1639)
    rule4481 = ReplacementRule(pattern4481, replacement4481)

    pattern4482 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1267, cons33, cons1322)
    rule4482 = ReplacementRule(pattern4482, replacement4482)

    pattern4483 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1267, cons33, cons1322)
    rule4483 = ReplacementRule(pattern4483, replacement4483)

    pattern4484 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1267, cons33, cons1322)
    rule4484 = ReplacementRule(pattern4484, replacement4484)

    pattern4485 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1267, cons33, cons1322)
    rule4485 = ReplacementRule(pattern4485, replacement4485)

    pattern4486 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1267, cons1323)
    rule4486 = ReplacementRule(pattern4486, replacement4486)

    pattern4487 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1267, cons1323)
    rule4487 = ReplacementRule(pattern4487, replacement4487)

    pattern4488 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1267, cons1323)
    rule4488 = ReplacementRule(pattern4488, replacement4488)

    pattern4489 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1267, cons1323)
    rule4489 = ReplacementRule(pattern4489, replacement4489)

    pattern4490 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons1640)
    rule4490 = ReplacementRule(pattern4490, replacement4490)

    pattern4491 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons1640)
    rule4491 = ReplacementRule(pattern4491, replacement4491)

    pattern4492 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons1640)
    rule4492 = ReplacementRule(pattern4492, replacement4492)

    pattern4493 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons1640)
    rule4493 = ReplacementRule(pattern4493, replacement4493)

    pattern4494 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4494 = ReplacementRule(pattern4494, replacement4494)

    pattern4495 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4495 = ReplacementRule(pattern4495, replacement4495)

    pattern4496 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269)
    rule4496 = ReplacementRule(pattern4496, replacement4496)

    pattern4497 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269)
    rule4497 = ReplacementRule(pattern4497, replacement4497)

    pattern4498 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons517, cons96)
    rule4498 = ReplacementRule(pattern4498, replacement4498)

    pattern4499 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons33, cons96)
    rule4499 = ReplacementRule(pattern4499, replacement4499)

    pattern4500 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons517, cons96)
    rule4500 = ReplacementRule(pattern4500, replacement4500)

    pattern4501 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons517, cons96)
    rule4501 = ReplacementRule(pattern4501, replacement4501)

    pattern4502 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons79)
    rule4502 = ReplacementRule(pattern4502, replacement4502)

    pattern4503 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons79)
    rule4503 = ReplacementRule(pattern4503, replacement4503)

    pattern4504 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1269, cons79)
    rule4504 = ReplacementRule(pattern4504, replacement4504)

    pattern4505 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1269, cons79)
    rule4505 = ReplacementRule(pattern4505, replacement4505)

    pattern4506 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons21)
    rule4506 = ReplacementRule(pattern4506, replacement4506)

    pattern4507 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons21)
    rule4507 = ReplacementRule(pattern4507, replacement4507)

    pattern4508 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons21)
    rule4508 = ReplacementRule(pattern4508, replacement4508)

    pattern4509 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons21)
    rule4509 = ReplacementRule(pattern4509, replacement4509)

    pattern4510 = Pattern(Integral(((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons5, cons21)
    rule4510 = ReplacementRule(pattern4510, replacement4510)

    pattern4511 = Pattern(Integral(((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons5, cons21)
    rule4511 = ReplacementRule(pattern4511, replacement4511)

    pattern4512 = Pattern(Integral(((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons5, cons21)
    rule4512 = ReplacementRule(pattern4512, replacement4512)

    pattern4513 = Pattern(Integral(((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons5, cons21)
    rule4513 = ReplacementRule(pattern4513, replacement4513)

    pattern4514 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons89, cons91)
    rule4514 = ReplacementRule(pattern4514, replacement4514)

    pattern4515 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons89, cons91)
    rule4515 = ReplacementRule(pattern4515, replacement4515)

    pattern4516 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons89, cons91)
    rule4516 = ReplacementRule(pattern4516, replacement4516)

    pattern4517 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons89, cons91)
    rule4517 = ReplacementRule(pattern4517, replacement4517)

    pattern4518 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons348)
    rule4518 = ReplacementRule(pattern4518, replacement4518)

    pattern4519 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons348)
    rule4519 = ReplacementRule(pattern4519, replacement4519)

    pattern4520 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons348)
    rule4520 = ReplacementRule(pattern4520, replacement4520)

    pattern4521 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons348)
    rule4521 = ReplacementRule(pattern4521, replacement4521)

    pattern4522 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1267)
    rule4522 = ReplacementRule(pattern4522, replacement4522)

    pattern4523 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1267)
    rule4523 = ReplacementRule(pattern4523, replacement4523)

    pattern4524 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1267)
    rule4524 = ReplacementRule(pattern4524, replacement4524)

    pattern4525 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1267)
    rule4525 = ReplacementRule(pattern4525, replacement4525)

    pattern4526 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1269)
    rule4526 = ReplacementRule(pattern4526, replacement4526)

    pattern4527 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1269)
    rule4527 = ReplacementRule(pattern4527, replacement4527)

    pattern4528 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1269)
    rule4528 = ReplacementRule(pattern4528, replacement4528)

    pattern4529 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1269)
    rule4529 = ReplacementRule(pattern4529, replacement4529)

    pattern4530 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons274)
    rule4530 = ReplacementRule(pattern4530, replacement4530)

    pattern4531 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons274)
    rule4531 = ReplacementRule(pattern4531, replacement4531)

    pattern4532 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule4532 = ReplacementRule(pattern4532, replacement4532)

    pattern4533 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule4533 = ReplacementRule(pattern4533, replacement4533)

    pattern4534 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1267, cons33, cons1322)
    rule4534 = ReplacementRule(pattern4534, replacement4534)

    pattern4535 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1267, cons33, cons1322)
    rule4535 = ReplacementRule(pattern4535, replacement4535)

    pattern4536 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1267, cons33, cons1322)
    rule4536 = ReplacementRule(pattern4536, replacement4536)

    pattern4537 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1267, cons33, cons1322)
    rule4537 = ReplacementRule(pattern4537, replacement4537)

    pattern4538 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1267, cons1323, cons1641)
    rule4538 = ReplacementRule(pattern4538, replacement4538)

    pattern4539 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1267, cons1323, cons1641)
    rule4539 = ReplacementRule(pattern4539, replacement4539)

    pattern4540 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1267, cons1323, cons1641)
    rule4540 = ReplacementRule(pattern4540, replacement4540)

    pattern4541 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1267, cons1323, cons1641)
    rule4541 = ReplacementRule(pattern4541, replacement4541)

    pattern4542 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons1267, cons1323, cons1642, cons685)
    rule4542 = ReplacementRule(pattern4542, replacement4542)

    pattern4543 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons1267, cons1323, cons1642, cons685)
    rule4543 = ReplacementRule(pattern4543, replacement4543)

    pattern4544 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1267, cons1323, cons1642, cons685)
    rule4544 = ReplacementRule(pattern4544, replacement4544)

    pattern4545 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1267, cons1323, cons1642, cons685)
    rule4545 = ReplacementRule(pattern4545, replacement4545)

    pattern4546 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons33, cons96)
    rule4546 = ReplacementRule(pattern4546, replacement4546)

    pattern4547 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1269, cons33, cons96)
    rule4547 = ReplacementRule(pattern4547, replacement4547)

    pattern4548 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons33, cons96)
    rule4548 = ReplacementRule(pattern4548, replacement4548)

    pattern4549 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1269, cons33, cons96)
    rule4549 = ReplacementRule(pattern4549, replacement4549)

    pattern4550 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons274)
    rule4550 = ReplacementRule(pattern4550, replacement4550)

    pattern4551 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons274)
    rule4551 = ReplacementRule(pattern4551, replacement4551)

    pattern4552 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1269, cons274)
    rule4552 = ReplacementRule(pattern4552, replacement4552)

    pattern4553 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1269, cons274)
    rule4553 = ReplacementRule(pattern4553, replacement4553)

    pattern4554 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269, cons95, cons170, cons1588)
    rule4554 = ReplacementRule(pattern4554, replacement4554)

    pattern4555 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269, cons95, cons170, cons1588)
    rule4555 = ReplacementRule(pattern4555, replacement4555)

    pattern4556 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269, cons95, cons170, cons1588)
    rule4556 = ReplacementRule(pattern4556, replacement4556)

    pattern4557 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269, cons95, cons170, cons1588)
    rule4557 = ReplacementRule(pattern4557, replacement4557)

    pattern4558 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1269, cons33, cons170, cons1571)
    rule4558 = ReplacementRule(pattern4558, replacement4558)

    pattern4559 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1269, cons33, cons170, cons1571)
    rule4559 = ReplacementRule(pattern4559, replacement4559)

    pattern4560 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1269, cons33, cons170, cons1571)
    rule4560 = ReplacementRule(pattern4560, replacement4560)

    pattern4561 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1269, cons33, cons170, cons1571)
    rule4561 = ReplacementRule(pattern4561, replacement4561)

    pattern4562 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269, cons95, cons96, cons90)
    rule4562 = ReplacementRule(pattern4562, replacement4562)

    pattern4563 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269, cons95, cons96, cons90)
    rule4563 = ReplacementRule(pattern4563, replacement4563)

    pattern4564 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269, cons95, cons96, cons90)
    rule4564 = ReplacementRule(pattern4564, replacement4564)

    pattern4565 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269, cons95, cons96, cons90)
    rule4565 = ReplacementRule(pattern4565, replacement4565)

    pattern4566 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1269, cons33, cons96, cons1633)
    rule4566 = ReplacementRule(pattern4566, replacement4566)

    pattern4567 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons1269, cons33, cons96, cons1633)
    rule4567 = ReplacementRule(pattern4567, replacement4567)

    pattern4568 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1269, cons33, cons96, cons1633)
    rule4568 = ReplacementRule(pattern4568, replacement4568)

    pattern4569 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons1269, cons33, cons96, cons1633)
    rule4569 = ReplacementRule(pattern4569, replacement4569)

    pattern4570 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons89, cons90)
    rule4570 = ReplacementRule(pattern4570, replacement4570)

    pattern4571 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons89, cons90)
    rule4571 = ReplacementRule(pattern4571, replacement4571)

    pattern4572 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1269, cons89, cons90)
    rule4572 = ReplacementRule(pattern4572, replacement4572)

    pattern4573 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1269, cons89, cons90)
    rule4573 = ReplacementRule(pattern4573, replacement4573)

    pattern4574 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons89, cons1588)
    rule4574 = ReplacementRule(pattern4574, replacement4574)

    pattern4575 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons1269, cons89, cons1588)
    rule4575 = ReplacementRule(pattern4575, replacement4575)

    pattern4576 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1269, cons89, cons1588)
    rule4576 = ReplacementRule(pattern4576, replacement4576)

    pattern4577 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons1269, cons89, cons1588)
    rule4577 = ReplacementRule(pattern4577, replacement4577)

    pattern4578 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4578 = ReplacementRule(pattern4578, replacement4578)

    pattern4579 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4579 = ReplacementRule(pattern4579, replacement4579)

    pattern4580 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule4580 = ReplacementRule(pattern4580, replacement4580)

    pattern4581 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule4581 = ReplacementRule(pattern4581, replacement4581)

    pattern4582 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4582 = ReplacementRule(pattern4582, replacement4582)

    pattern4583 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule4583 = ReplacementRule(pattern4583, replacement4583)

    pattern4584 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule4584 = ReplacementRule(pattern4584, replacement4584)

    pattern4585 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule4585 = ReplacementRule(pattern4585, replacement4585)

    pattern4586 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons1643)
    rule4586 = ReplacementRule(pattern4586, replacement4586)

    pattern4587 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons1643)
    rule4587 = ReplacementRule(pattern4587, replacement4587)

    pattern4588 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1644)
    rule4588 = ReplacementRule(pattern4588, replacement4588)

    pattern4589 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1644)
    rule4589 = ReplacementRule(pattern4589, replacement4589)

    pattern4590 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons25, cons20)
    rule4590 = ReplacementRule(pattern4590, replacement4590)

    pattern4591 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons25, cons20)
    rule4591 = ReplacementRule(pattern4591, replacement4591)

    pattern4592 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons25, cons20)
    rule4592 = ReplacementRule(pattern4592, replacement4592)

    pattern4593 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons4, cons25, cons20)
    rule4593 = ReplacementRule(pattern4593, replacement4593)

    pattern4594 = Pattern(Integral(((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons25)
    rule4594 = ReplacementRule(pattern4594, replacement4594)

    pattern4595 = Pattern(Integral(((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons25)
    rule4595 = ReplacementRule(pattern4595, replacement4595)

    pattern4596 = Pattern(Integral(((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons5, cons25)
    rule4596 = ReplacementRule(pattern4596, replacement4596)

    pattern4597 = Pattern(Integral(((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons5, cons25)
    rule4597 = ReplacementRule(pattern4597, replacement4597)

    pattern4598 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons3, cons8, cons29, cons4, cons1645)
    rule4598 = ReplacementRule(pattern4598, replacement4598)

    pattern4599 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**n_, x_), cons3, cons8, cons29, cons4, cons1645)
    rule4599 = ReplacementRule(pattern4599, replacement4599)

    pattern4600 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456)
    rule4600 = ReplacementRule(pattern4600, replacement4600)

    pattern4601 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456)
    rule4601 = ReplacementRule(pattern4601, replacement4601)

    pattern4602 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1480)
    rule4602 = ReplacementRule(pattern4602, replacement4602)

    pattern4603 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1480)
    rule4603 = ReplacementRule(pattern4603, replacement4603)

    pattern4604 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1480, cons56)
    rule4604 = ReplacementRule(pattern4604, replacement4604)

    pattern4605 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1480, cons56)
    rule4605 = ReplacementRule(pattern4605, replacement4605)

    pattern4606 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons5, cons1482, cons1481)
    rule4606 = ReplacementRule(pattern4606, With4606)

    pattern4607 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons5, cons1482, cons1481)
    rule4607 = ReplacementRule(pattern4607, With4607)

    pattern4608 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons378)
    rule4608 = ReplacementRule(pattern4608, With4608)

    pattern4609 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons378)
    rule4609 = ReplacementRule(pattern4609, With4609)

    pattern4610 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons5, cons1482, cons1481)
    rule4610 = ReplacementRule(pattern4610, With4610)

    pattern4611 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons5, cons1482, cons1481)
    rule4611 = ReplacementRule(pattern4611, With4611)

    pattern4612 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons1481, cons40)
    rule4612 = ReplacementRule(pattern4612, With4612)

    pattern4613 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons1481, cons40)
    rule4613 = ReplacementRule(pattern4613, With4613)

    pattern4614 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons377)
    rule4614 = ReplacementRule(pattern4614, replacement4614)

    pattern4615 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons377)
    rule4615 = ReplacementRule(pattern4615, replacement4615)

    pattern4616 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1483, cons87, cons40)
    rule4616 = ReplacementRule(pattern4616, With4616)

    pattern4617 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1483, cons87, cons40)
    rule4617 = ReplacementRule(pattern4617, With4617)

    pattern4618 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1482, cons1481)
    rule4618 = ReplacementRule(pattern4618, With4618)

    pattern4619 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1482, cons1481)
    rule4619 = ReplacementRule(pattern4619, With4619)

    pattern4620 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule4620 = ReplacementRule(pattern4620, replacement4620)

    pattern4621 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule4621 = ReplacementRule(pattern4621, replacement4621)

    pattern4622 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons47, cons149)
    rule4622 = ReplacementRule(pattern4622, replacement4622)

    pattern4623 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons47, cons149)
    rule4623 = ReplacementRule(pattern4623, replacement4623)

    pattern4624 = Pattern(Integral(S(1)/((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule4624 = ReplacementRule(pattern4624, With4624)

    pattern4625 = Pattern(Integral(S(1)/((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule4625 = ReplacementRule(pattern4625, With4625)

    pattern4626 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n2_*WC('c', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons378)
    rule4626 = ReplacementRule(pattern4626, With4626)

    pattern4627 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n2_*WC('c', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons378)
    rule4627 = ReplacementRule(pattern4627, With4627)

    pattern4628 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n2_*WC('c', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons48, cons1482, cons1481)
    rule4628 = ReplacementRule(pattern4628, With4628)

    pattern4629 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n2_*WC('c', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons5, cons48, cons1482, cons1481)
    rule4629 = ReplacementRule(pattern4629, With4629)

    pattern4630 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule4630 = ReplacementRule(pattern4630, replacement4630)

    pattern4631 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule4631 = ReplacementRule(pattern4631, replacement4631)

    pattern4632 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule4632 = ReplacementRule(pattern4632, replacement4632)

    pattern4633 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule4633 = ReplacementRule(pattern4633, replacement4633)

    pattern4634 = Pattern(Integral(((S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons377)
    rule4634 = ReplacementRule(pattern4634, replacement4634)

    pattern4635 = Pattern(Integral(((S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons377)
    rule4635 = ReplacementRule(pattern4635, replacement4635)

    pattern4636 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons1483, cons87, cons40)
    rule4636 = ReplacementRule(pattern4636, With4636)

    pattern4637 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons1483, cons87, cons40)
    rule4637 = ReplacementRule(pattern4637, With4637)

    pattern4638 = Pattern(Integral((a_ + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons1481)
    rule4638 = ReplacementRule(pattern4638, With4638)

    pattern4639 = Pattern(Integral((a_ + (S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons1481)
    rule4639 = ReplacementRule(pattern4639, With4639)

    pattern4640 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule4640 = ReplacementRule(pattern4640, replacement4640)

    pattern4641 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule4641 = ReplacementRule(pattern4641, replacement4641)

    pattern4642 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule4642 = ReplacementRule(pattern4642, replacement4642)

    pattern4643 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule4643 = ReplacementRule(pattern4643, replacement4643)

    pattern4644 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule4644 = ReplacementRule(pattern4644, With4644)

    pattern4645 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule4645 = ReplacementRule(pattern4645, With4645)

    pattern4646 = Pattern(Integral((A_ + WC('B', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule4646 = ReplacementRule(pattern4646, replacement4646)

    pattern4647 = Pattern(Integral((A_ + WC('B', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule4647 = ReplacementRule(pattern4647, replacement4647)

    pattern4648 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons64)
    rule4648 = ReplacementRule(pattern4648, replacement4648)

    pattern4649 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons64)
    rule4649 = ReplacementRule(pattern4649, replacement4649)

    pattern4650 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule4650 = ReplacementRule(pattern4650, replacement4650)

    pattern4651 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule4651 = ReplacementRule(pattern4651, replacement4651)

    pattern4652 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167, cons1646)
    rule4652 = ReplacementRule(pattern4652, replacement4652)

    pattern4653 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167, cons1646)
    rule4653 = ReplacementRule(pattern4653, replacement4653)

    pattern4654 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons1646, cons168)
    rule4654 = ReplacementRule(pattern4654, replacement4654)

    pattern4655 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons1646, cons168)
    rule4655 = ReplacementRule(pattern4655, replacement4655)

    pattern4656 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91)
    rule4656 = ReplacementRule(pattern4656, replacement4656)

    pattern4657 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91)
    rule4657 = ReplacementRule(pattern4657, replacement4657)

    pattern4658 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons168)
    rule4658 = ReplacementRule(pattern4658, replacement4658)

    pattern4659 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons168)
    rule4659 = ReplacementRule(pattern4659, replacement4659)

    pattern4660 = Pattern(Integral((WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons25)
    rule4660 = ReplacementRule(pattern4660, replacement4660)

    pattern4661 = Pattern(Integral((WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons25)
    rule4661 = ReplacementRule(pattern4661, replacement4661)

    pattern4662 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule4662 = ReplacementRule(pattern4662, replacement4662)

    pattern4663 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule4663 = ReplacementRule(pattern4663, replacement4663)

    pattern4664 = Pattern(Integral((a_ + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons198, cons64)
    rule4664 = ReplacementRule(pattern4664, replacement4664)

    pattern4665 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons198, cons64)
    rule4665 = ReplacementRule(pattern4665, replacement4665)

    pattern4666 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule4666 = ReplacementRule(pattern4666, replacement4666)

    pattern4667 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule4667 = ReplacementRule(pattern4667, replacement4667)

    pattern4668 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule4668 = ReplacementRule(pattern4668, replacement4668)

    pattern4669 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule4669 = ReplacementRule(pattern4669, replacement4669)

    pattern4670 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule4670 = ReplacementRule(pattern4670, replacement4670)

    pattern4671 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule4671 = ReplacementRule(pattern4671, replacement4671)

    pattern4672 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule4672 = ReplacementRule(pattern4672, replacement4672)

    pattern4673 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule4673 = ReplacementRule(pattern4673, replacement4673)

    pattern4674 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule4674 = ReplacementRule(pattern4674, replacement4674)

    pattern4675 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule4675 = ReplacementRule(pattern4675, replacement4675)

    pattern4676 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule4676 = ReplacementRule(pattern4676, replacement4676)

    pattern4677 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule4677 = ReplacementRule(pattern4677, replacement4677)

    pattern4678 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule4678 = ReplacementRule(pattern4678, replacement4678)

    pattern4679 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule4679 = ReplacementRule(pattern4679, replacement4679)

    pattern4680 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule4680 = ReplacementRule(pattern4680, replacement4680)

    pattern4681 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule4681 = ReplacementRule(pattern4681, replacement4681)

    pattern4682 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule4682 = ReplacementRule(pattern4682, replacement4682)

    pattern4683 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule4683 = ReplacementRule(pattern4683, replacement4683)

    pattern4684 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule4684 = ReplacementRule(pattern4684, replacement4684)

    pattern4685 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule4685 = ReplacementRule(pattern4685, replacement4685)

    pattern4686 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**p_*sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1647)
    rule4686 = ReplacementRule(pattern4686, replacement4686)

    pattern4687 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**p_*cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1647)
    rule4687 = ReplacementRule(pattern4687, replacement4687)
    return [rule3920, rule3921, rule3922, rule3923, rule3924, rule3925, rule3926, rule3927, rule3928, rule3929, rule3930, rule3931, rule3932, rule3933, rule3934, rule3935, rule3936, rule3937, rule3938, rule3939, rule3940, rule3941, rule3942, rule3943, rule3944, rule3945, rule3946, rule3947, rule3948, rule3949, rule3950, rule3951, rule3952, rule3953, rule3954, rule3955, rule3956, rule3957, rule3958, rule3959, rule3960, rule3961, rule3962, rule3963, rule3964, rule3965, rule3966, rule3967, rule3968, rule3969, rule3970, rule3971, rule3972, rule3973, rule3974, rule3975, rule3976, rule3977, rule3978, rule3979, rule3980, rule3981, rule3982, rule3983, rule3984, rule3985, rule3986, rule3987, rule3988, rule3989, rule3990, rule3991, rule3992, rule3993, rule3994, rule3995, rule3996, rule3997, rule3998, rule3999, rule4000, rule4001, rule4002, rule4003, rule4004, rule4005, rule4006, rule4007, rule4008, rule4009, rule4010, rule4011, rule4012, rule4013, rule4014, rule4015, rule4016, rule4017, rule4018, rule4019, rule4020, rule4021, rule4022, rule4023, rule4024, rule4025, rule4026, rule4027, rule4028, rule4029, rule4030, rule4031, rule4032, rule4033, rule4034, rule4035, rule4036, rule4037, rule4038, rule4039, rule4040, rule4041, rule4042, rule4043, rule4044, rule4045, rule4046, rule4047, rule4048, rule4049, rule4050, rule4051, rule4052, rule4053, rule4054, rule4055, rule4056, rule4057, rule4058, rule4059, rule4060, rule4061, rule4062, rule4063, rule4064, rule4065, rule4066, rule4067, rule4068, rule4069, rule4070, rule4071, rule4072, rule4073, rule4074, rule4075, rule4076, rule4077, rule4078, rule4079, rule4080, rule4081, rule4082, rule4083, rule4084, rule4085, rule4086, rule4087, rule4088, rule4089, rule4090, rule4091, rule4092, rule4093, rule4094, rule4095, rule4096, rule4097, rule4098, rule4099, rule4100, rule4101, rule4102, rule4103, rule4104, rule4105, rule4106, rule4107, rule4108, rule4109, rule4110, rule4111, rule4112, rule4113, rule4114, rule4115, rule4116, rule4117, rule4118, rule4119, rule4120, rule4121, rule4122, rule4123, rule4124, rule4125, rule4126, rule4127, rule4128, rule4129, rule4130, rule4131, rule4132, rule4133, rule4134, rule4135, rule4136, rule4137, rule4138, rule4139, rule4140, rule4141, rule4142, rule4143, rule4144, rule4145, rule4146, rule4147, rule4148, rule4149, rule4150, rule4151, rule4152, rule4153, rule4154, rule4155, rule4156, rule4157, rule4158, rule4159, rule4160, rule4161, rule4162, rule4163, rule4164, rule4165, rule4166, rule4167, rule4168, rule4169, rule4170, rule4171, rule4172, rule4173, rule4174, rule4175, rule4176, rule4177, rule4178, rule4179, rule4180, rule4181, rule4182, rule4183, rule4184, rule4185, rule4186, rule4187, rule4188, rule4189, rule4190, rule4191, rule4192, rule4193, rule4194, rule4195, rule4196, rule4197, rule4198, rule4199, rule4200, rule4201, rule4202, rule4203, rule4204, rule4205, rule4206, rule4207, rule4208, rule4209, rule4210, rule4211, rule4212, rule4213, rule4214, rule4215, rule4216, rule4217, rule4218, rule4219, rule4220, rule4221, rule4222, rule4223, rule4224, rule4225, rule4226, rule4227, rule4228, rule4229, rule4230, rule4231, rule4232, rule4233, rule4234, rule4235, rule4236, rule4237, rule4238, rule4239, rule4240, rule4241, rule4242, rule4243, rule4244, rule4245, rule4246, rule4247, rule4248, rule4249, rule4250, rule4251, rule4252, rule4253, rule4254, rule4255, rule4256, rule4257, rule4258, rule4259, rule4260, rule4261, rule4262, rule4263, rule4264, rule4265, rule4266, rule4267, rule4268, rule4269, rule4270, rule4271, rule4272, rule4273, rule4274, rule4275, rule4276, rule4277, rule4278, rule4279, rule4280, rule4281, rule4282, rule4283, rule4284, rule4285, rule4286, rule4287, rule4288, rule4289, rule4290, rule4291, rule4292, rule4293, rule4294, rule4295, rule4296, rule4297, rule4298, rule4299, rule4300, rule4301, rule4302, rule4303, rule4304, rule4305, rule4306, rule4307, rule4308, rule4309, rule4310, rule4311, rule4312, rule4313, rule4314, rule4315, rule4316, rule4317, rule4318, rule4319, rule4320, rule4321, rule4322, rule4323, rule4324, rule4325, rule4326, rule4327, rule4328, rule4329, rule4330, rule4331, rule4332, rule4333, rule4334, rule4335, rule4336, rule4337, rule4338, rule4339, rule4340, rule4341, rule4342, rule4343, rule4344, rule4345, rule4346, rule4347, rule4348, rule4349, rule4350, rule4351, rule4352, rule4353, rule4354, rule4355, rule4356, rule4357, rule4358, rule4359, rule4360, rule4361, rule4362, rule4363, rule4364, rule4365, rule4366, rule4367, rule4368, rule4369, rule4370, rule4371, rule4372, rule4373, rule4374, rule4375, rule4376, rule4377, rule4378, rule4379, rule4380, rule4381, rule4382, rule4383, rule4384, rule4385, rule4386, rule4387, rule4388, rule4389, rule4390, rule4391, rule4392, rule4393, rule4394, rule4395, rule4396, rule4397, rule4398, rule4399, rule4400, rule4401, rule4402, rule4403, rule4404, rule4405, rule4406, rule4407, rule4408, rule4409, rule4410, rule4411, rule4412, rule4413, rule4414, rule4415, rule4416, rule4417, rule4418, rule4419, rule4420, rule4421, rule4422, rule4423, rule4424, rule4425, rule4426, rule4427, rule4428, rule4429, rule4430, rule4431, rule4432, rule4433, rule4434, rule4435, rule4436, rule4437, rule4438, rule4439, rule4440, rule4441, rule4442, rule4443, rule4444, rule4445, rule4446, rule4447, rule4448, rule4449, rule4450, rule4451, rule4452, rule4453, rule4454, rule4455, rule4456, rule4457, rule4458, rule4459, rule4460, rule4461, rule4462, rule4463, rule4464, rule4465, rule4466, rule4467, rule4468, rule4469, rule4470, rule4471, rule4472, rule4473, rule4474, rule4475, rule4476, rule4477, rule4478, rule4479, rule4480, rule4481, rule4482, rule4483, rule4484, rule4485, rule4486, rule4487, rule4488, rule4489, rule4490, rule4491, rule4492, rule4493, rule4494, rule4495, rule4496, rule4497, rule4498, rule4499, rule4500, rule4501, rule4502, rule4503, rule4504, rule4505, rule4506, rule4507, rule4508, rule4509, rule4510, rule4511, rule4512, rule4513, rule4514, rule4515, rule4516, rule4517, rule4518, rule4519, rule4520, rule4521, rule4522, rule4523, rule4524, rule4525, rule4526, rule4527, rule4528, rule4529, rule4530, rule4531, rule4532, rule4533, rule4534, rule4535, rule4536, rule4537, rule4538, rule4539, rule4540, rule4541, rule4542, rule4543, rule4544, rule4545, rule4546, rule4547, rule4548, rule4549, rule4550, rule4551, rule4552, rule4553, rule4554, rule4555, rule4556, rule4557, rule4558, rule4559, rule4560, rule4561, rule4562, rule4563, rule4564, rule4565, rule4566, rule4567, rule4568, rule4569, rule4570, rule4571, rule4572, rule4573, rule4574, rule4575, rule4576, rule4577, rule4578, rule4579, rule4580, rule4581, rule4582, rule4583, rule4584, rule4585, rule4586, rule4587, rule4588, rule4589, rule4590, rule4591, rule4592, rule4593, rule4594, rule4595, rule4596, rule4597, rule4598, rule4599, rule4600, rule4601, rule4602, rule4603, rule4604, rule4605, rule4606, rule4607, rule4608, rule4609, rule4610, rule4611, rule4612, rule4613, rule4614, rule4615, rule4616, rule4617, rule4618, rule4619, rule4620, rule4621, rule4622, rule4623, rule4624, rule4625, rule4626, rule4627, rule4628, rule4629, rule4630, rule4631, rule4632, rule4633, rule4634, rule4635, rule4636, rule4637, rule4638, rule4639, rule4640, rule4641, rule4642, rule4643, rule4644, rule4645, rule4646, rule4647, rule4648, rule4649, rule4650, rule4651, rule4652, rule4653, rule4654, rule4655, rule4656, rule4657, rule4658, rule4659, rule4660, rule4661, rule4662, rule4663, rule4664, rule4665, rule4666, rule4667, rule4668, rule4669, rule4670, rule4671, rule4672, rule4673, rule4674, rule4675, rule4676, rule4677, rule4678, rule4679, rule4680, rule4681, rule4682, rule4683, rule4684, rule4685, rule4686, rule4687, ]





def replacement3920(a, b, e, f, m, n, x):
    return Simp(a*b*(a/sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)


def replacement3921(e, f, m, n, x):
    return Dist(S(1)/f, Subst(Int(x**(-m)*(x**S(2) + S(1))**(m/S(2) + n/S(2) + S(-1)), x), x, tan(e + f*x)), x)


def replacement3922(a, e, f, m, n, x):
    return -Dist(a**(S(1) - n)/f, Subst(Int((a*x)**(m + n + S(-1))*(x**S(2) + S(-1))**(-n/S(2) + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)


def replacement3923(a, e, f, m, n, x):
    return Dist(a**(S(1) - n)/f, Subst(Int((a*x)**(m + n + S(-1))*(x**S(2) + S(-1))**(-n/S(2) + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)


def replacement3924(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(n + S(1))/(b**S(2)*(m + S(-1))), Int((a/sin(e + f*x))**(m + S(-2))*(b/cos(e + f*x))**(n + S(2)), x), x) - Simp(a*(a/sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(1))/(b*f*(m + S(-1))), x)


def replacement3925(a, b, e, f, m, n, x):
    return Dist(b**S(2)*(m + S(1))/(a**S(2)*(n + S(-1))), Int((a/sin(e + f*x))**(m + S(2))*(b/cos(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a/sin(e + f*x))**(m + S(1))*(b/cos(e + f*x))**(n + S(-1))/(a*f*(n + S(-1))), x)


def replacement3926(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + n + S(-2))/(m + S(-1)), Int((a/sin(e + f*x))**(m + S(-2))*(b/cos(e + f*x))**n, x), x) - Simp(a*b*(a/sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(-1))/(f*(m + S(-1))), x)


def replacement3927(a, b, e, f, m, n, x):
    return Dist(b**S(2)*(m + n + S(-2))/(n + S(-1)), Int((a/sin(e + f*x))**m*(b/cos(e + f*x))**(n + S(-2)), x), x) + Simp(a*b*(a/sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)


def replacement3928(a, b, e, f, m, n, x):
    return Dist((m + S(1))/(a**S(2)*(m + n)), Int((a/sin(e + f*x))**(m + S(2))*(b/cos(e + f*x))**n, x), x) + Simp(b*(a/sin(e + f*x))**(m + S(1))*(b/cos(e + f*x))**(n + S(-1))/(a*f*(m + n)), x)


def replacement3929(a, b, e, f, m, n, x):
    return Dist((n + S(1))/(b**S(2)*(m + n)), Int((a/sin(e + f*x))**m*(b/cos(e + f*x))**(n + S(2)), x), x) - Simp(a*(a/sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(1))/(b*f*(m + n)), x)


def replacement3930(a, b, e, f, m, n, x):
    return Dist((a/sin(e + f*x))**m*(b/cos(e + f*x))**n*tan(e + f*x)**(-n), Int(tan(e + f*x)**n, x), x)


def replacement3931(a, b, e, f, m, n, x):
    return Dist((a/sin(e + f*x))**m*(a*sin(e + f*x))**m*(b/cos(e + f*x))**n*(b*cos(e + f*x))**n, Int((a*sin(e + f*x))**(-m)*(b*cos(e + f*x))**(-n), x), x)


def replacement3932(c, d, n, x):
    return Dist(S(1)/d, Subst(Int(ExpandIntegrand((x**S(2) + S(1))**(n/S(2) + S(-1)), x), x), x, tan(c + d*x)), x)


def replacement3933(c, d, n, x):
    return -Dist(S(1)/d, Subst(Int(ExpandIntegrand((x**S(2) + S(1))**(n/S(2) + S(-1)), x), x), x, S(1)/tan(c + d*x)), x)


def replacement3934(b, c, d, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/cos(c + d*x))**(n + S(-2)), x), x) + Simp(b*(b/cos(c + d*x))**(n + S(-1))*sin(c + d*x)/(d*(n + S(-1))), x)


def replacement3935(b, c, d, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/sin(c + d*x))**(n + S(-2)), x), x) - Simp(b*(b/sin(c + d*x))**(n + S(-1))*cos(c + d*x)/(d*(n + S(-1))), x)


def replacement3936(b, c, d, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/cos(c + d*x))**(n + S(2)), x), x) - Simp((b/cos(c + d*x))**(n + S(1))*sin(c + d*x)/(b*d*n), x)


def replacement3937(b, c, d, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/sin(c + d*x))**(n + S(2)), x), x) + Simp((b/sin(c + d*x))**(n + S(1))*cos(c + d*x)/(b*d*n), x)


def replacement3938(c, d, x):
    return Simp(atanh(sin(c + d*x))/d, x)


def replacement3939(c, d, x):
    return -Simp(atanh(cos(c + d*x))/d, x)


def replacement3940(b, c, d, n, x):
    return Dist((b/cos(c + d*x))**n*cos(c + d*x)**n, Int(cos(c + d*x)**(-n), x), x)


def replacement3941(b, c, d, n, x):
    return Dist((b/sin(c + d*x))**n*sin(c + d*x)**n, Int(sin(c + d*x)**(-n), x), x)


def replacement3942(b, c, d, n, x):
    return Simp((cos(c + d*x)/b)**(n + S(-1))*(b/cos(c + d*x))**(n + S(-1))*Int((cos(c + d*x)/b)**(-n), x), x)


def replacement3943(b, c, d, n, x):
    return Simp((sin(c + d*x)/b)**(n + S(-1))*(b/sin(c + d*x))**(n + S(-1))*Int((sin(c + d*x)/b)**(-n), x), x)


def replacement3944(a, b, c, d, x):
    return Dist(b**S(2), Int(cos(c + d*x)**(S(-2)), x), x) + Dist(S(2)*a*b, Int(S(1)/cos(c + d*x), x), x) + Simp(a**S(2)*x, x)


def replacement3945(a, b, c, d, x):
    return Dist(b**S(2), Int(sin(c + d*x)**(S(-2)), x), x) + Dist(S(2)*a*b, Int(S(1)/sin(c + d*x), x), x) + Simp(a**S(2)*x, x)


def replacement3946(a, b, c, d, x):
    return Dist(S(2)*b/d, Subst(Int(S(1)/(a + x**S(2)), x), x, b*tan(c + d*x)/sqrt(a + b/cos(c + d*x))), x)


def replacement3947(a, b, c, d, x):
    return Dist(-S(2)*b/d, Subst(Int(S(1)/(a + x**S(2)), x), x, b/(sqrt(a + b/sin(c + d*x))*tan(c + d*x))), x)


def replacement3948(a, b, c, d, n, x):
    return Dist(a/(n + S(-1)), Int((a + b/cos(c + d*x))**(n + S(-2))*(a*(n + S(-1)) + b*(S(3)*n + S(-4))/cos(c + d*x)), x), x) + Simp(b**S(2)*(a + b/cos(c + d*x))**(n + S(-2))*tan(c + d*x)/(d*(n + S(-1))), x)


def replacement3949(a, b, c, d, n, x):
    return Dist(a/(n + S(-1)), Int((a + b/sin(c + d*x))**(n + S(-2))*(a*(n + S(-1)) + b*(S(3)*n + S(-4))/sin(c + d*x)), x), x) - Simp(b**S(2)*(a + b/sin(c + d*x))**(n + S(-2))/(d*(n + S(-1))*tan(c + d*x)), x)


def replacement3950(a, b, c, d, x):
    return Dist(S(1)/a, Int(sqrt(a + b/cos(c + d*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b/cos(c + d*x))*cos(c + d*x)), x), x)


def replacement3951(a, b, c, d, x):
    return Dist(S(1)/a, Int(sqrt(a + b/sin(c + d*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b/sin(c + d*x))*sin(c + d*x)), x), x)


def replacement3952(a, b, c, d, n, x):
    return Dist(S(1)/(a**S(2)*(S(2)*n + S(1))), Int((a + b/cos(c + d*x))**(n + S(1))*(a*(S(2)*n + S(1)) - b*(n + S(1))/cos(c + d*x)), x), x) + Simp((a + b/cos(c + d*x))**n*tan(c + d*x)/(d*(S(2)*n + S(1))), x)


def replacement3953(a, b, c, d, n, x):
    return Dist(S(1)/(a**S(2)*(S(2)*n + S(1))), Int((a + b/sin(c + d*x))**(n + S(1))*(a*(S(2)*n + S(1)) - b*(n + S(1))/sin(c + d*x)), x), x) - Simp((a + b/sin(c + d*x))**n/(d*(S(2)*n + S(1))*tan(c + d*x)), x)


def replacement3954(a, b, c, d, n, x):
    return -Dist(a**n*tan(c + d*x)/(d*sqrt(S(1) - S(1)/cos(c + d*x))*sqrt(S(1) + S(1)/cos(c + d*x))), Subst(Int((S(1) + b*x/a)**(n + S(-1)/2)/(x*sqrt(S(1) - b*x/a)), x), x, S(1)/cos(c + d*x)), x)


def replacement3955(a, b, c, d, n, x):
    return Dist(a**n/(d*sqrt(S(1) - S(1)/sin(c + d*x))*sqrt(S(1) + S(1)/sin(c + d*x))*tan(c + d*x)), Subst(Int((S(1) + b*x/a)**(n + S(-1)/2)/(x*sqrt(S(1) - b*x/a)), x), x, S(1)/sin(c + d*x)), x)


def replacement3956(a, b, c, d, n, x):
    return Dist(a**IntPart(n)*(S(1) + b/(a*cos(c + d*x)))**(-FracPart(n))*(a + b/cos(c + d*x))**FracPart(n), Int((S(1) + b/(a*cos(c + d*x)))**n, x), x)


def replacement3957(a, b, c, d, n, x):
    return Dist(a**IntPart(n)*(S(1) + b/(a*sin(c + d*x)))**(-FracPart(n))*(a + b/sin(c + d*x))**FracPart(n), Int((S(1) + b/(a*sin(c + d*x)))**n, x), x)


def replacement3958(a, b, c, d, x):
    return Simp(-S(2)*sqrt(b*(S(1) + S(1)/cos(c + d*x))/(a + b/cos(c + d*x)))*sqrt(-b*(S(1) - S(1)/cos(c + d*x))/(a + b/cos(c + d*x)))*(a + b/cos(c + d*x))*EllipticPi(a/(a + b), asin(Rt(a + b, S(2))/sqrt(a + b/cos(c + d*x))), (a - b)/(a + b))/(d*Rt(a + b, S(2))*tan(c + d*x)), x)


def replacement3959(a, b, c, d, x):
    return Simp(S(2)*sqrt(b*(S(1) + S(1)/sin(c + d*x))/(a + b/sin(c + d*x)))*sqrt(-b*(S(1) - S(1)/sin(c + d*x))/(a + b/sin(c + d*x)))*(a + b/sin(c + d*x))*EllipticPi(a/(a + b), asin(Rt(a + b, S(2))/sqrt(a + b/sin(c + d*x))), (a - b)/(a + b))*tan(c + d*x)/(d*Rt(a + b, S(2))), x)


def replacement3960(a, b, c, d, x):
    return Dist(b**S(2), Int((S(1) + S(1)/cos(c + d*x))/(sqrt(a + b/cos(c + d*x))*cos(c + d*x)), x), x) + Int((a**S(2) + b*(S(2)*a - b)/cos(c + d*x))/sqrt(a + b/cos(c + d*x)), x)


def replacement3961(a, b, c, d, x):
    return Dist(b**S(2), Int((S(1) + S(1)/sin(c + d*x))/(sqrt(a + b/sin(c + d*x))*sin(c + d*x)), x), x) + Int((a**S(2) + b*(S(2)*a - b)/sin(c + d*x))/sqrt(a + b/sin(c + d*x)), x)


def replacement3962(a, b, c, d, n, x):
    return Dist(S(1)/(n + S(-1)), Int((a + b/cos(c + d*x))**(n + S(-3))*Simp(a**S(3)*(n + S(-1)) + a*b**S(2)*(S(3)*n + S(-4))/cos(c + d*x)**S(2) + b*(S(3)*a**S(2)*(n + S(-1)) + b**S(2)*(n + S(-2)))/cos(c + d*x), x), x), x) + Simp(b**S(2)*(a + b/cos(c + d*x))**(n + S(-2))*tan(c + d*x)/(d*(n + S(-1))), x)


def replacement3963(a, b, c, d, n, x):
    return Dist(S(1)/(n + S(-1)), Int((a + b/sin(c + d*x))**(n + S(-3))*Simp(a**S(3)*(n + S(-1)) + a*b**S(2)*(S(3)*n + S(-4))/sin(c + d*x)**S(2) + b*(S(3)*a**S(2)*(n + S(-1)) + b**S(2)*(n + S(-2)))/sin(c + d*x), x), x), x) - Simp(b**S(2)*(a + b/sin(c + d*x))**(n + S(-2))/(d*(n + S(-1))*tan(c + d*x)), x)


def replacement3964(a, b, c, d, x):
    return -Dist(S(1)/a, Int(S(1)/(a*cos(c + d*x)/b + S(1)), x), x) + Simp(x/a, x)


def replacement3965(a, b, c, d, x):
    return -Dist(S(1)/a, Int(S(1)/(a*sin(c + d*x)/b + S(1)), x), x) + Simp(x/a, x)


def replacement3966(a, b, c, d, x):
    return Simp(-S(2)*sqrt(b*(S(1) - S(1)/cos(c + d*x))/(a + b))*sqrt(-b*(S(1) + S(1)/cos(c + d*x))/(a - b))*EllipticPi((a + b)/a, asin(sqrt(a + b/cos(c + d*x))/Rt(a + b, S(2))), (a + b)/(a - b))*Rt(a + b, S(2))/(a*d*tan(c + d*x)), x)


def replacement3967(a, b, c, d, x):
    return Simp(S(2)*sqrt(b*(S(1) - S(1)/sin(c + d*x))/(a + b))*sqrt(-b*(S(1) + S(1)/sin(c + d*x))/(a - b))*EllipticPi((a + b)/a, asin(sqrt(a + b/sin(c + d*x))/Rt(a + b, S(2))), (a + b)/(a - b))*Rt(a + b, S(2))*tan(c + d*x)/(a*d), x)


def replacement3968(a, b, c, d, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(n + S(1))), Int((a + b/cos(c + d*x))**(n + S(1))*Simp(-a*b*(n + S(1))/cos(c + d*x) + b**S(2)*(n + S(2))/cos(c + d*x)**S(2) + (a**S(2) - b**S(2))*(n + S(1)), x), x), x) - Simp(b**S(2)*(a + b/cos(c + d*x))**(n + S(1))*tan(c + d*x)/(a*d*(a**S(2) - b**S(2))*(n + S(1))), x)


def replacement3969(a, b, c, d, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(n + S(1))), Int((a + b/sin(c + d*x))**(n + S(1))*Simp(-a*b*(n + S(1))/sin(c + d*x) + b**S(2)*(n + S(2))/sin(c + d*x)**S(2) + (a**S(2) - b**S(2))*(n + S(1)), x), x), x) + Simp(b**S(2)*(a + b/sin(c + d*x))**(n + S(1))/(a*d*(a**S(2) - b**S(2))*(n + S(1))*tan(c + d*x)), x)


def replacement3970(a, b, c, d, n, x):
    return Int((a + b/cos(c + d*x))**n, x)


def replacement3971(a, b, c, d, n, x):
    return Int((a + b/sin(c + d*x))**n, x)


def replacement3972(a, b, d, e, f, n, x):
    return Dist(a, Int((d/cos(e + f*x))**n, x), x) + Dist(b/d, Int((d/cos(e + f*x))**(n + S(1)), x), x)


def replacement3973(a, b, d, e, f, n, x):
    return Dist(a, Int((d/sin(e + f*x))**n, x), x) + Dist(b/d, Int((d/sin(e + f*x))**(n + S(1)), x), x)


def replacement3974(a, b, d, e, f, n, x):
    return Dist(S(2)*a*b/d, Int((d/cos(e + f*x))**(n + S(1)), x), x) + Int((d/cos(e + f*x))**n*(a**S(2) + b**S(2)/cos(e + f*x)**S(2)), x)


def replacement3975(a, b, d, e, f, n, x):
    return Dist(S(2)*a*b/d, Int((d/sin(e + f*x))**(n + S(1)), x), x) + Int((d/sin(e + f*x))**n*(a**S(2) + b**S(2)/sin(e + f*x)**S(2)), x)


def replacement3976(a, b, e, f, x):
    return Dist(S(1)/b, Int(S(1)/cos(e + f*x), x), x) - Dist(a/b, Int(S(1)/((a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement3977(a, b, e, f, x):
    return Dist(S(1)/b, Int(S(1)/sin(e + f*x), x), x) - Dist(a/b, Int(S(1)/((a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement3978(a, b, e, f, x):
    return -Dist(a/b, Int(S(1)/((a + b/cos(e + f*x))*cos(e + f*x)**S(2)), x), x) + Simp(tan(e + f*x)/(b*f), x)


def replacement3979(a, b, e, f, x):
    return -Dist(a/b, Int(S(1)/((a + b/sin(e + f*x))*sin(e + f*x)**S(2)), x), x) - Simp(S(1)/(b*f*tan(e + f*x)), x)


def replacement3980(a, b, d, e, f, m, n, x):
    return Int(ExpandTrig((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m, x), x)


def replacement3981(a, b, d, e, f, m, n, x):
    return Int(ExpandTrig((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m, x), x)


def replacement3982(a, b, e, f, x):
    return Simp(S(2)*b*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))), x)


def replacement3983(a, b, e, f, x):
    return Simp(-S(2)*b/(f*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement3984(a, b, e, f, m, x):
    return Dist(a*(S(2)*m + S(-1))/m, Int((a + b/cos(e + f*x))**(m + S(-1))/cos(e + f*x), x), x) + Simp(b*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*m), x)


def replacement3985(a, b, e, f, m, x):
    return Dist(a*(S(2)*m + S(-1))/m, Int((a + b/sin(e + f*x))**(m + S(-1))/sin(e + f*x), x), x) - Simp(b*(a + b/sin(e + f*x))**(m + S(-1))/(f*m*tan(e + f*x)), x)


def replacement3986(a, b, e, f, x):
    return Simp(tan(e + f*x)/(f*(a/cos(e + f*x) + b)), x)


def replacement3987(a, b, e, f, x):
    return -Simp(S(1)/(f*(a/sin(e + f*x) + b)*tan(e + f*x)), x)


def replacement3988(a, b, e, f, x):
    return Dist(S(2)/f, Subst(Int(S(1)/(S(2)*a + x**S(2)), x), x, b*tan(e + f*x)/sqrt(a + b/cos(e + f*x))), x)


def replacement3989(a, b, e, f, x):
    return Dist(-S(2)/f, Subst(Int(S(1)/(S(2)*a + x**S(2)), x), x, b/(sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement3990(a, b, e, f, m, x):
    return Dist((m + S(1))/(a*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x) - Simp(b*(a + b/cos(e + f*x))**m*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement3991(a, b, e, f, m, x):
    return Dist((m + S(1))/(a*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x) + Simp(b*(a + b/sin(e + f*x))**m/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement3992(a, b, e, f, m, x):
    return Dist(m/(b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement3993(a, b, e, f, m, x):
    return Dist(m/(b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement3994(a, b, e, f, m, x):
    return Dist(a*m/(b*(m + S(1))), Int((a + b/cos(e + f*x))**m/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement3995(a, b, e, f, m, x):
    return Dist(a*m/(b*(m + S(1))), Int((a + b/sin(e + f*x))**m/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement3996(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(1))/cos(e + f*x))/cos(e + f*x), x), x) - Simp(b*(a + b/cos(e + f*x))**m*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement3997(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(1))/sin(e + f*x))/sin(e + f*x), x), x) + Simp(b*(a + b/sin(e + f*x))**m/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement3998(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/cos(e + f*x))**m*(-a/cos(e + f*x) + b*(m + S(1)))/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(2))), x)


def replacement3999(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/sin(e + f*x))**m*(-a/sin(e + f*x) + b*(m + S(1)))/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(2))*tan(e + f*x)), x)


def replacement4000(a, b, d, e, f, x):
    return Dist(S(2)*a*sqrt(a*d/b)/(b*f), Subst(Int(S(1)/sqrt(S(1) + x**S(2)/a), x), x, b*tan(e + f*x)/sqrt(a + b/cos(e + f*x))), x)


def replacement4001(a, b, d, e, f, x):
    return Dist(-S(2)*a*sqrt(a*d/b)/(b*f), Subst(Int(S(1)/sqrt(S(1) + x**S(2)/a), x), x, b/(sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4002(a, b, d, e, f, x):
    return Dist(S(2)*b*d/f, Subst(Int(S(1)/(b - d*x**S(2)), x), x, b*tan(e + f*x)/(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x)))), x)


def replacement4003(a, b, d, e, f, x):
    return Dist(-S(2)*b*d/f, Subst(Int(S(1)/(b - d*x**S(2)), x), x, b/(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4004(a, b, d, e, f, n, x):
    return Dist(S(2)*a*d*(n + S(-1))/(b*(S(2)*n + S(-1))), Int((d/cos(e + f*x))**(n + S(-1))*sqrt(a + b/cos(e + f*x)), x), x) + Simp(S(2)*b*d*(d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(-1))), x)


def replacement4005(a, b, d, e, f, n, x):
    return Dist(S(2)*a*d*(n + S(-1))/(b*(S(2)*n + S(-1))), Int((d/sin(e + f*x))**(n + S(-1))*sqrt(a + b/sin(e + f*x)), x), x) + Simp(-S(2)*b*d*(d/sin(e + f*x))**(n + S(-1))/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(-1))*tan(e + f*x)), x)


def replacement4006(a, b, d, e, f, x):
    return Simp(S(2)*a*tan(e + f*x)/(f*sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), x)


def replacement4007(a, b, d, e, f, x):
    return Simp(-S(2)*a/(f*sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4008(a, b, d, e, f, n, x):
    return Dist(a*(S(2)*n + S(1))/(S(2)*b*d*n), Int((d/cos(e + f*x))**(n + S(1))*sqrt(a + b/cos(e + f*x)), x), x) - Simp(a*(d/cos(e + f*x))**n*tan(e + f*x)/(f*n*sqrt(a + b/cos(e + f*x))), x)


def replacement4009(a, b, d, e, f, n, x):
    return Dist(a*(S(2)*n + S(1))/(S(2)*b*d*n), Int((d/sin(e + f*x))**(n + S(1))*sqrt(a + b/sin(e + f*x)), x), x) + Simp(a*(d/sin(e + f*x))**n/(f*n*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4010(a, b, d, e, f, n, x):
    return -Dist(a**S(2)*d*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int((d*x)**(n + S(-1))/sqrt(a - b*x), x), x, S(1)/cos(e + f*x)), x)


def replacement4011(a, b, d, e, f, n, x):
    return Dist(a**S(2)*d/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), Subst(Int((d*x)**(n + S(-1))/sqrt(a - b*x), x), x, S(1)/sin(e + f*x)), x)


def replacement4012(a, b, d, e, f, x):
    return Dist(sqrt(S(2))*sqrt(a)/(b*f), Subst(Int(S(1)/sqrt(x**S(2) + S(1)), x), x, b*tan(e + f*x)/(a + b/cos(e + f*x))), x)


def replacement4013(a, b, d, e, f, x):
    return -Dist(sqrt(S(2))*sqrt(a)/(b*f), Subst(Int(S(1)/sqrt(x**S(2) + S(1)), x), x, b/((a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4014(a, b, d, e, f, x):
    return Dist(S(2)*b*d/(a*f), Subst(Int(S(1)/(S(2)*b - d*x**S(2)), x), x, b*tan(e + f*x)/(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x)))), x)


def replacement4015(a, b, d, e, f, x):
    return Dist(-S(2)*b*d/(a*f), Subst(Int(S(1)/(S(2)*b - d*x**S(2)), x), x, b/(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4016(a, b, d, e, f, m, n, x):
    return Dist(b*(S(2)*m + S(-1))/(d*m), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1)), x), x) + Simp(a*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*m), x)


def replacement4017(a, b, d, e, f, m, n, x):
    return Dist(b*(S(2)*m + S(-1))/(d*m), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1)), x), x) - Simp(a*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))/(f*m*tan(e + f*x)), x)


def replacement4018(a, b, d, e, f, m, n, x):
    return Dist(d*(m + S(1))/(b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1)), x), x) - Simp(b*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4019(a, b, d, e, f, m, n, x):
    return Dist(d*(m + S(1))/(b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1)), x), x) + Simp(b*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4020(a, b, d, e, f, m, n, x):
    return Dist(m/(a*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1)), x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement4021(a, b, d, e, f, m, n, x):
    return Dist(m/(a*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1)), x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4022(a, b, d, e, f, m, n, x):
    return Dist(a*m/(b*d*(m + S(1))), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m, x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4023(a, b, d, e, f, m, n, x):
    return Dist(a*m/(b*d*(m + S(1))), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m, x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4024(a, b, d, e, f, m, n, x):
    return -Dist(a/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-2))*(-a*(m + S(2)*n + S(-1))/cos(e + f*x) + b*(m - S(2)*n + S(-2))), x), x) - Simp(b**S(2)*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*tan(e + f*x)/(f*n), x)


def replacement4025(a, b, d, e, f, m, n, x):
    return -Dist(a/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-2))*(-a*(m + S(2)*n + S(-1))/sin(e + f*x) + b*(m - S(2)*n + S(-2))), x), x) + Simp(b**S(2)*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))/(f*n*tan(e + f*x)), x)


def replacement4026(a, b, d, e, f, m, n, x):
    return Dist(b/(m + n + S(-1)), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*(a*(S(3)*m + S(2)*n + S(-4))/cos(e + f*x) + b*(m + S(2)*n + S(-1))), x), x) + Simp(b**S(2)*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*tan(e + f*x)/(f*(m + n + S(-1))), x)


def replacement4027(a, b, d, e, f, m, n, x):
    return Dist(b/(m + n + S(-1)), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))*(a*(S(3)*m + S(2)*n + S(-4))/sin(e + f*x) + b*(m + S(2)*n + S(-1))), x), x) - Simp(b**S(2)*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))/(f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4028(a, b, d, e, f, m, n, x):
    return -Dist(d/(a*b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*(a*(n + S(-1)) - b*(m + n)/cos(e + f*x)), x), x) - Simp(b*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4029(a, b, d, e, f, m, n, x):
    return -Dist(d/(a*b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*(a*(n + S(-1)) - b*(m + n)/sin(e + f*x)), x), x) + Simp(b*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4030(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(a*b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*(a*(m - n + S(2))/cos(e + f*x) + b*(n + S(-2))), x), x) + Simp(d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement4031(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(a*b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))*(a*(m - n + S(2))/sin(e + f*x) + b*(n + S(-2))), x), x) - Simp(d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4032(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*(a*(S(2)*m + n + S(1)) - b*(m + n + S(1))/cos(e + f*x)), x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement4033(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*(a*(S(2)*m + n + S(1)) - b*(m + n + S(1))/sin(e + f*x)), x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4034(a, b, d, e, f, n, x):
    return -Dist(d**S(2)/(a*b), Int((d/cos(e + f*x))**(n + S(-2))*(-a*(n + S(-1))/cos(e + f*x) + b*(n + S(-2))), x), x) - Simp(d**S(2)*(d/cos(e + f*x))**(n + S(-2))*tan(e + f*x)/(f*(a + b/cos(e + f*x))), x)


def replacement4035(a, b, d, e, f, n, x):
    return -Dist(d**S(2)/(a*b), Int((d/sin(e + f*x))**(n + S(-2))*(-a*(n + S(-1))/sin(e + f*x) + b*(n + S(-2))), x), x) + Simp(d**S(2)*(d/sin(e + f*x))**(n + S(-2))/(f*(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4036(a, b, d, e, f, n, x):
    return -Dist(a**(S(-2)), Int((d/cos(e + f*x))**n*(a*(n + S(-1)) - b*n/cos(e + f*x)), x), x) - Simp((d/cos(e + f*x))**n*tan(e + f*x)/(f*(a + b/cos(e + f*x))), x)


def replacement4037(a, b, d, e, f, n, x):
    return -Dist(a**(S(-2)), Int((d/sin(e + f*x))**n*(a*(n + S(-1)) - b*n/sin(e + f*x)), x), x) + Simp((d/sin(e + f*x))**n/(f*(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4038(a, b, d, e, f, n, x):
    return Dist(d*(n + S(-1))/(a*b), Int((d/cos(e + f*x))**(n + S(-1))*(a - b/cos(e + f*x)), x), x) + Simp(b*d*(d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(a*f*(a + b/cos(e + f*x))), x)


def replacement4039(a, b, d, e, f, n, x):
    return Dist(d*(n + S(-1))/(a*b), Int((d/sin(e + f*x))**(n + S(-1))*(a - b/sin(e + f*x)), x), x) - Simp(b*d*(d/sin(e + f*x))**(n + S(-1))/(a*f*(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4040(a, b, d, e, f, x):
    return Dist(d/b, Int(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x)), x), x) - Dist(a*d/b, Int(sqrt(d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4041(a, b, d, e, f, x):
    return Dist(d/b, Int(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x)), x), x) - Dist(a*d/b, Int(sqrt(d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4042(a, b, d, e, f, n, x):
    return Dist(d**S(2)/(b*(S(2)*n + S(-3))), Int((d/cos(e + f*x))**(n + S(-2))*(-a/cos(e + f*x) + S(2)*b*(n + S(-2)))/sqrt(a + b/cos(e + f*x)), x), x) + Simp(S(2)*d**S(2)*(d/cos(e + f*x))**(n + S(-2))*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(-3))), x)


def replacement4043(a, b, d, e, f, n, x):
    return Dist(d**S(2)/(b*(S(2)*n + S(-3))), Int((d/sin(e + f*x))**(n + S(-2))*(-a/sin(e + f*x) + S(2)*b*(n + S(-2)))/sqrt(a + b/sin(e + f*x)), x), x) + Simp(-S(2)*d**S(2)*(d/sin(e + f*x))**(n + S(-2))/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(-3))*tan(e + f*x)), x)


def replacement4044(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*b*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b*(S(2)*n + S(1))/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) - Simp((d/cos(e + f*x))**n*tan(e + f*x)/(f*n*sqrt(a + b/cos(e + f*x))), x)


def replacement4045(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*b*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b*(S(2)*n + S(1))/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Simp((d/sin(e + f*x))**n/(f*n*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4046(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n + S(-1))), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**m*(a*m/cos(e + f*x) + b*(n + S(-2))), x), x) + Simp(d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(-1))), x)


def replacement4047(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n + S(-1))), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**m*(a*m/sin(e + f*x) + b*(n + S(-2))), x), x) - Simp(d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**m/(f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4048(a, b, d, e, f, m, n, x):
    return Dist(a**(S(2) - n)*(a*d/b)**n*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int((a - x)**(n + S(-1))*(S(2)*a - x)**(m + S(-1)/2)/sqrt(x), x), x, a - b/cos(e + f*x)), x)


def replacement4049(a, b, d, e, f, m, n, x):
    return -Dist(a**(S(2) - n)*(a*d/b)**n/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), Subst(Int((a - x)**(n + S(-1))*(S(2)*a - x)**(m + S(-1)/2)/sqrt(x), x), x, a - b/sin(e + f*x)), x)


def replacement4050(a, b, d, e, f, m, n, x):
    return Dist(a**(S(1) - n)*(-a*d/b)**n*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int(x**(m + S(-1)/2)*(a - x)**(n + S(-1))/sqrt(S(2)*a - x), x), x, a + b/cos(e + f*x)), x)


def replacement4051(a, b, d, e, f, m, n, x):
    return -Dist(a**(S(1) - n)*(-a*d/b)**n/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), Subst(Int(x**(m + S(-1)/2)*(a - x)**(n + S(-1))/sqrt(S(2)*a - x), x), x, a + b/sin(e + f*x)), x)


def replacement4052(a, b, d, e, f, m, n, x):
    return -Dist(a**S(2)*d*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int((d*x)**(n + S(-1))*(a + b*x)**(m + S(-1)/2)/sqrt(a - b*x), x), x, S(1)/cos(e + f*x)), x)


def replacement4053(a, b, d, e, f, m, n, x):
    return Dist(a**S(2)*d/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), Subst(Int((d*x)**(n + S(-1))*(a + b*x)**(m + S(-1)/2)/sqrt(a - b*x), x), x, S(1)/sin(e + f*x)), x)


def replacement4054(a, b, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*(S(1) + b/(a*cos(e + f*x)))**(-FracPart(m))*(a + b/cos(e + f*x))**FracPart(m), Int((d/cos(e + f*x))**n*(S(1) + b/(a*cos(e + f*x)))**m, x), x)


def replacement4055(a, b, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*(S(1) + b/(a*sin(e + f*x)))**(-FracPart(m))*(a + b/sin(e + f*x))**FracPart(m), Int((d/sin(e + f*x))**n*(S(1) + b/(a*sin(e + f*x)))**m, x), x)


def replacement4056(a, b, e, f, x):
    return Dist(b, Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Dist(a - b, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4057(a, b, e, f, x):
    return Dist(b, Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Dist(a - b, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4058(a, b, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/cos(e + f*x))**(m + S(-2))*(a**S(2)*m + a*b*(S(2)*m + S(-1))/cos(e + f*x) + b**S(2)*(m + S(-1)))/cos(e + f*x), x), x) + Simp(b*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*m), x)


def replacement4059(a, b, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/sin(e + f*x))**(m + S(-2))*(a**S(2)*m + a*b*(S(2)*m + S(-1))/sin(e + f*x) + b**S(2)*(m + S(-1)))/sin(e + f*x), x), x) - Simp(b*(a + b/sin(e + f*x))**(m + S(-1))/(f*m*tan(e + f*x)), x)


def replacement4060(a, b, e, f, x):
    return Dist(S(1)/b, Int(S(1)/(a*cos(e + f*x)/b + S(1)), x), x)


def replacement4061(a, b, e, f, x):
    return Dist(S(1)/b, Int(S(1)/(a*sin(e + f*x)/b + S(1)), x), x)


def replacement4062(a, b, e, f, x):
    return Simp(S(2)*sqrt(b*(S(1) - S(1)/cos(e + f*x))/(a + b))*sqrt(-b*(S(1) + S(1)/cos(e + f*x))/(a - b))*EllipticF(asin(sqrt(a + b/cos(e + f*x))/Rt(a + b, S(2))), (a + b)/(a - b))*Rt(a + b, S(2))/(b*f*tan(e + f*x)), x)


def replacement4063(a, b, e, f, x):
    return Simp(-S(2)*sqrt(b*(S(1) - S(1)/sin(e + f*x))/(a + b))*sqrt(-b*(S(1) + S(1)/sin(e + f*x))/(a - b))*EllipticF(asin(sqrt(a + b/sin(e + f*x))/Rt(a + b, S(2))), (a + b)/(a - b))*Rt(a + b, S(2))*tan(e + f*x)/(b*f), x)


def replacement4064(a, b, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(a*(m + S(1)) - b*(m + S(2))/cos(e + f*x))/cos(e + f*x), x), x) + Simp(b*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4065(a, b, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(a*(m + S(1)) - b*(m + S(2))/sin(e + f*x))/sin(e + f*x), x), x) - Simp(b*(a + b/sin(e + f*x))**(m + S(1))/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4066(a, b, e, f, m, x):
    return -Dist(tan(e + f*x)/(f*sqrt(S(1) - S(1)/cos(e + f*x))*sqrt(S(1) + S(1)/cos(e + f*x))), Subst(Int((a + b*x)**m/(sqrt(S(1) - x)*sqrt(x + S(1))), x), x, S(1)/cos(e + f*x)), x)


def replacement4067(a, b, e, f, m, x):
    return Dist(S(1)/(f*sqrt(S(1) - S(1)/sin(e + f*x))*sqrt(S(1) + S(1)/sin(e + f*x))*tan(e + f*x)), Subst(Int((a + b*x)**m/(sqrt(S(1) - x)*sqrt(x + S(1))), x), x, S(1)/sin(e + f*x)), x)


def replacement4068(a, b, e, f, m, x):
    return Dist(m/(m + S(1)), Int((a + b/cos(e + f*x))**(m + S(-1))*(a/cos(e + f*x) + b)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4069(a, b, e, f, m, x):
    return Dist(m/(m + S(1)), Int((a + b/sin(e + f*x))**(m + S(-1))*(a/sin(e + f*x) + b)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4070(a, b, e, f, m, x):
    return -Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(-a*(m + S(2))/cos(e + f*x) + b*(m + S(1)))/cos(e + f*x), x), x) - Simp(a*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4071(a, b, e, f, m, x):
    return -Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(-a*(m + S(2))/sin(e + f*x) + b*(m + S(1)))/sin(e + f*x), x), x) + Simp(a*(a + b/sin(e + f*x))**(m + S(1))/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4072(a, b, e, f, x):
    return -Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x) + Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x)


def replacement4073(a, b, e, f, x):
    return -Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x) + Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x)


def replacement4074(a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x) - Dist(a/b, Int((a + b/cos(e + f*x))**m/cos(e + f*x), x), x)


def replacement4075(a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x) - Dist(a/b, Int((a + b/sin(e + f*x))**m/sin(e + f*x), x), x)


def replacement4076(a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(a*b*(m + S(1)) - (a**S(2) + b**S(2)*(m + S(1)))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp(a**S(2)*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4077(a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(a*b*(m + S(1)) - (a**S(2) + b**S(2)*(m + S(1)))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp(a**S(2)*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4078(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/cos(e + f*x))**m*(-a/cos(e + f*x) + b*(m + S(1)))/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(2))), x)


def replacement4079(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/sin(e + f*x))**m*(-a/sin(e + f*x) + b*(m + S(1)))/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(2))*tan(e + f*x)), x)


def replacement4080(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-3))*Simp(a**S(2)*b*(m - S(2)*n + S(-2)) - a*(a**S(2)*(n + S(1)) + S(3)*b**S(2)*n)/cos(e + f*x) - b*(a**S(2)*(m + n + S(-1)) + b**S(2)*n)/cos(e + f*x)**S(2), x), x), x) - Simp(a**S(2)*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*tan(e + f*x)/(f*n), x)


def replacement4081(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-3))*Simp(a**S(2)*b*(m - S(2)*n + S(-2)) - a*(a**S(2)*(n + S(1)) + S(3)*b**S(2)*n)/sin(e + f*x) - b*(a**S(2)*(m + n + S(-1)) + b**S(2)*n)/sin(e + f*x)**S(2), x), x), x) + Simp(a**S(2)*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))/(f*n*tan(e + f*x)), x)


def replacement4082(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(-1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-3))*Simp(a**S(3)*d*(m + n + S(-1)) + a*b**S(2)*d*n + a*b**S(2)*d*(S(3)*m + S(2)*n + S(-4))/cos(e + f*x)**S(2) + b*(S(3)*a**S(2)*d*(m + n + S(-1)) + b**S(2)*d*(m + n + S(-2)))/cos(e + f*x), x), x), x) + Simp(b**S(2)*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*tan(e + f*x)/(f*(m + n + S(-1))), x)


def replacement4083(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(-1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-3))*Simp(a**S(3)*d*(m + n + S(-1)) + a*b**S(2)*d*n + a*b**S(2)*d*(S(3)*m + S(2)*n + S(-4))/sin(e + f*x)**S(2) + b*(S(3)*a**S(2)*d*(m + n + S(-1)) + b**S(2)*d*(m + n + S(-2)))/sin(e + f*x), x), x), x) - Simp(b**S(2)*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))/(f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4084(a, b, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*Simp(a*d*(m + S(1))/cos(e + f*x) + b*d*(n + S(-1)) - b*d*(m + n + S(1))/cos(e + f*x)**S(2), x), x), x) + Simp(b*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4085(a, b, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*Simp(a*d*(m + S(1))/sin(e + f*x) + b*d*(n + S(-1)) - b*d*(m + n + S(1))/sin(e + f*x)**S(2), x), x), x) - Simp(b*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4086(a, b, d, e, f, m, n, x):
    return -Dist(d**S(2)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*(-a*(m + n)/cos(e + f*x)**S(2) + a*(n + S(-2)) + b*(m + S(1))/cos(e + f*x)), x), x) - Simp(a*d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4087(a, b, d, e, f, m, n, x):
    return -Dist(d**S(2)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))*(-a*(m + n)/sin(e + f*x)**S(2) + a*(n + S(-2)) + b*(m + S(1))/sin(e + f*x)), x), x) + Simp(a*d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4088(a, b, d, e, f, m, n, x):
    return Dist(d**S(3)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-3))*(a + b/cos(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(-3)) + a*b*(m + S(1))/cos(e + f*x) - (a**S(2)*(n + S(-2)) + b**S(2)*(m + S(1)))/cos(e + f*x)**S(2), x), x), x) + Simp(a**S(2)*d**S(3)*(d/cos(e + f*x))**(n + S(-3))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4089(a, b, d, e, f, m, n, x):
    return Dist(d**S(3)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-3))*(a + b/sin(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(-3)) + a*b*(m + S(1))/sin(e + f*x) - (a**S(2)*(n + S(-2)) + b**S(2)*(m + S(1)))/sin(e + f*x)**S(2), x), x), x) - Simp(a**S(2)*d**S(3)*(d/sin(e + f*x))**(n + S(-3))*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4090(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(-a*(n + S(1))/cos(e + f*x) + b*(m + n + S(1)) - b*(m + n + S(2))/cos(e + f*x)**S(2), x), x), x) - Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f*n), x)


def replacement4091(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(-a*(n + S(1))/sin(e + f*x) + b*(m + n + S(1)) - b*(m + n + S(2))/sin(e + f*x)**S(2), x), x), x) + Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))/(a*f*n*tan(e + f*x)), x)


def replacement4092(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*(a**S(2)*(m + S(1)) - a*b*(m + S(1))/cos(e + f*x) - b**S(2)*(m + n + S(1)) + b**S(2)*(m + n + S(2))/cos(e + f*x)**S(2)), x), x) - Simp(b**S(2)*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4093(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*(a**S(2)*(m + S(1)) - a*b*(m + S(1))/sin(e + f*x) - b**S(2)*(m + n + S(1)) + b**S(2)*(m + n + S(2))/sin(e + f*x)**S(2)), x), x) + Simp(b**S(2)*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4094(a, b, d, e, f, x):
    return Dist(sqrt(d/cos(e + f*x))*sqrt(d*cos(e + f*x))/d, Int(sqrt(d*cos(e + f*x))/(a*cos(e + f*x) + b), x), x)


def replacement4095(a, b, d, e, f, x):
    return Dist(sqrt(d/sin(e + f*x))*sqrt(d*sin(e + f*x))/d, Int(sqrt(d*sin(e + f*x))/(a*sin(e + f*x) + b), x), x)


def replacement4096(a, b, d, e, f, x):
    return Dist(d*sqrt(d/cos(e + f*x))*sqrt(d*cos(e + f*x)), Int(S(1)/(sqrt(d*cos(e + f*x))*(a*cos(e + f*x) + b)), x), x)


def replacement4097(a, b, d, e, f, x):
    return Dist(d*sqrt(d/sin(e + f*x))*sqrt(d*sin(e + f*x)), Int(S(1)/(sqrt(d*sin(e + f*x))*(a*sin(e + f*x) + b)), x), x)


def replacement4098(a, b, d, e, f, x):
    return Dist(d/b, Int((d/cos(e + f*x))**(S(3)/2), x), x) - Dist(a*d/b, Int((d/cos(e + f*x))**(S(3)/2)/(a + b/cos(e + f*x)), x), x)


def replacement4099(a, b, d, e, f, x):
    return Dist(d/b, Int((d/sin(e + f*x))**(S(3)/2), x), x) - Dist(a*d/b, Int((d/sin(e + f*x))**(S(3)/2)/(a + b/sin(e + f*x)), x), x)


def replacement4100(a, b, d, e, f, n, x):
    return Dist(d**S(3)/(b*(n + S(-2))), Int((d/cos(e + f*x))**(n + S(-3))*Simp(a*(n + S(-3)) - a*(n + S(-2))/cos(e + f*x)**S(2) + b*(n + S(-3))/cos(e + f*x), x)/(a + b/cos(e + f*x)), x), x) + Simp(d**S(3)*(d/cos(e + f*x))**(n + S(-3))*tan(e + f*x)/(b*f*(n + S(-2))), x)


def replacement4101(a, b, d, e, f, n, x):
    return Dist(d**S(3)/(b*(n + S(-2))), Int((d/sin(e + f*x))**(n + S(-3))*Simp(a*(n + S(-3)) - a*(n + S(-2))/sin(e + f*x)**S(2) + b*(n + S(-3))/sin(e + f*x), x)/(a + b/sin(e + f*x)), x), x) - Simp(d**S(3)*(d/sin(e + f*x))**(n + S(-3))/(b*f*(n + S(-2))*tan(e + f*x)), x)


def replacement4102(a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((a - b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) + Dist(b**S(2)/(a**S(2)*d**S(2)), Int((d/cos(e + f*x))**(S(3)/2)/(a + b/cos(e + f*x)), x), x)


def replacement4103(a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((a - b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) + Dist(b**S(2)/(a**S(2)*d**S(2)), Int((d/sin(e + f*x))**(S(3)/2)/(a + b/sin(e + f*x)), x), x)


def replacement4104(a, b, d, e, f, n, x):
    return -Dist(S(1)/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(-a*(n + S(1))/cos(e + f*x) + b*n - b*(n + S(1))/cos(e + f*x)**S(2), x)/(a + b/cos(e + f*x)), x), x) - Simp((d/cos(e + f*x))**n*tan(e + f*x)/(a*f*n), x)


def replacement4105(a, b, d, e, f, n, x):
    return -Dist(S(1)/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(-a*(n + S(1))/sin(e + f*x) + b*n - b*(n + S(1))/sin(e + f*x)**S(2), x)/(a + b/sin(e + f*x)), x), x) + Simp((d/sin(e + f*x))**n/(a*f*n*tan(e + f*x)), x)


def replacement4106(a, b, d, e, f, x):
    return Dist(a, Int(sqrt(d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist(b/d, Int((d/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4107(a, b, d, e, f, x):
    return Dist(a, Int(sqrt(d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist(b/d, Int((d/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4108(a, b, d, e, f, n, x):
    return Dist(d**S(2)/(S(2)*n + S(-1)), Int((d/cos(e + f*x))**(n + S(-2))*Simp(S(2)*a*(n + S(-2)) + a/cos(e + f*x)**S(2) + b*(S(2)*n + S(-3))/cos(e + f*x), x)/sqrt(a + b/cos(e + f*x)), x), x) + Simp(S(2)*d*(d/cos(e + f*x))**(n + S(-1))*sqrt(a + b/cos(e + f*x))*sin(e + f*x)/(f*(S(2)*n + S(-1))), x)


def replacement4109(a, b, d, e, f, n, x):
    return Dist(d**S(2)/(S(2)*n + S(-1)), Int((d/sin(e + f*x))**(n + S(-2))*Simp(S(2)*a*(n + S(-2)) + a/sin(e + f*x)**S(2) + b*(S(2)*n + S(-3))/sin(e + f*x), x)/sqrt(a + b/sin(e + f*x)), x), x) + Simp(-S(2)*d*(d/sin(e + f*x))**(n + S(-1))*sqrt(a + b/sin(e + f*x))*cos(e + f*x)/(f*(S(2)*n + S(-1))), x)


def replacement4110(a, b, d, e, f, x):
    return Dist(sqrt(a + b/cos(e + f*x))/(sqrt(d/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)), Int(sqrt(a*cos(e + f*x) + b), x), x)


def replacement4111(a, b, d, e, f, x):
    return Dist(sqrt(a + b/sin(e + f*x))/(sqrt(d/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)), Int(sqrt(a*sin(e + f*x) + b), x), x)


def replacement4112(a, b, d, e, f, n, x):
    return -Dist(S(1)/(S(2)*d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(-S(2)*a*(n + S(1))/cos(e + f*x) - b*(S(2)*n + S(3))/cos(e + f*x)**S(2) + b, x)/sqrt(a + b/cos(e + f*x)), x), x) - Simp((d/cos(e + f*x))**n*sqrt(a + b/cos(e + f*x))*tan(e + f*x)/(f*n), x)


def replacement4113(a, b, d, e, f, n, x):
    return -Dist(S(1)/(S(2)*d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(-S(2)*a*(n + S(1))/sin(e + f*x) - b*(S(2)*n + S(3))/sin(e + f*x)**S(2) + b, x)/sqrt(a + b/sin(e + f*x)), x), x) + Simp((d/sin(e + f*x))**n*sqrt(a + b/sin(e + f*x))/(f*n*tan(e + f*x)), x)


def replacement4114(a, b, d, e, f, x):
    return Dist(sqrt(d/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)/sqrt(a + b/cos(e + f*x)), Int(S(1)/sqrt(a*cos(e + f*x) + b), x), x)


def replacement4115(a, b, d, e, f, x):
    return Dist(sqrt(d/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)/sqrt(a + b/sin(e + f*x)), Int(S(1)/sqrt(a*sin(e + f*x) + b), x), x)


def replacement4116(a, b, d, e, f, x):
    return Dist(d*sqrt(d/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)/sqrt(a + b/cos(e + f*x)), Int(S(1)/(sqrt(a*cos(e + f*x) + b)*cos(e + f*x)), x), x)


def replacement4117(a, b, d, e, f, x):
    return Dist(d*sqrt(d/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)/sqrt(a + b/sin(e + f*x)), Int(S(1)/(sqrt(a*sin(e + f*x) + b)*sin(e + f*x)), x), x)


def replacement4118(a, b, d, e, f, n, x):
    return Dist(d**S(3)/(b*(S(2)*n + S(-3))), Int((d/cos(e + f*x))**(n + S(-3))*Simp(S(2)*a*(n + S(-3)) - S(2)*a*(n + S(-2))/cos(e + f*x)**S(2) + b*(S(2)*n + S(-5))/cos(e + f*x), x)/sqrt(a + b/cos(e + f*x)), x), x) + Simp(S(2)*d**S(2)*(d/cos(e + f*x))**(n + S(-2))*sqrt(a + b/cos(e + f*x))*sin(e + f*x)/(b*f*(S(2)*n + S(-3))), x)


def replacement4119(a, b, d, e, f, n, x):
    return Dist(d**S(3)/(b*(S(2)*n + S(-3))), Int((d/sin(e + f*x))**(n + S(-3))*Simp(S(2)*a*(n + S(-3)) - S(2)*a*(n + S(-2))/sin(e + f*x)**S(2) + b*(S(2)*n + S(-5))/sin(e + f*x), x)/sqrt(a + b/sin(e + f*x)), x), x) + Simp(-S(2)*d**S(2)*(d/sin(e + f*x))**(n + S(-2))*sqrt(a + b/sin(e + f*x))*cos(e + f*x)/(b*f*(S(2)*n + S(-3))), x)


def replacement4120(a, b, e, f, x):
    return -Dist(b/(S(2)*a), Int((S(1) + cos(e + f*x)**(S(-2)))/sqrt(a + b/cos(e + f*x)), x), x) + Simp(sqrt(a + b/cos(e + f*x))*sin(e + f*x)/(a*f), x)


def replacement4121(a, b, e, f, x):
    return -Dist(b/(S(2)*a), Int((S(1) + sin(e + f*x)**(S(-2)))/sqrt(a + b/sin(e + f*x)), x), x) - Simp(sqrt(a + b/sin(e + f*x))*cos(e + f*x)/(a*f), x)


def replacement4122(a, b, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) - Dist(b/(a*d), Int(sqrt(d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4123(a, b, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) - Dist(b/(a*d), Int(sqrt(d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4124(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*a*d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(S(2)*a*(n + S(1))/cos(e + f*x) - b*(S(2)*n + S(1)) + b*(S(2)*n + S(3))/cos(e + f*x)**S(2), x)/sqrt(a + b/cos(e + f*x)), x), x) - Simp((d/cos(e + f*x))**(n + S(1))*sqrt(a + b/cos(e + f*x))*sin(e + f*x)/(a*d*f*n), x)


def replacement4125(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*a*d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(S(2)*a*(n + S(1))/sin(e + f*x) - b*(S(2)*n + S(1)) + b*(S(2)*n + S(3))/sin(e + f*x)**S(2), x)/sqrt(a + b/sin(e + f*x)), x), x) + Simp((d/sin(e + f*x))**(n + S(1))*sqrt(a + b/sin(e + f*x))*cos(e + f*x)/(a*d*f*n), x)


def replacement4126(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(a*b*(S(2)*n + S(-1)) + a*b*(S(2)*n + S(3))/cos(e + f*x)**S(2) + S(2)*(a**S(2)*(n + S(1)) + b**S(2)*n)/cos(e + f*x), x)/sqrt(a + b/cos(e + f*x)), x), x) - Simp(a*(d/cos(e + f*x))**n*sqrt(a + b/cos(e + f*x))*tan(e + f*x)/(f*n), x)


def replacement4127(a, b, d, e, f, n, x):
    return Dist(S(1)/(S(2)*d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(a*b*(S(2)*n + S(-1)) + a*b*(S(2)*n + S(3))/sin(e + f*x)**S(2) + S(2)*(a**S(2)*(n + S(1)) + b**S(2)*n)/sin(e + f*x), x)/sqrt(a + b/sin(e + f*x)), x), x) + Simp(a*(d/sin(e + f*x))**n*sqrt(a + b/sin(e + f*x))/(f*n*tan(e + f*x)), x)


def replacement4128(a, b, d, e, f, m, n, x):
    return Dist(d**S(3)/(b*(m + n + S(-1))), Int((d/cos(e + f*x))**(n + S(-3))*(a + b/cos(e + f*x))**m*Simp(a*(n + S(-3)) - a*(n + S(-2))/cos(e + f*x)**S(2) + b*(m + n + S(-2))/cos(e + f*x), x), x), x) + Simp(d**S(3)*(d/cos(e + f*x))**(n + S(-3))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + n + S(-1))), x)


def replacement4129(a, b, d, e, f, m, n, x):
    return Dist(d**S(3)/(b*(m + n + S(-1))), Int((d/sin(e + f*x))**(n + S(-3))*(a + b/sin(e + f*x))**m*Simp(a*(n + S(-3)) - a*(n + S(-2))/sin(e + f*x)**S(2) + b*(m + n + S(-2))/sin(e + f*x), x), x), x) - Simp(d**S(3)*(d/sin(e + f*x))**(n + S(-3))*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4130(a, b, d, e, f, m, n, x):
    return Dist(d/(m + n + S(-1)), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(-2))*Simp(a*b*(n + S(-1)) + a*b*(S(2)*m + n + S(-2))/cos(e + f*x)**S(2) + (a**S(2)*(m + n + S(-1)) + b**S(2)*(m + n + S(-2)))/cos(e + f*x), x), x), x) + Simp(b*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*(m + n + S(-1))), x)


def replacement4131(a, b, d, e, f, m, n, x):
    return Dist(d/(m + n + S(-1)), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(-2))*Simp(a*b*(n + S(-1)) + a*b*(S(2)*m + n + S(-2))/sin(e + f*x)**S(2) + (a**S(2)*(m + n + S(-1)) + b**S(2)*(m + n + S(-2)))/sin(e + f*x), x), x), x) - Simp(b*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(-1))/(f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4132(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n + S(-1))), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(a*b*m/cos(e + f*x)**S(2) + a*b*(n + S(-2)) + b**S(2)*(m + n + S(-2))/cos(e + f*x), x), x), x) + Simp(d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(-1))), x)


def replacement4133(a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n + S(-1))), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(a*b*m/sin(e + f*x)**S(2) + a*b*(n + S(-2)) + b**S(2)*(m + n + S(-2))/sin(e + f*x), x), x), x) - Simp(d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**m/(f*(m + n + S(-1))*tan(e + f*x)), x)


def replacement4134(a, b, d, e, f, x):
    return Dist(a, Int(sqrt(a + b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) + Dist(b/d, Int(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x)), x), x)


def replacement4135(a, b, d, e, f, x):
    return Dist(a, Int(sqrt(a + b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) + Dist(b/d, Int(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x)), x), x)


def replacement4136(a, b, d, e, f, m, n, x):
    return Dist(a, Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1)), x), x) + Dist(b/d, Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1)), x), x)


def replacement4137(a, b, d, e, f, m, n, x):
    return Dist(a, Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1)), x), x) + Dist(b/d, Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1)), x), x)


def replacement4138(a, b, d, e, f, m, n, x):
    return Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m, x)


def replacement4139(a, b, d, e, f, m, n, x):
    return Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m, x)


def replacement4140(a, b, e, f, g, m, p, x):
    return Int((g*sin(e + f*x))**p*(a*cos(e + f*x) + b)**m*cos(e + f*x)**(-m), x)


def replacement4141(a, b, e, f, g, m, p, x):
    return Int((g*cos(e + f*x))**p*(a*sin(e + f*x) + b)**m*sin(e + f*x)**(-m), x)


def replacement4142(a, b, e, f, m, p, x):
    return Dist(b**(S(1) - p)/f, Subst(Int(x**(-p + S(-1))*(-a + b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)


def replacement4143(a, b, e, f, m, p, x):
    return -Dist(b**(S(1) - p)/f, Subst(Int(x**(-p + S(-1))*(-a + b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)


def replacement4144(a, b, e, f, m, p, x):
    return Dist(S(1)/f, Subst(Int(x**(-p + S(-1))*(a + b*x)**m*(x + S(-1))**(p/S(2) + S(-1)/2)*(x + S(1))**(p/S(2) + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)


def replacement4145(a, b, e, f, m, p, x):
    return -Dist(S(1)/f, Subst(Int(x**(-p + S(-1))*(a + b*x)**m*(x + S(-1))**(p/S(2) + S(-1)/2)*(x + S(1))**(p/S(2) + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)


def replacement4146(a, b, e, f, m, x):
    return Dist(b*m, Int((a + b/cos(e + f*x))**(m + S(-1))/cos(e + f*x), x), x) - Simp((a + b/cos(e + f*x))**m/(f*tan(e + f*x)), x)


def replacement4147(a, b, e, f, m, x):
    return Dist(b*m, Int((a + b/sin(e + f*x))**(m + S(-1))/sin(e + f*x), x), x) + Simp((a + b/sin(e + f*x))**m*tan(e + f*x)/f, x)


def replacement4148(a, b, e, f, g, m, p, x):
    return Dist((a + b/cos(e + f*x))**FracPart(m)*(a*cos(e + f*x) + b)**(-FracPart(m))*cos(e + f*x)**FracPart(m), Int((g*sin(e + f*x))**p*(a*cos(e + f*x) + b)**m*cos(e + f*x)**(-m), x), x)


def replacement4149(a, b, e, f, g, m, p, x):
    return Dist((a + b/sin(e + f*x))**FracPart(m)*(a*sin(e + f*x) + b)**(-FracPart(m))*sin(e + f*x)**FracPart(m), Int((g*cos(e + f*x))**p*(a*sin(e + f*x) + b)**m*sin(e + f*x)**(-m), x), x)


def replacement4150(a, b, e, f, g, m, p, x):
    return Int((g*sin(e + f*x))**p*(a + b/cos(e + f*x))**m, x)


def replacement4151(a, b, e, f, g, m, p, x):
    return Int((g*cos(e + f*x))**p*(a + b/sin(e + f*x))**m, x)


def replacement4152(a, b, e, f, g, m, p, x):
    return Dist((g/sin(e + f*x))**p*(g*sin(e + f*x))**p, Int((g*sin(e + f*x))**(-p)*(a + b/cos(e + f*x))**m, x), x)


def replacement4153(a, b, e, f, g, m, p, x):
    return Dist((g/cos(e + f*x))**p*(g*cos(e + f*x))**p, Int((g*cos(e + f*x))**(-p)*(a + b/sin(e + f*x))**m, x), x)


def replacement4154(a, b, c, d, m, n, x):
    return -Dist(a**(-m + n + S(1))*b**(-n)/d, Subst(Int(x**(-m - n)*(a - b*x)**(m/S(2) + S(-1)/2)*(a + b*x)**(m/S(2) + n + S(-1)/2), x), x, cos(c + d*x)), x)


def replacement4155(a, b, c, d, m, n, x):
    return Dist(a**(-m + n + S(1))*b**(-n)/d, Subst(Int(x**(-m - n)*(a - b*x)**(m/S(2) + S(-1)/2)*(a + b*x)**(m/S(2) + n + S(-1)/2), x), x, sin(c + d*x)), x)


def replacement4156(a, b, c, d, m, n, x):
    return Dist(b**(S(1) - m)/d, Subst(Int((-a + b*x)**(m/S(2) + S(-1)/2)*(a + b*x)**(m/S(2) + n + S(-1)/2)/x, x), x, S(1)/cos(c + d*x)), x)


def replacement4157(a, b, c, d, m, n, x):
    return -Dist(b**(S(1) - m)/d, Subst(Int((-a + b*x)**(m/S(2) + S(-1)/2)*(a + b*x)**(m/S(2) + n + S(-1)/2)/x, x), x, S(1)/sin(c + d*x)), x)


def replacement4158(a, b, c, d, e, m, x):
    return -Dist(e**S(2)/m, Int((e*tan(c + d*x))**(m + S(-2))*(a*m + b*(m + S(-1))/cos(c + d*x)), x), x) + Simp(e*(e*tan(c + d*x))**(m + S(-1))*(a*m + b*(m + S(-1))/cos(c + d*x))/(d*m*(m + S(-1))), x)


def replacement4159(a, b, c, d, e, m, x):
    return -Dist(e**S(2)/m, Int((e/tan(c + d*x))**(m + S(-2))*(a*m + b*(m + S(-1))/sin(c + d*x)), x), x) - Simp(e*(e/tan(c + d*x))**(m + S(-1))*(a*m + b*(m + S(-1))/sin(c + d*x))/(d*m*(m + S(-1))), x)


def replacement4160(a, b, c, d, e, m, x):
    return -Dist(S(1)/(e**S(2)*(m + S(1))), Int((e*tan(c + d*x))**(m + S(2))*(a*(m + S(1)) + b*(m + S(2))/cos(c + d*x)), x), x) + Simp((e*tan(c + d*x))**(m + S(1))*(a + b/cos(c + d*x))/(d*e*(m + S(1))), x)


def replacement4161(a, b, c, d, e, m, x):
    return -Dist(S(1)/(e**S(2)*(m + S(1))), Int((e/tan(c + d*x))**(m + S(2))*(a*(m + S(1)) + b*(m + S(2))/sin(c + d*x)), x), x) - Simp((e/tan(c + d*x))**(m + S(1))*(a + b/sin(c + d*x))/(d*e*(m + S(1))), x)


def replacement4162(a, b, c, d, x):
    return Int((a*cos(c + d*x) + b)/sin(c + d*x), x)


def replacement4163(a, b, c, d, x):
    return Int((a*sin(c + d*x) + b)/cos(c + d*x), x)


def replacement4164(a, b, c, d, e, m, x):
    return Dist(a, Int((e*tan(c + d*x))**m, x), x) + Dist(b, Int((e*tan(c + d*x))**m/cos(c + d*x), x), x)


def replacement4165(a, b, c, d, e, m, x):
    return Dist(a, Int((e/tan(c + d*x))**m, x), x) + Dist(b, Int((e/tan(c + d*x))**m/sin(c + d*x), x), x)


def replacement4166(a, b, c, d, m, n, x):
    return Dist((S(-1))**(m/S(2) + S(-1)/2)*b**(S(1) - m)/d, Subst(Int((a + x)**n*(b**S(2) - x**S(2))**(m/S(2) + S(-1)/2)/x, x), x, b/cos(c + d*x)), x)


def replacement4167(a, b, c, d, m, n, x):
    return -Dist((S(-1))**(m/S(2) + S(-1)/2)*b**(S(1) - m)/d, Subst(Int((a + x)**n*(b**S(2) - x**S(2))**(m/S(2) + S(-1)/2)/x, x), x, b/sin(c + d*x)), x)


def replacement4168(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((e*tan(c + d*x))**m, (a + b/cos(c + d*x))**n, x), x)


def replacement4169(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((e/tan(c + d*x))**m, (a + b/sin(c + d*x))**n, x), x)


def replacement4170(a, b, c, d, m, n, x):
    return Dist(S(2)*a**(m/S(2) + n + S(1)/2)/d, Subst(Int(x**m*(a*x**S(2) + S(2))**(m/S(2) + n + S(-1)/2)/(a*x**S(2) + S(1)), x), x, tan(c + d*x)/sqrt(a + b/cos(c + d*x))), x)


def replacement4171(a, b, c, d, m, n, x):
    return Dist(-S(2)*a**(m/S(2) + n + S(1)/2)/d, Subst(Int(x**m*(a*x**S(2) + S(2))**(m/S(2) + n + S(-1)/2)/(a*x**S(2) + S(1)), x), x, S(1)/(sqrt(a + b/sin(c + d*x))*tan(c + d*x))), x)


def replacement4172(a, b, c, d, e, m, n, x):
    return Dist(a**(S(2)*n)*e**(-S(2)*n), Int((e*tan(c + d*x))**(m + S(2)*n)*(-a + b/cos(c + d*x))**(-n), x), x)


def replacement4173(a, b, c, d, e, m, n, x):
    return Dist(a**(S(2)*n)*e**(-S(2)*n), Int((e/tan(c + d*x))**(m + S(2)*n)*(-a + b/sin(c + d*x))**(-n), x), x)


def replacement4174(a, b, c, d, e, m, n, x):
    return Simp(S(2)**(m + n + S(1))*(a/(a + b/cos(c + d*x)))**(m + n + S(1))*(e*tan(c + d*x))**(m + S(1))*(a + b/cos(c + d*x))**n*AppellF1(m/S(2) + S(1)/2, m + n, S(1), m/S(2) + S(3)/2, -(a - b/cos(c + d*x))/(a + b/cos(c + d*x)), (a - b/cos(c + d*x))/(a + b/cos(c + d*x)))/(d*e*(m + S(1))), x)


def replacement4175(a, b, c, d, e, m, n, x):
    return -Simp(S(2)**(m + n + S(1))*(a/(a + b/sin(c + d*x)))**(m + n + S(1))*(e/tan(c + d*x))**(m + S(1))*(a + b/sin(c + d*x))**n*AppellF1(m/S(2) + S(1)/2, m + n, S(1), m/S(2) + S(3)/2, -(a - b/sin(c + d*x))/(a + b/sin(c + d*x)), (a - b/sin(c + d*x))/(a + b/sin(c + d*x)))/(d*e*(m + S(1))), x)


def replacement4176(a, b, c, d, e, x):
    return Dist(S(1)/a, Int(sqrt(e*tan(c + d*x)), x), x) - Dist(b/a, Int(sqrt(e*tan(c + d*x))/(a*cos(c + d*x) + b), x), x)


def replacement4177(a, b, c, d, e, x):
    return Dist(S(1)/a, Int(sqrt(e/tan(c + d*x)), x), x) - Dist(b/a, Int(sqrt(e/tan(c + d*x))/(a*sin(c + d*x) + b), x), x)


def replacement4178(a, b, c, d, e, m, x):
    return -Dist(e**S(2)/b**S(2), Int((e*tan(c + d*x))**(m + S(-2))*(a - b/cos(c + d*x)), x), x) + Dist(e**S(2)*(a**S(2) - b**S(2))/b**S(2), Int((e*tan(c + d*x))**(m + S(-2))/(a + b/cos(c + d*x)), x), x)


def replacement4179(a, b, c, d, e, m, x):
    return -Dist(e**S(2)/b**S(2), Int((e/tan(c + d*x))**(m + S(-2))*(a - b/sin(c + d*x)), x), x) + Dist(e**S(2)*(a**S(2) - b**S(2))/b**S(2), Int((e/tan(c + d*x))**(m + S(-2))/(a + b/sin(c + d*x)), x), x)


def replacement4180(a, b, c, d, e, x):
    return Dist(S(1)/a, Int(S(1)/sqrt(e*tan(c + d*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(e*tan(c + d*x))*(a*cos(c + d*x) + b)), x), x)


def replacement4181(a, b, c, d, e, x):
    return Dist(S(1)/a, Int(S(1)/sqrt(e/tan(c + d*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(e/tan(c + d*x))*(a*sin(c + d*x) + b)), x), x)


def replacement4182(a, b, c, d, e, m, x):
    return Dist(b**S(2)/(e**S(2)*(a**S(2) - b**S(2))), Int((e*tan(c + d*x))**(m + S(2))/(a + b/cos(c + d*x)), x), x) + Dist(S(1)/(a**S(2) - b**S(2)), Int((e*tan(c + d*x))**m*(a - b/cos(c + d*x)), x), x)


def replacement4183(a, b, c, d, e, m, x):
    return Dist(b**S(2)/(e**S(2)*(a**S(2) - b**S(2))), Int((e/tan(c + d*x))**(m + S(2))/(a + b/sin(c + d*x)), x), x) + Dist(S(1)/(a**S(2) - b**S(2)), Int((e/tan(c + d*x))**m*(a - b/sin(c + d*x)), x), x)


def replacement4184(a, b, c, d, n, x):
    return Int((S(-1) + cos(c + d*x)**(S(-2)))*(a + b/cos(c + d*x))**n, x)


def replacement4185(a, b, c, d, n, x):
    return Int((S(-1) + sin(c + d*x)**(S(-2)))*(a + b/sin(c + d*x))**n, x)


def replacement4186(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((e*tan(c + d*x))**m, (a + b/cos(c + d*x))**n, x), x)


def replacement4187(a, b, c, d, e, m, n, x):
    return Int(ExpandIntegrand((e/tan(c + d*x))**m, (a + b/sin(c + d*x))**n, x), x)


def replacement4188(a, b, c, d, m, n, x):
    return Int((a*cos(c + d*x) + b)**n*sin(c + d*x)**m*cos(c + d*x)**(-m - n), x)


def replacement4189(a, b, c, d, m, n, x):
    return Int((a*sin(c + d*x) + b)**n*sin(c + d*x)**(-m - n)*cos(c + d*x)**m, x)


def replacement4190(a, b, c, d, e, m, n, x):
    return Int((e*tan(c + d*x))**m*(a + b/cos(c + d*x))**n, x)


def replacement4191(a, b, c, d, e, m, n, x):
    return Int((e/tan(c + d*x))**m*(a + b/sin(c + d*x))**n, x)


def replacement4192(a, b, c, d, e, m, n, p, x):
    return Dist((e*tan(c + d*x))**(-m*p)*(e*tan(c + d*x)**p)**m, Int((e*tan(c + d*x))**(m*p)*(a + b/cos(c + d*x))**n, x), x)


def replacement4193(a, b, c, d, e, m, n, p, x):
    return Dist((e*(S(1)/tan(c + d*x))**p)**m*(e/tan(c + d*x))**(-m*p), Int((e/tan(c + d*x))**(m*p)*(a + b/sin(c + d*x))**n, x), x)


def replacement4194(a, b, c, d, e, f, m, n, x):
    return Dist(c**n, Int(ExpandTrig((S(1) + d/(c*cos(e + f*x)))**n, (a + b/cos(e + f*x))**m, x), x), x)


def replacement4195(a, b, c, d, e, f, m, n, x):
    return Dist(c**n, Int(ExpandTrig((S(1) + d/(c*sin(e + f*x)))**n, (a + b/sin(e + f*x))**m, x), x), x)


def replacement4196(a, b, c, d, e, f, m, n, x):
    return Dist((-a*c)**m, Int((c + d/cos(e + f*x))**(-m + n)*tan(e + f*x)**(S(2)*m), x), x)


def replacement4197(a, b, c, d, e, f, m, n, x):
    return Dist((-a*c)**m, Int((c + d/sin(e + f*x))**(-m + n)*(S(1)/tan(e + f*x))**(S(2)*m), x), x)


def replacement4198(a, b, c, d, e, f, m, x):
    return Dist((-a*c)**(m + S(1)/2)*tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Int(tan(e + f*x)**(S(2)*m), x), x)


def replacement4199(a, b, c, d, e, f, m, x):
    return Dist((-a*c)**(m + S(1)/2)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Int((S(1)/tan(e + f*x))**(S(2)*m), x), x)


def replacement4200(a, b, c, d, e, f, n, x):
    return Dist(c, Int(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**(n + S(-1)), x), x) + Simp(-S(2)*a*c*(c + d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(-1))), x)


def replacement4201(a, b, c, d, e, f, n, x):
    return Dist(c, Int(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**(n + S(-1)), x), x) + Simp(S(2)*a*c*(c + d/sin(e + f*x))**(n + S(-1))/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(-1))*tan(e + f*x)), x)


def replacement4202(a, b, c, d, e, f, n, x):
    return Dist(S(1)/c, Int(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**(n + S(1)), x), x) + Simp(S(2)*a*(c + d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4203(a, b, c, d, e, f, n, x):
    return Dist(S(1)/c, Int(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**(n + S(1)), x), x) + Simp(-S(2)*a*(c + d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4204(a, b, c, d, e, f, n, x):
    return Dist(a/c, Int(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**(n + S(1)), x), x) + Simp(S(4)*a**S(2)*(c + d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4205(a, b, c, d, e, f, n, x):
    return Dist(a/c, Int(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**(n + S(1)), x), x) + Simp(-S(4)*a**S(2)*(c + d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4206(a, b, c, d, e, f, n, x):
    return Dist(a, Int(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**n, x), x) + Simp(S(2)*a**S(2)*(c + d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4207(a, b, c, d, e, f, n, x):
    return Dist(a, Int(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**n, x), x) + Simp(-S(2)*a**S(2)*(c + d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4208(a, b, c, d, e, f, n, x):
    return Dist(a**S(2)/c**S(2), Int(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**(n + S(2)), x), x) + Simp(S(8)*a**S(3)*(c + d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4209(a, b, c, d, e, f, n, x):
    return Dist(a**S(2)/c**S(2), Int(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**(n + S(2)), x), x) + Simp(-S(8)*a**S(3)*(c + d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4210(a, b, c, d, e, f, m, n, x):
    return Dist(a*c*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Subst(Int(x**(-m - n)*(a*x + b)**(m + S(-1)/2)*(c*x + d)**(n + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement4211(a, b, c, d, e, f, m, n, x):
    return -Dist(a*c/(f*sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Subst(Int(x**(-m - n)*(a*x + b)**(m + S(-1)/2)*(c*x + d)**(n + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement4212(a, b, c, d, e, f, m, n, x):
    return -Dist(a*c*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2)/x, x), x, S(1)/cos(e + f*x)), x)


def replacement4213(a, b, c, d, e, f, m, n, x):
    return Dist(a*c/(f*sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2)/x, x), x, S(1)/sin(e + f*x)), x)


def replacement4214(a, b, c, d, e, f, x):
    return Dist(b*d, Int(cos(e + f*x)**(S(-2)), x), x) + Simp(a*c*x, x)


def replacement4215(a, b, c, d, e, f, x):
    return Dist(b*d, Int(sin(e + f*x)**(S(-2)), x), x) + Simp(a*c*x, x)


def replacement4216(a, b, c, d, e, f, x):
    return Dist(b*d, Int(cos(e + f*x)**(S(-2)), x), x) + Dist(a*d + b*c, Int(S(1)/cos(e + f*x), x), x) + Simp(a*c*x, x)


def replacement4217(a, b, c, d, e, f, x):
    return Dist(b*d, Int(sin(e + f*x)**(S(-2)), x), x) + Dist(a*d + b*c, Int(S(1)/sin(e + f*x), x), x) + Simp(a*c*x, x)


def replacement4218(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b/cos(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b/cos(e + f*x))/cos(e + f*x), x), x)


def replacement4219(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b/sin(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b/sin(e + f*x))/sin(e + f*x), x), x)


def replacement4220(a, b, c, d, e, f, x):
    return Dist(a*c, Int(S(1)/sqrt(a + b/cos(e + f*x)), x), x) + Int((a*d + b*c + b*d/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x)


def replacement4221(a, b, c, d, e, f, x):
    return Dist(a*c, Int(S(1)/sqrt(a + b/sin(e + f*x)), x), x) + Int((a*d + b*c + b*d/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x)


def replacement4222(a, b, c, d, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/cos(e + f*x))**(m + S(-1))*Simp(a*c*m + (a*d*(S(2)*m + S(-1)) + b*c*m)/cos(e + f*x), x), x), x) + Simp(b*d*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*m), x)


def replacement4223(a, b, c, d, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/sin(e + f*x))**(m + S(-1))*Simp(a*c*m + (a*d*(S(2)*m + S(-1)) + b*c*m)/sin(e + f*x), x), x), x) - Simp(b*d*(a + b/sin(e + f*x))**(m + S(-1))/(f*m*tan(e + f*x)), x)


def replacement4224(a, b, c, d, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/cos(e + f*x))**(m + S(-2))*Simp(a**S(2)*c*m + b*(a*d*(S(2)*m + S(-1)) + b*c*m)/cos(e + f*x)**S(2) + (a**S(2)*d*m + S(2)*a*b*c*m + b**S(2)*d*(m + S(-1)))/cos(e + f*x), x), x), x) + Simp(b*d*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*m), x)


def replacement4225(a, b, c, d, e, f, m, x):
    return Dist(S(1)/m, Int((a + b/sin(e + f*x))**(m + S(-2))*Simp(a**S(2)*c*m + b*(a*d*(S(2)*m + S(-1)) + b*c*m)/sin(e + f*x)**S(2) + (a**S(2)*d*m + S(2)*a*b*c*m + b**S(2)*d*(m + S(-1)))/sin(e + f*x), x), x), x) - Simp(b*d*(a + b/sin(e + f*x))**(m + S(-1))/(f*m*tan(e + f*x)), x)


def replacement4226(a, b, c, d, e, f, x):
    return -Dist((-a*d + b*c)/a, Int(S(1)/((a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Simp(c*x/a, x)


def replacement4227(a, b, c, d, e, f, x):
    return -Dist((-a*d + b*c)/a, Int(S(1)/((a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Simp(c*x/a, x)


def replacement4228(a, b, c, d, e, f, x):
    return Dist(c/a, Int(sqrt(a + b/cos(e + f*x)), x), x) - Dist((-a*d + b*c)/a, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4229(a, b, c, d, e, f, x):
    return Dist(c/a, Int(sqrt(a + b/sin(e + f*x)), x), x) - Dist((-a*d + b*c)/a, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4230(a, b, c, d, e, f, x):
    return Dist(c, Int(S(1)/sqrt(a + b/cos(e + f*x)), x), x) + Dist(d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4231(a, b, c, d, e, f, x):
    return Dist(c, Int(S(1)/sqrt(a + b/sin(e + f*x)), x), x) + Dist(d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4232(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(a*c*(S(2)*m + S(1)) - (m + S(1))*(-a*d + b*c)/cos(e + f*x), x), x), x) + Simp((a + b/cos(e + f*x))**m*(-a*d + b*c)*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4233(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(a*c*(S(2)*m + S(1)) - (m + S(1))*(-a*d + b*c)/sin(e + f*x), x), x), x) - Simp((a + b/sin(e + f*x))**m*(-a*d + b*c)/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4234(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-a*(m + S(1))*(-a*d + b*c)/cos(e + f*x) + b*(m + S(2))*(-a*d + b*c)/cos(e + f*x)**S(2) + c*(a**S(2) - b**S(2))*(m + S(1)), x), x), x) - Simp(b*(a + b/cos(e + f*x))**(m + S(1))*(-a*d + b*c)*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4235(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-a*(m + S(1))*(-a*d + b*c)/sin(e + f*x) + b*(m + S(2))*(-a*d + b*c)/sin(e + f*x)**S(2) + c*(a**S(2) - b**S(2))*(m + S(1)), x), x), x) + Simp(b*(a + b/sin(e + f*x))**(m + S(1))*(-a*d + b*c)/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4236(a, b, c, d, e, f, m, x):
    return Dist(c, Int((a + b/cos(e + f*x))**m, x), x) + Dist(d, Int((a + b/cos(e + f*x))**m/cos(e + f*x), x), x)


def replacement4237(a, b, c, d, e, f, m, x):
    return Dist(c, Int((a + b/sin(e + f*x))**m, x), x) + Dist(d, Int((a + b/sin(e + f*x))**m/sin(e + f*x), x), x)


def replacement4238(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/cos(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4239(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/sin(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4240(a, b, c, d, e, f, x):
    return Dist(a/c, Int(S(1)/sqrt(a + b/cos(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4241(a, b, c, d, e, f, x):
    return Dist(a/c, Int(S(1)/sqrt(a + b/sin(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4242(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(a + b/cos(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4243(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(a + b/sin(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4244(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*d), Int((a**S(2)*d + b**S(2)*c/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) - Dist((-a*d + b*c)**S(2)/(c*d), Int(S(1)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4245(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*d), Int((a**S(2)*d + b**S(2)*c/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) - Dist((-a*d + b*c)**S(2)/(c*d), Int(S(1)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4246(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*(-a*d + b*c)), Int((-a*d + b*c - b*d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist(d**S(2)/(c*(-a*d + b*c)), Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4247(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*(-a*d + b*c)), Int((-a*d + b*c - b*d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist(d**S(2)/(c*(-a*d + b*c)), Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4248(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(S(1)/sqrt(a + b/cos(e + f*x)), x), x) - Dist(d/c, Int(S(1)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4249(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(S(1)/sqrt(a + b/sin(e + f*x)), x), x) - Dist(d/c, Int(S(1)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4250(a, b, c, d, e, f, x):
    return Dist(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))/tan(e + f*x), Int(tan(e + f*x), x), x)


def replacement4251(a, b, c, d, e, f, x):
    return Dist(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x), Int(S(1)/tan(e + f*x), x), x)


def replacement4252(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b/cos(e + f*x))/sqrt(c + d/cos(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b/cos(e + f*x))/(sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4253(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b/sin(e + f*x))/sqrt(c + d/sin(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b/sin(e + f*x))/(sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4254(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/cos(e + f*x))/(sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4255(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/sin(e + f*x))/(sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4256(a, b, c, d, e, f, x):
    return Dist(S(2)*a/f, Subst(Int(S(1)/(a*c*x**S(2) + S(1)), x), x, tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x)))), x)


def replacement4257(a, b, c, d, e, f, x):
    return Dist(-S(2)*a/f, Subst(Int(S(1)/(a*c*x**S(2) + S(1)), x), x, S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x))), x)


def replacement4258(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(c + d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4259(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(c + d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4260(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((S(1) + S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b/cos(e + f*x))*(c - d)))*sqrt(-(S(1) - S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b/cos(e + f*x))*(c + d)))*(a + b/cos(e + f*x))*EllipticPi(a*(c + d)/(c*(a + b)), asin(sqrt(c + d/cos(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b/cos(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(c*f*Rt((a + b)/(c + d), S(2))*tan(e + f*x)), x)


def replacement4261(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((S(1) + S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b/sin(e + f*x))*(c - d)))*sqrt(-(S(1) - S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b/sin(e + f*x))*(c + d)))*(a + b/sin(e + f*x))*EllipticPi(a*(c + d)/(c*(a + b)), asin(sqrt(c + d/sin(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b/sin(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))*tan(e + f*x)/(c*f*Rt((a + b)/(c + d), S(2))), x)


def replacement4262(a, b, c, d, e, f, x):
    return Dist(tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Int(S(1)/tan(e + f*x), x), x)


def replacement4263(a, b, c, d, e, f, x):
    return Dist(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Int(tan(e + f*x), x), x)


def replacement4264(a, b, c, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b/cos(e + f*x))/sqrt(c + d/cos(e + f*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4265(a, b, c, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b/sin(e + f*x))/sqrt(c + d/sin(e + f*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4266(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/cos(e + f*x))/sqrt(c + d/cos(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))**(S(3)/2)*cos(e + f*x)), x), x)


def replacement4267(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b/sin(e + f*x))/sqrt(c + d/sin(e + f*x)), x), x) - Dist(d/c, Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))**(S(3)/2)*sin(e + f*x)), x), x)


def replacement4268(a, b, c, d, e, f, m, n, x):
    return -Dist(a**S(2)*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**n/(x*sqrt(a - b*x)), x), x, S(1)/cos(e + f*x)), x)


def replacement4269(a, b, c, d, e, f, m, n, x):
    return Dist(a**S(2)*cos(e + f*x)/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**n/(x*sqrt(a - b*x)), x), x, S(1)/sin(e + f*x)), x)


def replacement4270(a, b, c, d, e, f, m, n, x):
    return Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-m - n), x)


def replacement4271(a, b, c, d, e, f, m, n, x):
    return Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-m - n), x)


def replacement4272(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(a + b/cos(e + f*x))*sqrt(c*cos(e + f*x) + d)/(sqrt(c + d/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)), Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-m - n), x), x)


def replacement4273(a, b, c, d, e, f, m, n, x):
    return Dist(sqrt(a + b/sin(e + f*x))*sqrt(c*sin(e + f*x) + d)/(sqrt(c + d/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)), Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-m - n), x), x)


def replacement4274(a, b, c, d, e, f, m, n, x):
    return Dist((a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n*(a*cos(e + f*x) + b)**(-m)*(c*cos(e + f*x) + d)**(-n)*cos(e + f*x)**(m + n), Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-m - n), x), x)


def replacement4275(a, b, c, d, e, f, m, n, x):
    return Dist((a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n*(a*sin(e + f*x) + b)**(-m)*(c*sin(e + f*x) + d)**(-n)*sin(e + f*x)**(m + n), Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-m - n), x), x)


def replacement4276(a, b, c, d, e, f, m, n, x):
    return Int(ExpandTrig((a + b/cos(e + f*x))**m, (c + d/cos(e + f*x))**n, x), x)


def replacement4277(a, b, c, d, e, f, m, n, x):
    return Int(ExpandTrig((a + b/sin(e + f*x))**m, (c + d/sin(e + f*x))**n, x), x)


def replacement4278(a, b, c, d, e, f, m, n, x):
    return Int((a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n, x)


def replacement4279(a, b, c, d, e, f, m, n, x):
    return Int((a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n, x)


def replacement4280(a, b, d, e, f, m, n, x):
    return Dist(d**m, Int((d*cos(e + f*x))**(-m + n)*(a*cos(e + f*x) + b)**m, x), x)


def replacement4281(a, b, d, e, f, m, n, x):
    return Dist(d**m, Int((d*sin(e + f*x))**(-m + n)*(a*sin(e + f*x) + b)**m, x), x)


def replacement4282(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/cos(e + f*x))**p)**FracPart(n)*(d/cos(e + f*x))**(-p*FracPart(n)), Int((d/cos(e + f*x))**(n*p)*(a + b/cos(e + f*x))**m, x), x)


def replacement4283(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/sin(e + f*x))**p)**FracPart(n)*(d/sin(e + f*x))**(-p*FracPart(n)), Int((d*cos(e + f*x))**(n*p)*(a + b*cos(e + f*x))**m, x), x)


def replacement4284(a, b, c, d, e, f, m, n, x):
    return -Simp(b*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4285(a, b, c, d, e, f, m, n, x):
    return Simp(b*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4286(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(c + d/cos(e + f*x))**n/cos(e + f*x), x), x) - Simp(b*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4287(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(c + d/sin(e + f*x))**n/sin(e + f*x), x), x) + Simp(b*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4288(a, b, c, d, e, f, x):
    return -Simp(a*c*log(S(1) + b/(a*cos(e + f*x)))*tan(e + f*x)/(b*f*sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), x)


def replacement4289(a, b, c, d, e, f, x):
    return Simp(a*c*log(S(1) + b/(a*sin(e + f*x)))/(b*f*sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), x)


def replacement4290(a, b, c, d, e, f, m, x):
    return Simp(-S(2)*a*c*(a + b/cos(e + f*x))**m*tan(e + f*x)/(b*f*sqrt(c + d/cos(e + f*x))*(S(2)*m + S(1))), x)


def replacement4291(a, b, c, d, e, f, m, x):
    return Simp(S(2)*a*c*(a + b/sin(e + f*x))**m/(b*f*sqrt(c + d/sin(e + f*x))*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4292(a, b, c, d, e, f, m, n, x):
    return -Dist(d*(S(2)*n + S(-1))/(b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(c + d/cos(e + f*x))**(n + S(-1))/cos(e + f*x), x), x) + Simp(-S(2)*a*c*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4293(a, b, c, d, e, f, m, n, x):
    return -Dist(d*(S(2)*n + S(-1))/(b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(c + d/sin(e + f*x))**(n + S(-1))/sin(e + f*x), x), x) + Simp(S(2)*a*c*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**(n + S(-1))/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4294(a, b, c, d, e, f, m, n, x):
    return Dist(c*(S(2)*n + S(-1))/(m + n), Int((a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**(n + S(-1))/cos(e + f*x), x), x) + Simp(d*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(f*(m + n)), x)


def replacement4295(a, b, c, d, e, f, m, n, x):
    return Dist(c*(S(2)*n + S(-1))/(m + n), Int((a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**(n + S(-1))/sin(e + f*x), x), x) - Simp(d*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**(n + S(-1))/(f*(m + n)*tan(e + f*x)), x)


def replacement4296(a, b, c, d, e, f, n, x):
    return Dist(S(2)*c, Int((c + d/cos(e + f*x))**(n + S(-1))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Simp(S(2)*d*(c + d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(-1))), x)


def replacement4297(a, b, c, d, e, f, n, x):
    return Dist(S(2)*c, Int((c + d/sin(e + f*x))**(n + S(-1))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Simp(-S(2)*d*(c + d/sin(e + f*x))**(n + S(-1))/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(-1))*tan(e + f*x)), x)


def replacement4298(a, b, c, d, e, f, m, n, x):
    return -Dist(d*(S(2)*n + S(-1))/(b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(c + d/cos(e + f*x))**(n + S(-1))/cos(e + f*x), x), x) + Simp(-S(2)*a*c*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**(n + S(-1))*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4299(a, b, c, d, e, f, m, n, x):
    return -Dist(d*(S(2)*n + S(-1))/(b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(c + d/sin(e + f*x))**(n + S(-1))/sin(e + f*x), x), x) + Simp(S(2)*a*c*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**(n + S(-1))/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4300(a, b, c, d, e, f, m, n, x):
    return Dist((-a*c)**m, Int(ExpandTrig(tan(e + f*x)**(S(2)*m)/cos(e + f*x), (c + d/cos(e + f*x))**(-m + n), x), x), x)


def replacement4301(a, b, c, d, e, f, m, n, x):
    return Dist((-a*c)**m, Int(ExpandTrig((S(1)/tan(e + f*x))**(S(2)*m)/sin(e + f*x), (c + d/sin(e + f*x))**(-m + n), x), x), x)


def replacement4302(a, b, c, d, e, f, m, x):
    return Dist((-a*c)**(m + S(1)/2)*tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Int(tan(e + f*x)**(S(2)*m)/cos(e + f*x), x), x)


def replacement4303(a, b, c, d, e, f, m, x):
    return Dist((-a*c)**(m + S(1)/2)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Int((S(1)/tan(e + f*x))**(S(2)*m)/sin(e + f*x), x), x)


def replacement4304(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*(c + d/cos(e + f*x))**n/cos(e + f*x), x), x) - Simp(b*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4305(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*(c + d/sin(e + f*x))**n/sin(e + f*x), x), x) + Simp(b*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4306(a, b, c, d, e, f, m, n, x):
    return -Dist(a*c*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)


def replacement4307(a, b, c, d, e, f, m, n, x):
    return Dist(a*c/(f*sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)


def replacement4308(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((-a*c)**m, Int(ExpandTrig((g/cos(e + f*x))**p*tan(e + f*x)**(S(2)*m), (c + d/cos(e + f*x))**(-m + n), x), x), x)


def replacement4309(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((-a*c)**m, Int(ExpandTrig((g/sin(e + f*x))**p*(S(1)/tan(e + f*x))**(S(2)*m), (c + d/sin(e + f*x))**(-m + n), x), x), x)


def replacement4310(a, b, c, d, e, f, g, m, p, x):
    return Dist((-a*c)**(m + S(1)/2)*tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Int((g/cos(e + f*x))**p*tan(e + f*x)**(S(2)*m), x), x)


def replacement4311(a, b, c, d, e, f, g, m, p, x):
    return Dist((-a*c)**(m + S(1)/2)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Int((g/sin(e + f*x))**p*(S(1)/tan(e + f*x))**(S(2)*m), x), x)


def replacement4312(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(a*c*g*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))), Subst(Int((g*x)**(p + S(-1))*(a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)


def replacement4313(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a*c*g/(f*sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x)), Subst(Int((g*x)**(p + S(-1))*(a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)


def replacement4314(a, b, c, d, e, f, g, x):
    return Dist(S(2)*b*g/f, Subst(Int(S(1)/(a*d + b*c - c*g*x**S(2)), x), x, b*tan(e + f*x)/(sqrt(g/cos(e + f*x))*sqrt(a + b/cos(e + f*x)))), x)


def replacement4315(a, b, c, d, e, f, g, x):
    return Dist(-S(2)*b*g/f, Subst(Int(S(1)/(a*d + b*c - c*g*x**S(2)), x), x, b/(sqrt(g/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4316(a, b, c, d, e, f, g, x):
    return Dist(a/c, Int(sqrt(g/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist((-a*d + b*c)/(c*g), Int((g/cos(e + f*x))**(S(3)/2)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))), x), x)


def replacement4317(a, b, c, d, e, f, g, x):
    return Dist(a/c, Int(sqrt(g/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist((-a*d + b*c)/(c*g), Int((g/sin(e + f*x))**(S(3)/2)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))), x), x)


def replacement4318(a, b, c, d, e, f, x):
    return Dist(S(2)*b/f, Subst(Int(S(1)/(a*d + b*c + d*x**S(2)), x), x, b*tan(e + f*x)/sqrt(a + b/cos(e + f*x))), x)


def replacement4319(a, b, c, d, e, f, x):
    return Dist(-S(2)*b/f, Subst(Int(S(1)/(a*d + b*c + d*x**S(2)), x), x, b/(sqrt(a + b/sin(e + f*x))*tan(e + f*x))), x)


def replacement4320(a, b, c, d, e, f, x):
    return Simp(sqrt(c/(c + d/cos(e + f*x)))*sqrt(a + b/cos(e + f*x))*EllipticE(asin(c*tan(e + f*x)/(c + d/cos(e + f*x))), -(-a*d + b*c)/(a*d + b*c))/(d*f*sqrt(c*d*(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*(a*d + b*c)))), x)


def replacement4321(a, b, c, d, e, f, x):
    return -Simp(sqrt(c/(c + d/sin(e + f*x)))*sqrt(a + b/sin(e + f*x))*EllipticE(asin(c/((c + d/sin(e + f*x))*tan(e + f*x))), -(-a*d + b*c)/(a*d + b*c))/(d*f*sqrt(c*d*(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*(a*d + b*c)))), x)


def replacement4322(a, b, c, d, e, f, x):
    return Dist(b/d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4323(a, b, c, d, e, f, x):
    return Dist(b/d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4324(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int(sqrt(g/cos(e + f*x))*sqrt(a + b/cos(e + f*x)), x), x) - Dist(c*g/d, Int(sqrt(g/cos(e + f*x))*sqrt(a + b/cos(e + f*x))/(c + d/cos(e + f*x)), x), x)


def replacement4325(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int(sqrt(g/sin(e + f*x))*sqrt(a + b/sin(e + f*x)), x), x) - Dist(c*g/d, Int(sqrt(g/sin(e + f*x))*sqrt(a + b/sin(e + f*x))/(c + d/sin(e + f*x)), x), x)


def replacement4326(a, b, c, d, e, f, g, x):
    return Dist(b/d, Int((g/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int((g/cos(e + f*x))**(S(3)/2)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))), x), x)


def replacement4327(a, b, c, d, e, f, g, x):
    return Dist(b/d, Int((g/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int((g/sin(e + f*x))**(S(3)/2)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))), x), x)


def replacement4328(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4329(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4330(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((a + b/cos(e + f*x))/(a + b))*EllipticPi(S(2)*d/(c + d), asin(sqrt(S(2))*sqrt(S(1) - S(1)/cos(e + f*x))/S(2)), S(2)*b/(a + b))*tan(e + f*x)/(f*sqrt(-tan(e + f*x)**S(2))*sqrt(a + b/cos(e + f*x))*(c + d)), x)


def replacement4331(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((a + b/sin(e + f*x))/(a + b))*EllipticPi(S(2)*d/(c + d), asin(sqrt(S(2))*sqrt(S(1) - S(1)/sin(e + f*x))/S(2)), S(2)*b/(a + b))/(f*sqrt(-S(1)/tan(e + f*x)**S(2))*sqrt(a + b/sin(e + f*x))*(c + d)*tan(e + f*x)), x)


def replacement4332(a, b, c, d, e, f, g, x):
    return -Dist(a*g/(-a*d + b*c), Int(sqrt(g/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist(c*g/(-a*d + b*c), Int(sqrt(g/cos(e + f*x))*sqrt(a + b/cos(e + f*x))/(c + d/cos(e + f*x)), x), x)


def replacement4333(a, b, c, d, e, f, g, x):
    return -Dist(a*g/(-a*d + b*c), Int(sqrt(g/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist(c*g/(-a*d + b*c), Int(sqrt(g/sin(e + f*x))*sqrt(a + b/sin(e + f*x))/(c + d/sin(e + f*x)), x), x)


def replacement4334(a, b, c, d, e, f, g, x):
    return Dist(g*sqrt(g/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)/sqrt(a + b/cos(e + f*x)), Int(S(1)/(sqrt(a*cos(e + f*x) + b)*(c*cos(e + f*x) + d)), x), x)


def replacement4335(a, b, c, d, e, f, g, x):
    return Dist(g*sqrt(g/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)/sqrt(a + b/sin(e + f*x)), Int(S(1)/(sqrt(a*sin(e + f*x) + b)*(c*sin(e + f*x) + d)), x), x)


def replacement4336(a, b, c, d, e, f, x):
    return -Dist(a/(-a*d + b*c), Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Dist(c/(-a*d + b*c), Int(sqrt(a + b/cos(e + f*x))/((c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4337(a, b, c, d, e, f, x):
    return -Dist(a/(-a*d + b*c), Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Dist(c/(-a*d + b*c), Int(sqrt(a + b/sin(e + f*x))/((c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4338(a, b, c, d, e, f, x):
    return Dist(S(1)/d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) - Dist(c/d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4339(a, b, c, d, e, f, x):
    return Dist(S(1)/d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) - Dist(c/d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4340(a, b, c, d, e, f, g, x):
    return Dist(g**S(2)/(d*(-a*d + b*c)), Int(sqrt(g/cos(e + f*x))*(a*c + (-a*d + b*c)/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) - Dist(c**S(2)*g**S(2)/(d*(-a*d + b*c)), Int(sqrt(g/cos(e + f*x))*sqrt(a + b/cos(e + f*x))/(c + d/cos(e + f*x)), x), x)


def replacement4341(a, b, c, d, e, f, g, x):
    return Dist(g**S(2)/(d*(-a*d + b*c)), Int(sqrt(g/sin(e + f*x))*(a*c + (-a*d + b*c)/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) - Dist(c**S(2)*g**S(2)/(d*(-a*d + b*c)), Int(sqrt(g/sin(e + f*x))*sqrt(a + b/sin(e + f*x))/(c + d/sin(e + f*x)), x), x)


def replacement4342(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int((g/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x) - Dist(c*g/d, Int((g/cos(e + f*x))**(S(3)/2)/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))), x), x)


def replacement4343(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int((g/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x) - Dist(c*g/d, Int((g/sin(e + f*x))**(S(3)/2)/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))), x), x)


def replacement4344(a, b, c, d, e, f, x):
    return Dist(S(2)*b/f, Subst(Int(S(1)/(-b*d*x**S(2) + S(1)), x), x, tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x)))), x)


def replacement4345(a, b, c, d, e, f, x):
    return Dist(-S(2)*b/f, Subst(Int(S(1)/(-b*d*x**S(2) + S(1)), x), x, S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x))), x)


def replacement4346(a, b, c, d, e, f, x):
    return Dist(b/d, Int(sqrt(c + d/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4347(a, b, c, d, e, f, x):
    return Dist(b/d, Int(sqrt(c + d/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4348(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((S(1) + S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b/cos(e + f*x))*(c - d)))*sqrt(-(S(1) - S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b/cos(e + f*x))*(c + d)))*(a + b/cos(e + f*x))*EllipticPi(b*(c + d)/(d*(a + b)), asin(sqrt((a + b)/(c + d))*sqrt(c + d/cos(e + f*x))/sqrt(a + b/cos(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(d*f*sqrt((a + b)/(c + d))*tan(e + f*x)), x)


def replacement4349(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((S(1) + S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b/sin(e + f*x))*(c - d)))*sqrt(-(S(1) - S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b/sin(e + f*x))*(c + d)))*(a + b/sin(e + f*x))*EllipticPi(b*(c + d)/(d*(a + b)), asin(sqrt((a + b)/(c + d))*sqrt(c + d/sin(e + f*x))/sqrt(a + b/sin(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))*tan(e + f*x)/(d*f*sqrt((a + b)/(c + d))), x)


def replacement4350(a, b, c, d, e, f, x):
    return Dist(S(2)*a/(b*f), Subst(Int(S(1)/(x**S(2)*(a*c - b*d) + S(2)), x), x, tan(e + f*x)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x)))), x)


def replacement4351(a, b, c, d, e, f, x):
    return Dist(-S(2)*a/(b*f), Subst(Int(S(1)/(x**S(2)*(a*c - b*d) + S(2)), x), x, S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*tan(e + f*x))), x)


def replacement4352(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((S(1) - S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b)*(c + d/cos(e + f*x))))*sqrt(-(S(1) + S(1)/cos(e + f*x))*(-a*d + b*c)/((a - b)*(c + d/cos(e + f*x))))*(c + d/cos(e + f*x))*EllipticF(asin(sqrt(a + b/cos(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d/cos(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))/(f*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))*tan(e + f*x)), x)


def replacement4353(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((S(1) - S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d/sin(e + f*x))))*sqrt(-(S(1) + S(1)/sin(e + f*x))*(-a*d + b*c)/((a - b)*(c + d/sin(e + f*x))))*(c + d/sin(e + f*x))*EllipticF(asin(sqrt(a + b/sin(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d/sin(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))*tan(e + f*x)/(f*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))), x)


def replacement4354(a, b, c, d, e, f, x):
    return Dist(S(1)/b, Int(sqrt(a + b/cos(e + f*x))/(sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x) - Dist(a/b, Int(S(1)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4355(a, b, c, d, e, f, x):
    return Dist(S(1)/b, Int(sqrt(a + b/sin(e + f*x))/(sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x) - Dist(a/b, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4356(a, b, c, d, e, f, x):
    return Dist((a - b)/(c - d), Int(S(1)/(sqrt(a + b/cos(e + f*x))*sqrt(c + d/cos(e + f*x))*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/(c - d), Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*(c + d/cos(e + f*x))**(S(3)/2)*cos(e + f*x)), x), x)


def replacement4357(a, b, c, d, e, f, x):
    return Dist((a - b)/(c - d), Int(S(1)/(sqrt(a + b/sin(e + f*x))*sqrt(c + d/sin(e + f*x))*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/(c - d), Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*(c + d/sin(e + f*x))**(S(3)/2)*sin(e + f*x)), x), x)


def replacement4358(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(a**S(2)*g*tan(e + f*x)/(f*sqrt(a - b/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), Subst(Int((g*x)**(p + S(-1))*(a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, S(1)/cos(e + f*x)), x)


def replacement4359(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**S(2)*g/(f*sqrt(a - b/sin(e + f*x))*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), Subst(Int((g*x)**(p + S(-1))*(a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, S(1)/sin(e + f*x)), x)


def replacement4360(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(-m - n), Int((g/cos(e + f*x))**(m + n + p)*(a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement4361(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(-m - n), Int((g/sin(e + f*x))**(m + n + p)*(a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement4362(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(-m)*(g/cos(e + f*x))**(m + p)*(c + d/cos(e + f*x))**n*(c*cos(e + f*x) + d)**(-n), Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement4363(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(-m)*(g/sin(e + f*x))**(m + p)*(c + d/sin(e + f*x))**n*(c*sin(e + f*x) + d)**(-n), Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement4364(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/cos(e + f*x))**p*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n*(a*cos(e + f*x) + b)**(-m)*(c*cos(e + f*x) + d)**(-n), Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement4365(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/sin(e + f*x))**p*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n*(a*sin(e + f*x) + b)**(-m)*(c*sin(e + f*x) + d)**(-n), Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement4366(a, b, c, d, e, f, m, n, p, x):
    return Dist(sqrt(a + b/cos(e + f*x))*sqrt(c*cos(e + f*x) + d)/(sqrt(c + d/cos(e + f*x))*sqrt(a*cos(e + f*x) + b)), Int((a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-m - n - p), x), x)


def replacement4367(a, b, c, d, e, f, m, n, p, x):
    return Dist(sqrt(a + b/sin(e + f*x))*sqrt(c*sin(e + f*x) + d)/(sqrt(c + d/sin(e + f*x))*sqrt(a*sin(e + f*x) + b)), Int((a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-m - n - p), x), x)


def replacement4368(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g/cos(e + f*x))**p*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n, x), x)


def replacement4369(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g/sin(e + f*x))**p*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n, x), x)


def replacement4370(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g/cos(e + f*x))**p*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n, x)


def replacement4371(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g/sin(e + f*x))**p*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n, x)


def replacement4372(A, B, a, b, c, d, e, f, x):
    return Simp(S(2)*A*sqrt((S(1) - S(1)/cos(e + f*x))*(-a*d + b*c)/((a + b)*(c + d/cos(e + f*x))))*(S(1) + S(1)/cos(e + f*x))*EllipticE(asin(sqrt(a + b/cos(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d/cos(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))/(f*sqrt(-(S(1) + S(1)/cos(e + f*x))*(-a*d + b*c)/((a - b)*(c + d/cos(e + f*x))))*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))*tan(e + f*x)), x)


def replacement4373(A, B, a, b, c, d, e, f, x):
    return Simp(-S(2)*A*sqrt((S(1) - S(1)/sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d/sin(e + f*x))))*(S(1) + S(1)/sin(e + f*x))*EllipticE(asin(sqrt(a + b/sin(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d/sin(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))*tan(e + f*x)/(f*sqrt(-(S(1) + S(1)/sin(e + f*x))*(-a*d + b*c)/((a - b)*(c + d/sin(e + f*x))))*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))), x)


def replacement4374(A, B, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(n*(A*b + B*a) + (A*a*(n + S(1)) + B*b*n)/cos(e + f*x), x), x), x) - Simp(A*a*(d/cos(e + f*x))**n*tan(e + f*x)/(f*n), x)


def replacement4375(A, B, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(n*(A*b + B*a) + (A*a*(n + S(1)) + B*b*n)/sin(e + f*x), x), x), x) + Simp(A*a*(d/sin(e + f*x))**n/(f*n*tan(e + f*x)), x)


def replacement4376(A, B, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(1)), Int((d/cos(e + f*x))**n*Simp(A*a*(n + S(1)) + B*b*n + (n + S(1))*(A*b + B*a)/cos(e + f*x), x), x), x) + Simp(B*b*(d/cos(e + f*x))**n*tan(e + f*x)/(f*(n + S(1))), x)


def replacement4377(A, B, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(1)), Int((d/sin(e + f*x))**n*Simp(A*a*(n + S(1)) + B*b*n + (n + S(1))*(A*b + B*a)/sin(e + f*x), x), x), x) - Simp(B*b*(d/sin(e + f*x))**n/(f*(n + S(1))*tan(e + f*x)), x)


def replacement4378(A, B, a, b, e, f, x):
    return Dist(B/b, Int(S(1)/cos(e + f*x), x), x) + Dist((A*b - B*a)/b, Int(S(1)/((a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4379(A, B, a, b, e, f, x):
    return Dist(B/b, Int(S(1)/sin(e + f*x), x), x) + Dist((A*b - B*a)/b, Int(S(1)/((a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4380(A, B, a, b, e, f, m, x):
    return Simp(B*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4381(A, B, a, b, e, f, m, x):
    return -Simp(B*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4382(A, B, a, b, e, f, m, x):
    return Dist((A*b*(m + S(1)) + B*a*m)/(a*b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x) - Simp((a + b/cos(e + f*x))**m*(A*b - B*a)*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4383(A, B, a, b, e, f, m, x):
    return Dist((A*b*(m + S(1)) + B*a*m)/(a*b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x) + Simp((a + b/sin(e + f*x))**m*(A*b - B*a)/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4384(A, B, a, b, e, f, m, x):
    return Dist((A*b*(m + S(1)) + B*a*m)/(b*(m + S(1))), Int((a + b/cos(e + f*x))**m/cos(e + f*x), x), x) + Simp(B*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4385(A, B, a, b, e, f, m, x):
    return Dist((A*b*(m + S(1)) + B*a*m)/(b*(m + S(1))), Int((a + b/sin(e + f*x))**m/sin(e + f*x), x), x) - Simp(B*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4386(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + B*b*m + (A*b*(m + S(1)) + B*a*m)/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp(B*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4387(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + B*b*m + (A*b*(m + S(1)) + B*a*m)/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp(B*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4388(A, B, a, b, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp((m + S(1))*(A*a - B*b) - (m + S(2))*(A*b - B*a)/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**(m + S(1))*(A*b - B*a)*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4389(A, B, a, b, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp((m + S(1))*(A*a - B*b) - (m + S(2))*(A*b - B*a)/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**(m + S(1))*(A*b - B*a)/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4390(A, B, a, b, e, f, x):
    return Simp(S(2)*sqrt(b*(S(1) - S(1)/cos(e + f*x))/(a + b))*sqrt(-b*(S(1) + S(1)/cos(e + f*x))/(a - b))*(A*b - B*a)*EllipticE(asin(sqrt(a + b/cos(e + f*x))/Rt(a + B*b/A, S(2))), (A*a + B*b)/(A*a - B*b))*Rt(a + B*b/A, S(2))/(b**S(2)*f*tan(e + f*x)), x)


def replacement4391(A, B, a, b, e, f, x):
    return Simp(-S(2)*sqrt(b*(S(1) - S(1)/sin(e + f*x))/(a + b))*sqrt(-b*(S(1) + S(1)/sin(e + f*x))/(a - b))*(A*b - B*a)*EllipticE(asin(sqrt(a + b/sin(e + f*x))/Rt(a + B*b/A, S(2))), (A*a + B*b)/(A*a - B*b))*Rt(a + B*b/A, S(2))*tan(e + f*x)/(b**S(2)*f), x)


def replacement4392(A, B, a, b, e, f, x):
    return Dist(B, Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Dist(A - B, Int(S(1)/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x)


def replacement4393(A, B, a, b, e, f, x):
    return Dist(B, Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Dist(A - B, Int(S(1)/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x)


def replacement4394(A, B, a, b, e, f, m, x):
    return Simp(-S(2)*sqrt(S(2))*A*sqrt((A + B/cos(e + f*x))/A)*(A*(a + b/cos(e + f*x))/(A*a + B*b))**(-m)*(A - B/cos(e + f*x))*(a + b/cos(e + f*x))**m*AppellF1(S(1)/2, S(-1)/2, -m, S(3)/2, (A - B/cos(e + f*x))/(S(2)*A), b*(A - B/cos(e + f*x))/(A*b + B*a))/(B*f*tan(e + f*x)), x)


def replacement4395(A, B, a, b, e, f, m, x):
    return Simp(S(2)*sqrt(S(2))*A*sqrt((A + B/sin(e + f*x))/A)*(A*(a + b/sin(e + f*x))/(A*a + B*b))**(-m)*(A - B/sin(e + f*x))*(a + b/sin(e + f*x))**m*AppellF1(S(1)/2, S(-1)/2, -m, S(3)/2, (A - B/sin(e + f*x))/(S(2)*A), b*(A - B/sin(e + f*x))/(A*b + B*a))*tan(e + f*x)/(B*f), x)


def replacement4396(A, B, a, b, e, f, m, x):
    return Dist(B/b, Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x) + Dist((A*b - B*a)/b, Int((a + b/cos(e + f*x))**m/cos(e + f*x), x), x)


def replacement4397(A, B, a, b, e, f, m, x):
    return Dist(B/b, Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x) + Dist((A*b - B*a)/b, Int((a + b/sin(e + f*x))**m/sin(e + f*x), x), x)


def replacement4398(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(b**S(2)*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(A*b*m - B*a*m + B*b*(S(2)*m + S(1))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**m*(A*b - B*a)*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4399(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(b**S(2)*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(A*b*m - B*a*m + B*b*(S(2)*m + S(1))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**m*(A*b - B*a)/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4400(A, B, a, b, e, f, m, x):
    return -Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*b - B*a) - (A*a*b*(m + S(2)) - B*(a**S(2) + b**S(2)*(m + S(1))))/cos(e + f*x), x)/cos(e + f*x), x), x) - Simp(a*(a + b/cos(e + f*x))**(m + S(1))*(A*b - B*a)*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4401(A, B, a, b, e, f, m, x):
    return -Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*b - B*a) - (A*a*b*(m + S(2)) - B*(a**S(2) + b**S(2)*(m + S(1))))/sin(e + f*x), x)/sin(e + f*x), x), x) + Simp(a*(a + b/sin(e + f*x))**(m + S(1))*(A*b - B*a)/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4402(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/cos(e + f*x))**m*Simp(B*b*(m + S(1)) + (A*b*(m + S(2)) - B*a)/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp(B*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(2))), x)


def replacement4403(A, B, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/sin(e + f*x))**m*Simp(B*b*(m + S(1)) + (A*b*(m + S(2)) - B*a)/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp(B*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(2))*tan(e + f*x)), x)


def replacement4404(A, B, a, b, d, e, f, m, n, x):
    return -Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4405(A, B, a, b, d, e, f, m, n, x):
    return Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4406(A, B, a, b, d, e, f, m, n, x):
    return Dist((A*a*m + B*b*(m + S(1)))/(a**S(2)*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1)), x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*(A*b - B*a)*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4407(A, B, a, b, d, e, f, m, n, x):
    return Dist((A*a*m + B*b*(m + S(1)))/(a**S(2)*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1)), x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*(A*b - B*a)/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4408(A, B, a, b, d, e, f, m, n, x):
    return -Dist((A*a*m - B*b*n)/(b*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m, x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4409(A, B, a, b, d, e, f, m, n, x):
    return -Dist((A*a*m - B*b*n)/(b*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m, x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4410(A, B, a, b, d, e, f, n, x):
    return Simp(S(2)*B*b*(d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4411(A, B, a, b, d, e, f, n, x):
    return Simp(-S(2)*B*b*(d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4412(A, B, a, b, d, e, f, n, x):
    return Dist((A*b*(S(2)*n + S(1)) + S(2)*B*a*n)/(S(2)*a*d*n), Int((d/cos(e + f*x))**(n + S(1))*sqrt(a + b/cos(e + f*x)), x), x) - Simp(A*b**S(2)*(d/cos(e + f*x))**n*tan(e + f*x)/(a*f*n*sqrt(a + b/cos(e + f*x))), x)


def replacement4413(A, B, a, b, d, e, f, n, x):
    return Dist((A*b*(S(2)*n + S(1)) + S(2)*B*a*n)/(S(2)*a*d*n), Int((d/sin(e + f*x))**(n + S(1))*sqrt(a + b/sin(e + f*x)), x), x) + Simp(A*b**S(2)*(d/sin(e + f*x))**n/(a*f*n*sqrt(a + b/sin(e + f*x))*tan(e + f*x)), x)


def replacement4414(A, B, a, b, d, e, f, n, x):
    return Dist((A*b*(S(2)*n + S(1)) + S(2)*B*a*n)/(b*(S(2)*n + S(1))), Int((d/cos(e + f*x))**n*sqrt(a + b/cos(e + f*x)), x), x) + Simp(S(2)*B*b*(d/cos(e + f*x))**n*tan(e + f*x)/(f*sqrt(a + b/cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement4415(A, B, a, b, d, e, f, n, x):
    return Dist((A*b*(S(2)*n + S(1)) + S(2)*B*a*n)/(b*(S(2)*n + S(1))), Int((d/sin(e + f*x))**n*sqrt(a + b/sin(e + f*x)), x), x) + Simp(-S(2)*B*b*(d/sin(e + f*x))**n/(f*sqrt(a + b/sin(e + f*x))*(S(2)*n + S(1))*tan(e + f*x)), x)


def replacement4416(A, B, a, b, d, e, f, m, n, x):
    return -Dist(b/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m - n + S(-1)) - B*b*n - (A*b*(m + n) + B*a*n)/cos(e + f*x), x), x), x) - Simp(A*a*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*n), x)


def replacement4417(A, B, a, b, d, e, f, m, n, x):
    return -Dist(b/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m - n + S(-1)) - B*b*n - (A*b*(m + n) + B*a*n)/sin(e + f*x), x), x), x) + Simp(A*a*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))/(f*n*tan(e + f*x)), x)


def replacement4418(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*d*(m + n) + B*b*d*n + (A*b*d*(m + n) + B*a*d*(S(2)*m + n + S(-1)))/cos(e + f*x), x), x), x) + Simp(B*b*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*(m + n)), x)


def replacement4419(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*d*(m + n) + B*b*d*n + (A*b*d*(m + n) + B*a*d*(S(2)*m + n + S(-1)))/sin(e + f*x), x), x), x) - Simp(B*b*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))/(f*(m + n)*tan(e + f*x)), x)


def replacement4420(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*Simp(A*a*d*(n + S(-1)) - B*b*d*(n + S(-1)) - d*(A*b*(m + n) + B*a*(m - n + S(1)))/cos(e + f*x), x), x), x) - Simp(d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*(A*b - B*a)*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4421(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*Simp(A*a*d*(n + S(-1)) - B*b*d*(n + S(-1)) - d*(A*b*(m + n) + B*a*(m - n + S(1)))/sin(e + f*x), x), x), x) + Simp(d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m*(A*b - B*a)/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4422(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(-A*a*(S(2)*m + n + S(1)) + B*b*n + (A*b - B*a)*(m + n + S(1))/cos(e + f*x), x), x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*(A*b - B*a)*tan(e + f*x)/(b*f*(S(2)*m + S(1))), x)


def replacement4423(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(-A*a*(S(2)*m + n + S(1)) + B*b*n + (A*b - B*a)*(m + n + S(1))/sin(e + f*x), x), x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*(A*b - B*a)/(b*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4424(A, B, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n)), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*Simp(B*b*(n + S(-1)) + (A*b*(m + n) + B*a*m)/cos(e + f*x), x), x), x) + Simp(B*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n)), x)


def replacement4425(A, B, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n)), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m*Simp(B*b*(n + S(-1)) + (A*b*(m + n) + B*a*m)/sin(e + f*x), x), x), x) - Simp(B*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m/(f*(m + n)*tan(e + f*x)), x)


def replacement4426(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(A*a*m - A*b*(m + n + S(1))/cos(e + f*x) - B*b*n, x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4427(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(A*a*m - A*b*(m + n + S(1))/sin(e + f*x) - B*b*n, x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4428(A, B, a, b, d, e, f, m, n, x):
    return Dist(B/b, Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1)), x), x) + Dist((A*b - B*a)/b, Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m, x), x)


def replacement4429(A, B, a, b, d, e, f, m, n, x):
    return Dist(B/b, Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1)), x), x) + Dist((A*b - B*a)/b, Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m, x), x)


def replacement4430(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-2))*Simp(a*(-A*b*(m - n + S(-1)) + B*a*n) + b*(A*a*(m + n) + B*b*n)/cos(e + f*x)**S(2) + (A*(a**S(2)*(n + S(1)) + b**S(2)*n) + S(2)*B*a*b*n)/cos(e + f*x), x), x), x) - Simp(A*a*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*n), x)


def replacement4431(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-2))*Simp(a*(-A*b*(m - n + S(-1)) + B*a*n) + b*(A*a*(m + n) + B*b*n)/sin(e + f*x)**S(2) + (A*(a**S(2)*(n + S(1)) + b**S(2)*n) + S(2)*B*a*b*n)/sin(e + f*x), x), x), x) + Simp(A*a*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))/(f*n*tan(e + f*x)), x)


def replacement4432(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-2))*Simp(A*a**S(2)*(m + n) + B*a*b*n + b*(A*b*(m + n) + B*a*(S(2)*m + n + S(-1)))/cos(e + f*x)**S(2) + (B*b**S(2)*(m + n + S(-1)) + a*(m + n)*(S(2)*A*b + B*a))/cos(e + f*x), x), x), x) + Simp(B*b*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*tan(e + f*x)/(f*(m + n)), x)


def replacement4433(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-2))*Simp(A*a**S(2)*(m + n) + B*a*b*n + b*(A*b*(m + n) + B*a*(S(2)*m + n + S(-1)))/sin(e + f*x)**S(2) + (B*b**S(2)*(m + n + S(-1)) + a*(m + n)*(S(2)*A*b + B*a))/sin(e + f*x), x), x), x) - Simp(B*b*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))/(f*(m + n)*tan(e + f*x)), x)


def replacement4434(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*Simp(d*(m + S(1))*(A*a - B*b)/cos(e + f*x) + d*(n + S(-1))*(A*b - B*a) - d*(A*b - B*a)*(m + n + S(1))/cos(e + f*x)**S(2), x), x), x) + Simp(d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*(A*b - B*a)*tan(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4435(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*Simp(d*(m + S(1))*(A*a - B*b)/sin(e + f*x) + d*(n + S(-1))*(A*b - B*a) - d*(A*b - B*a)*(m + n + S(1))/sin(e + f*x)**S(2), x), x), x) - Simp(d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*(A*b - B*a)/(f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4436(A, B, a, b, d, e, f, m, n, x):
    return -Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*Simp(a*d*(n + S(-2))*(A*b - B*a) + b*d*(m + S(1))*(A*b - B*a)/cos(e + f*x) - (A*a*b*d*(m + n) - B*d*(a**S(2)*(n + S(-1)) + b**S(2)*(m + S(1))))/cos(e + f*x)**S(2), x), x), x) - Simp(a*d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*(A*b - B*a)*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4437(A, B, a, b, d, e, f, m, n, x):
    return -Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))*Simp(a*d*(n + S(-2))*(A*b - B*a) + b*d*(m + S(1))*(A*b - B*a)/sin(e + f*x) - (A*a*b*d*(m + n) - B*d*(a**S(2)*(n + S(-1)) + b**S(2)*(m + S(1))))/sin(e + f*x)**S(2), x), x), x) + Simp(a*d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))*(A*b - B*a)/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4438(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(A*(a**S(2)*(m + S(1)) - b**S(2)*(m + n + S(1))) + B*a*b*n - a*(m + S(1))*(A*b - B*a)/cos(e + f*x) + b*(A*b - B*a)*(m + n + S(2))/cos(e + f*x)**S(2), x), x), x) - Simp(b*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*(A*b - B*a)*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4439(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(A*(a**S(2)*(m + S(1)) - b**S(2)*(m + n + S(1))) + B*a*b*n - a*(m + S(1))*(A*b - B*a)/sin(e + f*x) + b*(A*b - B*a)*(m + n + S(2))/sin(e + f*x)**S(2), x), x), x) + Simp(b*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*(A*b - B*a)/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4440(A, B, a, b, d, e, f, m, n, x):
    return Dist(d/(m + n), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(B*a*(n + S(-1)) + (A*a*(m + n) + B*b*(m + n + S(-1)))/cos(e + f*x) + (A*b*(m + n) + B*a*m)/cos(e + f*x)**S(2), x), x), x) + Simp(B*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n)), x)


def replacement4441(A, B, a, b, d, e, f, m, n, x):
    return Dist(d/(m + n), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(B*a*(n + S(-1)) + (A*a*(m + n) + B*b*(m + n + S(-1)))/sin(e + f*x) + (A*b*(m + n) + B*a*m)/sin(e + f*x)**S(2), x), x), x) - Simp(B*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m/(f*(m + n)*tan(e + f*x)), x)


def replacement4442(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*b*m - A*b*(m + n + S(1))/cos(e + f*x)**S(2) - B*a*n - (A*a*(n + S(1)) + B*b*n)/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4443(A, B, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*b*m - A*b*(m + n + S(1))/sin(e + f*x)**S(2) - B*a*n - (A*a*(n + S(1)) + B*b*n)/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4444(A, B, a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n)), Int((d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**m*Simp(B*a*(n + S(-2)) + B*b*(m + n + S(-1))/cos(e + f*x) + (A*b*(m + n) - B*a*(n + S(-1)))/cos(e + f*x)**S(2), x), x), x) + Simp(B*d**S(2)*(d/cos(e + f*x))**(n + S(-2))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + n)), x)


def replacement4445(A, B, a, b, d, e, f, m, n, x):
    return Dist(d**S(2)/(b*(m + n)), Int((d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**m*Simp(B*a*(n + S(-2)) + B*b*(m + n + S(-1))/sin(e + f*x) + (A*b*(m + n) - B*a*(n + S(-1)))/sin(e + f*x)**S(2), x), x), x) - Simp(B*d**S(2)*(d/sin(e + f*x))**(n + S(-2))*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + n)*tan(e + f*x)), x)


def replacement4446(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(A*a*(n + S(1))/cos(e + f*x) - A*b*(m + n + S(1)) + A*b*(m + n + S(2))/cos(e + f*x)**S(2) + B*a*n, x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f*n), x)


def replacement4447(A, B, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(A*a*(n + S(1))/sin(e + f*x) - A*b*(m + n + S(1)) + A*b*(m + n + S(2))/sin(e + f*x)**S(2) + B*a*n, x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))/(a*f*n*tan(e + f*x)), x)


def replacement4448(A, B, a, b, d, e, f, x):
    return Dist(A/a, Int(sqrt(a + b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) - Dist((A*b - B*a)/(a*d), Int(sqrt(d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4449(A, B, a, b, d, e, f, x):
    return Dist(A/a, Int(sqrt(a + b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) - Dist((A*b - B*a)/(a*d), Int(sqrt(d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4450(A, B, a, b, d, e, f, x):
    return Dist(A, Int(sqrt(d/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x), x) + Dist(B/d, Int((d/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4451(A, B, a, b, d, e, f, x):
    return Dist(A, Int(sqrt(d/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x), x) + Dist(B/d, Int((d/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4452(A, B, a, b, d, e, f, x):
    return Dist(A, Int(sqrt(a + b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) + Dist(B/d, Int(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x)), x), x)


def replacement4453(A, B, a, b, d, e, f, x):
    return Dist(A, Int(sqrt(a + b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) + Dist(B/d, Int(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x)), x), x)


def replacement4454(A, B, a, b, d, e, f, n, x):
    return Dist(A/a, Int((d/cos(e + f*x))**n, x), x) - Dist((A*b - B*a)/(a*d), Int((d/cos(e + f*x))**(n + S(1))/(a + b/cos(e + f*x)), x), x)


def replacement4455(A, B, a, b, d, e, f, n, x):
    return Dist(A/a, Int((d/sin(e + f*x))**n, x), x) - Dist((A*b - B*a)/(a*d), Int((d/sin(e + f*x))**(n + S(1))/(a + b/sin(e + f*x)), x), x)


def replacement4456(A, B, a, b, d, e, f, m, n, x):
    return Int((d/cos(e + f*x))**n*(A + B/cos(e + f*x))*(a + b/cos(e + f*x))**m, x)


def replacement4457(A, B, a, b, d, e, f, m, n, x):
    return Int((d/sin(e + f*x))**n*(A + B/sin(e + f*x))*(a + b/sin(e + f*x))**m, x)


def replacement4458(A, B, a, b, c, d, e, f, m, n, p, x):
    return Dist((-a*c)**m, Int((A*cos(e + f*x) + B)**p*(c*cos(e + f*x) + d)**(-m + n)*sin(e + f*x)**(S(2)*m)*cos(e + f*x)**(-m - n - p), x), x)


def replacement4459(A, B, a, b, c, d, e, f, m, n, p, x):
    return Dist((-a*c)**m, Int((A*sin(e + f*x) + B)**p*(c*sin(e + f*x) + d)**(-m + n)*sin(e + f*x)**(-m - n - p)*cos(e + f*x)**(S(2)*m), x), x)


def replacement4460(A, B, C, a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b/cos(e + f*x), x), x), x)


def replacement4461(A, B, C, a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b/sin(e + f*x), x), x), x)


def replacement4462(A, C, a, b, e, f, m, x):
    return Dist(C/b**S(2), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-a + b/cos(e + f*x), x), x), x)


def replacement4463(A, C, a, b, e, f, m, x):
    return Dist(C/b**S(2), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-a + b/sin(e + f*x), x), x), x)


def replacement4464(A, C, b, e, f, m, x):
    return -Simp(A*(b/cos(e + f*x))**m*tan(e + f*x)/(f*m), x)


def replacement4465(A, C, b, e, f, m, x):
    return Simp(A*(b/sin(e + f*x))**m/(f*m*tan(e + f*x)), x)


def replacement4466(A, C, b, e, f, m, x):
    return Dist((A*(m + S(1)) + C*m)/(b**S(2)*m), Int((b/cos(e + f*x))**(m + S(2)), x), x) - Simp(A*(b/cos(e + f*x))**m*tan(e + f*x)/(f*m), x)


def replacement4467(A, C, b, e, f, m, x):
    return Dist((A*(m + S(1)) + C*m)/(b**S(2)*m), Int((b/sin(e + f*x))**(m + S(2)), x), x) + Simp(A*(b/sin(e + f*x))**m/(f*m*tan(e + f*x)), x)


def replacement4468(A, C, b, e, f, m, x):
    return Dist((A*(m + S(1)) + C*m)/(m + S(1)), Int((b/cos(e + f*x))**m, x), x) + Simp(C*(b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4469(A, C, b, e, f, m, x):
    return Dist((A*(m + S(1)) + C*m)/(m + S(1)), Int((b/sin(e + f*x))**m, x), x) - Simp(C*(b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4470(B, C, b, e, f, m, x):
    return Dist(B/b, Int((b/cos(e + f*x))**(m + S(1)), x), x) + Dist(C/b**S(2), Int((b/cos(e + f*x))**(m + S(2)), x), x)


def replacement4471(B, C, b, e, f, m, x):
    return Dist(B/b, Int((b/sin(e + f*x))**(m + S(1)), x), x) + Dist(C/b**S(2), Int((b/sin(e + f*x))**(m + S(2)), x), x)


def replacement4472(A, B, C, b, e, f, m, x):
    return Dist(B/b, Int((b/cos(e + f*x))**(m + S(1)), x), x) + Int((b/cos(e + f*x))**m*(A + C/cos(e + f*x)**S(2)), x)


def replacement4473(A, B, C, b, e, f, m, x):
    return Dist(B/b, Int((b/sin(e + f*x))**(m + S(1)), x), x) + Int((b/sin(e + f*x))**m*(A + C/sin(e + f*x)**S(2)), x)


def replacement4474(A, B, C, a, b, e, f, x):
    return Dist(S(1)/2, Int(Simp(S(2)*A*a + (S(2)*B*a + b*(S(2)*A + C))/cos(e + f*x) + S(2)*(B*b + C*a)/cos(e + f*x)**S(2), x), x), x) + Simp(C*b*tan(e + f*x)/(S(2)*f*cos(e + f*x)), x)


def replacement4475(A, B, C, a, b, e, f, x):
    return Dist(S(1)/2, Int(Simp(S(2)*A*a + (S(2)*B*a + b*(S(2)*A + C))/sin(e + f*x) + S(2)*(B*b + C*a)/sin(e + f*x)**S(2), x), x), x) - Simp(C*b/(S(2)*f*sin(e + f*x)*tan(e + f*x)), x)


def replacement4476(A, C, a, b, e, f, x):
    return Dist(S(1)/2, Int(Simp(S(2)*A*a + S(2)*C*a/cos(e + f*x)**S(2) + b*(S(2)*A + C)/cos(e + f*x), x), x), x) + Simp(C*b*tan(e + f*x)/(S(2)*f*cos(e + f*x)), x)


def replacement4477(A, C, a, b, e, f, x):
    return Dist(S(1)/2, Int(Simp(S(2)*A*a + S(2)*C*a/sin(e + f*x)**S(2) + b*(S(2)*A + C)/sin(e + f*x), x), x), x) - Simp(C*b/(S(2)*f*sin(e + f*x)*tan(e + f*x)), x)


def replacement4478(A, B, C, a, b, e, f, x):
    return Dist(S(1)/b, Int((A*b + (B*b - C*a)/cos(e + f*x))/(a + b/cos(e + f*x)), x), x) + Dist(C/b, Int(S(1)/cos(e + f*x), x), x)


def replacement4479(A, B, C, a, b, e, f, x):
    return Dist(S(1)/b, Int((A*b + (B*b - C*a)/sin(e + f*x))/(a + b/sin(e + f*x)), x), x) + Dist(C/b, Int(S(1)/sin(e + f*x), x), x)


def replacement4480(A, C, a, b, e, f, x):
    return Dist(S(1)/b, Int((A*b - C*a/cos(e + f*x))/(a + b/cos(e + f*x)), x), x) + Dist(C/b, Int(S(1)/cos(e + f*x), x), x)


def replacement4481(A, C, a, b, e, f, x):
    return Dist(S(1)/b, Int((A*b - C*a/sin(e + f*x))/(a + b/sin(e + f*x)), x), x) + Dist(C/b, Int(S(1)/sin(e + f*x), x), x)


def replacement4482(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + S(1)) + (B*b*(m + S(1)) - a*(A*(m + S(1)) - C*m))/cos(e + f*x), x), x), x) + Simp((a + b/cos(e + f*x))**m*(A*a - B*b + C*a)*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4483(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + S(1)) + (B*b*(m + S(1)) - a*(A*(m + S(1)) - C*m))/sin(e + f*x), x), x), x) - Simp((a + b/sin(e + f*x))**m*(A*a - B*b + C*a)/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4484(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + S(1)) - a*(A*(m + S(1)) - C*m)/cos(e + f*x), x), x), x) + Simp((A + C)*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement4485(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + S(1)) - a*(A*(m + S(1)) - C*m)/sin(e + f*x), x), x), x) - Simp((A + C)*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4486(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(1))), Int((a + b/cos(e + f*x))**m*Simp(A*b*(m + S(1)) + (B*b*(m + S(1)) + C*a*m)/cos(e + f*x), x), x), x) + Simp(C*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4487(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(1))), Int((a + b/sin(e + f*x))**m*Simp(A*b*(m + S(1)) + (B*b*(m + S(1)) + C*a*m)/sin(e + f*x), x), x), x) - Simp(C*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4488(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(1))), Int((a + b/cos(e + f*x))**m*Simp(A*b*(m + S(1)) + C*a*m/cos(e + f*x), x), x), x) + Simp(C*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4489(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(1))), Int((a + b/sin(e + f*x))**m*Simp(A*b*(m + S(1)) + C*a*m/sin(e + f*x), x), x), x) - Simp(C*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4490(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + (B*b*(m + S(1)) + C*a*m)/cos(e + f*x)**S(2) + (C*b*m + (m + S(1))*(A*b + B*a))/cos(e + f*x), x), x), x) + Simp(C*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4491(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + (B*b*(m + S(1)) + C*a*m)/sin(e + f*x)**S(2) + (C*b*m + (m + S(1))*(A*b + B*a))/sin(e + f*x), x), x), x) - Simp(C*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4492(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + C*a*m/cos(e + f*x)**S(2) + (A*b*(m + S(1)) + C*b*m)/cos(e + f*x), x), x), x) + Simp(C*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + S(1))), x)


def replacement4493(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m + S(1)) + C*a*m/sin(e + f*x)**S(2) + (A*b*(m + S(1)) + C*b*m)/sin(e + f*x), x), x), x) - Simp(C*(a + b/sin(e + f*x))**m/(f*(m + S(1))*tan(e + f*x)), x)


def replacement4494(A, B, C, a, b, e, f, x):
    return Dist(C, Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Int((A + (B - C)/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x)


def replacement4495(A, B, C, a, b, e, f, x):
    return Dist(C, Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Int((A + (B - C)/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x)


def replacement4496(A, C, a, b, e, f, x):
    return Dist(C, Int((S(1) + S(1)/cos(e + f*x))/(sqrt(a + b/cos(e + f*x))*cos(e + f*x)), x), x) + Int((A - C/cos(e + f*x))/sqrt(a + b/cos(e + f*x)), x)


def replacement4497(A, C, a, b, e, f, x):
    return Dist(C, Int((S(1) + S(1)/sin(e + f*x))/(sqrt(a + b/sin(e + f*x))*sin(e + f*x)), x), x) + Int((A - C/sin(e + f*x))/sqrt(a + b/sin(e + f*x)), x)


def replacement4498(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(A*(a**S(2) - b**S(2))*(m + S(1)) - a*(m + S(1))*(A*b - B*a + C*b)/cos(e + f*x) + (m + S(2))*(A*b**S(2) - B*a*b + C*a**S(2))/cos(e + f*x)**S(2), x), x), x) - Simp((a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4499(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(A*(a**S(2) - b**S(2))*(m + S(1)) - a*(m + S(1))*(A*b - B*a + C*b)/sin(e + f*x) + (m + S(2))*(A*b**S(2) - B*a*b + C*a**S(2))/sin(e + f*x)**S(2), x), x), x) + Simp((a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4500(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(A*(a**S(2) - b**S(2))*(m + S(1)) - a*b*(A + C)*(m + S(1))/cos(e + f*x) + (m + S(2))*(A*b**S(2) + C*a**S(2))/cos(e + f*x)**S(2), x), x), x) - Simp((a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4501(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(A*(a**S(2) - b**S(2))*(m + S(1)) - a*b*(A + C)*(m + S(1))/sin(e + f*x) + (m + S(2))*(A*b**S(2) + C*a**S(2))/sin(e + f*x)**S(2), x), x), x) + Simp((a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4502(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/cos(e + f*x))**m*(A*b + (B*b - C*a)/cos(e + f*x)), x), x) + Dist(C/b, Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x)


def replacement4503(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/sin(e + f*x))**m*(A*b + (B*b - C*a)/sin(e + f*x)), x), x) + Dist(C/b, Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x)


def replacement4504(A, C, a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/cos(e + f*x))**m*(A*b - C*a/cos(e + f*x)), x), x) + Dist(C/b, Int((a + b/cos(e + f*x))**(m + S(1))/cos(e + f*x), x), x)


def replacement4505(A, C, a, b, e, f, m, x):
    return Dist(S(1)/b, Int((a + b/sin(e + f*x))**m*(A*b - C*a/sin(e + f*x)), x), x) + Dist(C/b, Int((a + b/sin(e + f*x))**(m + S(1))/sin(e + f*x), x), x)


def replacement4506(A, B, C, b, e, f, m, x):
    return Dist(b**S(2), Int((b*cos(e + f*x))**(m + S(-2))*(A*cos(e + f*x)**S(2) + B*cos(e + f*x) + C), x), x)


def replacement4507(A, B, C, b, e, f, m, x):
    return Dist(b**S(2), Int((b*sin(e + f*x))**(m + S(-2))*(A*sin(e + f*x)**S(2) + B*sin(e + f*x) + C), x), x)


def replacement4508(A, C, b, e, f, m, x):
    return Dist(b**S(2), Int((b*cos(e + f*x))**(m + S(-2))*(A*cos(e + f*x)**S(2) + C), x), x)


def replacement4509(A, C, b, e, f, m, x):
    return Dist(b**S(2), Int((b*sin(e + f*x))**(m + S(-2))*(A*sin(e + f*x)**S(2) + C), x), x)


def replacement4510(A, B, C, a, b, e, f, m, p, x):
    return Dist(a**IntPart(m)*(a*(b/cos(e + f*x))**p)**FracPart(m)*(b/cos(e + f*x))**(-p*FracPart(m)), Int((b/cos(e + f*x))**(m*p)*(A + B/cos(e + f*x) + C/cos(e + f*x)**S(2)), x), x)


def replacement4511(A, B, C, a, b, e, f, m, p, x):
    return Dist(a**IntPart(m)*(a*(b/sin(e + f*x))**p)**FracPart(m)*(b/sin(e + f*x))**(-p*FracPart(m)), Int((b/sin(e + f*x))**(m*p)*(A + B/sin(e + f*x) + C/sin(e + f*x)**S(2)), x), x)


def replacement4512(A, C, a, b, e, f, m, p, x):
    return Dist(a**IntPart(m)*(a*(b/cos(e + f*x))**p)**FracPart(m)*(b/cos(e + f*x))**(-p*FracPart(m)), Int((b/cos(e + f*x))**(m*p)*(A + C/cos(e + f*x)**S(2)), x), x)


def replacement4513(A, C, a, b, e, f, m, p, x):
    return Dist(a**IntPart(m)*(a*(b/sin(e + f*x))**p)**FracPart(m)*(b/sin(e + f*x))**(-p*FracPart(m)), Int((b/sin(e + f*x))**(m*p)*(A + C/sin(e + f*x)**S(2)), x), x)


def replacement4514(A, B, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(C*b*n/cos(e + f*x)**S(2) + n*(A*b + B*a) + (A*a*(n + S(1)) + n*(B*b + C*a))/cos(e + f*x), x), x), x) - Simp(A*a*(d/cos(e + f*x))**n*tan(e + f*x)/(f*n), x)


def replacement4515(A, B, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(C*b*n/sin(e + f*x)**S(2) + n*(A*b + B*a) + (A*a*(n + S(1)) + n*(B*b + C*a))/sin(e + f*x), x), x), x) + Simp(A*a*(d/sin(e + f*x))**n/(f*n*tan(e + f*x)), x)


def replacement4516(A, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*Simp(A*b*n + C*b*n/cos(e + f*x)**S(2) + a*(A*(n + S(1)) + C*n)/cos(e + f*x), x), x), x) - Simp(A*a*(d/cos(e + f*x))**n*tan(e + f*x)/(f*n), x)


def replacement4517(A, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*Simp(A*b*n + C*b*n/sin(e + f*x)**S(2) + a*(A*(n + S(1)) + C*n)/sin(e + f*x), x), x), x) + Simp(A*a*(d/sin(e + f*x))**n/(f*n*tan(e + f*x)), x)


def replacement4518(A, B, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(2)), Int((d/cos(e + f*x))**n*Simp(A*a*(n + S(2)) + (n + S(2))*(B*b + C*a)/cos(e + f*x)**S(2) + (B*a*(n + S(2)) + b*(A*(n + S(2)) + C*(n + S(1))))/cos(e + f*x), x), x), x) + Simp(C*b*(d/cos(e + f*x))**n*tan(e + f*x)/(f*(n + S(2))*cos(e + f*x)), x)


def replacement4519(A, B, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(2)), Int((d/sin(e + f*x))**n*Simp(A*a*(n + S(2)) + (n + S(2))*(B*b + C*a)/sin(e + f*x)**S(2) + (B*a*(n + S(2)) + b*(A*(n + S(2)) + C*(n + S(1))))/sin(e + f*x), x), x), x) - Simp(C*b*(d/sin(e + f*x))**n/(f*(n + S(2))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4520(A, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(2)), Int((d/cos(e + f*x))**n*Simp(A*a*(n + S(2)) + C*a*(n + S(2))/cos(e + f*x)**S(2) + b*(A*(n + S(2)) + C*(n + S(1)))/cos(e + f*x), x), x), x) + Simp(C*b*(d/cos(e + f*x))**n*tan(e + f*x)/(f*(n + S(2))*cos(e + f*x)), x)


def replacement4521(A, C, a, b, d, e, f, n, x):
    return Dist(S(1)/(n + S(2)), Int((d/sin(e + f*x))**n*Simp(A*a*(n + S(2)) + C*a*(n + S(2))/sin(e + f*x)**S(2) + b*(A*(n + S(2)) + C*(n + S(1)))/sin(e + f*x), x), x), x) - Simp(C*b*(d/sin(e + f*x))**n/(f*(n + S(2))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4522(A, B, C, a, b, e, f, m, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-S(2)*A*b*(m + S(1)) + B*a - C*b - (B*b*(m + S(2)) - a*(A*(m + S(2)) - C*(m + S(-1))))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**m*(A*a - B*b + C*a)*tan(e + f*x)/(a*f*(S(2)*m + S(1))*cos(e + f*x)), x)


def replacement4523(A, B, C, a, b, e, f, m, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-S(2)*A*b*(m + S(1)) + B*a - C*b - (B*b*(m + S(2)) - a*(A*(m + S(2)) - C*(m + S(-1))))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**m*(A*a - B*b + C*a)/(a*f*(S(2)*m + S(1))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4524(A, C, a, b, e, f, m, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-S(2)*A*b*(m + S(1)) - C*b + a*(A*(m + S(2)) - C*(m + S(-1)))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((A + C)*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))*cos(e + f*x)), x)


def replacement4525(A, C, a, b, e, f, m, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-S(2)*A*b*(m + S(1)) - C*b + a*(A*(m + S(2)) - C*(m + S(-1)))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((A + C)*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4526(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*a - B*b + C*a) - (A*b**S(2) - B*a*b + C*a**S(2) + b*(m + S(1))*(A*b - B*a + C*b))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4527(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*a - B*b + C*a) - (A*b**S(2) - B*a*b + C*a**S(2) + b*(m + S(1))*(A*b - B*a + C*b))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4528(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(a*b*(A + C)*(m + S(1)) - (A*b**S(2) + C*a**S(2) + b*(m + S(1))*(A*b + C*b))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp((a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4529(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(a*b*(A + C)*(m + S(1)) - (A*b**S(2) + C*a**S(2) + b*(m + S(1))*(A*b + C*b))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp((a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4530(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/cos(e + f*x))**m*Simp(A*b*(m + S(2)) + C*b*(m + S(1)) + (B*b*(m + S(2)) - C*a)/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp(C*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(2))), x)


def replacement4531(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/sin(e + f*x))**m*Simp(A*b*(m + S(2)) + C*b*(m + S(1)) + (B*b*(m + S(2)) - C*a)/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp(C*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(2))*tan(e + f*x)), x)


def replacement4532(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/cos(e + f*x))**m*Simp(A*b*(m + S(2)) - C*a/cos(e + f*x) + C*b*(m + S(1)), x)/cos(e + f*x), x), x) + Simp(C*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(2))), x)


def replacement4533(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b/sin(e + f*x))**m*Simp(A*b*(m + S(2)) - C*a/sin(e + f*x) + C*b*(m + S(1)), x)/sin(e + f*x), x), x) - Simp(C*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(2))*tan(e + f*x)), x)


def replacement4534(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(-A*b*(S(2)*m + n + S(1)) + B*a*n - C*b*n - (B*b*(m + n + S(1)) - a*(A*(m + n + S(1)) - C*(m - n)))/cos(e + f*x), x), x), x) + Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*(A*a - B*b + C*a)*tan(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement4535(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(-A*b*(S(2)*m + n + S(1)) + B*a*n - C*b*n - (B*b*(m + n + S(1)) - a*(A*(m + n + S(1)) - C*(m - n)))/sin(e + f*x), x), x), x) - Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*(A*a - B*b + C*a)/(a*f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4536(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + n + S(1)) + C*b*n - a*(A*(m + n + S(1)) - C*(m - n))/cos(e + f*x), x), x), x) + Simp((d/cos(e + f*x))**n*(A + C)*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(S(2)*m + S(1))), x)


def replacement4537(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(A*b*(S(2)*m + n + S(1)) + C*b*n - a*(A*(m + n + S(1)) - C*(m - n))/sin(e + f*x), x), x), x) - Simp((d/sin(e + f*x))**n*(A + C)*(a + b/sin(e + f*x))**m/(f*(S(2)*m + S(1))*tan(e + f*x)), x)


def replacement4538(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(A*a*m - B*b*n - b*(A*(m + n + S(1)) + C*n)/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4539(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(A*a*m - B*b*n - b*(A*(m + n + S(1)) + C*n)/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4540(A, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(A*a*m - b*(A*(m + n + S(1)) + C*n)/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4541(A, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(A*a*m - b*(A*(m + n + S(1)) + C*n)/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4542(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*Simp(A*b*(m + n + S(1)) + C*b*n + (B*b*(m + n + S(1)) + C*a*m)/cos(e + f*x), x), x), x) + Simp(C*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(1))), x)


def replacement4543(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*Simp(A*b*(m + n + S(1)) + C*b*n + (B*b*(m + n + S(1)) + C*a*m)/sin(e + f*x), x), x), x) - Simp(C*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4544(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*Simp(A*b*(m + n + S(1)) + C*a*m/cos(e + f*x) + C*b*n, x), x), x) + Simp(C*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(1))), x)


def replacement4545(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*Simp(A*b*(m + n + S(1)) + C*a*m/sin(e + f*x) + C*b*n, x), x), x) - Simp(C*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4546(A, B, C, a, b, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-C*b*(a**S(2) - b**S(2))*(m + S(1))/cos(e + f*x)**S(2) + b*(m + S(1))*(A*b**S(2) - a*(B*b - C*a)) + (B*b*(a**S(2) + b**S(2)*(m + S(1))) - a*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1)))))/cos(e + f*x), x)/cos(e + f*x), x), x) - Simp(a*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*tan(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4547(A, B, C, a, b, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-C*b*(a**S(2) - b**S(2))*(m + S(1))/sin(e + f*x)**S(2) + b*(m + S(1))*(A*b**S(2) - a*(B*b - C*a)) + (B*b*(a**S(2) + b**S(2)*(m + S(1))) - a*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1)))))/sin(e + f*x), x)/sin(e + f*x), x), x) + Simp(a*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4548(A, C, a, b, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/cos(e + f*x))**(m + S(1))*Simp(-C*b*(a**S(2) - b**S(2))*(m + S(1))/cos(e + f*x)**S(2) - a*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1))))/cos(e + f*x) + b*(m + S(1))*(A*b**S(2) + C*a**S(2)), x)/cos(e + f*x), x), x) - Simp(a*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*tan(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4549(A, C, a, b, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b/sin(e + f*x))**(m + S(1))*Simp(-C*b*(a**S(2) - b**S(2))*(m + S(1))/sin(e + f*x)**S(2) - a*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1))))/sin(e + f*x) + b*(m + S(1))*(A*b**S(2) + C*a**S(2)), x)/sin(e + f*x), x), x) + Simp(a*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4550(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b/cos(e + f*x))**m*Simp(C*a + b*(A*(m + S(3)) + C*(m + S(2)))/cos(e + f*x) - (-B*b*(m + S(3)) + S(2)*C*a)/cos(e + f*x)**S(2), x)/cos(e + f*x), x), x) + Simp(C*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(3))*cos(e + f*x)), x)


def replacement4551(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b/sin(e + f*x))**m*Simp(C*a + b*(A*(m + S(3)) + C*(m + S(2)))/sin(e + f*x) - (-B*b*(m + S(3)) + S(2)*C*a)/sin(e + f*x)**S(2), x)/sin(e + f*x), x), x) - Simp(C*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(3))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4552(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b/cos(e + f*x))**m*Simp(C*a - S(2)*C*a/cos(e + f*x)**S(2) + b*(A*(m + S(3)) + C*(m + S(2)))/cos(e + f*x), x)/cos(e + f*x), x), x) + Simp(C*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + S(3))*cos(e + f*x)), x)


def replacement4553(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b/sin(e + f*x))**m*Simp(C*a - S(2)*C*a/sin(e + f*x)**S(2) + b*(A*(m + S(3)) + C*(m + S(2)))/sin(e + f*x), x)/sin(e + f*x), x), x) - Simp(C*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + S(3))*sin(e + f*x)*tan(e + f*x)), x)


def replacement4554(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*b*m - B*a*n - b*(A*(m + n + S(1)) + C*n)/cos(e + f*x)**S(2) - (B*b*n + a*(A*(n + S(1)) + C*n))/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4555(A, B, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*b*m - B*a*n - b*(A*(m + n + S(1)) + C*n)/sin(e + f*x)**S(2) - (B*b*n + a*(A*(n + S(1)) + C*n))/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4556(A, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*b*m - a*(A*(n + S(1)) + C*n)/cos(e + f*x) - b*(A*(m + n + S(1)) + C*n)/cos(e + f*x)**S(2), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*n), x)


def replacement4557(A, C, a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*b*m - a*(A*(n + S(1)) + C*n)/sin(e + f*x) - b*(A*(m + n + S(1)) + C*n)/sin(e + f*x)**S(2), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*n*tan(e + f*x)), x)


def replacement4558(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n + S(1)), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m + n + S(1)) + C*a*n + (B*b*(m + n + S(1)) + C*a*m)/cos(e + f*x)**S(2) + (C*b*(m + n) + (A*b + B*a)*(m + n + S(1)))/cos(e + f*x), x), x), x) + Simp(C*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(1))), x)


def replacement4559(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n + S(1)), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m + n + S(1)) + C*a*n + (B*b*(m + n + S(1)) + C*a*m)/sin(e + f*x)**S(2) + (C*b*(m + n) + (A*b + B*a)*(m + n + S(1)))/sin(e + f*x), x), x), x) - Simp(C*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4560(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n + S(1)), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(-1))*Simp(A*a*(m + n + S(1)) + C*a*m/cos(e + f*x)**S(2) + C*a*n + b*(A*(m + n + S(1)) + C*(m + n))/cos(e + f*x), x), x), x) + Simp(C*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*tan(e + f*x)/(f*(m + n + S(1))), x)


def replacement4561(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(m + n + S(1)), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(-1))*Simp(A*a*(m + n + S(1)) + C*a*m/sin(e + f*x)**S(2) + C*a*n + b*(A*(m + n + S(1)) + C*(m + n))/sin(e + f*x), x), x), x) - Simp(C*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m/(f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4562(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*Simp(A*b**S(2)*(n + S(-1)) - a*(n + S(-1))*(B*b - C*a) + b*(m + S(1))*(A*a - B*b + C*a)/cos(e + f*x) - (C*(a**S(2)*n + b**S(2)*(m + S(1))) + b*(A*b - B*a)*(m + n + S(1)))/cos(e + f*x)**S(2), x), x), x) + Simp(d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4563(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*Simp(A*b**S(2)*(n + S(-1)) - a*(n + S(-1))*(B*b - C*a) + b*(m + S(1))*(A*a - B*b + C*a)/sin(e + f*x) - (C*(a**S(2)*n + b**S(2)*(m + S(1))) + b*(A*b - B*a)*(m + n + S(1)))/sin(e + f*x)**S(2), x), x), x) - Simp(d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4564(A, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*Simp(A*b**S(2)*(n + S(-1)) + C*a**S(2)*(n + S(-1)) + a*b*(A + C)*(m + S(1))/cos(e + f*x) - (A*b**S(2)*(m + n + S(1)) + C*(a**S(2)*n + b**S(2)*(m + S(1))))/cos(e + f*x)**S(2), x), x), x) + Simp(d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*tan(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4565(A, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*Simp(A*b**S(2)*(n + S(-1)) + C*a**S(2)*(n + S(-1)) + a*b*(A + C)*(m + S(1))/sin(e + f*x) - (A*b**S(2)*(m + n + S(1)) + C*(a**S(2)*n + b**S(2)*(m + S(1))))/sin(e + f*x)**S(2), x), x), x) - Simp(d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(b*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4566(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(a*(m + S(1))*(A*a - B*b + C*a) - a*(m + S(1))*(A*b - B*a + C*b)/cos(e + f*x) - (m + n + S(1))*(A*b**S(2) - B*a*b + C*a**S(2)) + (m + n + S(2))*(A*b**S(2) - B*a*b + C*a**S(2))/cos(e + f*x)**S(2), x), x), x) - Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4567(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(a*(m + S(1))*(A*a - B*b + C*a) - a*(m + S(1))*(A*b - B*a + C*b)/sin(e + f*x) - (m + n + S(1))*(A*b**S(2) - B*a*b + C*a**S(2)) + (m + n + S(2))*(A*b**S(2) - B*a*b + C*a**S(2))/sin(e + f*x)**S(2), x), x), x) + Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4568(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*Simp(a**S(2)*(A + C)*(m + S(1)) - a*b*(A + C)*(m + S(1))/cos(e + f*x) - (A*b**S(2) + C*a**S(2))*(m + n + S(1)) + (A*b**S(2) + C*a**S(2))*(m + n + S(2))/cos(e + f*x)**S(2), x), x), x) - Simp((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*tan(e + f*x)/(a*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement4569(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*(a**S(2) - b**S(2))*(m + S(1))), Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*Simp(a**S(2)*(A + C)*(m + S(1)) - a*b*(A + C)*(m + S(1))/sin(e + f*x) - (A*b**S(2) + C*a**S(2))*(m + n + S(1)) + (A*b**S(2) + C*a**S(2))*(m + n + S(2))/sin(e + f*x)**S(2), x), x), x) + Simp((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(a*f*(a**S(2) - b**S(2))*(m + S(1))*tan(e + f*x)), x)


def replacement4570(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*Simp(C*a*(n + S(-1)) + (A*b*(m + n + S(1)) + C*b*(m + n))/cos(e + f*x) + (B*b*(m + n + S(1)) - C*a*n)/cos(e + f*x)**S(2), x), x), x) + Simp(C*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + n + S(1))), x)


def replacement4571(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m*Simp(C*a*(n + S(-1)) + (A*b*(m + n + S(1)) + C*b*(m + n))/sin(e + f*x) + (B*b*(m + n + S(1)) - C*a*n)/sin(e + f*x)**S(2), x), x), x) - Simp(C*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4572(A, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n + S(1))), Int((d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**m*Simp(-C*a*n/cos(e + f*x)**S(2) + C*a*(n + S(-1)) + (A*b*(m + n + S(1)) + C*b*(m + n))/cos(e + f*x), x), x), x) + Simp(C*d*(d/cos(e + f*x))**(n + S(-1))*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(b*f*(m + n + S(1))), x)


def replacement4573(A, C, a, b, d, e, f, m, n, x):
    return Dist(d/(b*(m + n + S(1))), Int((d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**m*Simp(-C*a*n/sin(e + f*x)**S(2) + C*a*(n + S(-1)) + (A*b*(m + n + S(1)) + C*b*(m + n))/sin(e + f*x), x), x), x) - Simp(C*d*(d/sin(e + f*x))**(n + S(-1))*(a + b/sin(e + f*x))**(m + S(1))/(b*f*(m + n + S(1))*tan(e + f*x)), x)


def replacement4574(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(-A*b*(m + n + S(1)) + A*b*(m + n + S(2))/cos(e + f*x)**S(2) + B*a*n + a*(A*n + A + C*n)/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f*n), x)


def replacement4575(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(-A*b*(m + n + S(1)) + A*b*(m + n + S(2))/sin(e + f*x)**S(2) + B*a*n + a*(A*n + A + C*n)/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))/(a*f*n*tan(e + f*x)), x)


def replacement4576(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/cos(e + f*x))**(n + S(1))*(a + b/cos(e + f*x))**m*Simp(-A*b*(m + n + S(1)) + A*b*(m + n + S(2))/cos(e + f*x)**S(2) + a*(A*n + A + C*n)/cos(e + f*x), x), x), x) - Simp(A*(d/cos(e + f*x))**n*(a + b/cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f*n), x)


def replacement4577(A, C, a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*d*n), Int((d/sin(e + f*x))**(n + S(1))*(a + b/sin(e + f*x))**m*Simp(-A*b*(m + n + S(1)) + A*b*(m + n + S(2))/sin(e + f*x)**S(2) + a*(A*n + A + C*n)/sin(e + f*x), x), x), x) + Simp(A*(d/sin(e + f*x))**n*(a + b/sin(e + f*x))**(m + S(1))/(a*f*n*tan(e + f*x)), x)


def replacement4578(A, B, C, a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((A*a - (A*b - B*a)/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) + Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2)*d**S(2)), Int((d/cos(e + f*x))**(S(3)/2)/(a + b/cos(e + f*x)), x), x)


def replacement4579(A, B, C, a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((A*a - (A*b - B*a)/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) + Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2)*d**S(2)), Int((d/sin(e + f*x))**(S(3)/2)/(a + b/sin(e + f*x)), x), x)


def replacement4580(A, C, a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((A*a - A*b/cos(e + f*x))/sqrt(d/cos(e + f*x)), x), x) + Dist((A*b**S(2) + C*a**S(2))/(a**S(2)*d**S(2)), Int((d/cos(e + f*x))**(S(3)/2)/(a + b/cos(e + f*x)), x), x)


def replacement4581(A, C, a, b, d, e, f, x):
    return Dist(a**(S(-2)), Int((A*a - A*b/sin(e + f*x))/sqrt(d/sin(e + f*x)), x), x) + Dist((A*b**S(2) + C*a**S(2))/(a**S(2)*d**S(2)), Int((d/sin(e + f*x))**(S(3)/2)/(a + b/sin(e + f*x)), x), x)


def replacement4582(A, B, C, a, b, d, e, f, x):
    return Dist(C/d**S(2), Int((d/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x) + Int((A + B/cos(e + f*x))/(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), x)


def replacement4583(A, B, C, a, b, d, e, f, x):
    return Dist(C/d**S(2), Int((d/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x) + Int((A + B/sin(e + f*x))/(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x))), x)


def replacement4584(A, C, a, b, d, e, f, x):
    return Dist(A, Int(S(1)/(sqrt(d/cos(e + f*x))*sqrt(a + b/cos(e + f*x))), x), x) + Dist(C/d**S(2), Int((d/cos(e + f*x))**(S(3)/2)/sqrt(a + b/cos(e + f*x)), x), x)


def replacement4585(A, C, a, b, d, e, f, x):
    return Dist(A, Int(S(1)/(sqrt(d/sin(e + f*x))*sqrt(a + b/sin(e + f*x))), x), x) + Dist(C/d**S(2), Int((d/sin(e + f*x))**(S(3)/2)/sqrt(a + b/sin(e + f*x)), x), x)


def replacement4586(A, B, C, a, b, d, e, f, m, n, x):
    return Int((d/cos(e + f*x))**n*(a + b/cos(e + f*x))**m*(A + B/cos(e + f*x) + C/cos(e + f*x)**S(2)), x)


def replacement4587(A, B, C, a, b, d, e, f, m, n, x):
    return Int((d/sin(e + f*x))**n*(a + b/sin(e + f*x))**m*(A + B/sin(e + f*x) + C/sin(e + f*x)**S(2)), x)


def replacement4588(A, C, a, b, d, e, f, m, n, x):
    return Int((d/cos(e + f*x))**n*(A + C/cos(e + f*x)**S(2))*(a + b/cos(e + f*x))**m, x)


def replacement4589(A, C, a, b, d, e, f, m, n, x):
    return Int((d/sin(e + f*x))**n*(A + C/sin(e + f*x)**S(2))*(a + b/sin(e + f*x))**m, x)


def replacement4590(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d**(m + S(2)), Int((d*cos(e + f*x))**(-m + n + S(-2))*(a*cos(e + f*x) + b)**m*(A*cos(e + f*x)**S(2) + B*cos(e + f*x) + C), x), x)


def replacement4591(A, B, C, a, b, d, e, f, m, n, x):
    return Dist(d**(m + S(2)), Int((d*sin(e + f*x))**(-m + n + S(-2))*(a*sin(e + f*x) + b)**m*(A*sin(e + f*x)**S(2) + B*sin(e + f*x) + C), x), x)


def replacement4592(A, C, a, b, d, e, f, m, n, x):
    return Dist(d**(m + S(2)), Int((d*cos(e + f*x))**(-m + n + S(-2))*(A*cos(e + f*x)**S(2) + C)*(a*cos(e + f*x) + b)**m, x), x)


def replacement4593(A, C, a, b, d, e, f, m, n, x):
    return Dist(d**(m + S(2)), Int((d*sin(e + f*x))**(-m + n + S(-2))*(A*sin(e + f*x)**S(2) + C)*(a*sin(e + f*x) + b)**m, x), x)


def replacement4594(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/cos(e + f*x))**p)**FracPart(n)*(d/cos(e + f*x))**(-p*FracPart(n)), Int((d/cos(e + f*x))**(n*p)*(a + b/cos(e + f*x))**m*(A + B/cos(e + f*x) + C/cos(e + f*x)**S(2)), x), x)


def replacement4595(A, B, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/sin(e + f*x))**p)**FracPart(n)*(d/sin(e + f*x))**(-p*FracPart(n)), Int((d/sin(e + f*x))**(n*p)*(a + b/sin(e + f*x))**m*(A + B/sin(e + f*x) + C/sin(e + f*x)**S(2)), x), x)


def replacement4596(A, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/cos(e + f*x))**p)**FracPart(n)*(d/cos(e + f*x))**(-p*FracPart(n)), Int((d/cos(e + f*x))**(n*p)*(A + C/cos(e + f*x)**S(2))*(a + b/cos(e + f*x))**m, x), x)


def replacement4597(A, C, a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d/sin(e + f*x))**p)**FracPart(n)*(d/sin(e + f*x))**(-p*FracPart(n)), Int((d/sin(e + f*x))**(n*p)*(A + C/sin(e + f*x)**S(2))*(a + b/sin(e + f*x))**m, x), x)


def replacement4598(b, c, d, n, x):
    return Dist(b/d, Subst(Int((b*x**S(2) + b)**(n + S(-1)), x), x, tan(c + d*x)), x)


def replacement4599(b, c, d, n, x):
    return -Dist(b/d, Subst(Int((b*x**S(2) + b)**(n + S(-1)), x), x, S(1)/tan(c + d*x)), x)


def replacement4600(a, b, c, d, p, x):
    return Int((-a*tan(c + d*x)**S(2))**p, x)


def replacement4601(a, b, c, d, p, x):
    return Int((-a/tan(c + d*x)**S(2))**p, x)


def replacement4602(a, b, c, d, x):
    return -Dist(b/a, Int(S(1)/(a*cos(c + d*x)**S(2) + b), x), x) + Simp(x/a, x)


def replacement4603(a, b, c, d, x):
    return -Dist(b/a, Int(S(1)/(a*sin(c + d*x)**S(2) + b), x), x) + Simp(x/a, x)


def replacement4604(a, b, c, d, p, x):
    return Dist(S(1)/d, Subst(Int((a + b*x**S(2) + b)**p/(x**S(2) + S(1)), x), x, tan(c + d*x)), x)


def replacement4605(a, b, c, d, p, x):
    return -Dist(S(1)/d, Subst(Int((a + b*x**S(2) + b)**p/(x**S(2) + S(1)), x), x, S(1)/tan(c + d*x)), x)


def With4606(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f**(m + S(1))/d, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, tan(c + d*x)/f), x)


def With4607(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f**(m + S(1))/d, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def With4608(a, b, c, d, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f/d, Subst(Int((f*x)**(-n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, cos(c + d*x)/f), x)


def With4609(a, b, c, d, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f/d, Subst(Int((f*x)**(-n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, sin(c + d*x)/f), x)


def With4610(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, tan(c + d*x)/f), x)


def With4611(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def With4612(a, b, c, d, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1)/2)*ExpandToSum(a*(-f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, sin(c + d*x)/f), x)


def With4613(a, b, c, d, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1)/2)*ExpandToSum(a*(-f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, cos(c + d*x)/f), x)


def replacement4614(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((a + b*(S(1)/cos(c + d*x))**n)**p*(S(1)/cos(c + d*x))**m, x), x)


def replacement4615(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((a + b*(S(1)/sin(c + d*x))**n)**p*(S(1)/sin(c + d*x))**m, x), x)


def With4616(a, b, c, d, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f**(-m - n*p + S(1))/d, Subst(Int(x**(-m - n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, cos(c + d*x)/f), x)


def With4617(a, b, c, d, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f**(-m - n*p + S(1))/d, Subst(Int(x**(-m - n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, sin(c + d*x)/f), x)


def With4618(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f**(m + S(1))/d, Subst(Int(x**m*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p/(f**S(2)*x**S(2) + S(1)), x), x, tan(c + d*x)/f), x)


def With4619(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f**(m + S(1))/d, Subst(Int(x**m*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p/(f**S(2)*x**S(2) + S(1)), x), x, S(1)/(f*tan(c + d*x))), x)


def replacement4620(a, b, c, d, e, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/cos(d + e*x))**n)**(S(2)*p), x), x)


def replacement4621(a, b, c, d, e, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/sin(d + e*x))**n)**(S(2)*p), x), x)


def replacement4622(a, b, c, d, e, n, n2, p, x):
    return Dist((b + S(2)*c*(S(1)/cos(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/cos(d + e*x))**n + c*(S(1)/cos(d + e*x))**(S(2)*n))**p, Int(u*(b + S(2)*c*(S(1)/cos(d + e*x))**n)**(S(2)*p), x), x)


def replacement4623(a, b, c, d, e, n, n2, p, x):
    return Dist((b + S(2)*c*(S(1)/sin(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/sin(d + e*x))**n + c*(S(1)/sin(d + e*x))**(S(2)*n))**p, Int(u*(b + S(2)*c*(S(1)/sin(d + e*x))**n)**(S(2)*p), x), x)


def With4624(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/cos(d + e*x))**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/cos(d + e*x))**n + q), x), x)


def With4625(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/sin(d + e*x))**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/sin(d + e*x))**n + q), x), x)


def With4626(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(cos(d + e*x), x)
    return -Dist(f/e, Subst(Int((f*x)**(-n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, cos(d + e*x)/f), x)


def With4627(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(sin(d + e*x), x)
    return Dist(f/e, Subst(Int((f*x)**(-n*p)*(a*(f*x)**n + b)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, sin(d + e*x)/f), x)


def With4628(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(f**S(2)*x**S(2) + S(1))**n, x)**p, x), x, tan(d + e*x)/f), x)


def With4629(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1))*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(f**S(2)*x**S(2) + S(1))**n, x)**p, x), x, S(1)/(f*tan(d + e*x))), x)


def replacement4630(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/cos(d + e*x))**n)**(S(2)*p)*(S(1)/cos(d + e*x))**m, x), x)


def replacement4631(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/sin(d + e*x))**n)**(S(2)*p)*(S(1)/sin(d + e*x))**m, x), x)


def replacement4632(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*(S(1)/cos(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/cos(d + e*x))**n + c*(S(1)/cos(d + e*x))**(S(2)*n))**p, Int((b + S(2)*c*(S(1)/cos(d + e*x))**n)**(S(2)*p)*(S(1)/cos(d + e*x))**m, x), x)


def replacement4633(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*(S(1)/sin(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/sin(d + e*x))**n + c*(S(1)/sin(d + e*x))**(S(2)*n))**p, Int((b + S(2)*c*(S(1)/sin(d + e*x))**n)**(S(2)*p)*(S(1)/sin(d + e*x))**m, x), x)


def replacement4634(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((a + b*(S(1)/cos(d + e*x))**n + c*(S(1)/cos(d + e*x))**(S(2)*n))**p*(S(1)/cos(d + e*x))**m, x), x)


def replacement4635(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((a + b*(S(1)/sin(d + e*x))**n + c*(S(1)/sin(d + e*x))**(S(2)*n))**p*(S(1)/sin(d + e*x))**m, x), x)


def With4636(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(cos(d + e*x), x)
    return -Dist(f**(-m - n*p + S(1))/e, Subst(Int(x**(-m - S(2)*n*p)*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(b*(f*x)**n + c*(f*x)**(S(2)*n) + c)**p, x), x, cos(d + e*x)/f), x)


def With4637(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(sin(d + e*x), x)
    return Dist(f**(-m - n*p + S(1))/e, Subst(Int(x**(-m - S(2)*n*p)*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(b*(f*x)**n + c*(f*x)**(S(2)*n) + c)**p, x), x, sin(d + e*x)/f), x)


def With4638(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f**(m + S(1))/e, Subst(Int(x**m*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(f**S(2)*x**S(2) + S(1))**n, x)**p/(f**S(2)*x**S(2) + S(1)), x), x, tan(d + e*x)/f), x)


def With4639(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f**(m + S(1))/e, Subst(Int(x**m*ExpandToSum(a + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(f**S(2)*x**S(2) + S(1))**n, x)**p/(f**S(2)*x**S(2) + S(1)), x), x, S(1)/(f*tan(d + e*x))), x)


def replacement4640(A, B, a, b, c, d, e, n, x):
    return Dist(S(4)**(-n)*c**(-n), Int((A + B/cos(d + e*x))*(b + S(2)*c/cos(d + e*x))**(S(2)*n), x), x)


def replacement4641(A, B, a, b, c, d, e, n, x):
    return Dist(S(4)**(-n)*c**(-n), Int((A + B/sin(d + e*x))*(b + S(2)*c/sin(d + e*x))**(S(2)*n), x), x)


def replacement4642(A, B, a, b, c, d, e, n, x):
    return Dist((b + S(2)*c/cos(d + e*x))**(-S(2)*n)*(a + b/cos(d + e*x) + c/cos(d + e*x)**S(2))**n, Int((A + B/cos(d + e*x))*(b + S(2)*c/cos(d + e*x))**(S(2)*n), x), x)


def replacement4643(A, B, a, b, c, d, e, n, x):
    return Dist((b + S(2)*c/sin(d + e*x))**(-S(2)*n)*(a + b/sin(d + e*x) + c/sin(d + e*x)**S(2))**n, Int((A + B/sin(d + e*x))*(b + S(2)*c/sin(d + e*x))**(S(2)*n), x), x)


def With4644(A, B, a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c/cos(d + e*x) - q), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c/cos(d + e*x) + q), x), x)


def With4645(A, B, a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c/sin(d + e*x) - q), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c/sin(d + e*x) + q), x), x)


def replacement4646(A, B, a, b, c, d, e, n, x):
    return Int(ExpandTrig((A + B/cos(d + e*x))*(a + b/cos(d + e*x) + c/cos(d + e*x)**S(2))**n, x), x)


def replacement4647(A, B, a, b, c, d, e, n, x):
    return Int(ExpandTrig((A + B/sin(d + e*x))*(a + b/sin(d + e*x) + c/sin(d + e*x)**S(2))**n, x), x)


def replacement4648(c, d, e, f, m, x):
    return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*log(-I*exp(I*(e + f*x)) + S(1)), x), x) + Dist(d*m/f, Int((c + d*x)**(m + S(-1))*log(I*exp(I*(e + f*x)) + S(1)), x), x) + Simp(-S(2)*I*(c + d*x)**m*ArcTan(exp(I*e + I*f*x))/f, x)


def replacement4649(c, d, e, f, m, x):
    return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*log(S(1) - exp(I*(e + f*x))), x), x) + Dist(d*m/f, Int((c + d*x)**(m + S(-1))*log(exp(I*(e + f*x)) + S(1)), x), x) + Simp(-S(2)*(c + d*x)**m*atanh(exp(I*e + I*f*x))/f, x)


def replacement4650(c, d, e, f, m, x):
    return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*tan(e + f*x), x), x) + Simp((c + d*x)**m*tan(e + f*x)/f, x)


def replacement4651(c, d, e, f, m, x):
    return Dist(d*m/f, Int((c + d*x)**(m + S(-1))/tan(e + f*x), x), x) - Simp((c + d*x)**m/(f*tan(e + f*x)), x)


def replacement4652(b, c, d, e, f, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/cos(e + f*x))**(n + S(-2))*(c + d*x), x), x) - Simp(b**S(2)*d*(b/cos(e + f*x))**(n + S(-2))/(f**S(2)*(n + S(-2))*(n + S(-1))), x) + Simp(b**S(2)*(b/cos(e + f*x))**(n + S(-2))*(c + d*x)*tan(e + f*x)/(f*(n + S(-1))), x)


def replacement4653(b, c, d, e, f, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/sin(e + f*x))**(n + S(-2))*(c + d*x), x), x) - Simp(b**S(2)*d*(b/sin(e + f*x))**(n + S(-2))/(f**S(2)*(n + S(-2))*(n + S(-1))), x) - Simp(b**S(2)*(b/sin(e + f*x))**(n + S(-2))*(c + d*x)/(f*(n + S(-1))*tan(e + f*x)), x)


def replacement4654(b, c, d, e, f, m, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/cos(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) + Dist(b**S(2)*d**S(2)*m*(m + S(-1))/(f**S(2)*(n + S(-2))*(n + S(-1))), Int((b/cos(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(-2)), x), x) + Simp(b**S(2)*(b/cos(e + f*x))**(n + S(-2))*(c + d*x)**m*tan(e + f*x)/(f*(n + S(-1))), x) - Simp(b**S(2)*d*m*(b/cos(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(-1))/(f**S(2)*(n + S(-2))*(n + S(-1))), x)


def replacement4655(b, c, d, e, f, m, n, x):
    return Dist(b**S(2)*(n + S(-2))/(n + S(-1)), Int((b/sin(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) + Dist(b**S(2)*d**S(2)*m*(m + S(-1))/(f**S(2)*(n + S(-2))*(n + S(-1))), Int((b/sin(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(-2)), x), x) - Simp(b**S(2)*(b/sin(e + f*x))**(n + S(-2))*(c + d*x)**m/(f*(n + S(-1))*tan(e + f*x)), x) - Simp(b**S(2)*d*m*(b/sin(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(-1))/(f**S(2)*(n + S(-2))*(n + S(-1))), x)


def replacement4656(b, c, d, e, f, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/cos(e + f*x))**(n + S(2))*(c + d*x), x), x) + Simp(d*(b/cos(e + f*x))**n/(f**S(2)*n**S(2)), x) - Simp((b/cos(e + f*x))**(n + S(1))*(c + d*x)*sin(e + f*x)/(b*f*n), x)


def replacement4657(b, c, d, e, f, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/sin(e + f*x))**(n + S(2))*(c + d*x), x), x) + Simp(d*(b/sin(e + f*x))**n/(f**S(2)*n**S(2)), x) + Simp((b/sin(e + f*x))**(n + S(1))*(c + d*x)*cos(e + f*x)/(b*f*n), x)


def replacement4658(b, c, d, e, f, m, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/cos(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) - Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b/cos(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) - Simp((b/cos(e + f*x))**(n + S(1))*(c + d*x)**m*sin(e + f*x)/(b*f*n), x) + Simp(d*m*(b/cos(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)


def replacement4659(b, c, d, e, f, m, n, x):
    return Dist((n + S(1))/(b**S(2)*n), Int((b/sin(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) - Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b/sin(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) + Simp((b/sin(e + f*x))**(n + S(1))*(c + d*x)**m*cos(e + f*x)/(b*f*n), x) + Simp(d*m*(b/sin(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)


def replacement4660(b, c, d, e, f, m, n, x):
    return Dist((b/cos(e + f*x))**n*(b*cos(e + f*x))**n, Int((b*cos(e + f*x))**(-n)*(c + d*x)**m, x), x)


def replacement4661(b, c, d, e, f, m, n, x):
    return Dist((b/sin(e + f*x))**n*(b*sin(e + f*x))**n, Int((b*sin(e + f*x))**(-n)*(c + d*x)**m, x), x)


def replacement4662(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a + b/cos(e + f*x))**n, x), x)


def replacement4663(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a + b/sin(e + f*x))**n, x), x)


def replacement4664(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a*cos(e + f*x) + b)**n*cos(e + f*x)**(-n), x), x)


def replacement4665(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a*sin(e + f*x) + b)**n*sin(e + f*x)**(-n), x), x)


def replacement4666(a, b, m, n, u, v, x):
    return Int((a + b/cos(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)


def replacement4667(a, b, m, n, u, v, x):
    return Int((a + b/sin(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)


def replacement4668(a, b, c, d, e, f, m, n, x):
    return Int((a + b/cos(e + f*x))**n*(c + d*x)**m, x)


def replacement4669(a, b, c, d, e, f, m, n, x):
    return Int((a + b/sin(e + f*x))**n*(c + d*x)**m, x)


def replacement4670(a, b, c, d, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/cos(c + d*x))**p, x), x, x**n), x)


def replacement4671(a, b, c, d, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/sin(c + d*x))**p, x), x, x**n), x)


def replacement4672(a, b, c, d, n, p, x):
    return Int((a + b/cos(c + d*x**n))**p, x)


def replacement4673(a, b, c, d, n, p, x):
    return Int((a + b/sin(c + d*x**n))**p, x)


def replacement4674(a, b, c, d, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/cos(c + d*x**n))**p, x), x, u), x)


def replacement4675(a, b, c, d, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/sin(c + d*x**n))**p, x), x, u), x)


def replacement4676(a, b, p, u, x):
    return Int((a + b/cos(ExpandToSum(u, x)))**p, x)


def replacement4677(a, b, p, u, x):
    return Int((a + b/sin(ExpandToSum(u, x)))**p, x)


def replacement4678(a, b, c, d, m, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/cos(c + d*x))**p, x), x, x**n), x)


def replacement4679(a, b, c, d, m, n, p, x):
    return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/sin(c + d*x))**p, x), x, x**n), x)


def replacement4680(a, b, c, d, m, n, p, x):
    return Int(x**m*(a + b/cos(c + d*x**n))**p, x)


def replacement4681(a, b, c, d, m, n, p, x):
    return Int(x**m*(a + b/sin(c + d*x**n))**p, x)


def replacement4682(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/cos(c + d*x**n))**p, x), x)


def replacement4683(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/sin(c + d*x**n))**p, x), x)


def replacement4684(a, b, e, m, p, u, x):
    return Int((e*x)**m*(a + b/cos(ExpandToSum(u, x)))**p, x)


def replacement4685(a, b, e, m, p, u, x):
    return Int((e*x)**m*(a + b/sin(ExpandToSum(u, x)))**p, x)


def replacement4686(a, b, m, n, p, x):
    return -Dist((m - n + S(1))/(b*n*(p + S(-1))), Int(x**(m - n)*(S(1)/cos(a + b*x**n))**(p + S(-1)), x), x) + Simp(x**(m - n + S(1))*(S(1)/cos(a + b*x**n))**(p + S(-1))/(b*n*(p + S(-1))), x)


def replacement4687(a, b, m, n, p, x):
    return Dist((m - n + S(1))/(b*n*(p + S(-1))), Int(x**(m - n)*(S(1)/sin(a + b*x**n))**(p + S(-1)), x), x) - Simp(x**(m - n + S(1))*(S(1)/sin(a + b*x**n))**(p + S(-1))/(b*n*(p + S(-1))), x)
