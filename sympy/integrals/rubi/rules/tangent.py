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
    i, ii , Pqq, Q, R, r, C, k, u = symbols('i ii Pqq Q R r C k u')
    _UseGamma = False
    ShowSteps = False
    StepCounter = None

def tangent(rubi):
    from sympy.integrals.rubi.constraints import cons1510, cons2, cons3, cons50, cons127, cons19, cons4, cons1511, cons20, cons25, cons95, cons167, cons96, cons1512, cons1172, cons89, cons1513, cons91, cons168, cons685, cons33, cons268, cons34, cons1514, cons1515, cons21, cons157, cons1252, cons1516, cons1517, cons1518, cons1519, cons1520, cons1521, cons1522, cons1523, cons1524, cons82, cons81, cons1525, cons1526, cons1527, cons1528, cons29, cons5, cons1529, cons8, cons1263, cons1266, cons1441, cons465, cons1442, cons1530, cons1257, cons1531, cons90, cons1532, cons1533, cons1534, cons517, cons1535, cons1536, cons1537, cons1538, cons1539, cons1540, cons68, cons1541, cons86, cons1230, cons150, cons198, cons1542, cons87, cons72, cons1543, cons73, cons1414, cons269, cons1305, cons170, cons1544, cons1545, cons1324, cons1546, cons1325, cons1547, cons1548, cons1549, cons1550, cons79, cons1551, cons1552, cons1553, cons1425, cons1387, cons113, cons274, cons1327, cons1329, cons1336, cons1554, cons1555, cons1335, cons1556, cons1557, cons1338, cons1558, cons1559, cons810, cons382, cons210, cons149, cons1419, cons52, cons40, cons36, cons37, cons348, cons1420, cons1428, cons1560, cons1561, cons1562, cons1258, cons1563, cons38, cons35, cons1435, cons1564, cons1565, cons1566, cons1567, cons1568, cons1569, cons1570, cons1571, cons1494, cons1458, cons1572, cons1483, cons1481, cons745, cons1499, cons48, cons47, cons228, cons1482, cons64, cons530, cons1573, cons1574, cons812, cons813, cons1362, cons1575, cons1497, cons70, cons71, cons825, cons826, cons1576, cons1577, cons1578, cons1579, cons1580, cons1581, cons49, cons241, cons1582

    pattern3401 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1510)
    rule3401 = ReplacementRule(pattern3401, replacement3401)
    pattern3402 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1510)
    rule3402 = ReplacementRule(pattern3402, replacement3402)
    pattern3403 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons50, cons127, cons1511)
    rule3403 = ReplacementRule(pattern3403, replacement3403)
    pattern3404 = Pattern(Integral((S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons50, cons127, cons1511)
    rule3404 = ReplacementRule(pattern3404, replacement3404)
    pattern3405 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons3, cons50, cons127, cons4, cons20, cons25)
    rule3405 = ReplacementRule(pattern3405, replacement3405)
    pattern3406 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons3, cons50, cons127, cons4, cons20, cons25)
    rule3406 = ReplacementRule(pattern3406, replacement3406)
    pattern3407 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons167, cons96, cons1512, cons1172)
    rule3407 = ReplacementRule(pattern3407, replacement3407)
    pattern3408 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons167, cons96, cons1512, cons1172)
    rule3408 = ReplacementRule(pattern3408, replacement3408)
    pattern3409 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons167, cons1172, cons1513)
    rule3409 = ReplacementRule(pattern3409, replacement3409)
    pattern3410 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons167, cons1172, cons1513)
    rule3410 = ReplacementRule(pattern3410, replacement3410)
    pattern3411 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons91, cons168, cons1172)
    rule3411 = ReplacementRule(pattern3411, replacement3411)
    pattern3412 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons91, cons168, cons1172)
    rule3412 = ReplacementRule(pattern3412, replacement3412)
    pattern3413 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons91, cons685, cons1172)
    rule3413 = ReplacementRule(pattern3413, replacement3413)
    pattern3414 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons91, cons685, cons1172)
    rule3414 = ReplacementRule(pattern3414, replacement3414)
    pattern3415 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons268, cons1172)
    rule3415 = ReplacementRule(pattern3415, replacement3415)
    pattern3416 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons268, cons1172)
    rule3416 = ReplacementRule(pattern3416, replacement3416)
    pattern3417 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons34, cons685, cons1172)
    rule3417 = ReplacementRule(pattern3417, replacement3417)
    pattern3418 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons4, cons33, cons34, cons685, cons1172)
    rule3418 = ReplacementRule(pattern3418, replacement3418)
    pattern3419 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1514)
    rule3419 = ReplacementRule(pattern3419, replacement3419)
    pattern3420 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1514)
    rule3420 = ReplacementRule(pattern3420, replacement3420)
    pattern3421 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1515)
    rule3421 = ReplacementRule(pattern3421, replacement3421)
    pattern3422 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1515)
    rule3422 = ReplacementRule(pattern3422, replacement3422)
    pattern3423 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule3423 = ReplacementRule(pattern3423, replacement3423)
    pattern3424 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule3424 = ReplacementRule(pattern3424, replacement3424)
    pattern3425 = Pattern(Integral((WC('a', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule3425 = ReplacementRule(pattern3425, replacement3425)
    pattern3426 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons157)
    rule3426 = ReplacementRule(pattern3426, replacement3426)
    pattern3427 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons157)
    rule3427 = ReplacementRule(pattern3427, replacement3427)
    pattern3428 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1252, cons1516)
    rule3428 = ReplacementRule(pattern3428, replacement3428)
    pattern3429 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1252, cons1516)
    rule3429 = ReplacementRule(pattern3429, replacement3429)
    pattern3430 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons3, cons50, cons127, cons4, cons1517, cons1518)
    rule3430 = ReplacementRule(pattern3430, replacement3430)
    pattern3431 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons3, cons50, cons127, cons4, cons1517, cons1518)
    rule3431 = ReplacementRule(pattern3431, replacement3431)
    pattern3432 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons91, cons1519, cons1172)
    rule3432 = ReplacementRule(pattern3432, replacement3432)
    pattern3433 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons91, cons1519, cons1172)
    rule3433 = ReplacementRule(pattern3433, replacement3433)
    pattern3434 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons91, cons1172)
    rule3434 = ReplacementRule(pattern3434, replacement3434)
    pattern3435 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons91, cons1172)
    rule3435 = ReplacementRule(pattern3435, replacement3435)
    pattern3436 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons167, cons1520, cons1172)
    rule3436 = ReplacementRule(pattern3436, replacement3436)
    pattern3437 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons167, cons1520, cons1172)
    rule3437 = ReplacementRule(pattern3437, replacement3437)
    pattern3438 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons167, cons1521, cons1172)
    rule3438 = ReplacementRule(pattern3438, replacement3438)
    pattern3439 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons89, cons167, cons1521, cons1172)
    rule3439 = ReplacementRule(pattern3439, replacement3439)
    pattern3440 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons1522, cons1172)
    rule3440 = ReplacementRule(pattern3440, replacement3440)
    pattern3441 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons1522, cons1172)
    rule3441 = ReplacementRule(pattern3441, replacement3441)
    pattern3442 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons1523, cons1521, cons1172)
    rule3442 = ReplacementRule(pattern3442, replacement3442)
    pattern3443 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons1523, cons1521, cons1172)
    rule3443 = ReplacementRule(pattern3443, replacement3443)
    pattern3444 = Pattern(Integral(S(1)/(sqrt(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons50, cons127, cons1524)
    rule3444 = ReplacementRule(pattern3444, replacement3444)
    pattern3445 = Pattern(Integral(S(1)/(sqrt(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons50, cons127, cons1524)
    rule3445 = ReplacementRule(pattern3445, replacement3445)
    pattern3446 = Pattern(Integral(sqrt(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons3, cons50, cons127, cons1524)
    rule3446 = ReplacementRule(pattern3446, replacement3446)
    pattern3447 = Pattern(Integral(sqrt(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons3, cons50, cons127, cons1524)
    rule3447 = ReplacementRule(pattern3447, replacement3447)
    pattern3448 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons82, cons81)
    rule3448 = ReplacementRule(pattern3448, replacement3448)
    pattern3449 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons82, cons81)
    rule3449 = ReplacementRule(pattern3449, replacement3449)
    pattern3450 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1525, cons1526)
    rule3450 = ReplacementRule(pattern3450, replacement3450)
    pattern3451 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1525, cons1526)
    rule3451 = ReplacementRule(pattern3451, replacement3451)
    pattern3452 = Pattern(Integral((WC('a', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule3452 = ReplacementRule(pattern3452, replacement3452)
    pattern3453 = Pattern(Integral((WC('a', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule3453 = ReplacementRule(pattern3453, replacement3453)
    pattern3454 = Pattern(Integral(((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons89, cons167, cons1527, cons1528)
    rule3454 = ReplacementRule(pattern3454, replacement3454)
    pattern3455 = Pattern(Integral(((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons89, cons167, cons1527, cons1528)
    rule3455 = ReplacementRule(pattern3455, replacement3455)
    pattern3456 = Pattern(Integral(((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**WC('m', S(1))*(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons89, cons91, cons1529, cons1528)
    rule3456 = ReplacementRule(pattern3456, replacement3456)
    pattern3457 = Pattern(Integral(((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('a', S(1)))**WC('m', S(1))*(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons5, cons89, cons91, cons1529, cons1528)
    rule3457 = ReplacementRule(pattern3457, replacement3457)
    pattern3458 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons89, cons167)
    rule3458 = ReplacementRule(pattern3458, replacement3458)
    pattern3459 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons89, cons167)
    rule3459 = ReplacementRule(pattern3459, replacement3459)
    pattern3460 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons89, cons91)
    rule3460 = ReplacementRule(pattern3460, replacement3460)
    pattern3461 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons89, cons91)
    rule3461 = ReplacementRule(pattern3461, replacement3461)
    pattern3462 = Pattern(Integral(tan(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons1263)
    rule3462 = ReplacementRule(pattern3462, replacement3462)
    pattern3463 = Pattern(Integral(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons1263)
    rule3463 = ReplacementRule(pattern3463, replacement3463)
    pattern3464 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons25)
    rule3464 = ReplacementRule(pattern3464, replacement3464)
    pattern3465 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons25)
    rule3465 = ReplacementRule(pattern3465, replacement3465)
    pattern3466 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule3466 = ReplacementRule(pattern3466, replacement3466)
    pattern3467 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule3467 = ReplacementRule(pattern3467, replacement3467)
    pattern3468 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1441, cons89, cons167)
    rule3468 = ReplacementRule(pattern3468, replacement3468)
    pattern3469 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1441, cons89, cons167)
    rule3469 = ReplacementRule(pattern3469, replacement3469)
    pattern3470 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1441, cons89, cons465)
    rule3470 = ReplacementRule(pattern3470, replacement3470)
    pattern3471 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1441, cons89, cons465)
    rule3471 = ReplacementRule(pattern3471, replacement3471)
    pattern3472 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1441)
    rule3472 = ReplacementRule(pattern3472, replacement3472)
    pattern3473 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1441)
    rule3473 = ReplacementRule(pattern3473, replacement3473)
    pattern3474 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1441)
    rule3474 = ReplacementRule(pattern3474, replacement3474)
    pattern3475 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1441)
    rule3475 = ReplacementRule(pattern3475, replacement3475)
    pattern3476 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons167)
    rule3476 = ReplacementRule(pattern3476, replacement3476)
    pattern3477 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons167)
    rule3477 = ReplacementRule(pattern3477, replacement3477)
    pattern3478 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons91)
    rule3478 = ReplacementRule(pattern3478, replacement3478)
    pattern3479 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons91)
    rule3479 = ReplacementRule(pattern3479, replacement3479)
    pattern3480 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3480 = ReplacementRule(pattern3480, replacement3480)
    pattern3481 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3481 = ReplacementRule(pattern3481, replacement3481)
    pattern3482 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1442)
    rule3482 = ReplacementRule(pattern3482, replacement3482)
    pattern3483 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1442)
    rule3483 = ReplacementRule(pattern3483, replacement3483)
    pattern3484 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1530)
    rule3484 = ReplacementRule(pattern3484, replacement3484)
    pattern3485 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1530)
    rule3485 = ReplacementRule(pattern3485, replacement3485)
    pattern3486 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1441, cons1517)
    rule3486 = ReplacementRule(pattern3486, replacement3486)
    pattern3487 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1441, cons1517)
    rule3487 = ReplacementRule(pattern3487, replacement3487)
    pattern3488 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1257)
    rule3488 = ReplacementRule(pattern3488, replacement3488)
    pattern3489 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1257)
    rule3489 = ReplacementRule(pattern3489, replacement3489)
    pattern3490 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1441)
    rule3490 = ReplacementRule(pattern3490, replacement3490)
    pattern3491 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1441)
    rule3491 = ReplacementRule(pattern3491, replacement3491)
    pattern3492 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons1531, cons95, cons90)
    rule3492 = ReplacementRule(pattern3492, replacement3492)
    pattern3493 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons1531, cons95, cons90)
    rule3493 = ReplacementRule(pattern3493, replacement3493)
    pattern3494 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons1531, cons95, cons91)
    rule3494 = ReplacementRule(pattern3494, replacement3494)
    pattern3495 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons1531, cons95, cons91)
    rule3495 = ReplacementRule(pattern3495, replacement3495)
    pattern3496 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1531)
    rule3496 = ReplacementRule(pattern3496, replacement3496)
    pattern3497 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1531)
    rule3497 = ReplacementRule(pattern3497, replacement3497)
    pattern3498 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1532)
    rule3498 = ReplacementRule(pattern3498, replacement3498)
    pattern3499 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1532)
    rule3499 = ReplacementRule(pattern3499, replacement3499)
    pattern3500 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1533, cons25)
    rule3500 = ReplacementRule(pattern3500, replacement3500)
    pattern3501 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1533, cons25)
    rule3501 = ReplacementRule(pattern3501, replacement3501)
    pattern3502 = Pattern(Integral(sqrt(WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3502 = ReplacementRule(pattern3502, replacement3502)
    pattern3503 = Pattern(Integral(sqrt(WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3503 = ReplacementRule(pattern3503, replacement3503)
    pattern3504 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons167, cons1534, cons517)
    rule3504 = ReplacementRule(pattern3504, replacement3504)
    pattern3505 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons167, cons1535, cons517)
    rule3505 = ReplacementRule(pattern3505, replacement3505)
    pattern3506 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons90, cons96, cons1172)
    rule3506 = ReplacementRule(pattern3506, replacement3506)
    pattern3507 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons90, cons96, cons1172)
    rule3507 = ReplacementRule(pattern3507, replacement3507)
    pattern3508 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons90, cons1521, cons1172)
    rule3508 = ReplacementRule(pattern3508, replacement3508)
    pattern3509 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons90, cons1521, cons1172)
    rule3509 = ReplacementRule(pattern3509, replacement3509)
    pattern3510 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3510 = ReplacementRule(pattern3510, replacement3510)
    pattern3511 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3511 = ReplacementRule(pattern3511, replacement3511)
    pattern3512 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons91, cons1536, cons517)
    rule3512 = ReplacementRule(pattern3512, replacement3512)
    pattern3513 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons91, cons1536, cons517)
    rule3513 = ReplacementRule(pattern3513, replacement3513)
    pattern3514 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons465, cons168, cons1537, cons1521, cons1172)
    rule3514 = ReplacementRule(pattern3514, replacement3514)
    pattern3515 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons1441, cons95, cons465, cons168, cons1537, cons1521, cons1172)
    rule3515 = ReplacementRule(pattern3515, replacement3515)
    pattern3516 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons465, cons1538, cons1172)
    rule3516 = ReplacementRule(pattern3516, replacement3516)
    pattern3517 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1441, cons89, cons465, cons1538, cons1172)
    rule3517 = ReplacementRule(pattern3517, replacement3517)
    pattern3518 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1539, cons89)
    rule3518 = ReplacementRule(pattern3518, replacement3518)
    pattern3519 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1539, cons89)
    rule3519 = ReplacementRule(pattern3519, replacement3519)
    pattern3520 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1540, cons1538)
    rule3520 = ReplacementRule(pattern3520, replacement3520)
    pattern3521 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441, cons1540, cons1538)
    rule3521 = ReplacementRule(pattern3521, replacement3521)
    pattern3522 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441)
    rule3522 = ReplacementRule(pattern3522, replacement3522)
    pattern3523 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1441)
    rule3523 = ReplacementRule(pattern3523, replacement3523)
    pattern3524 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1442, cons1517)
    rule3524 = ReplacementRule(pattern3524, replacement3524)
    pattern3525 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1442, cons1517)
    rule3525 = ReplacementRule(pattern3525, replacement3525)
    pattern3526 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1442)
    rule3526 = ReplacementRule(pattern3526, replacement3526)
    pattern3527 = Pattern(Integral((a_ + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)*tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons50, cons127, cons1442)
    rule3527 = ReplacementRule(pattern3527, replacement3527)
    pattern3528 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1442, cons68)
    rule3528 = ReplacementRule(pattern3528, replacement3528)
    pattern3529 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1442, cons68)
    rule3529 = ReplacementRule(pattern3529, replacement3529)
    pattern3530 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1442)
    rule3530 = ReplacementRule(pattern3530, replacement3530)
    pattern3531 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1442)
    rule3531 = ReplacementRule(pattern3531, replacement3531)
    pattern3532 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1442, cons1541)
    rule3532 = ReplacementRule(pattern3532, replacement3532)
    pattern3533 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1442, cons1541)
    rule3533 = ReplacementRule(pattern3533, replacement3533)
    pattern3534 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1442, cons86)
    rule3534 = ReplacementRule(pattern3534, replacement3534)
    pattern3535 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1442, cons86)
    rule3535 = ReplacementRule(pattern3535, replacement3535)
    pattern3536 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1442, cons1526)
    rule3536 = ReplacementRule(pattern3536, replacement3536)
    pattern3537 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1442, cons1526)
    rule3537 = ReplacementRule(pattern3537, replacement3537)
    pattern3538 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3538 = ReplacementRule(pattern3538, replacement3538)
    pattern3539 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3539 = ReplacementRule(pattern3539, replacement3539)
    pattern3540 = Pattern(Integral(S(1)/((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3540 = ReplacementRule(pattern3540, replacement3540)
    pattern3541 = Pattern(Integral(S(1)/((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1441)
    rule3541 = ReplacementRule(pattern3541, replacement3541)
    pattern3542 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons21)
    rule3542 = ReplacementRule(pattern3542, replacement3542)
    pattern3543 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons21)
    rule3543 = ReplacementRule(pattern3543, replacement3543)
    pattern3544 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1517)
    rule3544 = ReplacementRule(pattern3544, replacement3544)
    pattern3545 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**m_, x_), cons2, cons3, cons50, cons127, cons4, cons1517)
    rule3545 = ReplacementRule(pattern3545, replacement3545)
    pattern3546 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons1230, cons150)
    rule3546 = ReplacementRule(pattern3546, replacement3546)
    pattern3547 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons1230, cons150)
    rule3547 = ReplacementRule(pattern3547, replacement3547)
    pattern3548 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons1230, cons198, cons1542)
    rule3548 = ReplacementRule(pattern3548, replacement3548)
    pattern3549 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons1230, cons198, cons1542)
    rule3549 = ReplacementRule(pattern3549, replacement3549)
    pattern3550 = Pattern(Integral((WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons21)
    rule3550 = ReplacementRule(pattern3550, replacement3550)
    pattern3551 = Pattern(Integral((WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons21)
    rule3551 = ReplacementRule(pattern3551, replacement3551)
    pattern3552 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons5, cons87)
    rule3552 = ReplacementRule(pattern3552, replacement3552)
    pattern3553 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons5, cons87)
    rule3553 = ReplacementRule(pattern3553, replacement3553)
    pattern3554 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1441, cons20, cons1543)
    rule3554 = ReplacementRule(pattern3554, replacement3554)
    pattern3555 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1441, cons20, cons1543)
    rule3555 = ReplacementRule(pattern3555, replacement3555)
    pattern3556 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1441)
    rule3556 = ReplacementRule(pattern3556, replacement3556)
    pattern3557 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1441)
    rule3557 = ReplacementRule(pattern3557, replacement3557)
    pattern3558 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons72)
    rule3558 = ReplacementRule(pattern3558, replacement3558)
    pattern3559 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons72)
    rule3559 = ReplacementRule(pattern3559, replacement3559)
    pattern3560 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1414)
    rule3560 = ReplacementRule(pattern3560, replacement3560)
    pattern3561 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1414)
    rule3561 = ReplacementRule(pattern3561, replacement3561)
    pattern3562 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons33, cons269)
    rule3562 = ReplacementRule(pattern3562, replacement3562)
    pattern3563 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons33, cons269)
    rule3563 = ReplacementRule(pattern3563, replacement3563)
    pattern3564 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1305)
    rule3564 = ReplacementRule(pattern3564, replacement3564)
    pattern3565 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1305)
    rule3565 = ReplacementRule(pattern3565, replacement3565)
    pattern3566 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons33, cons170)
    rule3566 = ReplacementRule(pattern3566, replacement3566)
    pattern3567 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons33, cons170)
    rule3567 = ReplacementRule(pattern3567, replacement3567)
    pattern3568 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons33, cons96)
    rule3568 = ReplacementRule(pattern3568, replacement3568)
    pattern3569 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons33, cons96)
    rule3569 = ReplacementRule(pattern3569, replacement3569)
    pattern3570 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1544)
    rule3570 = ReplacementRule(pattern3570, replacement3570)
    pattern3571 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1544)
    rule3571 = ReplacementRule(pattern3571, replacement3571)
    pattern3572 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1545)
    rule3572 = ReplacementRule(pattern3572, replacement3572)
    pattern3573 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1545)
    rule3573 = ReplacementRule(pattern3573, replacement3573)
    pattern3574 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1324)
    rule3574 = ReplacementRule(pattern3574, replacement3574)
    pattern3575 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1324)
    rule3575 = ReplacementRule(pattern3575, replacement3575)
    pattern3576 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1546)
    rule3576 = ReplacementRule(pattern3576, replacement3576)
    pattern3577 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1546)
    rule3577 = ReplacementRule(pattern3577, replacement3577)
    pattern3578 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1547)
    rule3578 = ReplacementRule(pattern3578, replacement3578)
    pattern3579 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1547)
    rule3579 = ReplacementRule(pattern3579, replacement3579)
    pattern3580 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons1548)
    rule3580 = ReplacementRule(pattern3580, replacement3580)
    pattern3581 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons1548)
    rule3581 = ReplacementRule(pattern3581, replacement3581)

    pattern3582 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons1549, cons1550)
    rule3582 = ReplacementRule(pattern3582, With3582)

    pattern3583 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons1549, cons1550)
    rule3583 = ReplacementRule(pattern3583, With3583)
    pattern3584 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1546)
    rule3584 = ReplacementRule(pattern3584, replacement3584)
    pattern3585 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1546)
    rule3585 = ReplacementRule(pattern3585, replacement3585)
    pattern3586 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1547, cons79)
    rule3586 = ReplacementRule(pattern3586, replacement3586)
    pattern3587 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1547, cons79)
    rule3587 = ReplacementRule(pattern3587, replacement3587)
    pattern3588 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1547, cons21)
    rule3588 = ReplacementRule(pattern3588, replacement3588)
    pattern3589 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1547, cons21)
    rule3589 = ReplacementRule(pattern3589, replacement3589)
    pattern3590 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons34, cons1441)
    rule3590 = ReplacementRule(pattern3590, replacement3590)
    pattern3591 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons34, cons1441)
    rule3591 = ReplacementRule(pattern3591, replacement3591)
    pattern3592 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442)
    rule3592 = ReplacementRule(pattern3592, replacement3592)
    pattern3593 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442)
    rule3593 = ReplacementRule(pattern3593, replacement3593)
    pattern3594 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1442)
    rule3594 = ReplacementRule(pattern3594, replacement3594)
    pattern3595 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons33, cons96, cons1442)
    rule3595 = ReplacementRule(pattern3595, replacement3595)
    pattern3596 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1551, cons1552)
    rule3596 = ReplacementRule(pattern3596, replacement3596)
    pattern3597 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1551, cons1552)
    rule3597 = ReplacementRule(pattern3597, replacement3597)
    pattern3598 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3598 = ReplacementRule(pattern3598, replacement3598)
    pattern3599 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3599 = ReplacementRule(pattern3599, replacement3599)
    pattern3600 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons1553, cons1425)
    rule3600 = ReplacementRule(pattern3600, replacement3600)
    pattern3601 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons1553, cons1425)
    rule3601 = ReplacementRule(pattern3601, replacement3601)
    pattern3602 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons1553, cons1387)
    rule3602 = ReplacementRule(pattern3602, replacement3602)
    pattern3603 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons1553, cons1387)
    rule3603 = ReplacementRule(pattern3603, replacement3603)
    pattern3604 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons113, cons96)
    rule3604 = ReplacementRule(pattern3604, replacement3604)
    pattern3605 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons113, cons96)
    rule3605 = ReplacementRule(pattern3605, replacement3605)
    pattern3606 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1441, cons1547, cons157, cons274)
    rule3606 = ReplacementRule(pattern3606, replacement3606)
    pattern3607 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1441, cons1547, cons157, cons274)
    rule3607 = ReplacementRule(pattern3607, replacement3607)
    pattern3608 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons89, cons1327)
    rule3608 = ReplacementRule(pattern3608, replacement3608)
    pattern3609 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons89, cons1327)
    rule3609 = ReplacementRule(pattern3609, replacement3609)
    pattern3610 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons89, cons167)
    rule3610 = ReplacementRule(pattern3610, replacement3610)
    pattern3611 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons89, cons167)
    rule3611 = ReplacementRule(pattern3611, replacement3611)
    pattern3612 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3612 = ReplacementRule(pattern3612, replacement3612)
    pattern3613 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3613 = ReplacementRule(pattern3613, replacement3613)
    pattern3614 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons1329)
    rule3614 = ReplacementRule(pattern3614, replacement3614)
    pattern3615 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons1329)
    rule3615 = ReplacementRule(pattern3615, replacement3615)
    pattern3616 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons168, cons91, cons1336)
    rule3616 = ReplacementRule(pattern3616, replacement3616)
    pattern3617 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons168, cons91, cons1336)
    rule3617 = ReplacementRule(pattern3617, replacement3617)
    pattern3618 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3618 = ReplacementRule(pattern3618, replacement3618)
    pattern3619 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3619 = ReplacementRule(pattern3619, replacement3619)
    pattern3620 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3620 = ReplacementRule(pattern3620, replacement3620)
    pattern3621 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3621 = ReplacementRule(pattern3621, replacement3621)
    pattern3622 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons517, cons168, cons1521, cons1336)
    rule3622 = ReplacementRule(pattern3622, replacement3622)
    pattern3623 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons517, cons168, cons1521, cons1336)
    rule3623 = ReplacementRule(pattern3623, replacement3623)
    pattern3624 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons33, cons269, cons1554)
    rule3624 = ReplacementRule(pattern3624, replacement3624)
    pattern3625 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons33, cons269, cons1554)
    rule3625 = ReplacementRule(pattern3625, replacement3625)
    pattern3626 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons269, cons167, cons1336)
    rule3626 = ReplacementRule(pattern3626, replacement3626)
    pattern3627 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547, cons95, cons269, cons167, cons1336)
    rule3627 = ReplacementRule(pattern3627, replacement3627)
    pattern3628 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons33, cons269, cons1336)
    rule3628 = ReplacementRule(pattern3628, replacement3628)
    pattern3629 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1441, cons1547, cons33, cons269, cons1336)
    rule3629 = ReplacementRule(pattern3629, replacement3629)
    pattern3630 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547, cons89, cons167, cons1521, cons1555)
    rule3630 = ReplacementRule(pattern3630, replacement3630)
    pattern3631 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547, cons89, cons167, cons1521, cons1555)
    rule3631 = ReplacementRule(pattern3631, replacement3631)
    pattern3632 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547, cons89, cons91, cons1555)
    rule3632 = ReplacementRule(pattern3632, replacement3632)
    pattern3633 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547, cons89, cons91, cons1555)
    rule3633 = ReplacementRule(pattern3633, replacement3633)
    pattern3634 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547)
    rule3634 = ReplacementRule(pattern3634, replacement3634)
    pattern3635 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1441, cons1547)
    rule3635 = ReplacementRule(pattern3635, replacement3635)
    pattern3636 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3636 = ReplacementRule(pattern3636, replacement3636)
    pattern3637 = Pattern(Integral(sqrt(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1441, cons1547)
    rule3637 = ReplacementRule(pattern3637, replacement3637)
    pattern3638 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1441, cons1547)
    rule3638 = ReplacementRule(pattern3638, replacement3638)
    pattern3639 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1441, cons1547)
    rule3639 = ReplacementRule(pattern3639, replacement3639)
    pattern3640 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons1335, cons91, cons517)
    rule3640 = ReplacementRule(pattern3640, replacement3640)
    pattern3641 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons1335, cons91, cons517)
    rule3641 = ReplacementRule(pattern3641, replacement3641)
    pattern3642 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1442, cons1547, cons517, cons1335, cons1556, cons1557)
    rule3642 = ReplacementRule(pattern3642, replacement3642)
    pattern3643 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1442, cons1547, cons517, cons1335, cons1556, cons1557)
    rule3643 = ReplacementRule(pattern3643, replacement3643)
    pattern3644 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons96, cons1338, cons517)
    rule3644 = ReplacementRule(pattern3644, replacement3644)
    pattern3645 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons96, cons1338, cons517)
    rule3645 = ReplacementRule(pattern3645, replacement3645)
    pattern3646 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons96, cons90, cons517)
    rule3646 = ReplacementRule(pattern3646, replacement3646)
    pattern3647 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons96, cons90, cons517)
    rule3647 = ReplacementRule(pattern3647, replacement3647)
    pattern3648 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1442, cons1547, cons517, cons96, cons1558, cons1559)
    rule3648 = ReplacementRule(pattern3648, replacement3648)
    pattern3649 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1442, cons1547, cons517, cons96, cons1558, cons1559)
    rule3649 = ReplacementRule(pattern3649, replacement3649)
    pattern3650 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons168, cons90, cons810)
    rule3650 = ReplacementRule(pattern3650, replacement3650)
    pattern3651 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547, cons95, cons168, cons90, cons810)
    rule3651 = ReplacementRule(pattern3651, replacement3651)
    pattern3652 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547)
    rule3652 = ReplacementRule(pattern3652, replacement3652)
    pattern3653 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547)
    rule3653 = ReplacementRule(pattern3653, replacement3653)
    pattern3654 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons73, cons1442, cons1547)
    rule3654 = ReplacementRule(pattern3654, replacement3654)
    pattern3655 = Pattern(Integral(sqrt(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons73, cons1442, cons1547)
    rule3655 = ReplacementRule(pattern3655, replacement3655)
    pattern3656 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547)
    rule3656 = ReplacementRule(pattern3656, replacement3656)
    pattern3657 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1442, cons1547)
    rule3657 = ReplacementRule(pattern3657, replacement3657)
    pattern3658 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1547, cons21)
    rule3658 = ReplacementRule(pattern3658, replacement3658)
    pattern3659 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1442, cons1547, cons21)
    rule3659 = ReplacementRule(pattern3659, replacement3659)
    pattern3660 = Pattern(Integral((c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1442, cons1547)
    rule3660 = ReplacementRule(pattern3660, replacement3660)
    pattern3661 = Pattern(Integral((c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1442, cons1547)
    rule3661 = ReplacementRule(pattern3661, replacement3661)
    pattern3662 = Pattern(Integral((WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons25, cons20)
    rule3662 = ReplacementRule(pattern3662, replacement3662)
    pattern3663 = Pattern(Integral((WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons25, cons20)
    rule3663 = ReplacementRule(pattern3663, replacement3663)
    pattern3664 = Pattern(Integral(((WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25, cons21)
    rule3664 = ReplacementRule(pattern3664, replacement3664)
    pattern3665 = Pattern(Integral(((WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25, cons21)
    rule3665 = ReplacementRule(pattern3665, replacement3665)
    pattern3666 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule3666 = ReplacementRule(pattern3666, replacement3666)
    pattern3667 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons382)
    rule3667 = ReplacementRule(pattern3667, replacement3667)
    pattern3668 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons149, cons20, cons87)
    rule3668 = ReplacementRule(pattern3668, replacement3668)
    pattern3669 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons149, cons20, cons87)
    rule3669 = ReplacementRule(pattern3669, replacement3669)
    pattern3670 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**q_)**p_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons149, cons1419)
    rule3670 = ReplacementRule(pattern3670, replacement3670)
    pattern3671 = Pattern(Integral(((S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**q_*WC('g', S(1)))**p_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons52, cons149, cons1419)
    rule3671 = ReplacementRule(pattern3671, replacement3671)
    pattern3672 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons87)
    rule3672 = ReplacementRule(pattern3672, replacement3672)
    pattern3673 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons87)
    rule3673 = ReplacementRule(pattern3673, replacement3673)
    pattern3674 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20, cons40)
    rule3674 = ReplacementRule(pattern3674, replacement3674)
    pattern3675 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20, cons40)
    rule3675 = ReplacementRule(pattern3675, replacement3675)
    pattern3676 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons25, cons20, cons149)
    rule3676 = ReplacementRule(pattern3676, replacement3676)
    pattern3677 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons25, cons20, cons149)
    rule3677 = ReplacementRule(pattern3677, replacement3677)
    pattern3678 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons25, cons21)
    rule3678 = ReplacementRule(pattern3678, replacement3678)
    pattern3679 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons25, cons21)
    rule3679 = ReplacementRule(pattern3679, replacement3679)
    pattern3680 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1441)
    rule3680 = ReplacementRule(pattern3680, replacement3680)
    pattern3681 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1441)
    rule3681 = ReplacementRule(pattern3681, replacement3681)
    pattern3682 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73)
    rule3682 = ReplacementRule(pattern3682, replacement3682)
    pattern3683 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73)
    rule3683 = ReplacementRule(pattern3683, replacement3683)
    pattern3684 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons33, cons96, cons1441)
    rule3684 = ReplacementRule(pattern3684, replacement3684)
    pattern3685 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons33, cons96, cons1441)
    rule3685 = ReplacementRule(pattern3685, replacement3685)
    pattern3686 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons33, cons96, cons1442)
    rule3686 = ReplacementRule(pattern3686, replacement3686)
    pattern3687 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons33, cons96, cons1442)
    rule3687 = ReplacementRule(pattern3687, replacement3687)
    pattern3688 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1551)
    rule3688 = ReplacementRule(pattern3688, replacement3688)
    pattern3689 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1551)
    rule3689 = ReplacementRule(pattern3689, replacement3689)
    pattern3690 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1441, cons95, cons168, cons91)
    rule3690 = ReplacementRule(pattern3690, replacement3690)
    pattern3691 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1441, cons95, cons168, cons91)
    rule3691 = ReplacementRule(pattern3691, replacement3691)
    pattern3692 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1441, cons33, cons168, cons348)
    rule3692 = ReplacementRule(pattern3692, replacement3692)
    pattern3693 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1441, cons33, cons168, cons348)
    rule3693 = ReplacementRule(pattern3693, replacement3693)
    pattern3694 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1441, cons95, cons269, cons90)
    rule3694 = ReplacementRule(pattern3694, replacement3694)
    pattern3695 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1441, cons95, cons269, cons90)
    rule3695 = ReplacementRule(pattern3695, replacement3695)
    pattern3696 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1441, cons33, cons269, cons1329)
    rule3696 = ReplacementRule(pattern3696, replacement3696)
    pattern3697 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1441, cons33, cons269, cons1329)
    rule3697 = ReplacementRule(pattern3697, replacement3697)
    pattern3698 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons89, cons90)
    rule3698 = ReplacementRule(pattern3698, replacement3698)
    pattern3699 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons89, cons90)
    rule3699 = ReplacementRule(pattern3699, replacement3699)
    pattern3700 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons89, cons91)
    rule3700 = ReplacementRule(pattern3700, replacement3700)
    pattern3701 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons89, cons91)
    rule3701 = ReplacementRule(pattern3701, replacement3701)
    pattern3702 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1441, cons1420)
    rule3702 = ReplacementRule(pattern3702, replacement3702)
    pattern3703 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1441, cons1420)
    rule3703 = ReplacementRule(pattern3703, replacement3703)
    pattern3704 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons1428)
    rule3704 = ReplacementRule(pattern3704, replacement3704)
    pattern3705 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1441, cons1428)
    rule3705 = ReplacementRule(pattern3705, replacement3705)
    pattern3706 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1441, cons1428)
    rule3706 = ReplacementRule(pattern3706, replacement3706)
    pattern3707 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1441, cons1428)
    rule3707 = ReplacementRule(pattern3707, replacement3707)
    pattern3708 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons21, cons25, cons1515, cons1560)
    rule3708 = ReplacementRule(pattern3708, replacement3708)
    pattern3709 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons21, cons25, cons1515, cons1560)
    rule3709 = ReplacementRule(pattern3709, replacement3709)
    pattern3710 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons21, cons25, cons1515, cons1561)
    rule3710 = ReplacementRule(pattern3710, replacement3710)
    pattern3711 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons21, cons25, cons1515, cons1561)
    rule3711 = ReplacementRule(pattern3711, replacement3711)
    pattern3712 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons168, cons91, cons1336)
    rule3712 = ReplacementRule(pattern3712, replacement3712)
    pattern3713 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons168, cons91, cons1336)
    rule3713 = ReplacementRule(pattern3713, replacement3713)
    pattern3714 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547, cons33, cons168, cons1336, cons1562)
    rule3714 = ReplacementRule(pattern3714, replacement3714)
    pattern3715 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547, cons33, cons168, cons1336, cons1562)
    rule3715 = ReplacementRule(pattern3715, replacement3715)
    pattern3716 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons96, cons1327, cons1336)
    rule3716 = ReplacementRule(pattern3716, replacement3716)
    pattern3717 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons96, cons1327, cons1336)
    rule3717 = ReplacementRule(pattern3717, replacement3717)
    pattern3718 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1336, cons1559)
    rule3718 = ReplacementRule(pattern3718, replacement3718)
    pattern3719 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1336, cons1559)
    rule3719 = ReplacementRule(pattern3719, replacement3719)
    pattern3720 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons1258, cons1327)
    rule3720 = ReplacementRule(pattern3720, replacement3720)
    pattern3721 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547, cons95, cons1258, cons1327)
    rule3721 = ReplacementRule(pattern3721, replacement3721)
    pattern3722 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3722 = ReplacementRule(pattern3722, replacement3722)
    pattern3723 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3723 = ReplacementRule(pattern3723, replacement3723)
    pattern3724 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3724 = ReplacementRule(pattern3724, replacement3724)
    pattern3725 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3725 = ReplacementRule(pattern3725, replacement3725)
    pattern3726 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547)
    rule3726 = ReplacementRule(pattern3726, replacement3726)
    pattern3727 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1442, cons1547)
    rule3727 = ReplacementRule(pattern3727, replacement3727)
    pattern3728 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3728 = ReplacementRule(pattern3728, replacement3728)
    pattern3729 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1442, cons1547)
    rule3729 = ReplacementRule(pattern3729, replacement3729)
    pattern3730 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons1560)
    rule3730 = ReplacementRule(pattern3730, replacement3730)
    pattern3731 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons1560)
    rule3731 = ReplacementRule(pattern3731, replacement3731)
    pattern3732 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons1561)
    rule3732 = ReplacementRule(pattern3732, replacement3732)
    pattern3733 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1442, cons1561)
    rule3733 = ReplacementRule(pattern3733, replacement3733)
    pattern3734 = Pattern(Integral((A_ + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1563)
    rule3734 = ReplacementRule(pattern3734, replacement3734)
    pattern3735 = Pattern(Integral((A_ + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1563)
    rule3735 = ReplacementRule(pattern3735, replacement3735)
    pattern3736 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule3736 = ReplacementRule(pattern3736, replacement3736)
    pattern3737 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule3737 = ReplacementRule(pattern3737, replacement3737)
    pattern3738 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule3738 = ReplacementRule(pattern3738, replacement3738)
    pattern3739 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule3739 = ReplacementRule(pattern3739, replacement3739)
    pattern3740 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons33, cons34, cons1441)
    rule3740 = ReplacementRule(pattern3740, replacement3740)
    pattern3741 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons33, cons34, cons1441)
    rule3741 = ReplacementRule(pattern3741, replacement3741)
    pattern3742 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons33, cons34, cons1441)
    rule3742 = ReplacementRule(pattern3742, replacement3742)
    pattern3743 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons33, cons34, cons1441)
    rule3743 = ReplacementRule(pattern3743, replacement3743)
    pattern3744 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1442, cons1566)
    rule3744 = ReplacementRule(pattern3744, replacement3744)
    pattern3745 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1442, cons1566)
    rule3745 = ReplacementRule(pattern3745, replacement3745)
    pattern3746 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons50, cons127, cons36, cons37, cons38, cons1567)
    rule3746 = ReplacementRule(pattern3746, replacement3746)
    pattern3747 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons50, cons127, cons36, cons37, cons38, cons1567)
    rule3747 = ReplacementRule(pattern3747, replacement3747)
    pattern3748 = Pattern(Integral((A_ + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons50, cons127, cons36, cons38, cons1567)
    rule3748 = ReplacementRule(pattern3748, replacement3748)
    pattern3749 = Pattern(Integral((A_ + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons50, cons127, cons36, cons38, cons1567)
    rule3749 = ReplacementRule(pattern3749, replacement3749)
    pattern3750 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons1442, cons1568)
    rule3750 = ReplacementRule(pattern3750, replacement3750)
    pattern3751 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons1442, cons1568)
    rule3751 = ReplacementRule(pattern3751, replacement3751)
    pattern3752 = Pattern(Integral((A_ + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons1442, cons1567)
    rule3752 = ReplacementRule(pattern3752, replacement3752)
    pattern3753 = Pattern(Integral((A_ + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons1442, cons1567)
    rule3753 = ReplacementRule(pattern3753, replacement3753)
    pattern3754 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons33, cons96, cons1442)
    rule3754 = ReplacementRule(pattern3754, replacement3754)
    pattern3755 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons1564, cons33, cons96, cons1442)
    rule3755 = ReplacementRule(pattern3755, replacement3755)
    pattern3756 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons33, cons96, cons1442)
    rule3756 = ReplacementRule(pattern3756, replacement3756)
    pattern3757 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons36, cons38, cons1565, cons33, cons96, cons1442)
    rule3757 = ReplacementRule(pattern3757, replacement3757)
    pattern3758 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1564, cons1551)
    rule3758 = ReplacementRule(pattern3758, replacement3758)
    pattern3759 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1564, cons1551)
    rule3759 = ReplacementRule(pattern3759, replacement3759)
    pattern3760 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1565, cons1551)
    rule3760 = ReplacementRule(pattern3760, replacement3760)
    pattern3761 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1565, cons1551)
    rule3761 = ReplacementRule(pattern3761, replacement3761)
    pattern3762 = Pattern(Integral((A_ + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1563)
    rule3762 = ReplacementRule(pattern3762, replacement3762)
    pattern3763 = Pattern(Integral((A_ + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons1563)
    rule3763 = ReplacementRule(pattern3763, replacement3763)
    pattern3764 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1547, cons89, cons91)
    rule3764 = ReplacementRule(pattern3764, replacement3764)
    pattern3765 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1547, cons89, cons91)
    rule3765 = ReplacementRule(pattern3765, replacement3765)
    pattern3766 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1547, cons89, cons91)
    rule3766 = ReplacementRule(pattern3766, replacement3766)
    pattern3767 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1547, cons89, cons91)
    rule3767 = ReplacementRule(pattern3767, replacement3767)
    pattern3768 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1547, cons348)
    rule3768 = ReplacementRule(pattern3768, replacement3768)
    pattern3769 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1547, cons348)
    rule3769 = ReplacementRule(pattern3769, replacement3769)
    pattern3770 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1547, cons348)
    rule3770 = ReplacementRule(pattern3770, replacement3770)
    pattern3771 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1547, cons348)
    rule3771 = ReplacementRule(pattern3771, replacement3771)
    pattern3772 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1441, cons1569)
    rule3772 = ReplacementRule(pattern3772, replacement3772)
    pattern3773 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1441, cons1569)
    rule3773 = ReplacementRule(pattern3773, replacement3773)
    pattern3774 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1441, cons1569)
    rule3774 = ReplacementRule(pattern3774, replacement3774)
    pattern3775 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1441, cons1569)
    rule3775 = ReplacementRule(pattern3775, replacement3775)
    pattern3776 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1441, cons1305, cons89, cons91, cons1547)
    rule3776 = ReplacementRule(pattern3776, replacement3776)
    pattern3777 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1441, cons1305, cons89, cons91, cons1547)
    rule3777 = ReplacementRule(pattern3777, replacement3777)
    pattern3778 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1441, cons1305, cons89, cons91, cons1547)
    rule3778 = ReplacementRule(pattern3778, replacement3778)
    pattern3779 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1441, cons1305, cons89, cons91, cons1547)
    rule3779 = ReplacementRule(pattern3779, replacement3779)
    pattern3780 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1441, cons1305, cons685)
    rule3780 = ReplacementRule(pattern3780, replacement3780)
    pattern3781 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1441, cons1305, cons685)
    rule3781 = ReplacementRule(pattern3781, replacement3781)
    pattern3782 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1441, cons1305, cons685)
    rule3782 = ReplacementRule(pattern3782, replacement3782)
    pattern3783 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1441, cons1305, cons685)
    rule3783 = ReplacementRule(pattern3783, replacement3783)
    pattern3784 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1442, cons1547, cons95, cons170, cons91)
    rule3784 = ReplacementRule(pattern3784, replacement3784)
    pattern3785 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1442, cons1547, cons95, cons170, cons91)
    rule3785 = ReplacementRule(pattern3785, replacement3785)
    pattern3786 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1442, cons1547, cons95, cons170, cons91)
    rule3786 = ReplacementRule(pattern3786, replacement3786)
    pattern3787 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1442, cons1547, cons95, cons170, cons91)
    rule3787 = ReplacementRule(pattern3787, replacement3787)
    pattern3788 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons33, cons170, cons1570)
    rule3788 = ReplacementRule(pattern3788, replacement3788)
    pattern3789 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons33, cons170, cons1570)
    rule3789 = ReplacementRule(pattern3789, replacement3789)
    pattern3790 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons33, cons170, cons1570)
    rule3790 = ReplacementRule(pattern3790, replacement3790)
    pattern3791 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons33, cons170, cons1570)
    rule3791 = ReplacementRule(pattern3791, replacement3791)
    pattern3792 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1559)
    rule3792 = ReplacementRule(pattern3792, replacement3792)
    pattern3793 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1559)
    rule3793 = ReplacementRule(pattern3793, replacement3793)
    pattern3794 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1559)
    rule3794 = ReplacementRule(pattern3794, replacement3794)
    pattern3795 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons33, cons96, cons1559)
    rule3795 = ReplacementRule(pattern3795, replacement3795)
    pattern3796 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1442, cons1547)
    rule3796 = ReplacementRule(pattern3796, replacement3796)
    pattern3797 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1442, cons1547)
    rule3797 = ReplacementRule(pattern3797, replacement3797)
    pattern3798 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1442, cons1547)
    rule3798 = ReplacementRule(pattern3798, replacement3798)
    pattern3799 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1442, cons1547)
    rule3799 = ReplacementRule(pattern3799, replacement3799)
    pattern3800 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons1329, cons1571)
    rule3800 = ReplacementRule(pattern3800, replacement3800)
    pattern3801 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1442, cons1547, cons1329, cons1571)
    rule3801 = ReplacementRule(pattern3801, replacement3801)
    pattern3802 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons1329, cons1571)
    rule3802 = ReplacementRule(pattern3802, replacement3802)
    pattern3803 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1442, cons1547, cons1329, cons1571)
    rule3803 = ReplacementRule(pattern3803, replacement3803)
    pattern3804 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1442, cons1547)
    rule3804 = ReplacementRule(pattern3804, replacement3804)
    pattern3805 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1442, cons1547)
    rule3805 = ReplacementRule(pattern3805, replacement3805)
    pattern3806 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1442, cons1547)
    rule3806 = ReplacementRule(pattern3806, replacement3806)
    pattern3807 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1442, cons1547)
    rule3807 = ReplacementRule(pattern3807, replacement3807)
    pattern3808 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1494)
    rule3808 = ReplacementRule(pattern3808, replacement3808)
    pattern3809 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1494)
    rule3809 = ReplacementRule(pattern3809, replacement3809)
    pattern3810 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1458)
    rule3810 = ReplacementRule(pattern3810, replacement3810)
    pattern3811 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons1458)
    rule3811 = ReplacementRule(pattern3811, replacement3811)
    pattern3812 = Pattern(Integral((a_ + (WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule3812 = ReplacementRule(pattern3812, replacement3812)
    pattern3813 = Pattern(Integral((a_ + (WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1572)
    rule3813 = ReplacementRule(pattern3813, replacement3813)

    pattern3814 = Pattern(Integral((a_ + (WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1517)
    rule3814 = ReplacementRule(pattern3814, With3814)

    pattern3815 = Pattern(Integral((a_ + (WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1517)
    rule3815 = ReplacementRule(pattern3815, With3815)

    pattern3816 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons1481, cons40)
    rule3816 = ReplacementRule(pattern3816, With3816)

    pattern3817 = Pattern(Integral((a_ + (S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1483, cons1481, cons40)
    rule3817 = ReplacementRule(pattern3817, With3817)

    pattern3818 = Pattern(Integral((a_ + (WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1517)
    rule3818 = ReplacementRule(pattern3818, With3818)

    pattern3819 = Pattern(Integral((a_ + (WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1517)
    rule3819 = ReplacementRule(pattern3819, With3819)

    pattern3820 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*(S(1)/cos(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1230, cons745, cons40)
    rule3820 = ReplacementRule(pattern3820, With3820)

    pattern3821 = Pattern(Integral((a_ + (S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1230, cons745, cons40)
    rule3821 = ReplacementRule(pattern3821, With3821)
    pattern3822 = Pattern(Integral((a_ + (WC('e', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3822 = ReplacementRule(pattern3822, replacement3822)
    pattern3823 = Pattern(Integral((a_ + (WC('e', S(1))/tan(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**p_*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3823 = ReplacementRule(pattern3823, replacement3823)
    pattern3824 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule3824 = ReplacementRule(pattern3824, replacement3824)
    pattern3825 = Pattern(Integral((a_ + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule3825 = ReplacementRule(pattern3825, replacement3825)
    pattern3826 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons149)
    rule3826 = ReplacementRule(pattern3826, replacement3826)
    pattern3827 = Pattern(Integral((a_ + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons149)
    rule3827 = ReplacementRule(pattern3827, replacement3827)

    pattern3828 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule3828 = ReplacementRule(pattern3828, With3828)

    pattern3829 = Pattern(Integral(S(1)/((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule3829 = ReplacementRule(pattern3829, With3829)
    pattern3830 = Pattern(Integral(((WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1517)
    rule3830 = ReplacementRule(pattern3830, replacement3830)
    pattern3831 = Pattern(Integral(((WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1517)
    rule3831 = ReplacementRule(pattern3831, replacement3831)

    pattern3832 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons1481, cons40)
    rule3832 = ReplacementRule(pattern3832, With3832)

    pattern3833 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons1481, cons40)
    rule3833 = ReplacementRule(pattern3833, With3833)
    pattern3834 = Pattern(Integral(((WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1482)
    rule3834 = ReplacementRule(pattern3834, replacement3834)
    pattern3835 = Pattern(Integral(((WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1482)
    rule3835 = ReplacementRule(pattern3835, replacement3835)

    pattern3836 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons1481, cons40)
    rule3836 = ReplacementRule(pattern3836, With3836)

    pattern3837 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1483, cons1481, cons40)
    rule3837 = ReplacementRule(pattern3837, With3837)
    pattern3838 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3838 = ReplacementRule(pattern3838, replacement3838)
    pattern3839 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3839 = ReplacementRule(pattern3839, replacement3839)
    pattern3840 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3840 = ReplacementRule(pattern3840, replacement3840)
    pattern3841 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3841 = ReplacementRule(pattern3841, replacement3841)
    pattern3842 = Pattern(Integral(((WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons48, cons228)
    rule3842 = ReplacementRule(pattern3842, replacement3842)
    pattern3843 = Pattern(Integral(((WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons48, cons228)
    rule3843 = ReplacementRule(pattern3843, replacement3843)
    pattern3844 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3844 = ReplacementRule(pattern3844, replacement3844)
    pattern3845 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3845 = ReplacementRule(pattern3845, replacement3845)
    pattern3846 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3846 = ReplacementRule(pattern3846, replacement3846)
    pattern3847 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3847 = ReplacementRule(pattern3847, replacement3847)

    pattern3848 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons48, cons228, cons1481)
    rule3848 = ReplacementRule(pattern3848, With3848)

    pattern3849 = Pattern(Integral(((S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**n2_*WC('c', S(1)) + (S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons5, cons48, cons228, cons1481)
    rule3849 = ReplacementRule(pattern3849, With3849)
    pattern3850 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule3850 = ReplacementRule(pattern3850, replacement3850)
    pattern3851 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule3851 = ReplacementRule(pattern3851, replacement3851)
    pattern3852 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule3852 = ReplacementRule(pattern3852, replacement3852)
    pattern3853 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule3853 = ReplacementRule(pattern3853, replacement3853)

    pattern3854 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule3854 = ReplacementRule(pattern3854, With3854)

    pattern3855 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule3855 = ReplacementRule(pattern3855, With3855)
    pattern3856 = Pattern(Integral((A_ + WC('B', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule3856 = ReplacementRule(pattern3856, replacement3856)
    pattern3857 = Pattern(Integral((A_ + WC('B', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule3857 = ReplacementRule(pattern3857, replacement3857)
    pattern3858 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons64)
    rule3858 = ReplacementRule(pattern3858, replacement3858)
    pattern3859 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons64)
    rule3859 = ReplacementRule(pattern3859, replacement3859)
    pattern3860 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons170)
    rule3860 = ReplacementRule(pattern3860, replacement3860)
    pattern3861 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons170)
    rule3861 = ReplacementRule(pattern3861, replacement3861)
    pattern3862 = Pattern(Integral((WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons170)
    rule3862 = ReplacementRule(pattern3862, replacement3862)
    pattern3863 = Pattern(Integral((WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons170)
    rule3863 = ReplacementRule(pattern3863, replacement3863)
    pattern3864 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule3864 = ReplacementRule(pattern3864, replacement3864)
    pattern3865 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons530)
    rule3865 = ReplacementRule(pattern3865, replacement3865)
    pattern3866 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons33, cons170)
    rule3866 = ReplacementRule(pattern3866, replacement3866)
    pattern3867 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons33, cons170)
    rule3867 = ReplacementRule(pattern3867, replacement3867)
    pattern3868 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441)
    rule3868 = ReplacementRule(pattern3868, replacement3868)
    pattern3869 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441)
    rule3869 = ReplacementRule(pattern3869, replacement3869)
    pattern3870 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons33, cons96, cons1512)
    rule3870 = ReplacementRule(pattern3870, replacement3870)
    pattern3871 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons33, cons96, cons1512)
    rule3871 = ReplacementRule(pattern3871, replacement3871)
    pattern3872 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441)
    rule3872 = ReplacementRule(pattern3872, replacement3872)
    pattern3873 = Pattern(Integral(S(1)/((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441)
    rule3873 = ReplacementRule(pattern3873, replacement3873)
    pattern3874 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons21)
    rule3874 = ReplacementRule(pattern3874, replacement3874)
    pattern3875 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons21)
    rule3875 = ReplacementRule(pattern3875, replacement3875)
    pattern3876 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons1573)
    rule3876 = ReplacementRule(pattern3876, replacement3876)
    pattern3877 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons1573)
    rule3877 = ReplacementRule(pattern3877, replacement3877)
    pattern3878 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons198)
    rule3878 = ReplacementRule(pattern3878, replacement3878)
    pattern3879 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1441, cons198)
    rule3879 = ReplacementRule(pattern3879, replacement3879)

    pattern3880 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons1574, cons33, cons170)
    rule3880 = ReplacementRule(pattern3880, With3880)

    pattern3881 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1441, cons1574, cons33, cons170)
    rule3881 = ReplacementRule(pattern3881, With3881)
    pattern3882 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons64)
    rule3882 = ReplacementRule(pattern3882, replacement3882)
    pattern3883 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons64)
    rule3883 = ReplacementRule(pattern3883, replacement3883)
    pattern3884 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442)
    rule3884 = ReplacementRule(pattern3884, replacement3884)
    pattern3885 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))/(a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442)
    rule3885 = ReplacementRule(pattern3885, replacement3885)
    pattern3886 = Pattern(Integral((a_ + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons198, cons64)
    rule3886 = ReplacementRule(pattern3886, replacement3886)
    pattern3887 = Pattern(Integral((a_ + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1442, cons198, cons64)
    rule3887 = ReplacementRule(pattern3887, replacement3887)
    pattern3888 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule3888 = ReplacementRule(pattern3888, replacement3888)
    pattern3889 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule3889 = ReplacementRule(pattern3889, replacement3889)
    pattern3890 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule3890 = ReplacementRule(pattern3890, replacement3890)
    pattern3891 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule3891 = ReplacementRule(pattern3891, replacement3891)
    pattern3892 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule3892 = ReplacementRule(pattern3892, replacement3892)
    pattern3893 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1575, cons40)
    rule3893 = ReplacementRule(pattern3893, replacement3893)
    pattern3894 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule3894 = ReplacementRule(pattern3894, replacement3894)
    pattern3895 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule3895 = ReplacementRule(pattern3895, replacement3895)
    pattern3896 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule3896 = ReplacementRule(pattern3896, replacement3896)
    pattern3897 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71)
    rule3897 = ReplacementRule(pattern3897, replacement3897)
    pattern3898 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*tan(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule3898 = ReplacementRule(pattern3898, replacement3898)
    pattern3899 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/tan(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule3899 = ReplacementRule(pattern3899, replacement3899)
    pattern3900 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule3900 = ReplacementRule(pattern3900, replacement3900)
    pattern3901 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1576, cons40)
    rule3901 = ReplacementRule(pattern3901, replacement3901)
    pattern3902 = Pattern(Integral(x_**WC('m', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons8, cons29, cons19, cons4, cons1577)
    rule3902 = ReplacementRule(pattern3902, replacement3902)
    pattern3903 = Pattern(Integral(x_**WC('m', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0)))**S(2), x_), cons8, cons29, cons19, cons4, cons1577)
    rule3903 = ReplacementRule(pattern3903, replacement3903)
    pattern3904 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule3904 = ReplacementRule(pattern3904, replacement3904)
    pattern3905 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons5, cons1578)
    rule3905 = ReplacementRule(pattern3905, replacement3905)
    pattern3906 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3906 = ReplacementRule(pattern3906, replacement3906)
    pattern3907 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3907 = ReplacementRule(pattern3907, replacement3907)
    pattern3908 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*tan(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule3908 = ReplacementRule(pattern3908, replacement3908)
    pattern3909 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))/tan(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule3909 = ReplacementRule(pattern3909, replacement3909)
    pattern3910 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*tan(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1580)
    rule3910 = ReplacementRule(pattern3910, replacement3910)
    pattern3911 = Pattern(Integral(x_**WC('m', S(1))*(S(1)/sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('p', S(1))*(S(1)/tan(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))))**WC('q', S(1)), x_), cons2, cons3, cons5, cons33, cons87, cons1579, cons1580)
    rule3911 = ReplacementRule(pattern3911, replacement3911)
    pattern3912 = Pattern(Integral(tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule3912 = ReplacementRule(pattern3912, replacement3912)
    pattern3913 = Pattern(Integral((S(1)/tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons4, cons1581)
    rule3913 = ReplacementRule(pattern3913, replacement3913)
    pattern3914 = Pattern(Integral((d_ + x_*WC('e', S(1)))*tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons49)
    rule3914 = ReplacementRule(pattern3914, replacement3914)
    pattern3915 = Pattern(Integral((d_ + x_*WC('e', S(1)))/tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons49)
    rule3915 = ReplacementRule(pattern3915, replacement3915)
    pattern3916 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons241)
    rule3916 = ReplacementRule(pattern3916, replacement3916)
    pattern3917 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons241)
    rule3917 = ReplacementRule(pattern3917, replacement3917)
    pattern3918 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule3918 = ReplacementRule(pattern3918, replacement3918)
    pattern3919 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(S(1)/tan(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons1582)
    rule3919 = ReplacementRule(pattern3919, replacement3919)
    return [rule3401, rule3402, rule3403, rule3404, rule3405, rule3406, rule3407, rule3408, rule3409, rule3410, rule3411, rule3412, rule3413, rule3414, rule3415, rule3416, rule3417, rule3418, rule3419, rule3420, rule3421, rule3422, rule3423, rule3424, rule3425, rule3426, rule3427, rule3428, rule3429, rule3430, rule3431, rule3432, rule3433, rule3434, rule3435, rule3436, rule3437, rule3438, rule3439, rule3440, rule3441, rule3442, rule3443, rule3444, rule3445, rule3446, rule3447, rule3448, rule3449, rule3450, rule3451, rule3452, rule3453, rule3454, rule3455, rule3456, rule3457, rule3458, rule3459, rule3460, rule3461, rule3462, rule3463, rule3464, rule3465, rule3466, rule3467, rule3468, rule3469, rule3470, rule3471, rule3472, rule3473, rule3474, rule3475, rule3476, rule3477, rule3478, rule3479, rule3480, rule3481, rule3482, rule3483, rule3484, rule3485, rule3486, rule3487, rule3488, rule3489, rule3490, rule3491, rule3492, rule3493, rule3494, rule3495, rule3496, rule3497, rule3498, rule3499, rule3500, rule3501, rule3502, rule3503, rule3504, rule3505, rule3506, rule3507, rule3508, rule3509, rule3510, rule3511, rule3512, rule3513, rule3514, rule3515, rule3516, rule3517, rule3518, rule3519, rule3520, rule3521, rule3522, rule3523, rule3524, rule3525, rule3526, rule3527, rule3528, rule3529, rule3530, rule3531, rule3532, rule3533, rule3534, rule3535, rule3536, rule3537, rule3538, rule3539, rule3540, rule3541, rule3542, rule3543, rule3544, rule3545, rule3546, rule3547, rule3548, rule3549, rule3550, rule3551, rule3552, rule3553, rule3554, rule3555, rule3556, rule3557, rule3558, rule3559, rule3560, rule3561, rule3562, rule3563, rule3564, rule3565, rule3566, rule3567, rule3568, rule3569, rule3570, rule3571, rule3572, rule3573, rule3574, rule3575, rule3576, rule3577, rule3578, rule3579, rule3580, rule3581, rule3582, rule3583, rule3584, rule3585, rule3586, rule3587, rule3588, rule3589, rule3590, rule3591, rule3592, rule3593, rule3594, rule3595, rule3596, rule3597, rule3598, rule3599, rule3600, rule3601, rule3602, rule3603, rule3604, rule3605, rule3606, rule3607, rule3608, rule3609, rule3610, rule3611, rule3612, rule3613, rule3614, rule3615, rule3616, rule3617, rule3618, rule3619, rule3620, rule3621, rule3622, rule3623, rule3624, rule3625, rule3626, rule3627, rule3628, rule3629, rule3630, rule3631, rule3632, rule3633, rule3634, rule3635, rule3636, rule3637, rule3638, rule3639, rule3640, rule3641, rule3642, rule3643, rule3644, rule3645, rule3646, rule3647, rule3648, rule3649, rule3650, rule3651, rule3652, rule3653, rule3654, rule3655, rule3656, rule3657, rule3658, rule3659, rule3660, rule3661, rule3662, rule3663, rule3664, rule3665, rule3666, rule3667, rule3668, rule3669, rule3670, rule3671, rule3672, rule3673, rule3674, rule3675, rule3676, rule3677, rule3678, rule3679, rule3680, rule3681, rule3682, rule3683, rule3684, rule3685, rule3686, rule3687, rule3688, rule3689, rule3690, rule3691, rule3692, rule3693, rule3694, rule3695, rule3696, rule3697, rule3698, rule3699, rule3700, rule3701, rule3702, rule3703, rule3704, rule3705, rule3706, rule3707, rule3708, rule3709, rule3710, rule3711, rule3712, rule3713, rule3714, rule3715, rule3716, rule3717, rule3718, rule3719, rule3720, rule3721, rule3722, rule3723, rule3724, rule3725, rule3726, rule3727, rule3728, rule3729, rule3730, rule3731, rule3732, rule3733, rule3734, rule3735, rule3736, rule3737, rule3738, rule3739, rule3740, rule3741, rule3742, rule3743, rule3744, rule3745, rule3746, rule3747, rule3748, rule3749, rule3750, rule3751, rule3752, rule3753, rule3754, rule3755, rule3756, rule3757, rule3758, rule3759, rule3760, rule3761, rule3762, rule3763, rule3764, rule3765, rule3766, rule3767, rule3768, rule3769, rule3770, rule3771, rule3772, rule3773, rule3774, rule3775, rule3776, rule3777, rule3778, rule3779, rule3780, rule3781, rule3782, rule3783, rule3784, rule3785, rule3786, rule3787, rule3788, rule3789, rule3790, rule3791, rule3792, rule3793, rule3794, rule3795, rule3796, rule3797, rule3798, rule3799, rule3800, rule3801, rule3802, rule3803, rule3804, rule3805, rule3806, rule3807, rule3808, rule3809, rule3810, rule3811, rule3812, rule3813, rule3814, rule3815, rule3816, rule3817, rule3818, rule3819, rule3820, rule3821, rule3822, rule3823, rule3824, rule3825, rule3826, rule3827, rule3828, rule3829, rule3830, rule3831, rule3832, rule3833, rule3834, rule3835, rule3836, rule3837, rule3838, rule3839, rule3840, rule3841, rule3842, rule3843, rule3844, rule3845, rule3846, rule3847, rule3848, rule3849, rule3850, rule3851, rule3852, rule3853, rule3854, rule3855, rule3856, rule3857, rule3858, rule3859, rule3860, rule3861, rule3862, rule3863, rule3864, rule3865, rule3866, rule3867, rule3868, rule3869, rule3870, rule3871, rule3872, rule3873, rule3874, rule3875, rule3876, rule3877, rule3878, rule3879, rule3880, rule3881, rule3882, rule3883, rule3884, rule3885, rule3886, rule3887, rule3888, rule3889, rule3890, rule3891, rule3892, rule3893, rule3894, rule3895, rule3896, rule3897, rule3898, rule3899, rule3900, rule3901, rule3902, rule3903, rule3904, rule3905, rule3906, rule3907, rule3908, rule3909, rule3910, rule3911, rule3912, rule3913, rule3914, rule3915, rule3916, rule3917, rule3918, rule3919, ]



def replacement3401(m, n, x, f, a, e, b):
        # rubi.append(3401)
        return -Simp(b*(a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3402(m, n, x, f, a, e, b):
        # rubi.append(3402)
        return Simp(b*(a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3403(m, n, x, f, e):
        # rubi.append(3403)
        return -Dist(S(1)/f, Subst(Int(x**(-n)*(S(1) - x**S(2))**(m/S(2) + n/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)
def replacement3404(m, n, x, f, e):
        # rubi.append(3404)
        return Dist(S(1)/f, Subst(Int(x**(-n)*(S(1) - x**S(2))**(m/S(2) + n/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)
def replacement3405(m, n, x, f, e, b):
        # rubi.append(3405)
        return Dist(b**(-m), Int((b*tan(e + f*x))**(m + n)*(S(1)/cos(e + f*x))**(-m), x), x)
def replacement3406(m, n, x, f, e, b):
        # rubi.append(3406)
        return Dist(b**(-m), Int((b/tan(e + f*x))**(m + n)*(S(1)/sin(e + f*x))**(-m), x), x)
def replacement3407(m, n, x, f, a, e, b):
        # rubi.append(3407)
        return -Dist(b**S(2)*(m + S(2))/(a**S(2)*(n + S(-1))), Int((a*sin(e + f*x))**(m + S(2))*(b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a*sin(e + f*x))**(m + S(2))*(b*tan(e + f*x))**(n + S(-1))/(a**S(2)*f*(n + S(-1))), x)
def replacement3408(m, n, x, f, a, e, b):
        # rubi.append(3408)
        return -Dist(b**S(2)*(m + S(2))/(a**S(2)*(n + S(-1))), Int((a*cos(e + f*x))**(m + S(2))*(b/tan(e + f*x))**(n + S(-2)), x), x) - Simp(b*(a*cos(e + f*x))**(m + S(2))*(b/tan(e + f*x))**(n + S(-1))/(a**S(2)*f*(n + S(-1))), x)
def replacement3409(m, n, x, f, a, e, b):
        # rubi.append(3409)
        return -Dist(b**S(2)*(m + n + S(-1))/(n + S(-1)), Int((a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)
def replacement3410(m, n, x, f, a, e, b):
        # rubi.append(3410)
        return -Dist(b**S(2)*(m + n + S(-1))/(n + S(-1)), Int((a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(-2)), x), x) - Simp(b*(a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)
def replacement3411(m, n, x, f, a, e, b):
        # rubi.append(3411)
        return -Dist(a**S(2)*(n + S(1))/(b**S(2)*m), Int((a*sin(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**(n + S(2)), x), x) + Simp((a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3412(m, n, x, f, a, e, b):
        # rubi.append(3412)
        return -Dist(a**S(2)*(n + S(1))/(b**S(2)*m), Int((a*cos(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**(n + S(2)), x), x) - Simp((a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3413(m, n, x, f, a, e, b):
        # rubi.append(3413)
        return -Dist((n + S(1))/(b**S(2)*(m + n + S(1))), Int((a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(2)), x), x) + Simp((a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(1))), x)
def replacement3414(m, n, x, f, a, e, b):
        # rubi.append(3414)
        return -Dist((n + S(1))/(b**S(2)*(m + n + S(1))), Int((a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(2)), x), x) - Simp((a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(1))), x)
def replacement3415(m, n, x, f, a, e, b):
        # rubi.append(3415)
        return Dist(a**S(2)*(m + n + S(-1))/m, Int((a*sin(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**n, x), x) - Simp(b*(a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3416(m, n, x, f, a, e, b):
        # rubi.append(3416)
        return Dist(a**S(2)*(m + n + S(-1))/m, Int((a*cos(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**n, x), x) + Simp(b*(a*cos(e + f*x))**m*(b/tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3417(m, n, x, f, a, e, b):
        # rubi.append(3417)
        return Dist((m + S(2))/(a**S(2)*(m + n + S(1))), Int((a*sin(e + f*x))**(m + S(2))*(b*tan(e + f*x))**n, x), x) + Simp(b*(a*sin(e + f*x))**(m + S(2))*(b*tan(e + f*x))**(n + S(-1))/(a**S(2)*f*(m + n + S(1))), x)
def replacement3418(m, n, x, f, a, e, b):
        # rubi.append(3418)
        return Dist((m + S(2))/(a**S(2)*(m + n + S(1))), Int((a*cos(e + f*x))**(m + S(2))*(b/tan(e + f*x))**n, x), x) - Simp(b*(a*cos(e + f*x))**(m + S(2))*(b/tan(e + f*x))**(n + S(-1))/(a**S(2)*f*(m + n + S(1))), x)
def replacement3419(m, n, x, f, a, e):
        # rubi.append(3419)
        return Dist(S(1)/f, Subst(Int(x**(m + n)*(a**S(2) - x**S(2))**(-n/S(2) + S(-1)/2), x), x, a*sin(e + f*x)), x)
def replacement3420(m, n, x, f, a, e):
        # rubi.append(3420)
        return -Dist(S(1)/f, Subst(Int(x**(m + n)*(a**S(2) - x**S(2))**(-n/S(2) + S(-1)/2), x), x, a*cos(e + f*x)), x)
def replacement3421(m, n, x, f, a, e, b):
        # rubi.append(3421)
        return Dist(a**(S(1) - S(2)*IntPart(n/S(2) + S(1)/2))*b**(S(2)*IntPart(n/S(2) + S(1)/2) + S(-1))*(a*sin(e + f*x))**(-S(2)*FracPart(n/S(2) + S(1)/2))*(b*tan(e + f*x))**(S(2)*FracPart(n/S(2) + S(1)/2))*(cos(e + f*x)**S(2))**FracPart(n/S(2) + S(1)/2)/f, Subst(Int((a*x)**(m + n)*(S(1) - x**S(2))**(-n/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)
def replacement3422(m, n, x, f, a, e, b):
        # rubi.append(3422)
        return -Dist(a**(S(1) - S(2)*IntPart(n/S(2) + S(1)/2))*b**(S(2)*IntPart(n/S(2) + S(1)/2) + S(-1))*(a*cos(e + f*x))**(-S(2)*FracPart(n/S(2) + S(1)/2))*(b/tan(e + f*x))**(S(2)*FracPart(n/S(2) + S(1)/2))*(sin(e + f*x)**S(2))**FracPart(n/S(2) + S(1)/2)/f, Subst(Int((a*x)**(m + n)*(S(1) - x**S(2))**(-n/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)
def replacement3423(m, n, x, f, a, e, b):
        # rubi.append(3423)
        return Dist((S(1)/(a*cos(e + f*x)))**FracPart(m)*(a*cos(e + f*x))**FracPart(m), Int((S(1)/(a*cos(e + f*x)))**(-m)*(b*tan(e + f*x))**n, x), x)
def replacement3424(m, n, x, f, a, e, b):
        # rubi.append(3424)
        return Dist((S(1)/(a*sin(e + f*x)))**FracPart(m)*(a*sin(e + f*x))**FracPart(m), Int((S(1)/(a*sin(e + f*x)))**(-m)*(b/tan(e + f*x))**n, x), x)
def replacement3425(m, n, x, f, a, e, b):
        # rubi.append(3425)
        return Dist((a/tan(e + f*x))**m*(b*tan(e + f*x))**m, Int((b*tan(e + f*x))**(-m + n), x), x)
def replacement3426(m, n, x, f, a, e, b):
        # rubi.append(3426)
        return -Simp((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3427(m, n, x, f, a, e, b):
        # rubi.append(3427)
        return Simp((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3428(m, n, x, f, a, e, b):
        # rubi.append(3428)
        return Dist(a/f, Subst(Int((a*x)**(m + S(-1))*(x**S(2) + S(-1))**(n/S(2) + S(-1)/2), x), x, S(1)/cos(e + f*x)), x)
def replacement3429(m, n, x, f, a, e, b):
        # rubi.append(3429)
        return -Dist(a/f, Subst(Int((a*x)**(m + S(-1))*(x**S(2) + S(-1))**(n/S(2) + S(-1)/2), x), x, S(1)/sin(e + f*x)), x)
def replacement3430(m, n, x, f, e, b):
        # rubi.append(3430)
        return Dist(S(1)/f, Subst(Int((b*x)**n*(x**S(2) + S(1))**(m/S(2) + S(-1)), x), x, tan(e + f*x)), x)
def replacement3431(m, n, x, f, e, b):
        # rubi.append(3431)
        return -Dist(S(1)/f, Subst(Int((b*x)**n*(x**S(2) + S(1))**(m/S(2) + S(-1)), x), x, S(1)/tan(e + f*x)), x)
def replacement3432(m, n, x, f, a, e, b):
        # rubi.append(3432)
        return -Dist(a**S(2)*(m + S(-2))/(b**S(2)*(n + S(1))), Int((a/cos(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**(n + S(2)), x), x) + Simp(a**S(2)*(a/cos(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3433(m, n, x, f, a, e, b):
        # rubi.append(3433)
        return -Dist(a**S(2)*(m + S(-2))/(b**S(2)*(n + S(1))), Int((a/sin(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**(n + S(2)), x), x) - Simp(a**S(2)*(a/sin(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3434(m, n, x, f, a, e, b):
        # rubi.append(3434)
        return -Dist((m + n + S(1))/(b**S(2)*(n + S(1))), Int((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(2)), x), x) + Simp((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3435(m, n, x, f, a, e, b):
        # rubi.append(3435)
        return -Dist((m + n + S(1))/(b**S(2)*(n + S(1))), Int((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(2)), x), x) - Simp((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3436(m, n, x, f, a, e, b):
        # rubi.append(3436)
        return -Dist(b**S(2)*(n + S(-1))/(a**S(2)*m), Int((a/cos(e + f*x))**(m + S(2))*(b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3437(m, n, x, f, a, e, b):
        # rubi.append(3437)
        return -Dist(b**S(2)*(n + S(-1))/(a**S(2)*m), Int((a/sin(e + f*x))**(m + S(2))*(b/tan(e + f*x))**(n + S(-2)), x), x) - Simp(b*(a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3438(m, n, x, f, a, e, b):
        # rubi.append(3438)
        return -Dist(b**S(2)*(n + S(-1))/(m + n + S(-1)), Int((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3439(m, n, x, f, a, e, b):
        # rubi.append(3439)
        return -Dist(b**S(2)*(n + S(-1))/(m + n + S(-1)), Int((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(-2)), x), x) - Simp(b*(a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3440(m, n, x, f, a, e, b):
        # rubi.append(3440)
        return Dist((m + n + S(1))/(a**S(2)*m), Int((a/cos(e + f*x))**(m + S(2))*(b*tan(e + f*x))**n, x), x) - Simp((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3441(m, n, x, f, a, e, b):
        # rubi.append(3441)
        return Dist((m + n + S(1))/(a**S(2)*m), Int((a/sin(e + f*x))**(m + S(2))*(b/tan(e + f*x))**n, x), x) + Simp((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))/(b*f*m), x)
def replacement3442(m, n, x, f, a, e, b):
        # rubi.append(3442)
        return Dist(a**S(2)*(m + S(-2))/(m + n + S(-1)), Int((a/cos(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**n, x), x) + Simp(a**S(2)*(a/cos(e + f*x))**(m + S(-2))*(b*tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(-1))), x)
def replacement3443(m, n, x, f, a, e, b):
        # rubi.append(3443)
        return Dist(a**S(2)*(m + S(-2))/(m + n + S(-1)), Int((a/sin(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**n, x), x) - Simp(a**S(2)*(a/sin(e + f*x))**(m + S(-2))*(b/tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(-1))), x)
def replacement3444(x, e, f, b):
        # rubi.append(3444)
        return Dist(sqrt(sin(e + f*x))/(sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))), Int(S(1)/(sqrt(sin(e + f*x))*sqrt(cos(e + f*x))), x), x)
def replacement3445(x, e, f, b):
        # rubi.append(3445)
        return Dist(sqrt(cos(e + f*x))/(sqrt(b/tan(e + f*x))*sqrt(sin(e + f*x))), Int(S(1)/(sqrt(sin(e + f*x))*sqrt(cos(e + f*x))), x), x)
def replacement3446(x, e, f, b):
        # rubi.append(3446)
        return Dist(sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))/sqrt(sin(e + f*x)), Int(sqrt(sin(e + f*x))*sqrt(cos(e + f*x)), x), x)
def replacement3447(x, e, f, b):
        # rubi.append(3447)
        return Dist(sqrt(b/tan(e + f*x))*sqrt(sin(e + f*x))/sqrt(cos(e + f*x)), Int(sqrt(sin(e + f*x))*sqrt(cos(e + f*x)), x), x)
def replacement3448(m, n, x, f, a, e, b):
        # rubi.append(3448)
        return Dist(a**(m + n)*(a/cos(e + f*x))**(-n)*(b*sin(e + f*x))**(-n)*(b*tan(e + f*x))**n, Int((b*sin(e + f*x))**n*cos(e + f*x)**(-m - n), x), x)
def replacement3449(m, n, x, f, a, e, b):
        # rubi.append(3449)
        return Dist(a**(m + n)*(a/sin(e + f*x))**(-n)*(b*cos(e + f*x))**(-n)*(b/tan(e + f*x))**n, Int((b*cos(e + f*x))**n*sin(e + f*x)**(-m - n), x), x)
def replacement3450(m, n, x, f, a, e, b):
        # rubi.append(3450)
        return Simp((a/cos(e + f*x))**m*(b*tan(e + f*x))**(n + S(1))*(cos(e + f*x)**S(2))**(m/S(2) + n/S(2) + S(1)/2)*Hypergeometric2F1(n/S(2) + S(1)/2, m/S(2) + n/S(2) + S(1)/2, n/S(2) + S(3)/2, sin(e + f*x)**S(2))/(b*f*(n + S(1))), x)
def replacement3451(m, n, x, f, a, e, b):
        # rubi.append(3451)
        return -Simp((a/sin(e + f*x))**m*(b/tan(e + f*x))**(n + S(1))*(sin(e + f*x)**S(2))**(m/S(2) + n/S(2) + S(1)/2)*Hypergeometric2F1(n/S(2) + S(1)/2, m/S(2) + n/S(2) + S(1)/2, n/S(2) + S(3)/2, cos(e + f*x)**S(2))/(b*f*(n + S(1))), x)
def replacement3452(m, n, x, f, a, e, b):
        # rubi.append(3452)
        return Dist((sin(e + f*x)/a)**FracPart(m)*(a/sin(e + f*x))**FracPart(m), Int((sin(e + f*x)/a)**(-m)*(b*tan(e + f*x))**n, x), x)
def replacement3453(m, n, x, f, a, e, b):
        # rubi.append(3453)
        return Dist((cos(e + f*x)/a)**FracPart(m)*(a/cos(e + f*x))**FracPart(m), Int((cos(e + f*x)/a)**(-m)*(b/tan(e + f*x))**n, x), x)
def replacement3454(m, n, x, d, f, a, p, e, b):
        # rubi.append(3454)
        return -Dist(b**S(2)*(n + S(-1))/(m*p + n + S(-1)), Int((a*(d/cos(e + f*x))**p)**m*(b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(b*(a*(d/cos(e + f*x))**p)**m*(b*tan(e + f*x))**(n + S(-1))/(f*(m*p + n + S(-1))), x)
def replacement3455(m, n, x, d, f, a, p, e, b):
        # rubi.append(3455)
        return -Dist(b**S(2)*(n + S(-1))/(m*p + n + S(-1)), Int((a*(d/sin(e + f*x))**p)**m*(b/tan(e + f*x))**(n + S(-2)), x), x) - Simp(b*(a*(d/sin(e + f*x))**p)**m*(b/tan(e + f*x))**(n + S(-1))/(f*(m*p + n + S(-1))), x)
def replacement3456(m, n, x, d, f, a, p, e, b):
        # rubi.append(3456)
        return -Dist((m*p + n + S(1))/(b**S(2)*(n + S(1))), Int((a*(d/cos(e + f*x))**p)**m*(b*tan(e + f*x))**(n + S(2)), x), x) + Simp((a*(d/cos(e + f*x))**p)**m*(b*tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3457(m, n, x, d, f, a, p, e, b):
        # rubi.append(3457)
        return -Dist((m*p + n + S(1))/(b**S(2)*(n + S(1))), Int((a*(d/sin(e + f*x))**p)**m*(b/tan(e + f*x))**(n + S(2)), x), x) + Simp((a*(d/sin(e + f*x))**p)**m*(b/tan(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)
def replacement3458(c, n, x, d, b):
        # rubi.append(3458)
        return -Dist(b**S(2), Int((b*tan(c + d*x))**(n + S(-2)), x), x) + Simp(b*(b*tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3459(c, n, x, d, b):
        # rubi.append(3459)
        return -Dist(b**S(2), Int((b/tan(c + d*x))**(n + S(-2)), x), x) - Simp(b*(b/tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3460(c, n, x, d, b):
        # rubi.append(3460)
        return -Dist(b**(S(-2)), Int((b*tan(c + d*x))**(n + S(2)), x), x) + Simp((b*tan(c + d*x))**(n + S(1))/(b*d*(n + S(1))), x)
def replacement3461(c, n, x, d, b):
        # rubi.append(3461)
        return -Dist(b**(S(-2)), Int((b/tan(c + d*x))**(n + S(2)), x), x) - Simp((b/tan(c + d*x))**(n + S(1))/(b*d*(n + S(1))), x)
def replacement3462(c, d, x):
        # rubi.append(3462)
        return -Simp(log(RemoveContent(cos(c + d*x), x))/d, x)
def replacement3463(c, d, x):
        # rubi.append(3463)
        return Simp(log(RemoveContent(sin(c + d*x), x))/d, x)
def replacement3464(c, n, x, d, b):
        # rubi.append(3464)
        return Dist(b/d, Subst(Int(x**n/(b**S(2) + x**S(2)), x), x, b*tan(c + d*x)), x)
def replacement3465(c, n, x, d, b):
        # rubi.append(3465)
        return -Dist(b/d, Subst(Int(x**n/(b**S(2) + x**S(2)), x), x, b/tan(c + d*x)), x)
def replacement3466(c, x, d, a, b):
        # rubi.append(3466)
        return Dist(S(2)*a*b, Int(tan(c + d*x), x), x) + Simp(x*(a**S(2) - b**S(2)), x) + Simp(b**S(2)*tan(c + d*x)/d, x)
def replacement3467(c, x, d, a, b):
        # rubi.append(3467)
        return Dist(S(2)*a*b, Int(S(1)/tan(c + d*x), x), x) + Simp(x*(a**S(2) - b**S(2)), x) - Simp(b**S(2)/(d*tan(c + d*x)), x)
def replacement3468(c, n, x, d, a, b):
        # rubi.append(3468)
        return Dist(S(2)*a, Int((a + b*tan(c + d*x))**(n + S(-1)), x), x) + Simp(b*(a + b*tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3469(c, n, x, d, a, b):
        # rubi.append(3469)
        return Dist(S(2)*a, Int((a + b/tan(c + d*x))**(n + S(-1)), x), x) - Simp(b*(a + b/tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3470(c, n, x, d, a, b):
        # rubi.append(3470)
        return Dist(S(1)/(S(2)*a), Int((a + b*tan(c + d*x))**(n + S(1)), x), x) + Simp(a*(a + b*tan(c + d*x))**n/(S(2)*b*d*n), x)
def replacement3471(c, n, x, d, a, b):
        # rubi.append(3471)
        return Dist(S(1)/(S(2)*a), Int((a + b/tan(c + d*x))**(n + S(1)), x), x) - Simp(a*(a + b/tan(c + d*x))**n/(S(2)*b*d*n), x)
def replacement3472(c, x, d, a, b):
        # rubi.append(3472)
        return Dist(-S(2)*b/d, Subst(Int(S(1)/(S(2)*a - x**S(2)), x), x, sqrt(a + b*tan(c + d*x))), x)
def replacement3473(c, x, d, a, b):
        # rubi.append(3473)
        return Dist(S(2)*b/d, Subst(Int(S(1)/(S(2)*a - x**S(2)), x), x, sqrt(a + b/tan(c + d*x))), x)
def replacement3474(c, n, x, d, a, b):
        # rubi.append(3474)
        return -Dist(b/d, Subst(Int((a + x)**(n + S(-1))/(a - x), x), x, b*tan(c + d*x)), x)
def replacement3475(c, n, x, d, a, b):
        # rubi.append(3475)
        return Dist(b/d, Subst(Int((a + x)**(n + S(-1))/(a - x), x), x, b/tan(c + d*x)), x)
def replacement3476(c, n, x, d, a, b):
        # rubi.append(3476)
        return Int((a + b*tan(c + d*x))**(n + S(-2))*(a**S(2) + S(2)*a*b*tan(c + d*x) - b**S(2)), x) + Simp(b*(a + b*tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3477(c, n, x, d, a, b):
        # rubi.append(3477)
        return Int((a + b/tan(c + d*x))**(n + S(-2))*(a**S(2) + S(2)*a*b/tan(c + d*x) - b**S(2)), x) - Simp(b*(a + b/tan(c + d*x))**(n + S(-1))/(d*(n + S(-1))), x)
def replacement3478(c, n, x, d, a, b):
        # rubi.append(3478)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a - b*tan(c + d*x))*(a + b*tan(c + d*x))**(n + S(1)), x), x) + Simp(b*(a + b*tan(c + d*x))**(n + S(1))/(d*(a**S(2) + b**S(2))*(n + S(1))), x)
def replacement3479(c, n, x, d, a, b):
        # rubi.append(3479)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a - b/tan(c + d*x))*(a + b/tan(c + d*x))**(n + S(1)), x), x) - Simp(b*(a + b/tan(c + d*x))**(n + S(1))/(d*(a**S(2) + b**S(2))*(n + S(1))), x)
def replacement3480(c, x, d, a, b):
        # rubi.append(3480)
        return Dist(b/(a**S(2) + b**S(2)), Int((-a*tan(c + d*x) + b)/(a + b*tan(c + d*x)), x), x) + Simp(a*x/(a**S(2) + b**S(2)), x)
def replacement3481(c, x, d, a, b):
        # rubi.append(3481)
        return Dist(b/(a**S(2) + b**S(2)), Int((-a/tan(c + d*x) + b)/(a + b/tan(c + d*x)), x), x) + Simp(a*x/(a**S(2) + b**S(2)), x)
def replacement3482(c, n, x, d, a, b):
        # rubi.append(3482)
        return Dist(b/d, Subst(Int((a + x)**n/(b**S(2) + x**S(2)), x), x, b*tan(c + d*x)), x)
def replacement3483(c, n, x, d, a, b):
        # rubi.append(3483)
        return -Dist(b/d, Subst(Int((a + x)**n/(b**S(2) + x**S(2)), x), x, b/tan(c + d*x)), x)
def replacement3484(m, x, d, f, a, e, b):
        # rubi.append(3484)
        return Dist(a, Int((d/cos(e + f*x))**m, x), x) + Simp(b*(d/cos(e + f*x))**m/(f*m), x)
def replacement3485(m, x, d, f, a, e, b):
        # rubi.append(3485)
        return Dist(a, Int((d/sin(e + f*x))**m, x), x) - Simp(b*(d/sin(e + f*x))**m/(f*m), x)
def replacement3486(m, n, x, f, a, e, b):
        # rubi.append(3486)
        return Dist(a**(S(2) - m)/(b*f), Subst(Int((a - x)**(m/S(2) + S(-1))*(a + x)**(m/S(2) + n + S(-1)), x), x, b*tan(e + f*x)), x)
def replacement3487(m, n, x, f, a, e, b):
        # rubi.append(3487)
        return -Dist(a**(S(2) - m)/(b*f), Subst(Int((a - x)**(m/S(2) + S(-1))*(a + x)**(m/S(2) + n + S(-1)), x), x, b/tan(e + f*x)), x)
def replacement3488(m, n, x, d, f, a, e, b):
        # rubi.append(3488)
        return Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**n/(a*f*m), x)
def replacement3489(m, n, x, d, f, a, e, b):
        # rubi.append(3489)
        return -Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**n/(a*f*m), x)
def replacement3490(x, f, a, e, b):
        # rubi.append(3490)
        return Dist(-S(2)*a/(b*f), Subst(Int(S(1)/(-a*x**S(2) + S(2)), x), x, S(1)/(sqrt(a + b*tan(e + f*x))*cos(e + f*x))), x)
def replacement3491(x, f, a, e, b):
        # rubi.append(3491)
        return Dist(S(2)*a/(b*f), Subst(Int(S(1)/(-a*x**S(2) + S(2)), x), x, S(1)/(sqrt(a + b/tan(e + f*x))*sin(e + f*x))), x)
def replacement3492(m, n, x, d, f, a, e, b):
        # rubi.append(3492)
        return Dist(a/(S(2)*d**S(2)), Int((d/cos(e + f*x))**(m + S(2))*(a + b*tan(e + f*x))**(n + S(-1)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**n/(a*f*m), x)
def replacement3493(m, n, x, d, f, a, e, b):
        # rubi.append(3493)
        return Dist(a/(S(2)*d**S(2)), Int((d/sin(e + f*x))**(m + S(2))*(a + b/tan(e + f*x))**(n + S(-1)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**n/(a*f*m), x)
def replacement3494(m, n, x, d, f, a, e, b):
        # rubi.append(3494)
        return Dist(S(2)*d**S(2)/a, Int((d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(1)), x), x) + Simp(S(2)*d**S(2)*(d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(1))/(b*f*(m + S(-2))), x)
def replacement3495(m, n, x, d, f, a, e, b):
        # rubi.append(3495)
        return Dist(S(2)*d**S(2)/a, Int((d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(1)), x), x) + Simp(-S(2)*d**S(2)*(d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(1))/(b*f*(m + S(-2))), x)
def replacement3496(m, n, x, d, f, a, e, b):
        # rubi.append(3496)
        return Dist((a/d)**(S(2)*IntPart(n))*(d/cos(e + f*x))**(-S(2)*FracPart(n))*(a - b*tan(e + f*x))**FracPart(n)*(a + b*tan(e + f*x))**FracPart(n), Int((a - b*tan(e + f*x))**(-n), x), x)
def replacement3497(m, n, x, d, f, a, e, b):
        # rubi.append(3497)
        return Dist((a/d)**(S(2)*IntPart(n))*(d/sin(e + f*x))**(-S(2)*FracPart(n))*(a - b/tan(e + f*x))**FracPart(n)*(a + b/tan(e + f*x))**FracPart(n), Int((a - b/tan(e + f*x))**(-n), x), x)
def replacement3498(m, n, x, d, f, a, e, b):
        # rubi.append(3498)
        return Simp(S(2)*b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3499(m, n, x, d, f, a, e, b):
        # rubi.append(3499)
        return Simp(-S(2)*b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3500(m, n, x, d, f, a, e, b):
        # rubi.append(3500)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3501(m, n, x, d, f, a, e, b):
        # rubi.append(3501)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3502(x, d, f, a, e, b):
        # rubi.append(3502)
        return Dist(-S(4)*b*d**S(2)/f, Subst(Int(x**S(2)/(a**S(2) + d**S(2)*x**S(4)), x), x, sqrt(a + b*tan(e + f*x))/sqrt(d/cos(e + f*x))), x)
def replacement3503(x, d, f, a, e, b):
        # rubi.append(3503)
        return Dist(S(4)*b*d**S(2)/f, Subst(Int(x**S(2)/(a**S(2) + d**S(2)*x**S(4)), x), x, sqrt(a + b/tan(e + f*x))/sqrt(d/sin(e + f*x))), x)
def replacement3504(m, n, x, d, f, a, e, b):
        # rubi.append(3504)
        return -Dist(b**S(2)*(m + S(2)*n + S(-2))/(d**S(2)*m), Int((d/cos(e + f*x))**(m + S(2))*(a + b*tan(e + f*x))**(n + S(-2)), x), x) + Simp(S(2)*b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3505(m, n, x, d, f, a, e, b):
        # rubi.append(3505)
        return -Dist(b**S(2)*(m + S(2)*n + S(-2))/(d**S(2)*m), Int((d/sin(e + f*x))**(m + S(2))*(a + b/tan(e + f*x))**(n + S(-2)), x), x) + Simp(-S(2)*b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1))/(f*m), x)
def replacement3506(m, n, x, d, f, a, e, b):
        # rubi.append(3506)
        return Dist(a*(m + n)/(d**S(2)*m), Int((d/cos(e + f*x))**(m + S(2))*(a + b*tan(e + f*x))**(n + S(-1)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**n/(a*f*m), x)
def replacement3507(m, n, x, d, f, a, e, b):
        # rubi.append(3507)
        return Dist(a*(m + n)/(d**S(2)*m), Int((d/sin(e + f*x))**(m + S(2))*(a + b/tan(e + f*x))**(n + S(-1)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**n/(a*f*m), x)
def replacement3508(m, n, x, d, f, a, e, b):
        # rubi.append(3508)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3509(m, n, x, d, f, a, e, b):
        # rubi.append(3509)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3510(x, d, f, a, e, b):
        # rubi.append(3510)
        return Dist(d/(sqrt(a - b*tan(e + f*x))*sqrt(a + b*tan(e + f*x))*cos(e + f*x)), Int(sqrt(d/cos(e + f*x))*sqrt(a - b*tan(e + f*x)), x), x)
def replacement3511(x, d, f, a, e, b):
        # rubi.append(3511)
        return Dist(d/(sqrt(a - b/tan(e + f*x))*sqrt(a + b/tan(e + f*x))*sin(e + f*x)), Int(sqrt(d/sin(e + f*x))*sqrt(a - b/tan(e + f*x)), x), x)
def replacement3512(m, n, x, d, f, a, e, b):
        # rubi.append(3512)
        return -Dist(d**S(2)*(m + S(-2))/(b**S(2)*(m + S(2)*n)), Int((d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(2)), x), x) + Simp(S(2)*d**S(2)*(d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(1))/(b*f*(m + S(2)*n)), x)
def replacement3513(m, n, x, d, f, a, e, b):
        # rubi.append(3513)
        return -Dist(d**S(2)*(m + S(-2))/(b**S(2)*(m + S(2)*n)), Int((d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(2)), x), x) + Simp(-S(2)*d**S(2)*(d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(1))/(b*f*(m + S(2)*n)), x)
def replacement3514(m, n, x, d, f, a, e, b):
        # rubi.append(3514)
        return Dist(d**S(2)*(m + S(-2))/(a*(m + n + S(-1))), Int((d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(1)), x), x) + Simp(d**S(2)*(d/cos(e + f*x))**(m + S(-2))*(a + b*tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(-1))), x)
def replacement3515(m, n, x, d, f, a, e, b):
        # rubi.append(3515)
        return Dist(d**S(2)*(m + S(-2))/(a*(m + n + S(-1))), Int((d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(1)), x), x) - Simp(d**S(2)*(d/sin(e + f*x))**(m + S(-2))*(a + b/tan(e + f*x))**(n + S(1))/(b*f*(m + n + S(-1))), x)
def replacement3516(m, n, x, d, f, a, e, b):
        # rubi.append(3516)
        return Dist((m + n)/(a*(m + S(2)*n)), Int((d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(1)), x), x) + Simp(a*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**n/(b*f*(m + S(2)*n)), x)
def replacement3517(m, n, x, d, f, a, e, b):
        # rubi.append(3517)
        return Dist((m + n)/(a*(m + S(2)*n)), Int((d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(1)), x), x) - Simp(a*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**n/(b*f*(m + S(2)*n)), x)
def replacement3518(m, n, x, d, f, a, e, b):
        # rubi.append(3518)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3519(m, n, x, d, f, a, e, b):
        # rubi.append(3519)
        return Dist(a*(m + S(2)*n + S(-2))/(m + n + S(-1)), Int((d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3520(m, n, x, d, f, a, e, b):
        # rubi.append(3520)
        return Dist((m + n)/(a*(m + S(2)*n)), Int((d/cos(e + f*x))**m*(a + b*tan(e + f*x))**(n + S(1)), x), x) + Simp(a*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))**n/(b*f*(m + S(2)*n)), x)
def replacement3521(m, n, x, d, f, a, e, b):
        # rubi.append(3521)
        return Dist((m + n)/(a*(m + S(2)*n)), Int((d/sin(e + f*x))**m*(a + b/tan(e + f*x))**(n + S(1)), x), x) - Simp(a*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))**n/(b*f*(m + S(2)*n)), x)
def replacement3522(m, n, x, d, f, a, e, b):
        # rubi.append(3522)
        return Dist((d/cos(e + f*x))**m*(a - b*tan(e + f*x))**(-m/S(2))*(a + b*tan(e + f*x))**(-m/S(2)), Int((a - b*tan(e + f*x))**(m/S(2))*(a + b*tan(e + f*x))**(m/S(2) + n), x), x)
def replacement3523(m, n, x, d, f, a, e, b):
        # rubi.append(3523)
        return Dist((d/sin(e + f*x))**m*(a - b/tan(e + f*x))**(-m/S(2))*(a + b/tan(e + f*x))**(-m/S(2)), Int((a - b/tan(e + f*x))**(m/S(2))*(a + b/tan(e + f*x))**(m/S(2) + n), x), x)
def replacement3524(m, n, x, f, a, e, b):
        # rubi.append(3524)
        return Dist(S(1)/(b*f), Subst(Int((S(1) + x**S(2)/b**S(2))**(m/S(2) + S(-1))*(a + x)**n, x), x, b*tan(e + f*x)), x)
def replacement3525(m, n, x, f, a, e, b):
        # rubi.append(3525)
        return -Dist(S(1)/(b*f), Subst(Int((S(1) + x**S(2)/b**S(2))**(m/S(2) + S(-1))*(a + x)**n, x), x, b/tan(e + f*x)), x)
def replacement3526(x, f, a, e, b):
        # rubi.append(3526)
        return Simp(b**S(2)*atanh(sin(e + f*x))/f, x) + Simp((a**S(2) - b**S(2))*sin(e + f*x)/f, x) - Simp(S(2)*a*b*cos(e + f*x)/f, x)
def replacement3527(x, f, a, e, b):
        # rubi.append(3527)
        return -Simp(b**S(2)*atanh(cos(e + f*x))/f, x) - Simp((a**S(2) - b**S(2))*cos(e + f*x)/f, x) + Simp(S(2)*a*b*sin(e + f*x)/f, x)
def replacement3528(m, x, d, f, a, e, b):
        # rubi.append(3528)
        return Dist(S(1)/(m + S(1)), Int((d/cos(e + f*x))**m*(a**S(2)*(m + S(1)) + a*b*(m + S(2))*tan(e + f*x) - b**S(2)), x), x) + Simp(b*(d/cos(e + f*x))**m*(a + b*tan(e + f*x))/(f*(m + S(1))), x)
def replacement3529(m, x, d, f, a, e, b):
        # rubi.append(3529)
        return Dist(S(1)/(m + S(1)), Int((d/sin(e + f*x))**m*(a**S(2)*(m + S(1)) + a*b*(m + S(2))/tan(e + f*x) - b**S(2)), x), x) - Simp(b*(d/sin(e + f*x))**m*(a + b/tan(e + f*x))/(f*(m + S(1))), x)
def replacement3530(x, f, a, e, b):
        # rubi.append(3530)
        return -Dist(S(1)/f, Subst(Int(S(1)/(a**S(2) + b**S(2) - x**S(2)), x), x, (-a*tan(e + f*x) + b)*cos(e + f*x)), x)
def replacement3531(x, f, a, e, b):
        # rubi.append(3531)
        return Dist(S(1)/f, Subst(Int(S(1)/(a**S(2) + b**S(2) - x**S(2)), x), x, (-a/tan(e + f*x) + b)*sin(e + f*x)), x)
def replacement3532(m, x, d, f, a, e, b):
        # rubi.append(3532)
        return -Dist(d**S(2)/b**S(2), Int((d/cos(e + f*x))**(m + S(-2))*(a - b*tan(e + f*x)), x), x) + Dist(d**S(2)*(a**S(2) + b**S(2))/b**S(2), Int((d/cos(e + f*x))**(m + S(-2))/(a + b*tan(e + f*x)), x), x)
def replacement3533(m, x, d, f, a, e, b):
        # rubi.append(3533)
        return -Dist(d**S(2)/b**S(2), Int((d/sin(e + f*x))**(m + S(-2))*(a - b/tan(e + f*x)), x), x) + Dist(d**S(2)*(a**S(2) + b**S(2))/b**S(2), Int((d/sin(e + f*x))**(m + S(-2))/(a + b/tan(e + f*x)), x), x)
def replacement3534(m, x, d, f, a, e, b):
        # rubi.append(3534)
        return Dist(b**S(2)/(d**S(2)*(a**S(2) + b**S(2))), Int((d/cos(e + f*x))**(m + S(2))/(a + b*tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((d/cos(e + f*x))**m*(a - b*tan(e + f*x)), x), x)
def replacement3535(m, x, d, f, a, e, b):
        # rubi.append(3535)
        return Dist(b**S(2)/(d**S(2)*(a**S(2) + b**S(2))), Int((d/sin(e + f*x))**(m + S(2))/(a + b/tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((d/sin(e + f*x))**m*(a - b/tan(e + f*x)), x), x)
def replacement3536(m, n, x, d, f, a, e, b):
        # rubi.append(3536)
        return Dist(d**(S(2)*IntPart(m/S(2)))*(d/cos(e + f*x))**(S(2)*FracPart(m/S(2)))*(cos(e + f*x)**(S(-2)))**(-FracPart(m/S(2)))/(b*f), Subst(Int((S(1) + x**S(2)/b**S(2))**(m/S(2) + S(-1))*(a + x)**n, x), x, b*tan(e + f*x)), x)
def replacement3537(m, n, x, d, f, a, e, b):
        # rubi.append(3537)
        return -Dist(d**(S(2)*IntPart(m/S(2)))*(d/sin(e + f*x))**(S(2)*FracPart(m/S(2)))*(sin(e + f*x)**(S(-2)))**(-FracPart(m/S(2)))/(b*f), Subst(Int((S(1) + x**S(2)/b**S(2))**(m/S(2) + S(-1))*(a + x)**n, x), x, b/tan(e + f*x)), x)
def replacement3538(x, d, f, a, e, b):
        # rubi.append(3538)
        return Dist(-S(4)*b/f, Subst(Int(x**S(2)/(a**S(2)*d**S(2) + x**S(4)), x), x, sqrt(d*cos(e + f*x))*sqrt(a + b*tan(e + f*x))), x)
def replacement3539(x, d, f, a, e, b):
        # rubi.append(3539)
        return Dist(S(4)*b/f, Subst(Int(x**S(2)/(a**S(2)*d**S(2) + x**S(4)), x), x, sqrt(d*sin(e + f*x))*sqrt(a + b/tan(e + f*x))), x)
def replacement3540(x, d, f, a, e, b):
        # rubi.append(3540)
        return Dist(S(1)/(d*sqrt(a - b*tan(e + f*x))*sqrt(a + b*tan(e + f*x))*cos(e + f*x)), Int(sqrt(a - b*tan(e + f*x))/sqrt(d*cos(e + f*x)), x), x)
def replacement3541(x, d, f, a, e, b):
        # rubi.append(3541)
        return Dist(S(1)/(d*sqrt(a - b/tan(e + f*x))*sqrt(a + b/tan(e + f*x))*sin(e + f*x)), Int(sqrt(a - b/tan(e + f*x))/sqrt(d*sin(e + f*x)), x), x)
def replacement3542(m, n, x, d, f, a, e, b):
        # rubi.append(3542)
        return Dist((d/cos(e + f*x))**m*(d*cos(e + f*x))**m, Int((d/cos(e + f*x))**(-m)*(a + b*tan(e + f*x))**n, x), x)
def replacement3543(m, n, x, d, f, a, e, b):
        # rubi.append(3543)
        return Dist((d/sin(e + f*x))**m*(d*sin(e + f*x))**m, Int((d/sin(e + f*x))**(-m)*(a + b/tan(e + f*x))**n, x), x)
def replacement3544(m, n, x, f, a, e, b):
        # rubi.append(3544)
        return Dist(b/f, Subst(Int(x**m*(a + x)**n*(b**S(2) + x**S(2))**(-m/S(2) + S(-1)), x), x, b*tan(e + f*x)), x)
def replacement3545(m, n, x, f, a, e, b):
        # rubi.append(3545)
        return -Dist(b/f, Subst(Int(x**m*(a + x)**n*(b**S(2) + x**S(2))**(-m/S(2) + S(-1)), x), x, b/tan(e + f*x)), x)
def replacement3546(m, n, x, f, a, e, b):
        # rubi.append(3546)
        return Int((a + b*tan(e + f*x))**n*sin(e + f*x)**m, x)
def replacement3547(m, n, x, f, a, e, b):
        # rubi.append(3547)
        return Int((a + b/tan(e + f*x))**n*cos(e + f*x)**m, x)
def replacement3548(m, n, x, f, a, e, b):
        # rubi.append(3548)
        return Int((a*cos(e + f*x) + b*sin(e + f*x))**n*sin(e + f*x)**m*cos(e + f*x)**(-n), x)
def replacement3549(m, n, x, f, a, e, b):
        # rubi.append(3549)
        return Int((a*sin(e + f*x) + b*cos(e + f*x))**n*sin(e + f*x)**(-n)*cos(e + f*x)**m, x)
def replacement3550(m, n, x, d, f, a, e, b):
        # rubi.append(3550)
        return Dist((sin(e + f*x)/d)**FracPart(m)*(d/sin(e + f*x))**FracPart(m), Int((sin(e + f*x)/d)**(-m)*(a + b*tan(e + f*x))**n, x), x)
def replacement3551(m, n, x, d, f, a, e, b):
        # rubi.append(3551)
        return Dist((cos(e + f*x)/d)**FracPart(m)*(d/cos(e + f*x))**FracPart(m), Int((cos(e + f*x)/d)**(-m)*(a + b/tan(e + f*x))**n, x), x)
def replacement3552(m, n, x, f, a, p, e, b):
        # rubi.append(3552)
        return Int((a*cos(e + f*x) + b*sin(e + f*x))**n*sin(e + f*x)**p*cos(e + f*x)**(m - n), x)
def replacement3553(m, n, x, f, a, p, e, b):
        # rubi.append(3553)
        return Int((a*sin(e + f*x) + b*cos(e + f*x))**n*sin(e + f*x)**(m - n)*cos(e + f*x)**p, x)
def replacement3554(c, m, n, x, d, f, a, e, b):
        # rubi.append(3554)
        return Dist(a**m*c**m, Int((c + d*tan(e + f*x))**(-m + n)*(S(1)/cos(e + f*x))**(S(2)*m), x), x)
def replacement3555(c, m, n, x, d, f, a, e, b):
        # rubi.append(3555)
        return Dist(a**m*c**m, Int((c + d/tan(e + f*x))**(-m + n)*(S(1)/sin(e + f*x))**(S(2)*m), x), x)
def replacement3556(c, m, n, x, d, f, a, e, b):
        # rubi.append(3556)
        return Dist(a*c/f, Subst(Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1)), x), x, tan(e + f*x)), x)
def replacement3557(c, m, n, x, d, f, a, e, b):
        # rubi.append(3557)
        return -Dist(a*c/f, Subst(Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1)), x), x, S(1)/tan(e + f*x)), x)
def replacement3558(c, x, d, f, a, e, b):
        # rubi.append(3558)
        return Simp(x*(a*c - b*d), x) + Simp(b*d*tan(e + f*x)/f, x)
def replacement3559(c, x, d, f, a, e, b):
        # rubi.append(3559)
        return Simp(x*(a*c - b*d), x) - Simp(b*d/(f*tan(e + f*x)), x)
def replacement3560(c, x, d, f, a, e, b):
        # rubi.append(3560)
        return Dist(a*d + b*c, Int(tan(e + f*x), x), x) + Simp(x*(a*c - b*d), x) + Simp(b*d*tan(e + f*x)/f, x)
def replacement3561(c, x, d, f, a, e, b):
        # rubi.append(3561)
        return Dist(a*d + b*c, Int(S(1)/tan(e + f*x), x), x) + Simp(x*(a*c - b*d), x) - Simp(b*d/(f*tan(e + f*x)), x)
def replacement3562(c, m, x, d, f, a, e, b):
        # rubi.append(3562)
        return Dist((a*d + b*c)/(S(2)*a*b), Int((a + b*tan(e + f*x))**(m + S(1)), x), x) - Simp((a + b*tan(e + f*x))**m*(-a*d + b*c)/(S(2)*a*f*m), x)
def replacement3563(c, m, x, d, f, a, e, b):
        # rubi.append(3563)
        return Dist((a*d + b*c)/(S(2)*a*b), Int((a + b/tan(e + f*x))**(m + S(1)), x), x) + Simp((a + b/tan(e + f*x))**m*(-a*d + b*c)/(S(2)*a*f*m), x)
def replacement3564(c, m, x, d, f, a, e, b):
        # rubi.append(3564)
        return Dist((a*d + b*c)/b, Int((a + b*tan(e + f*x))**m, x), x) + Simp(d*(a + b*tan(e + f*x))**m/(f*m), x)
def replacement3565(c, m, x, d, f, a, e, b):
        # rubi.append(3565)
        return Dist((a*d + b*c)/b, Int((a + b/tan(e + f*x))**m, x), x) - Simp(d*(a + b/tan(e + f*x))**m/(f*m), x)
def replacement3566(c, m, x, d, f, a, e, b):
        # rubi.append(3566)
        return Int((a + b*tan(e + f*x))**(m + S(-1))*Simp(a*c - b*d + (a*d + b*c)*tan(e + f*x), x), x) + Simp(d*(a + b*tan(e + f*x))**m/(f*m), x)
def replacement3567(c, m, x, d, f, a, e, b):
        # rubi.append(3567)
        return Int((a + b/tan(e + f*x))**(m + S(-1))*Simp(a*c - b*d + (a*d + b*c)/tan(e + f*x), x), x) - Simp(d*(a + b/tan(e + f*x))**m/(f*m), x)
def replacement3568(c, m, x, d, f, a, e, b):
        # rubi.append(3568)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(a*c + b*d - (-a*d + b*c)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(-a*d + b*c)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3569(c, m, x, d, f, a, e, b):
        # rubi.append(3569)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(a*c + b*d - (-a*d + b*c)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(-a*d + b*c)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3570(c, x, d, f, a, e, b):
        # rubi.append(3570)
        return Simp(c*log(RemoveContent(a*cos(e + f*x) + b*sin(e + f*x), x))/(b*f), x)
def replacement3571(c, x, d, f, a, e, b):
        # rubi.append(3571)
        return -Simp(c*log(RemoveContent(a*sin(e + f*x) + b*cos(e + f*x), x))/(b*f), x)
def replacement3572(c, x, d, f, a, e, b):
        # rubi.append(3572)
        return Dist((-a*d + b*c)/(a**S(2) + b**S(2)), Int((-a*tan(e + f*x) + b)/(a + b*tan(e + f*x)), x), x) + Simp(x*(a*c + b*d)/(a**S(2) + b**S(2)), x)
def replacement3573(c, x, d, f, a, e, b):
        # rubi.append(3573)
        return Dist((-a*d + b*c)/(a**S(2) + b**S(2)), Int((-a/tan(e + f*x) + b)/(a + b/tan(e + f*x)), x), x) + Simp(x*(a*c + b*d)/(a**S(2) + b**S(2)), x)
def replacement3574(c, x, d, f, e, b):
        # rubi.append(3574)
        return Dist(-S(2)*d**S(2)/f, Subst(Int(S(1)/(b*x**S(2) + S(2)*c*d), x), x, (c - d*tan(e + f*x))/sqrt(b*tan(e + f*x))), x)
def replacement3575(c, x, d, f, e, b):
        # rubi.append(3575)
        return Dist(S(2)*d**S(2)/f, Subst(Int(S(1)/(b*x**S(2) + S(2)*c*d), x), x, (c - d/tan(e + f*x))/sqrt(b/tan(e + f*x))), x)
def replacement3576(c, x, d, f, e, b):
        # rubi.append(3576)
        return Dist(S(2)*c**S(2)/f, Subst(Int(S(1)/(b*c - d*x**S(2)), x), x, sqrt(b*tan(e + f*x))), x)
def replacement3577(c, x, d, f, e, b):
        # rubi.append(3577)
        return Dist(-S(2)*c**S(2)/f, Subst(Int(S(1)/(b*c - d*x**S(2)), x), x, sqrt(b/tan(e + f*x))), x)
def replacement3578(c, x, d, f, e, b):
        # rubi.append(3578)
        return Dist(S(2)/f, Subst(Int((b*c + d*x**S(2))/(b**S(2) + x**S(4)), x), x, sqrt(b*tan(e + f*x))), x)
def replacement3579(c, x, d, f, e, b):
        # rubi.append(3579)
        return Dist(-S(2)/f, Subst(Int((b*c + d*x**S(2))/(b**S(2) + x**S(4)), x), x, sqrt(b/tan(e + f*x))), x)
def replacement3580(c, x, d, f, a, e, b):
        # rubi.append(3580)
        return Dist(-S(2)*d**S(2)/f, Subst(Int(S(1)/(-S(4)*a*d**S(2) + S(2)*b*c*d + x**S(2)), x), x, (-S(2)*a*d + b*c - b*d*tan(e + f*x))/sqrt(a + b*tan(e + f*x))), x)
def replacement3581(c, x, d, f, a, e, b):
        # rubi.append(3581)
        return Dist(S(2)*d**S(2)/f, Subst(Int(S(1)/(-S(4)*a*d**S(2) + S(2)*b*c*d + x**S(2)), x), x, (-S(2)*a*d + b*c - b*d/tan(e + f*x))/sqrt(a + b/tan(e + f*x))), x)

def With3582(c, x, d, f, a, e, b):
        q = Rt(a**S(2) + b**S(2), S(2))
        # rubi.append(3582)
        return -Dist(S(1)/(S(2)*q), Int((a*c + b*d - c*q + (-a*d + b*c - d*q)*tan(e + f*x))/sqrt(a + b*tan(e + f*x)), x), x) + Dist(S(1)/(S(2)*q), Int((a*c + b*d + c*q + (-a*d + b*c + d*q)*tan(e + f*x))/sqrt(a + b*tan(e + f*x)), x), x)

def With3583(c, x, d, f, a, e, b):
        q = Rt(a**S(2) + b**S(2), S(2))
        # rubi.append(3583)
        return -Dist(S(1)/(S(2)*q), Int((a*c + b*d - c*q + (-a*d + b*c - d*q)/tan(e + f*x))/sqrt(a + b/tan(e + f*x)), x), x) + Dist(S(1)/(S(2)*q), Int((a*c + b*d + c*q + (-a*d + b*c + d*q)/tan(e + f*x))/sqrt(a + b/tan(e + f*x)), x), x)
def replacement3584(c, m, x, d, f, a, e, b):
        # rubi.append(3584)
        return Dist(c*d/f, Subst(Int((a + b*x/d)**m/(c*x + d**S(2)), x), x, d*tan(e + f*x)), x)
def replacement3585(c, m, x, d, f, a, e, b):
        # rubi.append(3585)
        return -Dist(c*d/f, Subst(Int((a + b*x/d)**m/(c*x + d**S(2)), x), x, d/tan(e + f*x)), x)
def replacement3586(c, m, x, d, f, e, b):
        # rubi.append(3586)
        return Dist(c, Int((b*tan(e + f*x))**m, x), x) + Dist(d/b, Int((b*tan(e + f*x))**(m + S(1)), x), x)
def replacement3587(c, m, x, d, f, e, b):
        # rubi.append(3587)
        return Dist(c, Int((b/tan(e + f*x))**m, x), x) + Dist(d/b, Int((b/tan(e + f*x))**(m + S(1)), x), x)
def replacement3588(c, m, x, d, f, a, e, b):
        # rubi.append(3588)
        return Dist(c/S(2) - I*d/S(2), Int((a + b*tan(e + f*x))**m*(I*tan(e + f*x) + S(1)), x), x) + Dist(c/S(2) + I*d/S(2), Int((a + b*tan(e + f*x))**m*(-I*tan(e + f*x) + S(1)), x), x)
def replacement3589(c, m, x, d, f, a, e, b):
        # rubi.append(3589)
        return Dist(c/S(2) - I*d/S(2), Int((S(1) + I/tan(e + f*x))*(a + b/tan(e + f*x))**m, x), x) + Dist(c/S(2) + I*d/S(2), Int((S(1) - I/tan(e + f*x))*(a + b/tan(e + f*x))**m, x), x)
def replacement3590(c, m, x, d, f, a, e, b):
        # rubi.append(3590)
        return Dist(S(1)/(S(2)*a**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(a*c**S(2) + a*d**S(2) - S(2)*b*c*d - S(2)*b*d**S(2)*tan(e + f*x), x), x), x) - Simp(b*(a + b*tan(e + f*x))**m*(a*c + b*d)**S(2)/(S(2)*a**S(3)*f*m), x)
def replacement3591(c, m, x, d, f, a, e, b):
        # rubi.append(3591)
        return Dist(S(1)/(S(2)*a**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(a*c**S(2) + a*d**S(2) - S(2)*b*c*d - S(2)*b*d**S(2)/tan(e + f*x), x), x), x) + Simp(b*(a + b/tan(e + f*x))**m*(a*c + b*d)**S(2)/(S(2)*a**S(3)*f*m), x)
def replacement3592(c, x, d, f, a, e, b):
        # rubi.append(3592)
        return Dist((-a*d + b*c)**S(2)/b**S(2), Int(S(1)/(a + b*tan(e + f*x)), x), x) + Dist(d**S(2)/b, Int(tan(e + f*x), x), x) + Simp(d*x*(-a*d + S(2)*b*c)/b**S(2), x)
def replacement3593(c, x, d, f, a, e, b):
        # rubi.append(3593)
        return Dist((-a*d + b*c)**S(2)/b**S(2), Int(S(1)/(a + b/tan(e + f*x)), x), x) + Dist(d**S(2)/b, Int(S(1)/tan(e + f*x), x), x) + Simp(d*x*(-a*d + S(2)*b*c)/b**S(2), x)
def replacement3594(c, m, x, d, f, a, e, b):
        # rubi.append(3594)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(a*c**S(2) - a*d**S(2) + S(2)*b*c*d - (-S(2)*a*c*d + b*c**S(2) - b*d**S(2))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(-a*d + b*c)**S(2)/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3595(c, m, x, d, f, a, e, b):
        # rubi.append(3595)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(a*c**S(2) - a*d**S(2) + S(2)*b*c*d - (-S(2)*a*c*d + b*c**S(2) - b*d**S(2))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(-a*d + b*c)**S(2)/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3596(c, m, x, d, f, a, e, b):
        # rubi.append(3596)
        return Int((a + b*tan(e + f*x))**m*Simp(c**S(2) + S(2)*c*d*tan(e + f*x) - d**S(2), x), x) + Simp(d**S(2)*(a + b*tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3597(c, m, x, d, f, a, e, b):
        # rubi.append(3597)
        return Int((a + b/tan(e + f*x))**m*Simp(c**S(2) + S(2)*c*d/tan(e + f*x) - d**S(2), x), x) - Simp(d**S(2)*(a + b/tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3598(c, x, d, f, a, e, b):
        # rubi.append(3598)
        return Dist(-S(2)*a*b/f, Subst(Int(S(1)/(-S(2)*a**S(2)*x**S(2) + a*c - b*d), x), x, sqrt(c + d*tan(e + f*x))/sqrt(a + b*tan(e + f*x))), x)
def replacement3599(c, x, d, f, a, e, b):
        # rubi.append(3599)
        return Dist(S(2)*a*b/f, Subst(Int(S(1)/(-S(2)*a**S(2)*x**S(2) + a*c - b*d), x), x, sqrt(c + d/tan(e + f*x))/sqrt(a + b/tan(e + f*x))), x)
def replacement3600(c, m, n, x, d, f, a, e, b):
        # rubi.append(3600)
        return Dist(S(2)*a**S(2)/(a*c - b*d), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1)), x), x) + Simp(a*b*(a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))/(f*(m + S(-1))*(a*c - b*d)), x)
def replacement3601(c, m, n, x, d, f, a, e, b):
        # rubi.append(3601)
        return Dist(S(2)*a**S(2)/(a*c - b*d), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1)), x), x) - Simp(a*b*(a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))/(f*(m + S(-1))*(a*c - b*d)), x)
def replacement3602(c, m, n, x, d, f, a, e, b):
        # rubi.append(3602)
        return -Dist((a*c - b*d)/(S(2)*b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-1)), x), x) + Simp(a*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n/(S(2)*b*f*m), x)
def replacement3603(c, m, n, x, d, f, a, e, b):
        # rubi.append(3603)
        return -Dist((a*c - b*d)/(S(2)*b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-1)), x), x) - Simp(a*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n/(S(2)*b*f*m), x)
def replacement3604(c, m, n, x, d, f, a, e, b):
        # rubi.append(3604)
        return Dist(S(1)/(S(2)*a), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n, x), x) + Simp(a*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3605(c, m, n, x, d, f, a, e, b):
        # rubi.append(3605)
        return Dist(S(1)/(S(2)*a), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n, x), x) - Simp(a*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3606(c, m, n, x, d, f, a, e, b):
        # rubi.append(3606)
        return Dist(a/(a*c - b*d), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1)), x), x) - Simp(d*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(f*m*(c**S(2) + d**S(2))), x)
def replacement3607(c, m, n, x, d, f, a, e, b):
        # rubi.append(3607)
        return Dist(a/(a*c - b*d), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1)), x), x) + Simp(d*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(f*m*(c**S(2) + d**S(2))), x)
def replacement3608(c, n, x, d, f, a, e, b):
        # rubi.append(3608)
        return Dist(S(1)/(S(2)*a*(-a*d + b*c)), Int((c + d*tan(e + f*x))**(n + S(-1))*Simp(a*c*d*(n + S(-1)) + b*c**S(2) + b*d**S(2)*n - d*(n + S(-1))*(-a*d + b*c)*tan(e + f*x), x), x), x) - Simp((c + d*tan(e + f*x))**n*(a*c + b*d)/(S(2)*f*(a + b*tan(e + f*x))*(-a*d + b*c)), x)
def replacement3609(c, n, x, d, f, a, e, b):
        # rubi.append(3609)
        return Dist(S(1)/(S(2)*a*(-a*d + b*c)), Int((c + d/tan(e + f*x))**(n + S(-1))*Simp(a*c*d*(n + S(-1)) + b*c**S(2) + b*d**S(2)*n - d*(n + S(-1))*(-a*d + b*c)/tan(e + f*x), x), x), x) + Simp((c + d/tan(e + f*x))**n*(a*c + b*d)/(S(2)*f*(a + b/tan(e + f*x))*(-a*d + b*c)), x)
def replacement3610(c, n, x, d, f, a, e, b):
        # rubi.append(3610)
        return Dist(S(1)/(S(2)*a**S(2)), Int((c + d*tan(e + f*x))**(n + S(-2))*Simp(a*c**S(2) + a*d**S(2)*(n + S(-1)) - b*c*d*n - d*(a*c*(n + S(-2)) + b*d*n)*tan(e + f*x), x), x), x) + Simp((c + d*tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(S(2)*a*f*(a + b*tan(e + f*x))), x)
def replacement3611(c, n, x, d, f, a, e, b):
        # rubi.append(3611)
        return Dist(S(1)/(S(2)*a**S(2)), Int((c + d/tan(e + f*x))**(n + S(-2))*Simp(a*c**S(2) + a*d**S(2)*(n + S(-1)) - b*c*d*n - d*(a*c*(n + S(-2)) + b*d*n)/tan(e + f*x), x), x), x) - Simp((c + d/tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(S(2)*a*f*(a + b/tan(e + f*x))), x)
def replacement3612(c, x, d, f, a, e, b):
        # rubi.append(3612)
        return Dist(b/(-a*d + b*c), Int(S(1)/(a + b*tan(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d*tan(e + f*x)), x), x)
def replacement3613(c, x, d, f, a, e, b):
        # rubi.append(3613)
        return Dist(b/(-a*d + b*c), Int(S(1)/(a + b/tan(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d/tan(e + f*x)), x), x)
def replacement3614(c, n, x, d, f, a, e, b):
        # rubi.append(3614)
        return Dist(S(1)/(S(2)*a*(-a*d + b*c)), Int((c + d*tan(e + f*x))**n*Simp(a*d*(n + S(-1)) + b*c - b*d*n*tan(e + f*x), x), x), x) - Simp(a*(c + d*tan(e + f*x))**(n + S(1))/(S(2)*f*(a + b*tan(e + f*x))*(-a*d + b*c)), x)
def replacement3615(c, n, x, d, f, a, e, b):
        # rubi.append(3615)
        return Dist(S(1)/(S(2)*a*(-a*d + b*c)), Int((c + d/tan(e + f*x))**n*Simp(a*d*(n + S(-1)) + b*c - b*d*n/tan(e + f*x), x), x), x) + Simp(a*(c + d/tan(e + f*x))**(n + S(1))/(S(2)*f*(a + b/tan(e + f*x))*(-a*d + b*c)), x)
def replacement3616(c, m, n, x, d, f, a, e, b):
        # rubi.append(3616)
        return Dist(a/(d*(n + S(1))*(a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))*Simp(b*(-a*d*(m - S(2)*n + S(-4)) + b*c*(m + S(-2))) + (-a**S(2)*d*(m + n + S(-1)) + a*b*c*(m + S(-2)) + b**S(2)*d*(n + S(1)))*tan(e + f*x), x), x), x) - Simp(a**S(2)*(a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))*(-a*d + b*c)/(d*f*(n + S(1))*(a*d + b*c)), x)
def replacement3617(c, m, n, x, d, f, a, e, b):
        # rubi.append(3617)
        return Dist(a/(d*(n + S(1))*(a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))*Simp(b*(-a*d*(m - S(2)*n + S(-4)) + b*c*(m + S(-2))) + (-a**S(2)*d*(m + n + S(-1)) + a*b*c*(m + S(-2)) + b**S(2)*d*(n + S(1)))/tan(e + f*x), x), x), x) + Simp(a**S(2)*(a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))*(-a*d + b*c)/(d*f*(n + S(1))*(a*d + b*c)), x)
def replacement3618(c, x, d, f, a, e, b):
        # rubi.append(3618)
        return Dist(S(2)*a**S(2)/(a*c - b*d), Int(sqrt(a + b*tan(e + f*x)), x), x) - Dist((a*(c**S(2) - d**S(2)) + S(2)*b*c*d)/(a*(c**S(2) + d**S(2))), Int((a - b*tan(e + f*x))*sqrt(a + b*tan(e + f*x))/(c + d*tan(e + f*x)), x), x)
def replacement3619(c, x, d, f, a, e, b):
        # rubi.append(3619)
        return Dist(S(2)*a**S(2)/(a*c - b*d), Int(sqrt(a + b/tan(e + f*x)), x), x) - Dist((a*(c**S(2) - d**S(2)) + S(2)*b*c*d)/(a*(c**S(2) + d**S(2))), Int((a - b/tan(e + f*x))*sqrt(a + b/tan(e + f*x))/(c + d/tan(e + f*x)), x), x)
def replacement3620(c, x, d, f, a, e, b):
        # rubi.append(3620)
        return Dist(S(2)*a, Int(sqrt(a + b*tan(e + f*x))/sqrt(c + d*tan(e + f*x)), x), x) + Dist(b/a, Int(sqrt(a + b*tan(e + f*x))*(a*tan(e + f*x) + b)/sqrt(c + d*tan(e + f*x)), x), x)
def replacement3621(c, x, d, f, a, e, b):
        # rubi.append(3621)
        return Dist(S(2)*a, Int(sqrt(a + b/tan(e + f*x))/sqrt(c + d/tan(e + f*x)), x), x) + Dist(b/a, Int(sqrt(a + b/tan(e + f*x))*(a/tan(e + f*x) + b)/sqrt(c + d/tan(e + f*x)), x), x)
def replacement3622(c, m, n, x, d, f, a, e, b):
        # rubi.append(3622)
        return Dist(a/(d*(m + n + S(-1))), Int((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**n*Simp(a*d*(m + S(2)*n) + b*c*(m + S(-2)) + (a*c*(m + S(-2)) + b*d*(S(3)*m + S(2)*n + S(-4)))*tan(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(-1))), x)
def replacement3623(c, m, n, x, d, f, a, e, b):
        # rubi.append(3623)
        return Dist(a/(d*(m + n + S(-1))), Int((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**n*Simp(a*d*(m + S(2)*n) + b*c*(m + S(-2)) + (a*c*(m + S(-2)) + b*d*(S(3)*m + S(2)*n + S(-4)))/tan(e + f*x), x), x), x) - Simp(b**S(2)*(a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(-1))), x)
def replacement3624(c, m, x, d, f, a, e, b):
        # rubi.append(3624)
        return Dist(S(1)/(S(4)*a**S(2)*m), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(S(2)*a*c*m + a*d*(S(2)*m + S(1))*tan(e + f*x) + b*d, x)/sqrt(c + d*tan(e + f*x)), x), x) - Simp(b*(a + b*tan(e + f*x))**m*sqrt(c + d*tan(e + f*x))/(S(2)*a*f*m), x)
def replacement3625(c, m, x, d, f, a, e, b):
        # rubi.append(3625)
        return Dist(S(1)/(S(4)*a**S(2)*m), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(S(2)*a*c*m + a*d*(S(2)*m + S(1))/tan(e + f*x) + b*d, x)/sqrt(c + d/tan(e + f*x)), x), x) + Simp(b*(a + b/tan(e + f*x))**m*sqrt(c + d/tan(e + f*x))/(S(2)*a*f*m), x)
def replacement3626(c, m, n, x, d, f, a, e, b):
        # rubi.append(3626)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-2))*Simp(c*(a*c*m + b*d*(n + S(-1))) - d*(-a*c*(m + n + S(-1)) + b*d*(m - n + S(1)))*tan(e + f*x) - d*(a*d*(n + S(-1)) + b*c*m), x), x), x) - Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(S(2)*a*f*m), x)
def replacement3627(c, m, n, x, d, f, a, e, b):
        # rubi.append(3627)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-2))*Simp(c*(a*c*m + b*d*(n + S(-1))) - d*(-a*c*(m + n + S(-1)) + b*d*(m - n + S(1)))/tan(e + f*x) - d*(a*d*(n + S(-1)) + b*c*m), x), x), x) + Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(S(2)*a*f*m), x)
def replacement3628(c, m, n, x, d, f, a, e, b):
        # rubi.append(3628)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(-a*d*(S(2)*m + n + S(1)) + b*c*m + b*d*(m + n + S(1))*tan(e + f*x), x), x), x) + Simp(a*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3629(c, m, n, x, d, f, a, e, b):
        # rubi.append(3629)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(-a*d*(S(2)*m + n + S(1)) + b*c*m + b*d*(m + n + S(1))/tan(e + f*x), x), x), x) - Simp(a*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3630(c, m, n, x, d, f, a, e, b):
        # rubi.append(3630)
        return -Dist(S(1)/(a*(m + n + S(-1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(-2))*Simp(-a*c**S(2)*(m + n + S(-1)) + d*(-a*c*(m + S(2)*n + S(-2)) + b*d*m)*tan(e + f*x) + d*(a*d*(n + S(-1)) + b*c*m), x), x), x) + Simp(d*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3631(c, m, n, x, d, f, a, e, b):
        # rubi.append(3631)
        return -Dist(S(1)/(a*(m + n + S(-1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(-2))*Simp(-a*c**S(2)*(m + n + S(-1)) + d*(-a*c*(m + S(2)*n + S(-2)) + b*d*m)/tan(e + f*x) + d*(a*d*(n + S(-1)) + b*c*m), x), x), x) - Simp(d*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(-1))/(f*(m + n + S(-1))), x)
def replacement3632(c, m, n, x, d, f, a, e, b):
        # rubi.append(3632)
        return -Dist(S(1)/(a*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*Simp(-a*c*(n + S(1)) + a*d*(m + n + S(1))*tan(e + f*x) + b*d*m, x), x), x) + Simp(d*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3633(c, m, n, x, d, f, a, e, b):
        # rubi.append(3633)
        return -Dist(S(1)/(a*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*Simp(-a*c*(n + S(1)) + a*d*(m + n + S(1))/tan(e + f*x) + b*d*m, x), x), x) - Simp(d*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3634(c, m, x, d, f, a, e, b):
        # rubi.append(3634)
        return Dist(a/(a*c - b*d), Int((a + b*tan(e + f*x))**m, x), x) - Dist(d/(a*c - b*d), Int((a + b*tan(e + f*x))**m*(a*tan(e + f*x) + b)/(c + d*tan(e + f*x)), x), x)
def replacement3635(c, m, x, d, f, a, e, b):
        # rubi.append(3635)
        return Dist(a/(a*c - b*d), Int((a + b/tan(e + f*x))**m, x), x) - Dist(d/(a*c - b*d), Int((a + b/tan(e + f*x))**m*(a/tan(e + f*x) + b)/(c + d/tan(e + f*x)), x), x)
def replacement3636(c, x, d, f, a, e, b):
        # rubi.append(3636)
        return Dist(d/a, Int(sqrt(a + b*tan(e + f*x))*(a*tan(e + f*x) + b)/sqrt(c + d*tan(e + f*x)), x), x) + Dist((a*c - b*d)/a, Int(sqrt(a + b*tan(e + f*x))/sqrt(c + d*tan(e + f*x)), x), x)
def replacement3637(c, x, d, f, a, e, b):
        # rubi.append(3637)
        return Dist(d/a, Int(sqrt(a + b/tan(e + f*x))*(a/tan(e + f*x) + b)/sqrt(c + d/tan(e + f*x)), x), x) + Dist((a*c - b*d)/a, Int(sqrt(a + b/tan(e + f*x))/sqrt(c + d/tan(e + f*x)), x), x)
def replacement3638(c, m, n, x, d, f, a, e, b):
        # rubi.append(3638)
        return Dist(a*b/f, Subst(Int((a + x)**(m + S(-1))*(c + d*x/b)**n/(a*x + b**S(2)), x), x, b*tan(e + f*x)), x)
def replacement3639(c, m, n, x, d, f, a, e, b):
        # rubi.append(3639)
        return -Dist(a*b/f, Subst(Int((a + x)**(m + S(-1))*(c + d*x/b)**n/(a*x + b**S(2)), x), x, b/tan(e + f*x)), x)
def replacement3640(c, m, n, x, d, f, a, e, b):
        # rubi.append(3640)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-3))*(c + d*tan(e + f*x))**(n + S(1))*Simp(a**S(2)*d*(-a*c*(n + S(1)) + b*d*(m + S(-2))) + b*(-S(2)*a*d + b*c)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b*(a*d*(-a*d + S(2)*b*c)*(m + n + S(-1)) - b**S(2)*(c**S(2)*(m + S(-2)) - d**S(2)*(n + S(1))))*tan(e + f*x)**S(2) - d*(n + S(1))*(-a**S(3)*d + S(3)*a**S(2)*b*c + S(3)*a*b**S(2)*d - b**S(3)*c)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))*(-a*d + b*c)**S(2)/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3641(c, m, n, x, d, f, a, e, b):
        # rubi.append(3641)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-3))*(c + d/tan(e + f*x))**(n + S(1))*Simp(a**S(2)*d*(-a*c*(n + S(1)) + b*d*(m + S(-2))) + b*(-S(2)*a*d + b*c)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b*(a*d*(-a*d + S(2)*b*c)*(m + n + S(-1)) - b**S(2)*(c**S(2)*(m + S(-2)) - d**S(2)*(n + S(1))))/tan(e + f*x)**S(2) - d*(n + S(1))*(-a**S(3)*d + S(3)*a**S(2)*b*c + S(3)*a*b**S(2)*d - b**S(3)*c)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))*(-a*d + b*c)**S(2)/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3642(c, m, n, x, d, f, a, e, b):
        # rubi.append(3642)
        return Dist(S(1)/(d*(m + n + S(-1))), Int((a + b*tan(e + f*x))**(m + S(-3))*(c + d*tan(e + f*x))**n*Simp(a**S(3)*d*(m + n + S(-1)) - b**S(2)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b**S(2)*(-a*d*(S(3)*m + S(2)*n + S(-4)) + b*c*(m + S(-2)))*tan(e + f*x)**S(2) + b*d*(S(3)*a**S(2) - b**S(2))*(m + n + S(-1))*tan(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(-1))), x)
def replacement3643(c, m, n, x, d, f, a, e, b):
        # rubi.append(3643)
        return Dist(S(1)/(d*(m + n + S(-1))), Int((a + b/tan(e + f*x))**(m + S(-3))*(c + d/tan(e + f*x))**n*Simp(a**S(3)*d*(m + n + S(-1)) - b**S(2)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b**S(2)*(-a*d*(S(3)*m + S(2)*n + S(-4)) + b*c*(m + S(-2)))/tan(e + f*x)**S(2) + b*d*(S(3)*a**S(2) - b**S(2))*(m + n + S(-1))/tan(e + f*x), x), x), x) - Simp(b**S(2)*(a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(-1))), x)
def replacement3644(c, m, n, x, d, f, a, e, b):
        # rubi.append(3644)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-2))*Simp(a*c**S(2)*(m + S(1)) + a*d**S(2)*(n + S(-1)) + b*c*d*(m - n + S(2)) - d*(m + n)*(-a*d + b*c)*tan(e + f*x)**S(2) - (m + S(1))*(-S(2)*a*c*d + b*c**S(2) - b*d**S(2))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3645(c, m, n, x, d, f, a, e, b):
        # rubi.append(3645)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-2))*Simp(a*c**S(2)*(m + S(1)) + a*d**S(2)*(n + S(-1)) + b*c*d*(m - n + S(2)) - d*(m + n)*(-a*d + b*c)/tan(e + f*x)**S(2) - (m + S(1))*(-S(2)*a*c*d + b*c**S(2) - b*d**S(2))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-1))*(-a*d + b*c)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3646(c, m, n, x, d, f, a, e, b):
        # rubi.append(3646)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-1))*Simp(a*c*(m + S(1)) - b*d*n - b*d*(m + n + S(1))*tan(e + f*x)**S(2) - (m + S(1))*(-a*d + b*c)*tan(e + f*x), x), x), x) + Simp(b*(a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3647(c, m, n, x, d, f, a, e, b):
        # rubi.append(3647)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-1))*Simp(a*c*(m + S(1)) - b*d*n - b*d*(m + n + S(1))/tan(e + f*x)**S(2) - (m + S(1))*(-a*d + b*c)/tan(e + f*x), x), x), x) - Simp(b*(a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3648(c, m, n, x, d, f, a, e, b):
        # rubi.append(3648)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))*tan(e + f*x)**S(2) - b**S(2)*d*(m + n + S(2)) - b*(m + S(1))*(-a*d + b*c)*tan(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(1))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3649(c, m, n, x, d, f, a, e, b):
        # rubi.append(3649)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2)) - b**S(2)*d*(m + n + S(2))/tan(e + f*x)**S(2) - b*(m + S(1))*(-a*d + b*c)/tan(e + f*x), x), x), x) - Simp(b**S(2)*(a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(1))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3650(c, m, n, x, d, f, a, e, b):
        # rubi.append(3650)
        return Dist(S(1)/(m + n + S(-1)), Int((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(-1))*Simp(a**S(2)*c*(m + n + S(-1)) - b*(a*d*n + b*c*(m + S(-1))) + b*(a*d*(S(2)*m + n + S(-2)) + b*c*n)*tan(e + f*x)**S(2) + (m + n + S(-1))*(a**S(2)*d + S(2)*a*b*c - b**S(2)*d)*tan(e + f*x), x), x), x) + Simp(b*(a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**n/(f*(m + n + S(-1))), x)
def replacement3651(c, m, n, x, d, f, a, e, b):
        # rubi.append(3651)
        return Dist(S(1)/(m + n + S(-1)), Int((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(-1))*Simp(a**S(2)*c*(m + n + S(-1)) - b*(a*d*n + b*c*(m + S(-1))) + b*(a*d*(S(2)*m + n + S(-2)) + b*c*n)/tan(e + f*x)**S(2) + (m + n + S(-1))*(a**S(2)*d + S(2)*a*b*c - b**S(2)*d)/tan(e + f*x), x), x), x) - Simp(b*(a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**n/(f*(m + n + S(-1))), x)
def replacement3652(c, x, d, f, a, e, b):
        # rubi.append(3652)
        return Dist(b**S(2)/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a*tan(e + f*x) + b)/(a + b*tan(e + f*x)), x), x) - Dist(d**S(2)/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c*tan(e + f*x) + d)/(c + d*tan(e + f*x)), x), x) + Simp(x*(a*c - b*d)/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3653(c, x, d, f, a, e, b):
        # rubi.append(3653)
        return Dist(b**S(2)/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a/tan(e + f*x) + b)/(a + b/tan(e + f*x)), x), x) - Dist(d**S(2)/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c/tan(e + f*x) + d)/(c + d/tan(e + f*x)), x), x) + Simp(x*(a*c - b*d)/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3654(c, x, d, f, a, e, b):
        # rubi.append(3654)
        return -Dist(d*(-a*d + b*c)/(c**S(2) + d**S(2)), Int((tan(e + f*x)**S(2) + S(1))/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int(Simp(a*c + b*d + (-a*d + b*c)*tan(e + f*x), x)/sqrt(a + b*tan(e + f*x)), x), x)
def replacement3655(c, x, d, f, a, e, b):
        # rubi.append(3655)
        return -Dist(d*(-a*d + b*c)/(c**S(2) + d**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/(sqrt(a + b/tan(e + f*x))*(c + d/tan(e + f*x))), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int(Simp(a*c + b*d + (-a*d + b*c)/tan(e + f*x), x)/sqrt(a + b/tan(e + f*x)), x), x)
def replacement3656(c, x, d, f, a, e, b):
        # rubi.append(3656)
        return Dist((-a*d + b*c)**S(2)/(c**S(2) + d**S(2)), Int((tan(e + f*x)**S(2) + S(1))/(sqrt(a + b*tan(e + f*x))*(c + d*tan(e + f*x))), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int(Simp(a**S(2)*c + S(2)*a*b*d - b**S(2)*c + (-a**S(2)*d + S(2)*a*b*c + b**S(2)*d)*tan(e + f*x), x)/sqrt(a + b*tan(e + f*x)), x), x)
def replacement3657(c, x, d, f, a, e, b):
        # rubi.append(3657)
        return Dist((-a*d + b*c)**S(2)/(c**S(2) + d**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/(sqrt(a + b/tan(e + f*x))*(c + d/tan(e + f*x))), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int(Simp(a**S(2)*c + S(2)*a*b*d - b**S(2)*c + (-a**S(2)*d + S(2)*a*b*c + b**S(2)*d)/tan(e + f*x), x)/sqrt(a + b/tan(e + f*x)), x), x)
def replacement3658(c, m, x, d, f, a, e, b):
        # rubi.append(3658)
        return Dist(d**S(2)/(c**S(2) + d**S(2)), Int((a + b*tan(e + f*x))**m*(tan(e + f*x)**S(2) + S(1))/(c + d*tan(e + f*x)), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int((a + b*tan(e + f*x))**m*(c - d*tan(e + f*x)), x), x)
def replacement3659(c, m, x, d, f, a, e, b):
        # rubi.append(3659)
        return Dist(d**S(2)/(c**S(2) + d**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))*(a + b/tan(e + f*x))**m/(c + d/tan(e + f*x)), x), x) + Dist(S(1)/(c**S(2) + d**S(2)), Int((a + b/tan(e + f*x))**m*(c - d/tan(e + f*x)), x), x)
def replacement3660(c, m, n, x, d, f, a, e, b):
        # rubi.append(3660)
        return Dist(b/f, Subst(Int((a + x)**m*(c + d*x/b)**n/(b**S(2) + x**S(2)), x), x, b*tan(e + f*x)), x)
def replacement3661(c, m, n, x, d, f, a, e, b):
        # rubi.append(3661)
        return -Dist(b/f, Subst(Int((a + x)**m*(c + d*x/b)**n/(b**S(2) + x**S(2)), x), x, b/tan(e + f*x)), x)
def replacement3662(m, n, x, d, f, a, e, b):
        # rubi.append(3662)
        return Dist(d**m, Int((d/tan(e + f*x))**(-m + n)*(a/tan(e + f*x) + b)**m, x), x)
def replacement3663(m, n, x, d, f, a, e, b):
        # rubi.append(3663)
        return Dist(d**m, Int((d*tan(e + f*x))**(-m + n)*(a*tan(e + f*x) + b)**m, x), x)
def replacement3664(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(3664)
        return Dist(c**IntPart(n)*(c*(d*tan(e + f*x))**p)**FracPart(n)*(d*tan(e + f*x))**(-p*FracPart(n)), Int((d*tan(e + f*x))**(n*p)*(a + b*tan(e + f*x))**m, x), x)
def replacement3665(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(3665)
        return Dist(c**IntPart(n)*(c*(d/tan(e + f*x))**p)**FracPart(n)*(d/tan(e + f*x))**(-p*FracPart(n)), Int((d/tan(e + f*x))**(n*p)*(a + b/tan(e + f*x))**m, x), x)
def replacement3666(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3666)
        return Int((g*tan(e + f*x))**p*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n, x)
def replacement3667(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3667)
        return Int((g/tan(e + f*x))**p*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x)
def replacement3668(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3668)
        return Dist(g**(m + n), Int((g/tan(e + f*x))**(-m - n + p)*(a/tan(e + f*x) + b)**m*(c/tan(e + f*x) + d)**n, x), x)
def replacement3669(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3669)
        return Dist(g**(m + n), Int((g*tan(e + f*x))**(-m - n + p)*(a*tan(e + f*x) + b)**m*(c*tan(e + f*x) + d)**n, x), x)
def replacement3670(c, m, q, n, x, d, f, g, a, p, e, b):
        # rubi.append(3670)
        return Dist((g*tan(e + f*x))**(-p*q)*(g*tan(e + f*x)**q)**p, Int((g*tan(e + f*x))**(p*q)*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n, x), x)
def replacement3671(c, m, q, n, x, d, f, g, a, p, e, b):
        # rubi.append(3671)
        return Dist((g*(S(1)/tan(e + f*x))**q)**p*(g/tan(e + f*x))**(-p*q), Int((g/tan(e + f*x))**(p*q)*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x)
def replacement3672(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3672)
        return Dist(g**n, Int((g*tan(e + f*x))**(-n + p)*(a + b*tan(e + f*x))**m*(c*tan(e + f*x) + d)**n, x), x)
def replacement3673(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3673)
        return Dist(g**n, Int((g/tan(e + f*x))**(-n + p)*(a + b/tan(e + f*x))**m*(c/tan(e + f*x) + d)**n, x), x)
def replacement3674(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(3674)
        return Int((c + d/tan(e + f*x))**n*(a/tan(e + f*x) + b)**m*(S(1)/tan(e + f*x))**(-m - p), x)
def replacement3675(c, m, n, x, d, f, a, p, e, b):
        # rubi.append(3675)
        return Int((c + d*tan(e + f*x))**n*(a*tan(e + f*x) + b)**m*tan(e + f*x)**(-m - p), x)
def replacement3676(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3676)
        return Dist((g*tan(e + f*x))**p*(S(1)/tan(e + f*x))**p, Int((c + d/tan(e + f*x))**n*(a/tan(e + f*x) + b)**m*(S(1)/tan(e + f*x))**(-m - p), x), x)
def replacement3677(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3677)
        return Dist((g/tan(e + f*x))**p*tan(e + f*x)**p, Int((c + d*tan(e + f*x))**n*(a*tan(e + f*x) + b)**m*tan(e + f*x)**(-m - p), x), x)
def replacement3678(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3678)
        return Dist((g*tan(e + f*x))**n*(c + d/tan(e + f*x))**n*(c*tan(e + f*x) + d)**(-n), Int((g*tan(e + f*x))**(-n + p)*(a + b*tan(e + f*x))**m*(c*tan(e + f*x) + d)**n, x), x)
def replacement3679(c, m, n, x, d, f, g, a, p, e, b):
        # rubi.append(3679)
        return Dist((g/tan(e + f*x))**n*(c + d*tan(e + f*x))**n*(c/tan(e + f*x) + d)**(-n), Int((g/tan(e + f*x))**(-n + p)*(a + b/tan(e + f*x))**m*(c/tan(e + f*x) + d)**n, x), x)
def replacement3680(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3680)
        return Dist(a*c/f, Subst(Int((A + B*x)*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1)), x), x, tan(e + f*x)), x)
def replacement3681(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3681)
        return -Dist(a*c/f, Subst(Int((A + B*x)*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1)), x), x, S(1)/tan(e + f*x)), x)
def replacement3682(c, x, d, B, A, f, a, e, b):
        # rubi.append(3682)
        return Dist(S(1)/b, Int(Simp(A*b*c + (A*b*d + B*(-a*d + b*c))*tan(e + f*x), x)/(a + b*tan(e + f*x)), x), x) + Dist(B*d/b, Int(tan(e + f*x), x), x)
def replacement3683(c, x, d, B, A, f, a, e, b):
        # rubi.append(3683)
        return Dist(S(1)/b, Int(Simp(A*b*c + (A*b*d + B*(-a*d + b*c))/tan(e + f*x), x)/(a + b/tan(e + f*x)), x), x) + Dist(B*d/b, Int(S(1)/tan(e + f*x), x), x)
def replacement3684(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3684)
        return Dist(S(1)/(S(2)*a*b), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(A*a*d + A*b*c + B*a*c + S(2)*B*a*d*tan(e + f*x) + B*b*d, x), x), x) - Simp((a + b*tan(e + f*x))**m*(A*b - B*a)*(a*c + b*d)/(S(2)*a**S(2)*f*m), x)
def replacement3685(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3685)
        return Dist(S(1)/(S(2)*a*b), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(A*a*d + A*b*c + B*a*c + S(2)*B*a*d/tan(e + f*x) + B*b*d, x), x), x) + Simp((a + b/tan(e + f*x))**m*(A*b - B*a)*(a*c + b*d)/(S(2)*a**S(2)*f*m), x)
def replacement3686(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3686)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(A*a*c + A*b*d - B*a*d + B*b*c - (-A*a*d + A*b*c - B*a*c - B*b*d)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(A*b - B*a)*(-a*d + b*c)/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3687(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3687)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(A*a*c + A*b*d - B*a*d + B*b*c - (-A*a*d + A*b*c - B*a*c - B*b*d)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(A*b - B*a)*(-a*d + b*c)/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3688(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3688)
        return Int((a + b*tan(e + f*x))**m*Simp(A*c - B*d + (A*d + B*c)*tan(e + f*x), x), x) + Simp(B*d*(a + b*tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3689(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3689)
        return Int((a + b/tan(e + f*x))**m*Simp(A*c - B*d + (A*d + B*c)/tan(e + f*x), x), x) - Simp(B*d*(a + b/tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3690(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3690)
        return -Dist(a/(d*(n + S(1))*(a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))*Simp(A*b*d*(m - n + S(-2)) - B*(a*d*(n + S(1)) + b*c*(m + S(-1))) + (A*a*d*(m + n) - B*(a*c*(m + S(-1)) + b*d*(n + S(1))))*tan(e + f*x), x), x), x) - Simp(a**S(2)*(a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))*(-A*d + B*c)/(d*f*(n + S(1))*(a*d + b*c)), x)
def replacement3691(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3691)
        return -Dist(a/(d*(n + S(1))*(a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))*Simp(A*b*d*(m - n + S(-2)) - B*(a*d*(n + S(1)) + b*c*(m + S(-1))) + (A*a*d*(m + n) - B*(a*c*(m + S(-1)) + b*d*(n + S(1))))/tan(e + f*x), x), x), x) + Simp(a**S(2)*(a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))*(-A*d + B*c)/(d*f*(n + S(1))*(a*d + b*c)), x)
def replacement3692(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3692)
        return Dist(S(1)/(d*(m + n)), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**n*Simp(A*a*d*(m + n) + B*(a*c*(m + S(-1)) - b*d*(n + S(1))) - (B*(m + S(-1))*(-a*d + b*c) - d*(m + n)*(A*b + B*a))*tan(e + f*x), x), x), x) + Simp(B*b*(a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n)), x)
def replacement3693(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3693)
        return Dist(S(1)/(d*(m + n)), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**n*Simp(A*a*d*(m + n) + B*(a*c*(m + S(-1)) - b*d*(n + S(1))) - (B*(m + S(-1))*(-a*d + b*c) - d*(m + n)*(A*b + B*a))/tan(e + f*x), x), x), x) - Simp(B*b*(a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n)), x)
def replacement3694(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3694)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-1))*Simp(A*(a*c*m + b*d*n) - B*(a*d*n + b*c*m) - d*(-A*a*(m + n) + B*b*(m - n))*tan(e + f*x), x), x), x) - Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(A*b - B*a)/(S(2)*a*f*m), x)
def replacement3695(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3695)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-1))*Simp(A*(a*c*m + b*d*n) - B*(a*d*n + b*c*m) - d*(-A*a*(m + n) + B*b*(m - n))/tan(e + f*x), x), x), x) + Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n*(A*b - B*a)/(S(2)*a*f*m), x)
def replacement3696(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3696)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(A*(-a*d*(S(2)*m + n + S(1)) + b*c*m) + B*(a*c*m - b*d*(n + S(1))) + d*(A*b - B*a)*(m + n + S(1))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*a + B*b)/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3697(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3697)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(A*(-a*d*(S(2)*m + n + S(1)) + b*c*m) + B*(a*c*m - b*d*(n + S(1))) + d*(A*b - B*a)*(m + n + S(1))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*a + B*b)/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3698(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3698)
        return Dist(S(1)/(a*(m + n)), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(-1))*Simp(A*a*c*(m + n) - B*(a*d*n + b*c*m) + (A*a*d*(m + n) - B*(-a*c*n + b*d*m))*tan(e + f*x), x), x), x) + Simp(B*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n/(f*(m + n)), x)
def replacement3699(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3699)
        return Dist(S(1)/(a*(m + n)), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(-1))*Simp(A*a*c*(m + n) - B*(a*d*n + b*c*m) + (A*a*d*(m + n) - B*(-a*c*n + b*d*m))/tan(e + f*x), x), x), x) - Simp(B*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n/(f*(m + n)), x)
def replacement3700(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3700)
        return -Dist(S(1)/(a*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*Simp(A*(-a*c*(n + S(1)) + b*d*m) - B*(a*d*(n + S(1)) + b*c*m) - a*(-A*d + B*c)*(m + n + S(1))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*d - B*c)/(f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3701(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3701)
        return -Dist(S(1)/(a*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*Simp(A*(-a*c*(n + S(1)) + b*d*m) - B*(a*d*(n + S(1)) + b*c*m) - a*(-A*d + B*c)*(m + n + S(1))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*d - B*c)/(f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3702(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3702)
        return Dist(B*b/f, Subst(Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x), x, tan(e + f*x)), x)
def replacement3703(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3703)
        return -Dist(B*b/f, Subst(Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x), x, S(1)/tan(e + f*x)), x)
def replacement3704(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3704)
        return Dist((A*b + B*a)/(a*d + b*c), Int((a + b*tan(e + f*x))**m, x), x) - Dist((-A*d + B*c)/(a*d + b*c), Int((a - b*tan(e + f*x))*(a + b*tan(e + f*x))**m/(c + d*tan(e + f*x)), x), x)
def replacement3705(c, m, x, d, B, A, f, a, e, b):
        # rubi.append(3705)
        return Dist((A*b + B*a)/(a*d + b*c), Int((a + b/tan(e + f*x))**m, x), x) - Dist((-A*d + B*c)/(a*d + b*c), Int((a - b/tan(e + f*x))*(a + b/tan(e + f*x))**m/(c + d/tan(e + f*x)), x), x)
def replacement3706(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3706)
        return -Dist(B/b, Int((a - b*tan(e + f*x))*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n, x), x) + Dist((A*b + B*a)/b, Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n, x), x)
def replacement3707(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3707)
        return -Dist(B/b, Int((a - b/tan(e + f*x))*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x) + Dist((A*b + B*a)/b, Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x)
def replacement3708(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3708)
        return Dist(A**S(2)/f, Subst(Int((a + b*x)**m*(c + d*x)**n/(A - B*x), x), x, tan(e + f*x)), x)
def replacement3709(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3709)
        return -Dist(A**S(2)/f, Subst(Int((a + b*x)**m*(c + d*x)**n/(A - B*x), x), x, S(1)/tan(e + f*x)), x)
def replacement3710(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3710)
        return Dist(A/S(2) - I*B/S(2), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(I*tan(e + f*x) + S(1)), x), x) + Dist(A/S(2) + I*B/S(2), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(-I*tan(e + f*x) + S(1)), x), x)
def replacement3711(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3711)
        return Dist(A/S(2) - I*B/S(2), Int((S(1) + I/tan(e + f*x))*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x) + Dist(A/S(2) + I*B/S(2), Int((S(1) - I/tan(e + f*x))*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x)
def replacement3712(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3712)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**(n + S(1))*Simp(A*a*d*(-a*c*(n + S(1)) + b*d*(m + S(-1))) - b*(-B*b*(c**S(2)*(m + S(-1)) - d**S(2)*(n + S(1))) + d*(m + n)*(-A*a*d + A*b*c + B*a*c))*tan(e + f*x)**S(2) - d*(n + S(1))*((A*a - B*b)*(-a*d + b*c) + (A*b + B*a)*(a*c + b*d))*tan(e + f*x) + (B*b*c - d*(A*b + B*a))*(a*d*(n + S(1)) + b*c*(m + S(-1))), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))*(-A*d + B*c)*(-a*d + b*c)/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3713(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3713)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**(n + S(1))*Simp(A*a*d*(-a*c*(n + S(1)) + b*d*(m + S(-1))) - b*(-B*b*(c**S(2)*(m + S(-1)) - d**S(2)*(n + S(1))) + d*(m + n)*(-A*a*d + A*b*c + B*a*c))/tan(e + f*x)**S(2) - d*(n + S(1))*((A*a - B*b)*(-a*d + b*c) + (A*b + B*a)*(a*c + b*d))/tan(e + f*x) + (B*b*c - d*(A*b + B*a))*(a*d*(n + S(1)) + b*c*(m + S(-1))), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))*(-A*d + B*c)*(-a*d + b*c)/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3714(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3714)
        return Dist(S(1)/(d*(m + n)), Int((a + b*tan(e + f*x))**(m + S(-2))*(c + d*tan(e + f*x))**n*Simp(A*a**S(2)*d*(m + n) - B*b*(a*d*(n + S(1)) + b*c*(m + S(-1))) + d*(m + n)*(S(2)*A*a*b + B*(a**S(2) - b**S(2)))*tan(e + f*x) - (B*b*(m + S(-1))*(-a*d + b*c) - b*d*(m + n)*(A*b + B*a))*tan(e + f*x)**S(2), x), x), x) + Simp(B*b*(a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n)), x)
def replacement3715(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3715)
        return Dist(S(1)/(d*(m + n)), Int((a + b/tan(e + f*x))**(m + S(-2))*(c + d/tan(e + f*x))**n*Simp(A*a**S(2)*d*(m + n) - B*b*(a*d*(n + S(1)) + b*c*(m + S(-1))) + d*(m + n)*(S(2)*A*a*b + B*(a**S(2) - b**S(2)))/tan(e + f*x) - (B*b*(m + S(-1))*(-a*d + b*c) - b*d*(m + n)*(A*b + B*a))/tan(e + f*x)**S(2), x), x), x) - Simp(B*b*(a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n)), x)
def replacement3716(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3716)
        return Dist(S(1)/(b*(a**S(2) + b**S(2))*(m + S(1))), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(-1))*Simp(A*b*(a*c*(m + S(1)) - b*d*n) + B*b*(a*d*n + b*c*(m + S(1))) - b*d*(A*b - B*a)*(m + n + S(1))*tan(e + f*x)**S(2) - b*(m + S(1))*(A*(-a*d + b*c) - B*(a*c + b*d))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*(A*b - B*a)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3717(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3717)
        return Dist(S(1)/(b*(a**S(2) + b**S(2))*(m + S(1))), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(-1))*Simp(A*b*(a*c*(m + S(1)) - b*d*n) + B*b*(a*d*n + b*c*(m + S(1))) - b*d*(A*b - B*a)*(m + n + S(1))/tan(e + f*x)**S(2) - b*(m + S(1))*(A*(-a*d + b*c) - B*(a*c + b*d))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*(A*b - B*a)/(f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3718(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3718)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) + B*b*(a*d*(n + S(1)) + b*c*(m + S(1))) - b*d*(A*b - B*a)*(m + n + S(2))*tan(e + f*x)**S(2) - (m + S(1))*(A*b - B*a)*(-a*d + b*c)*tan(e + f*x), x), x), x) + Simp(b*(a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(1))*(A*b - B*a)/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3719(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3719)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) + B*b*(a*d*(n + S(1)) + b*c*(m + S(1))) - b*d*(A*b - B*a)*(m + n + S(2))/tan(e + f*x)**S(2) - (m + S(1))*(A*b - B*a)*(-a*d + b*c)/tan(e + f*x), x), x), x) - Simp(b*(a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(1))*(A*b - B*a)/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3720(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3720)
        return Dist(S(1)/(m + n), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(-1))*Simp(A*a*c*(m + n) - B*(a*d*n + b*c*m) + (m + n)*(A*a*d + A*b*c + B*a*c - B*b*d)*tan(e + f*x) + (A*b*d*(m + n) + B*(a*d*m + b*c*n))*tan(e + f*x)**S(2), x), x), x) + Simp(B*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n/(f*(m + n)), x)
def replacement3721(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3721)
        return Dist(S(1)/(m + n), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(-1))*Simp(A*a*c*(m + n) - B*(a*d*n + b*c*m) + (m + n)*(A*a*d + A*b*c + B*a*c - B*b*d)/tan(e + f*x) + (A*b*d*(m + n) + B*(a*d*m + b*c*n))/tan(e + f*x)**S(2), x), x), x) - Simp(B*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n/(f*(m + n)), x)
def replacement3722(c, x, d, B, A, f, a, e, b):
        # rubi.append(3722)
        return Dist(b*(A*b - B*a)/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a*tan(e + f*x) + b)/(a + b*tan(e + f*x)), x), x) + Dist(d*(-A*d + B*c)/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c*tan(e + f*x) + d)/(c + d*tan(e + f*x)), x), x) + Simp(x*(A*(a*c - b*d) + B*(a*d + b*c))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3723(c, x, d, B, A, f, a, e, b):
        # rubi.append(3723)
        return Dist(b*(A*b - B*a)/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a/tan(e + f*x) + b)/(a + b/tan(e + f*x)), x), x) + Dist(d*(-A*d + B*c)/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c/tan(e + f*x) + d)/(c + d/tan(e + f*x)), x), x) + Simp(x*(A*(a*c - b*d) + B*(a*d + b*c))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3724(c, x, d, B, A, f, a, e, b):
        # rubi.append(3724)
        return -Dist((-A*b + B*a)*(-a*d + b*c)/(a**S(2) + b**S(2)), Int((tan(e + f*x)**S(2) + S(1))/((a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int(Simp(A*(a*c + b*d) + B*(-a*d + b*c) - (A*(-a*d + b*c) - B*(a*c + b*d))*tan(e + f*x), x)/sqrt(c + d*tan(e + f*x)), x), x)
def replacement3725(c, x, d, B, A, f, a, e, b):
        # rubi.append(3725)
        return -Dist((-A*b + B*a)*(-a*d + b*c)/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/((a + b/tan(e + f*x))*sqrt(c + d/tan(e + f*x))), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int(Simp(A*(a*c + b*d) + B*(-a*d + b*c) - (A*(-a*d + b*c) - B*(a*c + b*d))/tan(e + f*x), x)/sqrt(c + d/tan(e + f*x)), x), x)
def replacement3726(c, n, x, d, B, A, f, a, e, b):
        # rubi.append(3726)
        return Dist(b*(A*b - B*a)/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*(tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*Simp(A*a + B*b - (A*b - B*a)*tan(e + f*x), x), x), x)
def replacement3727(c, n, x, d, B, A, f, a, e, b):
        # rubi.append(3727)
        return Dist(b*(A*b - B*a)/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))*(c + d/tan(e + f*x))**n/(a + b/tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d/tan(e + f*x))**n*Simp(A*a + B*b - (A*b - B*a)/tan(e + f*x), x), x), x)
def replacement3728(c, x, d, B, A, f, a, e, b):
        # rubi.append(3728)
        return Dist(B*b, Int((tan(e + f*x)**S(2) + S(1))/(sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))), x), x) + Int(Simp(A*a - B*b + (A*b + B*a)*tan(e + f*x), x)/(sqrt(a + b*tan(e + f*x))*sqrt(c + d*tan(e + f*x))), x)
def replacement3729(c, x, d, B, A, f, a, e, b):
        # rubi.append(3729)
        return Dist(B*b, Int((S(1) + tan(e + f*x)**(S(-2)))/(sqrt(a + b/tan(e + f*x))*sqrt(c + d/tan(e + f*x))), x), x) + Int(Simp(A*a - B*b + (A*b + B*a)/tan(e + f*x), x)/(sqrt(a + b/tan(e + f*x))*sqrt(c + d/tan(e + f*x))), x)
def replacement3730(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3730)
        return Dist(A**S(2)/f, Subst(Int((a + b*x)**m*(c + d*x)**n/(A - B*x), x), x, tan(e + f*x)), x)
def replacement3731(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3731)
        return -Dist(A**S(2)/f, Subst(Int((a + b*x)**m*(c + d*x)**n/(A - B*x), x), x, S(1)/tan(e + f*x)), x)
def replacement3732(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3732)
        return Dist(A/S(2) - I*B/S(2), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(I*tan(e + f*x) + S(1)), x), x) + Dist(A/S(2) + I*B/S(2), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*(-I*tan(e + f*x) + S(1)), x), x)
def replacement3733(c, m, n, x, d, B, A, f, a, e, b):
        # rubi.append(3733)
        return Dist(A/S(2) - I*B/S(2), Int((S(1) + I/tan(e + f*x))*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x) + Dist(A/S(2) + I*B/S(2), Int((S(1) - I/tan(e + f*x))*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n, x), x)
def replacement3734(m, x, A, f, a, C, e, b):
        # rubi.append(3734)
        return Dist(A/(b*f), Subst(Int((a + x)**m, x), x, b*tan(e + f*x)), x)
def replacement3735(m, x, A, f, a, C, e, b):
        # rubi.append(3735)
        return -Dist(A/(b*f), Subst(Int((a + x)**m, x), x, b/tan(e + f*x)), x)
def replacement3736(m, x, B, A, f, a, C, e, b):
        # rubi.append(3736)
        return Dist(b**(S(-2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b*tan(e + f*x), x), x), x)
def replacement3737(m, x, B, A, f, a, C, e, b):
        # rubi.append(3737)
        return Dist(b**(S(-2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b/tan(e + f*x), x), x), x)
def replacement3738(m, x, A, f, a, C, e, b):
        # rubi.append(3738)
        return -Dist(C/b**S(2), Int((a - b*tan(e + f*x))*(a + b*tan(e + f*x))**(m + S(1)), x), x)
def replacement3739(m, x, A, f, a, C, e, b):
        # rubi.append(3739)
        return -Dist(C/b**S(2), Int((a - b/tan(e + f*x))*(a + b/tan(e + f*x))**(m + S(1)), x), x)
def replacement3740(m, x, B, A, f, a, C, e, b):
        # rubi.append(3740)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(A*a*(S(2)*m + S(1)) + B*b - C*a - (C*b*(m + S(-1)) + (m + S(1))*(A*b - B*a))*tan(e + f*x), x), x), x) - Simp((a + b*tan(e + f*x))**m*(A*a + B*b - C*a)*tan(e + f*x)/(S(2)*a*f*m), x)
def replacement3741(m, x, B, A, f, a, C, e, b):
        # rubi.append(3741)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(A*a*(S(2)*m + S(1)) + B*b - C*a - (C*b*(m + S(-1)) + (m + S(1))*(A*b - B*a))/tan(e + f*x), x), x), x) + Simp((a + b/tan(e + f*x))**m*(A*a + B*b - C*a)/(S(2)*a*f*m*tan(e + f*x)), x)
def replacement3742(m, x, A, f, a, C, e, b):
        # rubi.append(3742)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(A*a*(S(2)*m + S(1)) - C*a - (A*b*(m + S(1)) + C*b*(m + S(-1)))*tan(e + f*x), x), x), x) - Simp((a + b*tan(e + f*x))**m*(A*a - C*a)*tan(e + f*x)/(S(2)*a*f*m), x)
def replacement3743(m, x, A, f, a, C, e, b):
        # rubi.append(3743)
        return Dist(S(1)/(S(2)*a**S(2)*m), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(A*a*(S(2)*m + S(1)) - C*a - (A*b*(m + S(1)) + C*b*(m + S(-1)))/tan(e + f*x), x), x), x) + Simp((a + b/tan(e + f*x))**m*(A*a - C*a)/(S(2)*a*f*m*tan(e + f*x)), x)
def replacement3744(x, B, A, f, a, C, e, b):
        # rubi.append(3744)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) + Simp(x*(A*a + B*b - C*a)/(a**S(2) + b**S(2)), x)
def replacement3745(x, B, A, f, a, C, e, b):
        # rubi.append(3745)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/(a + b/tan(e + f*x)), x), x) + Simp(x*(A*a + B*b - C*a)/(a**S(2) + b**S(2)), x)
def replacement3746(x, B, A, f, C, e):
        # rubi.append(3746)
        return Dist(A, Int(S(1)/tan(e + f*x), x), x) + Dist(C, Int(tan(e + f*x), x), x) + Simp(B*x, x)
def replacement3747(x, B, A, f, C, e):
        # rubi.append(3747)
        return Dist(A, Int(tan(e + f*x), x), x) + Dist(C, Int(S(1)/tan(e + f*x), x), x) + Simp(B*x, x)
def replacement3748(x, A, f, C, e):
        # rubi.append(3748)
        return Dist(A, Int(S(1)/tan(e + f*x), x), x) + Dist(C, Int(tan(e + f*x), x), x)
def replacement3749(x, A, f, C, e):
        # rubi.append(3749)
        return Dist(A, Int(tan(e + f*x), x), x) + Dist(C, Int(S(1)/tan(e + f*x), x), x)
def replacement3750(x, B, A, f, a, C, e, b):
        # rubi.append(3750)
        return -Dist((A*b - B*a - C*b)/(a**S(2) + b**S(2)), Int(tan(e + f*x), x), x) + Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) + Simp(x*(A*a + B*b - C*a)/(a**S(2) + b**S(2)), x)
def replacement3751(x, B, A, f, a, C, e, b):
        # rubi.append(3751)
        return -Dist((A*b - B*a - C*b)/(a**S(2) + b**S(2)), Int(S(1)/tan(e + f*x), x), x) + Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/(a + b/tan(e + f*x)), x), x) + Simp(x*(A*a + B*b - C*a)/(a**S(2) + b**S(2)), x)
def replacement3752(x, A, f, a, C, e, b):
        # rubi.append(3752)
        return Dist((A*b**S(2) + C*a**S(2))/(a**S(2) + b**S(2)), Int((tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) - Dist(b*(A - C)/(a**S(2) + b**S(2)), Int(tan(e + f*x), x), x) + Simp(a*x*(A - C)/(a**S(2) + b**S(2)), x)
def replacement3753(x, A, f, a, C, e, b):
        # rubi.append(3753)
        return Dist((A*b**S(2) + C*a**S(2))/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))/(a + b/tan(e + f*x)), x), x) - Dist(b*(A - C)/(a**S(2) + b**S(2)), Int(S(1)/tan(e + f*x), x), x) + Simp(a*x*(A - C)/(a**S(2) + b**S(2)), x)
def replacement3754(m, x, B, A, f, a, C, e, b):
        # rubi.append(3754)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(B*b + a*(A - C) - (A*b - B*a - C*b)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3755(m, x, B, A, f, a, C, e, b):
        # rubi.append(3755)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(B*b + a*(A - C) - (A*b - B*a - C*b)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3756(m, x, A, f, a, C, e, b):
        # rubi.append(3756)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b*tan(e + f*x))**(m + S(1))*Simp(a*(A - C) - (A*b - C*b)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3757(m, x, A, f, a, C, e, b):
        # rubi.append(3757)
        return Dist(S(1)/(a**S(2) + b**S(2)), Int((a + b/tan(e + f*x))**(m + S(1))*Simp(a*(A - C) - (A*b - C*b)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))/(b*f*(a**S(2) + b**S(2))*(m + S(1))), x)
def replacement3758(m, x, B, A, f, a, C, e, b):
        # rubi.append(3758)
        return Int((a + b*tan(e + f*x))**m*Simp(A + B*tan(e + f*x) - C, x), x) + Simp(C*(a + b*tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3759(m, x, B, A, f, a, C, e, b):
        # rubi.append(3759)
        return Int((a + b/tan(e + f*x))**m*Simp(A + B/tan(e + f*x) - C, x), x) - Simp(C*(a + b/tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3760(m, x, A, f, a, C, e, b):
        # rubi.append(3760)
        return Dist(A - C, Int((a + b*tan(e + f*x))**m, x), x) + Simp(C*(a + b*tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3761(m, x, A, f, a, C, e, b):
        # rubi.append(3761)
        return Dist(A - C, Int((a + b/tan(e + f*x))**m, x), x) - Simp(C*(a + b/tan(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)
def replacement3762(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3762)
        return Dist(A/f, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, tan(e + f*x)), x)
def replacement3763(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3763)
        return -Dist(A/f, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, S(1)/tan(e + f*x)), x)
def replacement3764(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3764)
        return Dist(S(1)/(d*(c**S(2) + d**S(2))), Int((c + d*tan(e + f*x))**(n + S(1))*Simp(C*b*(c**S(2) + d**S(2))*tan(e + f*x)**S(2) + a*d*(A*c + B*d - C*c) + b*(A*d**S(2) - B*c*d + C*c**S(2)) + d*(-A*a*d + A*b*c + B*a*c + B*b*d + C*a*d - C*b*c)*tan(e + f*x), x), x), x) - Simp((c + d*tan(e + f*x))**(n + S(1))*(-a*d + b*c)*(A*d**S(2) - B*c*d + C*c**S(2))/(d**S(2)*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3765(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3765)
        return Dist(S(1)/(d*(c**S(2) + d**S(2))), Int((c + d/tan(e + f*x))**(n + S(1))*Simp(C*b*(c**S(2) + d**S(2))/tan(e + f*x)**S(2) + a*d*(A*c + B*d - C*c) + b*(A*d**S(2) - B*c*d + C*c**S(2)) + d*(-A*a*d + A*b*c + B*a*c + B*b*d + C*a*d - C*b*c)/tan(e + f*x), x), x), x) + Simp((c + d/tan(e + f*x))**(n + S(1))*(-a*d + b*c)*(A*d**S(2) - B*c*d + C*c**S(2))/(d**S(2)*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3766(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3766)
        return Dist(S(1)/(d*(c**S(2) + d**S(2))), Int((c + d*tan(e + f*x))**(n + S(1))*Simp(C*b*(c**S(2) + d**S(2))*tan(e + f*x)**S(2) + a*d*(A*c - C*c) + b*(A*d**S(2) + C*c**S(2)) + d*(-A*a*d + A*b*c + C*a*d - C*b*c)*tan(e + f*x), x), x), x) - Simp((c + d*tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*(-a*d + b*c)/(d**S(2)*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3767(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3767)
        return Dist(S(1)/(d*(c**S(2) + d**S(2))), Int((c + d/tan(e + f*x))**(n + S(1))*Simp(C*b*(c**S(2) + d**S(2))/tan(e + f*x)**S(2) + a*d*(A*c - C*c) + b*(A*d**S(2) + C*c**S(2)) + d*(-A*a*d + A*b*c + C*a*d - C*b*c)/tan(e + f*x), x), x), x) + Simp((c + d/tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*(-a*d + b*c)/(d**S(2)*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3768(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3768)
        return -Dist(S(1)/(d*(n + S(2))), Int((c + d*tan(e + f*x))**n*Simp(-A*a*d*(n + S(2)) + C*b*c - d*(n + S(2))*(A*b + B*a - C*b)*tan(e + f*x) - (C*a*d*(n + S(2)) - b*(-B*d*(n + S(2)) + C*c))*tan(e + f*x)**S(2), x), x), x) + Simp(C*b*(c + d*tan(e + f*x))**(n + S(1))*tan(e + f*x)/(d*f*(n + S(2))), x)
def replacement3769(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3769)
        return -Dist(S(1)/(d*(n + S(2))), Int((c + d/tan(e + f*x))**n*Simp(-A*a*d*(n + S(2)) + C*b*c - d*(n + S(2))*(A*b + B*a - C*b)/tan(e + f*x) - (C*a*d*(n + S(2)) - b*(-B*d*(n + S(2)) + C*c))/tan(e + f*x)**S(2), x), x), x) - Simp(C*b*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(n + S(2))*tan(e + f*x)), x)
def replacement3770(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3770)
        return -Dist(S(1)/(d*(n + S(2))), Int((c + d*tan(e + f*x))**n*Simp(-A*a*d*(n + S(2)) + C*b*c - d*(n + S(2))*(A*b - C*b)*tan(e + f*x) - (C*a*d*(n + S(2)) - C*b*c)*tan(e + f*x)**S(2), x), x), x) + Simp(C*b*(c + d*tan(e + f*x))**(n + S(1))*tan(e + f*x)/(d*f*(n + S(2))), x)
def replacement3771(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3771)
        return -Dist(S(1)/(d*(n + S(2))), Int((c + d/tan(e + f*x))**n*Simp(-A*a*d*(n + S(2)) + C*b*c - d*(n + S(2))*(A*b - C*b)/tan(e + f*x) - (C*a*d*(n + S(2)) - C*b*c)/tan(e + f*x)**S(2), x), x), x) - Simp(C*b*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(n + S(2))*tan(e + f*x)), x)
def replacement3772(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3772)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(a*(-A*d*(S(2)*m + n + S(1)) + B*c*m + C*d*(n + S(1))) + b*(-B*d*(n + S(1)) + c*m*(A + C)) + (A*b*d*(m + n + S(1)) + C*b*d*(m - n + S(-1)) + a*(-B*d*(m + n + S(1)) + S(2)*C*c*m))*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*a + B*b - C*a)/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3773(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3773)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(a*(-A*d*(S(2)*m + n + S(1)) + B*c*m + C*d*(n + S(1))) + b*(-B*d*(n + S(1)) + c*m*(A + C)) + (A*b*d*(m + n + S(1)) + C*b*d*(m - n + S(-1)) + a*(-B*d*(m + n + S(1)) + S(2)*C*c*m))/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*a + B*b - C*a)/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3774(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3774)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(a*(-A*d*(S(2)*m + n + S(1)) + C*d*(n + S(1))) + b*c*m*(A + C) + (A*b*d*(m + n + S(1)) + S(2)*C*a*c*m + C*b*d*(m - n + S(-1)))*tan(e + f*x), x), x), x) + Simp(a*(A - C)*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3775(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3775)
        return Dist(S(1)/(S(2)*a*m*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(a*(-A*d*(S(2)*m + n + S(1)) + C*d*(n + S(1))) + b*c*m*(A + C) + (A*b*d*(m + n + S(1)) + S(2)*C*a*c*m + C*b*d*(m - n + S(-1)))/tan(e + f*x), x), x), x) - Simp(a*(A - C)*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(S(2)*f*m*(-a*d + b*c)), x)
def replacement3776(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3776)
        return -Dist(S(1)/(a*d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*Simp(-a*d*(n + S(1))*(A*c + B*d - C*c) - a*(-C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(1)))*tan(e + f*x) + b*m*(A*d**S(2) - B*c*d + C*c**S(2)), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3777(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3777)
        return -Dist(S(1)/(a*d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*Simp(-a*d*(n + S(1))*(A*c + B*d - C*c) - a*(-C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(1)))/tan(e + f*x) + b*m*(A*d**S(2) - B*c*d + C*c**S(2)), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3778(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3778)
        return -Dist(S(1)/(a*d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*Simp(-a*d*(n + S(1))*(A*c - C*c) - a*(-A*d**S(2)*(m + n + S(1)) - C*(c**S(2)*m - d**S(2)*(n + S(1))))*tan(e + f*x) + b*m*(A*d**S(2) + C*c**S(2)), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3779(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3779)
        return -Dist(S(1)/(a*d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*Simp(-a*d*(n + S(1))*(A*c - C*c) - a*(-A*d**S(2)*(m + n + S(1)) - C*(c**S(2)*m - d**S(2)*(n + S(1))))/tan(e + f*x) + b*m*(A*d**S(2) + C*c**S(2)), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3780(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3780)
        return Dist(S(1)/(b*d*(m + n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*Simp(A*b*d*(m + n + S(1)) + C*(a*c*m - b*d*(n + S(1))) - (-B*b*d*(m + n + S(1)) + C*m*(-a*d + b*c))*tan(e + f*x), x), x), x) + Simp(C*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3781(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3781)
        return Dist(S(1)/(b*d*(m + n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n*Simp(A*b*d*(m + n + S(1)) + C*(a*c*m - b*d*(n + S(1))) - (-B*b*d*(m + n + S(1)) + C*m*(-a*d + b*c))/tan(e + f*x), x), x), x) - Simp(C*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3782(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3782)
        return Dist(S(1)/(b*d*(m + n + S(1))), Int((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**n*Simp(A*b*d*(m + n + S(1)) - C*m*(-a*d + b*c)*tan(e + f*x) + C*(a*c*m - b*d*(n + S(1))), x), x), x) + Simp(C*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3783(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3783)
        return Dist(S(1)/(b*d*(m + n + S(1))), Int((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**n*Simp(A*b*d*(m + n + S(1)) - C*m*(-a*d + b*c)/tan(e + f*x) + C*(a*c*m - b*d*(n + S(1))), x), x), x) - Simp(C*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3784(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3784)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))*Simp(A*d*(-a*c*(n + S(1)) + b*d*m) - b*(-C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(1)))*tan(e + f*x)**S(2) - d*(n + S(1))*(B*(a*c + b*d) + (A - C)*(-a*d + b*c))*tan(e + f*x) + (-B*d + C*c)*(a*d*(n + S(1)) + b*c*m), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*d**S(2) + c*(-B*d + C*c))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3785(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3785)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))*Simp(A*d*(-a*c*(n + S(1)) + b*d*m) - b*(-C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(1)))/tan(e + f*x)**S(2) - d*(n + S(1))*(B*(a*c + b*d) + (A - C)*(-a*d + b*c))/tan(e + f*x) + (-B*d + C*c)*(a*d*(n + S(1)) + b*c*m), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*d**S(2) + c*(-B*d + C*c))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3786(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3786)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**(n + S(1))*Simp(A*d*(-a*c*(n + S(1)) + b*d*m) + C*c*(a*d*(n + S(1)) + b*c*m) + b*(A*d**S(2)*(m + n + S(1)) + C*(c**S(2)*m - d**S(2)*(n + S(1))))*tan(e + f*x)**S(2) - d*(A - C)*(n + S(1))*(-a*d + b*c)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3787(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3787)
        return -Dist(S(1)/(d*(c**S(2) + d**S(2))*(n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**(n + S(1))*Simp(A*d*(-a*c*(n + S(1)) + b*d*m) + C*c*(a*d*(n + S(1)) + b*c*m) + b*(A*d**S(2)*(m + n + S(1)) + C*(c**S(2)*m - d**S(2)*(n + S(1))))/tan(e + f*x)**S(2) - d*(A - C)*(n + S(1))*(-a*d + b*c)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))/(d*f*(c**S(2) + d**S(2))*(n + S(1))), x)
def replacement3788(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3788)
        return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) - C*(a*d*(n + S(1)) + b*c*m) + d*(m + n + S(1))*(A*b + B*a - C*b)*tan(e + f*x) - (-B*b*d*(m + n + S(1)) + C*m*(-a*d + b*c))*tan(e + f*x)**S(2), x), x), x) + Simp(C*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3789(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3789)
        return Dist(S(1)/(d*(m + n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) - C*(a*d*(n + S(1)) + b*c*m) + d*(m + n + S(1))*(A*b + B*a - C*b)/tan(e + f*x) - (-B*b*d*(m + n + S(1)) + C*m*(-a*d + b*c))/tan(e + f*x)**S(2), x), x), x) - Simp(C*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3790(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3790)
        return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*tan(e + f*x))**(m + S(-1))*(c + d*tan(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) - C*m*(-a*d + b*c)*tan(e + f*x)**S(2) - C*(a*d*(n + S(1)) + b*c*m) + d*(A*b - C*b)*(m + n + S(1))*tan(e + f*x), x), x), x) + Simp(C*(a + b*tan(e + f*x))**m*(c + d*tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3791(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3791)
        return Dist(S(1)/(d*(m + n + S(1))), Int((a + b/tan(e + f*x))**(m + S(-1))*(c + d/tan(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) - C*m*(-a*d + b*c)/tan(e + f*x)**S(2) - C*(a*d*(n + S(1)) + b*c*m) + d*(A*b - C*b)*(m + n + S(1))/tan(e + f*x), x), x), x) - Simp(C*(a + b/tan(e + f*x))**m*(c + d/tan(e + f*x))**(n + S(1))/(d*f*(m + n + S(1))), x)
def replacement3792(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3792)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) - d*(A*b**S(2) - a*(B*b - C*a))*(m + n + S(2))*tan(e + f*x)**S(2) - (m + S(1))*(-a*d + b*c)*(A*b - B*a - C*b)*tan(e + f*x) + (B*b - C*a)*(a*d*(n + S(1)) + b*c*(m + S(1))), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(1))*(A*b**S(2) - a*(B*b - C*a))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3793(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3793)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) - d*(A*b**S(2) - a*(B*b - C*a))*(m + n + S(2))/tan(e + f*x)**S(2) - (m + S(1))*(-a*d + b*c)*(A*b - B*a - C*b)/tan(e + f*x) + (B*b - C*a)*(a*d*(n + S(1)) + b*c*(m + S(1))), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(1))*(A*b**S(2) - a*(B*b - C*a))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3794(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3794)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) - C*a*(a*d*(n + S(1)) + b*c*(m + S(1))) - d*(A*b**S(2) + C*a**S(2))*(m + n + S(2))*tan(e + f*x)**S(2) - (m + S(1))*(A*b - C*b)*(-a*d + b*c)*tan(e + f*x), x), x), x) + Simp((a + b*tan(e + f*x))**(m + S(1))*(c + d*tan(e + f*x))**(n + S(1))*(A*b**S(2) + C*a**S(2))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3795(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3795)
        return Dist(S(1)/((a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**n*Simp(A*(a*(m + S(1))*(-a*d + b*c) - b**S(2)*d*(m + n + S(2))) - C*a*(a*d*(n + S(1)) + b*c*(m + S(1))) - d*(A*b**S(2) + C*a**S(2))*(m + n + S(2))/tan(e + f*x)**S(2) - (m + S(1))*(A*b - C*b)*(-a*d + b*c)/tan(e + f*x), x), x), x) - Simp((a + b/tan(e + f*x))**(m + S(1))*(c + d/tan(e + f*x))**(n + S(1))*(A*b**S(2) + C*a**S(2))/(f*(a**S(2) + b**S(2))*(m + S(1))*(-a*d + b*c)), x)
def replacement3796(c, x, d, B, A, f, a, C, e, b):
        # rubi.append(3796)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a*tan(e + f*x) + b)/(a + b*tan(e + f*x)), x), x) - Dist((A*d**S(2) - B*c*d + C*c**S(2))/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c*tan(e + f*x) + d)/(c + d*tan(e + f*x)), x), x) + Simp(x*(a*(A*c + B*d - C*c) + b*(-A*d + B*c + C*d))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3797(c, x, d, B, A, f, a, C, e, b):
        # rubi.append(3797)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a/tan(e + f*x) + b)/(a + b/tan(e + f*x)), x), x) - Dist((A*d**S(2) - B*c*d + C*c**S(2))/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c/tan(e + f*x) + d)/(c + d/tan(e + f*x)), x), x) + Simp(x*(a*(A*c + B*d - C*c) + b*(-A*d + B*c + C*d))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3798(c, x, d, A, f, a, C, e, b):
        # rubi.append(3798)
        return Dist((A*b**S(2) + C*a**S(2))/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a*tan(e + f*x) + b)/(a + b*tan(e + f*x)), x), x) - Dist((A*d**S(2) + C*c**S(2))/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c*tan(e + f*x) + d)/(c + d*tan(e + f*x)), x), x) + Simp(x*(a*(A*c - C*c) - b*(A*d - C*d))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3799(c, x, d, A, f, a, C, e, b):
        # rubi.append(3799)
        return Dist((A*b**S(2) + C*a**S(2))/((a**S(2) + b**S(2))*(-a*d + b*c)), Int((-a/tan(e + f*x) + b)/(a + b/tan(e + f*x)), x), x) - Dist((A*d**S(2) + C*c**S(2))/((c**S(2) + d**S(2))*(-a*d + b*c)), Int((-c/tan(e + f*x) + d)/(c + d/tan(e + f*x)), x), x) + Simp(x*(a*(A*c - C*c) - b*(A*d - C*d))/((a**S(2) + b**S(2))*(c**S(2) + d**S(2))), x)
def replacement3800(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3800)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*(tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*Simp(B*b + a*(A - C) + (B*a - b*(A - C))*tan(e + f*x), x), x), x)
def replacement3801(c, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3801)
        return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))*(c + d/tan(e + f*x))**n/(a + b/tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d/tan(e + f*x))**n*Simp(B*b + a*(A - C) + (B*a - b*(A - C))/tan(e + f*x), x), x), x)
def replacement3802(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3802)
        return Dist((A*b**S(2) + C*a**S(2))/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*(tan(e + f*x)**S(2) + S(1))/(a + b*tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d*tan(e + f*x))**n*Simp(a*(A - C) - (A*b - C*b)*tan(e + f*x), x), x), x)
def replacement3803(c, n, x, d, A, f, a, C, e, b):
        # rubi.append(3803)
        return Dist((A*b**S(2) + C*a**S(2))/(a**S(2) + b**S(2)), Int((S(1) + tan(e + f*x)**(S(-2)))*(c + d/tan(e + f*x))**n/(a + b/tan(e + f*x)), x), x) + Dist(S(1)/(a**S(2) + b**S(2)), Int((c + d/tan(e + f*x))**n*Simp(a*(A - C) - (A*b - C*b)/tan(e + f*x), x), x), x)
def replacement3804(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3804)
        return Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n*(A*b**S(2) + B*b*x + C*x**S(2))/(b**S(2) + x**S(2)), x), x, b*tan(e + f*x)), x)
def replacement3805(c, m, n, x, d, B, A, f, a, C, e, b):
        # rubi.append(3805)
        return -Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n*(A*b**S(2) + B*b*x + C*x**S(2))/(b**S(2) + x**S(2)), x), x, b/tan(e + f*x)), x)
def replacement3806(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3806)
        return Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n*(A*b**S(2) + C*x**S(2))/(b**S(2) + x**S(2)), x), x, b*tan(e + f*x)), x)
def replacement3807(c, m, n, x, d, A, f, a, C, e, b):
        # rubi.append(3807)
        return -Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n*(A*b**S(2) + C*x**S(2))/(b**S(2) + x**S(2)), x), x, b/tan(e + f*x)), x)
def replacement3808(c, x, d, a, b):
        # rubi.append(3808)
        return Dist(S(1)/a, Int(cos(c + d*x)**S(2), x), x)
def replacement3809(c, x, d, a, b):
        # rubi.append(3809)
        return Dist(S(1)/a, Int(sin(c + d*x)**S(2), x), x)
def replacement3810(c, x, d, a, b):
        # rubi.append(3810)
        return -Dist(b/(a - b), Int(S(1)/((a + b*tan(c + d*x)**S(2))*cos(c + d*x)**S(2)), x), x) + Simp(x/(a - b), x)
def replacement3811(c, x, d, a, b):
        # rubi.append(3811)
        return -Dist(b/(a - b), Int(S(1)/((a + b/tan(c + d*x)**S(2))*sin(c + d*x)**S(2)), x), x) + Simp(x/(a - b), x)
def replacement3812(c, n, x, d, a, p, e, b):
        # rubi.append(3812)
        return Dist(e/d, Subst(Int((a + b*x**n)**p/(e**S(2) + x**S(2)), x), x, e*tan(c + d*x)), x)
def replacement3813(c, n, x, d, a, p, e, b):
        # rubi.append(3813)
        return -Dist(e/d, Subst(Int((a + b*x**n)**p/(e**S(2) + x**S(2)), x), x, e/tan(c + d*x)), x)

def With3814(c, m, n, x, d, a, p, e, b):
        f = FreeFactors(tan(c + d*x), x)
        # rubi.append(3814)
        return Dist(f**(m + S(1))/d, Subst(Int(x**m*(a + b*(e*f*x)**n)**p*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)), x), x, tan(c + d*x)/f), x)

def With3815(c, m, n, x, d, a, p, e, b):
        f = FreeFactors(S(1)/tan(c + d*x), x)
        # rubi.append(3815)
        return -Dist(f**(m + S(1))/d, Subst(Int(x**m*(a + b*(e*f*x)**n)**p*(f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)), x), x, S(1)/(f*tan(c + d*x))), x)

def With3816(c, m, n, x, d, a, p, b):
        f = FreeFactors(cos(c + d*x), x)
        # rubi.append(3816)
        return -Dist(f/d, Subst(Int((f*x)**(-n*p)*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a*(f*x)**n + b*(-f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, cos(c + d*x)/f), x)

def With3817(c, m, n, x, d, a, p, b):
        f = FreeFactors(sin(c + d*x), x)
        # rubi.append(3817)
        return Dist(f/d, Subst(Int((f*x)**(-n*p)*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a*(f*x)**n + b*(-f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, sin(c + d*x)/f), x)

def With3818(c, m, n, x, d, a, p, e, b):
        f = FreeFactors(tan(c + d*x), x)
        # rubi.append(3818)
        return Dist(f/d, Subst(Int((a + b*(e*f*x)**n)**p*(f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)), x), x, tan(c + d*x)/f), x)

def With3819(c, m, n, x, d, a, p, e, b):
        f = FreeFactors(S(1)/tan(c + d*x), x)
        # rubi.append(3819)
        return -Dist(f/d, Subst(Int((a + b*(e*f*x)**n)**p*(f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)), x), x, S(1)/(f*tan(c + d*x))), x)

def With3820(c, m, n, x, d, a, p, b):
        f = FreeFactors(sin(c + d*x), x)
        # rubi.append(3820)
        return Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1)/2)*ExpandToSum(a*(-f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*(f*x)**n, x)**p, x), x, sin(c + d*x)/f), x)

def With3821(c, m, n, x, d, a, p, b):
        f = FreeFactors(cos(c + d*x), x)
        # rubi.append(3821)
        return -Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1)/2)*ExpandToSum(a*(-f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*(f*x)**n, x)**p, x), x, cos(c + d*x)/f), x)
def replacement3822(c, m, n, x, d, a, p, e, b):
        # rubi.append(3822)
        return Dist(e/d, Subst(Int((x/e)**m*(a + b*x**n)**p/(e**S(2) + x**S(2)), x), x, e*tan(c + d*x)), x)
def replacement3823(c, m, n, x, d, a, p, e, b):
        # rubi.append(3823)
        return -Dist(e/d, Subst(Int((x/e)**m*(a + b*x**n)**p/(e**S(2) + x**S(2)), x), x, e/tan(c + d*x)), x)
def replacement3824(c, n, x, d, n2, a, p, e, b):
        # rubi.append(3824)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p), x), x)
def replacement3825(c, n, x, d, n2, a, p, e, b):
        # rubi.append(3825)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p), x), x)
def replacement3826(c, n, x, d, n2, a, p, e, b):
        # rubi.append(3826)
        return Dist((b + S(2)*c*tan(d + e*x)**n)**(-S(2)*p)*(a + b*tan(d + e*x)**n + c*tan(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p), x), x)
def replacement3827(c, n, x, d, n2, a, p, e, b):
        # rubi.append(3827)
        return Dist((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/tan(d + e*x))**n + c*(S(1)/tan(d + e*x))**(S(2)*n))**p, Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p), x), x)

def With3828(c, n, x, d, n2, a, e, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(3828)
        return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*tan(d + e*x)**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*tan(d + e*x)**n + q), x), x)

def With3829(c, n, x, d, n2, a, e, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(3829)
        return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/tan(d + e*x))**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*(S(1)/tan(d + e*x))**n + q), x), x)
def replacement3830(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3830)
        return Dist(f/e, Subst(Int(x**m*(f**S(2) + x**S(2))**(-m/S(2) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, f*tan(d + e*x)), x)
def replacement3831(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3831)
        return -Dist(f/e, Subst(Int(x**m*(f**S(2) + x**S(2))**(-m/S(2) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, f/tan(d + e*x)), x)

def With3832(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(cos(d + e*x), x)
        # rubi.append(3832)
        return -Dist(g/e, Subst(Int((g*x)**(-S(2)*n*p)*(-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a*(g*x)**(S(2)*n) + b*(g*x)**n*(-g**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(-g**S(2)*x**S(2) + S(1))**n, x)**p, x), x, cos(d + e*x)/g), x)

def With3833(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(sin(d + e*x), x)
        # rubi.append(3833)
        return Dist(g/e, Subst(Int((g*x)**(-S(2)*n*p)*(-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a*(g*x)**(S(2)*n) + b*(g*x)**n*(-g**S(2)*x**S(2) + S(1))**(n/S(2)) + c*(-g**S(2)*x**S(2) + S(1))**n, x)**p, x), x, sin(d + e*x)/g), x)
def replacement3834(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3834)
        return Dist(f**(m + S(1))/e, Subst(Int((f**S(2) + x**S(2))**(-m/S(2) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, f*tan(d + e*x)), x)
def replacement3835(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3835)
        return -Dist(f**(m + S(1))/e, Subst(Int((f**S(2) + x**S(2))**(-m/S(2) + S(-1))*(a + b*x**n + c*x**(S(2)*n))**p, x), x, f/tan(d + e*x)), x)

def With3836(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(sin(d + e*x), x)
        # rubi.append(3836)
        return Dist(g/e, Subst(Int((-g**S(2)*x**S(2) + S(1))**(m/S(2) - n*p + S(-1)/2)*ExpandToSum(a*(S(1) - x**S(2))**n + b*x**n*(S(1) - x**S(2))**(n/S(2)) + c*x**(S(2)*n), x)**p, x), x, sin(d + e*x)/g), x)

def With3837(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(cos(d + e*x), x)
        # rubi.append(3837)
        return -Dist(g/e, Subst(Int((-g**S(2)*x**S(2) + S(1))**(m/S(2) - n*p + S(-1)/2)*ExpandToSum(a*(S(1) - x**S(2))**n + b*x**n*(S(1) - x**S(2))**(n/S(2)) + c*x**(S(2)*n), x)**p, x), x, cos(d + e*x)/g), x)
def replacement3838(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3838)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)
def replacement3839(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3839)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)
def replacement3840(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3840)
        return Dist((b + S(2)*c*tan(d + e*x)**n)**(-S(2)*p)*(a + b*tan(d + e*x)**n + c*tan(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)
def replacement3841(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3841)
        return Dist((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/tan(d + e*x))**n + c*(S(1)/tan(d + e*x))**(S(2)*n))**p, Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)
def replacement3842(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3842)
        return Dist(f/e, Subst(Int((x/f)**m*(a + b*x**n + c*x**(S(2)*n))**p/(f**S(2) + x**S(2)), x), x, f*tan(d + e*x)), x)
def replacement3843(c, m, n, x, d, f, n2, a, p, e, b):
        # rubi.append(3843)
        return -Dist(f/e, Subst(Int((x/f)**m*(a + b*x**n + c*x**(S(2)*n))**p/(f**S(2) + x**S(2)), x), x, f/tan(d + e*x)), x)
def replacement3844(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3844)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)
def replacement3845(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3845)
        return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p)*tan(d + e*x)**m, x), x)
def replacement3846(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3846)
        return Dist((b + S(2)*c*tan(d + e*x)**n)**(-S(2)*p)*(a + b*tan(d + e*x)**n + c*tan(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*tan(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)
def replacement3847(c, m, n, x, d, n2, a, p, e, b):
        # rubi.append(3847)
        return Dist((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(-S(2)*p)*(a + b*(S(1)/tan(d + e*x))**n + c*(S(1)/tan(d + e*x))**(S(2)*n))**p, Int((b + S(2)*c*(S(1)/tan(d + e*x))**n)**(S(2)*p)*tan(d + e*x)**m, x), x)

def With3848(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(S(1)/tan(d + e*x), x)
        # rubi.append(3848)
        return Dist(g/e, Subst(Int((g*x)**(m - S(2)*n*p)*(a*(g*x)**(S(2)*n) + b*(g*x)**n + c)**p/(g**S(2)*x**S(2) + S(1)), x), x, S(1)/(g*tan(d + e*x))), x)

def With3849(c, m, n, x, d, n2, a, p, e, b):
        g = FreeFactors(tan(d + e*x), x)
        # rubi.append(3849)
        return -Dist(g/e, Subst(Int((g*x)**(m - S(2)*n*p)*(a*(g*x)**(S(2)*n) + b*(g*x)**n + c)**p/(g**S(2)*x**S(2) + S(1)), x), x, tan(d + e*x)/g), x)
def replacement3850(c, n, x, d, B, A, a, e, b):
        # rubi.append(3850)
        return Dist(S(4)**(-n)*c**(-n), Int((A + B*tan(d + e*x))*(b + S(2)*c*tan(d + e*x))**(S(2)*n), x), x)
def replacement3851(c, n, x, d, B, A, a, e, b):
        # rubi.append(3851)
        return Dist(S(4)**(-n)*c**(-n), Int((A + B/tan(d + e*x))*(b + S(2)*c/tan(d + e*x))**(S(2)*n), x), x)
def replacement3852(c, n, x, d, B, A, a, e, b):
        # rubi.append(3852)
        return Dist((b + S(2)*c*tan(d + e*x))**(-S(2)*n)*(a + b*tan(d + e*x) + c*tan(d + e*x)**S(2))**n, Int((A + B*tan(d + e*x))*(b + S(2)*c*tan(d + e*x))**(S(2)*n), x), x)
def replacement3853(c, n, x, d, B, A, a, e, b):
        # rubi.append(3853)
        return Dist((b + S(2)*c/tan(d + e*x))**(-S(2)*n)*(a + b/tan(d + e*x) + c/tan(d + e*x)**S(2))**n, Int((A + B/tan(d + e*x))*(b + S(2)*c/tan(d + e*x))**(S(2)*n), x), x)

def With3854(c, x, d, B, A, a, e, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(3854)
        return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/Simp(b + S(2)*c*tan(d + e*x) - q, x), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/Simp(b + S(2)*c*tan(d + e*x) + q, x), x), x)

def With3855(c, x, d, B, A, a, e, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        # rubi.append(3855)
        return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/Simp(b + S(2)*c/tan(d + e*x) - q, x), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/Simp(b + S(2)*c/tan(d + e*x) + q, x), x), x)
def replacement3856(c, n, x, d, B, A, a, e, b):
        # rubi.append(3856)
        return Int(ExpandTrig((A + B*tan(d + e*x))*(a + b*tan(d + e*x) + c*tan(d + e*x)**S(2))**n, x), x)
def replacement3857(c, n, x, d, B, A, a, e, b):
        # rubi.append(3857)
        return Int(ExpandTrig((A + B/tan(d + e*x))*(a + b/tan(d + e*x) + c/tan(d + e*x)**S(2))**n, x), x)
def replacement3858(c, m, x, d, f, e):
        # rubi.append(3858)
        return -Dist(S(2)*I, Int((c + d*x)**m*exp(S(2)*I*(e + f*x))/(exp(S(2)*I*(e + f*x)) + S(1)), x), x) + Simp(I*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)
def replacement3859(c, m, x, d, f, e):
        # rubi.append(3859)
        return -Dist(S(2)*I, Int((c + d*x)**m*exp(S(2)*I*(e + f*x))/(S(1) - exp(S(2)*I*(e + f*x))), x), x) - Simp(I*(c + d*x)**(m + S(1))/(d*(m + S(1))), x)
def replacement3860(c, m, n, x, d, f, e, b):
        # rubi.append(3860)
        return -Dist(b**S(2), Int((b*tan(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) - Dist(b*d*m/(f*(n + S(-1))), Int((b*tan(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(-1)), x), x) + Simp(b*(b*tan(e + f*x))**(n + S(-1))*(c + d*x)**m/(f*(n + S(-1))), x)
def replacement3861(c, m, n, x, d, f, e, b):
        # rubi.append(3861)
        return -Dist(b**S(2), Int((b/tan(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) + Dist(b*d*m/(f*(n + S(-1))), Int((b/tan(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(-1)), x), x) - Simp(b*(b/tan(e + f*x))**(n + S(-1))*(c + d*x)**m/(f*(n + S(-1))), x)
def replacement3862(c, m, n, x, d, f, e, b):
        # rubi.append(3862)
        return -Dist(b**(S(-2)), Int((b*tan(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) - Dist(d*m/(b*f*(n + S(1))), Int((b*tan(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1)), x), x) + Simp((b*tan(e + f*x))**(n + S(1))*(c + d*x)**m/(b*f*(n + S(1))), x)
def replacement3863(c, m, n, x, d, f, e, b):
        # rubi.append(3863)
        return -Dist(b**(S(-2)), Int((b/tan(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) + Dist(d*m/(b*f*(n + S(1))), Int((b/tan(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1)), x), x) - Simp((b/tan(e + f*x))**(n + S(1))*(c + d*x)**m/(b*f*(n + S(1))), x)
def replacement3864(c, m, n, x, d, f, a, e, b):
        # rubi.append(3864)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b*tan(e + f*x))**n, x), x)
def replacement3865(c, m, n, x, d, f, a, e, b):
        # rubi.append(3865)
        return Int(ExpandIntegrand((c + d*x)**m, (a + b/tan(e + f*x))**n, x), x)
def replacement3866(c, m, x, d, f, a, e, b):
        # rubi.append(3866)
        return Dist(a*d*m/(S(2)*b*f), Int((c + d*x)**(m + S(-1))/(a + b*tan(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x) - Simp(a*(c + d*x)**m/(S(2)*b*f*(a + b*tan(e + f*x))), x)
def replacement3867(c, m, x, d, f, a, e, b):
        # rubi.append(3867)
        return -Dist(a*d*m/(S(2)*b*f), Int((c + d*x)**(m + S(-1))/(a + b/tan(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x) + Simp(a*(c + d*x)**m/(S(2)*b*f*(a + b/tan(e + f*x))), x)
def replacement3868(c, x, d, f, a, e, b):
        # rubi.append(3868)
        return -Dist(f/(a*d), Int(sin(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(f/(b*d), Int(cos(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Simp(S(1)/(d*(a + b*tan(e + f*x))*(c + d*x)), x)
def replacement3869(c, x, d, f, a, e, b):
        # rubi.append(3869)
        return Dist(f/(a*d), Int(sin(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(f/(b*d), Int(cos(S(2)*e + S(2)*f*x)/(c + d*x), x), x) - Simp(S(1)/(d*(a + b/tan(e + f*x))*(c + d*x)), x)
def replacement3870(c, m, x, d, f, a, e, b):
        # rubi.append(3870)
        return Dist(S(2)*b*f/(a*d*(m + S(1))), Int((c + d*x)**(m + S(1))/(a + b*tan(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + b*tan(e + f*x))*(m + S(1))), x) + Simp(f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement3871(c, m, x, d, f, a, e, b):
        # rubi.append(3871)
        return -Dist(S(2)*b*f/(a*d*(m + S(1))), Int((c + d*x)**(m + S(1))/(a + b/tan(e + f*x)), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + b/tan(e + f*x))*(m + S(1))), x) - Simp(f*(c + d*x)**(m + S(2))/(b*d**S(2)*(m + S(1))*(m + S(2))), x)
def replacement3872(c, x, d, f, a, e, b):
        # rubi.append(3872)
        return Dist(S(1)/(S(2)*a), Int(cos(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(S(1)/(S(2)*b), Int(sin(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Simp(log(c + d*x)/(S(2)*a*d), x)
def replacement3873(c, x, d, f, a, e, b):
        # rubi.append(3873)
        return -Dist(S(1)/(S(2)*a), Int(cos(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Dist(S(1)/(S(2)*b), Int(sin(S(2)*e + S(2)*f*x)/(c + d*x), x), x) + Simp(log(c + d*x)/(S(2)*a*d), x)
def replacement3874(c, m, x, d, f, a, e, b):
        # rubi.append(3874)
        return Dist(S(1)/(S(2)*a), Int((c + d*x)**m*exp(S(2)*a*(e + f*x)/b), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x)
def replacement3875(c, m, x, d, f, a, e, b):
        # rubi.append(3875)
        return -Dist(S(1)/(S(2)*a), Int((c + d*x)**m*exp(-S(2)*a*(e + f*x)/b), x), x) + Simp((c + d*x)**(m + S(1))/(S(2)*a*d*(m + S(1))), x)
def replacement3876(c, m, n, x, d, f, a, e, b):
        # rubi.append(3876)
        return Int(ExpandIntegrand((c + d*x)**m, (sin(S(2)*e + S(2)*f*x)/(S(2)*b) + cos(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x)
def replacement3877(c, m, n, x, d, f, a, e, b):
        # rubi.append(3877)
        return Int(ExpandIntegrand((c + d*x)**m, (sin(S(2)*e + S(2)*f*x)/(S(2)*b) - cos(S(2)*e + S(2)*f*x)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x)
def replacement3878(c, m, n, x, d, f, a, e, b):
        # rubi.append(3878)
        return Int(ExpandIntegrand((c + d*x)**m, (exp(S(2)*a*(e + f*x)/b)/(S(2)*a) + S(1)/(S(2)*a))**(-n), x), x)
def replacement3879(c, m, n, x, d, f, a, e, b):
        # rubi.append(3879)
        return Int(ExpandIntegrand((c + d*x)**m, (S(1)/(S(2)*a) - exp(-S(2)*a*(e + f*x)/b)/(S(2)*a))**(-n), x), x)

def With3880(c, m, n, x, d, f, a, e, b):
        u = IntHide((a + b*tan(e + f*x))**n, x)
        # rubi.append(3880)
        return -Dist(d*m, Int(Dist((c + d*x)**(m + S(-1)), u, x), x), x) + Dist((c + d*x)**m, u, x)

def With3881(c, m, n, x, d, f, a, e, b):
        u = IntHide((a + b/tan(e + f*x))**n, x)
        # rubi.append(3881)
        return -Dist(d*m, Int(Dist((c + d*x)**(m + S(-1)), u, x), x), x) + Dist((c + d*x)**m, u, x)
def replacement3882(c, m, x, d, f, a, e, b):
        # rubi.append(3882)
        return -Dist(S(2)*I*b, Int((c + d*x)**m/(a**S(2) + b**S(2) + (a - I*b)**S(2)*exp(S(2)*I*(e + f*x))), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a - I*b)*(m + S(1))), x)
def replacement3883(c, m, x, d, f, a, e, b):
        # rubi.append(3883)
        return Dist(S(2)*I*b, Int((c + d*x)**m/(a**S(2) + b**S(2) - (a + I*b)**S(2)*exp(S(2)*I*(e + f*x))), x), x) + Simp((c + d*x)**(m + S(1))/(d*(a + I*b)*(m + S(1))), x)
def replacement3884(c, x, d, f, a, e, b):
        # rubi.append(3884)
        return Dist(S(1)/(f*(a**S(2) + b**S(2))), Int((S(2)*a*c*f + S(2)*a*d*f*x + b*d)/(a + b*tan(e + f*x)), x), x) - Simp((c + d*x)**S(2)/(S(2)*d*(a**S(2) + b**S(2))), x) - Simp(b*(c + d*x)/(f*(a + b*tan(e + f*x))*(a**S(2) + b**S(2))), x)
def replacement3885(c, x, d, f, a, e, b):
        # rubi.append(3885)
        return -Dist(S(1)/(f*(a**S(2) + b**S(2))), Int((-S(2)*a*c*f - S(2)*a*d*f*x + b*d)/(a + b/tan(e + f*x)), x), x) - Simp((c + d*x)**S(2)/(S(2)*d*(a**S(2) + b**S(2))), x) + Simp(b*(c + d*x)/(f*(a + b/tan(e + f*x))*(a**S(2) + b**S(2))), x)
def replacement3886(c, m, n, x, d, f, a, e, b):
        # rubi.append(3886)
        return Int(ExpandIntegrand((c + d*x)**m, (-S(2)*I*b/(a**S(2) + b**S(2) + (a - I*b)**S(2)*exp(S(2)*I*(e + f*x))) + S(1)/(a - I*b))**(-n), x), x)
def replacement3887(c, m, n, x, d, f, a, e, b):
        # rubi.append(3887)
        return Int(ExpandIntegrand((c + d*x)**m, (S(2)*I*b/(a**S(2) + b**S(2) - (a + I*b)**S(2)*exp(S(2)*I*(e + f*x))) + S(1)/(a + I*b))**(-n), x), x)
def replacement3888(m, n, x, v, u, a, b):
        # rubi.append(3888)
        return Int((a + b*tan(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement3889(m, n, x, v, u, a, b):
        # rubi.append(3889)
        return Int((a + b/tan(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)
def replacement3890(c, m, n, x, d, f, a, e, b):
        # rubi.append(3890)
        return Int((a + b*tan(e + f*x))**n*(c + d*x)**m, x)
def replacement3891(c, m, n, x, d, f, a, e, b):
        # rubi.append(3891)
        return Int((a + b/tan(e + f*x))**n*(c + d*x)**m, x)
def replacement3892(c, n, x, d, a, p, b):
        # rubi.append(3892)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b*tan(c + d*x))**p, x), x, x**n), x)
def replacement3893(c, n, x, d, a, p, b):
        # rubi.append(3893)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + S(1)/n)*(a + b/tan(c + d*x))**p, x), x, x**n), x)
def replacement3894(c, n, x, d, a, p, b):
        # rubi.append(3894)
        return Int((a + b*tan(c + d*x**n))**p, x)
def replacement3895(c, n, x, d, a, p, b):
        # rubi.append(3895)
        return Int((a + b/tan(c + d*x**n))**p, x)
def replacement3896(c, n, x, d, u, a, p, b):
        # rubi.append(3896)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*tan(c + d*x**n))**p, x), x, u), x)
def replacement3897(c, n, x, d, u, a, p, b):
        # rubi.append(3897)
        return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b/tan(c + d*x**n))**p, x), x, u), x)
def replacement3898(x, u, a, p, b):
        # rubi.append(3898)
        return Int((a + b*tan(ExpandToSum(u, x)))**p, x)
def replacement3899(x, u, a, p, b):
        # rubi.append(3899)
        return Int((a + b/tan(ExpandToSum(u, x)))**p, x)
def replacement3900(c, m, n, x, d, a, p, b):
        # rubi.append(3900)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*tan(c + d*x))**p, x), x, x**n), x)
def replacement3901(c, m, n, x, d, a, p, b):
        # rubi.append(3901)
        return Dist(S(1)/n, Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b/tan(c + d*x))**p, x), x, x**n), x)
def replacement3902(c, m, n, x, d):
        # rubi.append(3902)
        return -Dist((m - n + S(1))/(d*n), Int(x**(m - n)*tan(c + d*x**n), x), x) - Int(x**m, x) + Simp(x**(m - n + S(1))*tan(c + d*x**n)/(d*n), x)
def replacement3903(c, m, n, x, d):
        # rubi.append(3903)
        return Dist((m - n + S(1))/(d*n), Int(x**(m - n)/tan(c + d*x**n), x), x) - Int(x**m, x) - Simp(x**(m - n + S(1))/(d*n*tan(c + d*x**n)), x)
def replacement3904(c, m, n, x, d, a, p, b):
        # rubi.append(3904)
        return Int(x**m*(a + b*tan(c + d*x**n))**p, x)
def replacement3905(c, m, n, x, d, a, p, b):
        # rubi.append(3905)
        return Int(x**m*(a + b/tan(c + d*x**n))**p, x)
def replacement3906(c, m, n, x, d, a, p, e, b):
        # rubi.append(3906)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*tan(c + d*x**n))**p, x), x)
def replacement3907(c, m, n, x, d, a, p, e, b):
        # rubi.append(3907)
        return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b/tan(c + d*x**n))**p, x), x)
def replacement3908(m, x, u, a, p, e, b):
        # rubi.append(3908)
        return Int((e*x)**m*(a + b*tan(ExpandToSum(u, x)))**p, x)
def replacement3909(m, x, u, a, p, e, b):
        # rubi.append(3909)
        return Int((e*x)**m*(a + b/tan(ExpandToSum(u, x)))**p, x)
def replacement3910(m, q, n, x, a, p, b):
        # rubi.append(3910)
        return -Dist((m - n + S(1))/(b*n*p), Int(x**(m - n)*(S(1)/cos(a + b*x**n))**p, x), x) + Simp(x**(m - n + S(1))*(S(1)/cos(a + b*x**n))**p/(b*n*p), x)
def replacement3911(m, q, n, x, a, p, b):
        # rubi.append(3911)
        return Dist((m - n + S(1))/(b*n*p), Int(x**(m - n)*(S(1)/sin(a + b*x**n))**p, x), x) - Simp(x**(m - n + S(1))*(S(1)/sin(a + b*x**n))**p/(b*n*p), x)
def replacement3912(c, n, x, a, b):
        # rubi.append(3912)
        return Int(tan(a + b*x + c*x**S(2))**n, x)
def replacement3913(c, n, x, a, b):
        # rubi.append(3913)
        return Int((S(1)/tan(a + b*x + c*x**S(2)))**n, x)
def replacement3914(c, x, d, a, e, b):
        # rubi.append(3914)
        return -Simp(e*log(cos(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement3915(c, x, d, a, e, b):
        # rubi.append(3915)
        return Simp(e*log(sin(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement3916(c, x, d, a, e, b):
        # rubi.append(3916)
        return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(tan(a + b*x + c*x**S(2)), x), x) - Simp(e*log(cos(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement3917(c, x, d, a, e, b):
        # rubi.append(3917)
        return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(S(1)/tan(a + b*x + c*x**S(2)), x), x) + Simp(e*log(sin(a + b*x + c*x**S(2)))/(S(2)*c), x)
def replacement3918(c, m, n, x, d, a, e, b):
        # rubi.append(3918)
        return Int((d + e*x)**m*tan(a + b*x + c*x**S(2))**n, x)
def replacement3919(c, m, n, x, d, a, e, b):
        # rubi.append(3919)
        return Int((d + e*x)**m*(S(1)/tan(a + b*x + c*x**S(2)))**n, x)
