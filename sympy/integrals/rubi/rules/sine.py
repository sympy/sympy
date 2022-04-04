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


def sine():
    from sympy.integrals.rubi.constraints import cons1251, cons74, cons68, cons2, cons3, cons50, cons127, cons19, cons4, cons1252, cons1253, cons1254, cons95, cons168, cons91, cons1255, cons33, cons1172, cons96, cons1256, cons1257, cons1258, cons1259, cons1260, cons1261, cons167, cons1262, cons21, cons25, cons523, cons8, cons29, cons810, cons1263, cons1264, cons1265, cons545, cons1266, cons1267, cons150, cons1268, cons89, cons45, cons450, cons1269, cons1270, cons1271, cons1272, cons1273, cons483, cons484, cons1274, cons1275, cons1276, cons1277, cons1278, cons210, cons5, cons20, cons13, cons139, cons1279, cons1280, cons1281, cons1282, cons1283, cons145, cons1284, cons1285, cons1286, cons1287, cons246, cons170, cons1288, cons1289, cons1290, cons148, cons1291, cons1292, cons1293, cons248, cons1294, cons1295, cons1296, cons1297, cons1298, cons86, cons1299, cons149, cons56, cons1300, cons1301, cons1302, cons1303, cons1304, cons64, cons269, cons1305, cons1306, cons1307, cons1308, cons517, cons274, cons1309, cons1310, cons73, cons72, cons1311, cons1312, cons1313, cons1314, cons348, cons1315, cons157, cons1316, cons1317, cons116, cons1318, cons1319, cons1320, cons1321, cons1322, cons1323, cons79, cons1324, cons1325, cons1326, cons1327, cons1328, cons1329, cons1330, cons465, cons90, cons1331, cons1332, cons87, cons1333, cons1334, cons1335, cons1336, cons1337, cons1338, cons1339, cons1340, cons1341, cons1342, cons1343, cons1344, cons1345, cons1346, cons1347, cons1348, cons1349, cons1350, cons1351, cons1352, cons1353, cons1354, cons1355, cons1356, cons1357, cons1358, cons1359, cons1360, cons1361, cons1362, cons1363, cons1364, cons1365, cons1366, cons144, cons337, cons1367, cons1368, cons1369, cons1370, cons1371, cons1372, cons172, cons1373, cons1374, cons1375, cons1376, cons1377, cons255, cons1378, cons1379, cons1380, cons1381, cons1382, cons360, cons1383, cons1384, cons1385, cons1386, cons1387, cons1388, cons1389, cons1390, cons1391, cons1392, cons1393, cons1394, cons1395, cons1396, cons215, cons586, cons1397, cons1398, cons1399, cons1400, cons1401, cons1402, cons1403, cons1404, cons1405, cons1406, cons1407, cons1408, cons1409, cons1410, cons1411, cons1412, cons1413, cons107, cons1414, cons1415, cons1416, cons1417, cons1418, cons1419, cons40, cons1420, cons36, cons37, cons1421, cons1422, cons1423, cons685, cons1424, cons1425, cons1426, cons1427, cons1428, cons1429, cons1430, cons1431, cons1432, cons1433, cons38, cons1230, cons1434, cons35, cons1435, cons1436, cons1437, cons1438, cons216, cons1439, cons1440, cons1441, cons1442, cons1443, cons1444, cons1445, cons1446, cons1447, cons1448, cons1154, cons1449, cons198, cons130, cons65, cons152, cons377, cons324, cons1450, cons1451, cons1452, cons1453, cons1454, cons1455, cons1456, cons78, cons1457, cons1458, cons1459, cons1460, cons1461, cons1462, cons1463, cons1464, cons1465, cons1466, cons1467, cons1468, cons1469, cons1470, cons1471, cons1472, cons1473, cons1247, cons1474, cons1475, cons1476, cons1477, cons1478, cons1479, cons1045, cons1480, cons1481, cons1482, cons1483, cons1484, cons1485, cons1486, cons1487, cons1488, cons1489, cons48, cons47, cons228, cons378, cons1118, cons178, cons1490, cons1491, cons247, cons249, cons1492, cons1493, cons1494, cons1495, cons812, cons813, cons746, cons1496, cons1497, cons55, cons598, cons1498, cons1499, cons491, cons1500, cons70, cons71, cons825, cons826, cons1501, cons1502, cons1503, cons1504, cons58, cons1505, cons1506, cons369, cons1507, cons358, cons856, cons1508, cons820, cons1133, cons49, cons241, cons1134, cons1135, cons1509, cons821


    pattern2167 = Pattern(Integral(u_, x_), cons1251)
    rule2167 = ReplacementRule(pattern2167, replacement2167)

    pattern2168 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons74, cons68)
    rule2168 = ReplacementRule(pattern2168, replacement2168)

    pattern2169 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1252, cons1253)
    rule2169 = ReplacementRule(pattern2169, replacement2169)

    pattern2170 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons2, cons50, cons127, cons19, cons1252, cons1254)
    rule2170 = ReplacementRule(pattern2170, replacement2170)

    pattern2171 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons168, cons91, cons1255)
    rule2171 = ReplacementRule(pattern2171, replacement2171)

    pattern2172 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons168, cons91, cons1255)
    rule2172 = ReplacementRule(pattern2172, replacement2172)

    pattern2173 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons168, cons1172)
    rule2173 = ReplacementRule(pattern2173, replacement2173)

    pattern2174 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons168, cons1172)
    rule2174 = ReplacementRule(pattern2174, replacement2174)

    pattern2175 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons96, cons1172)
    rule2175 = ReplacementRule(pattern2175, replacement2175)

    pattern2176 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons96, cons1172)
    rule2176 = ReplacementRule(pattern2176, replacement2176)

    pattern2177 = Pattern(Integral(sqrt(WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons1256)
    rule2177 = ReplacementRule(pattern2177, replacement2177)

    pattern2178 = Pattern(Integral(S(1)/(sqrt(WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons1256)
    rule2178 = ReplacementRule(pattern2178, replacement2178)

    pattern2179 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons1257, cons33, cons1258)
    rule2179 = ReplacementRule(pattern2179, With2179)

    pattern2180 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons1257, cons33, cons1258)
    rule2180 = ReplacementRule(pattern2180, With2180)

    pattern2181 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1259)
    rule2181 = ReplacementRule(pattern2181, replacement2181)

    pattern2182 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1260)
    rule2182 = ReplacementRule(pattern2182, replacement2182)

    pattern2183 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1261, cons68)
    rule2183 = ReplacementRule(pattern2183, replacement2183)

    pattern2184 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons4, cons1261, cons68)
    rule2184 = ReplacementRule(pattern2184, replacement2184)

    pattern2185 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons168, cons167, cons1172)
    rule2185 = ReplacementRule(pattern2185, replacement2185)

    pattern2186 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons95, cons168, cons167, cons1172)
    rule2186 = ReplacementRule(pattern2186, replacement2186)

    pattern2187 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons168, cons1262, cons1172)
    rule2187 = ReplacementRule(pattern2187, replacement2187)

    pattern2188 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons168, cons1262, cons1172)
    rule2188 = ReplacementRule(pattern2188, replacement2188)

    pattern2189 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons96, cons1172)
    rule2189 = ReplacementRule(pattern2189, replacement2189)

    pattern2190 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons4, cons33, cons96, cons1172)
    rule2190 = ReplacementRule(pattern2190, replacement2190)

    pattern2191 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule2191 = ReplacementRule(pattern2191, replacement2191)

    pattern2192 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule2192 = ReplacementRule(pattern2192, replacement2192)

    pattern2193 = Pattern(Integral((WC('a', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule2193 = ReplacementRule(pattern2193, replacement2193)

    pattern2194 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons50, cons127, cons19, cons4, cons21, cons25)
    rule2194 = ReplacementRule(pattern2194, replacement2194)

    pattern2195 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons8, cons29, cons523)
    rule2195 = ReplacementRule(pattern2195, replacement2195)

    pattern2196 = Pattern(Integral(cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons8, cons29, cons523)
    rule2196 = ReplacementRule(pattern2196, replacement2196)

    pattern2197 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons810, cons167)
    rule2197 = ReplacementRule(pattern2197, replacement2197)

    pattern2198 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons810, cons167)
    rule2198 = ReplacementRule(pattern2198, replacement2198)

    pattern2199 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons810, cons91)
    rule2199 = ReplacementRule(pattern2199, replacement2199)

    pattern2200 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons810, cons91)
    rule2200 = ReplacementRule(pattern2200, replacement2200)

    pattern2201 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons1263)
    rule2201 = ReplacementRule(pattern2201, replacement2201)

    pattern2202 = Pattern(Integral(cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons8, cons1264)
    rule2202 = ReplacementRule(pattern2202, replacement2202)

    pattern2203 = Pattern(Integral(sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons8, cons29, cons1263)
    rule2203 = ReplacementRule(pattern2203, replacement2203)

    pattern2204 = Pattern(Integral(sqrt(cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons8, cons29, cons1263)
    rule2204 = ReplacementRule(pattern2204, replacement2204)

    pattern2205 = Pattern(Integral(sqrt(b_*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons3, cons8, cons29, cons1265)
    rule2205 = ReplacementRule(pattern2205, replacement2205)

    pattern2206 = Pattern(Integral(sqrt(b_*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons3, cons8, cons29, cons1265)
    rule2206 = ReplacementRule(pattern2206, replacement2206)

    pattern2207 = Pattern(Integral(S(1)/sqrt(sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons8, cons29, cons1263)
    rule2207 = ReplacementRule(pattern2207, replacement2207)

    pattern2208 = Pattern(Integral(S(1)/sqrt(cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons8, cons29, cons1263)
    rule2208 = ReplacementRule(pattern2208, replacement2208)

    pattern2209 = Pattern(Integral(S(1)/sqrt(b_*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons3, cons8, cons29, cons1265)
    rule2209 = ReplacementRule(pattern2209, replacement2209)

    pattern2210 = Pattern(Integral(S(1)/sqrt(b_*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons3, cons8, cons29, cons1265)
    rule2210 = ReplacementRule(pattern2210, replacement2210)

    pattern2211 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons545)
    rule2211 = ReplacementRule(pattern2211, replacement2211)

    pattern2212 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons3, cons8, cons29, cons4, cons545)
    rule2212 = ReplacementRule(pattern2212, replacement2212)

    pattern2213 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule2213 = ReplacementRule(pattern2213, replacement2213)

    pattern2214 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons1266)
    rule2214 = ReplacementRule(pattern2214, replacement2214)

    pattern2215 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons150)
    rule2215 = ReplacementRule(pattern2215, replacement2215)

    pattern2216 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons150)
    rule2216 = ReplacementRule(pattern2216, replacement2216)

    pattern2217 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2217 = ReplacementRule(pattern2217, replacement2217)

    pattern2218 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2218 = ReplacementRule(pattern2218, replacement2218)

    pattern2219 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons1268)
    rule2219 = ReplacementRule(pattern2219, replacement2219)

    pattern2220 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons1268)
    rule2220 = ReplacementRule(pattern2220, replacement2220)

    pattern2221 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2221 = ReplacementRule(pattern2221, replacement2221)

    pattern2222 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2222 = ReplacementRule(pattern2222, replacement2222)

    pattern2223 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2223 = ReplacementRule(pattern2223, replacement2223)

    pattern2224 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1267)
    rule2224 = ReplacementRule(pattern2224, replacement2224)

    pattern2225 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons91, cons810)
    rule2225 = ReplacementRule(pattern2225, replacement2225)

    pattern2226 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1267, cons89, cons91, cons810)
    rule2226 = ReplacementRule(pattern2226, replacement2226)

    pattern2227 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons45)
    rule2227 = ReplacementRule(pattern2227, replacement2227)

    pattern2228 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons45)
    rule2228 = ReplacementRule(pattern2228, replacement2228)

    pattern2229 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons450)
    rule2229 = ReplacementRule(pattern2229, replacement2229)

    pattern2230 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1267, cons545, cons450)
    rule2230 = ReplacementRule(pattern2230, replacement2230)

    pattern2231 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1270)
    rule2231 = ReplacementRule(pattern2231, replacement2231)

    pattern2232 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1270)
    rule2232 = ReplacementRule(pattern2232, replacement2232)

    pattern2233 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1271)
    rule2233 = ReplacementRule(pattern2233, replacement2233)

    pattern2234 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1271)
    rule2234 = ReplacementRule(pattern2234, replacement2234)

    pattern2235 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1272)
    rule2235 = ReplacementRule(pattern2235, replacement2235)

    pattern2236 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1272)
    rule2236 = ReplacementRule(pattern2236, replacement2236)

    pattern2237 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons167, cons810)
    rule2237 = ReplacementRule(pattern2237, replacement2237)

    pattern2238 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons167, cons810)
    rule2238 = ReplacementRule(pattern2238, replacement2238)

    pattern2239 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1273, cons483)
    rule2239 = ReplacementRule(pattern2239, With2239)

    pattern2240 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1273, cons483)
    rule2240 = ReplacementRule(pattern2240, With2240)

    pattern2241 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1273, cons484)
    rule2241 = ReplacementRule(pattern2241, With2241)

    pattern2242 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1273, cons484)
    rule2242 = ReplacementRule(pattern2242, With2242)

    pattern2243 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1274)
    rule2243 = ReplacementRule(pattern2243, With2243)

    pattern2244 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule2244 = ReplacementRule(pattern2244, With2244)

    pattern2245 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269)
    rule2245 = ReplacementRule(pattern2245, With2245)

    pattern2246 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1270)
    rule2246 = ReplacementRule(pattern2246, replacement2246)

    pattern2247 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1270)
    rule2247 = ReplacementRule(pattern2247, replacement2247)

    pattern2248 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1271)
    rule2248 = ReplacementRule(pattern2248, replacement2248)

    pattern2249 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1271)
    rule2249 = ReplacementRule(pattern2249, replacement2249)

    pattern2250 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1272)
    rule2250 = ReplacementRule(pattern2250, replacement2250)

    pattern2251 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1269, cons1272)
    rule2251 = ReplacementRule(pattern2251, replacement2251)

    pattern2252 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons91, cons810)
    rule2252 = ReplacementRule(pattern2252, replacement2252)

    pattern2253 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1269, cons89, cons91, cons810)
    rule2253 = ReplacementRule(pattern2253, replacement2253)

    pattern2254 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1269, cons545)
    rule2254 = ReplacementRule(pattern2254, replacement2254)

    pattern2255 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1269, cons545)
    rule2255 = ReplacementRule(pattern2255, replacement2255)

    pattern2256 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1275)
    rule2256 = ReplacementRule(pattern2256, replacement2256)

    pattern2257 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1267, cons1277)
    rule2257 = ReplacementRule(pattern2257, replacement2257)

    pattern2258 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1267, cons1277)
    rule2258 = ReplacementRule(pattern2258, replacement2258)

    pattern2259 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1269)
    rule2259 = ReplacementRule(pattern2259, replacement2259)

    pattern2260 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1276, cons1269)
    rule2260 = ReplacementRule(pattern2260, replacement2260)

    pattern2261 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1278)
    rule2261 = ReplacementRule(pattern2261, replacement2261)

    pattern2262 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1278)
    rule2262 = ReplacementRule(pattern2262, replacement2262)

    pattern2263 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons20, cons13, cons139, cons1279)
    rule2263 = ReplacementRule(pattern2263, replacement2263)

    pattern2264 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons20, cons13, cons139, cons1279)
    rule2264 = ReplacementRule(pattern2264, replacement2264)

    pattern2265 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1280, cons1281)
    rule2265 = ReplacementRule(pattern2265, replacement2265)

    pattern2266 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1280, cons1281)
    rule2266 = ReplacementRule(pattern2266, replacement2266)

    pattern2267 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1282, cons1283, cons145)
    rule2267 = ReplacementRule(pattern2267, replacement2267)

    pattern2268 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1282, cons1283, cons145)
    rule2268 = ReplacementRule(pattern2268, replacement2268)

    pattern2269 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1284, cons1285)
    rule2269 = ReplacementRule(pattern2269, replacement2269)

    pattern2270 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1284, cons1285)
    rule2270 = ReplacementRule(pattern2270, replacement2270)

    pattern2271 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1286, cons1287)
    rule2271 = ReplacementRule(pattern2271, replacement2271)

    pattern2272 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1286, cons1287)
    rule2272 = ReplacementRule(pattern2272, replacement2272)

    pattern2273 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons170, cons1288, cons1289)
    rule2273 = ReplacementRule(pattern2273, replacement2273)

    pattern2274 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons170, cons1288, cons1289)
    rule2274 = ReplacementRule(pattern2274, replacement2274)

    pattern2275 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons168, cons139, cons1290)
    rule2275 = ReplacementRule(pattern2275, replacement2275)

    pattern2276 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons168, cons139, cons1290)
    rule2276 = ReplacementRule(pattern2276, replacement2276)

    pattern2277 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2277 = ReplacementRule(pattern2277, replacement2277)

    pattern2278 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2278 = ReplacementRule(pattern2278, replacement2278)

    pattern2279 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons33, cons170, cons1287, cons1290)
    rule2279 = ReplacementRule(pattern2279, replacement2279)

    pattern2280 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons33, cons170, cons1287, cons1290)
    rule2280 = ReplacementRule(pattern2280, replacement2280)

    pattern2281 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons96, cons148, cons1291, cons1287, cons1290)
    rule2281 = ReplacementRule(pattern2281, replacement2281)

    pattern2282 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons96, cons148, cons1291, cons1287, cons1290)
    rule2282 = ReplacementRule(pattern2282, replacement2282)

    pattern2283 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons1292, cons148, cons1283, cons1293, cons1290)
    rule2283 = ReplacementRule(pattern2283, replacement2283)

    pattern2284 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons246, cons1292, cons148, cons1283, cons1293, cons1290)
    rule2284 = ReplacementRule(pattern2284, replacement2284)

    pattern2285 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons33, cons96, cons1283, cons1290)
    rule2285 = ReplacementRule(pattern2285, replacement2285)

    pattern2286 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons33, cons96, cons1283, cons1290)
    rule2286 = ReplacementRule(pattern2286, replacement2286)

    pattern2287 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons148, cons248)
    rule2287 = ReplacementRule(pattern2287, replacement2287)

    pattern2288 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons148, cons248)
    rule2288 = ReplacementRule(pattern2288, replacement2288)

    pattern2289 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons1294, cons248)
    rule2289 = ReplacementRule(pattern2289, replacement2289)

    pattern2290 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons1294, cons248)
    rule2290 = ReplacementRule(pattern2290, replacement2290)

    pattern2291 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2291 = ReplacementRule(pattern2291, replacement2291)

    pattern2292 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2292 = ReplacementRule(pattern2292, replacement2292)

    pattern2293 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2293 = ReplacementRule(pattern2293, replacement2293)

    pattern2294 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267)
    rule2294 = ReplacementRule(pattern2294, replacement2294)

    pattern2295 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons1295, cons248)
    rule2295 = ReplacementRule(pattern2295, replacement2295)

    pattern2296 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons1295, cons248)
    rule2296 = ReplacementRule(pattern2296, replacement2296)

    pattern2297 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons139, cons248)
    rule2297 = ReplacementRule(pattern2297, replacement2297)

    pattern2298 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1267, cons13, cons139, cons248)
    rule2298 = ReplacementRule(pattern2298, replacement2298)

    pattern2299 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons20)
    rule2299 = ReplacementRule(pattern2299, replacement2299)

    pattern2300 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons20)
    rule2300 = ReplacementRule(pattern2300, replacement2300)

    pattern2301 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons21)
    rule2301 = ReplacementRule(pattern2301, replacement2301)

    pattern2302 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons21)
    rule2302 = ReplacementRule(pattern2302, replacement2302)

    pattern2303 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons1258, cons139, cons1296)
    rule2303 = ReplacementRule(pattern2303, replacement2303)

    pattern2304 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons1258, cons139, cons1296)
    rule2304 = ReplacementRule(pattern2304, replacement2304)

    pattern2305 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons168, cons139, cons1296)
    rule2305 = ReplacementRule(pattern2305, replacement2305)

    pattern2306 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons168, cons139, cons1296)
    rule2306 = ReplacementRule(pattern2306, replacement2306)

    pattern2307 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons33, cons168, cons1287, cons1296)
    rule2307 = ReplacementRule(pattern2307, replacement2307)

    pattern2308 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons33, cons168, cons1287, cons1296)
    rule2308 = ReplacementRule(pattern2308, replacement2308)

    pattern2309 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons96, cons148, cons1290)
    rule2309 = ReplacementRule(pattern2309, replacement2309)

    pattern2310 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons246, cons96, cons148, cons1290)
    rule2310 = ReplacementRule(pattern2310, replacement2310)

    pattern2311 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons33, cons96, cons1290)
    rule2311 = ReplacementRule(pattern2311, replacement2311)

    pattern2312 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons33, cons96, cons1290)
    rule2312 = ReplacementRule(pattern2312, replacement2312)

    pattern2313 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons1269, cons13, cons148, cons1287, cons1290)
    rule2313 = ReplacementRule(pattern2313, replacement2313)

    pattern2314 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons1269, cons13, cons148, cons1287, cons1290)
    rule2314 = ReplacementRule(pattern2314, replacement2314)

    pattern2315 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons1269, cons13, cons139, cons1290)
    rule2315 = ReplacementRule(pattern2315, replacement2315)

    pattern2316 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons1269, cons13, cons139, cons1290)
    rule2316 = ReplacementRule(pattern2316, replacement2316)

    pattern2317 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2317 = ReplacementRule(pattern2317, replacement2317)

    pattern2318 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2318 = ReplacementRule(pattern2318, replacement2318)

    pattern2319 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1280)
    rule2319 = ReplacementRule(pattern2319, replacement2319)

    pattern2320 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1280)
    rule2320 = ReplacementRule(pattern2320, replacement2320)

    pattern2321 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1297)
    rule2321 = ReplacementRule(pattern2321, replacement2321)

    pattern2322 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1297)
    rule2322 = ReplacementRule(pattern2322, replacement2322)

    pattern2323 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1298)
    rule2323 = ReplacementRule(pattern2323, replacement2323)

    pattern2324 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons1298)
    rule2324 = ReplacementRule(pattern2324, replacement2324)

    pattern2325 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2325 = ReplacementRule(pattern2325, With2325)

    pattern2326 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2326 = ReplacementRule(pattern2326, With2326)

    pattern2327 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2327 = ReplacementRule(pattern2327, With2327)

    pattern2328 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2328 = ReplacementRule(pattern2328, With2328)

    pattern2329 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons86, cons1299)
    rule2329 = ReplacementRule(pattern2329, replacement2329)

    pattern2330 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons86, cons1299)
    rule2330 = ReplacementRule(pattern2330, replacement2330)

    pattern2331 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons145)
    rule2331 = ReplacementRule(pattern2331, replacement2331)

    pattern2332 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1269, cons145)
    rule2332 = ReplacementRule(pattern2332, replacement2332)

    pattern2333 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule2333 = ReplacementRule(pattern2333, replacement2333)

    pattern2334 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule2334 = ReplacementRule(pattern2334, replacement2334)

    pattern2335 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons56)
    rule2335 = ReplacementRule(pattern2335, replacement2335)

    pattern2336 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons56)
    rule2336 = ReplacementRule(pattern2336, replacement2336)

    pattern2337 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1300)
    rule2337 = ReplacementRule(pattern2337, replacement2337)

    pattern2338 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1300)
    rule2338 = ReplacementRule(pattern2338, replacement2338)

    pattern2339 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons1267, cons1301, cons1302)
    rule2339 = ReplacementRule(pattern2339, replacement2339)

    pattern2340 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons50, cons127, cons1267, cons1301, cons1302)
    rule2340 = ReplacementRule(pattern2340, replacement2340)

    pattern2341 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons1267, cons1303, cons1304)
    rule2341 = ReplacementRule(pattern2341, replacement2341)

    pattern2342 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons50, cons127, cons1267, cons1303, cons1304)
    rule2342 = ReplacementRule(pattern2342, replacement2342)

    pattern2343 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons64)
    rule2343 = ReplacementRule(pattern2343, replacement2343)

    pattern2344 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons64)
    rule2344 = ReplacementRule(pattern2344, replacement2344)

    pattern2345 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons86)
    rule2345 = ReplacementRule(pattern2345, replacement2345)

    pattern2346 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons86)
    rule2346 = ReplacementRule(pattern2346, replacement2346)

    pattern2347 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons21, cons33, cons269)
    rule2347 = ReplacementRule(pattern2347, replacement2347)

    pattern2348 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons21, cons33, cons269)
    rule2348 = ReplacementRule(pattern2348, replacement2348)

    pattern2349 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons21, cons1305)
    rule2349 = ReplacementRule(pattern2349, replacement2349)

    pattern2350 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons21, cons1305)
    rule2350 = ReplacementRule(pattern2350, replacement2350)

    pattern2351 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306)
    rule2351 = ReplacementRule(pattern2351, replacement2351)

    pattern2352 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306)
    rule2352 = ReplacementRule(pattern2352, replacement2352)

    pattern2353 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons1306, cons96)
    rule2353 = ReplacementRule(pattern2353, replacement2353)

    pattern2354 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons1306, cons96)
    rule2354 = ReplacementRule(pattern2354, replacement2354)

    pattern2355 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306, cons1307)
    rule2355 = ReplacementRule(pattern2355, replacement2355)

    pattern2356 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306, cons1307)
    rule2356 = ReplacementRule(pattern2356, replacement2356)

    pattern2357 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons1267, cons1306, cons96)
    rule2357 = ReplacementRule(pattern2357, replacement2357)

    pattern2358 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons1267, cons1306, cons96)
    rule2358 = ReplacementRule(pattern2358, replacement2358)

    pattern2359 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306, cons1307)
    rule2359 = ReplacementRule(pattern2359, replacement2359)

    pattern2360 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1306, cons1307)
    rule2360 = ReplacementRule(pattern2360, replacement2360)

    pattern2361 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons21, cons1308)
    rule2361 = ReplacementRule(pattern2361, replacement2361)

    pattern2362 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons21, cons1308)
    rule2362 = ReplacementRule(pattern2362, replacement2362)

    pattern2363 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons21, cons149)
    rule2363 = ReplacementRule(pattern2363, replacement2363)

    pattern2364 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons21, cons149)
    rule2364 = ReplacementRule(pattern2364, replacement2364)

    pattern2365 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons1300)
    rule2365 = ReplacementRule(pattern2365, replacement2365)

    pattern2366 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons1300)
    rule2366 = ReplacementRule(pattern2366, replacement2366)

    pattern2367 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons64)
    rule2367 = ReplacementRule(pattern2367, replacement2367)

    pattern2368 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons64)
    rule2368 = ReplacementRule(pattern2368, replacement2368)

    pattern2369 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1269)
    rule2369 = ReplacementRule(pattern2369, replacement2369)

    pattern2370 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1269)
    rule2370 = ReplacementRule(pattern2370, replacement2370)

    pattern2371 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96, cons517)
    rule2371 = ReplacementRule(pattern2371, replacement2371)

    pattern2372 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons1269, cons33, cons96, cons517)
    rule2372 = ReplacementRule(pattern2372, replacement2372)

    pattern2373 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons274, cons517)
    rule2373 = ReplacementRule(pattern2373, replacement2373)

    pattern2374 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons274, cons517)
    rule2374 = ReplacementRule(pattern2374, replacement2374)

    pattern2375 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(6), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons1285, cons517)
    rule2375 = ReplacementRule(pattern2375, replacement2375)

    pattern2376 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**S(6), x_), cons2, cons3, cons50, cons127, cons19, cons1269, cons1285, cons517)
    rule2376 = ReplacementRule(pattern2376, replacement2376)

    pattern2377 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons1309, cons148)
    rule2377 = ReplacementRule(pattern2377, replacement2377)

    pattern2378 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons1309, cons148)
    rule2378 = ReplacementRule(pattern2378, replacement2378)

    pattern2379 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons1309, cons139)
    rule2379 = ReplacementRule(pattern2379, replacement2379)

    pattern2380 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons1309, cons139)
    rule2380 = ReplacementRule(pattern2380, replacement2380)

    pattern2381 = Pattern(Integral(sqrt(WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2381 = ReplacementRule(pattern2381, replacement2381)

    pattern2382 = Pattern(Integral(sqrt(WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2382 = ReplacementRule(pattern2382, replacement2382)

    pattern2383 = Pattern(Integral(S(1)/(sqrt(g_*tan(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2383 = ReplacementRule(pattern2383, replacement2383)

    pattern2384 = Pattern(Integral(S(1)/(sqrt(g_/tan(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2384 = ReplacementRule(pattern2384, replacement2384)

    pattern2385 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*tan(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons1269, cons1303)
    rule2385 = ReplacementRule(pattern2385, replacement2385)

    pattern2386 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(S(1)/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons50, cons127, cons1269, cons1303)
    rule2386 = ReplacementRule(pattern2386, replacement2386)

    pattern2387 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1310)
    rule2387 = ReplacementRule(pattern2387, replacement2387)

    pattern2388 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1310)
    rule2388 = ReplacementRule(pattern2388, replacement2388)

    pattern2389 = Pattern(Integral((WC('g', S(1))/tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule2389 = ReplacementRule(pattern2389, replacement2389)

    pattern2390 = Pattern(Integral((WC('g', S(1))*tan(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons149)
    rule2390 = ReplacementRule(pattern2390, replacement2390)

    pattern2391 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2391 = ReplacementRule(pattern2391, replacement2391)

    pattern2392 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2392 = ReplacementRule(pattern2392, replacement2392)

    pattern2393 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2393 = ReplacementRule(pattern2393, replacement2393)

    pattern2394 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2394 = ReplacementRule(pattern2394, replacement2394)

    pattern2395 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons20, cons1311)
    rule2395 = ReplacementRule(pattern2395, replacement2395)

    pattern2396 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons20, cons1311)
    rule2396 = ReplacementRule(pattern2396, replacement2396)

    pattern2397 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule2397 = ReplacementRule(pattern2397, replacement2397)

    pattern2398 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule2398 = ReplacementRule(pattern2398, replacement2398)

    pattern2399 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1312)
    rule2399 = ReplacementRule(pattern2399, replacement2399)

    pattern2400 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1312)
    rule2400 = ReplacementRule(pattern2400, replacement2400)

    pattern2401 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1313, cons89, cons91, cons1314)
    rule2401 = ReplacementRule(pattern2401, replacement2401)

    pattern2402 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267, cons1313, cons89, cons91, cons1314)
    rule2402 = ReplacementRule(pattern2402, replacement2402)

    pattern2403 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1313, cons348, cons1315, cons1314)
    rule2403 = ReplacementRule(pattern2403, replacement2403)

    pattern2404 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons1313, cons348, cons1315, cons1314)
    rule2404 = ReplacementRule(pattern2404, replacement2404)

    pattern2405 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule2405 = ReplacementRule(pattern2405, replacement2405)

    pattern2406 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons72, cons1267)
    rule2406 = ReplacementRule(pattern2406, replacement2406)

    pattern2407 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons157, cons1316)
    rule2407 = ReplacementRule(pattern2407, replacement2407)

    pattern2408 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons157, cons1316)
    rule2408 = ReplacementRule(pattern2408, replacement2408)

    pattern2409 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1317, cons1316, cons116)
    rule2409 = ReplacementRule(pattern2409, replacement2409)

    pattern2410 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1317, cons1316, cons116)
    rule2410 = ReplacementRule(pattern2410, replacement2410)

    pattern2411 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons33, cons96, cons1318, cons1172)
    rule2411 = ReplacementRule(pattern2411, replacement2411)

    pattern2412 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons72, cons1267, cons33, cons96, cons1318, cons1172)
    rule2412 = ReplacementRule(pattern2412, replacement2412)

    pattern2413 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1319)
    rule2413 = ReplacementRule(pattern2413, replacement2413)

    pattern2414 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons72, cons1267, cons1319)
    rule2414 = ReplacementRule(pattern2414, replacement2414)

    pattern2415 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2415 = ReplacementRule(pattern2415, replacement2415)

    pattern2416 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2)/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2416 = ReplacementRule(pattern2416, replacement2416)

    pattern2417 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2417 = ReplacementRule(pattern2417, replacement2417)

    pattern2418 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73)
    rule2418 = ReplacementRule(pattern2418, replacement2418)

    pattern2419 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1320)
    rule2419 = ReplacementRule(pattern2419, replacement2419)

    pattern2420 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1320)
    rule2420 = ReplacementRule(pattern2420, replacement2420)

    pattern2421 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1321)
    rule2421 = ReplacementRule(pattern2421, replacement2421)

    pattern2422 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1321)
    rule2422 = ReplacementRule(pattern2422, replacement2422)

    pattern2423 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons33, cons1322)
    rule2423 = ReplacementRule(pattern2423, replacement2423)

    pattern2424 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons33, cons1322)
    rule2424 = ReplacementRule(pattern2424, replacement2424)

    pattern2425 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1323)
    rule2425 = ReplacementRule(pattern2425, replacement2425)

    pattern2426 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1323)
    rule2426 = ReplacementRule(pattern2426, replacement2426)

    pattern2427 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2427 = ReplacementRule(pattern2427, replacement2427)

    pattern2428 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2428 = ReplacementRule(pattern2428, replacement2428)

    pattern2429 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons170, cons517)
    rule2429 = ReplacementRule(pattern2429, replacement2429)

    pattern2430 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons170, cons517)
    rule2430 = ReplacementRule(pattern2430, replacement2430)

    pattern2431 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons96, cons517)
    rule2431 = ReplacementRule(pattern2431, replacement2431)

    pattern2432 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons96, cons517)
    rule2432 = ReplacementRule(pattern2432, replacement2432)

    pattern2433 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269, cons79, cons1324)
    rule2433 = ReplacementRule(pattern2433, replacement2433)

    pattern2434 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269, cons79, cons1324)
    rule2434 = ReplacementRule(pattern2434, replacement2434)

    pattern2435 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269)
    rule2435 = ReplacementRule(pattern2435, replacement2435)

    pattern2436 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269)
    rule2436 = ReplacementRule(pattern2436, replacement2436)

    pattern2437 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons64, cons89)
    rule2437 = ReplacementRule(pattern2437, replacement2437)

    pattern2438 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons64, cons89)
    rule2438 = ReplacementRule(pattern2438, replacement2438)

    pattern2439 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule2439 = ReplacementRule(pattern2439, replacement2439)

    pattern2440 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons1267, cons33, cons1322)
    rule2440 = ReplacementRule(pattern2440, replacement2440)

    pattern2441 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule2441 = ReplacementRule(pattern2441, replacement2441)

    pattern2442 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons19, cons1267, cons1323)
    rule2442 = ReplacementRule(pattern2442, replacement2442)

    pattern2443 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons33, cons96)
    rule2443 = ReplacementRule(pattern2443, replacement2443)

    pattern2444 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons33, cons96)
    rule2444 = ReplacementRule(pattern2444, replacement2444)

    pattern2445 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons274)
    rule2445 = ReplacementRule(pattern2445, replacement2445)

    pattern2446 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons274)
    rule2446 = ReplacementRule(pattern2446, replacement2446)

    pattern2447 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons168, cons91, cons1326)
    rule2447 = ReplacementRule(pattern2447, replacement2447)

    pattern2448 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons168, cons91, cons1326)
    rule2448 = ReplacementRule(pattern2448, replacement2448)

    pattern2449 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons33, cons168, cons348, cons1326)
    rule2449 = ReplacementRule(pattern2449, replacement2449)

    pattern2450 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons33, cons168, cons348, cons1326)
    rule2450 = ReplacementRule(pattern2450, replacement2450)

    pattern2451 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons96, cons1327, cons1328)
    rule2451 = ReplacementRule(pattern2451, replacement2451)

    pattern2452 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons96, cons1327, cons1328)
    rule2452 = ReplacementRule(pattern2452, replacement2452)

    pattern2453 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons96, cons167, cons1328)
    rule2453 = ReplacementRule(pattern2453, replacement2453)

    pattern2454 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons95, cons96, cons167, cons1328)
    rule2454 = ReplacementRule(pattern2454, replacement2454)

    pattern2455 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons33, cons96, cons1329, cons1328)
    rule2455 = ReplacementRule(pattern2455, replacement2455)

    pattern2456 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons33, cons96, cons1329, cons1328)
    rule2456 = ReplacementRule(pattern2456, replacement2456)

    pattern2457 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons167, cons1330)
    rule2457 = ReplacementRule(pattern2457, replacement2457)

    pattern2458 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons167, cons1330)
    rule2458 = ReplacementRule(pattern2458, replacement2458)

    pattern2459 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons465, cons1330)
    rule2459 = ReplacementRule(pattern2459, replacement2459)

    pattern2460 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons465, cons1330)
    rule2460 = ReplacementRule(pattern2460, replacement2460)

    pattern2461 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons1330)
    rule2461 = ReplacementRule(pattern2461, replacement2461)

    pattern2462 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons1330)
    rule2462 = ReplacementRule(pattern2462, replacement2462)

    pattern2463 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons90, cons810)
    rule2463 = ReplacementRule(pattern2463, replacement2463)

    pattern2464 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons90, cons810)
    rule2464 = ReplacementRule(pattern2464, replacement2464)

    pattern2465 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2465 = ReplacementRule(pattern2465, replacement2465)

    pattern2466 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2466 = ReplacementRule(pattern2466, replacement2466)

    pattern2467 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons91, cons1331, cons810)
    rule2467 = ReplacementRule(pattern2467, replacement2467)

    pattern2468 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons91, cons810)
    rule2468 = ReplacementRule(pattern2468, replacement2468)

    pattern2469 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2469 = ReplacementRule(pattern2469, replacement2469)

    pattern2470 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2470 = ReplacementRule(pattern2470, replacement2470)

    pattern2471 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332)
    rule2471 = ReplacementRule(pattern2471, replacement2471)

    pattern2472 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332)
    rule2472 = ReplacementRule(pattern2472, replacement2472)

    pattern2473 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2473 = ReplacementRule(pattern2473, replacement2473)

    pattern2474 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2474 = ReplacementRule(pattern2474, replacement2474)

    pattern2475 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons545)
    rule2475 = ReplacementRule(pattern2475, replacement2475)

    pattern2476 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons545)
    rule2476 = ReplacementRule(pattern2476, replacement2476)

    pattern2477 = Pattern(Integral(sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2477 = ReplacementRule(pattern2477, replacement2477)

    pattern2478 = Pattern(Integral(sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2478 = ReplacementRule(pattern2478, replacement2478)

    pattern2479 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons167, cons810)
    rule2479 = ReplacementRule(pattern2479, replacement2479)

    pattern2480 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons167, cons810)
    rule2480 = ReplacementRule(pattern2480, replacement2480)

    pattern2481 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons91, cons810)
    rule2481 = ReplacementRule(pattern2481, replacement2481)

    pattern2482 = Pattern(Integral((WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_/sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325, cons89, cons91, cons810)
    rule2482 = ReplacementRule(pattern2482, replacement2482)

    pattern2483 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2483 = ReplacementRule(pattern2483, replacement2483)

    pattern2484 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2484 = ReplacementRule(pattern2484, replacement2484)

    pattern2485 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332, cons45)
    rule2485 = ReplacementRule(pattern2485, replacement2485)

    pattern2486 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1332, cons45)
    rule2486 = ReplacementRule(pattern2486, replacement2486)

    pattern2487 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2487 = ReplacementRule(pattern2487, replacement2487)

    pattern2488 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1325)
    rule2488 = ReplacementRule(pattern2488, replacement2488)

    pattern2489 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1325, cons89, cons167, cons87)
    rule2489 = ReplacementRule(pattern2489, replacement2489)

    pattern2490 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1267, cons1325, cons89, cons167, cons87)
    rule2490 = ReplacementRule(pattern2490, replacement2490)

    pattern2491 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons20)
    rule2491 = ReplacementRule(pattern2491, replacement2491)

    pattern2492 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1267, cons1325, cons20)
    rule2492 = ReplacementRule(pattern2492, replacement2492)

    pattern2493 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons1333)
    rule2493 = ReplacementRule(pattern2493, replacement2493)

    pattern2494 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons1333)
    rule2494 = ReplacementRule(pattern2494, replacement2494)

    pattern2495 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons1334)
    rule2495 = ReplacementRule(pattern2495, replacement2495)

    pattern2496 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons45, cons1334)
    rule2496 = ReplacementRule(pattern2496, replacement2496)

    pattern2497 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons450)
    rule2497 = ReplacementRule(pattern2497, replacement2497)

    pattern2498 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons21, cons450)
    rule2498 = ReplacementRule(pattern2498, replacement2498)

    pattern2499 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1267, cons1325, cons21)
    rule2499 = ReplacementRule(pattern2499, replacement2499)

    pattern2500 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1267, cons1325, cons21)
    rule2500 = ReplacementRule(pattern2500, replacement2500)

    pattern2501 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1320)
    rule2501 = ReplacementRule(pattern2501, replacement2501)

    pattern2502 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons3, cons8, cons29, cons50, cons127, cons19, cons1320)
    rule2502 = ReplacementRule(pattern2502, replacement2502)

    pattern2503 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons96)
    rule2503 = ReplacementRule(pattern2503, replacement2503)

    pattern2504 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons33, cons96)
    rule2504 = ReplacementRule(pattern2504, replacement2504)

    pattern2505 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269, cons274)
    rule2505 = ReplacementRule(pattern2505, replacement2505)

    pattern2506 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons73, cons1269, cons274)
    rule2506 = ReplacementRule(pattern2506, replacement2506)

    pattern2507 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons1335, cons91, cons1336)
    rule2507 = ReplacementRule(pattern2507, replacement2507)

    pattern2508 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons1335, cons91, cons1336)
    rule2508 = ReplacementRule(pattern2508, replacement2508)

    pattern2509 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1269, cons1325, cons33, cons1335, cons1336, cons1337)
    rule2509 = ReplacementRule(pattern2509, replacement2509)

    pattern2510 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1269, cons1325, cons33, cons1335, cons1336, cons1337)
    rule2510 = ReplacementRule(pattern2510, replacement2510)

    pattern2511 = Pattern(Integral(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2511 = ReplacementRule(pattern2511, replacement2511)

    pattern2512 = Pattern(Integral(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2512 = ReplacementRule(pattern2512, replacement2512)

    pattern2513 = Pattern(Integral(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2513 = ReplacementRule(pattern2513, replacement2513)

    pattern2514 = Pattern(Integral(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2514 = ReplacementRule(pattern2514, replacement2514)

    pattern2515 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons96, cons1327, cons1172)
    rule2515 = ReplacementRule(pattern2515, replacement2515)

    pattern2516 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons96, cons1327, cons1172)
    rule2516 = ReplacementRule(pattern2516, replacement2516)

    pattern2517 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2517 = ReplacementRule(pattern2517, replacement2517)

    pattern2518 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2518 = ReplacementRule(pattern2518, replacement2518)

    pattern2519 = Pattern(Integral((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2519 = ReplacementRule(pattern2519, replacement2519)

    pattern2520 = Pattern(Integral((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2520 = ReplacementRule(pattern2520, replacement2520)

    pattern2521 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons96, cons1338, cons1172)
    rule2521 = ReplacementRule(pattern2521, replacement2521)

    pattern2522 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons96, cons1338, cons1172)
    rule2522 = ReplacementRule(pattern2522, replacement2522)

    pattern2523 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2523 = ReplacementRule(pattern2523, replacement2523)

    pattern2524 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2524 = ReplacementRule(pattern2524, replacement2524)

    pattern2525 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2525 = ReplacementRule(pattern2525, replacement2525)

    pattern2526 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2526 = ReplacementRule(pattern2526, replacement2526)

    pattern2527 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1172, cons1339)
    rule2527 = ReplacementRule(pattern2527, replacement2527)

    pattern2528 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1172, cons1339)
    rule2528 = ReplacementRule(pattern2528, replacement2528)

    pattern2529 = Pattern(Integral(sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2529 = ReplacementRule(pattern2529, replacement2529)

    pattern2530 = Pattern(Integral(sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2530 = ReplacementRule(pattern2530, replacement2530)

    pattern2531 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2531 = ReplacementRule(pattern2531, replacement2531)

    pattern2532 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2532 = ReplacementRule(pattern2532, replacement2532)

    pattern2533 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1340)
    rule2533 = ReplacementRule(pattern2533, replacement2533)

    pattern2534 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1340)
    rule2534 = ReplacementRule(pattern2534, replacement2534)

    pattern2535 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1341)
    rule2535 = ReplacementRule(pattern2535, replacement2535)

    pattern2536 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1341)
    rule2536 = ReplacementRule(pattern2536, replacement2536)

    pattern2537 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1342)
    rule2537 = ReplacementRule(pattern2537, replacement2537)

    pattern2538 = Pattern(Integral(S(1)/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1342)
    rule2538 = ReplacementRule(pattern2538, replacement2538)

    pattern2539 = Pattern(Integral(sqrt(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1343, cons1344, cons1345)
    rule2539 = ReplacementRule(pattern2539, replacement2539)

    pattern2540 = Pattern(Integral(sqrt(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1343, cons1344, cons1345)
    rule2540 = ReplacementRule(pattern2540, replacement2540)

    pattern2541 = Pattern(Integral(sqrt(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1344)
    rule2541 = ReplacementRule(pattern2541, replacement2541)

    pattern2542 = Pattern(Integral(sqrt(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1344)
    rule2542 = ReplacementRule(pattern2542, replacement2542)

    pattern2543 = Pattern(Integral(sqrt(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1346)
    rule2543 = ReplacementRule(pattern2543, replacement2543)

    pattern2544 = Pattern(Integral(sqrt(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons3, cons8, cons29, cons50, cons127, cons1325, cons1346)
    rule2544 = ReplacementRule(pattern2544, replacement2544)

    pattern2545 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1347)
    rule2545 = ReplacementRule(pattern2545, replacement2545)

    pattern2546 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1347)
    rule2546 = ReplacementRule(pattern2546, replacement2546)

    pattern2547 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1348)
    rule2547 = ReplacementRule(pattern2547, replacement2547)

    pattern2548 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1348)
    rule2548 = ReplacementRule(pattern2548, replacement2548)

    pattern2549 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1349, cons1350, cons1351)
    rule2549 = ReplacementRule(pattern2549, replacement2549)

    pattern2550 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1349, cons1350, cons1351)
    rule2550 = ReplacementRule(pattern2550, replacement2550)

    pattern2551 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1349, cons1352, cons1353)
    rule2551 = ReplacementRule(pattern2551, replacement2551)

    pattern2552 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1349, cons1352, cons1353)
    rule2552 = ReplacementRule(pattern2552, replacement2552)

    pattern2553 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1273, cons1354, cons1355)
    rule2553 = ReplacementRule(pattern2553, replacement2553)

    pattern2554 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1273, cons1354, cons1355)
    rule2554 = ReplacementRule(pattern2554, replacement2554)

    pattern2555 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1354)
    rule2555 = ReplacementRule(pattern2555, replacement2555)

    pattern2556 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1354)
    rule2556 = ReplacementRule(pattern2556, replacement2556)

    pattern2557 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1356)
    rule2557 = ReplacementRule(pattern2557, replacement2557)

    pattern2558 = Pattern(Integral(S(1)/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1356)
    rule2558 = ReplacementRule(pattern2558, replacement2558)

    pattern2559 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1357)
    rule2559 = ReplacementRule(pattern2559, replacement2559)

    pattern2560 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1357)
    rule2560 = ReplacementRule(pattern2560, replacement2560)

    pattern2561 = Pattern(Integral(S(1)/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1358)
    rule2561 = ReplacementRule(pattern2561, replacement2561)

    pattern2562 = Pattern(Integral(S(1)/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons1358)
    rule2562 = ReplacementRule(pattern2562, replacement2562)

    pattern2563 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2563 = ReplacementRule(pattern2563, replacement2563)

    pattern2564 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)/sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2564 = ReplacementRule(pattern2564, replacement2564)

    pattern2565 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons1359, cons1360, cons1361, cons1336)
    rule2565 = ReplacementRule(pattern2565, replacement2565)

    pattern2566 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325, cons95, cons1359, cons1360, cons1361, cons1336)
    rule2566 = ReplacementRule(pattern2566, replacement2566)

    pattern2567 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons64)
    rule2567 = ReplacementRule(pattern2567, replacement2567)

    pattern2568 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons64)
    rule2568 = ReplacementRule(pattern2568, replacement2568)

    pattern2569 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1269, cons1325)
    rule2569 = ReplacementRule(pattern2569, replacement2569)

    pattern2570 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons73, cons1269, cons1325)
    rule2570 = ReplacementRule(pattern2570, replacement2570)

    pattern2571 = Pattern(Integral(((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25)
    rule2571 = ReplacementRule(pattern2571, replacement2571)

    pattern2572 = Pattern(Integral(((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*WC('c', S(1)))**n_*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons25)
    rule2572 = ReplacementRule(pattern2572, replacement2572)

    pattern2573 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons87)
    rule2573 = ReplacementRule(pattern2573, replacement2573)

    pattern2574 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons87)
    rule2574 = ReplacementRule(pattern2574, replacement2574)

    pattern2575 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20)
    rule2575 = ReplacementRule(pattern2575, replacement2575)

    pattern2576 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20)
    rule2576 = ReplacementRule(pattern2576, replacement2576)

    pattern2577 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons25, cons21)
    rule2577 = ReplacementRule(pattern2577, replacement2577)

    pattern2578 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons25, cons21)
    rule2578 = ReplacementRule(pattern2578, replacement2578)

    pattern2579 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule2579 = ReplacementRule(pattern2579, replacement2579)

    pattern2580 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule2580 = ReplacementRule(pattern2580, replacement2580)

    pattern2581 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1276, cons87, cons1363)
    rule2581 = ReplacementRule(pattern2581, replacement2581)

    pattern2582 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1276, cons87, cons1363)
    rule2582 = ReplacementRule(pattern2582, replacement2582)

    pattern2583 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1276, cons1267, cons87, cons1364)
    rule2583 = ReplacementRule(pattern2583, replacement2583)

    pattern2584 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1276, cons1267, cons87, cons1364)
    rule2584 = ReplacementRule(pattern2584, replacement2584)

    pattern2585 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons8, cons29, cons19, cons4, cons1276, cons1267)
    rule2585 = ReplacementRule(pattern2585, replacement2585)

    pattern2586 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons50, cons127, cons8, cons29, cons19, cons4, cons1276, cons1267)
    rule2586 = ReplacementRule(pattern2586, replacement2586)

    pattern2587 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1276, cons1269)
    rule2587 = ReplacementRule(pattern2587, replacement2587)

    pattern2588 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1276, cons1269)
    rule2588 = ReplacementRule(pattern2588, replacement2588)

    pattern2589 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1365)
    rule2589 = ReplacementRule(pattern2589, replacement2589)

    pattern2590 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1365)
    rule2590 = ReplacementRule(pattern2590, replacement2590)

    pattern2591 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267)
    rule2591 = ReplacementRule(pattern2591, replacement2591)

    pattern2592 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267)
    rule2592 = ReplacementRule(pattern2592, replacement2592)

    pattern2593 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons20, cons1366)
    rule2593 = ReplacementRule(pattern2593, replacement2593)

    pattern2594 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons20, cons1366)
    rule2594 = ReplacementRule(pattern2594, replacement2594)

    pattern2595 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons72, cons1267, cons1308)
    rule2595 = ReplacementRule(pattern2595, replacement2595)

    pattern2596 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons72, cons1267, cons1308)
    rule2596 = ReplacementRule(pattern2596, replacement2596)

    pattern2597 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267)
    rule2597 = ReplacementRule(pattern2597, replacement2597)

    pattern2598 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267)
    rule2598 = ReplacementRule(pattern2598, replacement2598)

    pattern2599 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1284, cons144)
    rule2599 = ReplacementRule(pattern2599, replacement2599)

    pattern2600 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1284, cons144)
    rule2600 = ReplacementRule(pattern2600, replacement2600)

    pattern2601 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1284, cons337)
    rule2601 = ReplacementRule(pattern2601, replacement2601)

    pattern2602 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1284, cons337)
    rule2602 = ReplacementRule(pattern2602, replacement2602)

    pattern2603 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons1286, cons89, cons91, cons1367, cons1368)
    rule2603 = ReplacementRule(pattern2603, replacement2603)

    pattern2604 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons1286, cons89, cons91, cons1367, cons1368)
    rule2604 = ReplacementRule(pattern2604, replacement2604)

    pattern2605 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons1286, cons348, cons1369, cons1368)
    rule2605 = ReplacementRule(pattern2605, replacement2605)

    pattern2606 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons1286, cons348, cons1369, cons1368)
    rule2606 = ReplacementRule(pattern2606, replacement2606)

    pattern2607 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons72, cons1267, cons1370)
    rule2607 = ReplacementRule(pattern2607, replacement2607)

    pattern2608 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons72, cons1267, cons1370)
    rule2608 = ReplacementRule(pattern2608, replacement2608)

    pattern2609 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1371, cons1262)
    rule2609 = ReplacementRule(pattern2609, replacement2609)

    pattern2610 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1371, cons1262)
    rule2610 = ReplacementRule(pattern2610, replacement2610)

    pattern2611 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1372, cons1283, cons116)
    rule2611 = ReplacementRule(pattern2611, replacement2611)

    pattern2612 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1372, cons1283, cons116)
    rule2612 = ReplacementRule(pattern2612, replacement2612)

    pattern2613 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons95, cons170, cons91, cons1367, cons172)
    rule2613 = ReplacementRule(pattern2613, replacement2613)

    pattern2614 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons72, cons1267, cons95, cons170, cons91, cons1367, cons172)
    rule2614 = ReplacementRule(pattern2614, replacement2614)

    pattern2615 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons33, cons170, cons1373, cons1374, cons172)
    rule2615 = ReplacementRule(pattern2615, replacement2615)

    pattern2616 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons33, cons170, cons1373, cons1374, cons172)
    rule2616 = ReplacementRule(pattern2616, replacement2616)

    pattern2617 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons33, cons96, cons1283, cons1318, cons172)
    rule2617 = ReplacementRule(pattern2617, replacement2617)

    pattern2618 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons72, cons1267, cons33, cons96, cons1283, cons1318, cons172)
    rule2618 = ReplacementRule(pattern2618, replacement2618)

    pattern2619 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1319)
    rule2619 = ReplacementRule(pattern2619, replacement2619)

    pattern2620 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons72, cons1267, cons1319)
    rule2620 = ReplacementRule(pattern2620, replacement2620)

    pattern2621 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1375)
    rule2621 = ReplacementRule(pattern2621, replacement2621)

    pattern2622 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1375)
    rule2622 = ReplacementRule(pattern2622, replacement2622)

    pattern2623 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1267, cons246, cons1376, cons139)
    rule2623 = ReplacementRule(pattern2623, replacement2623)

    pattern2624 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1267, cons246, cons1376, cons139)
    rule2624 = ReplacementRule(pattern2624, replacement2624)

    pattern2625 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1377, cons255)
    rule2625 = ReplacementRule(pattern2625, replacement2625)

    pattern2626 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1377, cons255)
    rule2626 = ReplacementRule(pattern2626, replacement2626)

    pattern2627 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons1378)
    rule2627 = ReplacementRule(pattern2627, replacement2627)

    pattern2628 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons1378)
    rule2628 = ReplacementRule(pattern2628, replacement2628)

    pattern2629 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons1379)
    rule2629 = ReplacementRule(pattern2629, replacement2629)

    pattern2630 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1267, cons33, cons1379)
    rule2630 = ReplacementRule(pattern2630, replacement2630)

    pattern2631 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1380, cons1283)
    rule2631 = ReplacementRule(pattern2631, replacement2631)

    pattern2632 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons1380, cons1283)
    rule2632 = ReplacementRule(pattern2632, replacement2632)

    pattern2633 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons255)
    rule2633 = ReplacementRule(pattern2633, replacement2633)

    pattern2634 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons1267, cons255)
    rule2634 = ReplacementRule(pattern2634, replacement2634)

    pattern2635 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269, cons246, cons170, cons139, cons517)
    rule2635 = ReplacementRule(pattern2635, replacement2635)

    pattern2636 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269, cons246, cons170, cons139, cons517)
    rule2636 = ReplacementRule(pattern2636, replacement2636)

    pattern2637 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons33, cons170, cons1381, cons517)
    rule2637 = ReplacementRule(pattern2637, replacement2637)

    pattern2638 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons33, cons170, cons1381, cons517)
    rule2638 = ReplacementRule(pattern2638, replacement2638)

    pattern2639 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269, cons246, cons96, cons148, cons255, cons517)
    rule2639 = ReplacementRule(pattern2639, replacement2639)

    pattern2640 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269, cons246, cons96, cons148, cons255, cons517)
    rule2640 = ReplacementRule(pattern2640, replacement2640)

    pattern2641 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons33, cons96, cons517)
    rule2641 = ReplacementRule(pattern2641, replacement2641)

    pattern2642 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons33, cons96, cons517)
    rule2642 = ReplacementRule(pattern2642, replacement2642)

    pattern2643 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons1269, cons13, cons148, cons1287, cons255, cons517)
    rule2643 = ReplacementRule(pattern2643, replacement2643)

    pattern2644 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons1269, cons13, cons148, cons1287, cons255, cons517)
    rule2644 = ReplacementRule(pattern2644, replacement2644)

    pattern2645 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons1269, cons13, cons139, cons517)
    rule2645 = ReplacementRule(pattern2645, replacement2645)

    pattern2646 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons1269, cons13, cons139, cons517)
    rule2646 = ReplacementRule(pattern2646, replacement2646)

    pattern2647 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269)
    rule2647 = ReplacementRule(pattern2647, replacement2647)

    pattern2648 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons1269)
    rule2648 = ReplacementRule(pattern2648, replacement2648)

    pattern2649 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1269, cons1324)
    rule2649 = ReplacementRule(pattern2649, replacement2649)

    pattern2650 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons1269, cons1324)
    rule2650 = ReplacementRule(pattern2650, replacement2650)

    pattern2651 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons1301, cons1382)
    rule2651 = ReplacementRule(pattern2651, replacement2651)

    pattern2652 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons1301, cons1382)
    rule2652 = ReplacementRule(pattern2652, replacement2652)

    pattern2653 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons360)
    rule2653 = ReplacementRule(pattern2653, replacement2653)

    pattern2654 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons360)
    rule2654 = ReplacementRule(pattern2654, replacement2654)

    pattern2655 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1280)
    rule2655 = ReplacementRule(pattern2655, replacement2655)

    pattern2656 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1280)
    rule2656 = ReplacementRule(pattern2656, replacement2656)

    pattern2657 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1383, cons1384)
    rule2657 = ReplacementRule(pattern2657, replacement2657)

    pattern2658 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons1267, cons1383, cons1384)
    rule2658 = ReplacementRule(pattern2658, replacement2658)

    pattern2659 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64)
    rule2659 = ReplacementRule(pattern2659, replacement2659)

    pattern2660 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64)
    rule2660 = ReplacementRule(pattern2660, replacement2660)

    pattern2661 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1385)
    rule2661 = ReplacementRule(pattern2661, replacement2661)

    pattern2662 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1385)
    rule2662 = ReplacementRule(pattern2662, replacement2662)

    pattern2663 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons86)
    rule2663 = ReplacementRule(pattern2663, replacement2663)

    pattern2664 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons86)
    rule2664 = ReplacementRule(pattern2664, replacement2664)

    pattern2665 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons1267, cons20, cons13, cons1386)
    rule2665 = ReplacementRule(pattern2665, replacement2665)

    pattern2666 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons1267, cons20, cons13, cons1386)
    rule2666 = ReplacementRule(pattern2666, replacement2666)

    pattern2667 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons33, cons1387, cons1283)
    rule2667 = ReplacementRule(pattern2667, replacement2667)

    pattern2668 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1267, cons33, cons1387, cons1283)
    rule2668 = ReplacementRule(pattern2668, replacement2668)

    pattern2669 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1388)
    rule2669 = ReplacementRule(pattern2669, replacement2669)

    pattern2670 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons50, cons127, cons210, cons19, cons5, cons1267, cons1388)
    rule2670 = ReplacementRule(pattern2670, replacement2670)

    pattern2671 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1172)
    rule2671 = ReplacementRule(pattern2671, replacement2671)

    pattern2672 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1172)
    rule2672 = ReplacementRule(pattern2672, replacement2672)

    pattern2673 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons96)
    rule2673 = ReplacementRule(pattern2673, replacement2673)

    pattern2674 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons33, cons96)
    rule2674 = ReplacementRule(pattern2674, replacement2674)

    pattern2675 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons145)
    rule2675 = ReplacementRule(pattern2675, replacement2675)

    pattern2676 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons145)
    rule2676 = ReplacementRule(pattern2676, replacement2676)

    pattern2677 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons1308, cons20)
    rule2677 = ReplacementRule(pattern2677, replacement2677)

    pattern2678 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1267, cons1308, cons20)
    rule2678 = ReplacementRule(pattern2678, replacement2678)

    pattern2679 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1308, cons21)
    rule2679 = ReplacementRule(pattern2679, replacement2679)

    pattern2680 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1267, cons1308, cons21)
    rule2680 = ReplacementRule(pattern2680, replacement2680)

    pattern2681 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64, cons1389)
    rule2681 = ReplacementRule(pattern2681, replacement2681)

    pattern2682 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64, cons1389)
    rule2682 = ReplacementRule(pattern2682, replacement2682)

    pattern2683 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1267, cons20)
    rule2683 = ReplacementRule(pattern2683, replacement2683)

    pattern2684 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons4, cons5, cons1267, cons20)
    rule2684 = ReplacementRule(pattern2684, replacement2684)

    pattern2685 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons5, cons1267, cons21)
    rule2685 = ReplacementRule(pattern2685, replacement2685)

    pattern2686 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons5, cons1267, cons21)
    rule2686 = ReplacementRule(pattern2686, replacement2686)

    pattern2687 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons33, cons96, cons1390)
    rule2687 = ReplacementRule(pattern2687, replacement2687)

    pattern2688 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons33, cons96, cons1390)
    rule2688 = ReplacementRule(pattern2688, replacement2688)

    pattern2689 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons33, cons170, cons1391)
    rule2689 = ReplacementRule(pattern2689, replacement2689)

    pattern2690 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons1269, cons33, cons170, cons1391)
    rule2690 = ReplacementRule(pattern2690, replacement2690)

    pattern2691 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392)
    rule2691 = ReplacementRule(pattern2691, replacement2691)

    pattern2692 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392)
    rule2692 = ReplacementRule(pattern2692, replacement2692)

    pattern2693 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1172, cons96, cons91)
    rule2693 = ReplacementRule(pattern2693, replacement2693)

    pattern2694 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1172, cons96, cons91)
    rule2694 = ReplacementRule(pattern2694, replacement2694)

    pattern2695 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons1172, cons96, cons1393, cons1394)
    rule2695 = ReplacementRule(pattern2695, replacement2695)

    pattern2696 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons1172, cons96, cons1393, cons1394)
    rule2696 = ReplacementRule(pattern2696, replacement2696)

    pattern2697 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons1172, cons96, cons1393, cons1395)
    rule2697 = ReplacementRule(pattern2697, replacement2697)

    pattern2698 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons1269, cons1172, cons96, cons1393, cons1395)
    rule2698 = ReplacementRule(pattern2698, replacement2698)

    pattern2699 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons1392, cons1307, cons89, cons91, cons1396)
    rule2699 = ReplacementRule(pattern2699, replacement2699)

    pattern2700 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons1392, cons1307, cons89, cons91, cons1396)
    rule2700 = ReplacementRule(pattern2700, replacement2700)

    pattern2701 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons1392, cons1307, cons89, cons91, cons1395)
    rule2701 = ReplacementRule(pattern2701, replacement2701)

    pattern2702 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons1269, cons1392, cons1307, cons89, cons91, cons1395)
    rule2702 = ReplacementRule(pattern2702, replacement2702)

    pattern2703 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392, cons1307, cons348, cons215, cons1395)
    rule2703 = ReplacementRule(pattern2703, replacement2703)

    pattern2704 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(4), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392, cons1307, cons348, cons215, cons1395)
    rule2704 = ReplacementRule(pattern2704, replacement2704)

    pattern2705 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(6), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1172, cons586, cons1397, cons1398, cons1399, cons145)
    rule2705 = ReplacementRule(pattern2705, replacement2705)

    pattern2706 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(6), x_), cons2, cons3, cons29, cons50, cons127, cons19, cons4, cons1269, cons1172, cons586, cons1397, cons1398, cons1399, cons145)
    rule2706 = ReplacementRule(pattern2706, replacement2706)

    pattern2707 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1400, cons1401)
    rule2707 = ReplacementRule(pattern2707, replacement2707)

    pattern2708 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons29, cons50, cons127, cons1269, cons1400, cons1401)
    rule2708 = ReplacementRule(pattern2708, replacement2708)

    pattern2709 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**n_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons87, cons1402)
    rule2709 = ReplacementRule(pattern2709, replacement2709)

    pattern2710 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**n_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons50, cons127, cons210, cons5, cons1269, cons87, cons1402)
    rule2710 = ReplacementRule(pattern2710, replacement2710)

    pattern2711 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148, cons1404)
    rule2711 = ReplacementRule(pattern2711, replacement2711)

    pattern2712 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148, cons1404)
    rule2712 = ReplacementRule(pattern2712, replacement2712)

    pattern2713 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148, cons1405)
    rule2713 = ReplacementRule(pattern2713, replacement2713)

    pattern2714 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148, cons1405)
    rule2714 = ReplacementRule(pattern2714, replacement2714)

    pattern2715 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148)
    rule2715 = ReplacementRule(pattern2715, replacement2715)

    pattern2716 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons148)
    rule2716 = ReplacementRule(pattern2716, replacement2716)

    pattern2717 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139, cons167)
    rule2717 = ReplacementRule(pattern2717, replacement2717)

    pattern2718 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139, cons167)
    rule2718 = ReplacementRule(pattern2718, replacement2718)

    pattern2719 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139, cons90)
    rule2719 = ReplacementRule(pattern2719, replacement2719)

    pattern2720 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139, cons90)
    rule2720 = ReplacementRule(pattern2720, replacement2720)

    pattern2721 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139)
    rule2721 = ReplacementRule(pattern2721, replacement2721)

    pattern2722 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons139)
    rule2722 = ReplacementRule(pattern2722, replacement2722)

    pattern2723 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2723 = ReplacementRule(pattern2723, replacement2723)

    pattern2724 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons210, cons1269)
    rule2724 = ReplacementRule(pattern2724, replacement2724)

    pattern2725 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(d_*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269)
    rule2725 = ReplacementRule(pattern2725, replacement2725)

    pattern2726 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(d_*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269)
    rule2726 = ReplacementRule(pattern2726, replacement2726)

    pattern2727 = Pattern(Integral(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2727 = ReplacementRule(pattern2727, With2727)

    pattern2728 = Pattern(Integral(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons1269)
    rule2728 = ReplacementRule(pattern2728, With2728)

    pattern2729 = Pattern(Integral(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269)
    rule2729 = ReplacementRule(pattern2729, replacement2729)

    pattern2730 = Pattern(Integral(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269)
    rule2730 = ReplacementRule(pattern2730, replacement2730)

    pattern2731 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons1406, cons90)
    rule2731 = ReplacementRule(pattern2731, replacement2731)

    pattern2732 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons1406, cons90)
    rule2732 = ReplacementRule(pattern2732, replacement2732)

    pattern2733 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons1406, cons465)
    rule2733 = ReplacementRule(pattern2733, replacement2733)

    pattern2734 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1403, cons1406, cons465)
    rule2734 = ReplacementRule(pattern2734, replacement2734)

    pattern2735 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1269, cons20, cons1407)
    rule2735 = ReplacementRule(pattern2735, replacement2735)

    pattern2736 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons4, cons5, cons1269, cons20, cons1407)
    rule2736 = ReplacementRule(pattern2736, replacement2736)

    pattern2737 = Pattern(Integral((WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1408, cons269, cons148, cons1409)
    rule2737 = ReplacementRule(pattern2737, replacement2737)

    pattern2738 = Pattern(Integral((WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_, x_), cons2, cons3, cons29, cons50, cons127, cons210, cons1269, cons1408, cons269, cons148, cons1409)
    rule2738 = ReplacementRule(pattern2738, replacement2738)

    pattern2739 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1267, cons1301, cons1382)
    rule2739 = ReplacementRule(pattern2739, replacement2739)

    pattern2740 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1267, cons1301, cons1382)
    rule2740 = ReplacementRule(pattern2740, replacement2740)

    pattern2741 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1267, cons20, cons13, cons1386)
    rule2741 = ReplacementRule(pattern2741, replacement2741)

    pattern2742 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons1267, cons20, cons13, cons1386)
    rule2742 = ReplacementRule(pattern2742, replacement2742)

    pattern2743 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1267, cons1172)
    rule2743 = ReplacementRule(pattern2743, replacement2743)

    pattern2744 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1267, cons1172)
    rule2744 = ReplacementRule(pattern2744, replacement2744)

    pattern2745 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1267, cons1308, cons20)
    rule2745 = ReplacementRule(pattern2745, replacement2745)

    pattern2746 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1267, cons1308, cons20)
    rule2746 = ReplacementRule(pattern2746, replacement2746)

    pattern2747 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1267, cons1308, cons21)
    rule2747 = ReplacementRule(pattern2747, replacement2747)

    pattern2748 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1267, cons1308, cons21)
    rule2748 = ReplacementRule(pattern2748, replacement2748)

    pattern2749 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64, cons1389)
    rule2749 = ReplacementRule(pattern2749, replacement2749)

    pattern2750 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons1267, cons64, cons1389)
    rule2750 = ReplacementRule(pattern2750, replacement2750)

    pattern2751 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons1267, cons20)
    rule2751 = ReplacementRule(pattern2751, replacement2751)

    pattern2752 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons1267, cons20)
    rule2752 = ReplacementRule(pattern2752, replacement2752)

    pattern2753 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1267, cons21)
    rule2753 = ReplacementRule(pattern2753, replacement2753)

    pattern2754 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons1267, cons21)
    rule2754 = ReplacementRule(pattern2754, replacement2754)

    pattern2755 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392)
    rule2755 = ReplacementRule(pattern2755, replacement2755)

    pattern2756 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1269, cons1392)
    rule2756 = ReplacementRule(pattern2756, replacement2756)

    pattern2757 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1269, cons1410, cons1392)
    rule2757 = ReplacementRule(pattern2757, replacement2757)

    pattern2758 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1269, cons1410, cons1392)
    rule2758 = ReplacementRule(pattern2758, replacement2758)

    pattern2759 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons1172)
    rule2759 = ReplacementRule(pattern2759, replacement2759)

    pattern2760 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons1269, cons1172)
    rule2760 = ReplacementRule(pattern2760, replacement2760)

    pattern2761 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons1269)
    rule2761 = ReplacementRule(pattern2761, replacement2761)

    pattern2762 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons1269)
    rule2762 = ReplacementRule(pattern2762, replacement2762)

    pattern2763 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons149)
    rule2763 = ReplacementRule(pattern2763, replacement2763)

    pattern2764 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons149)
    rule2764 = ReplacementRule(pattern2764, replacement2764)

    pattern2765 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2765 = ReplacementRule(pattern2765, replacement2765)

    pattern2766 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2766 = ReplacementRule(pattern2766, replacement2766)

    pattern2767 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2767 = ReplacementRule(pattern2767, replacement2767)

    pattern2768 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2768 = ReplacementRule(pattern2768, replacement2768)

    pattern2769 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule2769 = ReplacementRule(pattern2769, replacement2769)

    pattern2770 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1267)
    rule2770 = ReplacementRule(pattern2770, replacement2770)

    pattern2771 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1412, cons1413, cons107)
    rule2771 = ReplacementRule(pattern2771, replacement2771)

    pattern2772 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1412, cons1413, cons107)
    rule2772 = ReplacementRule(pattern2772, replacement2772)

    pattern2773 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1324)
    rule2773 = ReplacementRule(pattern2773, replacement2773)

    pattern2774 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1324)
    rule2774 = ReplacementRule(pattern2774, replacement2774)

    pattern2775 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2775 = ReplacementRule(pattern2775, replacement2775)

    pattern2776 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2776 = ReplacementRule(pattern2776, replacement2776)

    pattern2777 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule2777 = ReplacementRule(pattern2777, replacement2777)

    pattern2778 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule2778 = ReplacementRule(pattern2778, replacement2778)

    pattern2779 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2779 = ReplacementRule(pattern2779, replacement2779)

    pattern2780 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2780 = ReplacementRule(pattern2780, replacement2780)

    pattern2781 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2781 = ReplacementRule(pattern2781, replacement2781)

    pattern2782 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2782 = ReplacementRule(pattern2782, replacement2782)

    pattern2783 = Pattern(Integral(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2783 = ReplacementRule(pattern2783, replacement2783)

    pattern2784 = Pattern(Integral(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2784 = ReplacementRule(pattern2784, replacement2784)

    pattern2785 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2785 = ReplacementRule(pattern2785, replacement2785)

    pattern2786 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1411)
    rule2786 = ReplacementRule(pattern2786, replacement2786)

    pattern2787 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2787 = ReplacementRule(pattern2787, replacement2787)

    pattern2788 = Pattern(Integral(S(1)/(sqrt(WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons73, cons1269, cons1325)
    rule2788 = ReplacementRule(pattern2788, replacement2788)

    pattern2789 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule2789 = ReplacementRule(pattern2789, replacement2789)

    pattern2790 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267)
    rule2790 = ReplacementRule(pattern2790, replacement2790)

    pattern2791 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2791 = ReplacementRule(pattern2791, replacement2791)

    pattern2792 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269)
    rule2792 = ReplacementRule(pattern2792, replacement2792)

    pattern2793 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons72)
    rule2793 = ReplacementRule(pattern2793, replacement2793)

    pattern2794 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons72)
    rule2794 = ReplacementRule(pattern2794, replacement2794)

    pattern2795 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1414)
    rule2795 = ReplacementRule(pattern2795, replacement2795)

    pattern2796 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1414)
    rule2796 = ReplacementRule(pattern2796, replacement2796)

    pattern2797 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule2797 = ReplacementRule(pattern2797, replacement2797)

    pattern2798 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1324)
    rule2798 = ReplacementRule(pattern2798, replacement2798)

    pattern2799 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2799 = ReplacementRule(pattern2799, replacement2799)

    pattern2800 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1269, cons1325)
    rule2800 = ReplacementRule(pattern2800, replacement2800)

    pattern2801 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule2801 = ReplacementRule(pattern2801, replacement2801)

    pattern2802 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule2802 = ReplacementRule(pattern2802, replacement2802)

    pattern2803 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1415)
    rule2803 = ReplacementRule(pattern2803, replacement2803)

    pattern2804 = Pattern(Integral(S(1)/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1415)
    rule2804 = ReplacementRule(pattern2804, replacement2804)

    pattern2805 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule2805 = ReplacementRule(pattern2805, replacement2805)

    pattern2806 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1267, cons1324)
    rule2806 = ReplacementRule(pattern2806, replacement2806)

    pattern2807 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1415)
    rule2807 = ReplacementRule(pattern2807, replacement2807)

    pattern2808 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons73, cons1415)
    rule2808 = ReplacementRule(pattern2808, replacement2808)

    pattern2809 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1416, cons87)
    rule2809 = ReplacementRule(pattern2809, replacement2809)

    pattern2810 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons72, cons1267, cons1416, cons87)
    rule2810 = ReplacementRule(pattern2810, replacement2810)

    pattern2811 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1267, cons1325, cons1306)
    rule2811 = ReplacementRule(pattern2811, replacement2811)

    pattern2812 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons1267, cons1325, cons1306)
    rule2812 = ReplacementRule(pattern2812, replacement2812)

    pattern2813 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons73, cons1417, cons1418)
    rule2813 = ReplacementRule(pattern2813, replacement2813)

    pattern2814 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons73, cons1417, cons1418)
    rule2814 = ReplacementRule(pattern2814, replacement2814)

    pattern2815 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons1418)
    rule2815 = ReplacementRule(pattern2815, replacement2815)

    pattern2816 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons1418)
    rule2816 = ReplacementRule(pattern2816, replacement2816)

    pattern2817 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons149, cons20, cons87)
    rule2817 = ReplacementRule(pattern2817, replacement2817)

    pattern2818 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons149, cons20, cons87)
    rule2818 = ReplacementRule(pattern2818, replacement2818)

    pattern2819 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons149, cons1419)
    rule2819 = ReplacementRule(pattern2819, replacement2819)

    pattern2820 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons149, cons1419)
    rule2820 = ReplacementRule(pattern2820, replacement2820)

    pattern2821 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons87)
    rule2821 = ReplacementRule(pattern2821, replacement2821)

    pattern2822 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons87)
    rule2822 = ReplacementRule(pattern2822, replacement2822)

    pattern2823 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20, cons40)
    rule2823 = ReplacementRule(pattern2823, replacement2823)

    pattern2824 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons25, cons20, cons40)
    rule2824 = ReplacementRule(pattern2824, replacement2824)

    pattern2825 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons25, cons20, cons149)
    rule2825 = ReplacementRule(pattern2825, replacement2825)

    pattern2826 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons25, cons20, cons149)
    rule2826 = ReplacementRule(pattern2826, replacement2826)

    pattern2827 = Pattern(Integral((WC('g', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons25, cons21)
    rule2827 = ReplacementRule(pattern2827, replacement2827)

    pattern2828 = Pattern(Integral((WC('g', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons25, cons21)
    rule2828 = ReplacementRule(pattern2828, replacement2828)

    pattern2829 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons149, cons20, cons87)
    rule2829 = ReplacementRule(pattern2829, replacement2829)

    pattern2830 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons5, cons73, cons149, cons20, cons87)
    rule2830 = ReplacementRule(pattern2830, replacement2830)

    pattern2831 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons149, cons1419)
    rule2831 = ReplacementRule(pattern2831, replacement2831)

    pattern2832 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons73, cons149, cons1419)
    rule2832 = ReplacementRule(pattern2832, replacement2832)

    pattern2833 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons20)
    rule2833 = ReplacementRule(pattern2833, replacement2833)

    pattern2834 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons4, cons5, cons20)
    rule2834 = ReplacementRule(pattern2834, replacement2834)

    pattern2835 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(S(1)/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons21, cons87, cons40)
    rule2835 = ReplacementRule(pattern2835, replacement2835)

    pattern2836 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(S(1)/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons21, cons87, cons40)
    rule2836 = ReplacementRule(pattern2836, replacement2836)

    pattern2837 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons21, cons87, cons149)
    rule2837 = ReplacementRule(pattern2837, replacement2837)

    pattern2838 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons5, cons21, cons87, cons149)
    rule2838 = ReplacementRule(pattern2838, replacement2838)

    pattern2839 = Pattern(Integral((WC('g', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons21, cons25)
    rule2839 = ReplacementRule(pattern2839, replacement2839)

    pattern2840 = Pattern(Integral((WC('g', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('p', S(1))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))/cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons19, cons4, cons5, cons21, cons25)
    rule2840 = ReplacementRule(pattern2840, replacement2840)

    pattern2841 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1420, cons1267, cons20, cons87)
    rule2841 = ReplacementRule(pattern2841, replacement2841)

    pattern2842 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons1420, cons1267, cons20, cons87)
    rule2842 = ReplacementRule(pattern2842, replacement2842)

    pattern2843 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons72, cons1267, cons20, cons1311)
    rule2843 = ReplacementRule(pattern2843, replacement2843)

    pattern2844 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons72, cons1267, cons20, cons1311)
    rule2844 = ReplacementRule(pattern2844, replacement2844)

    pattern2845 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73)
    rule2845 = ReplacementRule(pattern2845, replacement2845)

    pattern2846 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73)
    rule2846 = ReplacementRule(pattern2846, replacement2846)

    pattern2847 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons72, cons1267)
    rule2847 = ReplacementRule(pattern2847, replacement2847)

    pattern2848 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons72, cons1267)
    rule2848 = ReplacementRule(pattern2848, replacement2848)

    pattern2849 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1421, cons1316)
    rule2849 = ReplacementRule(pattern2849, replacement2849)

    pattern2850 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1421, cons1316)
    rule2850 = ReplacementRule(pattern2850, replacement2850)

    pattern2851 = Pattern(Integral((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons72, cons1267)
    rule2851 = ReplacementRule(pattern2851, replacement2851)

    pattern2852 = Pattern(Integral((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons72, cons1267)
    rule2852 = ReplacementRule(pattern2852, replacement2852)

    pattern2853 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1422, cons1423)
    rule2853 = ReplacementRule(pattern2853, replacement2853)

    pattern2854 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1422, cons1423)
    rule2854 = ReplacementRule(pattern2854, replacement2854)

    pattern2855 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1323, cons685)
    rule2855 = ReplacementRule(pattern2855, replacement2855)

    pattern2856 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons72, cons1267, cons1323, cons685)
    rule2856 = ReplacementRule(pattern2856, replacement2856)

    pattern2857 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1267, cons1325, cons74, cons1424)
    rule2857 = ReplacementRule(pattern2857, replacement2857)

    pattern2858 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1267, cons1325, cons74, cons1424)
    rule2858 = ReplacementRule(pattern2858, replacement2858)

    pattern2859 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons95, cons1425, cons91, cons517, cons1330)
    rule2859 = ReplacementRule(pattern2859, replacement2859)

    pattern2860 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons95, cons1425, cons91, cons517, cons1330)
    rule2860 = ReplacementRule(pattern2860, replacement2860)

    pattern2861 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons33, cons1425, cons348, cons517, cons1330)
    rule2861 = ReplacementRule(pattern2861, replacement2861)

    pattern2862 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons33, cons1425, cons348, cons517, cons1330)
    rule2862 = ReplacementRule(pattern2862, replacement2862)

    pattern2863 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons95, cons1322, cons90, cons517, cons1330)
    rule2863 = ReplacementRule(pattern2863, replacement2863)

    pattern2864 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons95, cons1322, cons90, cons517, cons1330)
    rule2864 = ReplacementRule(pattern2864, replacement2864)

    pattern2865 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons33, cons1322, cons1329, cons517, cons1330)
    rule2865 = ReplacementRule(pattern2865, replacement2865)

    pattern2866 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons33, cons1322, cons1329, cons517, cons1330)
    rule2866 = ReplacementRule(pattern2866, replacement2866)

    pattern2867 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons1426)
    rule2867 = ReplacementRule(pattern2867, replacement2867)

    pattern2868 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons1426)
    rule2868 = ReplacementRule(pattern2868, replacement2868)

    pattern2869 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons89, cons91)
    rule2869 = ReplacementRule(pattern2869, replacement2869)

    pattern2870 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325, cons89, cons91)
    rule2870 = ReplacementRule(pattern2870, replacement2870)

    pattern2871 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons348)
    rule2871 = ReplacementRule(pattern2871, replacement2871)

    pattern2872 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1267, cons1325, cons348)
    rule2872 = ReplacementRule(pattern2872, replacement2872)

    pattern2873 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325)
    rule2873 = ReplacementRule(pattern2873, replacement2873)

    pattern2874 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325)
    rule2874 = ReplacementRule(pattern2874, replacement2874)

    pattern2875 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons89, cons90, cons1427)
    rule2875 = ReplacementRule(pattern2875, replacement2875)

    pattern2876 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons89, cons90, cons1427)
    rule2876 = ReplacementRule(pattern2876, replacement2876)

    pattern2877 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons89, cons91, cons1427)
    rule2877 = ReplacementRule(pattern2877, replacement2877)

    pattern2878 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons89, cons91, cons1427)
    rule2878 = ReplacementRule(pattern2878, replacement2878)

    pattern2879 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325)
    rule2879 = ReplacementRule(pattern2879, replacement2879)

    pattern2880 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1267, cons1325)
    rule2880 = ReplacementRule(pattern2880, replacement2880)

    pattern2881 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons1316)
    rule2881 = ReplacementRule(pattern2881, replacement2881)

    pattern2882 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1267, cons1325, cons1316)
    rule2882 = ReplacementRule(pattern2882, replacement2882)

    pattern2883 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1267, cons1325, cons1428)
    rule2883 = ReplacementRule(pattern2883, replacement2883)

    pattern2884 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1267, cons1325, cons1428)
    rule2884 = ReplacementRule(pattern2884, replacement2884)

    pattern2885 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons95, cons168, cons91)
    rule2885 = ReplacementRule(pattern2885, replacement2885)

    pattern2886 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons95, cons168, cons91)
    rule2886 = ReplacementRule(pattern2886, replacement2886)

    pattern2887 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1269, cons1325, cons33, cons168, cons1429)
    rule2887 = ReplacementRule(pattern2887, replacement2887)

    pattern2888 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1269, cons1325, cons33, cons168, cons1429)
    rule2888 = ReplacementRule(pattern2888, replacement2888)

    pattern2889 = Pattern(Integral(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325)
    rule2889 = ReplacementRule(pattern2889, replacement2889)

    pattern2890 = Pattern(Integral(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325)
    rule2890 = ReplacementRule(pattern2890, replacement2890)

    pattern2891 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2891 = ReplacementRule(pattern2891, replacement2891)

    pattern2892 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2892 = ReplacementRule(pattern2892, replacement2892)

    pattern2893 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1269)
    rule2893 = ReplacementRule(pattern2893, replacement2893)

    pattern2894 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons1269)
    rule2894 = ReplacementRule(pattern2894, replacement2894)

    pattern2895 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325, cons1430, cons1344)
    rule2895 = ReplacementRule(pattern2895, replacement2895)

    pattern2896 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325, cons1430, cons1344)
    rule2896 = ReplacementRule(pattern2896, replacement2896)

    pattern2897 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325, cons1430, cons1346)
    rule2897 = ReplacementRule(pattern2897, replacement2897)

    pattern2898 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons1325, cons1430, cons1346)
    rule2898 = ReplacementRule(pattern2898, replacement2898)

    pattern2899 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430, cons1347)
    rule2899 = ReplacementRule(pattern2899, replacement2899)

    pattern2900 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430, cons1347)
    rule2900 = ReplacementRule(pattern2900, replacement2900)

    pattern2901 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430, cons1348)
    rule2901 = ReplacementRule(pattern2901, replacement2901)

    pattern2902 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1430, cons1348)
    rule2902 = ReplacementRule(pattern2902, replacement2902)

    pattern2903 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1431)
    rule2903 = ReplacementRule(pattern2903, replacement2903)

    pattern2904 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons1431)
    rule2904 = ReplacementRule(pattern2904, replacement2904)

    pattern2905 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons95, cons96, cons90)
    rule2905 = ReplacementRule(pattern2905, replacement2905)

    pattern2906 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons95, cons96, cons90)
    rule2906 = ReplacementRule(pattern2906, replacement2906)

    pattern2907 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule2907 = ReplacementRule(pattern2907, replacement2907)

    pattern2908 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule2908 = ReplacementRule(pattern2908, replacement2908)

    pattern2909 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2909 = ReplacementRule(pattern2909, replacement2909)

    pattern2910 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2910 = ReplacementRule(pattern2910, replacement2910)

    pattern2911 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1269, cons1325)
    rule2911 = ReplacementRule(pattern2911, replacement2911)

    pattern2912 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_/(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons73, cons1269, cons1325)
    rule2912 = ReplacementRule(pattern2912, replacement2912)

    pattern2913 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons89, cons1432)
    rule2913 = ReplacementRule(pattern2913, replacement2913)

    pattern2914 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325, cons89, cons1432)
    rule2914 = ReplacementRule(pattern2914, replacement2914)

    pattern2915 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons107, cons1413, cons1430)
    rule2915 = ReplacementRule(pattern2915, replacement2915)

    pattern2916 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons107, cons1413, cons1430)
    rule2916 = ReplacementRule(pattern2916, replacement2916)

    pattern2917 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(d_*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons29, cons36, cons37, cons107, cons1413, cons1430)
    rule2917 = ReplacementRule(pattern2917, replacement2917)

    pattern2918 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(d_*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons50, cons127, cons29, cons36, cons37, cons107, cons1413, cons1430)
    rule2918 = ReplacementRule(pattern2918, replacement2918)

    pattern2919 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2919 = ReplacementRule(pattern2919, replacement2919)

    pattern2920 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))/(sqrt(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons73, cons1269, cons1325)
    rule2920 = ReplacementRule(pattern2920, replacement2920)

    pattern2921 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1269, cons1325)
    rule2921 = ReplacementRule(pattern2921, replacement2921)

    pattern2922 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons73, cons1269, cons1325)
    rule2922 = ReplacementRule(pattern2922, replacement2922)

    pattern2923 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons5, cons72, cons1267)
    rule2923 = ReplacementRule(pattern2923, replacement2923)

    pattern2924 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons19, cons4, cons5, cons72, cons1267)
    rule2924 = ReplacementRule(pattern2924, replacement2924)

    pattern2925 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons37, cons38, cons19, cons1433)
    rule2925 = ReplacementRule(pattern2925, replacement2925)

    pattern2926 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons37, cons38, cons19, cons1433)
    rule2926 = ReplacementRule(pattern2926, replacement2926)

    pattern2927 = Pattern(Integral((A_ + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons50, cons127, cons36, cons38, cons1230)
    rule2927 = ReplacementRule(pattern2927, replacement2927)

    pattern2928 = Pattern(Integral((A_ + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**WC('m', S(1)), x_), cons50, cons127, cons36, cons38, cons1230)
    rule2928 = ReplacementRule(pattern2928, replacement2928)

    pattern2929 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(A_ + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1434)
    rule2929 = ReplacementRule(pattern2929, replacement2929)

    pattern2930 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(A_ + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons1434)
    rule2930 = ReplacementRule(pattern2930, replacement2930)

    pattern2931 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(A_ + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons33, cons96)
    rule2931 = ReplacementRule(pattern2931, replacement2931)

    pattern2932 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(A_ + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons33, cons96)
    rule2932 = ReplacementRule(pattern2932, replacement2932)

    pattern2933 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule2933 = ReplacementRule(pattern2933, replacement2933)

    pattern2934 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(A_ + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule2934 = ReplacementRule(pattern2934, replacement2934)

    pattern2935 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule2935 = ReplacementRule(pattern2935, replacement2935)

    pattern2936 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons35)
    rule2936 = ReplacementRule(pattern2936, replacement2936)

    pattern2937 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule2937 = ReplacementRule(pattern2937, replacement2937)

    pattern2938 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1435)
    rule2938 = ReplacementRule(pattern2938, replacement2938)

    pattern2939 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1436, cons79)
    rule2939 = ReplacementRule(pattern2939, replacement2939)

    pattern2940 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons1436, cons79)
    rule2940 = ReplacementRule(pattern2940, replacement2940)

    pattern2941 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1437, cons79)
    rule2941 = ReplacementRule(pattern2941, replacement2941)

    pattern2942 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons1437, cons79)
    rule2942 = ReplacementRule(pattern2942, replacement2942)

    pattern2943 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1267)
    rule2943 = ReplacementRule(pattern2943, replacement2943)

    pattern2944 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1267)
    rule2944 = ReplacementRule(pattern2944, replacement2944)

    pattern2945 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1267)
    rule2945 = ReplacementRule(pattern2945, replacement2945)

    pattern2946 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1267)
    rule2946 = ReplacementRule(pattern2946, replacement2946)

    pattern2947 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1269)
    rule2947 = ReplacementRule(pattern2947, replacement2947)

    pattern2948 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons33, cons96, cons1269)
    rule2948 = ReplacementRule(pattern2948, replacement2948)

    pattern2949 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1269)
    rule2949 = ReplacementRule(pattern2949, replacement2949)

    pattern2950 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons33, cons96, cons1269)
    rule2950 = ReplacementRule(pattern2950, replacement2950)

    pattern2951 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons274)
    rule2951 = ReplacementRule(pattern2951, replacement2951)

    pattern2952 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons274)
    rule2952 = ReplacementRule(pattern2952, replacement2952)

    pattern2953 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule2953 = ReplacementRule(pattern2953, replacement2953)

    pattern2954 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons50, cons127, cons36, cons38, cons19, cons274)
    rule2954 = ReplacementRule(pattern2954, replacement2954)

    pattern2955 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons5, cons21)
    rule2955 = ReplacementRule(pattern2955, replacement2955)

    pattern2956 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons37, cons38, cons19, cons5, cons21)
    rule2956 = ReplacementRule(pattern2956, replacement2956)

    pattern2957 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons5, cons21)
    rule2957 = ReplacementRule(pattern2957, replacement2957)

    pattern2958 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons50, cons127, cons36, cons38, cons19, cons5, cons21)
    rule2958 = ReplacementRule(pattern2958, replacement2958)

    pattern2959 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons33, cons96)
    rule2959 = ReplacementRule(pattern2959, replacement2959)

    pattern2960 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons33, cons96)
    rule2960 = ReplacementRule(pattern2960, replacement2960)

    pattern2961 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons33, cons96)
    rule2961 = ReplacementRule(pattern2961, replacement2961)

    pattern2962 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons33, cons96)
    rule2962 = ReplacementRule(pattern2962, replacement2962)

    pattern2963 = Pattern(Integral((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1269, cons274)
    rule2963 = ReplacementRule(pattern2963, replacement2963)

    pattern2964 = Pattern(Integral((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1269, cons274)
    rule2964 = ReplacementRule(pattern2964, replacement2964)

    pattern2965 = Pattern(Integral((c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1269, cons274)
    rule2965 = ReplacementRule(pattern2965, replacement2965)

    pattern2966 = Pattern(Integral((c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1269, cons274)
    rule2966 = ReplacementRule(pattern2966, replacement2966)

    pattern2967 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons72, cons1267, cons1438)
    rule2967 = ReplacementRule(pattern2967, replacement2967)

    pattern2968 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons72, cons1267, cons1438)
    rule2968 = ReplacementRule(pattern2968, replacement2968)

    pattern2969 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons72, cons1267, cons1438)
    rule2969 = ReplacementRule(pattern2969, replacement2969)

    pattern2970 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons72, cons1267, cons1438)
    rule2970 = ReplacementRule(pattern2970, replacement2970)

    pattern2971 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons72, cons1267, cons1323)
    rule2971 = ReplacementRule(pattern2971, replacement2971)

    pattern2972 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons72, cons1267, cons1323)
    rule2972 = ReplacementRule(pattern2972, replacement2972)

    pattern2973 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))/sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons72, cons1267, cons1323)
    rule2973 = ReplacementRule(pattern2973, replacement2973)

    pattern2974 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))/sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons72, cons1267, cons1323)
    rule2974 = ReplacementRule(pattern2974, replacement2974)

    pattern2975 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons72, cons1267, cons1323, cons216)
    rule2975 = ReplacementRule(pattern2975, replacement2975)

    pattern2976 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons72, cons1267, cons1323, cons216)
    rule2976 = ReplacementRule(pattern2976, replacement2976)

    pattern2977 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons72, cons1267, cons1323, cons216)
    rule2977 = ReplacementRule(pattern2977, replacement2977)

    pattern2978 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons72, cons1267, cons1323, cons216)
    rule2978 = ReplacementRule(pattern2978, replacement2978)

    pattern2979 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1267, cons1325, cons33, cons1322)
    rule2979 = ReplacementRule(pattern2979, replacement2979)

    pattern2980 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1267, cons1325, cons33, cons1322)
    rule2980 = ReplacementRule(pattern2980, replacement2980)

    pattern2981 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1267, cons1325, cons33, cons1322)
    rule2981 = ReplacementRule(pattern2981, replacement2981)

    pattern2982 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1267, cons1325, cons33, cons1322)
    rule2982 = ReplacementRule(pattern2982, replacement2982)

    pattern2983 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1267, cons1325, cons1323, cons1439)
    rule2983 = ReplacementRule(pattern2983, replacement2983)

    pattern2984 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons73, cons1267, cons1325, cons1323, cons1439)
    rule2984 = ReplacementRule(pattern2984, replacement2984)

    pattern2985 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1267, cons1325, cons1323, cons1439)
    rule2985 = ReplacementRule(pattern2985, replacement2985)

    pattern2986 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons73, cons1267, cons1325, cons1323, cons1439)
    rule2986 = ReplacementRule(pattern2986, replacement2986)

    pattern2987 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1267, cons1325, cons1323, cons216)
    rule2987 = ReplacementRule(pattern2987, replacement2987)

    pattern2988 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1267, cons1325, cons1323, cons216)
    rule2988 = ReplacementRule(pattern2988, replacement2988)

    pattern2989 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1267, cons1325, cons1323, cons216)
    rule2989 = ReplacementRule(pattern2989, replacement2989)

    pattern2990 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1267, cons1325, cons1323, cons216)
    rule2990 = ReplacementRule(pattern2990, replacement2990)

    pattern2991 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325, cons95, cons170, cons91)
    rule2991 = ReplacementRule(pattern2991, replacement2991)

    pattern2992 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325, cons95, cons170, cons91)
    rule2992 = ReplacementRule(pattern2992, replacement2992)

    pattern2993 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325, cons95, cons170, cons91)
    rule2993 = ReplacementRule(pattern2993, replacement2993)

    pattern2994 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325, cons95, cons170, cons91)
    rule2994 = ReplacementRule(pattern2994, replacement2994)

    pattern2995 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1269, cons1325, cons33, cons170, cons1440)
    rule2995 = ReplacementRule(pattern2995, replacement2995)

    pattern2996 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1269, cons1325, cons33, cons170, cons1440)
    rule2996 = ReplacementRule(pattern2996, replacement2996)

    pattern2997 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1269, cons1325, cons33, cons170, cons1440)
    rule2997 = ReplacementRule(pattern2997, replacement2997)

    pattern2998 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('m', S(1))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1269, cons1325, cons33, cons170, cons1440)
    rule2998 = ReplacementRule(pattern2998, replacement2998)

    pattern2999 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule2999 = ReplacementRule(pattern2999, replacement2999)

    pattern3000 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons37, cons38, cons1269)
    rule3000 = ReplacementRule(pattern3000, replacement3000)

    pattern3001 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule3001 = ReplacementRule(pattern3001, replacement3001)

    pattern3002 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)), x_), cons2, cons3, cons29, cons50, cons127, cons36, cons38, cons1269)
    rule3002 = ReplacementRule(pattern3002, replacement3002)

    pattern3003 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3003 = ReplacementRule(pattern3003, replacement3003)

    pattern3004 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3004 = ReplacementRule(pattern3004, replacement3004)

    pattern3005 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3005 = ReplacementRule(pattern3005, replacement3005)

    pattern3006 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**(S(3)/2)*sqrt(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3006 = ReplacementRule(pattern3006, replacement3006)

    pattern3007 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule3007 = ReplacementRule(pattern3007, replacement3007)

    pattern3008 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule3008 = ReplacementRule(pattern3008, replacement3008)

    pattern3009 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule3009 = ReplacementRule(pattern3009, replacement3009)

    pattern3010 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons4, cons73, cons1269, cons1325, cons33, cons96, cons1339)
    rule3010 = ReplacementRule(pattern3010, replacement3010)

    pattern3011 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3011 = ReplacementRule(pattern3011, replacement3011)

    pattern3012 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3012 = ReplacementRule(pattern3012, replacement3012)

    pattern3013 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3013 = ReplacementRule(pattern3013, replacement3013)

    pattern3014 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3014 = ReplacementRule(pattern3014, replacement3014)

    pattern3015 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3015 = ReplacementRule(pattern3015, replacement3015)

    pattern3016 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3016 = ReplacementRule(pattern3016, replacement3016)

    pattern3017 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3017 = ReplacementRule(pattern3017, replacement3017)

    pattern3018 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3018 = ReplacementRule(pattern3018, replacement3018)

    pattern3019 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3019 = ReplacementRule(pattern3019, replacement3019)

    pattern3020 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons73, cons1269, cons1325)
    rule3020 = ReplacementRule(pattern3020, replacement3020)

    pattern3021 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(c_ + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3021 = ReplacementRule(pattern3021, replacement3021)

    pattern3022 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))/(sqrt(c_ + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))*sqrt(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons73, cons1269, cons1325)
    rule3022 = ReplacementRule(pattern3022, replacement3022)

    pattern3023 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1269, cons1325)
    rule3023 = ReplacementRule(pattern3023, replacement3023)

    pattern3024 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons73, cons1269, cons1325)
    rule3024 = ReplacementRule(pattern3024, replacement3024)

    pattern3025 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1269, cons1325)
    rule3025 = ReplacementRule(pattern3025, replacement3025)

    pattern3026 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons73, cons1269, cons1325)
    rule3026 = ReplacementRule(pattern3026, replacement3026)

    pattern3027 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons21)
    rule3027 = ReplacementRule(pattern3027, replacement3027)

    pattern3028 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2)), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons37, cons38, cons19, cons4, cons5, cons21)
    rule3028 = ReplacementRule(pattern3028, replacement3028)

    pattern3029 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons5, cons21)
    rule3029 = ReplacementRule(pattern3029, replacement3029)

    pattern3030 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**p_)**m_*(WC('A', S(0)) + WC('C', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))**S(2))*(WC('c', S(0)) + WC('d', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons36, cons38, cons19, cons4, cons5, cons21)
    rule3030 = ReplacementRule(pattern3030, replacement3030)

    pattern3031 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1441)
    rule3031 = ReplacementRule(pattern3031, replacement3031)

    pattern3032 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons523)
    rule3032 = ReplacementRule(pattern3032, replacement3032)

    pattern3033 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons1443, cons89, cons167)
    rule3033 = ReplacementRule(pattern3033, replacement3033)

    pattern3034 = Pattern(Integral(S(1)/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3034 = ReplacementRule(pattern3034, replacement3034)

    pattern3035 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**(S(-2)), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3035 = ReplacementRule(pattern3035, replacement3035)

    pattern3036 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons91, cons1444)
    rule3036 = ReplacementRule(pattern3036, replacement3036)

    pattern3037 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1445, cons1446)
    rule3037 = ReplacementRule(pattern3037, replacement3037)

    pattern3038 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons4, cons1445, cons1447)
    rule3038 = ReplacementRule(pattern3038, replacement3038)

    pattern3039 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1257, cons1441, cons89, cons167)
    rule3039 = ReplacementRule(pattern3039, replacement3039)

    pattern3040 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1257, cons1441, cons89, cons167)
    rule3040 = ReplacementRule(pattern3040, replacement3040)

    pattern3041 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1257, cons1441, cons89, cons465)
    rule3041 = ReplacementRule(pattern3041, replacement3041)

    pattern3042 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1257, cons1441, cons89, cons465)
    rule3042 = ReplacementRule(pattern3042, replacement3042)

    pattern3043 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1257, cons1441, cons25)
    rule3043 = ReplacementRule(pattern3043, replacement3043)

    pattern3044 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons1257, cons1441, cons25)
    rule3044 = ReplacementRule(pattern3044, replacement3044)

    pattern3045 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1257, cons87, cons1442)
    rule3045 = ReplacementRule(pattern3045, replacement3045)

    pattern3046 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1257, cons87, cons1442)
    rule3046 = ReplacementRule(pattern3046, replacement3046)

    pattern3047 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons87, cons1448, cons1154, cons1449)
    rule3047 = ReplacementRule(pattern3047, replacement3047)

    pattern3048 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons87, cons1448, cons1154, cons1449)
    rule3048 = ReplacementRule(pattern3048, replacement3048)

    pattern3049 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons20, cons150)
    rule3049 = ReplacementRule(pattern3049, replacement3049)

    pattern3050 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons20, cons150)
    rule3050 = ReplacementRule(pattern3050, replacement3050)

    pattern3051 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons1441, cons198)
    rule3051 = ReplacementRule(pattern3051, replacement3051)

    pattern3052 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons1441, cons198)
    rule3052 = ReplacementRule(pattern3052, replacement3052)

    pattern3053 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_/sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons91)
    rule3053 = ReplacementRule(pattern3053, replacement3053)

    pattern3054 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_/cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons1442, cons89, cons91)
    rule3054 = ReplacementRule(pattern3054, replacement3054)

    pattern3055 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1442, cons95, cons167, cons96)
    rule3055 = ReplacementRule(pattern3055, replacement3055)

    pattern3056 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1442, cons95, cons167, cons96)
    rule3056 = ReplacementRule(pattern3056, replacement3056)

    pattern3057 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3057 = ReplacementRule(pattern3057, replacement3057)

    pattern3058 = Pattern(Integral(cos(x_*WC('d', S(1)) + WC('c', S(0)))/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3058 = ReplacementRule(pattern3058, replacement3058)

    pattern3059 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442, cons33, cons168)
    rule3059 = ReplacementRule(pattern3059, replacement3059)

    pattern3060 = Pattern(Integral(cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442, cons33, cons168)
    rule3060 = ReplacementRule(pattern3060, replacement3060)

    pattern3061 = Pattern(Integral(S(1)/((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3061 = ReplacementRule(pattern3061, replacement3061)

    pattern3062 = Pattern(Integral(S(1)/((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))*cos(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442)
    rule3062 = ReplacementRule(pattern3062, replacement3062)

    pattern3063 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442, cons33, cons96)
    rule3063 = ReplacementRule(pattern3063, replacement3063)

    pattern3064 = Pattern(Integral(cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442, cons33, cons96)
    rule3064 = ReplacementRule(pattern3064, replacement3064)

    pattern3065 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1442, cons95, cons91, cons96)
    rule3065 = ReplacementRule(pattern3065, replacement3065)

    pattern3066 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1442, cons95, cons91, cons96)
    rule3066 = ReplacementRule(pattern3066, replacement3066)

    pattern3067 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons130)
    rule3067 = ReplacementRule(pattern3067, replacement3067)

    pattern3068 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons1441, cons65)
    rule3068 = ReplacementRule(pattern3068, replacement3068)

    pattern3069 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons1442, cons152, cons170, cons90)
    rule3069 = ReplacementRule(pattern3069, replacement3069)

    pattern3070 = Pattern(Integral(sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons152)
    rule3070 = ReplacementRule(pattern3070, replacement3070)

    pattern3071 = Pattern(Integral((WC('a', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))) + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons1442, cons377, cons170, cons90, cons324)
    rule3071 = ReplacementRule(pattern3071, replacement3071)

    pattern3072 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1450)
    rule3072 = ReplacementRule(pattern3072, replacement3072)

    pattern3073 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1450, cons89, cons90)
    rule3073 = ReplacementRule(pattern3073, replacement3073)

    pattern3074 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1450)
    rule3074 = ReplacementRule(pattern3074, replacement3074)

    pattern3075 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1450)
    rule3075 = ReplacementRule(pattern3075, replacement3075)

    pattern3076 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1450, cons89, cons91)
    rule3076 = ReplacementRule(pattern3076, replacement3076)

    pattern3077 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1451)
    rule3077 = ReplacementRule(pattern3077, replacement3077)

    pattern3078 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1452, cons1453)
    rule3078 = ReplacementRule(pattern3078, replacement3078)

    pattern3079 = Pattern(Integral(sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1454, cons1452, cons1455)
    rule3079 = ReplacementRule(pattern3079, replacement3079)

    pattern3080 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1454, cons89, cons167)
    rule3080 = ReplacementRule(pattern3080, replacement3080)

    pattern3081 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1456)
    rule3081 = ReplacementRule(pattern3081, With3081)

    pattern3082 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons78)
    rule3082 = ReplacementRule(pattern3082, With3082)

    pattern3083 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1457, cons1458)
    rule3083 = ReplacementRule(pattern3083, With3083)

    pattern3084 = Pattern(Integral(S(1)/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1454)
    rule3084 = ReplacementRule(pattern3084, With3084)

    pattern3085 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1451)
    rule3085 = ReplacementRule(pattern3085, replacement3085)

    pattern3086 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1452, cons1453)
    rule3086 = ReplacementRule(pattern3086, replacement3086)

    pattern3087 = Pattern(Integral(S(1)/sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1454, cons1452, cons1455)
    rule3087 = ReplacementRule(pattern3087, replacement3087)

    pattern3088 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**(S(-3)/2), x_), cons2, cons3, cons8, cons29, cons50, cons1454)
    rule3088 = ReplacementRule(pattern3088, replacement3088)

    pattern3089 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons1454, cons89, cons91, cons1459)
    rule3089 = ReplacementRule(pattern3089, replacement3089)

    pattern3090 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1451)
    rule3090 = ReplacementRule(pattern3090, replacement3090)

    pattern3091 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1451)
    rule3091 = ReplacementRule(pattern3091, replacement3091)

    pattern3092 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1451)
    rule3092 = ReplacementRule(pattern3092, replacement3092)

    pattern3093 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1452, cons1460)
    rule3093 = ReplacementRule(pattern3093, replacement3093)

    pattern3094 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1452, cons1461)
    rule3094 = ReplacementRule(pattern3094, replacement3094)

    pattern3095 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1452, cons1462)
    rule3095 = ReplacementRule(pattern3095, replacement3095)

    pattern3096 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1452, cons1463)
    rule3096 = ReplacementRule(pattern3096, replacement3096)

    pattern3097 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1452, cons1464)
    rule3097 = ReplacementRule(pattern3097, replacement3097)

    pattern3098 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1452, cons1465)
    rule3098 = ReplacementRule(pattern3098, replacement3098)

    pattern3099 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons4, cons586, cons1450, cons1466)
    rule3099 = ReplacementRule(pattern3099, replacement3099)

    pattern3100 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons4, cons586, cons1450, cons1467)
    rule3100 = ReplacementRule(pattern3100, replacement3100)

    pattern3101 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons4, cons586, cons1450, cons1468)
    rule3101 = ReplacementRule(pattern3101, replacement3101)

    pattern3102 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons4, cons586, cons1450, cons1469)
    rule3102 = ReplacementRule(pattern3102, replacement3102)

    pattern3103 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons4, cons586, cons1450, cons1470)
    rule3103 = ReplacementRule(pattern3103, replacement3103)

    pattern3104 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons4, cons586, cons1450, cons1471)
    rule3104 = ReplacementRule(pattern3104, replacement3104)

    pattern3105 = Pattern(Integral((WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons3, cons8, cons29, cons50, cons37, cons38, cons586, cons1452, cons1472)
    rule3105 = ReplacementRule(pattern3105, replacement3105)

    pattern3106 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*(WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons89, cons90, cons1454)
    rule3106 = ReplacementRule(pattern3106, replacement3106)

    pattern3107 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons89, cons90, cons1454)
    rule3107 = ReplacementRule(pattern3107, replacement3107)

    pattern3108 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons89, cons90, cons1454)
    rule3108 = ReplacementRule(pattern3108, replacement3108)

    pattern3109 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/sqrt(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1473, cons1247)
    rule3109 = ReplacementRule(pattern3109, replacement3109)

    pattern3110 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1454, cons1474)
    rule3110 = ReplacementRule(pattern3110, replacement3110)

    pattern3111 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1454, cons1475)
    rule3111 = ReplacementRule(pattern3111, replacement3111)

    pattern3112 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1454, cons1476)
    rule3112 = ReplacementRule(pattern3112, replacement3112)

    pattern3113 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons1454, cons1477)
    rule3113 = ReplacementRule(pattern3113, replacement3113)

    pattern3114 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons1454, cons1478)
    rule3114 = ReplacementRule(pattern3114, replacement3114)

    pattern3115 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons1454, cons1479)
    rule3115 = ReplacementRule(pattern3115, replacement3115)

    pattern3116 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons38, cons89, cons91, cons1454, cons1444)
    rule3116 = ReplacementRule(pattern3116, replacement3116)

    pattern3117 = Pattern(Integral((WC('A', S(0)) + WC('C', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons38, cons89, cons91, cons1454, cons1444)
    rule3117 = ReplacementRule(pattern3117, replacement3117)

    pattern3118 = Pattern(Integral((WC('A', S(0)) + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons89, cons91, cons1454, cons1444)
    rule3118 = ReplacementRule(pattern3118, replacement3118)

    pattern3119 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule3119 = ReplacementRule(pattern3119, replacement3119)

    pattern3120 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons1045)
    rule3120 = ReplacementRule(pattern3120, replacement3120)

    pattern3121 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons87)
    rule3121 = ReplacementRule(pattern3121, replacement3121)

    pattern3122 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons87)
    rule3122 = ReplacementRule(pattern3122, replacement3122)

    pattern3123 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**n_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons25)
    rule3123 = ReplacementRule(pattern3123, replacement3123)

    pattern3124 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**n_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons25)
    rule3124 = ReplacementRule(pattern3124, replacement3124)

    pattern3125 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1257, cons87)
    rule3125 = ReplacementRule(pattern3125, replacement3125)

    pattern3126 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1257, cons87)
    rule3126 = ReplacementRule(pattern3126, replacement3126)

    pattern3127 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_*(S(1)/cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1257, cons25)
    rule3127 = ReplacementRule(pattern3127, replacement3127)

    pattern3128 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))/sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_*(S(1)/sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons1257, cons25)
    rule3128 = ReplacementRule(pattern3128, replacement3128)

    pattern3129 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons3, cons8, cons29, cons13, cons148)
    rule3129 = ReplacementRule(pattern3129, replacement3129)

    pattern3130 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons3, cons8, cons29, cons13, cons148)
    rule3130 = ReplacementRule(pattern3130, replacement3130)

    pattern3131 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons3, cons8, cons29, cons13, cons139)
    rule3131 = ReplacementRule(pattern3131, replacement3131)

    pattern3132 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons3, cons8, cons29, cons13, cons139)
    rule3132 = ReplacementRule(pattern3132, replacement3132)

    pattern3133 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons1456, cons40)
    rule3133 = ReplacementRule(pattern3133, replacement3133)

    pattern3134 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons1456, cons40)
    rule3134 = ReplacementRule(pattern3134, replacement3134)

    pattern3135 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456, cons149)
    rule3135 = ReplacementRule(pattern3135, replacement3135)

    pattern3136 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456, cons149)
    rule3136 = ReplacementRule(pattern3136, replacement3136)

    pattern3137 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons40)
    rule3137 = ReplacementRule(pattern3137, With3137)

    pattern3138 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons40)
    rule3138 = ReplacementRule(pattern3138, With3138)

    pattern3139 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1480, cons149)
    rule3139 = ReplacementRule(pattern3139, replacement3139)

    pattern3140 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1480, cons149)
    rule3140 = ReplacementRule(pattern3140, replacement3140)

    pattern3141 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_, x_), cons2, cons3, cons8, cons29, cons87, cons130)
    rule3141 = ReplacementRule(pattern3141, replacement3141)

    pattern3142 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_, x_), cons2, cons3, cons8, cons29, cons87, cons130)
    rule3142 = ReplacementRule(pattern3142, replacement3142)

    pattern3143 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_, x_), cons2, cons3, cons8, cons29, cons1481, cons40, cons139)
    rule3143 = ReplacementRule(pattern3143, With3143)

    pattern3144 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_, x_), cons2, cons3, cons8, cons29, cons1481, cons40, cons139)
    rule3144 = ReplacementRule(pattern3144, With3144)

    pattern3145 = Pattern(Integral(u_*(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons1456, cons40)
    rule3145 = ReplacementRule(pattern3145, replacement3145)

    pattern3146 = Pattern(Integral(u_*(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons1456, cons40)
    rule3146 = ReplacementRule(pattern3146, replacement3146)

    pattern3147 = Pattern(Integral(u_*(a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456, cons149)
    rule3147 = ReplacementRule(pattern3147, replacement3147)

    pattern3148 = Pattern(Integral(u_*(a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**S(2))**p_, x_), cons2, cons3, cons8, cons29, cons5, cons1456, cons149)
    rule3148 = ReplacementRule(pattern3148, replacement3148)

    pattern3149 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1481, cons1482, cons40)
    rule3149 = ReplacementRule(pattern3149, With3149)

    pattern3150 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1481, cons1482, cons40)
    rule3150 = ReplacementRule(pattern3150, With3150)

    pattern3151 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1481, cons1483)
    rule3151 = ReplacementRule(pattern3151, With3151)

    pattern3152 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons5, cons1481, cons1483)
    rule3152 = ReplacementRule(pattern3152, With3152)

    pattern3153 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons377)
    rule3153 = ReplacementRule(pattern3153, replacement3153)

    pattern3154 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons377)
    rule3154 = ReplacementRule(pattern3154, replacement3154)

    pattern3155 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1481, cons40)
    rule3155 = ReplacementRule(pattern3155, With3155)

    pattern3156 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1481, cons40)
    rule3156 = ReplacementRule(pattern3156, With3156)

    pattern3157 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1484, cons40, cons170)
    rule3157 = ReplacementRule(pattern3157, replacement3157)

    pattern3158 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1484, cons40, cons170)
    rule3158 = ReplacementRule(pattern3158, replacement3158)

    pattern3159 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*cos(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1484, cons40, cons269, cons139)
    rule3159 = ReplacementRule(pattern3159, replacement3159)

    pattern3160 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**p_*sin(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons1482, cons1484, cons40, cons269, cons139)
    rule3160 = ReplacementRule(pattern3160, replacement3160)

    pattern3161 = Pattern(Integral((a_ + (WC('e', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1483)
    rule3161 = ReplacementRule(pattern3161, With3161)

    pattern3162 = Pattern(Integral((a_ + (WC('e', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons1483)
    rule3162 = ReplacementRule(pattern3162, With3162)

    pattern3163 = Pattern(Integral((a_ + (WC('e', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1483, cons248)
    rule3163 = ReplacementRule(pattern3163, With3163)

    pattern3164 = Pattern(Integral((a_ + (WC('e', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))))**n_*WC('b', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons1483, cons248)
    rule3164 = ReplacementRule(pattern3164, With3164)

    pattern3165 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons19, cons1485, cons1481, cons40)
    rule3165 = ReplacementRule(pattern3165, With3165)

    pattern3166 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons19, cons1485, cons1481, cons40)
    rule3166 = ReplacementRule(pattern3166, With3166)

    pattern3167 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*tan(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons1485, cons1481, cons149)
    rule3167 = ReplacementRule(pattern3167, With3167)

    pattern3168 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0)))**n_)**WC('p', S(1))*(S(1)/tan(x_*WC('d', S(1)) + WC('c', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons19, cons5, cons1485, cons1481, cons149)
    rule3168 = ReplacementRule(pattern3168, With3168)

    pattern3169 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**p_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**q_)**n_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons1482, cons1486, cons1487, cons87, cons1488)
    rule3169 = ReplacementRule(pattern3169, With3169)

    pattern3170 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**p_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**q_)**n_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons1482, cons1486, cons1487, cons87, cons1488)
    rule3170 = ReplacementRule(pattern3170, With3170)

    pattern3171 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**p_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**q_)**n_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons1482, cons1486, cons1487, cons87, cons1489)
    rule3171 = ReplacementRule(pattern3171, With3171)

    pattern3172 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**p_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**q_)**n_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons1482, cons1486, cons1487, cons87, cons1489)
    rule3172 = ReplacementRule(pattern3172, With3172)

    pattern3173 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule3173 = ReplacementRule(pattern3173, replacement3173)

    pattern3174 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons47, cons40)
    rule3174 = ReplacementRule(pattern3174, replacement3174)

    pattern3175 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons47, cons149)
    rule3175 = ReplacementRule(pattern3175, replacement3175)

    pattern3176 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons4, cons5, cons48, cons47, cons149)
    rule3176 = ReplacementRule(pattern3176, replacement3176)

    pattern3177 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule3177 = ReplacementRule(pattern3177, With3177)

    pattern3178 = Pattern(Integral(S(1)/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1))), x_), cons2, cons3, cons8, cons29, cons50, cons4, cons48, cons228)
    rule3178 = ReplacementRule(pattern3178, With3178)

    pattern3179 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3179 = ReplacementRule(pattern3179, replacement3179)

    pattern3180 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons47, cons40)
    rule3180 = ReplacementRule(pattern3180, replacement3180)

    pattern3181 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3181 = ReplacementRule(pattern3181, replacement3181)

    pattern3182 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons47, cons149)
    rule3182 = ReplacementRule(pattern3182, replacement3182)

    pattern3183 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons1481, cons40)
    rule3183 = ReplacementRule(pattern3183, With3183)

    pattern3184 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons1481, cons40)
    rule3184 = ReplacementRule(pattern3184, With3184)

    pattern3185 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons377)
    rule3185 = ReplacementRule(pattern3185, replacement3185)

    pattern3186 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons228, cons377)
    rule3186 = ReplacementRule(pattern3186, replacement3186)

    pattern3187 = Pattern(Integral(((WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1483)
    rule3187 = ReplacementRule(pattern3187, With3187)

    pattern3188 = Pattern(Integral(((WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n', S(1))*WC('b', S(1)) + (WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons48, cons1483)
    rule3188 = ReplacementRule(pattern3188, With3188)

    pattern3189 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3189 = ReplacementRule(pattern3189, replacement3189)

    pattern3190 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3190 = ReplacementRule(pattern3190, replacement3190)

    pattern3191 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3191 = ReplacementRule(pattern3191, replacement3191)

    pattern3192 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*sin(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3192 = ReplacementRule(pattern3192, replacement3192)

    pattern3193 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons1481, cons40)
    rule3193 = ReplacementRule(pattern3193, With3193)

    pattern3194 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons1481, cons40)
    rule3194 = ReplacementRule(pattern3194, With3194)

    pattern3195 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3195 = ReplacementRule(pattern3195, replacement3195)

    pattern3196 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3196 = ReplacementRule(pattern3196, replacement3196)

    pattern3197 = Pattern(Integral((a_ + (WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1483, cons248)
    rule3197 = ReplacementRule(pattern3197, With3197)

    pattern3198 = Pattern(Integral((a_ + (WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1483, cons248)
    rule3198 = ReplacementRule(pattern3198, With3198)

    pattern3199 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3199 = ReplacementRule(pattern3199, replacement3199)

    pattern3200 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3200 = ReplacementRule(pattern3200, replacement3200)

    pattern3201 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3201 = ReplacementRule(pattern3201, replacement3201)

    pattern3202 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3202 = ReplacementRule(pattern3202, replacement3202)

    pattern3203 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons48, cons1485, cons228, cons1481, cons40)
    rule3203 = ReplacementRule(pattern3203, With3203)

    pattern3204 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons48, cons1485, cons228, cons1481, cons40)
    rule3204 = ReplacementRule(pattern3204, With3204)

    pattern3205 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3205 = ReplacementRule(pattern3205, replacement3205)

    pattern3206 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3206 = ReplacementRule(pattern3206, replacement3206)

    pattern3207 = Pattern(Integral((a_ + (WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (WC('f', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1483, cons248)
    rule3207 = ReplacementRule(pattern3207, With3207)

    pattern3208 = Pattern(Integral((a_ + (WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**n_*WC('b', S(1)) + (WC('f', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))**WC('n2', S(1))*WC('c', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons1483, cons248)
    rule3208 = ReplacementRule(pattern3208, With3208)

    pattern3209 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3209 = ReplacementRule(pattern3209, replacement3209)

    pattern3210 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons48, cons1485, cons47, cons40)
    rule3210 = ReplacementRule(pattern3210, replacement3210)

    pattern3211 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3211 = ReplacementRule(pattern3211, replacement3211)

    pattern3212 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**p_*tan(x_*WC('e', S(1)) + WC('d', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons48, cons1485, cons47, cons149)
    rule3212 = ReplacementRule(pattern3212, replacement3212)

    pattern3213 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons48, cons1481, cons40)
    rule3213 = ReplacementRule(pattern3213, With3213)

    pattern3214 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n_ + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**n2_)**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons48, cons1481, cons40)
    rule3214 = ReplacementRule(pattern3214, With3214)

    pattern3215 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*(S(1)/tan(x_*WC('e', S(1)) + WC('d', S(0))))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3215 = ReplacementRule(pattern3215, replacement3215)

    pattern3216 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n', S(1)) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**WC('n2', S(1)))**WC('p', S(1))*tan(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons48, cons1482, cons228, cons378)
    rule3216 = ReplacementRule(pattern3216, replacement3216)

    pattern3217 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule3217 = ReplacementRule(pattern3217, replacement3217)

    pattern3218 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons87)
    rule3218 = ReplacementRule(pattern3218, replacement3218)

    pattern3219 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule3219 = ReplacementRule(pattern3219, replacement3219)

    pattern3220 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(a_ + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons47, cons25)
    rule3220 = ReplacementRule(pattern3220, replacement3220)

    pattern3221 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule3221 = ReplacementRule(pattern3221, With3221)

    pattern3222 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))/(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2)), x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228)
    rule3222 = ReplacementRule(pattern3222, With3222)

    pattern3223 = Pattern(Integral((A_ + WC('B', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*sin(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule3223 = ReplacementRule(pattern3223, replacement3223)

    pattern3224 = Pattern(Integral((A_ + WC('B', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0))) + WC('c', S(1))*cos(x_*WC('e', S(1)) + WC('d', S(0)))**S(2))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons36, cons37, cons228, cons87)
    rule3224 = ReplacementRule(pattern3224, replacement3224)

    pattern3225 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule3225 = ReplacementRule(pattern3225, replacement3225)

    pattern3226 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons170)
    rule3226 = ReplacementRule(pattern3226, replacement3226)

    pattern3227 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons96)
    rule3227 = ReplacementRule(pattern3227, replacement3227)

    pattern3228 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons33, cons96)
    rule3228 = ReplacementRule(pattern3228, replacement3228)

    pattern3229 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule3229 = ReplacementRule(pattern3229, replacement3229)

    pattern3230 = Pattern(Integral(cos(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule3230 = ReplacementRule(pattern3230, replacement3230)

    pattern3231 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule3231 = ReplacementRule(pattern3231, replacement3231)

    pattern3232 = Pattern(Integral(cos(x_*WC('f', S(1)) + WC('e', S(0)))/(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule3232 = ReplacementRule(pattern3232, replacement3232)

    pattern3233 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule3233 = ReplacementRule(pattern3233, replacement3233)

    pattern3234 = Pattern(Integral(cos(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons1118)
    rule3234 = ReplacementRule(pattern3234, replacement3234)

    pattern3235 = Pattern(Integral(sin(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule3235 = ReplacementRule(pattern3235, replacement3235)

    pattern3236 = Pattern(Integral(cos(x_*WC('f', S(1)) + WC('e', S(0)))/sqrt(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons127, cons178)
    rule3236 = ReplacementRule(pattern3236, replacement3236)

    pattern3237 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons19, cons1490)
    rule3237 = ReplacementRule(pattern3237, replacement3237)

    pattern3238 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons8, cons29, cons50, cons127, cons19, cons1490)
    rule3238 = ReplacementRule(pattern3238, replacement3238)

    pattern3239 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167)
    rule3239 = ReplacementRule(pattern3239, replacement3239)

    pattern3240 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons167)
    rule3240 = ReplacementRule(pattern3240, replacement3240)

    pattern3241 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons168)
    rule3241 = ReplacementRule(pattern3241, replacement3241)

    pattern3242 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons168)
    rule3242 = ReplacementRule(pattern3242, replacement3242)

    pattern3243 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons1491)
    rule3243 = ReplacementRule(pattern3243, replacement3243)

    pattern3244 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons1491)
    rule3244 = ReplacementRule(pattern3244, replacement3244)

    pattern3245 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*sin(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons33, cons247)
    rule3245 = ReplacementRule(pattern3245, replacement3245)

    pattern3246 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**m_*cos(x_*WC('f', S(1)) + WC('e', S(0)))**n_, x_), cons8, cons29, cons50, cons127, cons19, cons87, cons167, cons33, cons247)
    rule3246 = ReplacementRule(pattern3246, replacement3246)

    pattern3247 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons249)
    rule3247 = ReplacementRule(pattern3247, replacement3247)

    pattern3248 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**m_, x_), cons3, cons8, cons29, cons50, cons127, cons95, cons167, cons249)
    rule3248 = ReplacementRule(pattern3248, replacement3248)

    pattern3249 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91, cons1444)
    rule3249 = ReplacementRule(pattern3249, replacement3249)

    pattern3250 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons3, cons8, cons29, cons50, cons127, cons89, cons91, cons1444)
    rule3250 = ReplacementRule(pattern3250, replacement3250)

    pattern3251 = Pattern(Integral((WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons1444, cons168)
    rule3251 = ReplacementRule(pattern3251, replacement3251)

    pattern3252 = Pattern(Integral((WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons3, cons8, cons29, cons50, cons127, cons95, cons91, cons1444, cons168)
    rule3252 = ReplacementRule(pattern3252, replacement3252)

    pattern3253 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons150, cons1492)
    rule3253 = ReplacementRule(pattern3253, replacement3253)

    pattern3254 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons150, cons1492)
    rule3254 = ReplacementRule(pattern3254, replacement3254)

    pattern3255 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons87)
    rule3255 = ReplacementRule(pattern3255, replacement3255)

    pattern3256 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1267, cons810, cons1493)
    rule3256 = ReplacementRule(pattern3256, replacement3256)

    pattern3257 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1494, cons87)
    rule3257 = ReplacementRule(pattern3257, replacement3257)

    pattern3258 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1456, cons87)
    rule3258 = ReplacementRule(pattern3258, replacement3258)

    pattern3259 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1494, cons810, cons1493)
    rule3259 = ReplacementRule(pattern3259, replacement3259)

    pattern3260 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons1456, cons810, cons1493)
    rule3260 = ReplacementRule(pattern3260, replacement3260)

    pattern3261 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule3261 = ReplacementRule(pattern3261, replacement3261)

    pattern3262 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule3262 = ReplacementRule(pattern3262, replacement3262)

    pattern3263 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule3263 = ReplacementRule(pattern3263, replacement3263)

    pattern3264 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))/(a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**S(2), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons64)
    rule3264 = ReplacementRule(pattern3264, replacement3264)

    pattern3265 = Pattern(Integral((a_ + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons1495, cons64)
    rule3265 = ReplacementRule(pattern3265, replacement3265)

    pattern3266 = Pattern(Integral((a_ + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**n_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons1269, cons1495, cons64)
    rule3266 = ReplacementRule(pattern3266, replacement3266)

    pattern3267 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule3267 = ReplacementRule(pattern3267, replacement3267)

    pattern3268 = Pattern(Integral(u_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(v_))**WC('n', S(1)), x_), cons2, cons3, cons19, cons4, cons812, cons813)
    rule3268 = ReplacementRule(pattern3268, replacement3268)

    pattern3269 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule3269 = ReplacementRule(pattern3269, replacement3269)

    pattern3270 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_*WC('f', S(1)) + WC('e', S(0))))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons1362)
    rule3270 = ReplacementRule(pattern3270, replacement3270)

    pattern3271 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule3271 = ReplacementRule(pattern3271, replacement3271)

    pattern3272 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule3272 = ReplacementRule(pattern3272, replacement3272)

    pattern3273 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons139, cons746)
    rule3273 = ReplacementRule(pattern3273, replacement3273)

    pattern3274 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons139, cons746)
    rule3274 = ReplacementRule(pattern3274, replacement3274)

    pattern3275 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons150, cons1496)
    rule3275 = ReplacementRule(pattern3275, replacement3275)

    pattern3276 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons150, cons1496)
    rule3276 = ReplacementRule(pattern3276, replacement3276)

    pattern3277 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons198)
    rule3277 = ReplacementRule(pattern3277, replacement3277)

    pattern3278 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons198)
    rule3278 = ReplacementRule(pattern3278, replacement3278)

    pattern3279 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule3279 = ReplacementRule(pattern3279, replacement3279)

    pattern3280 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons1497)
    rule3280 = ReplacementRule(pattern3280, replacement3280)

    pattern3281 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule3281 = ReplacementRule(pattern3281, replacement3281)

    pattern3282 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule3282 = ReplacementRule(pattern3282, replacement3282)

    pattern3283 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons55, cons13, cons139, cons598)
    rule3283 = ReplacementRule(pattern3283, replacement3283)

    pattern3284 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons55, cons13, cons139, cons598)
    rule3284 = ReplacementRule(pattern3284, replacement3284)

    pattern3285 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons33, cons139, cons1498)
    rule3285 = ReplacementRule(pattern3285, replacement3285)

    pattern3286 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons40, cons150, cons33, cons139, cons1498)
    rule3286 = ReplacementRule(pattern3286, replacement3286)

    pattern3287 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons20, cons150, cons1496)
    rule3287 = ReplacementRule(pattern3287, replacement3287)

    pattern3288 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons65, cons20, cons150, cons1496)
    rule3288 = ReplacementRule(pattern3288, replacement3288)

    pattern3289 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons65, cons198)
    rule3289 = ReplacementRule(pattern3289, replacement3289)

    pattern3290 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons19, cons65, cons198)
    rule3290 = ReplacementRule(pattern3290, replacement3290)

    pattern3291 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*sin(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3291 = ReplacementRule(pattern3291, replacement3291)

    pattern3292 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*cos(x_*WC('d', S(1)) + WC('c', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons1499)
    rule3292 = ReplacementRule(pattern3292, replacement3292)

    pattern3293 = Pattern(Integral(sin(x_**S(2)*WC('d', S(1))), x_), cons29, cons29)
    rule3293 = ReplacementRule(pattern3293, replacement3293)

    pattern3294 = Pattern(Integral(cos(x_**S(2)*WC('d', S(1))), x_), cons29, cons29)
    rule3294 = ReplacementRule(pattern3294, replacement3294)

    pattern3295 = Pattern(Integral(sin(c_ + x_**S(2)*WC('d', S(1))), x_), cons8, cons29, cons1263)
    rule3295 = ReplacementRule(pattern3295, replacement3295)

    pattern3296 = Pattern(Integral(cos(c_ + x_**S(2)*WC('d', S(1))), x_), cons8, cons29, cons1263)
    rule3296 = ReplacementRule(pattern3296, replacement3296)

    pattern3297 = Pattern(Integral(sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons87, cons746)
    rule3297 = ReplacementRule(pattern3297, replacement3297)

    pattern3298 = Pattern(Integral(cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons87, cons746)
    rule3298 = ReplacementRule(pattern3298, replacement3298)

    pattern3299 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons378, cons167, cons148)
    rule3299 = ReplacementRule(pattern3299, replacement3299)

    pattern3300 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons378, cons167, cons148)
    rule3300 = ReplacementRule(pattern3300, replacement3300)

    pattern3301 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198)
    rule3301 = ReplacementRule(pattern3301, replacement3301)

    pattern3302 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198)
    rule3302 = ReplacementRule(pattern3302, replacement3302)

    pattern3303 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons491)
    rule3303 = ReplacementRule(pattern3303, With3303)

    pattern3304 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons491)
    rule3304 = ReplacementRule(pattern3304, With3304)

    pattern3305 = Pattern(Integral(sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons4, cons1500)
    rule3305 = ReplacementRule(pattern3305, replacement3305)

    pattern3306 = Pattern(Integral(cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons4, cons1500)
    rule3306 = ReplacementRule(pattern3306, replacement3306)

    pattern3307 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule3307 = ReplacementRule(pattern3307, replacement3307)

    pattern3308 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons130)
    rule3308 = ReplacementRule(pattern3308, replacement3308)

    pattern3309 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons40, cons70, cons71)
    rule3309 = ReplacementRule(pattern3309, replacement3309)

    pattern3310 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons40, cons70, cons71)
    rule3310 = ReplacementRule(pattern3310, replacement3310)

    pattern3311 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(u_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70)
    rule3311 = ReplacementRule(pattern3311, replacement3311)

    pattern3312 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(u_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70)
    rule3312 = ReplacementRule(pattern3312, replacement3312)

    pattern3313 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*sin(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule3313 = ReplacementRule(pattern3313, replacement3313)

    pattern3314 = Pattern(Integral((WC('a', S(0)) + WC('b', S(1))*cos(u_))**WC('p', S(1)), x_), cons2, cons3, cons5, cons825, cons826)
    rule3314 = ReplacementRule(pattern3314, replacement3314)

    pattern3315 = Pattern(Integral(sin(x_**n_*WC('d', S(1)))/x_, x_), cons29, cons4, cons1501)
    rule3315 = ReplacementRule(pattern3315, replacement3315)

    pattern3316 = Pattern(Integral(cos(x_**n_*WC('d', S(1)))/x_, x_), cons29, cons4, cons1501)
    rule3316 = ReplacementRule(pattern3316, replacement3316)

    pattern3317 = Pattern(Integral(sin(c_ + x_**n_*WC('d', S(1)))/x_, x_), cons8, cons29, cons4, cons1500)
    rule3317 = ReplacementRule(pattern3317, replacement3317)

    pattern3318 = Pattern(Integral(cos(c_ + x_**n_*WC('d', S(1)))/x_, x_), cons8, cons29, cons4, cons1500)
    rule3318 = ReplacementRule(pattern3318, replacement3318)

    pattern3319 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, CustomConstraint(With3319))
    rule3319 = ReplacementRule(pattern3319, replacement3319)

    pattern3320 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, CustomConstraint(With3320))
    rule3320 = ReplacementRule(pattern3320, replacement3320)

    pattern3321 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, CustomConstraint(With3321))
    rule3321 = ReplacementRule(pattern3321, replacement3321)

    pattern3322 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, CustomConstraint(With3322))
    rule3322 = ReplacementRule(pattern3322, replacement3322)

    pattern3323 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons4, cons1502)
    rule3323 = ReplacementRule(pattern3323, replacement3323)

    pattern3324 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons4, cons1502)
    rule3324 = ReplacementRule(pattern3324, replacement3324)

    pattern3325 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons1503)
    rule3325 = ReplacementRule(pattern3325, replacement3325)

    pattern3326 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons1503)
    rule3326 = ReplacementRule(pattern3326, replacement3326)

    pattern3327 = Pattern(Integral((x_*WC('e', S(1)))**m_*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons96)
    rule3327 = ReplacementRule(pattern3327, replacement3327)

    pattern3328 = Pattern(Integral((x_*WC('e', S(1)))**m_*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons150, cons33, cons96)
    rule3328 = ReplacementRule(pattern3328, replacement3328)

    pattern3329 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons150)
    rule3329 = ReplacementRule(pattern3329, replacement3329)

    pattern3330 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons150)
    rule3330 = ReplacementRule(pattern3330, replacement3330)

    pattern3331 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons378, cons1257, cons148, cons1504)
    rule3331 = ReplacementRule(pattern3331, replacement3331)

    pattern3332 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons378, cons1257, cons148, cons1504)
    rule3332 = ReplacementRule(pattern3332, replacement3332)

    pattern3333 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons148)
    rule3333 = ReplacementRule(pattern3333, replacement3333)

    pattern3334 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons148)
    rule3334 = ReplacementRule(pattern3334, replacement3334)

    pattern3335 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1505)
    rule3335 = ReplacementRule(pattern3335, replacement3335)

    pattern3336 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1505)
    rule3336 = ReplacementRule(pattern3336, replacement3336)

    pattern3337 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1506, cons685)
    rule3337 = ReplacementRule(pattern3337, replacement3337)

    pattern3338 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons148, cons1506, cons685)
    rule3338 = ReplacementRule(pattern3338, replacement3338)

    pattern3339 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150, cons369)
    rule3339 = ReplacementRule(pattern3339, With3339)

    pattern3340 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons150, cons369)
    rule3340 = ReplacementRule(pattern3340, With3340)

    pattern3341 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons148)
    rule3341 = ReplacementRule(pattern3341, replacement3341)

    pattern3342 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons150, cons148)
    rule3342 = ReplacementRule(pattern3342, replacement3342)

    pattern3343 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons139, cons1507)
    rule3343 = ReplacementRule(pattern3343, replacement3343)

    pattern3344 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons19, cons4, cons58, cons13, cons139, cons1507)
    rule3344 = ReplacementRule(pattern3344, replacement3344)

    pattern3345 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons139, cons1507, cons1505)
    rule3345 = ReplacementRule(pattern3345, replacement3345)

    pattern3346 = Pattern(Integral(x_**WC('m', S(1))*cos(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons152, cons13, cons139, cons1507, cons1505)
    rule3346 = ReplacementRule(pattern3346, replacement3346)

    pattern3347 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198, cons20)
    rule3347 = ReplacementRule(pattern3347, replacement3347)

    pattern3348 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons40, cons198, cons20)
    rule3348 = ReplacementRule(pattern3348, replacement3348)

    pattern3349 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons198, cons369)
    rule3349 = ReplacementRule(pattern3349, With3349)

    pattern3350 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons40, cons198, cons369)
    rule3350 = ReplacementRule(pattern3350, With3350)

    pattern3351 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons198, cons358)
    rule3351 = ReplacementRule(pattern3351, replacement3351)

    pattern3352 = Pattern(Integral((x_*WC('e', S(1)))**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons198, cons358)
    rule3352 = ReplacementRule(pattern3352, replacement3352)

    pattern3353 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons40, cons491)
    rule3353 = ReplacementRule(pattern3353, With3353)

    pattern3354 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons40, cons491)
    rule3354 = ReplacementRule(pattern3354, With3354)

    pattern3355 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons491)
    rule3355 = ReplacementRule(pattern3355, replacement3355)

    pattern3356 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons40, cons491)
    rule3356 = ReplacementRule(pattern3356, replacement3356)

    pattern3357 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, cons68, cons856, cons25)
    rule3357 = ReplacementRule(pattern3357, replacement3357)

    pattern3358 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons40, cons68, cons856, cons25)
    rule3358 = ReplacementRule(pattern3358, replacement3358)

    pattern3359 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons68, cons856, cons25)
    rule3359 = ReplacementRule(pattern3359, replacement3359)

    pattern3360 = Pattern(Integral((e_*x_)**m_*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons40, cons68, cons856, cons25)
    rule3360 = ReplacementRule(pattern3360, replacement3360)

    pattern3361 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons4, cons1508)
    rule3361 = ReplacementRule(pattern3361, replacement3361)

    pattern3362 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))), x_), cons8, cons29, cons50, cons19, cons4, cons1508)
    rule3362 = ReplacementRule(pattern3362, replacement3362)

    pattern3363 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule3363 = ReplacementRule(pattern3363, replacement3363)

    pattern3364 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(x_**n_*WC('d', S(1)) + WC('c', S(0))))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons130)
    rule3364 = ReplacementRule(pattern3364, replacement3364)

    pattern3365 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71, cons20)
    rule3365 = ReplacementRule(pattern3365, replacement3365)

    pattern3366 = Pattern(Integral(x_**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons5, cons70, cons71, cons20)
    rule3366 = ReplacementRule(pattern3366, replacement3366)

    pattern3367 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons70)
    rule3367 = ReplacementRule(pattern3367, replacement3367)

    pattern3368 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(u_**n_*WC('d', S(1)) + WC('c', S(0))))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons4, cons5, cons70)
    rule3368 = ReplacementRule(pattern3368, replacement3368)

    pattern3369 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*sin(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule3369 = ReplacementRule(pattern3369, replacement3369)

    pattern3370 = Pattern(Integral((e_*x_)**WC('m', S(1))*(WC('a', S(0)) + WC('b', S(1))*cos(u_))**WC('p', S(1)), x_), cons2, cons3, cons50, cons19, cons5, cons825, cons826)
    rule3370 = ReplacementRule(pattern3370, replacement3370)

    pattern3371 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons19, cons4, cons5, cons55, cons56)
    rule3371 = ReplacementRule(pattern3371, replacement3371)

    pattern3372 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons19, cons4, cons5, cons55, cons56)
    rule3372 = ReplacementRule(pattern3372, replacement3372)

    pattern3373 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons5, cons95, cons1503, cons56)
    rule3373 = ReplacementRule(pattern3373, replacement3373)

    pattern3374 = Pattern(Integral(x_**WC('m', S(1))*sin(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*cos(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons5, cons95, cons1503, cons56)
    rule3374 = ReplacementRule(pattern3374, replacement3374)

    pattern3375 = Pattern(Integral(sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons47)
    rule3375 = ReplacementRule(pattern3375, replacement3375)

    pattern3376 = Pattern(Integral(cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons47)
    rule3376 = ReplacementRule(pattern3376, replacement3376)

    pattern3377 = Pattern(Integral(sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons228)
    rule3377 = ReplacementRule(pattern3377, replacement3377)

    pattern3378 = Pattern(Integral(cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons228)
    rule3378 = ReplacementRule(pattern3378, replacement3378)

    pattern3379 = Pattern(Integral(sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons87, cons167)
    rule3379 = ReplacementRule(pattern3379, replacement3379)

    pattern3380 = Pattern(Integral(cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons87, cons167)
    rule3380 = ReplacementRule(pattern3380, replacement3380)

    pattern3381 = Pattern(Integral(sin(v_)**WC('n', S(1)), x_), cons150, cons820, cons1133)
    rule3381 = ReplacementRule(pattern3381, replacement3381)

    pattern3382 = Pattern(Integral(cos(v_)**WC('n', S(1)), x_), cons150, cons820, cons1133)
    rule3382 = ReplacementRule(pattern3382, replacement3382)

    pattern3383 = Pattern(Integral((d_ + x_*WC('e', S(1)))*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons49)
    rule3383 = ReplacementRule(pattern3383, replacement3383)

    pattern3384 = Pattern(Integral((d_ + x_*WC('e', S(1)))*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons49)
    rule3384 = ReplacementRule(pattern3384, replacement3384)

    pattern3385 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons241)
    rule3385 = ReplacementRule(pattern3385, replacement3385)

    pattern3386 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons241)
    rule3386 = ReplacementRule(pattern3386, replacement3386)

    pattern3387 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1134)
    rule3387 = ReplacementRule(pattern3387, replacement3387)

    pattern3388 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1134)
    rule3388 = ReplacementRule(pattern3388, replacement3388)

    pattern3389 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1135)
    rule3389 = ReplacementRule(pattern3389, replacement3389)

    pattern3390 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons168, cons1135)
    rule3390 = ReplacementRule(pattern3390, replacement3390)

    pattern3391 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1134)
    rule3391 = ReplacementRule(pattern3391, replacement3391)

    pattern3392 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1134)
    rule3392 = ReplacementRule(pattern3392, replacement3392)

    pattern3393 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1135)
    rule3393 = ReplacementRule(pattern3393, replacement3393)

    pattern3394 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons33, cons96, cons1135)
    rule3394 = ReplacementRule(pattern3394, replacement3394)

    pattern3395 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule3395 = ReplacementRule(pattern3395, replacement3395)

    pattern3396 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons19, cons1509)
    rule3396 = ReplacementRule(pattern3396, replacement3396)

    pattern3397 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*sin(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons87, cons167)
    rule3397 = ReplacementRule(pattern3397, replacement3397)

    pattern3398 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*cos(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons50, cons19, cons87, cons167)
    rule3398 = ReplacementRule(pattern3398, replacement3398)

    pattern3399 = Pattern(Integral(u_**WC('m', S(1))*sin(v_)**WC('n', S(1)), x_), cons19, cons150, cons70, cons820, cons821)
    rule3399 = ReplacementRule(pattern3399, replacement3399)

    pattern3400 = Pattern(Integral(u_**WC('m', S(1))*cos(v_)**WC('n', S(1)), x_), cons19, cons150, cons70, cons820, cons821)
    rule3400 = ReplacementRule(pattern3400, replacement3400)
    return [rule2167, rule2168, rule2169, rule2170, rule2171, rule2172, rule2173, rule2174, rule2175, rule2176, rule2177, rule2178, rule2179, rule2180, rule2181, rule2182, rule2183, rule2184, rule2185, rule2186, rule2187, rule2188, rule2189, rule2190, rule2191, rule2192, rule2193, rule2194, rule2195, rule2196, rule2197, rule2198, rule2199, rule2200, rule2201, rule2202, rule2203, rule2204, rule2205, rule2206, rule2207, rule2208, rule2209, rule2210, rule2211, rule2212, rule2213, rule2214, rule2215, rule2216, rule2217, rule2218, rule2219, rule2220, rule2221, rule2222, rule2223, rule2224, rule2225, rule2226, rule2227, rule2228, rule2229, rule2230, rule2231, rule2232, rule2233, rule2234, rule2235, rule2236, rule2237, rule2238, rule2239, rule2240, rule2241, rule2242, rule2243, rule2244, rule2245, rule2246, rule2247, rule2248, rule2249, rule2250, rule2251, rule2252, rule2253, rule2254, rule2255, rule2256, rule2257, rule2258, rule2259, rule2260, rule2261, rule2262, rule2263, rule2264, rule2265, rule2266, rule2267, rule2268, rule2269, rule2270, rule2271, rule2272, rule2273, rule2274, rule2275, rule2276, rule2277, rule2278, rule2279, rule2280, rule2281, rule2282, rule2283, rule2284, rule2285, rule2286, rule2287, rule2288, rule2289, rule2290, rule2291, rule2292, rule2293, rule2294, rule2295, rule2296, rule2297, rule2298, rule2299, rule2300, rule2301, rule2302, rule2303, rule2304, rule2305, rule2306, rule2307, rule2308, rule2309, rule2310, rule2311, rule2312, rule2313, rule2314, rule2315, rule2316, rule2317, rule2318, rule2319, rule2320, rule2321, rule2322, rule2323, rule2324, rule2325, rule2326, rule2327, rule2328, rule2329, rule2330, rule2331, rule2332, rule2333, rule2334, rule2335, rule2336, rule2337, rule2338, rule2339, rule2340, rule2341, rule2342, rule2343, rule2344, rule2345, rule2346, rule2347, rule2348, rule2349, rule2350, rule2351, rule2352, rule2353, rule2354, rule2355, rule2356, rule2357, rule2358, rule2359, rule2360, rule2361, rule2362, rule2363, rule2364, rule2365, rule2366, rule2367, rule2368, rule2369, rule2370, rule2371, rule2372, rule2373, rule2374, rule2375, rule2376, rule2377, rule2378, rule2379, rule2380, rule2381, rule2382, rule2383, rule2384, rule2385, rule2386, rule2387, rule2388, rule2389, rule2390, rule2391, rule2392, rule2393, rule2394, rule2395, rule2396, rule2397, rule2398, rule2399, rule2400, rule2401, rule2402, rule2403, rule2404, rule2405, rule2406, rule2407, rule2408, rule2409, rule2410, rule2411, rule2412, rule2413, rule2414, rule2415, rule2416, rule2417, rule2418, rule2419, rule2420, rule2421, rule2422, rule2423, rule2424, rule2425, rule2426, rule2427, rule2428, rule2429, rule2430, rule2431, rule2432, rule2433, rule2434, rule2435, rule2436, rule2437, rule2438, rule2439, rule2440, rule2441, rule2442, rule2443, rule2444, rule2445, rule2446, rule2447, rule2448, rule2449, rule2450, rule2451, rule2452, rule2453, rule2454, rule2455, rule2456, rule2457, rule2458, rule2459, rule2460, rule2461, rule2462, rule2463, rule2464, rule2465, rule2466, rule2467, rule2468, rule2469, rule2470, rule2471, rule2472, rule2473, rule2474, rule2475, rule2476, rule2477, rule2478, rule2479, rule2480, rule2481, rule2482, rule2483, rule2484, rule2485, rule2486, rule2487, rule2488, rule2489, rule2490, rule2491, rule2492, rule2493, rule2494, rule2495, rule2496, rule2497, rule2498, rule2499, rule2500, rule2501, rule2502, rule2503, rule2504, rule2505, rule2506, rule2507, rule2508, rule2509, rule2510, rule2511, rule2512, rule2513, rule2514, rule2515, rule2516, rule2517, rule2518, rule2519, rule2520, rule2521, rule2522, rule2523, rule2524, rule2525, rule2526, rule2527, rule2528, rule2529, rule2530, rule2531, rule2532, rule2533, rule2534, rule2535, rule2536, rule2537, rule2538, rule2539, rule2540, rule2541, rule2542, rule2543, rule2544, rule2545, rule2546, rule2547, rule2548, rule2549, rule2550, rule2551, rule2552, rule2553, rule2554, rule2555, rule2556, rule2557, rule2558, rule2559, rule2560, rule2561, rule2562, rule2563, rule2564, rule2565, rule2566, rule2567, rule2568, rule2569, rule2570, rule2571, rule2572, rule2573, rule2574, rule2575, rule2576, rule2577, rule2578, rule2579, rule2580, rule2581, rule2582, rule2583, rule2584, rule2585, rule2586, rule2587, rule2588, rule2589, rule2590, rule2591, rule2592, rule2593, rule2594, rule2595, rule2596, rule2597, rule2598, rule2599, rule2600, rule2601, rule2602, rule2603, rule2604, rule2605, rule2606, rule2607, rule2608, rule2609, rule2610, rule2611, rule2612, rule2613, rule2614, rule2615, rule2616, rule2617, rule2618, rule2619, rule2620, rule2621, rule2622, rule2623, rule2624, rule2625, rule2626, rule2627, rule2628, rule2629, rule2630, rule2631, rule2632, rule2633, rule2634, rule2635, rule2636, rule2637, rule2638, rule2639, rule2640, rule2641, rule2642, rule2643, rule2644, rule2645, rule2646, rule2647, rule2648, rule2649, rule2650, rule2651, rule2652, rule2653, rule2654, rule2655, rule2656, rule2657, rule2658, rule2659, rule2660, rule2661, rule2662, rule2663, rule2664, rule2665, rule2666, rule2667, rule2668, rule2669, rule2670, rule2671, rule2672, rule2673, rule2674, rule2675, rule2676, rule2677, rule2678, rule2679, rule2680, rule2681, rule2682, rule2683, rule2684, rule2685, rule2686, rule2687, rule2688, rule2689, rule2690, rule2691, rule2692, rule2693, rule2694, rule2695, rule2696, rule2697, rule2698, rule2699, rule2700, rule2701, rule2702, rule2703, rule2704, rule2705, rule2706, rule2707, rule2708, rule2709, rule2710, rule2711, rule2712, rule2713, rule2714, rule2715, rule2716, rule2717, rule2718, rule2719, rule2720, rule2721, rule2722, rule2723, rule2724, rule2725, rule2726, rule2727, rule2728, rule2729, rule2730, rule2731, rule2732, rule2733, rule2734, rule2735, rule2736, rule2737, rule2738, rule2739, rule2740, rule2741, rule2742, rule2743, rule2744, rule2745, rule2746, rule2747, rule2748, rule2749, rule2750, rule2751, rule2752, rule2753, rule2754, rule2755, rule2756, rule2757, rule2758, rule2759, rule2760, rule2761, rule2762, rule2763, rule2764, rule2765, rule2766, rule2767, rule2768, rule2769, rule2770, rule2771, rule2772, rule2773, rule2774, rule2775, rule2776, rule2777, rule2778, rule2779, rule2780, rule2781, rule2782, rule2783, rule2784, rule2785, rule2786, rule2787, rule2788, rule2789, rule2790, rule2791, rule2792, rule2793, rule2794, rule2795, rule2796, rule2797, rule2798, rule2799, rule2800, rule2801, rule2802, rule2803, rule2804, rule2805, rule2806, rule2807, rule2808, rule2809, rule2810, rule2811, rule2812, rule2813, rule2814, rule2815, rule2816, rule2817, rule2818, rule2819, rule2820, rule2821, rule2822, rule2823, rule2824, rule2825, rule2826, rule2827, rule2828, rule2829, rule2830, rule2831, rule2832, rule2833, rule2834, rule2835, rule2836, rule2837, rule2838, rule2839, rule2840, rule2841, rule2842, rule2843, rule2844, rule2845, rule2846, rule2847, rule2848, rule2849, rule2850, rule2851, rule2852, rule2853, rule2854, rule2855, rule2856, rule2857, rule2858, rule2859, rule2860, rule2861, rule2862, rule2863, rule2864, rule2865, rule2866, rule2867, rule2868, rule2869, rule2870, rule2871, rule2872, rule2873, rule2874, rule2875, rule2876, rule2877, rule2878, rule2879, rule2880, rule2881, rule2882, rule2883, rule2884, rule2885, rule2886, rule2887, rule2888, rule2889, rule2890, rule2891, rule2892, rule2893, rule2894, rule2895, rule2896, rule2897, rule2898, rule2899, rule2900, rule2901, rule2902, rule2903, rule2904, rule2905, rule2906, rule2907, rule2908, rule2909, rule2910, rule2911, rule2912, rule2913, rule2914, rule2915, rule2916, rule2917, rule2918, rule2919, rule2920, rule2921, rule2922, rule2923, rule2924, rule2925, rule2926, rule2927, rule2928, rule2929, rule2930, rule2931, rule2932, rule2933, rule2934, rule2935, rule2936, rule2937, rule2938, rule2939, rule2940, rule2941, rule2942, rule2943, rule2944, rule2945, rule2946, rule2947, rule2948, rule2949, rule2950, rule2951, rule2952, rule2953, rule2954, rule2955, rule2956, rule2957, rule2958, rule2959, rule2960, rule2961, rule2962, rule2963, rule2964, rule2965, rule2966, rule2967, rule2968, rule2969, rule2970, rule2971, rule2972, rule2973, rule2974, rule2975, rule2976, rule2977, rule2978, rule2979, rule2980, rule2981, rule2982, rule2983, rule2984, rule2985, rule2986, rule2987, rule2988, rule2989, rule2990, rule2991, rule2992, rule2993, rule2994, rule2995, rule2996, rule2997, rule2998, rule2999, rule3000, rule3001, rule3002, rule3003, rule3004, rule3005, rule3006, rule3007, rule3008, rule3009, rule3010, rule3011, rule3012, rule3013, rule3014, rule3015, rule3016, rule3017, rule3018, rule3019, rule3020, rule3021, rule3022, rule3023, rule3024, rule3025, rule3026, rule3027, rule3028, rule3029, rule3030, rule3031, rule3032, rule3033, rule3034, rule3035, rule3036, rule3037, rule3038, rule3039, rule3040, rule3041, rule3042, rule3043, rule3044, rule3045, rule3046, rule3047, rule3048, rule3049, rule3050, rule3051, rule3052, rule3053, rule3054, rule3055, rule3056, rule3057, rule3058, rule3059, rule3060, rule3061, rule3062, rule3063, rule3064, rule3065, rule3066, rule3067, rule3068, rule3069, rule3070, rule3071, rule3072, rule3073, rule3074, rule3075, rule3076, rule3077, rule3078, rule3079, rule3080, rule3081, rule3082, rule3083, rule3084, rule3085, rule3086, rule3087, rule3088, rule3089, rule3090, rule3091, rule3092, rule3093, rule3094, rule3095, rule3096, rule3097, rule3098, rule3099, rule3100, rule3101, rule3102, rule3103, rule3104, rule3105, rule3106, rule3107, rule3108, rule3109, rule3110, rule3111, rule3112, rule3113, rule3114, rule3115, rule3116, rule3117, rule3118, rule3119, rule3120, rule3121, rule3122, rule3123, rule3124, rule3125, rule3126, rule3127, rule3128, rule3129, rule3130, rule3131, rule3132, rule3133, rule3134, rule3135, rule3136, rule3137, rule3138, rule3139, rule3140, rule3141, rule3142, rule3143, rule3144, rule3145, rule3146, rule3147, rule3148, rule3149, rule3150, rule3151, rule3152, rule3153, rule3154, rule3155, rule3156, rule3157, rule3158, rule3159, rule3160, rule3161, rule3162, rule3163, rule3164, rule3165, rule3166, rule3167, rule3168, rule3169, rule3170, rule3171, rule3172, rule3173, rule3174, rule3175, rule3176, rule3177, rule3178, rule3179, rule3180, rule3181, rule3182, rule3183, rule3184, rule3185, rule3186, rule3187, rule3188, rule3189, rule3190, rule3191, rule3192, rule3193, rule3194, rule3195, rule3196, rule3197, rule3198, rule3199, rule3200, rule3201, rule3202, rule3203, rule3204, rule3205, rule3206, rule3207, rule3208, rule3209, rule3210, rule3211, rule3212, rule3213, rule3214, rule3215, rule3216, rule3217, rule3218, rule3219, rule3220, rule3221, rule3222, rule3223, rule3224, rule3225, rule3226, rule3227, rule3228, rule3229, rule3230, rule3231, rule3232, rule3233, rule3234, rule3235, rule3236, rule3237, rule3238, rule3239, rule3240, rule3241, rule3242, rule3243, rule3244, rule3245, rule3246, rule3247, rule3248, rule3249, rule3250, rule3251, rule3252, rule3253, rule3254, rule3255, rule3256, rule3257, rule3258, rule3259, rule3260, rule3261, rule3262, rule3263, rule3264, rule3265, rule3266, rule3267, rule3268, rule3269, rule3270, rule3271, rule3272, rule3273, rule3274, rule3275, rule3276, rule3277, rule3278, rule3279, rule3280, rule3281, rule3282, rule3283, rule3284, rule3285, rule3286, rule3287, rule3288, rule3289, rule3290, rule3291, rule3292, rule3293, rule3294, rule3295, rule3296, rule3297, rule3298, rule3299, rule3300, rule3301, rule3302, rule3303, rule3304, rule3305, rule3306, rule3307, rule3308, rule3309, rule3310, rule3311, rule3312, rule3313, rule3314, rule3315, rule3316, rule3317, rule3318, rule3319, rule3320, rule3321, rule3322, rule3323, rule3324, rule3325, rule3326, rule3327, rule3328, rule3329, rule3330, rule3331, rule3332, rule3333, rule3334, rule3335, rule3336, rule3337, rule3338, rule3339, rule3340, rule3341, rule3342, rule3343, rule3344, rule3345, rule3346, rule3347, rule3348, rule3349, rule3350, rule3351, rule3352, rule3353, rule3354, rule3355, rule3356, rule3357, rule3358, rule3359, rule3360, rule3361, rule3362, rule3363, rule3364, rule3365, rule3366, rule3367, rule3368, rule3369, rule3370, rule3371, rule3372, rule3373, rule3374, rule3375, rule3376, rule3377, rule3378, rule3379, rule3380, rule3381, rule3382, rule3383, rule3384, rule3385, rule3386, rule3387, rule3388, rule3389, rule3390, rule3391, rule3392, rule3393, rule3394, rule3395, rule3396, rule3397, rule3398, rule3399, rule3400, ]





def replacement2167(u, x):
    return Int(DeactivateTrig(u, x), x)


def replacement2168(a, b, e, f, m, n, x):
    return Simp((a*sin(e + f*x))**(m + S(1))*(b*cos(e + f*x))**(n + S(1))/(a*b*f*(m + S(1))), x)


def replacement2169(a, e, f, m, n, x):
    return Dist(S(1)/(a*f), Subst(Int(x**m*(S(1) - x**S(2)/a**S(2))**(n/S(2) + S(-1)/2), x), x, a*sin(e + f*x)), x)


def replacement2170(a, e, f, m, n, x):
    return -Dist(S(1)/(a*f), Subst(Int(x**m*(S(1) - x**S(2)/a**S(2))**(n/S(2) + S(-1)/2), x), x, a*cos(e + f*x)), x)


def replacement2171(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(b**S(2)*(n + S(1))), Int((a*sin(e + f*x))**(m + S(-2))*(b*cos(e + f*x))**(n + S(2)), x), x) - Simp(a*(a*sin(e + f*x))**(m + S(-1))*(b*cos(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)


def replacement2172(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(b**S(2)*(n + S(1))), Int((a*cos(e + f*x))**(m + S(-2))*(b*sin(e + f*x))**(n + S(2)), x), x) + Simp(a*(a*cos(e + f*x))**(m + S(-1))*(b*sin(e + f*x))**(n + S(1))/(b*f*(n + S(1))), x)


def replacement2173(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(m + n), Int((a*sin(e + f*x))**(m + S(-2))*(b*cos(e + f*x))**n, x), x) - Simp(a*(a*sin(e + f*x))**(m + S(-1))*(b*cos(e + f*x))**(n + S(1))/(b*f*(m + n)), x)


def replacement2174(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(m + n), Int((a*cos(e + f*x))**(m + S(-2))*(b*sin(e + f*x))**n, x), x) + Simp(a*(a*cos(e + f*x))**(m + S(-1))*(b*sin(e + f*x))**(n + S(1))/(b*f*(m + n)), x)


def replacement2175(a, b, e, f, m, n, x):
    return Dist((m + n + S(2))/(a**S(2)*(m + S(1))), Int((a*sin(e + f*x))**(m + S(2))*(b*cos(e + f*x))**n, x), x) + Simp((a*sin(e + f*x))**(m + S(1))*(b*cos(e + f*x))**(n + S(1))/(a*b*f*(m + S(1))), x)


def replacement2176(a, b, e, f, m, n, x):
    return Dist((m + n + S(2))/(a**S(2)*(m + S(1))), Int((a*cos(e + f*x))**(m + S(2))*(b*sin(e + f*x))**n, x), x) - Simp((a*cos(e + f*x))**(m + S(1))*(b*sin(e + f*x))**(n + S(1))/(a*b*f*(m + S(1))), x)


def replacement2177(a, b, e, f, x):
    return Dist(sqrt(a*sin(e + f*x))*sqrt(b*cos(e + f*x))/sqrt(sin(S(2)*e + S(2)*f*x)), Int(sqrt(sin(S(2)*e + S(2)*f*x)), x), x)


def replacement2178(a, b, e, f, x):
    return Dist(sqrt(sin(S(2)*e + S(2)*f*x))/(sqrt(a*sin(e + f*x))*sqrt(b*cos(e + f*x))), Int(S(1)/sqrt(sin(S(2)*e + S(2)*f*x)), x), x)


def With2179(a, b, e, f, m, n, x):
    k = Denominator(m)
    return Dist(a*b*k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))/(a**S(2) + b**S(2)*x**(S(2)*k)), x), x, (a*sin(e + f*x))**(S(1)/k)*(b*cos(e + f*x))**(-S(1)/k)), x)


def With2180(a, b, e, f, m, n, x):
    k = Denominator(m)
    return -Dist(a*b*k/f, Subst(Int(x**(k*(m + S(1)) + S(-1))/(a**S(2) + b**S(2)*x**(S(2)*k)), x), x, (a*cos(e + f*x))**(S(1)/k)*(b*sin(e + f*x))**(-S(1)/k)), x)


def replacement2181(a, b, e, f, m, n, x):
    return Simp(b**(S(2)*IntPart(n/S(2) + S(-1)/2) + S(1))*(a*sin(e + f*x))**(m + S(1))*(b*cos(e + f*x))**(S(2)*FracPart(n/S(2) + S(-1)/2))*(cos(e + f*x)**S(2))**(-FracPart(n/S(2) + S(-1)/2))*Hypergeometric2F1(m/S(2) + S(1)/2, S(1)/2 - n/S(2), m/S(2) + S(3)/2, sin(e + f*x)**S(2))/(a*f*(m + S(1))), x)


def replacement2182(a, b, e, f, m, n, x):
    return -Simp(b**(S(2)*IntPart(n/S(2) + S(-1)/2) + S(1))*(a*cos(e + f*x))**(m + S(1))*(b*sin(e + f*x))**(S(2)*FracPart(n/S(2) + S(-1)/2))*(sin(e + f*x)**S(2))**(-FracPart(n/S(2) + S(-1)/2))*Hypergeometric2F1(m/S(2) + S(1)/2, S(1)/2 - n/S(2), m/S(2) + S(3)/2, cos(e + f*x)**S(2))/(a*f*(m + S(1))), x)


def replacement2183(a, b, e, f, m, n, x):
    return Simp(b*(a*sin(e + f*x))**(m + S(1))*(b/cos(e + f*x))**(n + S(-1))/(a*f*(m + S(1))), x)


def replacement2184(a, b, e, f, m, n, x):
    return -Simp(b*(a*cos(e + f*x))**(m + S(1))*(b/sin(e + f*x))**(n + S(-1))/(a*f*(m + S(1))), x)


def replacement2185(a, b, e, f, m, n, x):
    return -Dist(a**S(2)*b**S(2)*(m + S(-1))/(n + S(-1)), Int((a*sin(e + f*x))**(m + S(-2))*(b/cos(e + f*x))**(n + S(-2)), x), x) + Simp(a*b*(a*sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)


def replacement2186(a, b, e, f, m, n, x):
    return -Dist(a**S(2)*b**S(2)*(m + S(-1))/(n + S(-1)), Int((a*cos(e + f*x))**(m + S(-2))*(b/sin(e + f*x))**(n + S(-2)), x), x) - Simp(a*b*(a*cos(e + f*x))**(m + S(-1))*(b/sin(e + f*x))**(n + S(-1))/(f*(n + S(-1))), x)


def replacement2187(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(m - n), Int((a*sin(e + f*x))**(m + S(-2))*(b/cos(e + f*x))**n, x), x) - Simp(a*b*(a*sin(e + f*x))**(m + S(-1))*(b/cos(e + f*x))**(n + S(-1))/(f*(m - n)), x)


def replacement2188(a, b, e, f, m, n, x):
    return Dist(a**S(2)*(m + S(-1))/(m - n), Int((a*cos(e + f*x))**(m + S(-2))*(b/sin(e + f*x))**n, x), x) + Simp(a*b*(a*cos(e + f*x))**(m + S(-1))*(b/sin(e + f*x))**(n + S(-1))/(f*(m - n)), x)


def replacement2189(a, b, e, f, m, n, x):
    return Dist((m - n + S(2))/(a**S(2)*(m + S(1))), Int((a*sin(e + f*x))**(m + S(2))*(b/cos(e + f*x))**n, x), x) + Simp(b*(a*sin(e + f*x))**(m + S(1))*(b/cos(e + f*x))**(n + S(-1))/(a*f*(m + S(1))), x)


def replacement2190(a, b, e, f, m, n, x):
    return Dist((m - n + S(2))/(a**S(2)*(m + S(1))), Int((a*cos(e + f*x))**(m + S(2))*(b/sin(e + f*x))**n, x), x) - Simp(b*(a*cos(e + f*x))**(m + S(1))*(b/sin(e + f*x))**(n + S(-1))/(a*f*(m + S(1))), x)


def replacement2191(a, b, e, f, m, n, x):
    return Dist((cos(e + f*x)/b)**(FracPart(n) + S(1))*(b/cos(e + f*x))**(FracPart(n) + S(1)), Int((a*sin(e + f*x))**m*(cos(e + f*x)/b)**(-n), x), x)


def replacement2192(a, b, e, f, m, n, x):
    return Dist((sin(e + f*x)/b)**(FracPart(n) + S(1))*(b/sin(e + f*x))**(FracPart(n) + S(1)), Int((a*cos(e + f*x))**m*(sin(e + f*x)/b)**(-n), x), x)


def replacement2193(a, b, e, f, m, n, x):
    return Dist((a*b)**IntPart(n)*(a*sin(e + f*x))**FracPart(n)*(b/sin(e + f*x))**FracPart(n), Int((a*sin(e + f*x))**(m - n), x), x)


def replacement2194(a, b, e, f, m, n, x):
    return Dist((a*b)**IntPart(n)*(a*cos(e + f*x))**FracPart(n)*(b/cos(e + f*x))**FracPart(n), Int((a*cos(e + f*x))**(m - n), x), x)


def replacement2195(c, d, n, x):
    return -Dist(S(1)/d, Subst(Int((S(1) - x**S(2))**(n/S(2))/sqrt(S(1) - x**S(2)), x), x, cos(c + d*x)), x)


def replacement2196(c, d, n, x):
    return Dist(S(1)/d, Subst(Int((S(1) - x**S(2))**(n/S(2))/sqrt(S(1) - x**S(2)), x), x, sin(c + d*x)), x)


def replacement2197(b, c, d, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*sin(c + d*x))**(n + S(-2)), x), x) - Simp(b*(b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)/(d*n), x)


def replacement2198(b, c, d, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*cos(c + d*x))**(n + S(-2)), x), x) + Simp(b*(b*cos(c + d*x))**(n + S(-1))*sin(c + d*x)/(d*n), x)


def replacement2199(b, c, d, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*sin(c + d*x))**(n + S(2)), x), x) + Simp((b*sin(c + d*x))**(n + S(1))*cos(c + d*x)/(b*d*(n + S(1))), x)


def replacement2200(b, c, d, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*cos(c + d*x))**(n + S(2)), x), x) - Simp((b*cos(c + d*x))**(n + S(1))*sin(c + d*x)/(b*d*(n + S(1))), x)


def replacement2201(c, d, x):
    return -Simp(cos(c + d*x)/d, x)


def replacement2202(c, d, x):
    return Simp(sin(c + d*x)/d, x)


def replacement2203(c, d, x):
    return Simp(S(2)*EllipticE(-Pi/S(4) + c/S(2) + d*x/S(2), S(2))/d, x)


def replacement2204(c, d, x):
    return Simp(S(2)*EllipticE(c/S(2) + d*x/S(2), S(2))/d, x)


def replacement2205(b, c, d, x):
    return Dist(sqrt(b*sin(c + d*x))/sqrt(sin(c + d*x)), Int(sqrt(sin(c + d*x)), x), x)


def replacement2206(b, c, d, x):
    return Dist(sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)), Int(sqrt(cos(c + d*x)), x), x)


def replacement2207(c, d, x):
    return Simp(S(2)*EllipticF(-Pi/S(4) + c/S(2) + d*x/S(2), S(2))/d, x)


def replacement2208(c, d, x):
    return Simp(S(2)*EllipticF(c/S(2) + d*x/S(2), S(2))/d, x)


def replacement2209(b, c, d, x):
    return Dist(sqrt(sin(c + d*x))/sqrt(b*sin(c + d*x)), Int(S(1)/sqrt(sin(c + d*x)), x), x)


def replacement2210(b, c, d, x):
    return Dist(sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x)), Int(S(1)/sqrt(cos(c + d*x)), x), x)


def replacement2211(b, c, d, n, x):
    return Simp((b*sin(c + d*x))**(n + S(1))*Hypergeometric2F1(S(1)/2, n/S(2) + S(1)/2, n/S(2) + S(3)/2, sin(c + d*x)**S(2))*cos(c + d*x)/(b*d*(n + S(1))*sqrt(cos(c + d*x)**S(2))), x)


def replacement2212(b, c, d, n, x):
    return -Simp((b*cos(c + d*x))**(n + S(1))*Hypergeometric2F1(S(1)/2, n/S(2) + S(1)/2, n/S(2) + S(3)/2, cos(c + d*x)**S(2))*sin(c + d*x)/(b*d*(n + S(1))*sqrt(sin(c + d*x)**S(2))), x)


def replacement2213(a, b, c, d, x):
    return Dist(S(2)*a*b, Int(sin(c + d*x), x), x) + Simp(x*(S(2)*a**S(2) + b**S(2))/S(2), x) - Simp(b**S(2)*sin(c + d*x)*cos(c + d*x)/(S(2)*d), x)


def replacement2214(a, b, c, d, x):
    return Dist(S(2)*a*b, Int(cos(c + d*x), x), x) + Simp(x*(S(2)*a**S(2) + b**S(2))/S(2), x) + Simp(b**S(2)*sin(c + d*x)*cos(c + d*x)/(S(2)*d), x)


def replacement2215(a, b, c, d, n, x):
    return Int(ExpandTrig((a + b*sin(c + d*x))**n, x), x)


def replacement2216(a, b, c, d, n, x):
    return Int(ExpandTrig((a + b*cos(c + d*x))**n, x), x)


def replacement2217(a, b, c, d, x):
    return Simp(-S(2)*b*cos(c + d*x)/(d*sqrt(a + b*sin(c + d*x))), x)


def replacement2218(a, b, c, d, x):
    return Simp(S(2)*b*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))), x)


def replacement2219(a, b, c, d, n, x):
    return Dist(a*(S(2)*n + S(-1))/n, Int((a + b*sin(c + d*x))**(n + S(-1)), x), x) - Simp(b*(a + b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)/(d*n), x)


def replacement2220(a, b, c, d, n, x):
    return Dist(a*(S(2)*n + S(-1))/n, Int((a + b*cos(c + d*x))**(n + S(-1)), x), x) + Simp(b*(a + b*cos(c + d*x))**(n + S(-1))*sin(c + d*x)/(d*n), x)


def replacement2221(a, b, c, d, x):
    return -Simp(cos(c + d*x)/(d*(a*sin(c + d*x) + b)), x)


def replacement2222(a, b, c, d, x):
    return Simp(sin(c + d*x)/(d*(a*cos(c + d*x) + b)), x)


def replacement2223(a, b, c, d, x):
    return Dist(-S(2)/d, Subst(Int(S(1)/(S(2)*a - x**S(2)), x), x, b*cos(c + d*x)/sqrt(a + b*sin(c + d*x))), x)


def replacement2224(a, b, c, d, x):
    return Dist(S(2)/d, Subst(Int(S(1)/(S(2)*a - x**S(2)), x), x, b*sin(c + d*x)/sqrt(a + b*cos(c + d*x))), x)


def replacement2225(a, b, c, d, n, x):
    return Dist((n + S(1))/(a*(S(2)*n + S(1))), Int((a + b*sin(c + d*x))**(n + S(1)), x), x) + Simp(b*(a + b*sin(c + d*x))**n*cos(c + d*x)/(a*d*(S(2)*n + S(1))), x)


def replacement2226(a, b, c, d, n, x):
    return Dist((n + S(1))/(a*(S(2)*n + S(1))), Int((a + b*cos(c + d*x))**(n + S(1)), x), x) - Simp(b*(a + b*cos(c + d*x))**n*sin(c + d*x)/(a*d*(S(2)*n + S(1))), x)


def replacement2227(a, b, c, d, n, x):
    return -Simp(S(2)**(n + S(1)/2)*a**(n + S(-1)/2)*b*Hypergeometric2F1(S(1)/2, S(1)/2 - n, S(3)/2, S(1)/2 - b*sin(c + d*x)/(S(2)*a))*cos(c + d*x)/(d*sqrt(a + b*sin(c + d*x))), x)


def replacement2228(a, b, c, d, n, x):
    return Simp(S(2)**(n + S(1)/2)*a**(n + S(-1)/2)*b*Hypergeometric2F1(S(1)/2, S(1)/2 - n, S(3)/2, S(1)/2 - b*cos(c + d*x)/(S(2)*a))*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))), x)


def replacement2229(a, b, c, d, n, x):
    return Dist(a**IntPart(n)*(S(1) + b*sin(c + d*x)/a)**(-FracPart(n))*(a + b*sin(c + d*x))**FracPart(n), Int((S(1) + b*sin(c + d*x)/a)**n, x), x)


def replacement2230(a, b, c, d, n, x):
    return Dist(a**IntPart(n)*(S(1) + b*cos(c + d*x)/a)**(-FracPart(n))*(a + b*cos(c + d*x))**FracPart(n), Int((S(1) + b*cos(c + d*x)/a)**n, x), x)


def replacement2231(a, b, c, d, x):
    return Simp(S(2)*sqrt(a + b)*EllipticE(-Pi/S(4) + c/S(2) + d*x/S(2), S(2)*b/(a + b))/d, x)


def replacement2232(a, b, c, d, x):
    return Simp(S(2)*sqrt(a + b)*EllipticE(c/S(2) + d*x/S(2), S(2)*b/(a + b))/d, x)


def replacement2233(a, b, c, d, x):
    return Simp(S(2)*sqrt(a - b)*EllipticE(Pi/S(4) + c/S(2) + d*x/S(2), -S(2)*b/(a - b))/d, x)


def replacement2234(a, b, c, d, x):
    return Simp(S(2)*sqrt(a - b)*EllipticE(Pi/S(2) + c/S(2) + d*x/S(2), -S(2)*b/(a - b))/d, x)


def replacement2235(a, b, c, d, x):
    return Dist(sqrt(a + b*sin(c + d*x))/sqrt((a + b*sin(c + d*x))/(a + b)), Int(sqrt(a/(a + b) + b*sin(c + d*x)/(a + b)), x), x)


def replacement2236(a, b, c, d, x):
    return Dist(sqrt(a + b*cos(c + d*x))/sqrt((a + b*cos(c + d*x))/(a + b)), Int(sqrt(a/(a + b) + b*cos(c + d*x)/(a + b)), x), x)


def replacement2237(a, b, c, d, n, x):
    return Dist(S(1)/n, Int((a + b*sin(c + d*x))**(n + S(-2))*Simp(a**S(2)*n + a*b*(S(2)*n + S(-1))*sin(c + d*x) + b**S(2)*(n + S(-1)), x), x), x) - Simp(b*(a + b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)/(d*n), x)


def replacement2238(a, b, c, d, n, x):
    return Dist(S(1)/n, Int((a + b*cos(c + d*x))**(n + S(-2))*Simp(a**S(2)*n + a*b*(S(2)*n + S(-1))*cos(c + d*x) + b**S(2)*(n + S(-1)), x), x), x) + Simp(b*(a + b*cos(c + d*x))**(n + S(-1))*sin(c + d*x)/(d*n), x)


def With2239(a, b, c, d, x):
    q = Rt(a**S(2) - b**S(2), S(2))
    return Simp(x/q, x) + Simp(S(2)*ArcTan(b*cos(c + d*x)/(a + b*sin(c + d*x) + q))/(d*q), x)


def With2240(a, b, c, d, x):
    q = Rt(a**S(2) - b**S(2), S(2))
    return Simp(x/q, x) - Simp(S(2)*ArcTan(b*sin(c + d*x)/(a + b*cos(c + d*x) + q))/(d*q), x)


def With2241(a, b, c, d, x):
    q = Rt(a**S(2) - b**S(2), S(2))
    return -Simp(x/q, x) - Simp(S(2)*ArcTan(b*cos(c + d*x)/(a + b*sin(c + d*x) - q))/(d*q), x)


def With2242(a, b, c, d, x):
    q = Rt(a**S(2) - b**S(2), S(2))
    return -Simp(x/q, x) + Simp(S(2)*ArcTan(b*sin(c + d*x)/(a + b*cos(c + d*x) - q))/(d*q), x)


def With2243(a, b, c, d, x):
    e = FreeFactors(tan(-Pi/S(4) + c/S(2) + d*x/S(2)), x)
    return Dist(S(2)*e/d, Subst(Int(S(1)/(a + b + e**S(2)*x**S(2)*(a - b)), x), x, tan(-Pi/S(4) + c/S(2) + d*x/S(2))/e), x)


def With2244(a, b, c, d, x):
    e = FreeFactors(tan(c/S(2) + d*x/S(2)), x)
    return Dist(S(2)*e/d, Subst(Int(S(1)/(a*e**S(2)*x**S(2) + a + S(2)*b*e*x), x), x, tan(c/S(2) + d*x/S(2))/e), x)


def With2245(a, b, c, d, x):
    e = FreeFactors(tan(c/S(2) + d*x/S(2)), x)
    return Dist(S(2)*e/d, Subst(Int(S(1)/(a + b + e**S(2)*x**S(2)*(a - b)), x), x, tan(c/S(2) + d*x/S(2))/e), x)


def replacement2246(a, b, c, d, x):
    return Simp(S(2)*EllipticF(-Pi/S(4) + c/S(2) + d*x/S(2), S(2)*b/(a + b))/(d*sqrt(a + b)), x)


def replacement2247(a, b, c, d, x):
    return Simp(S(2)*EllipticF(c/S(2) + d*x/S(2), S(2)*b/(a + b))/(d*sqrt(a + b)), x)


def replacement2248(a, b, c, d, x):
    return Simp(S(2)*EllipticF(Pi/S(4) + c/S(2) + d*x/S(2), -S(2)*b/(a - b))/(d*sqrt(a - b)), x)


def replacement2249(a, b, c, d, x):
    return Simp(S(2)*EllipticF(Pi/S(2) + c/S(2) + d*x/S(2), -S(2)*b/(a - b))/(d*sqrt(a - b)), x)


def replacement2250(a, b, c, d, x):
    return Dist(sqrt((a + b*sin(c + d*x))/(a + b))/sqrt(a + b*sin(c + d*x)), Int(S(1)/sqrt(a/(a + b) + b*sin(c + d*x)/(a + b)), x), x)


def replacement2251(a, b, c, d, x):
    return Dist(sqrt((a + b*cos(c + d*x))/(a + b))/sqrt(a + b*cos(c + d*x)), Int(S(1)/sqrt(a/(a + b) + b*cos(c + d*x)/(a + b)), x), x)


def replacement2252(a, b, c, d, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(n + S(1))), Int((a + b*sin(c + d*x))**(n + S(1))*Simp(a*(n + S(1)) - b*(n + S(2))*sin(c + d*x), x), x), x) - Simp(b*(a + b*sin(c + d*x))**(n + S(1))*cos(c + d*x)/(d*(a**S(2) - b**S(2))*(n + S(1))), x)


def replacement2253(a, b, c, d, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(n + S(1))), Int((a + b*cos(c + d*x))**(n + S(1))*Simp(a*(n + S(1)) - b*(n + S(2))*cos(c + d*x), x), x), x) + Simp(b*(a + b*cos(c + d*x))**(n + S(1))*sin(c + d*x)/(d*(a**S(2) - b**S(2))*(n + S(1))), x)


def replacement2254(a, b, c, d, n, x):
    return Dist(cos(c + d*x)/(d*sqrt(S(1) - sin(c + d*x))*sqrt(sin(c + d*x) + S(1))), Subst(Int((a + b*x)**n/(sqrt(S(1) - x)*sqrt(x + S(1))), x), x, sin(c + d*x)), x)


def replacement2255(a, b, c, d, n, x):
    return -Dist(sin(c + d*x)/(d*sqrt(S(1) - cos(c + d*x))*sqrt(cos(c + d*x) + S(1))), Subst(Int((a + b*x)**n/(sqrt(S(1) - x)*sqrt(x + S(1))), x), x, cos(c + d*x)), x)


def replacement2256(a, b, c, d, n, x):
    return Int((a + b*sin(S(2)*c + S(2)*d*x)/S(2))**n, x)


def replacement2257(a, b, e, f, m, p, x):
    return Dist(b**(-p)/f, Subst(Int((a - x)**(p/S(2) + S(-1)/2)*(a + x)**(m + p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2258(a, b, e, f, m, p, x):
    return -Dist(b**(-p)/f, Subst(Int((a - x)**(p/S(2) + S(-1)/2)*(a + x)**(m + p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2259(a, b, e, f, m, p, x):
    return Dist(b**(-p)/f, Subst(Int((a + x)**m*(b**S(2) - x**S(2))**(p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2260(a, b, e, f, m, p, x):
    return -Dist(b**(-p)/f, Subst(Int((a + x)**m*(b**S(2) - x**S(2))**(p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2261(a, b, e, f, g, p, x):
    return Dist(a, Int((g*cos(e + f*x))**p, x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))/(f*g*(p + S(1))), x)


def replacement2262(a, b, e, f, g, p, x):
    return Dist(a, Int((g*sin(e + f*x))**p, x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))/(f*g*(p + S(1))), x)


def replacement2263(a, b, e, f, g, m, p, x):
    return Dist((a/g)**(S(2)*m), Int((g*cos(e + f*x))**(S(2)*m + p)*(a - b*sin(e + f*x))**(-m), x), x)


def replacement2264(a, b, e, f, g, m, p, x):
    return Dist((a/g)**(S(2)*m), Int((g*sin(e + f*x))**(S(2)*m + p)*(a - b*cos(e + f*x))**(-m), x), x)


def replacement2265(a, b, e, f, g, m, p, x):
    return Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*m), x)


def replacement2266(a, b, e, f, g, m, p, x):
    return -Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*m), x)


def replacement2267(a, b, e, f, g, m, p, x):
    return Dist((m + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1)), x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2268(a, b, e, f, g, m, p, x):
    return Dist((m + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2269(a, b, e, f, g, m, p, x):
    return Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))/(f*g*(m + S(-1))), x)


def replacement2270(a, b, e, f, g, m, p, x):
    return -Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))/(f*g*(m + S(-1))), x)


def replacement2271(a, b, e, f, g, m, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + p), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1)), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2272(a, b, e, f, g, m, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + p), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1)), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2273(a, b, e, f, g, m, p, x):
    return Dist(a*(m + p + S(1))/(g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1)), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*(p + S(1))), x)


def replacement2274(a, b, e, f, g, m, p, x):
    return Dist(a*(m + p + S(1))/(g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1)), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*(p + S(1))), x)


def replacement2275(a, b, e, f, g, m, p, x):
    return Dist(b**S(2)*(S(2)*m + p + S(-1))/(g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-2)), x), x) + Simp(-S(2)*b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))/(f*g*(p + S(1))), x)


def replacement2276(a, b, e, f, g, m, p, x):
    return Dist(b**S(2)*(S(2)*m + p + S(-1))/(g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-2)), x), x) + Simp(S(2)*b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))/(f*g*(p + S(1))), x)


def replacement2277(a, b, e, f, g, x):
    return Dist(a*sqrt(a + b*sin(e + f*x))*sqrt(cos(e + f*x) + S(1))/(a*cos(e + f*x) + a + b*sin(e + f*x)), Int(sqrt(cos(e + f*x) + S(1))/sqrt(g*cos(e + f*x)), x), x) + Dist(b*sqrt(a + b*sin(e + f*x))*sqrt(cos(e + f*x) + S(1))/(a*cos(e + f*x) + a + b*sin(e + f*x)), Int(sin(e + f*x)/(sqrt(g*cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), x), x)


def replacement2278(a, b, e, f, g, x):
    return Dist(a*sqrt(a + b*cos(e + f*x))*sqrt(sin(e + f*x) + S(1))/(a*sin(e + f*x) + a + b*cos(e + f*x)), Int(sqrt(sin(e + f*x) + S(1))/sqrt(g*sin(e + f*x)), x), x) + Dist(b*sqrt(a + b*cos(e + f*x))*sqrt(sin(e + f*x) + S(1))/(a*sin(e + f*x) + a + b*cos(e + f*x)), Int(cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), x), x)


def replacement2279(a, b, e, f, g, m, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + p), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1)), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2280(a, b, e, f, g, m, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + p), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1)), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2281(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(a*(m + p)), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1)), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))/(b*f*(m + p)), x)


def replacement2282(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(a*(m + p)), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))/(b*f*(m + p)), x)


def replacement2283(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(2)), x), x) + Simp(S(2)*g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))/(b*f*(S(2)*m + p + S(1))), x)


def replacement2284(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(2)), x), x) + Simp(-S(2)*g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))/(b*f*(S(2)*m + p + S(1))), x)


def replacement2285(a, b, e, f, g, m, p, x):
    return Dist((m + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1)), x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2286(a, b, e, f, g, m, p, x):
    return Dist((m + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2287(a, b, e, f, g, p, x):
    return Dist(g**S(2)/a, Int((g*cos(e + f*x))**(p + S(-2)), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))/(b*f*(p + S(-1))), x)


def replacement2288(a, b, e, f, g, p, x):
    return Dist(g**S(2)/a, Int((g*sin(e + f*x))**(p + S(-2)), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))/(b*f*(p + S(-1))), x)


def replacement2289(a, b, e, f, g, p, x):
    return Dist(p/(a*(p + S(-1))), Int((g*cos(e + f*x))**p, x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))/(a*f*g*(a + b*sin(e + f*x))*(p + S(-1))), x)


def replacement2290(a, b, e, f, g, p, x):
    return Dist(p/(a*(p + S(-1))), Int((g*sin(e + f*x))**p, x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))/(a*f*g*(a + b*cos(e + f*x))*(p + S(-1))), x)


def replacement2291(a, b, e, f, g, x):
    return -Dist(g*sqrt(a + b*sin(e + f*x))*sqrt(cos(e + f*x) + S(1))/(a*sin(e + f*x) + b*cos(e + f*x) + b), Int(sin(e + f*x)/(sqrt(g*cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), x), x) + Dist(g*sqrt(a + b*sin(e + f*x))*sqrt(cos(e + f*x) + S(1))/(a*cos(e + f*x) + a + b*sin(e + f*x)), Int(sqrt(cos(e + f*x) + S(1))/sqrt(g*cos(e + f*x)), x), x)


def replacement2292(a, b, e, f, g, x):
    return Dist(g*sqrt(a + b*cos(e + f*x))*sqrt(sin(e + f*x) + S(1))/(a*sin(e + f*x) + a + b*cos(e + f*x)), Int(sqrt(sin(e + f*x) + S(1))/sqrt(g*sin(e + f*x)), x), x) - Dist(g*sqrt(a + b*cos(e + f*x))*sqrt(sin(e + f*x) + S(1))/(a*cos(e + f*x) + b*sin(e + f*x) + b), Int(cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), x), x)


def replacement2293(a, b, e, f, g, x):
    return Dist(g**S(2)/(S(2)*a), Int(sqrt(a + b*sin(e + f*x))/sqrt(g*cos(e + f*x)), x), x) + Simp(g*sqrt(g*cos(e + f*x))*sqrt(a + b*sin(e + f*x))/(b*f), x)


def replacement2294(a, b, e, f, g, x):
    return Dist(g**S(2)/(S(2)*a), Int(sqrt(a + b*cos(e + f*x))/sqrt(g*sin(e + f*x)), x), x) - Simp(g*sqrt(g*sin(e + f*x))*sqrt(a + b*cos(e + f*x))/(b*f), x)


def replacement2295(a, b, e, f, g, p, x):
    return Dist(S(2)*a*(p + S(-2))/(S(2)*p + S(-1)), Int((g*cos(e + f*x))**p/(a + b*sin(e + f*x))**(S(3)/2), x), x) + Simp(-S(2)*b*(g*cos(e + f*x))**(p + S(1))/(f*g*(a + b*sin(e + f*x))**(S(3)/2)*(S(2)*p + S(-1))), x)


def replacement2296(a, b, e, f, g, p, x):
    return Dist(S(2)*a*(p + S(-2))/(S(2)*p + S(-1)), Int((g*sin(e + f*x))**p/(a + b*cos(e + f*x))**(S(3)/2), x), x) + Simp(S(2)*b*(g*sin(e + f*x))**(p + S(1))/(f*g*(a + b*cos(e + f*x))**(S(3)/2)*(S(2)*p + S(-1))), x)


def replacement2297(a, b, e, f, g, p, x):
    return Dist(a*(S(2)*p + S(1))/(S(2)*g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))/(a + b*sin(e + f*x))**(S(3)/2), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))/(a*f*g*sqrt(a + b*sin(e + f*x))*(p + S(1))), x)


def replacement2298(a, b, e, f, g, p, x):
    return Dist(a*(S(2)*p + S(1))/(S(2)*g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))/(a + b*cos(e + f*x))**(S(3)/2), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))/(a*f*g*sqrt(a + b*cos(e + f*x))*(p + S(1))), x)


def replacement2299(a, b, e, f, g, m, p, x):
    return Dist(a**m*(g*cos(e + f*x))**(p + S(1))*(S(1) - sin(e + f*x))**(-p/S(2) + S(-1)/2)*(sin(e + f*x) + S(1))**(-p/S(2) + S(-1)/2)/(f*g), Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2300(a, b, e, f, g, m, p, x):
    return -Dist(a**m*(g*sin(e + f*x))**(p + S(1))*(S(1) - cos(e + f*x))**(-p/S(2) + S(-1)/2)*(cos(e + f*x) + S(1))**(-p/S(2) + S(-1)/2)/(f*g), Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2301(a, b, e, f, g, m, p, x):
    return Dist(a**S(2)*(g*cos(e + f*x))**(p + S(1))*(a - b*sin(e + f*x))**(-p/S(2) + S(-1)/2)*(a + b*sin(e + f*x))**(-p/S(2) + S(-1)/2)/(f*g), Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2302(a, b, e, f, g, m, p, x):
    return -Dist(a**S(2)*(g*sin(e + f*x))**(p + S(1))*(a - b*cos(e + f*x))**(-p/S(2) + S(-1)/2)*(a + b*cos(e + f*x))**(-p/S(2) + S(-1)/2)/(f*g), Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2303(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1))*(a*(p + S(2)) + b*(m + p + S(2))*sin(e + f*x)), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*sin(e + f*x)/(f*g*(p + S(1))), x)


def replacement2304(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1))*(a*(p + S(2)) + b*(m + p + S(2))*cos(e + f*x)), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*cos(e + f*x)/(f*g*(p + S(1))), x)


def replacement2305(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-2))*(a**S(2)*(p + S(2)) + a*b*(m + p + S(1))*sin(e + f*x) + b**S(2)*(m + S(-1))), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(a*sin(e + f*x) + b)/(f*g*(p + S(1))), x)


def replacement2306(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-2))*(a**S(2)*(p + S(2)) + a*b*(m + p + S(1))*cos(e + f*x) + b**S(2)*(m + S(-1))), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(a*cos(e + f*x) + b)/(f*g*(p + S(1))), x)


def replacement2307(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(m + p), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-2))*(a**S(2)*(m + p) + a*b*(S(2)*m + p + S(-1))*sin(e + f*x) + b**S(2)*(m + S(-1))), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2308(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(m + p), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-2))*(a**S(2)*(m + p) + a*b*(S(2)*m + p + S(-1))*cos(e + f*x) + b**S(2)*(m + S(-1))), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))/(f*g*(m + p)), x)


def replacement2309(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b*(m + S(1))), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1))*sin(e + f*x), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)


def replacement2310(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b*(m + S(1))), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1))*cos(e + f*x), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))/(b*f*(m + S(1))), x)


def replacement2311(a, b, e, f, g, m, p, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1))*(a*(m + S(1)) - b*(m + p + S(2))*sin(e + f*x)), x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))/(f*g*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2312(a, b, e, f, g, m, p, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1))*(a*(m + S(1)) - b*(m + p + S(2))*cos(e + f*x)), x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))/(f*g*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2313(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b*(m + p)), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**m*(a*sin(e + f*x) + b), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))/(b*f*(m + p)), x)


def replacement2314(a, b, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b*(m + p)), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**m*(a*cos(e + f*x) + b), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))/(b*f*(m + p)), x)


def replacement2315(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(a**S(2) - b**S(2))*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m*(a**S(2)*(p + S(2)) + a*b*(m + p + S(3))*sin(e + f*x) - b**S(2)*(m + p + S(2))), x), x) + Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(-a*sin(e + f*x) + b)/(f*g*(a**S(2) - b**S(2))*(p + S(1))), x)


def replacement2316(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(a**S(2) - b**S(2))*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m*(a**S(2)*(p + S(2)) + a*b*(m + p + S(3))*cos(e + f*x) - b**S(2)*(m + p + S(2))), x), x) - Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(-a*cos(e + f*x) + b)/(f*g*(a**S(2) - b**S(2))*(p + S(1))), x)


def replacement2317(a, b, e, f, g, x):
    return Dist(S(2)*sqrt(S(2))*sqrt(g*cos(e + f*x))*sqrt((a + b*sin(e + f*x))/((S(1) - sin(e + f*x))*(a - b)))/(f*g*sqrt((sin(e + f*x) + cos(e + f*x) + S(1))/(-sin(e + f*x) + cos(e + f*x) + S(1)))*sqrt(a + b*sin(e + f*x))), Subst(Int(S(1)/sqrt(x**S(4)*(a + b)/(a - b) + S(1)), x), x, sqrt((sin(e + f*x) + cos(e + f*x) + S(1))/(-sin(e + f*x) + cos(e + f*x) + S(1)))), x)


def replacement2318(a, b, e, f, g, x):
    return Dist(-S(2)*sqrt(S(2))*sqrt(g*sin(e + f*x))*sqrt((a + b*cos(e + f*x))/((S(1) - cos(e + f*x))*(a - b)))/(f*g*sqrt((sin(e + f*x) + cos(e + f*x) + S(1))/(sin(e + f*x) - cos(e + f*x) + S(1)))*sqrt(a + b*cos(e + f*x))), Subst(Int(S(1)/sqrt(x**S(4)*(a + b)/(a - b) + S(1)), x), x, sqrt((sin(e + f*x) + cos(e + f*x) + S(1))/(sin(e + f*x) - cos(e + f*x) + S(1)))), x)


def replacement2319(a, b, e, f, g, m, p, x):
    return Simp(g*(g*cos(e + f*x))**(p + S(-1))*(-(S(1) - sin(e + f*x))*(a - b)/((a + b)*(sin(e + f*x) + S(1))))**(m/S(2))*(S(1) - sin(e + f*x))*(a + b*sin(e + f*x))**(m + S(1))*Hypergeometric2F1(m + S(1), m/S(2) + S(1), m + S(2), S(2)*(a + b*sin(e + f*x))/((a + b)*(sin(e + f*x) + S(1))))/(f*(a + b)*(m + S(1))), x)


def replacement2320(a, b, e, f, g, m, p, x):
    return -Simp(g*(g*sin(e + f*x))**(p + S(-1))*(-(S(1) - cos(e + f*x))*(a - b)/((a + b)*(cos(e + f*x) + S(1))))**(m/S(2))*(S(1) - cos(e + f*x))*(a + b*cos(e + f*x))**(m + S(1))*Hypergeometric2F1(m + S(1), m/S(2) + S(1), m + S(2), S(2)*(a + b*cos(e + f*x))/((a + b)*(cos(e + f*x) + S(1))))/(f*(a + b)*(m + S(1))), x)


def replacement2321(a, b, e, f, g, m, p, x):
    return Dist(a/(g**S(2)*(a - b)), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m/(S(1) - sin(e + f*x)), x), x) + Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))/(f*g*(a - b)*(p + S(1))), x)


def replacement2322(a, b, e, f, g, m, p, x):
    return Dist(a/(g**S(2)*(a - b)), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m/(S(1) - cos(e + f*x)), x), x) - Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))/(f*g*(a - b)*(p + S(1))), x)


def replacement2323(a, b, e, f, g, m, p, x):
    return Dist(a/(g**S(2)*(a - b)), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m/(S(1) - sin(e + f*x)), x), x) - Dist(b*(m + p + S(2))/(g**S(2)*(a - b)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m, x), x) + Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))/(f*g*(a - b)*(p + S(1))), x)


def replacement2324(a, b, e, f, g, m, p, x):
    return Dist(a/(g**S(2)*(a - b)), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m/(S(1) - cos(e + f*x)), x), x) - Dist(b*(m + p + S(2))/(g**S(2)*(a - b)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m, x), x) - Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))/(f*g*(a - b)*(p + S(1))), x)


def With2325(a, b, e, f, g, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return -Dist(a*g/(S(2)*b), Int(S(1)/(sqrt(g*cos(e + f*x))*(-b*cos(e + f*x) + q)), x), x) + Dist(a*g/(S(2)*b), Int(S(1)/(sqrt(g*cos(e + f*x))*(b*cos(e + f*x) + q)), x), x) + Dist(b*g/f, Subst(Int(sqrt(x)/(b**S(2)*x**S(2) + g**S(2)*(a**S(2) - b**S(2))), x), x, g*cos(e + f*x)), x)


def With2326(a, b, e, f, g, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return -Dist(a*g/(S(2)*b), Int(S(1)/(sqrt(g*sin(e + f*x))*(-b*sin(e + f*x) + q)), x), x) + Dist(a*g/(S(2)*b), Int(S(1)/(sqrt(g*sin(e + f*x))*(b*sin(e + f*x) + q)), x), x) - Dist(b*g/f, Subst(Int(sqrt(x)/(b**S(2)*x**S(2) + g**S(2)*(a**S(2) - b**S(2))), x), x, g*sin(e + f*x)), x)


def With2327(a, b, e, f, g, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return -Dist(a/(S(2)*q), Int(S(1)/(sqrt(g*cos(e + f*x))*(-b*cos(e + f*x) + q)), x), x) - Dist(a/(S(2)*q), Int(S(1)/(sqrt(g*cos(e + f*x))*(b*cos(e + f*x) + q)), x), x) + Dist(b*g/f, Subst(Int(S(1)/(sqrt(x)*(b**S(2)*x**S(2) + g**S(2)*(a**S(2) - b**S(2)))), x), x, g*cos(e + f*x)), x)


def With2328(a, b, e, f, g, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return -Dist(a/(S(2)*q), Int(S(1)/(sqrt(g*sin(e + f*x))*(-b*sin(e + f*x) + q)), x), x) - Dist(a/(S(2)*q), Int(S(1)/(sqrt(g*sin(e + f*x))*(b*sin(e + f*x) + q)), x), x) - Dist(b*g/f, Subst(Int(S(1)/(sqrt(x)*(b**S(2)*x**S(2) + g**S(2)*(a**S(2) - b**S(2)))), x), x, g*sin(e + f*x)), x)


def replacement2329(a, b, e, f, g, m, p, x):
    return Simp(g*(g*cos(e + f*x))**(p + S(-1))*(b*(sin(e + f*x) + S(1))/(a + b*sin(e + f*x)))**(S(1)/2 - p/S(2))*(-b*(S(1) - sin(e + f*x))/(a + b*sin(e + f*x)))**(S(1)/2 - p/S(2))*(a + b*sin(e + f*x))**(m + S(1))*AppellF1(-m - p, S(1)/2 - p/S(2), S(1)/2 - p/S(2), -m - p + S(1), (a + b)/(a + b*sin(e + f*x)), (a - b)/(a + b*sin(e + f*x)))/(b*f*(m + p)), x)


def replacement2330(a, b, e, f, g, m, p, x):
    return -Simp(g*(g*sin(e + f*x))**(p + S(-1))*(b*(cos(e + f*x) + S(1))/(a + b*cos(e + f*x)))**(S(1)/2 - p/S(2))*(-b*(S(1) - cos(e + f*x))/(a + b*cos(e + f*x)))**(S(1)/2 - p/S(2))*(a + b*cos(e + f*x))**(m + S(1))*AppellF1(-m - p, S(1)/2 - p/S(2), S(1)/2 - p/S(2), -m - p + S(1), (a + b)/(a + b*cos(e + f*x)), (a - b)/(a + b*cos(e + f*x)))/(b*f*(m + p)), x)


def replacement2331(a, b, e, f, g, m, p, x):
    return Dist(g*(g*cos(e + f*x))**(p + S(-1))*(S(1) - (a + b*sin(e + f*x))/(a - b))**(S(1)/2 - p/S(2))*(S(1) - (a + b*sin(e + f*x))/(a + b))**(S(1)/2 - p/S(2))/f, Subst(Int((a + b*x)**m*(-b*x/(a - b) - b/(a - b))**(p/S(2) + S(-1)/2)*(-b*x/(a + b) + b/(a + b))**(p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2332(a, b, e, f, g, m, p, x):
    return -Dist(g*(g*sin(e + f*x))**(p + S(-1))*(S(1) - (a + b*cos(e + f*x))/(a - b))**(S(1)/2 - p/S(2))*(S(1) - (a + b*cos(e + f*x))/(a + b))**(S(1)/2 - p/S(2))/f, Subst(Int((a + b*x)**m*(-b*x/(a - b) - b/(a - b))**(p/S(2) + S(-1)/2)*(-b*x/(a + b) + b/(a + b))**(p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2333(a, b, e, f, g, m, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/cos(e + f*x))**FracPart(p)*(g*cos(e + f*x))**FracPart(p), Int((g*cos(e + f*x))**(-p)*(a + b*sin(e + f*x))**m, x), x)


def replacement2334(a, b, e, f, g, m, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/sin(e + f*x))**FracPart(p)*(g*sin(e + f*x))**FracPart(p), Int((g*sin(e + f*x))**(-p)*(a + b*cos(e + f*x))**m, x), x)


def replacement2335(a, b, e, f, g, p, x):
    return Dist(S(1)/a, Int((g*tan(e + f*x))**p/cos(e + f*x)**S(2), x), x) - Dist(S(1)/(b*g), Int((g*tan(e + f*x))**(p + S(1))/cos(e + f*x), x), x)


def replacement2336(a, b, e, f, g, p, x):
    return Dist(S(1)/a, Int((g/tan(e + f*x))**p/sin(e + f*x)**S(2), x), x) - Dist(S(1)/(b*g), Int((g/tan(e + f*x))**(p + S(1))/sin(e + f*x), x), x)


def replacement2337(a, b, e, f, m, p, x):
    return Dist(S(1)/f, Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2338(a, b, e, f, m, p, x):
    return -Dist(S(1)/f, Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2339(a, b, e, f, m, p, x):
    return Dist(a**p, Int((a - b*sin(e + f*x))**(-m)*sin(e + f*x)**p, x), x)


def replacement2340(a, b, e, f, m, p, x):
    return Dist(a**p, Int((a - b*cos(e + f*x))**(-m)*cos(e + f*x)**p, x), x)


def replacement2341(a, b, e, f, m, p, x):
    return Dist(a**p, Int(ExpandIntegrand((a - b*sin(e + f*x))**(-p/S(2))*(a + b*sin(e + f*x))**(m - p/S(2))*sin(e + f*x)**p, x), x), x)


def replacement2342(a, b, e, f, m, p, x):
    return Dist(a**p, Int(ExpandIntegrand((a - b*cos(e + f*x))**(-p/S(2))*(a + b*cos(e + f*x))**(m - p/S(2))*cos(e + f*x)**p, x), x), x)


def replacement2343(a, b, e, f, g, m, p, x):
    return Int(ExpandIntegrand((g*tan(e + f*x))**p, (a + b*sin(e + f*x))**m, x), x)


def replacement2344(a, b, e, f, g, m, p, x):
    return Int(ExpandIntegrand((g/tan(e + f*x))**p, (a + b*cos(e + f*x))**m, x), x)


def replacement2345(a, b, e, f, g, m, p, x):
    return Dist(a**(S(2)*m), Int(ExpandIntegrand((g*tan(e + f*x))**p*(S(1)/cos(e + f*x))**(-m), (a/cos(e + f*x) - b*tan(e + f*x))**(-m), x), x), x)


def replacement2346(a, b, e, f, g, m, p, x):
    return Dist(a**(S(2)*m), Int(ExpandIntegrand((g/tan(e + f*x))**p*(S(1)/sin(e + f*x))**(-m), (a/sin(e + f*x) - b/tan(e + f*x))**(-m), x), x), x)


def replacement2347(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(-1))), Int((a + b*sin(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(-1))*sin(e + f*x))/cos(e + f*x)**S(2), x), x) + Simp(b*(a + b*sin(e + f*x))**m/(a*f*(S(2)*m + S(-1))*cos(e + f*x)), x)


def replacement2348(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(-1))), Int((a + b*cos(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(-1))*cos(e + f*x))/sin(e + f*x)**S(2), x), x) - Simp(b*(a + b*cos(e + f*x))**m/(a*f*(S(2)*m + S(-1))*sin(e + f*x)), x)


def replacement2349(a, b, e, f, m, x):
    return Dist(S(1)/(b*m), Int((a + b*sin(e + f*x))**m*(a*sin(e + f*x) + b*(m + S(1)))/cos(e + f*x)**S(2), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))/(b*f*m*cos(e + f*x)), x)


def replacement2350(a, b, e, f, m, x):
    return Dist(S(1)/(b*m), Int((a + b*cos(e + f*x))**m*(a*cos(e + f*x) + b*(m + S(1)))/sin(e + f*x)**S(2), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))/(b*f*m*sin(e + f*x)), x)


def replacement2351(a, b, e, f, m, x):
    return -Int((S(1) - S(2)*sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m/cos(e + f*x)**S(4), x) + Int((a + b*sin(e + f*x))**m, x)


def replacement2352(a, b, e, f, m, x):
    return -Int((S(1) - S(2)*cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m/sin(e + f*x)**S(4), x) + Int((a + b*cos(e + f*x))**m, x)


def replacement2353(a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b*sin(e + f*x))**(m + S(1))*(-a*(m + S(1))*sin(e + f*x) + b*m)/sin(e + f*x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))/(a*f*tan(e + f*x)), x)


def replacement2354(a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b*cos(e + f*x))**(m + S(1))*(-a*(m + S(1))*cos(e + f*x) + b*m)/cos(e + f*x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*tan(e + f*x)/(a*f), x)


def replacement2355(a, b, e, f, m, x):
    return Dist(S(1)/a, Int((a + b*sin(e + f*x))**m*(-a*(m + S(1))*sin(e + f*x) + b*m)/sin(e + f*x), x), x) - Simp((a + b*sin(e + f*x))**m/(f*tan(e + f*x)), x)


def replacement2356(a, b, e, f, m, x):
    return Dist(S(1)/a, Int((a + b*cos(e + f*x))**m*(-a*(m + S(1))*cos(e + f*x) + b*m)/cos(e + f*x), x), x) + Simp((a + b*cos(e + f*x))**m*tan(e + f*x)/f, x)


def replacement2357(a, b, e, f, m, x):
    return Dist(a**(S(-2)), Int((a + b*sin(e + f*x))**(m + S(2))*(sin(e + f*x)**S(2) + S(1))/sin(e + f*x)**S(4), x), x) + Dist(-S(2)/(a*b), Int((a + b*sin(e + f*x))**(m + S(2))/sin(e + f*x)**S(3), x), x)


def replacement2358(a, b, e, f, m, x):
    return Dist(a**(S(-2)), Int((a + b*cos(e + f*x))**(m + S(2))*(cos(e + f*x)**S(2) + S(1))/cos(e + f*x)**S(4), x), x) + Dist(-S(2)/(a*b), Int((a + b*cos(e + f*x))**(m + S(2))/cos(e + f*x)**S(3), x), x)


def replacement2359(a, b, e, f, m, x):
    return Int((S(1) - S(2)*sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m/sin(e + f*x)**S(4), x) + Int((a + b*sin(e + f*x))**m, x)


def replacement2360(a, b, e, f, m, x):
    return Int((S(1) - S(2)*cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m/cos(e + f*x)**S(4), x) + Int((a + b*cos(e + f*x))**m, x)


def replacement2361(a, b, e, f, m, p, x):
    return Dist(sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))/(b*f*cos(e + f*x)), Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2362(a, b, e, f, m, p, x):
    return -Dist(sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))/(b*f*sin(e + f*x)), Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2363(a, b, e, f, g, m, p, x):
    return Dist((b*sin(e + f*x))**(-p + S(-1))*(g*tan(e + f*x))**(p + S(1))*(a - b*sin(e + f*x))**(p/S(2) + S(1)/2)*(a + b*sin(e + f*x))**(p/S(2) + S(1)/2)/(f*g), Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2364(a, b, e, f, g, m, p, x):
    return -Dist((b*cos(e + f*x))**(-p + S(-1))*(g/tan(e + f*x))**(p + S(1))*(a - b*cos(e + f*x))**(p/S(2) + S(1)/2)*(a + b*cos(e + f*x))**(p/S(2) + S(1)/2)/(f*g), Subst(Int(x**p*(a - x)**(-p/S(2) + S(-1)/2)*(a + x)**(m - p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2365(a, b, e, f, m, p, x):
    return Dist(S(1)/f, Subst(Int(x**p*(a + x)**m*(b**S(2) - x**S(2))**(-p/S(2) + S(-1)/2), x), x, b*sin(e + f*x)), x)


def replacement2366(a, b, e, f, m, p, x):
    return -Dist(S(1)/f, Subst(Int(x**p*(a + x)**m*(b**S(2) - x**S(2))**(-p/S(2) + S(-1)/2), x), x, b*cos(e + f*x)), x)


def replacement2367(a, b, e, f, g, m, p, x):
    return Int(ExpandIntegrand((g*tan(e + f*x))**p, (a + b*sin(e + f*x))**m, x), x)


def replacement2368(a, b, e, f, g, m, p, x):
    return Int(ExpandIntegrand((g/tan(e + f*x))**p, (a + b*cos(e + f*x))**m, x), x)


def replacement2369(a, b, e, f, m, x):
    return Int((S(1) - sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m/sin(e + f*x)**S(2), x)


def replacement2370(a, b, e, f, m, x):
    return Int((S(1) - cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m/cos(e + f*x)**S(2), x)


def replacement2371(a, b, e, f, m, x):
    return -Dist(S(1)/(S(3)*a**S(2)*b*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(S(6)*a**S(2) + a*b*(m + S(1))*sin(e + f*x) - b**S(2)*(m + S(-2))*(m + S(-1)) - (S(3)*a**S(2) - b**S(2)*m*(m + S(-2)))*sin(e + f*x)**S(2), x)/sin(e + f*x)**S(3), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(S(3)*a*f*sin(e + f*x)**S(3)), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(S(3)*a**S(2) + b**S(2)*(m + S(-2)))*cos(e + f*x)/(S(3)*a**S(2)*b*f*(m + S(1))*sin(e + f*x)**S(2)), x)


def replacement2372(a, b, e, f, m, x):
    return -Dist(S(1)/(S(3)*a**S(2)*b*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(S(6)*a**S(2) + a*b*(m + S(1))*cos(e + f*x) - b**S(2)*(m + S(-2))*(m + S(-1)) - (S(3)*a**S(2) - b**S(2)*m*(m + S(-2)))*cos(e + f*x)**S(2), x)/cos(e + f*x)**S(3), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(S(3)*a*f*cos(e + f*x)**S(3)), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(S(3)*a**S(2) + b**S(2)*(m + S(-2)))*sin(e + f*x)/(S(3)*a**S(2)*b*f*(m + S(1))*cos(e + f*x)**S(2)), x)


def replacement2373(a, b, e, f, m, x):
    return -Dist(S(1)/(S(6)*a**S(2)), Int((a + b*sin(e + f*x))**m*Simp(S(8)*a**S(2) + a*b*m*sin(e + f*x) - b**S(2)*(m + S(-2))*(m + S(-1)) - (S(6)*a**S(2) - b**S(2)*m*(m + S(-2)))*sin(e + f*x)**S(2), x)/sin(e + f*x)**S(2), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(S(3)*a*f*sin(e + f*x)**S(3)), x) - Simp(b*(a + b*sin(e + f*x))**(m + S(1))*(m + S(-2))*cos(e + f*x)/(S(6)*a**S(2)*f*sin(e + f*x)**S(2)), x)


def replacement2374(a, b, e, f, m, x):
    return -Dist(S(1)/(S(6)*a**S(2)), Int((a + b*cos(e + f*x))**m*Simp(S(8)*a**S(2) + a*b*m*cos(e + f*x) - b**S(2)*(m + S(-2))*(m + S(-1)) - (S(6)*a**S(2) - b**S(2)*m*(m + S(-2)))*cos(e + f*x)**S(2), x)/cos(e + f*x)**S(2), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(S(3)*a*f*cos(e + f*x)**S(3)), x) + Simp(b*(a + b*cos(e + f*x))**(m + S(1))*(m + S(-2))*sin(e + f*x)/(S(6)*a**S(2)*f*cos(e + f*x)**S(2)), x)


def replacement2375(a, b, e, f, m, x):
    return Dist(S(1)/(S(20)*a**S(2)*b**S(2)*m*(m + S(-1))), Int((a + b*sin(e + f*x))**m*Simp(S(60)*a**S(4) - S(44)*a**S(2)*b**S(2)*m*(m + S(-1)) + a*b*m*(S(20)*a**S(2) - b**S(2)*m*(m + S(-1)))*sin(e + f*x) + b**S(4)*m*(m + S(-4))*(m + S(-3))*(m + S(-1)) - (S(40)*a**S(4) - S(20)*a**S(2)*b**S(2)*(m + S(-1))*(S(2)*m + S(1)) + b**S(4)*m*(m + S(-4))*(m + S(-2))*(m + S(-1)))*sin(e + f*x)**S(2), x)/sin(e + f*x)**S(4), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(S(5)*a*f*sin(e + f*x)**S(5)), x) + Simp((a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*m*sin(e + f*x)**S(2)), x) - Simp(b*(a + b*sin(e + f*x))**(m + S(1))*(m + S(-4))*cos(e + f*x)/(S(20)*a**S(2)*f*sin(e + f*x)**S(4)), x) + Simp(a*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b**S(2)*f*m*(m + S(-1))*sin(e + f*x)**S(3)), x)


def replacement2376(a, b, e, f, m, x):
    return Dist(S(1)/(S(20)*a**S(2)*b**S(2)*m*(m + S(-1))), Int((a + b*cos(e + f*x))**m*Simp(S(60)*a**S(4) - S(44)*a**S(2)*b**S(2)*m*(m + S(-1)) + a*b*m*(S(20)*a**S(2) - b**S(2)*m*(m + S(-1)))*cos(e + f*x) + b**S(4)*m*(m + S(-4))*(m + S(-3))*(m + S(-1)) - (S(40)*a**S(4) - S(20)*a**S(2)*b**S(2)*(m + S(-1))*(S(2)*m + S(1)) + b**S(4)*m*(m + S(-4))*(m + S(-2))*(m + S(-1)))*cos(e + f*x)**S(2), x)/cos(e + f*x)**S(4), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(S(5)*a*f*cos(e + f*x)**S(5)), x) - Simp((a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*m*cos(e + f*x)**S(2)), x) + Simp(b*(a + b*cos(e + f*x))**(m + S(1))*(m + S(-4))*sin(e + f*x)/(S(20)*a**S(2)*f*cos(e + f*x)**S(4)), x) - Simp(a*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b**S(2)*f*m*(m + S(-1))*cos(e + f*x)**S(3)), x)


def replacement2377(a, b, e, f, g, p, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((g*tan(e + f*x))**p/sin(e + f*x)**S(2), x), x) - Dist(a**S(2)*g**S(2)/(a**S(2) - b**S(2)), Int((g*tan(e + f*x))**(p + S(-2))/(a + b*sin(e + f*x)), x), x) - Dist(b*g/(a**S(2) - b**S(2)), Int((g*tan(e + f*x))**(p + S(-1))/cos(e + f*x), x), x)


def replacement2378(a, b, e, f, g, p, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((g/tan(e + f*x))**p/cos(e + f*x)**S(2), x), x) - Dist(a**S(2)*g**S(2)/(a**S(2) - b**S(2)), Int((g/tan(e + f*x))**(p + S(-2))/(a + b*cos(e + f*x)), x), x) - Dist(b*g/(a**S(2) - b**S(2)), Int((g/tan(e + f*x))**(p + S(-1))/sin(e + f*x), x), x)


def replacement2379(a, b, e, f, g, p, x):
    return Dist(S(1)/a, Int((g*tan(e + f*x))**p/cos(e + f*x)**S(2), x), x) - Dist(b/(a**S(2)*g), Int((g*tan(e + f*x))**(p + S(1))/cos(e + f*x), x), x) - Dist((a**S(2) - b**S(2))/(a**S(2)*g**S(2)), Int((g*tan(e + f*x))**(p + S(2))/(a + b*sin(e + f*x)), x), x)


def replacement2380(a, b, e, f, g, p, x):
    return Dist(S(1)/a, Int((g/tan(e + f*x))**p/sin(e + f*x)**S(2), x), x) - Dist(b/(a**S(2)*g), Int((g/tan(e + f*x))**(p + S(1))/sin(e + f*x), x), x) - Dist((a**S(2) - b**S(2))/(a**S(2)*g**S(2)), Int((g/tan(e + f*x))**(p + S(2))/(a + b*cos(e + f*x)), x), x)


def replacement2381(a, b, e, f, g, x):
    return Dist(sqrt(g*tan(e + f*x))*sqrt(cos(e + f*x))/sqrt(sin(e + f*x)), Int(sqrt(sin(e + f*x))/((a + b*sin(e + f*x))*sqrt(cos(e + f*x))), x), x)


def replacement2382(a, b, e, f, g, x):
    return Dist(sqrt(g/tan(e + f*x))*sqrt(sin(e + f*x))/sqrt(cos(e + f*x)), Int(sqrt(cos(e + f*x))/((a + b*cos(e + f*x))*sqrt(sin(e + f*x))), x), x)


def replacement2383(a, b, e, f, g, x):
    return Dist(sqrt(sin(e + f*x))/(sqrt(g*tan(e + f*x))*sqrt(cos(e + f*x))), Int(sqrt(cos(e + f*x))/((a + b*sin(e + f*x))*sqrt(sin(e + f*x))), x), x)


def replacement2384(a, b, e, f, g, x):
    return Dist(sqrt(cos(e + f*x))/(sqrt(g/tan(e + f*x))*sqrt(sin(e + f*x))), Int(sqrt(sin(e + f*x))/((a + b*cos(e + f*x))*sqrt(cos(e + f*x))), x), x)


def replacement2385(a, b, e, f, m, p, x):
    return Int(ExpandIntegrand((S(1) - sin(e + f*x)**S(2))**(-p/S(2))*(a + b*sin(e + f*x))**m*sin(e + f*x)**p, x), x)


def replacement2386(a, b, e, f, m, p, x):
    return Int(ExpandIntegrand((S(1) - cos(e + f*x)**S(2))**(-p/S(2))*(a + b*cos(e + f*x))**m*cos(e + f*x)**p, x), x)


def replacement2387(a, b, e, f, g, m, p, x):
    return Int((g*tan(e + f*x))**p*(a + b*sin(e + f*x))**m, x)


def replacement2388(a, b, e, f, g, m, p, x):
    return Int((g/tan(e + f*x))**p*(a + b*cos(e + f*x))**m, x)


def replacement2389(a, b, e, f, g, m, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/tan(e + f*x))**FracPart(p)*(g*tan(e + f*x))**FracPart(p), Int((g*tan(e + f*x))**(-p)*(a + b*sin(e + f*x))**m, x), x)


def replacement2390(a, b, e, f, g, m, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/tan(e + f*x))**FracPart(p)*(g*tan(e + f*x))**FracPart(p), Int((g/tan(e + f*x))**(-p)*(a + b*cos(e + f*x))**m, x), x)


def replacement2391(a, b, c, d, e, f, x):
    return Simp(x*(S(2)*a*c + b*d)/S(2), x) - Simp((a*d + b*c)*cos(e + f*x)/f, x) - Simp(b*d*sin(e + f*x)*cos(e + f*x)/(S(2)*f), x)


def replacement2392(a, b, c, d, e, f, x):
    return Simp(x*(S(2)*a*c + b*d)/S(2), x) + Simp((a*d + b*c)*sin(e + f*x)/f, x) + Simp(b*d*sin(e + f*x)*cos(e + f*x)/(S(2)*f), x)


def replacement2393(a, b, c, d, e, f, x):
    return -Dist((-a*d + b*c)/d, Int(S(1)/(c + d*sin(e + f*x)), x), x) + Simp(b*x/d, x)


def replacement2394(a, b, c, d, e, f, x):
    return -Dist((-a*d + b*c)/d, Int(S(1)/(c + d*cos(e + f*x)), x), x) + Simp(b*x/d, x)


def replacement2395(a, b, c, d, e, f, m, n, x):
    return Dist(a**m*c**m, Int((c + d*sin(e + f*x))**(-m + n)*cos(e + f*x)**(S(2)*m), x), x)


def replacement2396(a, b, c, d, e, f, m, n, x):
    return Dist(a**m*c**m, Int((c + d*cos(e + f*x))**(-m + n)*sin(e + f*x)**(S(2)*m), x), x)


def replacement2397(a, b, c, d, e, f, x):
    return Dist(a*c*cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), Int(cos(e + f*x)/(c + d*sin(e + f*x)), x), x)


def replacement2398(a, b, c, d, e, f, x):
    return Dist(a*c*sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), Int(sin(e + f*x)/(c + d*cos(e + f*x)), x), x)


def replacement2399(a, b, c, d, e, f, n, x):
    return Simp(-S(2)*b*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(S(2)*n + S(1))), x)


def replacement2400(a, b, c, d, e, f, n, x):
    return Simp(S(2)*b*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement2401(a, b, c, d, e, f, m, n, x):
    return -Dist(b*(S(2)*m + S(-1))/(d*(S(2)*n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1)), x), x) + Simp(-S(2)*b*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(S(2)*n + S(1))), x)


def replacement2402(a, b, c, d, e, f, m, n, x):
    return -Dist(b*(S(2)*m + S(-1))/(d*(S(2)*n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1)), x), x) + Simp(S(2)*b*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(S(2)*n + S(1))), x)


def replacement2403(a, b, c, d, e, f, m, n, x):
    return Dist(a*(S(2)*m + S(-1))/(m + n), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n, x), x) - Simp(b*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(m + n)), x)


def replacement2404(a, b, c, d, e, f, m, n, x):
    return Dist(a*(S(2)*m + S(-1))/(m + n), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n, x), x) + Simp(b*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(m + n)), x)


def replacement2405(a, b, c, d, e, f, x):
    return Dist(cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), Int(S(1)/cos(e + f*x), x), x)


def replacement2406(a, b, c, d, e, f, x):
    return Dist(sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), Int(S(1)/sin(e + f*x), x), x)


def replacement2407(a, b, c, d, e, f, m, n, x):
    return Simp(b*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2408(a, b, c, d, e, f, m, n, x):
    return -Simp(b*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2409(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Simp(b*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2410(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) - Simp(b*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2411(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Simp(b*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2412(a, b, c, d, e, f, m, n, x):
    return Dist((m + n + S(1))/(a*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) - Simp(b*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2413(a, b, c, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*(a + b*sin(e + f*x))**FracPart(m)*(c + d*sin(e + f*x))**FracPart(m)*cos(e + f*x)**(-S(2)*FracPart(m)), Int((c + d*sin(e + f*x))**(-m + n)*cos(e + f*x)**(S(2)*m), x), x)


def replacement2414(a, b, c, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*(a + b*cos(e + f*x))**FracPart(m)*(c + d*cos(e + f*x))**FracPart(m)*sin(e + f*x)**(-S(2)*FracPart(m)), Int((c + d*cos(e + f*x))**(-m + n)*sin(e + f*x)**(S(2)*m), x), x)


def replacement2415(a, b, c, d, e, f, x):
    return Dist(S(1)/d, Int(Simp(a**S(2)*d - b*(-S(2)*a*d + b*c)*sin(e + f*x), x)/(c + d*sin(e + f*x)), x), x) - Simp(b**S(2)*cos(e + f*x)/(d*f), x)


def replacement2416(a, b, c, d, e, f, x):
    return Dist(S(1)/d, Int(Simp(a**S(2)*d - b*(-S(2)*a*d + b*c)*cos(e + f*x), x)/(c + d*cos(e + f*x)), x), x) + Simp(b**S(2)*sin(e + f*x)/(d*f), x)


def replacement2417(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(a + b*sin(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d*sin(e + f*x)), x), x)


def replacement2418(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(a + b*cos(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d*cos(e + f*x)), x), x)


def replacement2419(b, c, d, e, f, m, x):
    return Dist(c, Int((b*sin(e + f*x))**m, x), x) + Dist(d/b, Int((b*sin(e + f*x))**(m + S(1)), x), x)


def replacement2420(b, c, d, e, f, m, x):
    return Dist(c, Int((b*cos(e + f*x))**m, x), x) + Dist(d/b, Int((b*cos(e + f*x))**(m + S(1)), x), x)


def replacement2421(a, b, c, d, e, f, m, x):
    return -Simp(d*(a + b*sin(e + f*x))**m*cos(e + f*x)/(f*(m + S(1))), x)


def replacement2422(a, b, c, d, e, f, m, x):
    return Simp(d*(a + b*cos(e + f*x))**m*sin(e + f*x)/(f*(m + S(1))), x)


def replacement2423(a, b, c, d, e, f, m, x):
    return Dist((a*d*m + b*c*(m + S(1)))/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1)), x), x) + Simp((a + b*sin(e + f*x))**m*(-a*d + b*c)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2424(a, b, c, d, e, f, m, x):
    return Dist((a*d*m + b*c*(m + S(1)))/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1)), x), x) - Simp((a + b*cos(e + f*x))**m*(-a*d + b*c)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2425(a, b, c, d, e, f, m, x):
    return Dist((a*d*m + b*c*(m + S(1)))/(b*(m + S(1))), Int((a + b*sin(e + f*x))**m, x), x) - Simp(d*(a + b*sin(e + f*x))**m*cos(e + f*x)/(f*(m + S(1))), x)


def replacement2426(a, b, c, d, e, f, m, x):
    return Dist((a*d*m + b*c*(m + S(1)))/(b*(m + S(1))), Int((a + b*cos(e + f*x))**m, x), x) + Simp(d*(a + b*cos(e + f*x))**m*sin(e + f*x)/(f*(m + S(1))), x)


def replacement2427(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(a + b*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/sqrt(a + b*sin(e + f*x)), x), x)


def replacement2428(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(a + b*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/sqrt(a + b*cos(e + f*x)), x), x)


def replacement2429(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b*sin(e + f*x))**(m + S(-1))*Simp(a*c*(m + S(1)) + b*d*m + (a*d*m + b*c*(m + S(1)))*sin(e + f*x), x), x), x) - Simp(d*(a + b*sin(e + f*x))**m*cos(e + f*x)/(f*(m + S(1))), x)


def replacement2430(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(m + S(1)), Int((a + b*cos(e + f*x))**(m + S(-1))*Simp(a*c*(m + S(1)) + b*d*m + (a*d*m + b*c*(m + S(1)))*cos(e + f*x), x), x), x) + Simp(d*(a + b*cos(e + f*x))**m*sin(e + f*x)/(f*(m + S(1))), x)


def replacement2431(a, b, c, d, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp((m + S(1))*(a*c - b*d) - (m + S(2))*(-a*d + b*c)*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(-a*d + b*c)*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2432(a, b, c, d, e, f, m, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp((m + S(1))*(a*c - b*d) - (m + S(2))*(-a*d + b*c)*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(-a*d + b*c)*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2433(a, b, c, d, e, f, m, x):
    return Dist(c*cos(e + f*x)/(f*sqrt(S(1) - sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), Subst(Int(sqrt(S(1) + d*x/c)*(a + b*x)**m/sqrt(S(1) - d*x/c), x), x, sin(e + f*x)), x)


def replacement2434(a, b, c, d, e, f, m, x):
    return -Dist(c*sin(e + f*x)/(f*sqrt(S(1) - cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), Subst(Int(sqrt(S(1) + d*x/c)*(a + b*x)**m/sqrt(S(1) - d*x/c), x), x, cos(e + f*x)), x)


def replacement2435(a, b, c, d, e, f, m, x):
    return Dist(d/b, Int((a + b*sin(e + f*x))**(m + S(1)), x), x) + Dist((-a*d + b*c)/b, Int((a + b*sin(e + f*x))**m, x), x)


def replacement2436(a, b, c, d, e, f, m, x):
    return Dist(d/b, Int((a + b*cos(e + f*x))**(m + S(1)), x), x) + Dist((-a*d + b*c)/b, Int((a + b*cos(e + f*x))**m, x), x)


def replacement2437(a, b, d, e, f, m, n, x):
    return Int(ExpandTrig((d*sin(e + f*x))**n*(a + b*sin(e + f*x))**m, x), x)


def replacement2438(a, b, d, e, f, m, n, x):
    return Int(ExpandTrig((d*cos(e + f*x))**n*(a + b*cos(e + f*x))**m, x), x)


def replacement2439(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(1))*sin(e + f*x)), x), x) + Simp(b*(a + b*sin(e + f*x))**m*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2440(a, b, e, f, m, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + S(1))*cos(e + f*x)), x), x) - Simp(b*(a + b*cos(e + f*x))**m*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2441(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*sin(e + f*x))**m*(-a*sin(e + f*x) + b*(m + S(1))), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2442(a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*cos(e + f*x))**m*(-a*cos(e + f*x) + b*(m + S(1))), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2443(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(a*c*d*(m + S(-1)) + b*(c**S(2)*(m + S(1)) + d**S(2)) + d*(a*d*(m + S(-1)) + b*c*(m + S(2)))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))*(-a*d + b*c)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2444(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(a*c*d*(m + S(-1)) + b*(c**S(2)*(m + S(1)) + d**S(2)) + d*(a*d*(m + S(-1)) + b*c*(m + S(2)))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))*(-a*d + b*c)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2445(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*sin(e + f*x))**m*Simp(b*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))) - d*(a*d - S(2)*b*c*(m + S(2)))*sin(e + f*x), x), x), x) - Simp(d**S(2)*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2446(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*cos(e + f*x))**m*Simp(b*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))) - d*(a*d - S(2)*b*c*(m + S(2)))*cos(e + f*x), x), x), x) + Simp(d**S(2)*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2447(a, b, c, d, e, f, m, n, x):
    return Dist(b**S(2)/(d*(n + S(1))*(a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*Simp(a*c*(m + S(-2)) - b*d*(m - S(2)*n + S(-4)) - (-a*d*(m + S(2)*n + S(1)) + b*c*(m + S(-1)))*sin(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*(-a*d + b*c)*cos(e + f*x)/(d*f*(n + S(1))*(a*d + b*c)), x)


def replacement2448(a, b, c, d, e, f, m, n, x):
    return Dist(b**S(2)/(d*(n + S(1))*(a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*Simp(a*c*(m + S(-2)) - b*d*(m - S(2)*n + S(-4)) - (-a*d*(m + S(2)*n + S(1)) + b*c*(m + S(-1)))*cos(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*(-a*d + b*c)*sin(e + f*x)/(d*f*(n + S(1))*(a*d + b*c)), x)


def replacement2449(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**n*Simp(a**S(2)*d*(m + n) + a*b*c*(m + S(-2)) + b**S(2)*d*(n + S(1)) - b*(-a*d*(S(3)*m + S(2)*n + S(-2)) + b*c*(m + S(-1)))*sin(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n)), x)


def replacement2450(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**n*Simp(a**S(2)*d*(m + n) + a*b*c*(m + S(-2)) + b**S(2)*d*(n + S(1)) - b*(-a*d*(S(3)*m + S(2)*n + S(-2)) + b*c*(m + S(-1)))*cos(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n)), x)


def replacement2451(a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-1))*Simp(a*d*n - b*c*(m + S(1)) - b*d*(m + n + S(1))*sin(e + f*x), x), x), x) + Simp(b*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2452(a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-1))*Simp(a*d*n - b*c*(m + S(1)) - b*d*(m + n + S(1))*cos(e + f*x), x), x), x) - Simp(b*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2453(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-2))*Simp(a*c*d*(m - n + S(1)) + b*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(-1))) + d*(a*d*(m - n + S(1)) + b*c*(m + n))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(-1))*(-a*d + b*c)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2454(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-2))*Simp(a*c*d*(m - n + S(1)) + b*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(-1))) + d*(a*d*(m - n + S(1)) + b*c*(m + n))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(-1))*(-a*d + b*c)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2455(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(-a*d*(S(2)*m + n + S(2)) + b*c*(m + S(1)) + b*d*(m + n + S(2))*sin(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(a*f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2456(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(-a*d*(S(2)*m + n + S(2)) + b*c*(m + S(1)) + b*d*(m + n + S(2))*cos(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(a*f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2457(a, b, c, d, e, f, n, x):
    return -Dist(d/(a*b), Int((c + d*sin(e + f*x))**(n + S(-2))*Simp(-a*c*n + b*d*(n + S(-1)) + (-a*d*n + b*c*(n + S(-1)))*sin(e + f*x), x), x), x) - Simp((c + d*sin(e + f*x))**(n + S(-1))*(-a*d + b*c)*cos(e + f*x)/(a*f*(a + b*sin(e + f*x))), x)


def replacement2458(a, b, c, d, e, f, n, x):
    return -Dist(d/(a*b), Int((c + d*cos(e + f*x))**(n + S(-2))*Simp(-a*c*n + b*d*(n + S(-1)) + (-a*d*n + b*c*(n + S(-1)))*cos(e + f*x), x), x), x) + Simp((c + d*cos(e + f*x))**(n + S(-1))*(-a*d + b*c)*sin(e + f*x)/(a*f*(a + b*cos(e + f*x))), x)


def replacement2459(a, b, c, d, e, f, n, x):
    return Dist(d/(a*(-a*d + b*c)), Int((c + d*sin(e + f*x))**n*(a*n - b*(n + S(1))*sin(e + f*x)), x), x) - Simp(b**S(2)*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(a*f*(a + b*sin(e + f*x))*(-a*d + b*c)), x)


def replacement2460(a, b, c, d, e, f, n, x):
    return Dist(d/(a*(-a*d + b*c)), Int((c + d*cos(e + f*x))**n*(a*n - b*(n + S(1))*cos(e + f*x)), x), x) + Simp(b**S(2)*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(a*f*(a + b*cos(e + f*x))*(-a*d + b*c)), x)


def replacement2461(a, b, c, d, e, f, n, x):
    return Dist(d*n/(a*b), Int((a - b*sin(e + f*x))*(c + d*sin(e + f*x))**(n + S(-1)), x), x) - Simp(b*(c + d*sin(e + f*x))**n*cos(e + f*x)/(a*f*(a + b*sin(e + f*x))), x)


def replacement2462(a, b, c, d, e, f, n, x):
    return Dist(d*n/(a*b), Int((a - b*cos(e + f*x))*(c + d*cos(e + f*x))**(n + S(-1)), x), x) + Simp(b*(c + d*cos(e + f*x))**n*sin(e + f*x)/(a*f*(a + b*cos(e + f*x))), x)


def replacement2463(a, b, c, d, e, f, n, x):
    return Dist(S(2)*n*(a*d + b*c)/(b*(S(2)*n + S(1))), Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(n + S(-1)), x), x) + Simp(-S(2)*b*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(S(2)*n + S(1))), x)


def replacement2464(a, b, c, d, e, f, n, x):
    return Dist(S(2)*n*(a*d + b*c)/(b*(S(2)*n + S(1))), Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**(n + S(-1)), x), x) + Simp(S(2)*b*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*(S(2)*n + S(1))), x)


def replacement2465(a, b, c, d, e, f, x):
    return Simp(-S(2)*b**S(2)*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))*(a*d + b*c)), x)


def replacement2466(a, b, c, d, e, f, x):
    return Simp(S(2)*b**S(2)*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))*(a*d + b*c)), x)


def replacement2467(a, b, c, d, e, f, n, x):
    return Dist((S(2)*n + S(3))*(-a*d + b*c)/(S(2)*b*(c**S(2) - d**S(2))*(n + S(1))), Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(n + S(1)), x), x) + Simp((c + d*sin(e + f*x))**(n + S(1))*(-a*d + b*c)*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2468(a, b, c, d, e, f, n, x):
    return Dist((S(2)*n + S(3))*(-a*d + b*c)/(S(2)*b*(c**S(2) - d**S(2))*(n + S(1))), Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**(n + S(1)), x), x) - Simp((c + d*cos(e + f*x))**(n + S(1))*(-a*d + b*c)*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2469(a, b, c, d, e, f, x):
    return Dist(-S(2)*b/f, Subst(Int(S(1)/(a*d + b*c - d*x**S(2)), x), x, b*cos(e + f*x)/sqrt(a + b*sin(e + f*x))), x)


def replacement2470(a, b, c, d, e, f, x):
    return Dist(S(2)*b/f, Subst(Int(S(1)/(a*d + b*c - d*x**S(2)), x), x, b*sin(e + f*x)/sqrt(a + b*cos(e + f*x))), x)


def replacement2471(a, b, d, e, f, x):
    return Dist(-S(2)/f, Subst(Int(S(1)/sqrt(S(1) - x**S(2)/a), x), x, b*cos(e + f*x)/sqrt(a + b*sin(e + f*x))), x)


def replacement2472(a, b, d, e, f, x):
    return Dist(S(2)/f, Subst(Int(S(1)/sqrt(S(1) - x**S(2)/a), x), x, b*sin(e + f*x)/sqrt(a + b*cos(e + f*x))), x)


def replacement2473(a, b, c, d, e, f, x):
    return Dist(-S(2)*b/f, Subst(Int(S(1)/(b + d*x**S(2)), x), x, b*cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))), x)


def replacement2474(a, b, c, d, e, f, x):
    return Dist(S(2)*b/f, Subst(Int(S(1)/(b + d*x**S(2)), x), x, b*sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x)))), x)


def replacement2475(a, b, c, d, e, f, n, x):
    return Dist(a**S(2)*cos(e + f*x)/(f*sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), Subst(Int((c + d*x)**n/sqrt(a - b*x), x), x, sin(e + f*x)), x)


def replacement2476(a, b, c, d, e, f, n, x):
    return -Dist(a**S(2)*sin(e + f*x)/(f*sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), Subst(Int((c + d*x)**n/sqrt(a - b*x), x), x, cos(e + f*x)), x)


def replacement2477(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2478(a, b, c, d, e, f, x):
    return Dist(d/b, Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2479(a, b, c, d, e, f, n, x):
    return -Dist(S(1)/(b*(S(2)*n + S(-1))), Int((c + d*sin(e + f*x))**(n + S(-2))*Simp(a*c*d - b*(c**S(2)*(S(2)*n + S(-1)) + S(2)*d**S(2)*(n + S(-1))) + d*(a*d - b*c*(S(4)*n + S(-3)))*sin(e + f*x), x)/sqrt(a + b*sin(e + f*x)), x), x) + Simp(-S(2)*d*(c + d*sin(e + f*x))**(n + S(-1))*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(S(2)*n + S(-1))), x)


def replacement2480(a, b, c, d, e, f, n, x):
    return -Dist(S(1)/(b*(S(2)*n + S(-1))), Int((c + d*cos(e + f*x))**(n + S(-2))*Simp(a*c*d - b*(c**S(2)*(S(2)*n + S(-1)) + S(2)*d**S(2)*(n + S(-1))) + d*(a*d - b*c*(S(4)*n + S(-3)))*cos(e + f*x), x)/sqrt(a + b*cos(e + f*x)), x), x) + Simp(S(2)*d*(c + d*cos(e + f*x))**(n + S(-1))*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*(S(2)*n + S(-1))), x)


def replacement2481(a, b, c, d, e, f, n, x):
    return -Dist(S(1)/(S(2)*b*(c**S(2) - d**S(2))*(n + S(1))), Int((c + d*sin(e + f*x))**(n + S(1))*Simp(a*d - S(2)*b*c*(n + S(1)) + b*d*(S(2)*n + S(3))*sin(e + f*x), x)/sqrt(a + b*sin(e + f*x)), x), x) - Simp(d*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2482(a, b, c, d, e, f, n, x):
    return -Dist(S(1)/(S(2)*b*(c**S(2) - d**S(2))*(n + S(1))), Int((c + d*cos(e + f*x))**(n + S(1))*Simp(a*d - S(2)*b*c*(n + S(1)) + b*d*(S(2)*n + S(3))*cos(e + f*x), x)/sqrt(a + b*cos(e + f*x)), x), x) + Simp(d*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(f*sqrt(a + b*cos(e + f*x))*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2483(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/sqrt(a + b*sin(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x)), x), x)


def replacement2484(a, b, c, d, e, f, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/sqrt(a + b*cos(e + f*x)), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b*cos(e + f*x))/(c + d*cos(e + f*x)), x), x)


def replacement2485(a, b, d, e, f, x):
    return -Dist(sqrt(S(2))/(sqrt(a)*f), Subst(Int(S(1)/sqrt(S(1) - x**S(2)), x), x, b*cos(e + f*x)/(a + b*sin(e + f*x))), x)


def replacement2486(a, b, d, e, f, x):
    return Dist(sqrt(S(2))/(sqrt(a)*f), Subst(Int(S(1)/sqrt(S(1) - x**S(2)), x), x, b*sin(e + f*x)/(a + b*cos(e + f*x))), x)


def replacement2487(a, b, c, d, e, f, x):
    return Dist(-S(2)*a/f, Subst(Int(S(1)/(S(2)*b**S(2) - x**S(2)*(a*c - b*d)), x), x, b*cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))), x)


def replacement2488(a, b, c, d, e, f, x):
    return Dist(S(2)*a/f, Subst(Int(S(1)/(S(2)*b**S(2) - x**S(2)*(a*c - b*d)), x), x, b*sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x)))), x)


def replacement2489(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n)), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(-2))*Simp(b*c**S(2)*(m + n) + d*(a*c*m + b*d*(n + S(-1))) + d*(a*d*m + b*c*(m + S(2)*n + S(-1)))*sin(e + f*x), x), x), x) - Simp(d*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(-1))*cos(e + f*x)/(f*(m + n)), x)


def replacement2490(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n)), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(-2))*Simp(b*c**S(2)*(m + n) + d*(a*c*m + b*d*(n + S(-1))) + d*(a*d*m + b*c*(m + S(2)*n + S(-1)))*cos(e + f*x), x), x), x) + Simp(d*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(-1))*sin(e + f*x)/(f*(m + n)), x)


def replacement2491(a, b, c, d, e, f, m, n, x):
    return Dist(a**m*cos(e + f*x)/(f*sqrt(S(1) - sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), Subst(Int((S(1) + b*x/a)**(m + S(-1)/2)*(c + d*x)**n/sqrt(S(1) - b*x/a), x), x, sin(e + f*x)), x)


def replacement2492(a, b, c, d, e, f, m, n, x):
    return -Dist(a**m*sin(e + f*x)/(f*sqrt(S(1) - cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), Subst(Int((S(1) + b*x/a)**(m + S(-1)/2)*(c + d*x)**n/sqrt(S(1) - b*x/a), x), x, cos(e + f*x)), x)


def replacement2493(a, b, d, e, f, m, n, x):
    return -Dist(b*(d/b)**n*cos(e + f*x)/(f*sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), Subst(Int((a - x)**n*(S(2)*a - x)**(m + S(-1)/2)/sqrt(x), x), x, a - b*sin(e + f*x)), x)


def replacement2494(a, b, d, e, f, m, n, x):
    return Dist(b*(d/b)**n*sin(e + f*x)/(f*sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), Subst(Int((a - x)**n*(S(2)*a - x)**(m + S(-1)/2)/sqrt(x), x), x, a - b*cos(e + f*x)), x)


def replacement2495(a, b, d, e, f, m, n, x):
    return Dist((d/b)**IntPart(n)*(b*sin(e + f*x))**(-FracPart(n))*(d*sin(e + f*x))**FracPart(n), Int((b*sin(e + f*x))**n*(a + b*sin(e + f*x))**m, x), x)


def replacement2496(a, b, d, e, f, m, n, x):
    return Dist((d/b)**IntPart(n)*(b*cos(e + f*x))**(-FracPart(n))*(d*cos(e + f*x))**FracPart(n), Int((b*cos(e + f*x))**n*(a + b*cos(e + f*x))**m, x), x)


def replacement2497(a, b, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*(S(1) + b*sin(e + f*x)/a)**(-FracPart(m))*(a + b*sin(e + f*x))**FracPart(m), Int((d*sin(e + f*x))**n*(S(1) + b*sin(e + f*x)/a)**m, x), x)


def replacement2498(a, b, d, e, f, m, n, x):
    return Dist(a**IntPart(m)*(S(1) + b*cos(e + f*x)/a)**(-FracPart(m))*(a + b*cos(e + f*x))**FracPart(m), Int((d*cos(e + f*x))**n*(S(1) + b*cos(e + f*x)/a)**m, x), x)


def replacement2499(a, b, c, d, e, f, m, n, x):
    return Dist(a**S(2)*cos(e + f*x)/(f*sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, sin(e + f*x)), x)


def replacement2500(a, b, c, d, e, f, m, n, x):
    return -Dist(a**S(2)*sin(e + f*x)/(f*sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), Subst(Int((a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, cos(e + f*x)), x)


def replacement2501(b, c, d, e, f, m, x):
    return Dist(S(2)*c*d/b, Int((b*sin(e + f*x))**(m + S(1)), x), x) + Int((b*sin(e + f*x))**m*(c**S(2) + d**S(2)*sin(e + f*x)**S(2)), x)


def replacement2502(b, c, d, e, f, m, x):
    return Dist(S(2)*c*d/b, Int((b*cos(e + f*x))**(m + S(1)), x), x) + Int((b*cos(e + f*x))**m*(c**S(2) + d**S(2)*cos(e + f*x)**S(2)), x)


def replacement2503(a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(-a*(c**S(2) + d**S(2)) + S(2)*b*c*d) + (a**S(2)*d**S(2) - S(2)*a*b*c*d*(m + S(2)) + b**S(2)*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(a**S(2)*d**S(2) - S(2)*a*b*c*d + b**S(2)*c**S(2))*cos(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2504(a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(-a*(c**S(2) + d**S(2)) + S(2)*b*c*d) + (a**S(2)*d**S(2) - S(2)*a*b*c*d*(m + S(2)) + b**S(2)*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(a**S(2)*d**S(2) - S(2)*a*b*c*d + b**S(2)*c**S(2))*sin(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2505(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*sin(e + f*x))**m*Simp(b*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))) - d*(a*d - S(2)*b*c*(m + S(2)))*sin(e + f*x), x), x), x) - Simp(d**S(2)*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2506(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*cos(e + f*x))**m*Simp(b*(c**S(2)*(m + S(2)) + d**S(2)*(m + S(1))) - d*(a*d - S(2)*b*c*(m + S(2)))*cos(e + f*x), x), x), x) + Simp(d**S(2)*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2507(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-3))*(c + d*sin(e + f*x))**(n + S(1))*Simp(a*d*(n + S(1))*(-S(2)*a*b*d + c*(a**S(2) + b**S(2))) + b*(m + S(-2))*(-a*d + b*c)**S(2) + b*(b**S(2)*(c**S(2) - d**S(2)) + d*n*(S(2)*a*b*c - d*(a**S(2) + b**S(2))) - m*(-a*d + b*c)**S(2))*sin(e + f*x)**S(2) + (-a*(n + S(2))*(-a*d + b*c)**S(2) + b*(n + S(1))*(a*b*c**S(2) - S(3)*a*b*d**S(2) + c*d*(a**S(2) + b**S(2))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*(a**S(2)*d**S(2) - S(2)*a*b*c*d + b**S(2)*c**S(2))*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2508(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-3))*(c + d*cos(e + f*x))**(n + S(1))*Simp(a*d*(n + S(1))*(-S(2)*a*b*d + c*(a**S(2) + b**S(2))) + b*(m + S(-2))*(-a*d + b*c)**S(2) + b*(b**S(2)*(c**S(2) - d**S(2)) + d*n*(S(2)*a*b*c - d*(a**S(2) + b**S(2))) - m*(-a*d + b*c)**S(2))*cos(e + f*x)**S(2) + (-a*(n + S(2))*(-a*d + b*c)**S(2) + b*(n + S(1))*(a*b*c**S(2) - S(3)*a*b*d**S(2) + c*d*(a**S(2) + b**S(2))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*(a**S(2)*d**S(2) - S(2)*a*b*c*d + b**S(2)*c**S(2))*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2509(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*sin(e + f*x))**(m + S(-3))*(c + d*sin(e + f*x))**n*Simp(a**S(3)*d*(m + n) + b**S(2)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b**S(2)*(-a*d*(S(3)*m + S(2)*n + S(-2)) + b*c*(m + S(-1)))*sin(e + f*x)**S(2) - b*(-S(3)*a**S(2)*d*(m + n) + a*b*c - b**S(2)*d*(m + n + S(-1)))*sin(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n)), x)


def replacement2510(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*cos(e + f*x))**(m + S(-3))*(c + d*cos(e + f*x))**n*Simp(a**S(3)*d*(m + n) + b**S(2)*(a*d*(n + S(1)) + b*c*(m + S(-2))) - b**S(2)*(-a*d*(S(3)*m + S(2)*n + S(-2)) + b*c*(m + S(-1)))*cos(e + f*x)**S(2) - b*(-S(3)*a**S(2)*d*(m + n) + a*b*c - b**S(2)*d*(m + n + S(-1)))*cos(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n)), x)


def replacement2511(a, b, d, e, f, x):
    return -Dist(d**S(2)/(a**S(2) - b**S(2)), Int(sqrt(a + b*sin(e + f*x))/(d*sin(e + f*x))**(S(3)/2), x), x) + Simp(-S(2)*a*d*cos(e + f*x)/(f*sqrt(d*sin(e + f*x))*sqrt(a + b*sin(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2512(a, b, d, e, f, x):
    return -Dist(d**S(2)/(a**S(2) - b**S(2)), Int(sqrt(a + b*cos(e + f*x))/(d*cos(e + f*x))**(S(3)/2), x), x) + Simp(S(2)*a*d*sin(e + f*x)/(f*sqrt(d*cos(e + f*x))*sqrt(a + b*cos(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2513(a, b, c, d, e, f, x):
    return Dist((c - d)/(a - b), Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x) - Dist((-a*d + b*c)/(a - b), Int((sin(e + f*x) + S(1))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2514(a, b, c, d, e, f, x):
    return Dist((c - d)/(a - b), Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x) - Dist((-a*d + b*c)/(a - b), Int((cos(e + f*x) + S(1))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2515(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-1))*Simp(a*c*(m + S(1)) + b*d*n - b*d*(m + n + S(2))*sin(e + f*x)**S(2) + (a*d*(m + S(1)) - b*c*(m + S(2)))*sin(e + f*x), x), x), x) - Simp(b*(a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2516(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-1))*Simp(a*c*(m + S(1)) + b*d*n - b*d*(m + n + S(2))*cos(e + f*x)**S(2) + (a*d*(m + S(1)) - b*c*(m + S(2)))*cos(e + f*x), x), x), x) + Simp(b*(a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2517(a, b, d, e, f, x):
    return Dist(d/b, Int(sqrt(d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x) - Dist(a*d/b, Int(sqrt(d*sin(e + f*x))/(a + b*sin(e + f*x))**(S(3)/2), x), x)


def replacement2518(a, b, d, e, f, x):
    return Dist(d/b, Int(sqrt(d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x) - Dist(a*d/b, Int(sqrt(d*cos(e + f*x))/(a + b*cos(e + f*x))**(S(3)/2), x), x)


def replacement2519(a, b, c, d, e, f, x):
    return Dist(d**S(2)/b**S(2), Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/b**S(2), Int(Simp(a*d + b*c + S(2)*b*d*sin(e + f*x), x)/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2520(a, b, c, d, e, f, x):
    return Dist(d**S(2)/b**S(2), Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/b**S(2), Int(Simp(a*d + b*c + S(2)*b*d*cos(e + f*x), x)/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2521(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-2))*Simp(c*(m + S(1))*(a*c - b*d) + d*(n + S(-1))*(-a*d + b*c) - d*(-a*d + b*c)*(m + n + S(1))*sin(e + f*x)**S(2) + (-c*(m + S(2))*(-a*d + b*c) + d*(m + S(1))*(a*c - b*d))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-1))*(-a*d + b*c)*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2522(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-2))*Simp(c*(m + S(1))*(a*c - b*d) + d*(n + S(-1))*(-a*d + b*c) - d*(-a*d + b*c)*(m + n + S(1))*cos(e + f*x)**S(2) + (-c*(m + S(2))*(-a*d + b*c) + d*(m + S(1))*(a*c - b*d))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-1))*(-a*d + b*c)*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2523(a, b, d, e, f, x):
    return Dist(d/(a**S(2) - b**S(2)), Int((a*sin(e + f*x) + b)/((d*sin(e + f*x))**(S(3)/2)*sqrt(a + b*sin(e + f*x))), x), x) + Simp(S(2)*b*cos(e + f*x)/(f*sqrt(d*sin(e + f*x))*sqrt(a + b*sin(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2524(a, b, d, e, f, x):
    return Dist(d/(a**S(2) - b**S(2)), Int((a*cos(e + f*x) + b)/((d*cos(e + f*x))**(S(3)/2)*sqrt(a + b*cos(e + f*x))), x), x) + Simp(-S(2)*b*sin(e + f*x)/(f*sqrt(d*cos(e + f*x))*sqrt(a + b*cos(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2525(a, b, c, d, e, f, x):
    return -Dist(b/(a - b), Int((sin(e + f*x) + S(1))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x) + Dist(S(1)/(a - b), Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2526(a, b, c, d, e, f, x):
    return -Dist(b/(a - b), Int((cos(e + f*x) + S(1))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x) + Dist(S(1)/(a - b), Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2527(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(a*(m + S(1))*(-a*d + b*c) + b**S(2)*d*(m + n + S(2)) - b**S(2)*d*(m + n + S(3))*sin(e + f*x)**S(2) - (b**S(2)*c + b*(m + S(1))*(-a*d + b*c))*sin(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement2528(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(a*(m + S(1))*(-a*d + b*c) + b**S(2)*d*(m + n + S(2)) - b**S(2)*d*(m + n + S(3))*cos(e + f*x)**S(2) - (b**S(2)*c + b*(m + S(1))*(-a*d + b*c))*cos(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement2529(a, b, c, d, e, f, x):
    return Dist(d/b, Int(S(1)/sqrt(c + d*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/((a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2530(a, b, c, d, e, f, x):
    return Dist(d/b, Int(S(1)/sqrt(c + d*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/b, Int(S(1)/((a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2531(a, b, c, d, e, f, x):
    return Dist(b/d, Int(sqrt(a + b*sin(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x)), x), x)


def replacement2532(a, b, c, d, e, f, x):
    return Dist(b/d, Int(sqrt(a + b*cos(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(sqrt(a + b*cos(e + f*x))/(c + d*cos(e + f*x)), x), x)


def replacement2533(a, b, c, d, e, f, x):
    return Simp(S(2)*EllipticPi(S(2)*b/(a + b), -Pi/S(4) + e/S(2) + f*x/S(2), S(2)*d/(c + d))/(f*(a + b)*sqrt(c + d)), x)


def replacement2534(a, b, c, d, e, f, x):
    return Simp(S(2)*EllipticPi(S(2)*b/(a + b), e/S(2) + f*x/S(2), S(2)*d/(c + d))/(f*(a + b)*sqrt(c + d)), x)


def replacement2535(a, b, c, d, e, f, x):
    return Simp(S(2)*EllipticPi(-S(2)*b/(a - b), Pi/S(4) + e/S(2) + f*x/S(2), -S(2)*d/(c - d))/(f*(a - b)*sqrt(c - d)), x)


def replacement2536(a, b, c, d, e, f, x):
    return Simp(S(2)*EllipticPi(-S(2)*b/(a - b), Pi/S(2) + e/S(2) + f*x/S(2), -S(2)*d/(c - d))/(f*(a - b)*sqrt(c - d)), x)


def replacement2537(a, b, c, d, e, f, x):
    return Dist(sqrt((c + d*sin(e + f*x))/(c + d))/sqrt(c + d*sin(e + f*x)), Int(S(1)/((a + b*sin(e + f*x))*sqrt(c/(c + d) + d*sin(e + f*x)/(c + d))), x), x)


def replacement2538(a, b, c, d, e, f, x):
    return Dist(sqrt((c + d*cos(e + f*x))/(c + d))/sqrt(c + d*cos(e + f*x)), Int(S(1)/((a + b*cos(e + f*x))*sqrt(c/(c + d) + d*cos(e + f*x)/(c + d))), x), x)


def replacement2539(b, c, d, e, f, x):
    return Simp(S(2)*c*sqrt(S(1) - S(1)/sin(e + f*x))*sqrt(S(1) + S(1)/sin(e + f*x))*EllipticPi((c + d)/d, asin(sqrt(c + d*sin(e + f*x))/(sqrt(b*sin(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt(b*(c + d), S(2))*tan(e + f*x)/(d*f*sqrt(c**S(2) - d**S(2))), x)


def replacement2540(b, c, d, e, f, x):
    return Simp(-S(2)*c*sqrt(S(1) - S(1)/cos(e + f*x))*sqrt(S(1) + S(1)/cos(e + f*x))*EllipticPi((c + d)/d, asin(sqrt(c + d*cos(e + f*x))/(sqrt(b*cos(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt(b*(c + d), S(2))/(d*f*sqrt(c**S(2) - d**S(2))*tan(e + f*x)), x)


def replacement2541(b, c, d, e, f, x):
    return Simp(S(2)*b*sqrt(c*(S(1) - S(1)/sin(e + f*x))/(c + d))*sqrt(c*(S(1) + S(1)/sin(e + f*x))/(c - d))*EllipticPi((c + d)/d, asin(sqrt(c + d*sin(e + f*x))/(sqrt(b*sin(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt((c + d)/b, S(2))*tan(e + f*x)/(d*f), x)


def replacement2542(b, c, d, e, f, x):
    return Simp(-S(2)*b*sqrt(c*(S(1) - S(1)/cos(e + f*x))/(c + d))*sqrt(c*(S(1) + S(1)/cos(e + f*x))/(c - d))*EllipticPi((c + d)/d, asin(sqrt(c + d*cos(e + f*x))/(sqrt(b*cos(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt((c + d)/b, S(2))/(d*f*tan(e + f*x)), x)


def replacement2543(b, c, d, e, f, x):
    return Dist(sqrt(b*sin(e + f*x))/sqrt(-b*sin(e + f*x)), Int(sqrt(-b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x)


def replacement2544(b, c, d, e, f, x):
    return Dist(sqrt(b*cos(e + f*x))/sqrt(-b*cos(e + f*x)), Int(sqrt(-b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x)


def replacement2545(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((-a*d + b*c)*(sin(e + f*x) + S(1))/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(S(1) - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*EllipticPi(b*(c + d)/(d*(a + b)), asin(sqrt(c + d*sin(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*sin(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(d*f*Rt((a + b)/(c + d), S(2))*cos(e + f*x)), x)


def replacement2546(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((-a*d + b*c)*(cos(e + f*x) + S(1))/((a + b*cos(e + f*x))*(c - d)))*sqrt(-(S(1) - cos(e + f*x))*(-a*d + b*c)/((a + b*cos(e + f*x))*(c + d)))*(a + b*cos(e + f*x))*EllipticPi(b*(c + d)/(d*(a + b)), asin(sqrt(c + d*cos(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*cos(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(d*f*Rt((a + b)/(c + d), S(2))*sin(e + f*x)), x)


def replacement2547(a, b, c, d, e, f, x):
    return Dist(sqrt(-c - d*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), Int(sqrt(a + b*sin(e + f*x))/sqrt(-c - d*sin(e + f*x)), x), x)


def replacement2548(a, b, c, d, e, f, x):
    return Dist(sqrt(-c - d*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), Int(sqrt(a + b*cos(e + f*x))/sqrt(-c - d*cos(e + f*x)), x), x)


def replacement2549(a, b, d, e, f, x):
    return Simp(-S(2)*d*EllipticF(asin(cos(e + f*x)/(d*sin(e + f*x) + S(1))), -(a - b*d)/(a + b*d))/(f*sqrt(a + b*d)), x)


def replacement2550(a, b, d, e, f, x):
    return Simp(S(2)*d*EllipticF(asin(sin(e + f*x)/(d*cos(e + f*x) + S(1))), -(a - b*d)/(a + b*d))/(f*sqrt(a + b*d)), x)


def replacement2551(a, b, d, e, f, x):
    return Dist(sqrt(sin(e + f*x)*sign(b))/sqrt(d*sin(e + f*x)), Int(S(1)/(sqrt(sin(e + f*x)*sign(b))*sqrt(a + b*sin(e + f*x))), x), x)


def replacement2552(a, b, d, e, f, x):
    return Dist(sqrt(cos(e + f*x)*sign(b))/sqrt(d*cos(e + f*x)), Int(S(1)/(sqrt(cos(e + f*x)*sign(b))*sqrt(a + b*cos(e + f*x))), x), x)


def replacement2553(a, b, d, e, f, x):
    return Simp(-S(2)*sqrt(-S(1)/tan(e + f*x)**S(2))*sqrt(a**S(2))*EllipticF(asin(sqrt(a + b*sin(e + f*x))/(sqrt(d*sin(e + f*x))*Rt((a + b)/d, S(2)))), -(a + b)/(a - b))*Rt((a + b)/d, S(2))*tan(e + f*x)/(a*f*sqrt(a**S(2) - b**S(2))), x)


def replacement2554(a, b, d, e, f, x):
    return Simp(S(2)*sqrt(-tan(e + f*x)**S(2))*sqrt(a**S(2))*EllipticF(asin(sqrt(a + b*cos(e + f*x))/(sqrt(d*cos(e + f*x))*Rt((a + b)/d, S(2)))), -(a + b)/(a - b))*Rt((a + b)/d, S(2))/(a*f*sqrt(a**S(2) - b**S(2))*tan(e + f*x)), x)


def replacement2555(a, b, d, e, f, x):
    return Simp(-S(2)*sqrt(a*(S(1) - S(1)/sin(e + f*x))/(a + b))*sqrt(a*(S(1) + S(1)/sin(e + f*x))/(a - b))*EllipticF(asin(sqrt(a + b*sin(e + f*x))/(sqrt(d*sin(e + f*x))*Rt((a + b)/d, S(2)))), -(a + b)/(a - b))*Rt((a + b)/d, S(2))*tan(e + f*x)/(a*f), x)


def replacement2556(a, b, d, e, f, x):
    return Simp(S(2)*sqrt(a*(S(1) - S(1)/cos(e + f*x))/(a + b))*sqrt(a*(S(1) + S(1)/cos(e + f*x))/(a - b))*EllipticF(asin(sqrt(a + b*cos(e + f*x))/(sqrt(d*cos(e + f*x))*Rt((a + b)/d, S(2)))), -(a + b)/(a - b))*Rt((a + b)/d, S(2))/(a*f*tan(e + f*x)), x)


def replacement2557(a, b, d, e, f, x):
    return Dist(sqrt(-d*sin(e + f*x))/sqrt(d*sin(e + f*x)), Int(S(1)/(sqrt(-d*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), x), x)


def replacement2558(a, b, d, e, f, x):
    return Dist(sqrt(-d*cos(e + f*x))/sqrt(d*cos(e + f*x)), Int(S(1)/(sqrt(-d*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), x), x)


def replacement2559(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((S(1) - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + S(1))/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*EllipticF(asin(sqrt(a + b*sin(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d*sin(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))/(f*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))*cos(e + f*x)), x)


def replacement2560(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((S(1) - cos(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*cos(e + f*x))))*sqrt(-(-a*d + b*c)*(cos(e + f*x) + S(1))/((a - b)*(c + d*cos(e + f*x))))*(c + d*cos(e + f*x))*EllipticF(asin(sqrt(a + b*cos(e + f*x))*Rt((c + d)/(a + b), S(2))/sqrt(c + d*cos(e + f*x))), (a + b)*(c - d)/((a - b)*(c + d)))/(f*(-a*d + b*c)*Rt((c + d)/(a + b), S(2))*sin(e + f*x)), x)


def replacement2561(a, b, c, d, e, f, x):
    return Dist(sqrt(-a - b*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), Int(S(1)/(sqrt(-a - b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2562(a, b, c, d, e, f, x):
    return Dist(sqrt(-a - b*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), Int(S(1)/(sqrt(-a - b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2563(a, b, d, e, f, x):
    return Dist(d/(S(2)*b), Int(sqrt(d*sin(e + f*x))*(a + S(2)*b*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x) - Dist(a*d/(S(2)*b), Int(sqrt(d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x)


def replacement2564(a, b, d, e, f, x):
    return Dist(d/(S(2)*b), Int(sqrt(d*cos(e + f*x))*(a + S(2)*b*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x) - Dist(a*d/(S(2)*b), Int(sqrt(d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x)


def replacement2565(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(-1))*Simp(a**S(2)*c*d*(m + n) + b*d*(a*d*n + b*c*(m + S(-1))) + b*d*(a*d*(S(2)*m + n + S(-1)) + b*c*n)*sin(e + f*x)**S(2) + (a*d*(m + n)*(a*d + S(2)*b*c) - b*d*(a*c - b*d*(m + n + S(-1))))*sin(e + f*x), x), x), x) - Simp(b*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(m + n)), x)


def replacement2566(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n)), Int((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(-1))*Simp(a**S(2)*c*d*(m + n) + b*d*(a*d*n + b*c*(m + S(-1))) + b*d*(a*d*(S(2)*m + n + S(-1)) + b*c*n)*cos(e + f*x)**S(2) + (a*d*(m + n)*(a*d + S(2)*b*c) - b*d*(a*c - b*d*(m + n + S(-1))))*cos(e + f*x), x), x), x) + Simp(b*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(m + n)), x)


def replacement2567(a, b, c, d, e, f, m, n, x):
    return Dist(b/d, Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1)), x), x) - Dist((-a*d + b*c)/d, Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n, x), x)


def replacement2568(a, b, c, d, e, f, m, n, x):
    return Dist(b/d, Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1)), x), x) - Dist((-a*d + b*c)/d, Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n, x), x)


def replacement2569(a, b, c, d, e, f, m, n, x):
    return Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement2570(a, b, c, d, e, f, m, n, x):
    return Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement2571(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d*sin(e + f*x))**p)**FracPart(n)*(d*sin(e + f*x))**(-p*FracPart(n)), Int((d*sin(e + f*x))**(n*p)*(a + b*sin(e + f*x))**m, x), x)


def replacement2572(a, b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(c*(d*cos(e + f*x))**p)**FracPart(n)*(d*cos(e + f*x))**(-p*FracPart(n)), Int((d*cos(e + f*x))**(n*p)*(a + b*cos(e + f*x))**m, x), x)


def replacement2573(a, b, c, d, e, f, m, n, x):
    return Int((a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-n), x)


def replacement2574(a, b, c, d, e, f, m, n, x):
    return Int((a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-n), x)


def replacement2575(a, b, c, d, e, f, m, n, x):
    return Int((c + d/sin(e + f*x))**n*(a/sin(e + f*x) + b)**m*(S(1)/sin(e + f*x))**(-m), x)


def replacement2576(a, b, c, d, e, f, m, n, x):
    return Int((c + d/cos(e + f*x))**n*(a/cos(e + f*x) + b)**m*(S(1)/cos(e + f*x))**(-m), x)


def replacement2577(a, b, c, d, e, f, m, n, x):
    return Dist((c + d/sin(e + f*x))**n*(c*sin(e + f*x) + d)**(-n)*sin(e + f*x)**n, Int((a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-n), x), x)


def replacement2578(a, b, c, d, e, f, m, n, x):
    return Dist((c + d/cos(e + f*x))**n*(c*cos(e + f*x) + d)**(-n)*cos(e + f*x)**n, Int((a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-n), x), x)


def replacement2579(a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n, x), x, b*sin(e + f*x)), x)


def replacement2580(a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(b*f), Subst(Int((a + x)**m*(c + d*x/b)**n, x), x, b*cos(e + f*x)), x)


def replacement2581(a, b, d, e, f, n, p, x):
    return Dist(a, Int((d*sin(e + f*x))**n*cos(e + f*x)**p, x), x) + Dist(b/d, Int((d*sin(e + f*x))**(n + S(1))*cos(e + f*x)**p, x), x)


def replacement2582(a, b, d, e, f, n, p, x):
    return Dist(a, Int((d*cos(e + f*x))**n*sin(e + f*x)**p, x), x) + Dist(b/d, Int((d*cos(e + f*x))**(n + S(1))*sin(e + f*x)**p, x), x)


def replacement2583(a, b, d, e, f, n, p, x):
    return Dist(S(1)/a, Int((d*sin(e + f*x))**n*cos(e + f*x)**(p + S(-2)), x), x) - Dist(S(1)/(b*d), Int((d*sin(e + f*x))**(n + S(1))*cos(e + f*x)**(p + S(-2)), x), x)


def replacement2584(a, b, d, e, f, n, p, x):
    return Dist(S(1)/a, Int((d*cos(e + f*x))**n*sin(e + f*x)**(p + S(-2)), x), x) - Dist(S(1)/(b*d), Int((d*cos(e + f*x))**(n + S(1))*sin(e + f*x)**(p + S(-2)), x), x)


def replacement2585(a, b, c, d, e, f, m, n, p, x):
    return Dist(b**(-p)/f, Subst(Int((a - x)**(p/S(2) + S(-1)/2)*(a + x)**(m + p/S(2) + S(-1)/2)*(c + d*x/b)**n, x), x, b*sin(e + f*x)), x)


def replacement2586(a, b, c, d, e, f, m, n, p, x):
    return -Dist(b**(-p)/f, Subst(Int((a - x)**(p/S(2) + S(-1)/2)*(a + x)**(m + p/S(2) + S(-1)/2)*(c + d*x/b)**n, x), x, b*cos(e + f*x)), x)


def replacement2587(a, b, c, d, e, f, m, n, p, x):
    return Dist(b**(-p)/f, Subst(Int((a + x)**m*(b**S(2) - x**S(2))**(p/S(2) + S(-1)/2)*(c + d*x/b)**n, x), x, b*sin(e + f*x)), x)


def replacement2588(a, b, c, d, e, f, m, n, p, x):
    return -Dist(b**(-p)/f, Subst(Int((a + x)**m*(b**S(2) - x**S(2))**(p/S(2) + S(-1)/2)*(c + d*x/b)**n, x), x, b*cos(e + f*x)), x)


def replacement2589(a, b, d, e, f, g, n, p, x):
    return Dist(a, Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**p, x), x) + Dist(b/d, Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**p, x), x)


def replacement2590(a, b, d, e, f, g, n, p, x):
    return Dist(a, Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**p, x), x) + Dist(b/d, Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**p, x), x)


def replacement2591(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/a, Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2)), x), x) - Dist(g**S(2)/(b*d), Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**(p + S(-2)), x), x)


def replacement2592(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/a, Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2)), x), x) - Dist(g**S(2)/(b*d), Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**(p + S(-2)), x), x)


def replacement2593(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**m*c**m*g**(-S(2)*m), Int((g*cos(e + f*x))**(S(2)*m + p)*(c + d*sin(e + f*x))**(-m + n), x), x)


def replacement2594(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**m*c**m*g**(-S(2)*m), Int((g*sin(e + f*x))**(S(2)*m + p)*(c + d*cos(e + f*x))**(-m + n), x), x)


def replacement2595(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**(-p/S(2))*c**(-p/S(2)), Int((a + b*sin(e + f*x))**(m + p/S(2))*(c + d*sin(e + f*x))**(n + p/S(2)), x), x)


def replacement2596(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**(-p/S(2))*c**(-p/S(2)), Int((a + b*cos(e + f*x))**(m + p/S(2))*(c + d*cos(e + f*x))**(n + p/S(2)), x), x)


def replacement2597(a, b, c, d, e, f, g, p, x):
    return Dist(g*cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), Int((g*cos(e + f*x))**(p + S(-1)), x), x)


def replacement2598(a, b, c, d, e, f, g, p, x):
    return Dist(g*sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), Int((g*sin(e + f*x))**(p + S(-1)), x), x)


def replacement2599(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*cos(e + f*x))**(-S(2)*FracPart(m))*(a + b*sin(e + f*x))**FracPart(m)*(c + d*sin(e + f*x))**FracPart(m), Int((g*cos(e + f*x))**(S(2)*m + p)/(c + d*sin(e + f*x)), x), x)


def replacement2600(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*sin(e + f*x))**(-S(2)*FracPart(m))*(a + b*cos(e + f*x))**FracPart(m)*(c + d*cos(e + f*x))**FracPart(m), Int((g*sin(e + f*x))**(S(2)*m + p)/(c + d*cos(e + f*x)), x), x)


def replacement2601(a, b, c, d, e, f, g, m, n, p, x):
    return Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n/(f*g*(m - n + S(-1))), x)


def replacement2602(a, b, c, d, e, f, g, m, n, p, x):
    return -Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n/(f*g*(m - n + S(-1))), x)


def replacement2603(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(b*(S(2)*m + p + S(-1))/(d*(S(2)*n + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1)), x), x) + Simp(-S(2)*b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n/(f*g*(S(2)*n + p + S(1))), x)


def replacement2604(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(b*(S(2)*m + p + S(-1))/(d*(S(2)*n + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1)), x), x) + Simp(S(2)*b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n/(f*g*(S(2)*n + p + S(1))), x)


def replacement2605(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + n + p), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n, x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n/(f*g*(m + n + p)), x)


def replacement2606(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + n + p), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n, x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n/(f*g*(m + n + p)), x)


def replacement2607(a, b, c, d, e, f, g, m, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*cos(e + f*x))**(-S(2)*FracPart(m))*(a + b*sin(e + f*x))**FracPart(m)*(c + d*sin(e + f*x))**FracPart(m), Int((g*cos(e + f*x))**(S(2)*m + p), x), x)


def replacement2608(a, b, c, d, e, f, g, m, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*sin(e + f*x))**(-S(2)*FracPart(m))*(a + b*cos(e + f*x))**FracPart(m)*(c + d*cos(e + f*x))**FracPart(m), Int((g*sin(e + f*x))**(S(2)*m + p), x), x)


def replacement2609(a, b, c, d, e, f, g, m, n, p, x):
    return Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n/(a*f*g*(m - n)), x)


def replacement2610(a, b, c, d, e, f, g, m, n, p, x):
    return -Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n/(a*f*g*(m - n)), x)


def replacement2611(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((m + n + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2612(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((m + n + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2613(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(b*(S(2)*m + p + S(-1))/(d*(S(2)*n + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1)), x), x) + Simp(-S(2)*b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n/(f*g*(S(2)*n + p + S(1))), x)


def replacement2614(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(b*(S(2)*m + p + S(-1))/(d*(S(2)*n + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1)), x), x) + Simp(S(2)*b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n/(f*g*(S(2)*n + p + S(1))), x)


def replacement2615(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + n + p), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n, x), x) - Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n/(f*g*(m + n + p)), x)


def replacement2616(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a*(S(2)*m + p + S(-1))/(m + n + p), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n, x), x) + Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n/(f*g*(m + n + p)), x)


def replacement2617(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((m + n + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2618(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((m + n + p + S(1))/(a*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2619(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*cos(e + f*x))**(-S(2)*FracPart(m))*(a + b*sin(e + f*x))**FracPart(m)*(c + d*sin(e + f*x))**FracPart(m), Int((g*cos(e + f*x))**(S(2)*m + p)*(c + d*sin(e + f*x))**(-m + n), x), x)


def replacement2620(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**IntPart(m)*c**IntPart(m)*g**(-S(2)*IntPart(m))*(g*sin(e + f*x))**(-S(2)*FracPart(m))*(a + b*cos(e + f*x))**FracPart(m)*(c + d*cos(e + f*x))**FracPart(m), Int((g*sin(e + f*x))**(S(2)*m + p)*(c + d*cos(e + f*x))**(-m + n), x), x)


def replacement2621(a, b, c, d, e, f, g, m, p, x):
    return -Simp(d*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2622(a, b, c, d, e, f, g, m, p, x):
    return Simp(d*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2623(a, b, c, d, e, f, g, m, p, x):
    return Dist(b*(a*d*m + b*c*(m + p + S(1)))/(a*g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1)), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(a*d + b*c)/(a*f*g*(p + S(1))), x)


def replacement2624(a, b, c, d, e, f, g, m, p, x):
    return Dist(b*(a*d*m + b*c*(m + p + S(1)))/(a*g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1)), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(a*d + b*c)/(a*f*g*(p + S(1))), x)


def replacement2625(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(b*(m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**m, x), x) - Simp(d*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2626(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(b*(m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**m, x), x) + Simp(d*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2627(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b**S(3)*(S(2)*m + S(3))), Int((a + b*sin(e + f*x))**(m + S(2))*(S(2)*a*d*(m + S(1)) + b*c - b*d*(S(2)*m + S(3))*sin(e + f*x)), x), x) + Simp(S(2)*(a + b*sin(e + f*x))**(m + S(1))*(-a*d + b*c)*cos(e + f*x)/(b**S(2)*f*(S(2)*m + S(3))), x)


def replacement2628(a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b**S(3)*(S(2)*m + S(3))), Int((a + b*cos(e + f*x))**(m + S(2))*(S(2)*a*d*(m + S(1)) + b*c - b*d*(S(2)*m + S(3))*cos(e + f*x)), x), x) + Simp(-S(2)*(a + b*cos(e + f*x))**(m + S(1))*(-a*d + b*c)*sin(e + f*x)/(b**S(2)*f*(S(2)*m + S(3))), x)


def replacement2629(a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(m + S(3))), Int((a + b*sin(e + f*x))**(m + S(1))*(-a*c*(m + S(3)) + b*d*(m + S(2)) + (-a*d*(m + S(4)) + b*c*(m + S(3)))*sin(e + f*x)), x), x) + Simp(d*(a + b*sin(e + f*x))**(m + S(2))*cos(e + f*x)/(b**S(2)*f*(m + S(3))), x)


def replacement2630(a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(m + S(3))), Int((a + b*cos(e + f*x))**(m + S(1))*(-a*c*(m + S(3)) + b*d*(m + S(2)) + (-a*d*(m + S(4)) + b*c*(m + S(3)))*cos(e + f*x)), x), x) - Simp(d*(a + b*cos(e + f*x))**(m + S(2))*sin(e + f*x)/(b**S(2)*f*(m + S(3))), x)


def replacement2631(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(a*b*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1)), x), x) + Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(-a*d + b*c)/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2632(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(a*b*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(-a*d + b*c)/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2633(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(b*(m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**m, x), x) - Simp(d*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2634(a, b, c, d, e, f, g, m, p, x):
    return Dist((a*d*m + b*c*(m + p + S(1)))/(b*(m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**m, x), x) + Simp(d*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2635(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1))*Simp(a*c*(p + S(2)) + b*c*(m + p + S(2))*sin(e + f*x) + b*d*m, x), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)/(f*g*(p + S(1))), x)


def replacement2636(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1))*Simp(a*c*(p + S(2)) + b*c*(m + p + S(2))*cos(e + f*x) + b*d*m, x), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)/(f*g*(p + S(1))), x)


def replacement2637(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(m + p + S(1)), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(-1))*Simp(a*c*(m + p + S(1)) + b*d*m + (a*d*m + b*c*(m + p + S(1)))*sin(e + f*x), x), x), x) - Simp(d*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2638(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(m + p + S(1)), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(-1))*Simp(a*c*(m + p + S(1)) + b*d*m + (a*d*m + b*c*(m + p + S(1)))*cos(e + f*x), x), x), x) + Simp(d*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(f*g*(m + p + S(1))), x)


def replacement2639(a, b, c, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(m + S(1))*(m + p + S(1))), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1))*Simp(b*d*(m + S(1)) + (-a*d*p + b*c*(m + p + S(1)))*sin(e + f*x), x), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))*(-a*d*p + b*c*(m + p + S(1)) + b*d*(m + S(1))*sin(e + f*x))/(b**S(2)*f*(m + S(1))*(m + p + S(1))), x)


def replacement2640(a, b, c, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(m + S(1))*(m + p + S(1))), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1))*Simp(b*d*(m + S(1)) + (-a*d*p + b*c*(m + p + S(1)))*cos(e + f*x), x), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))*(-a*d*p + b*c*(m + p + S(1)) + b*d*(m + S(1))*cos(e + f*x))/(b**S(2)*f*(m + S(1))*(m + p + S(1))), x)


def replacement2641(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1))*Simp((m + S(1))*(a*c - b*d) - (-a*d + b*c)*(m + p + S(2))*sin(e + f*x), x), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(-a*d + b*c)/(f*g*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2642(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1))*Simp((m + S(1))*(a*c - b*d) - (-a*d + b*c)*(m + p + S(2))*cos(e + f*x), x), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(-a*d + b*c)/(f*g*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2643(a, b, c, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(m + p)*(m + p + S(1))), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**m*Simp(b*(a*d*m + b*c*(m + p + S(1))) + (a*b*c*(m + p + S(1)) - d*(a**S(2)*p - b**S(2)*(m + p)))*sin(e + f*x), x), x), x) + Simp(g*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))*(-a*d*p + b*c*(m + p + S(1)) + b*d*(m + p)*sin(e + f*x))/(b**S(2)*f*(m + p)*(m + p + S(1))), x)


def replacement2644(a, b, c, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(p + S(-1))/(b**S(2)*(m + p)*(m + p + S(1))), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**m*Simp(b*(a*d*m + b*c*(m + p + S(1))) + (a*b*c*(m + p + S(1)) - d*(a**S(2)*p - b**S(2)*(m + p)))*cos(e + f*x), x), x), x) - Simp(g*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))*(-a*d*p + b*c*(m + p + S(1)) + b*d*(m + p)*cos(e + f*x))/(b**S(2)*f*(m + p)*(m + p + S(1))), x)


def replacement2645(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(a**S(2) - b**S(2))*(p + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m*Simp(a*b*d*m + b*(a*c - b*d)*(m + p + S(3))*sin(e + f*x) + c*(a**S(2)*(p + S(2)) - b**S(2)*(m + p + S(2))), x), x), x) + Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(-a*d + b*c - (a*c - b*d)*sin(e + f*x))/(f*g*(a**S(2) - b**S(2))*(p + S(1))), x)


def replacement2646(a, b, c, d, e, f, g, m, p, x):
    return Dist(S(1)/(g**S(2)*(a**S(2) - b**S(2))*(p + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m*Simp(a*b*d*m + b*(a*c - b*d)*(m + p + S(3))*cos(e + f*x) + c*(a**S(2)*(p + S(2)) - b**S(2)*(m + p + S(2))), x), x), x) - Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(-a*d + b*c - (a*c - b*d)*cos(e + f*x))/(f*g*(a**S(2) - b**S(2))*(p + S(1))), x)


def replacement2647(a, b, c, d, e, f, g, p, x):
    return Dist(d/b, Int((g*cos(e + f*x))**p, x), x) + Dist((-a*d + b*c)/b, Int((g*cos(e + f*x))**p/(a + b*sin(e + f*x)), x), x)


def replacement2648(a, b, c, d, e, f, g, p, x):
    return Dist(d/b, Int((g*sin(e + f*x))**p, x), x) + Dist((-a*d + b*c)/b, Int((g*sin(e + f*x))**p/(a + b*cos(e + f*x)), x), x)


def replacement2649(a, b, c, d, e, f, g, m, p, x):
    return Dist(c*g*(g*cos(e + f*x))**(p + S(-1))*(S(1) - sin(e + f*x))**(S(1)/2 - p/S(2))*(sin(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((S(1) - d*x/c)**(p/S(2) + S(-1)/2)*(S(1) + d*x/c)**(p/S(2) + S(1)/2)*(a + b*x)**m, x), x, sin(e + f*x)), x)


def replacement2650(a, b, c, d, e, f, g, m, p, x):
    return -Dist(c*g*(g*sin(e + f*x))**(p + S(-1))*(S(1) - cos(e + f*x))**(S(1)/2 - p/S(2))*(cos(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((S(1) - d*x/c)**(p/S(2) + S(-1)/2)*(S(1) + d*x/c)**(p/S(2) + S(1)/2)*(a + b*x)**m, x), x, cos(e + f*x)), x)


def replacement2651(a, b, d, e, f, m, n, p, x):
    return Dist(a**(S(2)*m), Int((d*sin(e + f*x))**n*(a - b*sin(e + f*x))**(-m), x), x)


def replacement2652(a, b, d, e, f, m, n, p, x):
    return Dist(a**(S(2)*m), Int((d*cos(e + f*x))**n*(a - b*cos(e + f*x))**(-m), x), x)


def replacement2653(a, b, e, f, g, m, p, x):
    return Dist(a/(S(2)*g**S(2)), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1)), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))/(S(2)*b*f*g*(m + S(1))), x)


def replacement2654(a, b, e, f, g, m, p, x):
    return Dist(a/(S(2)*g**S(2)), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1)), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))/(S(2)*b*f*g*(m + S(1))), x)


def replacement2655(a, b, e, f, g, m, p, x):
    return -Dist(g**(S(-2)), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**m, x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*m), x)


def replacement2656(a, b, e, f, g, m, p, x):
    return -Dist(g**(S(-2)), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**m, x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*m), x)


def replacement2657(a, b, d, e, f, m, n, p, x):
    return Dist(a**(-p), Int(ExpandTrig((d*sin(e + f*x))**n*(a - b*sin(e + f*x))**(p/S(2))*(a + b*sin(e + f*x))**(m + p/S(2)), x), x), x)


def replacement2658(a, b, d, e, f, m, n, p, x):
    return Dist(a**(-p), Int(ExpandTrig((d*cos(e + f*x))**n*(a - b*cos(e + f*x))**(p/S(2))*(a + b*cos(e + f*x))**(m + p/S(2)), x), x), x)


def replacement2659(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p, (d*sin(e + f*x))**n*(a + b*sin(e + f*x))**m, x), x)


def replacement2660(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p, (d*cos(e + f*x))**n*(a + b*cos(e + f*x))**m, x), x)


def replacement2661(a, b, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((d*sin(e + f*x))**n*(a - b*sin(e + f*x))*(a + b*sin(e + f*x))**(m + S(1)), x), x)


def replacement2662(a, b, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((d*cos(e + f*x))**n*(a - b*cos(e + f*x))*(a + b*cos(e + f*x))**(m + S(1)), x), x)


def replacement2663(a, b, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(S(2)*m + p)*(a - b*sin(e + f*x))**(-m), x), x)


def replacement2664(a, b, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(S(2)*m + p)*(a - b*cos(e + f*x))**(-m), x), x)


def replacement2665(a, b, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(S(2)*m + p)*(a - b*sin(e + f*x))**(-m), x), x)


def replacement2666(a, b, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(S(2)*m + p)*(a - b*cos(e + f*x))**(-m), x), x)


def replacement2667(a, b, e, f, g, m, p, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + p + S(1))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + p + S(1))*sin(e + f*x)), x), x) + Simp(b*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2668(a, b, e, f, g, m, p, x):
    return -Dist(S(1)/(a**S(2)*(S(2)*m + p + S(1))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**(m + S(1))*(a*m - b*(S(2)*m + p + S(1))*cos(e + f*x)), x), x) - Simp(b*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(a*f*g*(S(2)*m + p + S(1))), x)


def replacement2669(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(b*(m + p + S(2))), Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**m*(-a*(p + S(1))*sin(e + f*x) + b*(m + S(1))), x), x) - Simp((g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**(m + S(1))/(b*f*g*(m + p + S(2))), x)


def replacement2670(a, b, e, f, g, m, p, x):
    return Dist(S(1)/(b*(m + p + S(2))), Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**m*(-a*(p + S(1))*cos(e + f*x) + b*(m + S(1))), x), x) + Simp((g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**(m + S(1))/(b*f*g*(m + p + S(2))), x)


def replacement2671(a, b, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((d*sin(e + f*x))**n*(a - b*sin(e + f*x))*(a + b*sin(e + f*x))**(m + S(1)), x), x)


def replacement2672(a, b, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((d*cos(e + f*x))**n*(a - b*cos(e + f*x))*(a + b*cos(e + f*x))**(m + S(1)), x), x)


def replacement2673(a, b, d, e, f, m, n, x):
    return Dist(a**(S(-2)), Int((d*sin(e + f*x))**n*(a + b*sin(e + f*x))**(m + S(2))*(sin(e + f*x)**S(2) + S(1)), x), x) + Dist(-S(2)/(a*b*d), Int((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(2)), x), x)


def replacement2674(a, b, d, e, f, m, n, x):
    return Dist(a**(S(-2)), Int((d*cos(e + f*x))**n*(a + b*cos(e + f*x))**(m + S(2))*(cos(e + f*x)**S(2) + S(1)), x), x) + Dist(-S(2)/(a*b*d), Int((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(2)), x), x)


def replacement2675(a, b, d, e, f, m, n, x):
    return Dist(d**(S(-4)), Int((d*sin(e + f*x))**(n + S(4))*(a + b*sin(e + f*x))**m, x), x) + Int((d*sin(e + f*x))**n*(S(1) - S(2)*sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m, x)


def replacement2676(a, b, d, e, f, m, n, x):
    return Dist(d**(S(-4)), Int((d*cos(e + f*x))**(n + S(4))*(a + b*cos(e + f*x))**m, x), x) + Int((d*cos(e + f*x))**n*(S(1) - S(2)*cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m, x)


def replacement2677(a, b, d, e, f, m, n, p, x):
    return Dist(a**m*cos(e + f*x)/(f*sqrt(S(1) - sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), Subst(Int((d*x)**n*(S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2678(a, b, d, e, f, m, n, p, x):
    return -Dist(a**m*sin(e + f*x)/(f*sqrt(S(1) - cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), Subst(Int((d*x)**n*(S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2679(a, b, d, e, f, m, n, p, x):
    return Dist(a**(S(2) - p)*cos(e + f*x)/(f*sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), Subst(Int((d*x)**n*(a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2680(a, b, d, e, f, m, n, p, x):
    return -Dist(a**(S(2) - p)*sin(e + f*x)/(f*sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), Subst(Int((d*x)**n*(a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2681(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p, (d*sin(e + f*x))**n*(a + b*sin(e + f*x))**m, x), x)


def replacement2682(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p, (d*cos(e + f*x))**n*(a + b*cos(e + f*x))**m, x), x)


def replacement2683(a, b, d, e, f, g, m, n, p, x):
    return Dist(a**m*g*(g*cos(e + f*x))**(p + S(-1))*(S(1) - sin(e + f*x))**(S(1)/2 - p/S(2))*(sin(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((d*x)**n*(S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2684(a, b, d, e, f, g, m, n, p, x):
    return -Dist(a**m*g*(g*sin(e + f*x))**(p + S(-1))*(S(1) - cos(e + f*x))**(S(1)/2 - p/S(2))*(cos(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((d*x)**n*(S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2685(a, b, d, e, f, g, m, n, p, x):
    return Dist(g*(g*cos(e + f*x))**(p + S(-1))*(a - b*sin(e + f*x))**(S(1)/2 - p/S(2))*(a + b*sin(e + f*x))**(S(1)/2 - p/S(2))/f, Subst(Int((d*x)**n*(a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2686(a, b, d, e, f, g, m, n, p, x):
    return -Dist(g*(g*sin(e + f*x))**(p + S(-1))*(a - b*cos(e + f*x))**(S(1)/2 - p/S(2))*(a + b*cos(e + f*x))**(S(1)/2 - p/S(2))/f, Subst(Int((d*x)**n*(a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2687(a, b, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(S(2)*m + S(3))/(S(2)*a*(m + S(1))), Int((g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1))/sqrt(d*sin(e + f*x)), x), x) - Simp(g*sqrt(d*sin(e + f*x))*(g*cos(e + f*x))**(p + S(-1))*(a + b*sin(e + f*x))**(m + S(1))/(a*d*f*(m + S(1))), x)


def replacement2688(a, b, d, e, f, g, m, p, x):
    return Dist(g**S(2)*(S(2)*m + S(3))/(S(2)*a*(m + S(1))), Int((g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1))/sqrt(d*cos(e + f*x)), x), x) + Simp(g*sqrt(d*cos(e + f*x))*(g*sin(e + f*x))**(p + S(-1))*(a + b*cos(e + f*x))**(m + S(1))/(a*d*f*(m + S(1))), x)


def replacement2689(a, b, d, e, f, g, m, p, x):
    return Dist(S(2)*a*m/(g**S(2)*(S(2)*m + S(1))), Int((g*cos(e + f*x))**(p + S(2))*(a + b*sin(e + f*x))**(m + S(-1))/sqrt(d*sin(e + f*x)), x), x) + Simp(S(2)*sqrt(d*sin(e + f*x))*(g*cos(e + f*x))**(p + S(1))*(a + b*sin(e + f*x))**m/(d*f*g*(S(2)*m + S(1))), x)


def replacement2690(a, b, d, e, f, g, m, p, x):
    return Dist(S(2)*a*m/(g**S(2)*(S(2)*m + S(1))), Int((g*sin(e + f*x))**(p + S(2))*(a + b*cos(e + f*x))**(m + S(-1))/sqrt(d*cos(e + f*x)), x), x) + Simp(-S(2)*sqrt(d*cos(e + f*x))*(g*sin(e + f*x))**(p + S(1))*(a + b*cos(e + f*x))**m/(d*f*g*(S(2)*m + S(1))), x)


def replacement2691(a, b, d, e, f, m, n, x):
    return Int((d*sin(e + f*x))**n*(S(1) - sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m, x)


def replacement2692(a, b, d, e, f, m, n, x):
    return Int((d*cos(e + f*x))**n*(S(1) - cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m, x)


def replacement2693(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*b*d*(m + S(1))*(n + S(1))), Int((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(1))*(n + S(2)) + a*b*(m + S(1))*sin(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(1))*(n + S(3)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*sin(e + f*x)**S(2), x), x), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(a*d*f*(n + S(1))), x) - Simp((d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**(m + S(1))*(a**S(2)*(n + S(1)) - b**S(2)*(m + n + S(2)))*cos(e + f*x)/(a**S(2)*b*d**S(2)*f*(m + S(1))*(n + S(1))), x)


def replacement2694(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*b*d*(m + S(1))*(n + S(1))), Int((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(1))*(n + S(2)) + a*b*(m + S(1))*cos(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(1))*(n + S(3)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*cos(e + f*x)**S(2), x), x), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(a*d*f*(n + S(1))), x) + Simp((d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**(m + S(1))*(a**S(2)*(n + S(1)) - b**S(2)*(m + n + S(2)))*sin(e + f*x)/(a**S(2)*b*d**S(2)*f*(m + S(1))*(n + S(1))), x)


def replacement2695(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*b**S(2)*(m + S(1))*(m + S(2))), Int((d*sin(e + f*x))**n*(a + b*sin(e + f*x))**(m + S(2))*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*(m + S(2))*sin(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*sin(e + f*x)**S(2), x), x), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(a**S(2) - b**S(2))*cos(e + f*x)/(a*b**S(2)*d*f*(m + S(1))), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(2))*(a**S(2)*(-m + n + S(1)) - b**S(2)*(m + n + S(2)))*cos(e + f*x)/(a**S(2)*b**S(2)*d*f*(m + S(1))*(m + S(2))), x)


def replacement2696(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*b**S(2)*(m + S(1))*(m + S(2))), Int((d*cos(e + f*x))**n*(a + b*cos(e + f*x))**(m + S(2))*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*(m + S(2))*cos(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*cos(e + f*x)**S(2), x), x), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(a**S(2) - b**S(2))*sin(e + f*x)/(a*b**S(2)*d*f*(m + S(1))), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(2))*(a**S(2)*(-m + n + S(1)) - b**S(2)*(m + n + S(2)))*sin(e + f*x)/(a**S(2)*b**S(2)*d*f*(m + S(1))*(m + S(2))), x)


def replacement2697(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b**S(2)*(m + S(1))*(m + n + S(4))), Int((d*sin(e + f*x))**n*(a + b*sin(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*(m + S(1))*sin(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(4)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(4)))*sin(e + f*x)**S(2), x), x), x) - Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(2))*cos(e + f*x)/(b**S(2)*d*f*(m + n + S(4))), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(a**S(2) - b**S(2))*cos(e + f*x)/(a*b**S(2)*d*f*(m + S(1))), x)


def replacement2698(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b**S(2)*(m + S(1))*(m + n + S(4))), Int((d*cos(e + f*x))**n*(a + b*cos(e + f*x))**(m + S(1))*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*(m + S(1))*cos(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(4)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(4)))*cos(e + f*x)**S(2), x), x), x) + Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(2))*sin(e + f*x)/(b**S(2)*d*f*(m + n + S(4))), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(a**S(2) - b**S(2))*sin(e + f*x)/(a*b**S(2)*d*f*(m + S(1))), x)


def replacement2699(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*d**S(2)*(n + S(1))*(n + S(2))), Int((d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**m*Simp(a**S(2)*n*(n + S(2)) + a*b*m*sin(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(1))*(n + S(2)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*sin(e + f*x)**S(2), x), x), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(a*d*f*(n + S(1))), x) - Simp(b*(d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**(m + S(1))*(m + n + S(2))*cos(e + f*x)/(a**S(2)*d**S(2)*f*(n + S(1))*(n + S(2))), x)


def replacement2700(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(a**S(2)*d**S(2)*(n + S(1))*(n + S(2))), Int((d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**m*Simp(a**S(2)*n*(n + S(2)) + a*b*m*cos(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(3)) - (a**S(2)*(n + S(1))*(n + S(2)) - b**S(2)*(m + n + S(2))*(m + n + S(4)))*cos(e + f*x)**S(2), x), x), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(a*d*f*(n + S(1))), x) + Simp(b*(d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**(m + S(1))*(m + n + S(2))*sin(e + f*x)/(a**S(2)*d**S(2)*f*(n + S(1))*(n + S(2))), x)


def replacement2701(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*d*(n + S(1))*(m + n + S(4))), Int((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**m*Simp(a**S(2)*(n + S(1))*(n + S(2)) + a*b*(m + S(3))*sin(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(4)) - (a**S(2)*(n + S(1))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(4)))*sin(e + f*x)**S(2), x), x), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(a*d*f*(n + S(1))), x) - Simp((d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*d**S(2)*f*(m + n + S(4))), x)


def replacement2702(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a*b*d*(n + S(1))*(m + n + S(4))), Int((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**m*Simp(a**S(2)*(n + S(1))*(n + S(2)) + a*b*(m + S(3))*cos(e + f*x) - b**S(2)*(m + n + S(2))*(m + n + S(4)) - (a**S(2)*(n + S(1))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(4)))*cos(e + f*x)**S(2), x), x), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(a*d*f*(n + S(1))), x) + Simp((d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*d**S(2)*f*(m + n + S(4))), x)


def replacement2703(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b**S(2)*(m + n + S(3))*(m + n + S(4))), Int((d*sin(e + f*x))**n*(a + b*sin(e + f*x))**m*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*m*sin(e + f*x) - b**S(2)*(m + n + S(3))*(m + n + S(4)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(5)))*sin(e + f*x)**S(2), x), x), x) - Simp((d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*d**S(2)*f*(m + n + S(4))), x) + Simp(a*(d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*(n + S(3))*cos(e + f*x)/(b**S(2)*d*f*(m + n + S(3))*(m + n + S(4))), x)


def replacement2704(a, b, d, e, f, m, n, x):
    return -Dist(S(1)/(b**S(2)*(m + n + S(3))*(m + n + S(4))), Int((d*cos(e + f*x))**n*(a + b*cos(e + f*x))**m*Simp(a**S(2)*(n + S(1))*(n + S(3)) + a*b*m*cos(e + f*x) - b**S(2)*(m + n + S(3))*(m + n + S(4)) - (a**S(2)*(n + S(2))*(n + S(3)) - b**S(2)*(m + n + S(3))*(m + n + S(5)))*cos(e + f*x)**S(2), x), x), x) + Simp((d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*d**S(2)*f*(m + n + S(4))), x) - Simp(a*(d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*(n + S(3))*sin(e + f*x)/(b**S(2)*d*f*(m + n + S(3))*(m + n + S(4))), x)


def replacement2705(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*b**S(2)*d**S(2)*(n + S(1))*(n + S(2))*(m + n + S(5))*(m + n + S(6))), Int((d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**m*Simp(a**S(4)*(n + S(1))*(n + S(2))*(n + S(3))*(n + S(5)) - a**S(2)*b**S(2)*(n + S(2))*(S(2)*n + S(1))*(m + n + S(5))*(m + n + S(6)) + a*b*m*(a**S(2)*(n + S(1))*(n + S(2)) - b**S(2)*(m + n + S(5))*(m + n + S(6)))*sin(e + f*x) + b**S(4)*(m + n + S(2))*(m + n + S(3))*(m + n + S(5))*(m + n + S(6)) - (a**S(4)*(n + S(1))*(n + S(2))*(n + S(4))*(n + S(5)) - a**S(2)*b**S(2)*(n + S(1))*(n + S(2))*(m + n + S(5))*(S(2)*m + S(2)*n + S(13)) + b**S(4)*(m + n + S(2))*(m + n + S(4))*(m + n + S(5))*(m + n + S(6)))*sin(e + f*x)**S(2), x), x), x) + Simp((d*sin(e + f*x))**(n + S(1))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(a*d*f*(n + S(1))), x) + Simp((d*sin(e + f*x))**(n + S(4))*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*d**S(4)*f*(m + n + S(6))), x) - Simp(b*(d*sin(e + f*x))**(n + S(2))*(a + b*sin(e + f*x))**(m + S(1))*(m + n + S(2))*cos(e + f*x)/(a**S(2)*d**S(2)*f*(n + S(1))*(n + S(2))), x) - Simp(a*(d*sin(e + f*x))**(n + S(3))*(a + b*sin(e + f*x))**(m + S(1))*(n + S(5))*cos(e + f*x)/(b**S(2)*d**S(3)*f*(m + n + S(5))*(m + n + S(6))), x)


def replacement2706(a, b, d, e, f, m, n, x):
    return Dist(S(1)/(a**S(2)*b**S(2)*d**S(2)*(n + S(1))*(n + S(2))*(m + n + S(5))*(m + n + S(6))), Int((d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**m*Simp(a**S(4)*(n + S(1))*(n + S(2))*(n + S(3))*(n + S(5)) - a**S(2)*b**S(2)*(n + S(2))*(S(2)*n + S(1))*(m + n + S(5))*(m + n + S(6)) + a*b*m*(a**S(2)*(n + S(1))*(n + S(2)) - b**S(2)*(m + n + S(5))*(m + n + S(6)))*cos(e + f*x) + b**S(4)*(m + n + S(2))*(m + n + S(3))*(m + n + S(5))*(m + n + S(6)) - (a**S(4)*(n + S(1))*(n + S(2))*(n + S(4))*(n + S(5)) - a**S(2)*b**S(2)*(n + S(1))*(n + S(2))*(m + n + S(5))*(S(2)*m + S(2)*n + S(13)) + b**S(4)*(m + n + S(2))*(m + n + S(4))*(m + n + S(5))*(m + n + S(6)))*cos(e + f*x)**S(2), x), x), x) - Simp((d*cos(e + f*x))**(n + S(1))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(a*d*f*(n + S(1))), x) - Simp((d*cos(e + f*x))**(n + S(4))*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*d**S(4)*f*(m + n + S(6))), x) + Simp(b*(d*cos(e + f*x))**(n + S(2))*(a + b*cos(e + f*x))**(m + S(1))*(m + n + S(2))*sin(e + f*x)/(a**S(2)*d**S(2)*f*(n + S(1))*(n + S(2))), x) + Simp(a*(d*cos(e + f*x))**(n + S(3))*(a + b*cos(e + f*x))**(m + S(1))*(n + S(5))*sin(e + f*x)/(b**S(2)*d**S(3)*f*(m + n + S(5))*(m + n + S(6))), x)


def replacement2707(a, b, d, e, f, m, n, p, x):
    return Int(ExpandTrig((d*sin(e + f*x))**n*(S(1) - sin(e + f*x)**S(2))**(p/S(2))*(a + b*sin(e + f*x))**m, x), x)


def replacement2708(a, b, d, e, f, m, n, p, x):
    return Int(ExpandTrig((d*cos(e + f*x))**n*(S(1) - cos(e + f*x)**S(2))**(p/S(2))*(a + b*cos(e + f*x))**m, x), x)


def replacement2709(a, b, e, f, g, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p, sin(e + f*x)**n/(a + b*sin(e + f*x)), x), x)


def replacement2710(a, b, e, f, g, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p, cos(e + f*x)**n/(a + b*cos(e + f*x)), x), x)


def replacement2711(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/a, Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2)), x), x) - Dist(b*g**S(2)/(a**S(2)*d), Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**(p + S(-2)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/(a**S(2)*d**S(2)), Int((d*sin(e + f*x))**(n + S(2))*(g*cos(e + f*x))**(p + S(-2))/(a + b*sin(e + f*x)), x), x)


def replacement2712(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/a, Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2)), x), x) - Dist(b*g**S(2)/(a**S(2)*d), Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**(p + S(-2)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/(a**S(2)*d**S(2)), Int((d*cos(e + f*x))**(n + S(2))*(g*sin(e + f*x))**(p + S(-2))/(a + b*cos(e + f*x)), x), x)


def replacement2713(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/(a*b), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2))*(-a*sin(e + f*x) + b), x), x) + Dist(g**S(2)*(a**S(2) - b**S(2))/(a*b*d), Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**(p + S(-2))/(a + b*sin(e + f*x)), x), x)


def replacement2714(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/(a*b), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2))*(-a*cos(e + f*x) + b), x), x) + Dist(g**S(2)*(a**S(2) - b**S(2))/(a*b*d), Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**(p + S(-2))/(a + b*cos(e + f*x)), x), x)


def replacement2715(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/b**S(2), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2))*(a - b*sin(e + f*x)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/b**S(2), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2))/(a + b*sin(e + f*x)), x), x)


def replacement2716(a, b, d, e, f, g, n, p, x):
    return Dist(g**S(2)/b**S(2), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2))*(a - b*cos(e + f*x)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/b**S(2), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2))/(a + b*cos(e + f*x)), x), x)


def replacement2717(a, b, d, e, f, g, n, p, x):
    return Dist(a*d**S(2)/(a**S(2) - b**S(2)), Int((d*sin(e + f*x))**(n + S(-2))*(g*cos(e + f*x))**p, x), x) - Dist(b*d/(a**S(2) - b**S(2)), Int((d*sin(e + f*x))**(n + S(-1))*(g*cos(e + f*x))**p, x), x) - Dist(a**S(2)*d**S(2)/(g**S(2)*(a**S(2) - b**S(2))), Int((d*sin(e + f*x))**(n + S(-2))*(g*cos(e + f*x))**(p + S(2))/(a + b*sin(e + f*x)), x), x)


def replacement2718(a, b, d, e, f, g, n, p, x):
    return Dist(a*d**S(2)/(a**S(2) - b**S(2)), Int((d*cos(e + f*x))**(n + S(-2))*(g*sin(e + f*x))**p, x), x) - Dist(b*d/(a**S(2) - b**S(2)), Int((d*cos(e + f*x))**(n + S(-1))*(g*sin(e + f*x))**p, x), x) - Dist(a**S(2)*d**S(2)/(g**S(2)*(a**S(2) - b**S(2))), Int((d*cos(e + f*x))**(n + S(-2))*(g*sin(e + f*x))**(p + S(2))/(a + b*cos(e + f*x)), x), x)


def replacement2719(a, b, d, e, f, g, n, p, x):
    return -Dist(d/(a**S(2) - b**S(2)), Int((d*sin(e + f*x))**(n + S(-1))*(g*cos(e + f*x))**p*(-a*sin(e + f*x) + b), x), x) + Dist(a*b*d/(g**S(2)*(a**S(2) - b**S(2))), Int((d*sin(e + f*x))**(n + S(-1))*(g*cos(e + f*x))**(p + S(2))/(a + b*sin(e + f*x)), x), x)


def replacement2720(a, b, d, e, f, g, n, p, x):
    return -Dist(d/(a**S(2) - b**S(2)), Int((d*cos(e + f*x))**(n + S(-1))*(g*sin(e + f*x))**p*(-a*cos(e + f*x) + b), x), x) + Dist(a*b*d/(g**S(2)*(a**S(2) - b**S(2))), Int((d*cos(e + f*x))**(n + S(-1))*(g*sin(e + f*x))**(p + S(2))/(a + b*cos(e + f*x)), x), x)


def replacement2721(a, b, d, e, f, g, n, p, x):
    return -Dist(b**S(2)/(g**S(2)*(a**S(2) - b**S(2))), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(2))/(a + b*sin(e + f*x)), x), x) + Dist(S(1)/(a**S(2) - b**S(2)), Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**p*(a - b*sin(e + f*x)), x), x)


def replacement2722(a, b, d, e, f, g, n, p, x):
    return -Dist(b**S(2)/(g**S(2)*(a**S(2) - b**S(2))), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(2))/(a + b*cos(e + f*x)), x), x) + Dist(S(1)/(a**S(2) - b**S(2)), Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**p*(a - b*cos(e + f*x)), x), x)


def replacement2723(a, b, e, f, g, x):
    return Dist(-S(4)*sqrt(S(2))*g/f, Subst(Int(x**S(2)/(sqrt(S(1) - x**S(4)/g**S(2))*(g**S(2)*(a + b) + x**S(4)*(a - b))), x), x, sqrt(g*cos(e + f*x))/sqrt(sin(e + f*x) + S(1))), x)


def replacement2724(a, b, e, f, g, x):
    return Dist(S(4)*sqrt(S(2))*g/f, Subst(Int(x**S(2)/(sqrt(S(1) - x**S(4)/g**S(2))*(g**S(2)*(a + b) + x**S(4)*(a - b))), x), x, sqrt(g*sin(e + f*x))/sqrt(cos(e + f*x) + S(1))), x)


def replacement2725(a, b, d, e, f, g, x):
    return Dist(sqrt(sin(e + f*x))/sqrt(d*sin(e + f*x)), Int(sqrt(g*cos(e + f*x))/((a + b*sin(e + f*x))*sqrt(sin(e + f*x))), x), x)


def replacement2726(a, b, d, e, f, g, x):
    return Dist(sqrt(cos(e + f*x))/sqrt(d*cos(e + f*x)), Int(sqrt(g*sin(e + f*x))/((a + b*cos(e + f*x))*sqrt(cos(e + f*x))), x), x)


def With2727(a, b, d, e, f, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return -Dist(S(2)*sqrt(S(2))*d*(b - q)/(f*q), Subst(Int(S(1)/(sqrt(S(1) - x**S(4)/d**S(2))*(a*x**S(2) + d*(b - q))), x), x, sqrt(d*sin(e + f*x))/sqrt(cos(e + f*x) + S(1))), x) + Dist(S(2)*sqrt(S(2))*d*(b + q)/(f*q), Subst(Int(S(1)/(sqrt(S(1) - x**S(4)/d**S(2))*(a*x**S(2) + d*(b + q))), x), x, sqrt(d*sin(e + f*x))/sqrt(cos(e + f*x) + S(1))), x)


def With2728(a, b, d, e, f, x):
    q = Rt(-a**S(2) + b**S(2), S(2))
    return Dist(S(2)*sqrt(S(2))*d*(b - q)/(f*q), Subst(Int(S(1)/(sqrt(S(1) - x**S(4)/d**S(2))*(a*x**S(2) + d*(b - q))), x), x, sqrt(d*cos(e + f*x))/sqrt(sin(e + f*x) + S(1))), x) + Dist(-S(2)*sqrt(S(2))*d*(b + q)/(f*q), Subst(Int(S(1)/(sqrt(S(1) - x**S(4)/d**S(2))*(a*x**S(2) + d*(b + q))), x), x, sqrt(d*cos(e + f*x))/sqrt(sin(e + f*x) + S(1))), x)


def replacement2729(a, b, d, e, f, g, x):
    return Dist(sqrt(cos(e + f*x))/sqrt(g*cos(e + f*x)), Int(sqrt(d*sin(e + f*x))/((a + b*sin(e + f*x))*sqrt(cos(e + f*x))), x), x)


def replacement2730(a, b, d, e, f, g, x):
    return Dist(sqrt(sin(e + f*x))/sqrt(g*sin(e + f*x)), Int(sqrt(d*cos(e + f*x))/((a + b*cos(e + f*x))*sqrt(sin(e + f*x))), x), x)


def replacement2731(a, b, d, e, f, g, n, p, x):
    return Dist(d/b, Int((d*sin(e + f*x))**(n + S(-1))*(g*cos(e + f*x))**p, x), x) - Dist(a*d/b, Int((d*sin(e + f*x))**(n + S(-1))*(g*cos(e + f*x))**p/(a + b*sin(e + f*x)), x), x)


def replacement2732(a, b, d, e, f, g, n, p, x):
    return Dist(d/b, Int((d*cos(e + f*x))**(n + S(-1))*(g*sin(e + f*x))**p, x), x) - Dist(a*d/b, Int((d*cos(e + f*x))**(n + S(-1))*(g*sin(e + f*x))**p/(a + b*cos(e + f*x)), x), x)


def replacement2733(a, b, d, e, f, g, n, p, x):
    return Dist(S(1)/a, Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**p, x), x) - Dist(b/(a*d), Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**p/(a + b*sin(e + f*x)), x), x)


def replacement2734(a, b, d, e, f, g, n, p, x):
    return Dist(S(1)/a, Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**p, x), x) - Dist(b/(a*d), Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**p/(a + b*cos(e + f*x)), x), x)


def replacement2735(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p, (d*sin(e + f*x))**n*(a + b*sin(e + f*x))**m, x), x)


def replacement2736(a, b, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p, (d*cos(e + f*x))**n*(a + b*cos(e + f*x))**m, x), x)


def replacement2737(a, b, d, e, f, g, m, n, p, x):
    return Dist(g**S(2)/a, Int((d*sin(e + f*x))**n*(g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1)), x), x) - Dist(b*g**S(2)/(a**S(2)*d), Int((d*sin(e + f*x))**(n + S(1))*(g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**(m + S(1)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/(a**S(2)*d**S(2)), Int((d*sin(e + f*x))**(n + S(2))*(g*cos(e + f*x))**(p + S(-2))*(a + b*sin(e + f*x))**m, x), x)


def replacement2738(a, b, d, e, f, g, m, n, p, x):
    return Dist(g**S(2)/a, Int((d*cos(e + f*x))**n*(g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Dist(b*g**S(2)/(a**S(2)*d), Int((d*cos(e + f*x))**(n + S(1))*(g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**(m + S(1)), x), x) - Dist(g**S(2)*(a**S(2) - b**S(2))/(a**S(2)*d**S(2)), Int((d*cos(e + f*x))**(n + S(2))*(g*sin(e + f*x))**(p + S(-2))*(a + b*cos(e + f*x))**m, x), x)


def replacement2739(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**(S(2)*m), Int((a - b*sin(e + f*x))**(-m)*(c + d*sin(e + f*x))**n, x), x)


def replacement2740(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**(S(2)*m), Int((a - b*cos(e + f*x))**(-m)*(c + d*cos(e + f*x))**n, x), x)


def replacement2741(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((g*cos(e + f*x))**(S(2)*m + p)*(a - b*sin(e + f*x))**(-m)*(c + d*sin(e + f*x))**n, x), x)


def replacement2742(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((a/g)**(S(2)*m), Int((g*sin(e + f*x))**(S(2)*m + p)*(a - b*cos(e + f*x))**(-m)*(c + d*cos(e + f*x))**n, x), x)


def replacement2743(a, b, c, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((a - b*sin(e + f*x))*(a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x)


def replacement2744(a, b, c, d, e, f, m, n, x):
    return Dist(b**(S(-2)), Int((a - b*cos(e + f*x))*(a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x)


def replacement2745(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**m*cos(e + f*x)/(f*sqrt(S(1) - sin(e + f*x))*sqrt(sin(e + f*x) + S(1))), Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, sin(e + f*x)), x)


def replacement2746(a, b, c, d, e, f, m, n, p, x):
    return -Dist(a**m*sin(e + f*x)/(f*sqrt(S(1) - cos(e + f*x))*sqrt(cos(e + f*x) + S(1))), Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, cos(e + f*x)), x)


def replacement2747(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**(S(2) - p)*cos(e + f*x)/(f*sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, sin(e + f*x)), x)


def replacement2748(a, b, c, d, e, f, m, n, p, x):
    return -Dist(a**(S(2) - p)*sin(e + f*x)/(f*sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, cos(e + f*x)), x)


def replacement2749(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p, (a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2750(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p, (a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2751(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(a**m*g*(g*cos(e + f*x))**(p + S(-1))*(S(1) - sin(e + f*x))**(S(1)/2 - p/S(2))*(sin(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, sin(e + f*x)), x)


def replacement2752(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(a**m*g*(g*sin(e + f*x))**(p + S(-1))*(S(1) - cos(e + f*x))**(S(1)/2 - p/S(2))*(cos(e + f*x) + S(1))**(S(1)/2 - p/S(2))/f, Subst(Int((S(1) - b*x/a)**(p/S(2) + S(-1)/2)*(S(1) + b*x/a)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, cos(e + f*x)), x)


def replacement2753(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g*(g*cos(e + f*x))**(p + S(-1))*(a - b*sin(e + f*x))**(S(1)/2 - p/S(2))*(a + b*sin(e + f*x))**(S(1)/2 - p/S(2))/f, Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, sin(e + f*x)), x)


def replacement2754(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(g*(g*sin(e + f*x))**(p + S(-1))*(a - b*cos(e + f*x))**(S(1)/2 - p/S(2))*(a + b*cos(e + f*x))**(S(1)/2 - p/S(2))/f, Subst(Int((a - b*x)**(p/S(2) + S(-1)/2)*(a + b*x)**(m + p/S(2) + S(-1)/2)*(c + d*x)**n, x), x, cos(e + f*x)), x)


def replacement2755(a, b, c, d, e, f, m, n, x):
    return Int((S(1) - sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement2756(a, b, c, d, e, f, m, n, x):
    return Int((S(1) - cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement2757(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandTrig((S(1) - sin(e + f*x)**S(2))**(p/S(2))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2758(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandTrig((S(1) - cos(e + f*x)**S(2))**(p/S(2))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2759(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2760(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2761(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g*cos(e + f*x))**p*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement2762(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g*sin(e + f*x))**p*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement2763(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/cos(e + f*x))**FracPart(p)*(g*cos(e + f*x))**FracPart(p), Int((g*cos(e + f*x))**(-p)*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2764(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(S(2)*IntPart(p))*(g/sin(e + f*x))**FracPart(p)*(g*sin(e + f*x))**FracPart(p), Int((g*sin(e + f*x))**(-p)*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2765(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int(sqrt(a + b*sin(e + f*x))/sqrt(g*sin(e + f*x)), x), x) - Dist(c*g/d, Int(sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2766(a, b, c, d, e, f, g, x):
    return Dist(g/d, Int(sqrt(a + b*cos(e + f*x))/sqrt(g*cos(e + f*x)), x), x) - Dist(c*g/d, Int(sqrt(a + b*cos(e + f*x))/(sqrt(g*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2767(a, b, c, d, e, f, g, x):
    return Dist(b/d, Int(sqrt(g*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(sqrt(g*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2768(a, b, c, d, e, f, g, x):
    return Dist(b/d, Int(sqrt(g*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x) - Dist((-a*d + b*c)/d, Int(sqrt(g*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2769(a, b, c, d, e, f, g, x):
    return Dist(-S(2)*b/f, Subst(Int(S(1)/(a*d + b*c + c*g*x**S(2)), x), x, b*cos(e + f*x)/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x)))), x)


def replacement2770(a, b, c, d, e, f, g, x):
    return Dist(S(2)*b/f, Subst(Int(S(1)/(a*d + b*c + c*g*x**S(2)), x), x, b*sin(e + f*x)/(sqrt(g*cos(e + f*x))*sqrt(a + b*cos(e + f*x)))), x)


def replacement2771(a, b, c, d, e, f, x):
    return -Simp(sqrt(a + b)*EllipticE(asin(cos(e + f*x)/(sin(e + f*x) + S(1))), -(a - b)/(a + b))/(c*f), x)


def replacement2772(a, b, c, d, e, f, x):
    return Simp(sqrt(a + b)*EllipticE(asin(sin(e + f*x)/(cos(e + f*x) + S(1))), -(a - b)/(a + b))/(c*f), x)


def replacement2773(a, b, c, d, e, f, g, x):
    return -Simp(sqrt(d*sin(e + f*x)/(c + d*sin(e + f*x)))*sqrt(a + b*sin(e + f*x))*EllipticE(asin(c*cos(e + f*x)/(c + d*sin(e + f*x))), (-a*d + b*c)/(a*d + b*c))/(d*f*sqrt(g*sin(e + f*x))*sqrt(c**S(2)*(a + b*sin(e + f*x))/((c + d*sin(e + f*x))*(a*c + b*d)))), x)


def replacement2774(a, b, c, d, e, f, g, x):
    return Simp(sqrt(d*cos(e + f*x)/(c + d*cos(e + f*x)))*sqrt(a + b*cos(e + f*x))*EllipticE(asin(c*sin(e + f*x)/(c + d*cos(e + f*x))), (-a*d + b*c)/(a*d + b*c))/(d*f*sqrt(g*cos(e + f*x))*sqrt(c**S(2)*(a + b*cos(e + f*x))/((c + d*cos(e + f*x))*(a*c + b*d)))), x)


def replacement2775(a, b, c, d, e, f, g, x):
    return Dist(a/c, Int(S(1)/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), x), x) + Dist((-a*d + b*c)/(c*g), Int(sqrt(g*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2776(a, b, c, d, e, f, g, x):
    return Dist(a/c, Int(S(1)/(sqrt(g*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), x), x) + Dist((-a*d + b*c)/(c*g), Int(sqrt(g*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2777(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b*sin(e + f*x))/sin(e + f*x), x), x) - Dist(d/c, Int(sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x)), x), x)


def replacement2778(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b*cos(e + f*x))/cos(e + f*x), x), x) - Dist(d/c, Int(sqrt(a + b*cos(e + f*x))/(c + d*cos(e + f*x)), x), x)


def replacement2779(a, b, c, d, e, f, x):
    return Dist(a/c, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2780(a, b, c, d, e, f, x):
    return Dist(a/c, Int(S(1)/(sqrt(a + b*cos(e + f*x))*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2781(a, b, c, d, e, f, g, x):
    return -Dist(a*g/(-a*d + b*c), Int(S(1)/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), x), x) + Dist(c*g/(-a*d + b*c), Int(sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2782(a, b, c, d, e, f, g, x):
    return -Dist(a*g/(-a*d + b*c), Int(S(1)/(sqrt(g*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), x), x) + Dist(c*g/(-a*d + b*c), Int(sqrt(a + b*cos(e + f*x))/(sqrt(g*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2783(a, b, c, d, e, f, g, x):
    return Simp(S(2)*sqrt(-S(1)/tan(e + f*x)**S(2))*sqrt(g*sin(e + f*x))*sqrt((a/sin(e + f*x) + b)/(a + b))*EllipticPi(S(2)*c/(c + d), asin(sqrt(S(2))*sqrt(S(1) - S(1)/sin(e + f*x))/S(2)), S(2)*a/(a + b))*tan(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(c + d)), x)


def replacement2784(a, b, c, d, e, f, g, x):
    return Simp(-S(2)*sqrt(-tan(e + f*x)**S(2))*sqrt(g*cos(e + f*x))*sqrt((a/cos(e + f*x) + b)/(a + b))*EllipticPi(S(2)*c/(c + d), asin(sqrt(S(2))*sqrt(S(1) - S(1)/cos(e + f*x))/S(2)), S(2)*a/(a + b))/(f*sqrt(a + b*cos(e + f*x))*(c + d)*tan(e + f*x)), x)


def replacement2785(a, b, c, d, e, f, g, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b*sin(e + f*x))/(sqrt(g*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2786(a, b, c, d, e, f, g, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(sqrt(g*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), x), x) - Dist(d/(-a*d + b*c), Int(sqrt(a + b*cos(e + f*x))/(sqrt(g*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2787(a, b, c, d, e, f, g, x):
    return Dist(S(1)/c, Int(S(1)/(sqrt(g*sin(e + f*x))*sqrt(a + b*sin(e + f*x))), x), x) - Dist(d/(c*g), Int(sqrt(g*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2788(a, b, c, d, e, f, g, x):
    return Dist(S(1)/c, Int(S(1)/(sqrt(g*cos(e + f*x))*sqrt(a + b*cos(e + f*x))), x), x) - Dist(d/(c*g), Int(sqrt(g*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2789(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*(-a*d + b*c)), Int((-a*d + b*c - b*d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sin(e + f*x)), x), x) + Dist(d**S(2)/(c*(-a*d + b*c)), Int(sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x)), x), x)


def replacement2790(a, b, c, d, e, f, x):
    return Dist(S(1)/(c*(-a*d + b*c)), Int((-a*d + b*c - b*d*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*cos(e + f*x)), x), x) + Dist(d**S(2)/(c*(-a*d + b*c)), Int(sqrt(a + b*cos(e + f*x))/(c + d*cos(e + f*x)), x), x)


def replacement2791(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sin(e + f*x)), x), x) - Dist(d/c, Int(S(1)/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x)


def replacement2792(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(S(1)/(sqrt(a + b*cos(e + f*x))*cos(e + f*x)), x), x) - Dist(d/c, Int(S(1)/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x)


def replacement2793(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/sin(e + f*x), x), x) - Dist(d/c, Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x)


def replacement2794(a, b, c, d, e, f, x):
    return Dist(S(1)/c, Int(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))/cos(e + f*x), x), x) - Dist(d/c, Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x)


def replacement2795(a, b, c, d, e, f, x):
    return Dist(-S(2)*a/f, Subst(Int(S(1)/(-a*c*x**S(2) + S(1)), x), x, cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x)))), x)


def replacement2796(a, b, c, d, e, f, x):
    return Dist(S(2)*a/f, Subst(Int(S(1)/(-a*c*x**S(2) + S(1)), x), x, sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x)))), x)


def replacement2797(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sin(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2798(a, b, c, d, e, f, x):
    return Dist(a/c, Int(sqrt(c + d*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*cos(e + f*x)), x), x) + Dist((-a*d + b*c)/c, Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2799(a, b, c, d, e, f, x):
    return Simp(-S(2)*sqrt((-a*d + b*c)*(sin(e + f*x) + S(1))/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(S(1) - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*EllipticPi(a*(c + d)/(c*(a + b)), asin(sqrt(c + d*sin(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*sin(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(c*f*Rt((a + b)/(c + d), S(2))*cos(e + f*x)), x)


def replacement2800(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt((-a*d + b*c)*(cos(e + f*x) + S(1))/((a + b*cos(e + f*x))*(c - d)))*sqrt(-(S(1) - cos(e + f*x))*(-a*d + b*c)/((a + b*cos(e + f*x))*(c + d)))*(a + b*cos(e + f*x))*EllipticPi(a*(c + d)/(c*(a + b)), asin(sqrt(c + d*cos(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*cos(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(c*f*Rt((a + b)/(c + d), S(2))*sin(e + f*x)), x)


def replacement2801(a, b, c, d, e, f, x):
    return Dist(cos(e + f*x)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), Int(S(1)/(sin(e + f*x)*cos(e + f*x)), x), x)


def replacement2802(a, b, c, d, e, f, x):
    return Dist(sin(e + f*x)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), Int(S(1)/(sin(e + f*x)*cos(e + f*x)), x), x)


def replacement2803(a, b, c, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b*sin(e + f*x))/(sqrt(c + d*sin(e + f*x))*sin(e + f*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2804(a, b, c, d, e, f, x):
    return Dist(S(1)/a, Int(sqrt(a + b*cos(e + f*x))/(sqrt(c + d*cos(e + f*x))*cos(e + f*x)), x), x) - Dist(b/a, Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2805(a, b, c, d, e, f, x):
    return Dist(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/cos(e + f*x), Int(S(1)/tan(e + f*x), x), x)


def replacement2806(a, b, c, d, e, f, x):
    return Dist(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))/sin(e + f*x), Int(tan(e + f*x), x), x)


def replacement2807(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b*sin(e + f*x))/(sqrt(c + d*sin(e + f*x))*sin(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x)


def replacement2808(a, b, c, d, e, f, x):
    return Dist(c, Int(sqrt(a + b*cos(e + f*x))/(sqrt(c + d*cos(e + f*x))*cos(e + f*x)), x), x) + Dist(d, Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x)


def replacement2809(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**n*c**n, Int((a + b*sin(e + f*x))**(m - n)*tan(e + f*x)**p, x), x)


def replacement2810(a, b, c, d, e, f, m, n, p, x):
    return Dist(a**n*c**n, Int((a + b*cos(e + f*x))**(m - n)*(S(1)/tan(e + f*x))**p, x), x)


def replacement2811(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(sqrt(a - b*sin(e + f*x))*sqrt(a + b*sin(e + f*x))/(f*cos(e + f*x)), Subst(Int((g*x)**p*(a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, sin(e + f*x)), x)


def replacement2812(a, b, c, d, e, f, g, m, n, p, x):
    return -Dist(sqrt(a - b*cos(e + f*x))*sqrt(a + b*cos(e + f*x))/(f*sin(e + f*x)), Subst(Int((g*x)**p*(a + b*x)**(m + S(-1)/2)*(c + d*x)**n/sqrt(a - b*x), x), x, cos(e + f*x)), x)


def replacement2813(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*sin(e + f*x))**p*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2814(a, b, c, d, e, f, g, m, n, p, x):
    return Int(ExpandTrig((g*cos(e + f*x))**p*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2815(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g*sin(e + f*x))**p*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement2816(a, b, c, d, e, f, g, m, n, p, x):
    return Int((g*cos(e + f*x))**p*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement2817(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(m + n), Int((g*sin(e + f*x))**(-m - n + p)*(a*sin(e + f*x) + b)**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement2818(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(m + n), Int((g*cos(e + f*x))**(-m - n + p)*(a*cos(e + f*x) + b)**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement2819(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/sin(e + f*x))**p*(g*sin(e + f*x))**p, Int((g/sin(e + f*x))**(-p)*(a + b/sin(e + f*x))**m*(c + d/sin(e + f*x))**n, x), x)


def replacement2820(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/cos(e + f*x))**p*(g*cos(e + f*x))**p, Int((g/cos(e + f*x))**(-p)*(a + b/cos(e + f*x))**m*(c + d/cos(e + f*x))**n, x), x)


def replacement2821(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**n, Int((g*sin(e + f*x))**(-n + p)*(a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement2822(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**n, Int((g*cos(e + f*x))**(-n + p)*(a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement2823(a, b, c, d, e, f, m, n, p, x):
    return Int((c + d/sin(e + f*x))**n*(a/sin(e + f*x) + b)**m*(S(1)/sin(e + f*x))**(-m - p), x)


def replacement2824(a, b, c, d, e, f, m, n, p, x):
    return Int((c + d/cos(e + f*x))**n*(a/cos(e + f*x) + b)**m*(S(1)/cos(e + f*x))**(-m - p), x)


def replacement2825(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g*sin(e + f*x))**p*(S(1)/sin(e + f*x))**p, Int((c + d/sin(e + f*x))**n*(a/sin(e + f*x) + b)**m*(S(1)/sin(e + f*x))**(-m - p), x), x)


def replacement2826(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g*cos(e + f*x))**p*(S(1)/cos(e + f*x))**p, Int((c + d/cos(e + f*x))**n*(a/cos(e + f*x) + b)**m*(S(1)/cos(e + f*x))**(-m - p), x), x)


def replacement2827(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g*sin(e + f*x))**n*(c + d/sin(e + f*x))**n*(c*sin(e + f*x) + d)**(-n), Int((g*sin(e + f*x))**(-n + p)*(a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n, x), x)


def replacement2828(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g*cos(e + f*x))**n*(c + d/cos(e + f*x))**n*(c*cos(e + f*x) + d)**(-n), Int((g*cos(e + f*x))**(-n + p)*(a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n, x), x)


def replacement2829(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(m + n), Int((g/sin(e + f*x))**(-m - n + p)*(a/sin(e + f*x) + b)**m*(c/sin(e + f*x) + d)**n, x), x)


def replacement2830(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**(m + n), Int((g/cos(e + f*x))**(-m - n + p)*(a/cos(e + f*x) + b)**m*(c/cos(e + f*x) + d)**n, x), x)


def replacement2831(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/sin(e + f*x))**p*(g*sin(e + f*x))**p, Int((g*sin(e + f*x))**(-p)*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2832(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/cos(e + f*x))**p*(g*cos(e + f*x))**p, Int((g*cos(e + f*x))**(-p)*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2833(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**m, Int((g/sin(e + f*x))**(-m + p)*(c + d/sin(e + f*x))**n*(a/sin(e + f*x) + b)**m, x), x)


def replacement2834(a, b, c, d, e, f, g, m, n, p, x):
    return Dist(g**m, Int((g/cos(e + f*x))**(-m + p)*(c + d/cos(e + f*x))**n*(a/cos(e + f*x) + b)**m, x), x)


def replacement2835(a, b, c, d, e, f, m, n, p, x):
    return Int((a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-n - p), x)


def replacement2836(a, b, c, d, e, f, m, n, p, x):
    return Int((a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-n - p), x)


def replacement2837(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/sin(e + f*x))**p*sin(e + f*x)**p, Int((a + b*sin(e + f*x))**m*(c*sin(e + f*x) + d)**n*sin(e + f*x)**(-n - p), x), x)


def replacement2838(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/cos(e + f*x))**p*cos(e + f*x)**p, Int((a + b*cos(e + f*x))**m*(c*cos(e + f*x) + d)**n*cos(e + f*x)**(-n - p), x), x)


def replacement2839(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/sin(e + f*x))**m*(a + b*sin(e + f*x))**m*(a/sin(e + f*x) + b)**(-m), Int((g/sin(e + f*x))**(-m + p)*(c + d/sin(e + f*x))**n*(a/sin(e + f*x) + b)**m, x), x)


def replacement2840(a, b, c, d, e, f, g, m, n, p, x):
    return Dist((g/cos(e + f*x))**m*(a + b*cos(e + f*x))**m*(a/cos(e + f*x) + b)**(-m), Int((g/cos(e + f*x))**(-m + p)*(c + d/cos(e + f*x))**n*(a/cos(e + f*x) + b)**m, x), x)


def replacement2841(A, B, a, b, e, f, m, n, x):
    return Int(ExpandTrig((A + B*sin(e + f*x))*(a + b*sin(e + f*x))**m*sin(e + f*x)**n, x), x)


def replacement2842(A, B, a, b, e, f, m, n, x):
    return Int(ExpandTrig((A + B*cos(e + f*x))*(a + b*cos(e + f*x))**m*cos(e + f*x)**n, x), x)


def replacement2843(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(a**m*c**m, Int((A + B*sin(e + f*x))*(c + d*sin(e + f*x))**(-m + n)*cos(e + f*x)**(S(2)*m), x), x)


def replacement2844(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(a**m*c**m, Int((A + B*cos(e + f*x))*(c + d*cos(e + f*x))**(-m + n)*sin(e + f*x)**(S(2)*m), x), x)


def replacement2845(A, B, a, b, c, d, e, f, m, x):
    return Int((a + b*sin(e + f*x))**m*(A*c + B*d*sin(e + f*x)**S(2) + (A*d + B*c)*sin(e + f*x)), x)


def replacement2846(A, B, a, b, c, d, e, f, m, x):
    return Int((a + b*cos(e + f*x))**m*(A*c + B*d*cos(e + f*x)**S(2) + (A*d + B*c)*cos(e + f*x)), x)


def replacement2847(A, B, a, b, c, d, e, f, x):
    return Dist((A*b + B*a)/(S(2)*a*b), Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x) + Dist((A*d + B*c)/(S(2)*c*d), Int(sqrt(c + d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x)


def replacement2848(A, B, a, b, c, d, e, f, x):
    return Dist((A*b + B*a)/(S(2)*a*b), Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x) + Dist((A*d + B*c)/(S(2)*c*d), Int(sqrt(c + d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x)


def replacement2849(A, B, a, b, c, d, e, f, m, n, x):
    return -Simp(B*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(m + n + S(1))), x)


def replacement2850(A, B, a, b, c, d, e, f, m, n, x):
    return Simp(B*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(m + n + S(1))), x)


def replacement2851(A, B, a, b, c, d, e, f, n, x):
    return Dist(B/d, Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(n + S(1)), x), x) - Dist((-A*d + B*c)/d, Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**n, x), x)


def replacement2852(A, B, a, b, c, d, e, f, n, x):
    return Dist(B/d, Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**(n + S(1)), x), x) - Dist((-A*d + B*c)/d, Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**n, x), x)


def replacement2853(A, B, a, b, c, d, e, f, m, n, x):
    return Dist((A*b*(m + n + S(1)) + B*a*(m - n))/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*(A*b - B*a)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2854(A, B, a, b, c, d, e, f, m, n, x):
    return Dist((A*b*(m + n + S(1)) + B*a*(m - n))/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*(A*b - B*a)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2855(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist((-A*d*(m + n + S(1)) + B*c*(m - n))/(d*(m + n + S(1))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x) - Simp(B*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(m + n + S(1))), x)


def replacement2856(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist((-A*d*(m + n + S(1)) + B*c*(m - n))/(d*(m + n + S(1))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x) + Simp(B*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(m + n + S(1))), x)


def replacement2857(A, B, a, b, c, d, e, f, m, n, x):
    return Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(-A*d + B*c)*cos(e + f*x)/(f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2858(A, B, a, b, c, d, e, f, m, n, x):
    return -Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(-A*d + B*c)*sin(e + f*x)/(f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2859(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist(b/(d*(n + S(1))*(a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*a*d*(m - n + S(-2)) - B*(a*c*(m + S(-1)) + b*d*(n + S(1))) - (A*b*d*(m + n + S(1)) - B*(-a*d*(n + S(1)) + b*c*m))*sin(e + f*x), x), x), x) - Simp(b**S(2)*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*(-A*d + B*c)*cos(e + f*x)/(d*f*(n + S(1))*(a*d + b*c)), x)


def replacement2860(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist(b/(d*(n + S(1))*(a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*a*d*(m - n + S(-2)) - B*(a*c*(m + S(-1)) + b*d*(n + S(1))) - (A*b*d*(m + n + S(1)) - B*(-a*d*(n + S(1)) + b*c*m))*cos(e + f*x), x), x), x) + Simp(b**S(2)*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*(-A*d + B*c)*sin(e + f*x)/(d*f*(n + S(1))*(a*d + b*c)), x)


def replacement2861(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) + B*(a*c*(m + S(-1)) + b*d*(n + S(1))) + (A*b*d*(m + n + S(1)) - B*(-a*d*(S(2)*m + n) + b*c*m))*sin(e + f*x), x), x), x) - Simp(B*b*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(1))), x)


def replacement2862(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*Simp(A*a*d*(m + n + S(1)) + B*(a*c*(m + S(-1)) + b*d*(n + S(1))) + (A*b*d*(m + n + S(1)) - B*(-a*d*(S(2)*m + n) + b*c*m))*cos(e + f*x), x), x), x) + Simp(B*b*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(1))), x)


def replacement2863(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-1))*Simp(A*(a*d*n - b*c*(m + S(1))) - B*(a*c*m + b*d*n) - d*(A*b*(m + n + S(1)) + B*a*(m - n))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*(A*b - B*a)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2864(A, B, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(a*b*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-1))*Simp(A*(a*d*n - b*c*(m + S(1))) - B*(a*c*m + b*d*n) - d*(A*b*(m + n + S(1)) + B*a*(m - n))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*(A*b - B*a)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2865(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(A*(-a*d*(S(2)*m + n + S(2)) + b*c*(m + S(1))) + B*(a*c*m + b*d*(n + S(1))) + d*(A*b - B*a)*(m + n + S(2))*sin(e + f*x), x), x), x) + Simp(b*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*b - B*a)*cos(e + f*x)/(a*f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2866(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(a*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(A*(-a*d*(S(2)*m + n + S(2)) + b*c*(m + S(1))) + B*(a*c*m + b*d*(n + S(1))) + d*(A*b - B*a)*(m + n + S(2))*cos(e + f*x), x), x), x) - Simp(b*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*b - B*a)*sin(e + f*x)/(a*f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2867(A, B, a, b, c, d, e, f, n, x):
    return Simp(-S(2)*B*b*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*sqrt(a + b*sin(e + f*x))*(S(2)*n + S(3))), x)


def replacement2868(A, B, a, b, c, d, e, f, n, x):
    return Simp(S(2)*B*b*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*sqrt(a + b*cos(e + f*x))*(S(2)*n + S(3))), x)


def replacement2869(A, B, a, b, c, d, e, f, n, x):
    return Dist((A*b*d*(S(2)*n + S(3)) - B*(-S(2)*a*d*(n + S(1)) + b*c))/(S(2)*d*(n + S(1))*(a*d + b*c)), Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**(n + S(1)), x), x) - Simp(b**S(2)*(c + d*sin(e + f*x))**(n + S(1))*(-A*d + B*c)*cos(e + f*x)/(d*f*sqrt(a + b*sin(e + f*x))*(n + S(1))*(a*d + b*c)), x)


def replacement2870(A, B, a, b, c, d, e, f, n, x):
    return Dist((A*b*d*(S(2)*n + S(3)) - B*(-S(2)*a*d*(n + S(1)) + b*c))/(S(2)*d*(n + S(1))*(a*d + b*c)), Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**(n + S(1)), x), x) + Simp(b**S(2)*(c + d*cos(e + f*x))**(n + S(1))*(-A*d + B*c)*sin(e + f*x)/(d*f*sqrt(a + b*cos(e + f*x))*(n + S(1))*(a*d + b*c)), x)


def replacement2871(A, B, a, b, c, d, e, f, n, x):
    return Dist((A*b*d*(S(2)*n + S(3)) - B*(-S(2)*a*d*(n + S(1)) + b*c))/(b*d*(S(2)*n + S(3))), Int(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**n, x), x) + Simp(-S(2)*B*b*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*sqrt(a + b*sin(e + f*x))*(S(2)*n + S(3))), x)


def replacement2872(A, B, a, b, c, d, e, f, n, x):
    return Dist((A*b*d*(S(2)*n + S(3)) - B*(-S(2)*a*d*(n + S(1)) + b*c))/(b*d*(S(2)*n + S(3))), Int(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**n, x), x) + Simp(S(2)*B*b*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*sqrt(a + b*cos(e + f*x))*(S(2)*n + S(3))), x)


def replacement2873(A, B, a, b, c, d, e, f, x):
    return Dist(B/b, Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x) + Dist((A*b - B*a)/b, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2874(A, B, a, b, c, d, e, f, x):
    return Dist(B/b, Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x) + Dist((A*b - B*a)/b, Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2875(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(-1))*Simp(A*b*c*(m + n + S(1)) + B*(a*c*m + b*d*n) + (A*b*d*(m + n + S(1)) + B*(a*d*m + b*c*n))*sin(e + f*x), x), x), x) - Simp(B*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(m + n + S(1))), x)


def replacement2876(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(m + n + S(1))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(-1))*Simp(A*b*c*(m + n + S(1)) + B*(a*c*m + b*d*n) + (A*b*d*(m + n + S(1)) + B*(a*d*m + b*c*n))*cos(e + f*x), x), x), x) + Simp(B*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(m + n + S(1))), x)


def replacement2877(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*(a*d*m + b*c*(n + S(1))) - B*(a*c*m + b*d*(n + S(1))) + b*(-A*d + B*c)*(m + n + S(2))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(-A*d + B*c)*cos(e + f*x)/(f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2878(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*(a*d*m + b*c*(n + S(1))) - B*(a*c*m + b*d*(n + S(1))) + b*(-A*d + B*c)*(m + n + S(2))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(-A*d + B*c)*sin(e + f*x)/(f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2879(A, B, a, b, c, d, e, f, x):
    return Dist((A*b - B*a)/(-a*d + b*c), Int(S(1)/sqrt(a + b*sin(e + f*x)), x), x) + Dist((-A*d + B*c)/(-a*d + b*c), Int(sqrt(a + b*sin(e + f*x))/(c + d*sin(e + f*x)), x), x)


def replacement2880(A, B, a, b, c, d, e, f, x):
    return Dist((A*b - B*a)/(-a*d + b*c), Int(S(1)/sqrt(a + b*cos(e + f*x)), x), x) + Dist((-A*d + B*c)/(-a*d + b*c), Int(sqrt(a + b*cos(e + f*x))/(c + d*cos(e + f*x)), x), x)


def replacement2881(A, B, a, b, c, d, e, f, m, x):
    return Dist(B/d, Int((a + b*sin(e + f*x))**m, x), x) - Dist((-A*d + B*c)/d, Int((a + b*sin(e + f*x))**m/(c + d*sin(e + f*x)), x), x)


def replacement2882(A, B, a, b, c, d, e, f, m, x):
    return Dist(B/d, Int((a + b*cos(e + f*x))**m, x), x) - Dist((-A*d + B*c)/d, Int((a + b*cos(e + f*x))**m/(c + d*cos(e + f*x)), x), x)


def replacement2883(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(B/b, Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n, x), x) + Dist((A*b - B*a)/b, Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x), x)


def replacement2884(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(B/b, Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n, x), x) + Dist((A*b - B*a)/b, Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x), x)


def replacement2885(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**(n + S(1))*Simp(a*d*(n + S(1))*(A*a*c + B*b*c - d*(A*b + B*a)) + b*(m + S(-1))*(-A*d + B*c)*(-a*d + b*c) + b*(-B*b*(c**S(2)*m + d**S(2)*(n + S(1))) + d*(m + n + S(1))*(-A*a*d + A*b*c + B*a*c))*sin(e + f*x)**S(2) + (-a*(n + S(2))*(-A*d + B*c)*(-a*d + b*c) + b*(n + S(1))*(a*(A*c*d + B*(c**S(2) - S(2)*d**S(2))) + b*d*(-A*d + B*c)))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*(-A*d + B*c)*(-a*d + b*c)*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2886(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**(n + S(1))*Simp(a*d*(n + S(1))*(A*a*c + B*b*c - d*(A*b + B*a)) + b*(m + S(-1))*(-A*d + B*c)*(-a*d + b*c) + b*(-B*b*(c**S(2)*m + d**S(2)*(n + S(1))) + d*(m + n + S(1))*(-A*a*d + A*b*c + B*a*c))*cos(e + f*x)**S(2) + (-a*(n + S(2))*(-A*d + B*c)*(-a*d + b*c) + b*(n + S(1))*(a*(A*c*d + B*(c**S(2) - S(2)*d**S(2))) + b*d*(-A*d + B*c)))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*(-A*d + B*c)*(-a*d + b*c)*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2887(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-2))*(c + d*sin(e + f*x))**n*Simp(A*a**S(2)*d*(m + n + S(1)) + B*b*(a*d*(n + S(1)) + b*c*(m + S(-1))) + b*(A*b*d*(m + n + S(1)) - B*(-a*d*(S(2)*m + n) + b*c*m))*sin(e + f*x)**S(2) + (-B*b*(a*c - b*d*(m + n)) + a*d*(S(2)*A*b + B*a)*(m + n + S(1)))*sin(e + f*x), x), x), x) - Simp(B*b*(a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(1))), x)


def replacement2888(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-2))*(c + d*cos(e + f*x))**n*Simp(A*a**S(2)*d*(m + n + S(1)) + B*b*(a*d*(n + S(1)) + b*c*(m + S(-1))) + b*(A*b*d*(m + n + S(1)) - B*(-a*d*(S(2)*m + n) + b*c*m))*cos(e + f*x)**S(2) + (-B*b*(a*c - b*d*(m + n)) + a*d*(S(2)*A*b + B*a)*(m + n + S(1)))*cos(e + f*x), x), x), x) + Simp(B*b*(a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(1))), x)


def replacement2889(A, B, b, c, d, e, f, x):
    return Dist(B*d/b**S(2), Int(sqrt(b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x) + Int((A*c + (A*d + B*c)*sin(e + f*x))/((b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x)


def replacement2890(A, B, b, c, d, e, f, x):
    return Dist(B*d/b**S(2), Int(sqrt(b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x) + Int((A*c + (A*d + B*c)*cos(e + f*x))/((b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x)


def replacement2891(A, B, a, b, c, d, e, f, x):
    return Dist(B/b, Int(sqrt(c + d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x) + Dist((A*b - B*a)/b, Int(sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**(S(3)/2), x), x)


def replacement2892(A, B, a, b, c, d, e, f, x):
    return Dist(B/b, Int(sqrt(c + d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x) + Dist((A*b - B*a)/b, Int(sqrt(c + d*cos(e + f*x))/(a + b*cos(e + f*x))**(S(3)/2), x), x)


def replacement2893(A, B, a, b, d, e, f, x):
    return Dist(d/(a**S(2) - b**S(2)), Int((A*b - B*a + (A*a - B*b)*sin(e + f*x))/((d*sin(e + f*x))**(S(3)/2)*sqrt(a + b*sin(e + f*x))), x), x) + Simp(S(2)*(A*b - B*a)*cos(e + f*x)/(f*sqrt(d*sin(e + f*x))*sqrt(a + b*sin(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2894(A, B, a, b, d, e, f, x):
    return Dist(d/(a**S(2) - b**S(2)), Int((A*b - B*a + (A*a - B*b)*cos(e + f*x))/((d*cos(e + f*x))**(S(3)/2)*sqrt(a + b*cos(e + f*x))), x), x) + Simp(-S(2)*(A*b - B*a)*sin(e + f*x)/(f*sqrt(d*cos(e + f*x))*sqrt(a + b*cos(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement2895(A, B, b, c, d, e, f, x):
    return Simp(-S(2)*A*sqrt(c*(S(1) - S(1)/sin(e + f*x))/(c + d))*sqrt(c*(S(1) + S(1)/sin(e + f*x))/(c - d))*(c - d)*EllipticE(asin(sqrt(c + d*sin(e + f*x))/(sqrt(b*sin(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt((c + d)/b, S(2))*tan(e + f*x)/(b*c**S(2)*f), x)


def replacement2896(A, B, b, c, d, e, f, x):
    return Simp(S(2)*A*sqrt(c*(S(1) - S(1)/cos(e + f*x))/(c + d))*sqrt(c*(S(1) + S(1)/cos(e + f*x))/(c - d))*(c - d)*EllipticE(asin(sqrt(c + d*cos(e + f*x))/(sqrt(b*cos(e + f*x))*Rt((c + d)/b, S(2)))), -(c + d)/(c - d))*Rt((c + d)/b, S(2))/(b*c**S(2)*f*tan(e + f*x)), x)


def replacement2897(A, B, b, c, d, e, f, x):
    return -Dist(sqrt(-b*sin(e + f*x))/sqrt(b*sin(e + f*x)), Int((A + B*sin(e + f*x))/((-b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2898(A, B, b, c, d, e, f, x):
    return -Dist(sqrt(-b*cos(e + f*x))/sqrt(b*cos(e + f*x)), Int((A + B*cos(e + f*x))/((-b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2899(A, B, a, b, c, d, e, f, x):
    return Simp(-S(2)*A*sqrt((-a*d + b*c)*(sin(e + f*x) + S(1))/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(S(1) - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*(c - d)*EllipticE(asin(sqrt(c + d*sin(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*sin(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(f*(-a*d + b*c)**S(2)*Rt((a + b)/(c + d), S(2))*cos(e + f*x)), x)


def replacement2900(A, B, a, b, c, d, e, f, x):
    return Simp(S(2)*A*sqrt((-a*d + b*c)*(cos(e + f*x) + S(1))/((a + b*cos(e + f*x))*(c - d)))*sqrt(-(S(1) - cos(e + f*x))*(-a*d + b*c)/((a + b*cos(e + f*x))*(c + d)))*(a + b*cos(e + f*x))*(c - d)*EllipticE(asin(sqrt(c + d*cos(e + f*x))*Rt((a + b)/(c + d), S(2))/sqrt(a + b*cos(e + f*x))), (a - b)*(c + d)/((a + b)*(c - d)))/(f*(-a*d + b*c)**S(2)*Rt((a + b)/(c + d), S(2))*sin(e + f*x)), x)


def replacement2901(A, B, a, b, c, d, e, f, x):
    return Dist(sqrt(-c - d*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), Int((A + B*sin(e + f*x))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(-c - d*sin(e + f*x))), x), x)


def replacement2902(A, B, a, b, c, d, e, f, x):
    return Dist(sqrt(-c - d*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), Int((A + B*cos(e + f*x))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(-c - d*cos(e + f*x))), x), x)


def replacement2903(A, B, a, b, c, d, e, f, x):
    return Dist((A - B)/(a - b), Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x) - Dist((A*b - B*a)/(a - b), Int((sin(e + f*x) + S(1))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2904(A, B, a, b, c, d, e, f, x):
    return Dist((A - B)/(a - b), Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x) - Dist((A*b - B*a)/(a - b), Int((cos(e + f*x) + S(1))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2905(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(-1))*Simp(c*(m + S(1))*(A*a - B*b) + d*n*(A*b - B*a) - d*(A*b - B*a)*(m + n + S(2))*sin(e + f*x)**S(2) + (-c*(m + S(2))*(A*b - B*a) + d*(m + S(1))*(A*a - B*b))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*(-A*b + B*a)*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2906(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(-1))*Simp(c*(m + S(1))*(A*a - B*b) + d*n*(A*b - B*a) - d*(A*b - B*a)*(m + n + S(2))*cos(e + f*x)**S(2) + (-c*(m + S(2))*(A*b - B*a) + d*(m + S(1))*(A*a - B*b))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*(-A*b + B*a)*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2907(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(b*d*(A*b - B*a)*(m + n + S(2)) - b*d*(A*b - B*a)*(m + n + S(3))*sin(e + f*x)**S(2) + (m + S(1))*(A*a - B*b)*(-a*d + b*c) + (A*b - B*a)*(a*d*(m + S(1)) - b*c*(m + S(2)))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(1))*(A*b**S(2) - B*a*b)*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement2908(A, B, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(b*d*(A*b - B*a)*(m + n + S(2)) - b*d*(A*b - B*a)*(m + n + S(3))*cos(e + f*x)**S(2) + (m + S(1))*(A*a - B*b)*(-a*d + b*c) + (A*b - B*a)*(a*d*(m + S(1)) - b*c*(m + S(2)))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(1))*(A*b**S(2) - B*a*b)*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement2909(A, B, a, b, c, d, e, f, x):
    return Dist((A*b - B*a)/(-a*d + b*c), Int(S(1)/(a + b*sin(e + f*x)), x), x) + Dist((-A*d + B*c)/(-a*d + b*c), Int(S(1)/(c + d*sin(e + f*x)), x), x)


def replacement2910(A, B, a, b, c, d, e, f, x):
    return Dist((A*b - B*a)/(-a*d + b*c), Int(S(1)/(a + b*cos(e + f*x)), x), x) + Dist((-A*d + B*c)/(-a*d + b*c), Int(S(1)/(c + d*cos(e + f*x)), x), x)


def replacement2911(A, B, a, b, c, d, e, f, m, x):
    return Dist(B/d, Int((a + b*sin(e + f*x))**m, x), x) - Dist((-A*d + B*c)/d, Int((a + b*sin(e + f*x))**m/(c + d*sin(e + f*x)), x), x)


def replacement2912(A, B, a, b, c, d, e, f, m, x):
    return Dist(B/d, Int((a + b*cos(e + f*x))**m, x), x) - Dist((-A*d + B*c)/d, Int((a + b*cos(e + f*x))**m/(c + d*cos(e + f*x)), x), x)


def replacement2913(A, B, a, b, c, d, e, f, n, x):
    return Dist(S(1)/(S(2)*n + S(3)), Int((c + d*sin(e + f*x))**(n + S(-1))*Simp(A*a*c*(S(2)*n + S(3)) + B*(S(2)*a*d*n + b*c) + (A*(S(2)*n + S(3))*(a*d + b*c) + B*(S(2)*n + S(1))*(a*c + b*d))*sin(e + f*x) + (A*b*d*(S(2)*n + S(3)) + B*(a*d + S(2)*b*c*n))*sin(e + f*x)**S(2), x)/sqrt(a + b*sin(e + f*x)), x), x) + Simp(-S(2)*B*sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))**n*cos(e + f*x)/(f*(S(2)*n + S(3))), x)


def replacement2914(A, B, a, b, c, d, e, f, n, x):
    return Dist(S(1)/(S(2)*n + S(3)), Int((c + d*cos(e + f*x))**(n + S(-1))*Simp(A*a*c*(S(2)*n + S(3)) + B*(S(2)*a*d*n + b*c) + (A*(S(2)*n + S(3))*(a*d + b*c) + B*(S(2)*n + S(1))*(a*c + b*d))*cos(e + f*x) + (A*b*d*(S(2)*n + S(3)) + B*(a*d + S(2)*b*c*n))*cos(e + f*x)**S(2), x)/sqrt(a + b*cos(e + f*x)), x), x) + Simp(S(2)*B*sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))**n*sin(e + f*x)/(f*(S(2)*n + S(3))), x)


def replacement2915(A, B, a, b, e, f, x):
    return Simp(S(4)*A*EllipticPi(S(-1), -asin(cos(e + f*x)/(sin(e + f*x) + S(1))), -(a - b)/(a + b))/(f*sqrt(a + b)), x)


def replacement2916(A, B, a, b, e, f, x):
    return Simp(S(4)*A*EllipticPi(S(-1), asin(sin(e + f*x)/(cos(e + f*x) + S(1))), -(a - b)/(a + b))/(f*sqrt(a + b)), x)


def replacement2917(A, B, a, b, d, e, f, x):
    return Dist(sqrt(sin(e + f*x))/sqrt(d*sin(e + f*x)), Int((A + B*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(sin(e + f*x))), x), x)


def replacement2918(A, B, a, b, d, e, f, x):
    return Dist(sqrt(cos(e + f*x))/sqrt(d*cos(e + f*x)), Int((A + B*cos(e + f*x))/(sqrt(a + b*cos(e + f*x))*sqrt(cos(e + f*x))), x), x)


def replacement2919(A, B, a, b, c, d, e, f, x):
    return Dist(B/d, Int(sqrt(c + d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x) - Dist((-A*d + B*c)/d, Int(S(1)/(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))), x), x)


def replacement2920(A, B, a, b, c, d, e, f, x):
    return Dist(B/d, Int(sqrt(c + d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x) - Dist((-A*d + B*c)/d, Int(S(1)/(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))), x), x)


def replacement2921(A, B, a, b, c, d, e, f, m, n, x):
    return Int((A + B*sin(e + f*x))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement2922(A, B, a, b, c, d, e, f, m, n, x):
    return Int((A + B*cos(e + f*x))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement2923(A, B, a, b, c, d, e, f, m, n, p, x):
    return Dist(sqrt(a + b*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/(f*cos(e + f*x)), Subst(Int((A + B*x)**p*(a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, sin(e + f*x)), x)


def replacement2924(A, B, a, b, c, d, e, f, m, n, p, x):
    return -Dist(sqrt(a + b*cos(e + f*x))*sqrt(c + d*cos(e + f*x))/(f*sin(e + f*x)), Subst(Int((A + B*x)**p*(a + b*x)**(m + S(-1)/2)*(c + d*x)**(n + S(-1)/2), x), x, cos(e + f*x)), x)


def replacement2925(B, C, b, e, f, m, x):
    return Dist(S(1)/b, Int((b*sin(e + f*x))**(m + S(1))*(B + C*sin(e + f*x)), x), x)


def replacement2926(B, C, b, e, f, m, x):
    return Dist(S(1)/b, Int((b*cos(e + f*x))**(m + S(1))*(B + C*cos(e + f*x)), x), x)


def replacement2927(A, C, e, f, m, x):
    return -Dist(S(1)/f, Subst(Int((S(1) - x**S(2))**(m/S(2) + S(-1)/2)*(A - C*x**S(2) + C), x), x, cos(e + f*x)), x)


def replacement2928(A, C, e, f, m, x):
    return Dist(S(1)/f, Subst(Int((S(1) - x**S(2))**(m/S(2) + S(-1)/2)*(A - C*x**S(2) + C), x), x, sin(e + f*x)), x)


def replacement2929(A, C, b, e, f, m, x):
    return Simp(A*(b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(1))), x)


def replacement2930(A, C, b, e, f, m, x):
    return -Simp(A*(b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(1))), x)


def replacement2931(A, C, b, e, f, m, x):
    return Dist((A*(m + S(2)) + C*(m + S(1)))/(b**S(2)*(m + S(1))), Int((b*sin(e + f*x))**(m + S(2)), x), x) + Simp(A*(b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(1))), x)


def replacement2932(A, C, b, e, f, m, x):
    return Dist((A*(m + S(2)) + C*(m + S(1)))/(b**S(2)*(m + S(1))), Int((b*cos(e + f*x))**(m + S(2)), x), x) - Simp(A*(b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(1))), x)


def replacement2933(A, C, b, e, f, m, x):
    return Dist((A*(m + S(2)) + C*(m + S(1)))/(m + S(2)), Int((b*sin(e + f*x))**m, x), x) - Simp(C*(b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2934(A, C, b, e, f, m, x):
    return Dist((A*(m + S(2)) + C*(m + S(1)))/(m + S(2)), Int((b*cos(e + f*x))**m, x), x) + Simp(C*(b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2935(A, B, C, a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b*sin(e + f*x), x), x), x)


def replacement2936(A, B, C, a, b, e, f, m, x):
    return Dist(b**(S(-2)), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(B*b - C*a + C*b*cos(e + f*x), x), x), x)


def replacement2937(A, C, a, b, e, f, m, x):
    return Dist(C/b**S(2), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(-a + b*sin(e + f*x), x), x), x)


def replacement2938(A, C, a, b, e, f, m, x):
    return Dist(C/b**S(2), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(-a + b*cos(e + f*x), x), x), x)


def replacement2939(A, B, C, a, b, e, f, m, x):
    return Dist(C, Int((a + b*sin(e + f*x))**m*(sin(e + f*x) + S(1))**S(2), x), x) + Dist(A - C, Int((a + b*sin(e + f*x))**m*(sin(e + f*x) + S(1)), x), x)


def replacement2940(A, B, C, a, b, e, f, m, x):
    return Dist(C, Int((a + b*cos(e + f*x))**m*(cos(e + f*x) + S(1))**S(2), x), x) + Dist(A - C, Int((a + b*cos(e + f*x))**m*(cos(e + f*x) + S(1)), x), x)


def replacement2941(A, C, a, b, e, f, m, x):
    return Dist(C, Int((a + b*sin(e + f*x))**m*(sin(e + f*x) + S(1))**S(2), x), x) + Dist(A - C, Int((a + b*sin(e + f*x))**m*(sin(e + f*x) + S(1)), x), x)


def replacement2942(A, C, a, b, e, f, m, x):
    return Dist(C, Int((a + b*cos(e + f*x))**m*(cos(e + f*x) + S(1))**S(2), x), x) + Dist(A - C, Int((a + b*cos(e + f*x))**m*(cos(e + f*x) + S(1)), x), x)


def replacement2943(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(A*a*(m + S(1)) + C*b*(S(2)*m + S(1))*sin(e + f*x) + m*(B*b - C*a), x), x), x) + Simp((a + b*sin(e + f*x))**m*(A*b - B*a + C*b)*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2944(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(A*a*(m + S(1)) + C*b*(S(2)*m + S(1))*cos(e + f*x) + m*(B*b - C*a), x), x), x) - Simp((a + b*cos(e + f*x))**m*(A*b - B*a + C*b)*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2945(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(A*a*(m + S(1)) - C*a*m + C*b*(S(2)*m + S(1))*sin(e + f*x), x), x), x) + Simp(b*(A + C)*(a + b*sin(e + f*x))**m*cos(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2946(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(a**S(2)*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(A*a*(m + S(1)) - C*a*m + C*b*(S(2)*m + S(1))*cos(e + f*x), x), x), x) - Simp(b*(A + C)*(a + b*cos(e + f*x))**m*sin(e + f*x)/(a*f*(S(2)*m + S(1))), x)


def replacement2947(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*a - B*b + C*a) - (A*b**S(2) - B*a*b + C*a**S(2) + b*(m + S(1))*(A*b - B*a + C*b))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*cos(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2948(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(b*(m + S(1))*(A*a - B*b + C*a) - (A*b**S(2) - B*a*b + C*a**S(2) + b*(m + S(1))*(A*b - B*a + C*b))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*sin(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2949(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(a*b*(A + C)*(m + S(1)) - (A*b**S(2) + C*a**S(2) + b**S(2)*(A + C)*(m + S(1)))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*cos(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2950(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(a*b*(A + C)*(m + S(1)) - (A*b**S(2) + C*a**S(2) + b**S(2)*(A + C)*(m + S(1)))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*sin(e + f*x)/(b*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2951(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*sin(e + f*x))**m*Simp(A*b*(m + S(2)) + C*b*(m + S(1)) + (B*b*(m + S(2)) - C*a)*sin(e + f*x), x), x), x) - Simp(C*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2952(A, B, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*cos(e + f*x))**m*Simp(A*b*(m + S(2)) + C*b*(m + S(1)) + (B*b*(m + S(2)) - C*a)*cos(e + f*x), x), x), x) + Simp(C*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2953(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*sin(e + f*x))**m*Simp(A*b*(m + S(2)) - C*a*sin(e + f*x) + C*b*(m + S(1)), x), x), x) - Simp(C*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*(m + S(2))), x)


def replacement2954(A, C, a, b, e, f, m, x):
    return Dist(S(1)/(b*(m + S(2))), Int((a + b*cos(e + f*x))**m*Simp(A*b*(m + S(2)) - C*a*cos(e + f*x) + C*b*(m + S(1)), x), x), x) + Simp(C*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*(m + S(2))), x)


def replacement2955(A, B, C, b, e, f, m, p, x):
    return Dist((b*sin(e + f*x))**(-m*p)*(b*sin(e + f*x)**p)**m, Int((b*sin(e + f*x))**(m*p)*(A + B*sin(e + f*x) + C*sin(e + f*x)**S(2)), x), x)


def replacement2956(A, B, C, b, e, f, m, p, x):
    return Dist((b*cos(e + f*x))**(-m*p)*(b*cos(e + f*x)**p)**m, Int((b*cos(e + f*x))**(m*p)*(A + B*cos(e + f*x) + C*cos(e + f*x)**S(2)), x), x)


def replacement2957(A, C, b, e, f, m, p, x):
    return Dist((b*sin(e + f*x))**(-m*p)*(b*sin(e + f*x)**p)**m, Int((b*sin(e + f*x))**(m*p)*(A + C*sin(e + f*x)**S(2)), x), x)


def replacement2958(A, C, b, e, f, m, p, x):
    return Dist((b*cos(e + f*x))**(-m*p)*(b*cos(e + f*x)**p)**m, Int((b*cos(e + f*x))**(m*p)*(A + C*cos(e + f*x)**S(2)), x), x)


def replacement2959(A, B, C, a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(-C*b*d*(a**S(2) - b**S(2))*(m + S(1))*sin(e + f*x)**S(2) + b*(m + S(1))*(-A*b*(a*c - b*d) + (B*b - C*a)*(-a*d + b*c)) + (B*b*(a**S(2)*d - a*b*c*(m + S(2)) + b**S(2)*d*(m + S(1))) + (-a*d + b*c)*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1)))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(-a*d + b*c)*(A*b**S(2) - B*a*b + C*a**S(2))*cos(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2960(A, B, C, a, b, c, d, e, f, m, x):
    return -Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(-C*b*d*(a**S(2) - b**S(2))*(m + S(1))*cos(e + f*x)**S(2) + b*(m + S(1))*(-A*b*(a*c - b*d) + (B*b - C*a)*(-a*d + b*c)) + (B*b*(a**S(2)*d - a*b*c*(m + S(2)) + b**S(2)*d*(m + S(1))) + (-a*d + b*c)*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1)))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(-a*d + b*c)*(A*b**S(2) - B*a*b + C*a**S(2))*sin(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2961(A, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*Simp(C*b*d*(a**S(2) - b**S(2))*(m + S(1))*sin(e + f*x)**S(2) + b*(m + S(1))*(A*b*(a*c - b*d) + C*a*(-a*d + b*c)) - (-a*d + b*c)*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*(-a*d + b*c)*cos(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2962(A, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b**S(2)*(a**S(2) - b**S(2))*(m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*Simp(C*b*d*(a**S(2) - b**S(2))*(m + S(1))*cos(e + f*x)**S(2) + b*(m + S(1))*(A*b*(a*c - b*d) + C*a*(-a*d + b*c)) - (-a*d + b*c)*(A*b**S(2)*(m + S(2)) + C*(a**S(2) + b**S(2)*(m + S(1))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(A*b**S(2) + C*a**S(2))*(-a*d + b*c)*sin(e + f*x)/(b**S(2)*f*(a**S(2) - b**S(2))*(m + S(1))), x)


def replacement2963(A, B, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b*sin(e + f*x))**m*Simp(A*b*c*(m + S(3)) + C*a*d + b*(B*c*(m + S(3)) + d*(A*(m + S(3)) + C*(m + S(2))))*sin(e + f*x) - (S(2)*C*a*d - b*(m + S(3))*(B*d + C*c))*sin(e + f*x)**S(2), x), x), x) - Simp(C*d*(a + b*sin(e + f*x))**(m + S(1))*sin(e + f*x)*cos(e + f*x)/(b*f*(m + S(3))), x)


def replacement2964(A, B, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b*cos(e + f*x))**m*Simp(A*b*c*(m + S(3)) + C*a*d + b*(B*c*(m + S(3)) + d*(A*(m + S(3)) + C*(m + S(2))))*cos(e + f*x) - (S(2)*C*a*d - b*(m + S(3))*(B*d + C*c))*cos(e + f*x)**S(2), x), x), x) + Simp(C*d*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)*cos(e + f*x)/(b*f*(m + S(3))), x)


def replacement2965(A, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b*sin(e + f*x))**m*Simp(A*b*c*(m + S(3)) + C*a*d + b*d*(A*(m + S(3)) + C*(m + S(2)))*sin(e + f*x) - (S(2)*C*a*d - C*b*c*(m + S(3)))*sin(e + f*x)**S(2), x), x), x) - Simp(C*d*(a + b*sin(e + f*x))**(m + S(1))*sin(e + f*x)*cos(e + f*x)/(b*f*(m + S(3))), x)


def replacement2966(A, C, a, b, c, d, e, f, m, x):
    return Dist(S(1)/(b*(m + S(3))), Int((a + b*cos(e + f*x))**m*Simp(A*b*c*(m + S(3)) + C*a*d + b*d*(A*(m + S(3)) + C*(m + S(2)))*cos(e + f*x) - (S(2)*C*a*d - C*b*c*(m + S(3)))*cos(e + f*x)**S(2), x), x), x) + Simp(C*d*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)*cos(e + f*x)/(b*f*(m + S(3))), x)


def replacement2967(A, B, C, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(S(2)*b*c*d*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(A*(c**S(2)*(m + S(1)) + d**S(2)*(S(2)*m + n + S(2))) - B*c*d*(m - n + S(-1)) - C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-C*c*(S(3)*m - n) + (A*c + B*d)*(m + n + S(2)))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*a - B*b + C*a)*cos(e + f*x)/(S(2)*b*c*f*(S(2)*m + S(1))), x)


def replacement2968(A, B, C, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(S(2)*b*c*d*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(A*(c**S(2)*(m + S(1)) + d**S(2)*(S(2)*m + n + S(2))) - B*c*d*(m - n + S(-1)) - C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(-C*c*(S(3)*m - n) + (A*c + B*d)*(m + n + S(2)))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*a - B*b + C*a)*sin(e + f*x)/(S(2)*b*c*f*(S(2)*m + S(1))), x)


def replacement2969(A, C, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(S(2)*b*c*d*(S(2)*m + S(1))), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(A*(c**S(2)*(m + S(1)) + d**S(2)*(S(2)*m + n + S(2))) - C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(A*c*(m + n + S(2)) - C*c*(S(3)*m - n))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*a + C*a)*cos(e + f*x)/(S(2)*b*c*f*(S(2)*m + S(1))), x)


def replacement2970(A, C, a, b, c, d, e, f, m, n, x):
    return -Dist(S(1)/(S(2)*b*c*d*(S(2)*m + S(1))), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(A*(c**S(2)*(m + S(1)) + d**S(2)*(S(2)*m + n + S(2))) - C*(c**S(2)*m - d**S(2)*(n + S(1))) + d*(A*c*(m + n + S(2)) - C*c*(S(3)*m - n))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*a + C*a)*sin(e + f*x)/(S(2)*b*c*f*(S(2)*m + S(1))), x)


def replacement2971(A, B, C, a, b, c, d, e, f, m, x):
    return Int((a + b*sin(e + f*x))**m*Simp(A + B*sin(e + f*x) + C, x)/sqrt(c + d*sin(e + f*x)), x) + Simp(-S(2)*C*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*sqrt(c + d*sin(e + f*x))*(S(2)*m + S(3))), x)


def replacement2972(A, B, C, a, b, c, d, e, f, m, x):
    return Int((a + b*cos(e + f*x))**m*Simp(A + B*cos(e + f*x) + C, x)/sqrt(c + d*cos(e + f*x)), x) + Simp(S(2)*C*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*sqrt(c + d*cos(e + f*x))*(S(2)*m + S(3))), x)


def replacement2973(A, C, a, b, c, d, e, f, m, x):
    return Dist(A + C, Int((a + b*sin(e + f*x))**m/sqrt(c + d*sin(e + f*x)), x), x) + Simp(-S(2)*C*(a + b*sin(e + f*x))**(m + S(1))*cos(e + f*x)/(b*f*sqrt(c + d*sin(e + f*x))*(S(2)*m + S(3))), x)


def replacement2974(A, C, a, b, c, d, e, f, m, x):
    return Dist(A + C, Int((a + b*cos(e + f*x))**m/sqrt(c + d*cos(e + f*x)), x), x) + Simp(S(2)*C*(a + b*cos(e + f*x))**(m + S(1))*sin(e + f*x)/(b*f*sqrt(c + d*cos(e + f*x))*(S(2)*m + S(3))), x)


def replacement2975(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + (B*b*d*(m + n + S(2)) - C*b*c*(S(2)*m + S(1)))*sin(e + f*x), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2976(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + (B*b*d*(m + n + S(2)) - C*b*c*(S(2)*m + S(1)))*cos(e + f*x), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2977(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) - C*b*c*(S(2)*m + S(1))*sin(e + f*x) + C*(a*c*m + b*d*(n + S(1))), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2978(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) - C*b*c*(S(2)*m + S(1))*cos(e + f*x) + C*(a*c*m + b*d*(n + S(1))), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2979(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(A*(a*c*(m + S(1)) - b*d*(S(2)*m + n + S(2))) + B*(a*d*(n + S(1)) + b*c*m) - C*(a*c*m + b*d*(n + S(1))) + (C*(-a*d*(m - n + S(-1)) + b*c*(S(2)*m + S(1))) + d*(A*a - B*b)*(m + n + S(2)))*sin(e + f*x), x), x), x) + Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*a - B*b + C*a)*cos(e + f*x)/(f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2980(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(A*(a*c*(m + S(1)) - b*d*(S(2)*m + n + S(2))) + B*(a*d*(n + S(1)) + b*c*m) - C*(a*c*m + b*d*(n + S(1))) + (C*(-a*d*(m - n + S(-1)) + b*c*(S(2)*m + S(1))) + d*(A*a - B*b)*(m + n + S(2)))*cos(e + f*x), x), x), x) - Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*a - B*b + C*a)*sin(e + f*x)/(f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2981(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(A*(a*c*(m + S(1)) - b*d*(S(2)*m + n + S(2))) - C*(a*c*m + b*d*(n + S(1))) + (A*a*d*(m + n + S(2)) + C*(-a*d*(m - n + S(-1)) + b*c*(S(2)*m + S(1))))*sin(e + f*x), x), x), x) + Simp(a*(A + C)*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2982(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*(S(2)*m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(A*(a*c*(m + S(1)) - b*d*(S(2)*m + n + S(2))) - C*(a*c*m + b*d*(n + S(1))) + (A*a*d*(m + n + S(2)) + C*(-a*d*(m - n + S(-1)) + b*c*(S(2)*m + S(1))))*cos(e + f*x), x), x), x) - Simp(a*(A + C)*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(f*(S(2)*m + S(1))*(-a*d + b*c)), x)


def replacement2983(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*d*(a*d*m + b*c*(n + S(1))) + b*(-C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(2)))*sin(e + f*x) + (-B*d + C*c)*(a*c*m + b*d*(n + S(1))), x), x), x) - Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2984(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*d*(a*d*m + b*c*(n + S(1))) + b*(-C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(2)))*cos(e + f*x) + (-B*d + C*c)*(a*c*m + b*d*(n + S(1))), x), x), x) + Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2985(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*d*(a*d*m + b*c*(n + S(1))) + C*c*(a*c*m + b*d*(n + S(1))) - b*(A*d**S(2)*(m + n + S(2)) + C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2986(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*d*(a*d*m + b*c*(n + S(1))) + C*c*(a*c*m + b*d*(n + S(1))) - b*(A*d**S(2)*(m + n + S(2)) + C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2987(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + (B*b*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1))))*sin(e + f*x), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2988(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + (B*b*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1))))*cos(e + f*x), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2989(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + C*(a*d*m - b*c*(m + S(1)))*sin(e + f*x), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2990(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(b*d*(m + n + S(2))), Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*Simp(A*b*d*(m + n + S(2)) + C*(a*c*m + b*d*(n + S(1))) + C*(a*d*m - b*c*(m + S(1)))*cos(e + f*x), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2991(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*d*(a*c*(n + S(1)) + b*d*m) + b*(-C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(2)))*sin(e + f*x)**S(2) + (-B*d + C*c)*(a*d*(n + S(1)) + b*c*m) - (-C*(-a*(c**S(2) + d**S(2)*(n + S(1))) + b*c*d*(n + S(1))) + d*(A*(a*d*(n + S(2)) - b*c*(n + S(1))) + B*(-a*c*(n + S(2)) + b*d*(n + S(1)))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2992(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*d*(a*c*(n + S(1)) + b*d*m) + b*(-C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))) + d*(-A*d + B*c)*(m + n + S(2)))*cos(e + f*x)**S(2) + (-B*d + C*c)*(a*d*(n + S(1)) + b*c*m) - (-C*(-a*(c**S(2) + d**S(2)*(n + S(1))) + b*c*d*(n + S(1))) + d*(A*(a*d*(n + S(2)) - b*c*(n + S(1))) + B*(-a*c*(n + S(2)) + b*d*(n + S(1)))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*d**S(2) - B*c*d + C*c**S(2))*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2993(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**(n + S(1))*Simp(A*d*(a*c*(n + S(1)) + b*d*m) + C*c*(a*d*(n + S(1)) + b*c*m) - b*(A*d**S(2)*(m + n + S(2)) + C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))))*sin(e + f*x)**S(2) - (A*d*(a*d*(n + S(2)) - b*c*(n + S(1))) - C*(-a*(c**S(2) + d**S(2)*(n + S(1))) + b*c*d*(n + S(1))))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*cos(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2994(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(c**S(2) - d**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**(n + S(1))*Simp(A*d*(a*c*(n + S(1)) + b*d*m) + C*c*(a*d*(n + S(1)) + b*c*m) - b*(A*d**S(2)*(m + n + S(2)) + C*(c**S(2)*(m + S(1)) + d**S(2)*(n + S(1))))*cos(e + f*x)**S(2) - (A*d*(a*d*(n + S(2)) - b*c*(n + S(1))) - C*(-a*(c**S(2) + d**S(2)*(n + S(1))) + b*c*d*(n + S(1))))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*(A*d**S(2) + C*c**S(2))*sin(e + f*x)/(d*f*(c**S(2) - d**S(2))*(n + S(1))), x)


def replacement2995(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(2))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*Simp(A*a*d*(m + n + S(2)) + C*(a*d*(n + S(1)) + b*c*m) + (-C*(a*c - b*d*(m + n + S(1))) + d*(A*b + B*a)*(m + n + S(2)))*sin(e + f*x) + (B*b*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1))))*sin(e + f*x)**S(2), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2996(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(2))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*Simp(A*a*d*(m + n + S(2)) + C*(a*d*(n + S(1)) + b*c*m) + (-C*(a*c - b*d*(m + n + S(1))) + d*(A*b + B*a)*(m + n + S(2)))*cos(e + f*x) + (B*b*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1))))*cos(e + f*x)**S(2), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2997(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(2))), Int((a + b*sin(e + f*x))**(m + S(-1))*(c + d*sin(e + f*x))**n*Simp(A*a*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1)))*sin(e + f*x)**S(2) + C*(a*d*(n + S(1)) + b*c*m) + (A*b*d*(m + n + S(2)) - C*(a*c - b*d*(m + n + S(1))))*sin(e + f*x), x), x), x) - Simp(C*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**(n + S(1))*cos(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2998(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/(d*(m + n + S(2))), Int((a + b*cos(e + f*x))**(m + S(-1))*(c + d*cos(e + f*x))**n*Simp(A*a*d*(m + n + S(2)) + C*(a*d*m - b*c*(m + S(1)))*cos(e + f*x)**S(2) + C*(a*d*(n + S(1)) + b*c*m) + (A*b*d*(m + n + S(2)) - C*(a*c - b*d*(m + n + S(1))))*cos(e + f*x), x), x), x) + Simp(C*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**(n + S(1))*sin(e + f*x)/(d*f*(m + n + S(2))), x)


def replacement2999(A, B, C, a, b, d, e, f, x):
    return Dist(S(1)/b, Int((A*b + (B*b - C*a)*sin(e + f*x))/(sqrt(d*sin(e + f*x))*(a + b*sin(e + f*x))**(S(3)/2)), x), x) + Dist(C/(b*d), Int(sqrt(d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x)


def replacement3000(A, B, C, a, b, d, e, f, x):
    return Dist(S(1)/b, Int((A*b + (B*b - C*a)*cos(e + f*x))/(sqrt(d*cos(e + f*x))*(a + b*cos(e + f*x))**(S(3)/2)), x), x) + Dist(C/(b*d), Int(sqrt(d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x)


def replacement3001(A, C, a, b, d, e, f, x):
    return Dist(S(1)/b, Int((A*b - C*a*sin(e + f*x))/(sqrt(d*sin(e + f*x))*(a + b*sin(e + f*x))**(S(3)/2)), x), x) + Dist(C/(b*d), Int(sqrt(d*sin(e + f*x))/sqrt(a + b*sin(e + f*x)), x), x)


def replacement3002(A, C, a, b, d, e, f, x):
    return Dist(S(1)/b, Int((A*b - C*a*cos(e + f*x))/(sqrt(d*cos(e + f*x))*(a + b*cos(e + f*x))**(S(3)/2)), x), x) + Dist(C/(b*d), Int(sqrt(d*cos(e + f*x))/sqrt(a + b*cos(e + f*x)), x), x)


def replacement3003(A, B, C, a, b, c, d, e, f, x):
    return Dist(b**(S(-2)), Int((A*b**S(2) - C*a**S(2) + b*(B*b - S(2)*C*a)*sin(e + f*x))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x) + Dist(C/b**S(2), Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x)


def replacement3004(A, B, C, a, b, c, d, e, f, x):
    return Dist(b**(S(-2)), Int((A*b**S(2) - C*a**S(2) + b*(B*b - S(2)*C*a)*cos(e + f*x))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x) + Dist(C/b**S(2), Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x)


def replacement3005(A, C, a, b, c, d, e, f, x):
    return Dist(b**(S(-2)), Int((A*b**S(2) - C*a**S(2) - S(2)*C*a*b*sin(e + f*x))/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x) + Dist(C/b**S(2), Int(sqrt(a + b*sin(e + f*x))/sqrt(c + d*sin(e + f*x)), x), x)


def replacement3006(A, C, a, b, c, d, e, f, x):
    return Dist(b**(S(-2)), Int((A*b**S(2) - C*a**S(2) - S(2)*C*a*b*cos(e + f*x))/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x) + Dist(C/b**S(2), Int(sqrt(a + b*cos(e + f*x))/sqrt(c + d*cos(e + f*x)), x), x)


def replacement3007(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(d*(m + n + S(2))*(A*b**S(2) - B*a*b + C*a**S(2)) - d*(m + n + S(3))*(A*b**S(2) - B*a*b + C*a**S(2))*sin(e + f*x)**S(2) + (m + S(1))*(-a*d + b*c)*(A*a - B*b + C*a) - (c*(A*b**S(2) - B*a*b + C*a**S(2)) + (m + S(1))*(-a*d + b*c)*(A*b - B*a + C*b))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement3008(A, B, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(d*(m + n + S(2))*(A*b**S(2) - B*a*b + C*a**S(2)) - d*(m + n + S(3))*(A*b**S(2) - B*a*b + C*a**S(2))*cos(e + f*x)**S(2) + (m + S(1))*(-a*d + b*c)*(A*a - B*b + C*a) - (c*(A*b**S(2) - B*a*b + C*a**S(2)) + (m + S(1))*(-a*d + b*c)*(A*b - B*a + C*b))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(1))*(A*b**S(2) - B*a*b + C*a**S(2))*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement3009(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**n*Simp(a*(A + C)*(m + S(1))*(-a*d + b*c) + d*(A*b**S(2) + C*a**S(2))*(m + n + S(2)) - d*(A*b**S(2) + C*a**S(2))*(m + n + S(3))*sin(e + f*x)**S(2) - (b*(A + C)*(m + S(1))*(-a*d + b*c) + c*(A*b**S(2) + C*a**S(2)))*sin(e + f*x), x), x), x) - Simp((a + b*sin(e + f*x))**(m + S(1))*(c + d*sin(e + f*x))**(n + S(1))*(A*b**S(2) + C*a**S(2))*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement3010(A, C, a, b, c, d, e, f, m, n, x):
    return Dist(S(1)/((a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), Int((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**n*Simp(a*(A + C)*(m + S(1))*(-a*d + b*c) + d*(A*b**S(2) + C*a**S(2))*(m + n + S(2)) - d*(A*b**S(2) + C*a**S(2))*(m + n + S(3))*cos(e + f*x)**S(2) - (b*(A + C)*(m + S(1))*(-a*d + b*c) + c*(A*b**S(2) + C*a**S(2)))*cos(e + f*x), x), x), x) + Simp((a + b*cos(e + f*x))**(m + S(1))*(c + d*cos(e + f*x))**(n + S(1))*(A*b**S(2) + C*a**S(2))*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(m + S(1))*(-a*d + b*c)), x)


def replacement3011(A, B, C, a, b, c, d, e, f, x):
    return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(b*(-a*d + b*c)), Int(S(1)/(a + b*sin(e + f*x)), x), x) - Dist((A*d**S(2) - B*c*d + C*c**S(2))/(d*(-a*d + b*c)), Int(S(1)/(c + d*sin(e + f*x)), x), x) + Simp(C*x/(b*d), x)


def replacement3012(A, B, C, a, b, c, d, e, f, x):
    return Dist((A*b**S(2) - B*a*b + C*a**S(2))/(b*(-a*d + b*c)), Int(S(1)/(a + b*cos(e + f*x)), x), x) - Dist((A*d**S(2) - B*c*d + C*c**S(2))/(d*(-a*d + b*c)), Int(S(1)/(c + d*cos(e + f*x)), x), x) + Simp(C*x/(b*d), x)


def replacement3013(A, C, a, b, c, d, e, f, x):
    return Dist((A*b**S(2) + C*a**S(2))/(b*(-a*d + b*c)), Int(S(1)/(a + b*sin(e + f*x)), x), x) - Dist((A*d**S(2) + C*c**S(2))/(d*(-a*d + b*c)), Int(S(1)/(c + d*sin(e + f*x)), x), x) + Simp(C*x/(b*d), x)


def replacement3014(A, C, a, b, c, d, e, f, x):
    return Dist((A*b**S(2) + C*a**S(2))/(b*(-a*d + b*c)), Int(S(1)/(a + b*cos(e + f*x)), x), x) - Dist((A*d**S(2) + C*c**S(2))/(d*(-a*d + b*c)), Int(S(1)/(c + d*cos(e + f*x)), x), x) + Simp(C*x/(b*d), x)


def replacement3015(A, B, C, a, b, c, d, e, f, x):
    return -Dist(S(1)/(b*d), Int(Simp(-A*b*d + C*a*c + (-B*b*d + C*a*d + C*b*c)*sin(e + f*x), x)/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x) + Dist(C/(b*d), Int(sqrt(a + b*sin(e + f*x)), x), x)


def replacement3016(A, B, C, a, b, c, d, e, f, x):
    return -Dist(S(1)/(b*d), Int(Simp(-A*b*d + C*a*c + (-B*b*d + C*a*d + C*b*c)*cos(e + f*x), x)/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x) + Dist(C/(b*d), Int(sqrt(a + b*cos(e + f*x)), x), x)


def replacement3017(A, C, a, b, c, d, e, f, x):
    return -Dist(S(1)/(b*d), Int(Simp(-A*b*d + C*a*c + (C*a*d + C*b*c)*sin(e + f*x), x)/(sqrt(a + b*sin(e + f*x))*(c + d*sin(e + f*x))), x), x) + Dist(C/(b*d), Int(sqrt(a + b*sin(e + f*x)), x), x)


def replacement3018(A, C, a, b, c, d, e, f, x):
    return -Dist(S(1)/(b*d), Int(Simp(-A*b*d + C*a*c + (C*a*d + C*b*c)*cos(e + f*x), x)/(sqrt(a + b*cos(e + f*x))*(c + d*cos(e + f*x))), x), x) + Dist(C/(b*d), Int(sqrt(a + b*cos(e + f*x)), x), x)


def replacement3019(A, B, C, a, b, c, d, e, f, x):
    return Dist(S(1)/(S(2)*d), Int(Simp(S(2)*A*a*d - C*(-a*d + b*c) + (S(2)*B*b*d - C*(a*d + b*c))*sin(e + f*x)**S(2) - S(2)*(C*a*c - d*(A*b + B*a))*sin(e + f*x), x)/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x) - Simp(C*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(d*f*sqrt(a + b*sin(e + f*x))), x)


def replacement3020(A, B, C, a, b, c, d, e, f, x):
    return Dist(S(1)/(S(2)*d), Int(Simp(S(2)*A*a*d - C*(-a*d + b*c) + (S(2)*B*b*d - C*(a*d + b*c))*cos(e + f*x)**S(2) - S(2)*(C*a*c - d*(A*b + B*a))*cos(e + f*x), x)/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x) + Simp(C*sqrt(c + d*cos(e + f*x))*sin(e + f*x)/(d*f*sqrt(a + b*cos(e + f*x))), x)


def replacement3021(A, C, a, b, c, d, e, f, x):
    return Dist(S(1)/(S(2)*d), Int(Simp(S(2)*A*a*d - C*(-a*d + b*c) - C*(a*d + b*c)*sin(e + f*x)**S(2) - S(2)*(-A*b*d + C*a*c)*sin(e + f*x), x)/((a + b*sin(e + f*x))**(S(3)/2)*sqrt(c + d*sin(e + f*x))), x), x) - Simp(C*sqrt(c + d*sin(e + f*x))*cos(e + f*x)/(d*f*sqrt(a + b*sin(e + f*x))), x)


def replacement3022(A, C, a, b, c, d, e, f, x):
    return Dist(S(1)/(S(2)*d), Int(Simp(S(2)*A*a*d - C*(-a*d + b*c) - C*(a*d + b*c)*cos(e + f*x)**S(2) - S(2)*(-A*b*d + C*a*c)*cos(e + f*x), x)/((a + b*cos(e + f*x))**(S(3)/2)*sqrt(c + d*cos(e + f*x))), x), x) + Simp(C*sqrt(c + d*cos(e + f*x))*sin(e + f*x)/(d*f*sqrt(a + b*cos(e + f*x))), x)


def replacement3023(A, B, C, a, b, c, d, e, f, m, n, x):
    return Int((a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n*(A + B*sin(e + f*x) + C*sin(e + f*x)**S(2)), x)


def replacement3024(A, B, C, a, b, c, d, e, f, m, n, x):
    return Int((a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n*(A + B*cos(e + f*x) + C*cos(e + f*x)**S(2)), x)


def replacement3025(A, C, a, b, c, d, e, f, m, n, x):
    return Int((A + C*sin(e + f*x)**S(2))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n, x)


def replacement3026(A, C, a, b, c, d, e, f, m, n, x):
    return Int((A + C*cos(e + f*x)**S(2))*(a + b*cos(e + f*x))**m*(c + d*cos(e + f*x))**n, x)


def replacement3027(A, B, C, b, c, d, e, f, m, n, p, x):
    return Dist((b*sin(e + f*x))**(-m*p)*(b*sin(e + f*x)**p)**m, Int((b*sin(e + f*x))**(m*p)*(c + d*sin(e + f*x))**n*(A + B*sin(e + f*x) + C*sin(e + f*x)**S(2)), x), x)


def replacement3028(A, B, C, b, c, d, e, f, m, n, p, x):
    return Dist((b*cos(e + f*x))**(-m*p)*(b*cos(e + f*x)**p)**m, Int((b*cos(e + f*x))**(m*p)*(c + d*cos(e + f*x))**n*(A + B*cos(e + f*x) + C*cos(e + f*x)**S(2)), x), x)


def replacement3029(A, C, b, c, d, e, f, m, n, p, x):
    return Dist((b*sin(e + f*x))**(-m*p)*(b*sin(e + f*x)**p)**m, Int((b*sin(e + f*x))**(m*p)*(A + C*sin(e + f*x)**S(2))*(c + d*sin(e + f*x))**n, x), x)


def replacement3030(A, C, b, c, d, e, f, m, n, p, x):
    return Dist((b*cos(e + f*x))**(-m*p)*(b*cos(e + f*x)**p)**m, Int((b*cos(e + f*x))**(m*p)*(A + C*cos(e + f*x)**S(2))*(c + d*cos(e + f*x))**n, x), x)


def replacement3031(a, b, c, d, n, x):
    return Simp(a*(a*cos(c + d*x) + b*sin(c + d*x))**n/(b*d*n), x)


def replacement3032(a, b, c, d, n, x):
    return -Dist(S(1)/d, Subst(Int((a**S(2) + b**S(2) - x**S(2))**(n/S(2) + S(-1)/2), x), x, -a*sin(c + d*x) + b*cos(c + d*x)), x)


def replacement3033(a, b, c, d, n, x):
    return Dist((a**S(2) + b**S(2))*(n + S(-1))/n, Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-2)), x), x) - Simp((-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))/(d*n), x)


def replacement3034(a, b, c, d, x):
    return -Dist(S(1)/d, Subst(Int(S(1)/(a**S(2) + b**S(2) - x**S(2)), x), x, -a*sin(c + d*x) + b*cos(c + d*x)), x)


def replacement3035(a, b, c, d, x):
    return Simp(sin(c + d*x)/(a*d*(a*cos(c + d*x) + b*sin(c + d*x))), x)


def replacement3036(a, b, c, d, n, x):
    return Dist((n + S(2))/((a**S(2) + b**S(2))*(n + S(1))), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(2)), x), x) + Simp((-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))/(d*(a**S(2) + b**S(2))*(n + S(1))), x)


def replacement3037(a, b, c, d, n, x):
    return Dist((a**S(2) + b**S(2))**(n/S(2)), Int(cos(c + d*x - ArcTan(a, b))**n, x), x)


def replacement3038(a, b, c, d, n, x):
    return Dist(((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**S(2) + b**S(2)))**(-n)*(a*cos(c + d*x) + b*sin(c + d*x))**n, Int(cos(c + d*x - ArcTan(a, b))**n, x), x)


def replacement3039(a, b, c, d, m, n, x):
    return Dist(S(2)*b, Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*sin(c + d*x)**(S(1) - n), x), x) - Simp(a*(a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*sin(c + d*x)**(S(1) - n)/(d*(n + S(-1))), x)


def replacement3040(a, b, c, d, m, n, x):
    return Dist(S(2)*a, Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)**(S(1) - n), x), x) + Simp(b*(a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)**(S(1) - n)/(d*(n + S(-1))), x)


def replacement3041(a, b, c, d, m, n, x):
    return Dist(S(1)/(S(2)*b), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))*sin(c + d*x)**(-n + S(-1)), x), x) + Simp(a*(a*cos(c + d*x) + b*sin(c + d*x))**n*sin(c + d*x)**(-n)/(S(2)*b*d*n), x)


def replacement3042(a, b, c, d, m, n, x):
    return Dist(S(1)/(S(2)*a), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))*cos(c + d*x)**(-n + S(-1)), x), x) - Simp(b*(a*cos(c + d*x) + b*sin(c + d*x))**n*cos(c + d*x)**(-n)/(S(2)*a*d*n), x)


def replacement3043(a, b, c, d, m, n, x):
    return Simp(a*(a*cos(c + d*x) + b*sin(c + d*x))**n*Hypergeometric2F1(S(1), n, n + S(1), (a/tan(c + d*x) + b)/(S(2)*b))*sin(c + d*x)**(-n)/(S(2)*b*d*n), x)


def replacement3044(a, b, c, d, m, n, x):
    return -Simp(b*(a*cos(c + d*x) + b*sin(c + d*x))**n*Hypergeometric2F1(S(1), n, n + S(1), (a + b*tan(c + d*x))/(S(2)*a))*cos(c + d*x)**(-n)/(S(2)*a*d*n), x)


def replacement3045(a, b, c, d, m, n, x):
    return Int((a/tan(c + d*x) + b)**n, x)


def replacement3046(a, b, c, d, m, n, x):
    return Int((a + b*tan(c + d*x))**n, x)


def replacement3047(a, b, c, d, m, n, x):
    return Dist(S(1)/d, Subst(Int(x**m*(a + b*x)**n*(x**S(2) + S(1))**(-m/S(2) - n/S(2) + S(-1)), x), x, tan(c + d*x)), x)


def replacement3048(a, b, c, d, m, n, x):
    return -Dist(S(1)/d, Subst(Int(x**m*(x**S(2) + S(1))**(-m/S(2) - n/S(2) + S(-1))*(a*x + b)**n, x), x, S(1)/tan(c + d*x)), x)


def replacement3049(a, b, c, d, m, n, x):
    return Int(ExpandTrig((a*cos(c + d*x) + b*sin(c + d*x))**n*sin(c + d*x)**m, x), x)


def replacement3050(a, b, c, d, m, n, x):
    return Int(ExpandTrig((a*cos(c + d*x) + b*sin(c + d*x))**n*cos(c + d*x)**m, x), x)


def replacement3051(a, b, c, d, m, n, x):
    return Dist(a**n*b**n, Int((a*sin(c + d*x) + b*cos(c + d*x))**(-n)*sin(c + d*x)**m, x), x)


def replacement3052(a, b, c, d, m, n, x):
    return Dist(a**n*b**n, Int((a*sin(c + d*x) + b*cos(c + d*x))**(-n)*cos(c + d*x)**m, x), x)


def replacement3053(a, b, c, d, n, x):
    return Dist(a**(S(-2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(2))/sin(c + d*x), x), x) - Dist(b/a**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1)), x), x) - Simp((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))/(a*d*(n + S(1))), x)


def replacement3054(a, b, c, d, n, x):
    return Dist(b**(S(-2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(2))/cos(c + d*x), x), x) - Dist(a/b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1)), x), x) + Simp((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))/(b*d*(n + S(1))), x)


def replacement3055(a, b, c, d, m, n, x):
    return Dist(a**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-2))*sin(c + d*x)**m, x), x) + Dist(S(2)*b, Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*sin(c + d*x)**(m + S(1)), x), x) - Dist(a**S(2) + b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-2))*sin(c + d*x)**(m + S(2)), x), x)


def replacement3056(a, b, c, d, m, n, x):
    return Dist(S(2)*a, Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-1))*cos(c + d*x)**(m + S(1)), x), x) + Dist(b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-2))*cos(c + d*x)**m, x), x) - Dist(a**S(2) + b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(-2))*cos(c + d*x)**(m + S(2)), x), x)


def replacement3057(a, b, c, d, x):
    return -Dist(a/(a**S(2) + b**S(2)), Int((-a*sin(c + d*x) + b*cos(c + d*x))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Simp(b*x/(a**S(2) + b**S(2)), x)


def replacement3058(a, b, c, d, x):
    return Dist(b/(a**S(2) + b**S(2)), Int((-a*sin(c + d*x) + b*cos(c + d*x))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Simp(a*x/(a**S(2) + b**S(2)), x)


def replacement3059(a, b, c, d, m, x):
    return Dist(a**S(2)/(a**S(2) + b**S(2)), Int(sin(c + d*x)**(m + S(-2))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Dist(b/(a**S(2) + b**S(2)), Int(sin(c + d*x)**(m + S(-1)), x), x) - Simp(a*sin(c + d*x)**(m + S(-1))/(d*(a**S(2) + b**S(2))*(m + S(-1))), x)


def replacement3060(a, b, c, d, m, x):
    return Dist(a/(a**S(2) + b**S(2)), Int(cos(c + d*x)**(m + S(-1)), x), x) + Dist(b**S(2)/(a**S(2) + b**S(2)), Int(cos(c + d*x)**(m + S(-2))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Simp(b*cos(c + d*x)**(m + S(-1))/(d*(a**S(2) + b**S(2))*(m + S(-1))), x)


def replacement3061(a, b, c, d, x):
    return -Dist(S(1)/a, Int((-a*sin(c + d*x) + b*cos(c + d*x))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Dist(S(1)/a, Int(S(1)/tan(c + d*x), x), x)


def replacement3062(a, b, c, d, x):
    return Dist(S(1)/b, Int((-a*sin(c + d*x) + b*cos(c + d*x))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Dist(S(1)/b, Int(tan(c + d*x), x), x)


def replacement3063(a, b, c, d, m, x):
    return -Dist(b/a**S(2), Int(sin(c + d*x)**(m + S(1)), x), x) + Dist((a**S(2) + b**S(2))/a**S(2), Int(sin(c + d*x)**(m + S(2))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) + Simp(sin(c + d*x)**(m + S(1))/(a*d*(m + S(1))), x)


def replacement3064(a, b, c, d, m, x):
    return -Dist(a/b**S(2), Int(cos(c + d*x)**(m + S(1)), x), x) + Dist((a**S(2) + b**S(2))/b**S(2), Int(cos(c + d*x)**(m + S(2))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x) - Simp(cos(c + d*x)**(m + S(1))/(b*d*(m + S(1))), x)


def replacement3065(a, b, c, d, m, n, x):
    return Dist(a**(S(-2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(2))*sin(c + d*x)**m, x), x) - Dist(S(2)*b/a**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))*sin(c + d*x)**(m + S(1)), x), x) + Dist((a**S(2) + b**S(2))/a**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**n*sin(c + d*x)**(m + S(2)), x), x)


def replacement3066(a, b, c, d, m, n, x):
    return Dist(b**(S(-2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(2))*cos(c + d*x)**m, x), x) - Dist(S(2)*a/b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**(n + S(1))*cos(c + d*x)**(m + S(1)), x), x) + Dist((a**S(2) + b**S(2))/b**S(2), Int((a*cos(c + d*x) + b*sin(c + d*x))**n*cos(c + d*x)**(m + S(2)), x), x)


def replacement3067(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((a*cos(c + d*x) + b*sin(c + d*x))**p*sin(c + d*x)**n*cos(c + d*x)**m, x), x)


def replacement3068(a, b, c, d, m, n, p, x):
    return Dist(a**p*b**p, Int((a*sin(c + d*x) + b*cos(c + d*x))**(-p)*sin(c + d*x)**n*cos(c + d*x)**m, x), x)


def replacement3069(a, b, c, d, m, n, x):
    return Dist(a/(a**S(2) + b**S(2)), Int(sin(c + d*x)**n*cos(c + d*x)**(m + S(-1)), x), x) + Dist(b/(a**S(2) + b**S(2)), Int(sin(c + d*x)**(n + S(-1))*cos(c + d*x)**m, x), x) - Dist(a*b/(a**S(2) + b**S(2)), Int(sin(c + d*x)**(n + S(-1))*cos(c + d*x)**(m + S(-1))/(a*cos(c + d*x) + b*sin(c + d*x)), x), x)


def replacement3070(a, b, c, d, m, n, x):
    return Int(ExpandTrig(sin(c + d*x)**n*cos(c + d*x)**m/(a*cos(c + d*x) + b*sin(c + d*x)), x), x)


def replacement3071(a, b, c, d, m, n, p, x):
    return Dist(a/(a**S(2) + b**S(2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(p + S(1))*sin(c + d*x)**n*cos(c + d*x)**(m + S(-1)), x), x) + Dist(b/(a**S(2) + b**S(2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**(p + S(1))*sin(c + d*x)**(n + S(-1))*cos(c + d*x)**m, x), x) - Dist(a*b/(a**S(2) + b**S(2)), Int((a*cos(c + d*x) + b*sin(c + d*x))**p*sin(c + d*x)**(n + S(-1))*cos(c + d*x)**(m + S(-1)), x), x)


def replacement3072(a, b, c, d, e, x):
    return Simp(-S(2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))), x)


def replacement3073(a, b, c, d, e, n, x):
    return Dist(a*(S(2)*n + S(-1))/n, Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1)), x), x) - Simp((-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1))/(e*n), x)


def replacement3074(a, b, c, d, e, x):
    return -Simp((-a*sin(d + e*x) + c)/(c*e*(-b*sin(d + e*x) + c*cos(d + e*x))), x)


def replacement3075(a, b, c, d, e, x):
    return Int(S(1)/sqrt(a + sqrt(b**S(2) + c**S(2))*cos(d + e*x - ArcTan(b, c))), x)


def replacement3076(a, b, c, d, e, n, x):
    return Dist((n + S(1))/(a*(S(2)*n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1)), x), x) + Simp((-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(S(2)*n + S(1))), x)


def replacement3077(a, b, c, d, e, x):
    return Dist(b/(c*e), Subst(Int(sqrt(a + x)/x, x), x, b*cos(d + e*x) + c*sin(d + e*x)), x)


def replacement3078(a, b, c, d, e, x):
    return Int(sqrt(a + sqrt(b**S(2) + c**S(2))*cos(d + e*x - ArcTan(b, c))), x)


def replacement3079(a, b, c, d, e, x):
    return Dist(sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))/sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**S(2) + c**S(2)))), Int(sqrt(a/(a + sqrt(b**S(2) + c**S(2))) + sqrt(b**S(2) + c**S(2))*cos(d + e*x - ArcTan(b, c))/(a + sqrt(b**S(2) + c**S(2)))), x), x)


def replacement3080(a, b, c, d, e, n, x):
    return Dist(S(1)/n, Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-2))*Simp(a**S(2)*n + a*b*(S(2)*n + S(-1))*cos(d + e*x) + a*c*(S(2)*n + S(-1))*sin(d + e*x) + (b**S(2) + c**S(2))*(n + S(-1)), x), x), x) - Simp((-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1))/(e*n), x)


def With3081(a, b, c, d, e, x):
    f = FreeFactors(S(1)/tan(d/S(2) + e*x/S(2)), x)
    return -Dist(f/e, Subst(Int(S(1)/(a + c*f*x), x), x, S(1)/(f*tan(d/S(2) + e*x/S(2)))), x)


def With3082(a, b, c, d, e, x):
    f = FreeFactors(tan(Pi/S(4) + d/S(2) + e*x/S(2)), x)
    return Dist(f/e, Subst(Int(S(1)/(a + b*f*x), x), x, tan(Pi/S(4) + d/S(2) + e*x/S(2))/f), x)


def With3083(a, b, c, d, e, x):
    f = FreeFactors(S(1)/tan(Pi/S(4) + d/S(2) + e*x/S(2)), x)
    return -Dist(f/e, Subst(Int(S(1)/(a + b*f*x), x), x, S(1)/(f*tan(Pi/S(4) + d/S(2) + e*x/S(2)))), x)


def With3084(a, b, c, d, e, x):
    f = FreeFactors(tan(d/S(2) + e*x/S(2)), x)
    return Dist(S(2)*f/e, Subst(Int(S(1)/(a + b + S(2)*c*f*x + f**S(2)*x**S(2)*(a - b)), x), x, tan(d/S(2) + e*x/S(2))/f), x)


def replacement3085(a, b, c, d, e, x):
    return Dist(b/(c*e), Subst(Int(S(1)/(x*sqrt(a + x)), x), x, b*cos(d + e*x) + c*sin(d + e*x)), x)


def replacement3086(a, b, c, d, e, x):
    return Int(S(1)/sqrt(a + sqrt(b**S(2) + c**S(2))*cos(d + e*x - ArcTan(b, c))), x)


def replacement3087(a, b, c, d, e, x):
    return Dist(sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**S(2) + c**S(2))))/sqrt(a + b*cos(d + e*x) + c*sin(d + e*x)), Int(S(1)/sqrt(a/(a + sqrt(b**S(2) + c**S(2))) + sqrt(b**S(2) + c**S(2))*cos(d + e*x - ArcTan(b, c))/(a + sqrt(b**S(2) + c**S(2)))), x), x)


def replacement3088(a, b, c, d, e, x):
    return Dist(S(1)/(a**S(2) - b**S(2) - c**S(2)), Int(sqrt(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Simp(S(2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3089(a, b, c, d, e, n, x):
    return Dist(S(1)/((n + S(1))*(a**S(2) - b**S(2) - c**S(2))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*(a*(n + S(1)) - b*(n + S(2))*cos(d + e*x) - c*(n + S(2))*sin(d + e*x)), x), x) + Simp((b*sin(d + e*x) - c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))/(e*(n + S(1))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3090(A, B, C, a, b, c, d, e, x):
    return Simp(x*(S(2)*A*a - B*b - C*c)/(S(2)*a**S(2)), x) + Simp((-S(2)*A*a*b**S(2) + a**S(2)*(B*b - C*c) + b**S(2)*(B*b + C*c))*log(RemoveContent(a + b*cos(d + e*x) + c*sin(d + e*x), x))/(S(2)*a**S(2)*b*c*e), x) - Simp((B*b + C*c)*(b*cos(d + e*x) - c*sin(d + e*x))/(S(2)*a*b*c*e), x)


def replacement3091(A, C, a, b, c, d, e, x):
    return Simp(x*(S(2)*A*a - C*c)/(S(2)*a**S(2)), x) - Simp(C*cos(d + e*x)/(S(2)*a*e), x) + Simp((S(2)*A*a*c - C*a**S(2) + C*b**S(2))*log(RemoveContent(a + b*cos(d + e*x) + c*sin(d + e*x), x))/(S(2)*a**S(2)*b*e), x) + Simp(C*c*sin(d + e*x)/(S(2)*a*b*e), x)


def replacement3092(A, B, a, b, c, d, e, x):
    return Simp(x*(S(2)*A*a - B*b)/(S(2)*a**S(2)), x) + Simp(B*sin(d + e*x)/(S(2)*a*e), x) + Simp((-S(2)*A*a*b + B*a**S(2) + B*b**S(2))*log(RemoveContent(a + b*cos(d + e*x) + c*sin(d + e*x), x))/(S(2)*a**S(2)*c*e), x) - Simp(B*b*cos(d + e*x)/(S(2)*a*c*e), x)


def replacement3093(A, B, C, a, b, c, d, e, x):
    return Simp(x*(B*b + C*c)/(b**S(2) + c**S(2)), x) + Simp((B*c - C*b)*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x)


def replacement3094(A, C, a, b, c, d, e, x):
    return Simp(C*c*x/(b**S(2) + c**S(2)), x) - Simp(C*b*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x)


def replacement3095(A, B, a, b, c, d, e, x):
    return Simp(B*b*x/(b**S(2) + c**S(2)), x) + Simp(B*c*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x)


def replacement3096(A, B, C, a, b, c, d, e, x):
    return Dist((A*(b**S(2) + c**S(2)) - a*(B*b + C*c))/(b**S(2) + c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Simp(x*(B*b + C*c)/(b**S(2) + c**S(2)), x) + Simp((B*c - C*b)*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x)


def replacement3097(A, C, a, b, c, d, e, x):
    return Dist((A*(b**S(2) + c**S(2)) - C*a*c)/(b**S(2) + c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) - Simp(C*b*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x) + Simp(C*c*(d + e*x)/(e*(b**S(2) + c**S(2))), x)


def replacement3098(A, B, a, b, c, d, e, x):
    return Dist((A*(b**S(2) + c**S(2)) - B*a*b)/(b**S(2) + c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Simp(B*b*(d + e*x)/(e*(b**S(2) + c**S(2))), x) + Simp(B*c*log(a + b*cos(d + e*x) + c*sin(d + e*x))/(e*(b**S(2) + c**S(2))), x)


def replacement3099(A, B, C, a, b, c, d, e, n, x):
    return Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**n*(B*a*sin(d + e*x) + B*c - C*a*cos(d + e*x) - C*b)/(a*e*(n + S(1))), x)


def replacement3100(A, C, a, b, c, d, e, n, x):
    return -Simp((C*a*cos(d + e*x) + C*b)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3101(A, B, a, b, c, d, e, n, x):
    return Simp((B*a*sin(d + e*x) + B*c)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3102(A, B, C, a, b, c, d, e, n, x):
    return Dist((A*a*(n + S(1)) + n*(B*b + C*c))/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**n, x), x) + Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**n*(B*a*sin(d + e*x) + B*c - C*a*cos(d + e*x) - C*b)/(a*e*(n + S(1))), x)


def replacement3103(A, C, a, b, c, d, e, n, x):
    return Dist((A*a*(n + S(1)) + C*c*n)/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**n, x), x) - Simp((C*a*cos(d + e*x) + C*b)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3104(A, B, a, b, c, d, e, n, x):
    return Dist((A*a*(n + S(1)) + B*b*n)/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**n, x), x) + Simp((B*a*sin(d + e*x) + B*c)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3105(B, C, b, c, d, e, n, x):
    return Simp((B*c - C*b)*(b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))/(e*(b**S(2) + c**S(2))*(n + S(1))), x)


def replacement3106(A, B, C, a, b, c, d, e, n, x):
    return Dist(S(1)/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1))*Simp(A*a**S(2)*(n + S(1)) + a*n*(B*b + C*c) + (A*a*b*(n + S(1)) + n*(B*a**S(2) - B*c**S(2) + C*b*c))*cos(d + e*x) + (A*a*c*(n + S(1)) + n*(B*b*c + C*a**S(2) - C*b**S(2)))*sin(d + e*x), x), x), x) + Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**n*(B*a*sin(d + e*x) + B*c - C*a*cos(d + e*x) - C*b)/(a*e*(n + S(1))), x)


def replacement3107(A, C, a, b, c, d, e, n, x):
    return Dist(S(1)/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1))*Simp(A*a**S(2)*(n + S(1)) + C*a*c*n + (A*a*b*(n + S(1)) + C*b*c*n)*cos(d + e*x) + (A*a*c*(n + S(1)) + C*a**S(2)*n - C*b**S(2)*n)*sin(d + e*x), x), x), x) - Simp((C*a*cos(d + e*x) + C*b)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3108(A, B, a, b, c, d, e, n, x):
    return Dist(S(1)/(a*(n + S(1))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(-1))*Simp(A*a**S(2)*(n + S(1)) + B*a*b*n + (A*a*c*(n + S(1)) + B*b*c*n)*sin(d + e*x) + (A*a*b*(n + S(1)) + B*a**S(2)*n - B*c**S(2)*n)*cos(d + e*x), x), x), x) + Simp((B*a*sin(d + e*x) + B*c)*(a + b*cos(d + e*x) + c*sin(d + e*x))**n/(a*e*(n + S(1))), x)


def replacement3109(A, B, C, a, b, c, d, e, x):
    return Dist(B/b, Int(sqrt(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Dist((A*b - B*a)/b, Int(S(1)/sqrt(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x)


def replacement3110(A, B, C, a, b, c, d, e, x):
    return Simp((B*c - C*b + (-A*b + B*a)*sin(d + e*x) - (-A*c + C*a)*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3111(A, C, a, b, c, d, e, x):
    return -Simp((A*b*sin(d + e*x) + C*b + (-A*c + C*a)*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3112(A, B, a, b, c, d, e, x):
    return Simp((A*c*cos(d + e*x) + B*c + (-A*b + B*a)*sin(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3113(A, B, C, a, b, c, d, e, x):
    return Dist((A*a - B*b - C*c)/(a**S(2) - b**S(2) - c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Simp((B*c - C*b + (-A*b + B*a)*sin(d + e*x) - (-A*c + C*a)*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3114(A, C, a, b, c, d, e, x):
    return Dist((A*a - C*c)/(a**S(2) - b**S(2) - c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) - Simp((A*b*sin(d + e*x) + C*b + (-A*c + C*a)*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3115(A, B, a, b, c, d, e, x):
    return Dist((A*a - B*b)/(a**S(2) - b**S(2) - c**S(2)), Int(S(1)/(a + b*cos(d + e*x) + c*sin(d + e*x)), x), x) + Simp((A*c*cos(d + e*x) + B*c + (-A*b + B*a)*sin(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3116(A, B, C, a, b, c, d, e, n, x):
    return Dist(S(1)/((n + S(1))*(a**S(2) - b**S(2) - c**S(2))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*Simp((n + S(1))*(A*a - B*b - C*c) + (n + S(2))*(-A*b + B*a)*cos(d + e*x) + (n + S(2))*(-A*c + C*a)*sin(d + e*x), x), x), x) - Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*(B*c - C*b + (-A*b + B*a)*sin(d + e*x) - (-A*c + C*a)*cos(d + e*x))/(e*(n + S(1))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3117(A, C, a, b, c, d, e, n, x):
    return Dist(S(1)/((n + S(1))*(a**S(2) - b**S(2) - c**S(2))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*Simp(-A*b*(n + S(2))*cos(d + e*x) + (n + S(1))*(A*a - C*c) + (n + S(2))*(-A*c + C*a)*sin(d + e*x), x), x), x) + Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*(A*b*sin(d + e*x) + C*b + (-A*c + C*a)*cos(d + e*x))/(e*(n + S(1))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3118(A, B, a, b, c, d, e, n, x):
    return Dist(S(1)/((n + S(1))*(a**S(2) - b**S(2) - c**S(2))), Int((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*Simp(-A*c*(n + S(2))*sin(d + e*x) + (n + S(1))*(A*a - B*b) + (n + S(2))*(-A*b + B*a)*cos(d + e*x), x), x), x) - Simp((a + b*cos(d + e*x) + c*sin(d + e*x))**(n + S(1))*(A*c*cos(d + e*x) + B*c + (-A*b + B*a)*sin(d + e*x))/(e*(n + S(1))*(a**S(2) - b**S(2) - c**S(2))), x)


def replacement3119(a, b, c, d, e, x):
    return Int(cos(d + e*x)/(a*cos(d + e*x) + b + c*sin(d + e*x)), x)


def replacement3120(a, b, c, d, e, x):
    return Int(sin(d + e*x)/(a*sin(d + e*x) + b + c*cos(d + e*x)), x)


def replacement3121(a, b, c, d, e, n, x):
    return Int((a*cos(d + e*x) + b + c*sin(d + e*x))**n, x)


def replacement3122(a, b, c, d, e, n, x):
    return Int((a*sin(d + e*x) + b + c*cos(d + e*x))**n, x)


def replacement3123(a, b, c, d, e, n, x):
    return Dist((a + b/cos(d + e*x) + c*tan(d + e*x))**n*(a*cos(d + e*x) + b + c*sin(d + e*x))**(-n)*cos(d + e*x)**n, Int((a*cos(d + e*x) + b + c*sin(d + e*x))**n, x), x)


def replacement3124(a, b, c, d, e, n, x):
    return Dist((a + b/sin(d + e*x) + c/tan(d + e*x))**n*(a*sin(d + e*x) + b + c*cos(d + e*x))**(-n)*sin(d + e*x)**n, Int((a*sin(d + e*x) + b + c*cos(d + e*x))**n, x), x)


def replacement3125(a, b, c, d, e, m, n, x):
    return Int((a*cos(d + e*x) + b + c*sin(d + e*x))**(-n), x)


def replacement3126(a, b, c, d, e, m, n, x):
    return Int((a*sin(d + e*x) + b + c*cos(d + e*x))**(-n), x)


def replacement3127(a, b, c, d, e, m, n, x):
    return Dist((a + b/cos(d + e*x) + c*tan(d + e*x))**(-n)*(a*cos(d + e*x) + b + c*sin(d + e*x))**n*(S(1)/cos(d + e*x))**n, Int((a*cos(d + e*x) + b + c*sin(d + e*x))**(-n), x), x)


def replacement3128(a, b, c, d, e, m, n, x):
    return Dist((a + b/sin(d + e*x) + c/tan(d + e*x))**(-n)*(a*sin(d + e*x) + b + c*cos(d + e*x))**n*(S(1)/sin(d + e*x))**n, Int((a*sin(d + e*x) + b + c*cos(d + e*x))**(-n), x), x)


def replacement3129(b, c, d, p, x):
    return Dist(b*(S(2)*p + S(-1))/(S(2)*p), Int((b*sin(c + d*x)**S(2))**(p + S(-1)), x), x) - Simp((b*sin(c + d*x)**S(2))**p/(S(2)*d*p*tan(c + d*x)), x)


def replacement3130(b, c, d, p, x):
    return Dist(b*(S(2)*p + S(-1))/(S(2)*p), Int((b*cos(c + d*x)**S(2))**(p + S(-1)), x), x) + Simp((b*cos(c + d*x)**S(2))**p*tan(c + d*x)/(S(2)*d*p), x)


def replacement3131(b, c, d, p, x):
    return Dist(S(2)*(p + S(1))/(b*(S(2)*p + S(1))), Int((b*sin(c + d*x)**S(2))**(p + S(1)), x), x) + Simp((b*sin(c + d*x)**S(2))**(p + S(1))/(b*d*(S(2)*p + S(1))*tan(c + d*x)), x)


def replacement3132(b, c, d, p, x):
    return Dist(S(2)*(p + S(1))/(b*(S(2)*p + S(1))), Int((b*cos(c + d*x)**S(2))**(p + S(1)), x), x) - Simp((b*cos(c + d*x)**S(2))**(p + S(1))*tan(c + d*x)/(b*d*(S(2)*p + S(1))), x)


def replacement3133(a, b, c, d, p, x):
    return Dist(a**p, Int(cos(c + d*x)**(S(2)*p), x), x)


def replacement3134(a, b, c, d, p, x):
    return Dist(a**p, Int(sin(c + d*x)**(S(2)*p), x), x)


def replacement3135(a, b, c, d, p, x):
    return Int((a*cos(c + d*x)**S(2))**p, x)


def replacement3136(a, b, c, d, p, x):
    return Int((a*sin(c + d*x)**S(2))**p, x)


def With3137(a, b, c, d, p, x):
    e = FreeFactors(tan(c + d*x), x)
    return Dist(e/d, Subst(Int((a + e**S(2)*x**S(2)*(a + b))**p*(e**S(2)*x**S(2) + S(1))**(-p + S(-1)), x), x, tan(c + d*x)/e), x)


def With3138(a, b, c, d, p, x):
    e = FreeFactors(tan(c + d*x), x)
    return Dist(e/d, Subst(Int((e**S(2)*x**S(2) + S(1))**(-p + S(-1))*(a*e**S(2)*x**S(2) + a + b)**p, x), x, tan(c + d*x)/e), x)


def replacement3139(a, b, c, d, p, x):
    return Dist(S(2)**(-p), Int((S(2)*a - b*cos(S(2)*c + S(2)*d*x) + b)**p, x), x)


def replacement3140(a, b, c, d, p, x):
    return Dist(S(2)**(-p), Int((S(2)*a + b*cos(S(2)*c + S(2)*d*x) + b)**p, x), x)


def replacement3141(a, b, c, d, n, p, x):
    return Int((a + b*sin(c + d*x)**n)**p, x)


def replacement3142(a, b, c, d, n, p, x):
    return Int((a + b*cos(c + d*x)**n)**p, x)


def With3143(a, b, c, d, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def With3144(a, b, c, d, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, tan(c + d*x)/f), x)


def replacement3145(a, b, c, d, p, u, x):
    return Dist(a**p, Int(ActivateTrig(u)*cos(c + d*x)**(S(2)*p), x), x)


def replacement3146(a, b, c, d, p, u, x):
    return Dist(a**p, Int(ActivateTrig(u)*sin(c + d*x)**(S(2)*p), x), x)


def replacement3147(a, b, c, d, p, u, x):
    return Dist((a*cos(c + d*x)**S(2))**p*cos(c + d*x)**(-S(2)*p), Int(ActivateTrig(u)*cos(c + d*x)**(S(2)*p), x), x)


def replacement3148(a, b, c, d, p, u, x):
    return Dist((a*sin(c + d*x)**S(2))**p*sin(c + d*x)**(-S(2)*p), Int(ActivateTrig(u)*sin(c + d*x)**(S(2)*p), x), x)


def With3149(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def With3150(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b, x)**p, x), x, tan(c + d*x)/f), x)


def With3151(a, b, c, d, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a + b*(-f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, cos(c + d*x)/f), x)


def With3152(a, b, c, d, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f/d, Subst(Int((-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*ExpandToSum(a + b*(-f**S(2)*x**S(2) + S(1))**(n/S(2)), x)**p, x), x, sin(c + d*x)/f), x)


def replacement3153(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((a + b*sin(c + d*x)**n)**p*sin(c + d*x)**m, x), x)


def replacement3154(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((a + b*cos(c + d*x)**n)**p*cos(c + d*x)**m, x), x)


def With3155(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x)**p, x), x, tan(c + d*x)/f), x)


def With3156(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f/d, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def replacement3157(a, b, c, d, m, n, p, x):
    return Int((S(1) - sin(c + d*x)**S(2))**(m/S(2))*(a + b*sin(c + d*x)**n)**p, x)


def replacement3158(a, b, c, d, m, n, p, x):
    return Int((S(1) - cos(c + d*x)**S(2))**(m/S(2))*(a + b*cos(c + d*x)**n)**p, x)


def replacement3159(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((S(1) - sin(c + d*x)**S(2))**(m/S(2))*(a + b*sin(c + d*x)**n)**p, x), x)


def replacement3160(a, b, c, d, m, n, p, x):
    return Int(ExpandTrig((S(1) - cos(c + d*x)**S(2))**(m/S(2))*(a + b*cos(c + d*x)**n)**p, x), x)


def With3161(a, b, c, d, e, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f/d, Subst(Int((a + b*(e*f*x)**n)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, sin(c + d*x)/f), x)


def With3162(a, b, c, d, e, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f/d, Subst(Int((a + b*(e*f*x)**n)**p*(-f**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2), x), x, cos(c + d*x)/f), x)


def With3163(a, b, c, d, e, m, n, p, x):
    f = FreeFactors(sin(c + d*x), x)
    return Dist(f**(m + S(1))/d, Subst(Int(x**m*(a + b*(e*f*x)**n)**p*(-f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)/2), x), x, sin(c + d*x)/f), x)


def With3164(a, b, c, d, e, m, n, p, x):
    f = FreeFactors(cos(c + d*x), x)
    return -Dist(f**(m + S(1))/d, Subst(Int(x**m*(a + b*(e*f*x)**n)**p*(-f**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)/2), x), x, cos(c + d*x)/f), x)


def With3165(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f**(m + S(1))/d, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x)**p, x), x, tan(c + d*x)/f), x)


def With3166(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f**(m + S(1))/d, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x)**p, x), x, S(1)/(f*tan(c + d*x))), x)


def With3167(a, b, c, d, m, n, p, x):
    f = FreeFactors(tan(c + d*x), x)
    return Dist(f**(m + S(1))/d, Subst(Int(x**m*((f**S(2)*x**S(2) + S(1))**(-n/S(2))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x))**p/(f**S(2)*x**S(2) + S(1)), x), x, tan(c + d*x)/f), x)


def With3168(a, b, c, d, m, n, p, x):
    f = FreeFactors(S(1)/tan(c + d*x), x)
    return -Dist(f**(m + S(1))/d, Subst(Int(x**m*((f**S(2)*x**S(2) + S(1))**(-n/S(2))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + b*f**n*x**n, x))**p/(f**S(2)*x**S(2) + S(1)), x), x, S(1)/(f*tan(c + d*x))), x)


def With3169(a, b, c, d, e, m, n, p, q, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*q/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(q/S(2)) + b*(f**S(2)*x**S(2) + S(1))**(-p/S(2) + q/S(2)) + c, x)**n, x), x, S(1)/(f*tan(d + e*x))), x)


def With3170(a, b, c, d, e, m, n, p, q, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*q/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(q/S(2)) + b*(f**S(2)*x**S(2) + S(1))**(-p/S(2) + q/S(2)) + c, x)**n, x), x, tan(d + e*x)/f), x)


def With3171(a, b, c, d, e, m, n, p, q, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(p/S(2)) + b*f**p*x**p + c*(f**S(2)*x**S(2) + S(1))**(p/S(2) - q/S(2)), x)**n, x), x, S(1)/(f*tan(d + e*x))), x)


def With3172(a, b, c, d, e, m, n, p, q, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p/S(2) + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**(p/S(2)) + b*f**p*x**p + c*(f**S(2)*x**S(2) + S(1))**(p/S(2) - q/S(2)), x)**n, x), x, tan(d + e*x)/f), x)


def replacement3173(a, b, c, d, e, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p), x), x)


def replacement3174(a, b, c, d, e, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p), x), x)


def replacement3175(a, b, c, d, e, n, n2, p, x):
    return Dist((b + S(2)*c*sin(d + e*x)**n)**(-S(2)*p)*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, Int(u*(b + S(2)*c*sin(d + e*x)**n)**(S(2)*p), x), x)


def replacement3176(a, b, c, d, e, n, n2, p, x):
    return Dist((b + S(2)*c*cos(d + e*x)**n)**(-S(2)*p)*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, Int(u*(b + S(2)*c*cos(d + e*x)**n)**(S(2)*p), x), x)


def With3177(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*sin(d + e*x)**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*sin(d + e*x)**n + q), x), x)


def With3178(a, b, c, d, e, n, n2, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*cos(d + e*x)**n - q), x), x) - Dist(S(2)*c/q, Int(S(1)/(b + S(2)*c*cos(d + e*x)**n + q), x), x)


def replacement3179(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*sin(d + e*x)**m, x), x)


def replacement3180(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*cos(d + e*x)**m, x), x)


def replacement3181(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*sin(d + e*x)**n)**(-S(2)*p)*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*sin(d + e*x)**m, x), x)


def replacement3182(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*cos(d + e*x)**n)**(-S(2)*p)*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*cos(d + e*x)**m, x), x)


def With3183(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*(x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, S(1)/(f*tan(d + e*x))), x)


def With3184(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f/e, Subst(Int((f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*(x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, tan(d + e*x)/f), x)


def replacement3185(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p*sin(d + e*x)**m, x), x)


def replacement3186(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p*cos(d + e*x)**m, x), x)


def With3187(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(sin(d + e*x), x)
    return Dist(g/e, Subst(Int((-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, sin(d + e*x)/g), x)


def With3188(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(cos(d + e*x), x)
    return -Dist(g/e, Subst(Int((-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, cos(d + e*x)/g), x)


def replacement3189(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*cos(d + e*x)**m, x), x)


def replacement3190(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*sin(d + e*x)**m, x), x)


def replacement3191(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*sin(d + e*x)**n)**(-S(2)*p)*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*cos(d + e*x)**m, x), x)


def replacement3192(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*cos(d + e*x)**n)**(-S(2)*p)*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*sin(d + e*x)**m, x), x)


def With3193(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*(x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, S(1)/(f*tan(d + e*x))), x)


def With3194(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-m/S(2) - n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*(x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, tan(d + e*x)/f), x)


def replacement3195(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - sin(d + e*x)**S(2))**(m/S(2))*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, x), x)


def replacement3196(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - cos(d + e*x)**S(2))**(m/S(2))*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, x), x)


def With3197(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(sin(d + e*x), x)
    return Dist(g**(m + S(1))/e, Subst(Int(x**m*(-g**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, sin(d + e*x)/g), x)


def With3198(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(cos(d + e*x), x)
    return -Dist(g**(m + S(1))/e, Subst(Int(x**m*(-g**S(2)*x**S(2) + S(1))**(-m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, cos(d + e*x)/g), x)


def replacement3199(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)


def replacement3200(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)


def replacement3201(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*sin(d + e*x)**n)**(-S(2)*p)*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)


def replacement3202(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*cos(d + e*x)**n)**(-S(2)*p)*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)


def With3203(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*x**n*(x**S(2) + S(1))**(n/S(2)) + c*x**(S(2)*n), x)**p, x), x, tan(d + e*x)/f), x)


def With3204(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p + S(-1))*ExpandToSum(a*(x**S(2) + S(1))**n + b*x**n*(x**S(2) + S(1))**(n/S(2)) + c*x**(S(2)*n), x)**p, x), x, S(1)/(f*tan(d + e*x))), x)


def replacement3205(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - sin(d + e*x)**S(2))**(-m/S(2))*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p*sin(d + e*x)**m, x), x)


def replacement3206(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - cos(d + e*x)**S(2))**(-m/S(2))*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p*cos(d + e*x)**m, x), x)


def With3207(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(sin(d + e*x), x)
    return Dist(g**(m + S(1))/e, Subst(Int(x**(-m)*(-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, sin(d + e*x)/g), x)


def With3208(a, b, c, d, e, f, m, n, n2, p, x):
    g = FreeFactors(cos(d + e*x), x)
    return -Dist(g**(m + S(1))/e, Subst(Int(x**(-m)*(-g**S(2)*x**S(2) + S(1))**(m/S(2) + S(-1)/2)*(a + b*(f*g*x)**n + c*(f*g*x)**(S(2)*n))**p, x), x, cos(d + e*x)/g), x)


def replacement3209(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)


def replacement3210(a, b, c, d, e, m, n, n2, p, x):
    return Dist(S(4)**(-p)*c**(-p), Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)


def replacement3211(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*sin(d + e*x)**n)**(-S(2)*p)*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*sin(d + e*x)**n)**(S(2)*p)*(S(1)/tan(d + e*x))**m, x), x)


def replacement3212(a, b, c, d, e, m, n, n2, p, x):
    return Dist((b + S(2)*c*cos(d + e*x)**n)**(-S(2)*p)*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p, Int((b + S(2)*c*cos(d + e*x)**n)**(S(2)*p)*tan(d + e*x)**m, x), x)


def With3213(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(S(1)/tan(d + e*x), x)
    return -Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**n + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, S(1)/(f*tan(d + e*x))), x)


def With3214(a, b, c, d, e, m, n, n2, p, x):
    f = FreeFactors(tan(d + e*x), x)
    return Dist(f**(m + S(1))/e, Subst(Int(x**m*(f**S(2)*x**S(2) + S(1))**(-n*p + S(-1))*ExpandToSum(a*(f**S(2)*x**S(2) + S(1))**n + b*(f**S(2)*x**S(2) + S(1))**(n/S(2)) + c, x)**p, x), x, tan(d + e*x)/f), x)


def replacement3215(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - sin(d + e*x)**S(2))**(m/S(2))*(a + b*sin(d + e*x)**n + c*sin(d + e*x)**(S(2)*n))**p*sin(d + e*x)**(-m), x), x)


def replacement3216(a, b, c, d, e, m, n, n2, p, x):
    return Int(ExpandTrig((S(1) - cos(d + e*x)**S(2))**(m/S(2))*(a + b*cos(d + e*x)**n + c*cos(d + e*x)**(S(2)*n))**p*cos(d + e*x)**(-m), x), x)


def replacement3217(A, B, a, b, c, d, e, n, x):
    return Dist(S(4)**(-n)*c**(-n), Int((A + B*sin(d + e*x))*(b + S(2)*c*sin(d + e*x))**(S(2)*n), x), x)


def replacement3218(A, B, a, b, c, d, e, n, x):
    return Dist(S(4)**(-n)*c**(-n), Int((A + B*cos(d + e*x))*(b + S(2)*c*cos(d + e*x))**(S(2)*n), x), x)


def replacement3219(A, B, a, b, c, d, e, n, x):
    return Dist((b + S(2)*c*sin(d + e*x))**(-S(2)*n)*(a + b*sin(d + e*x) + c*sin(d + e*x)**S(2))**n, Int((A + B*sin(d + e*x))*(b + S(2)*c*sin(d + e*x))**(S(2)*n), x), x)


def replacement3220(A, B, a, b, c, d, e, n, x):
    return Dist((b + S(2)*c*cos(d + e*x))**(-S(2)*n)*(a + b*cos(d + e*x) + c*cos(d + e*x)**S(2))**n, Int((A + B*cos(d + e*x))*(b + S(2)*c*cos(d + e*x))**(S(2)*n), x), x)


def With3221(A, B, a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c*sin(d + e*x) - q), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c*sin(d + e*x) + q), x), x)


def With3222(A, B, a, b, c, d, e, x):
    q = Rt(-S(4)*a*c + b**S(2), S(2))
    return Dist(B - (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c*cos(d + e*x) - q), x), x) + Dist(B + (-S(2)*A*c + B*b)/q, Int(S(1)/(b + S(2)*c*cos(d + e*x) + q), x), x)


def replacement3223(A, B, a, b, c, d, e, n, x):
    return Int(ExpandTrig((A + B*sin(d + e*x))*(a + b*sin(d + e*x) + c*sin(d + e*x)**S(2))**n, x), x)


def replacement3224(A, B, a, b, c, d, e, n, x):
    return Int(ExpandTrig((A + B*cos(d + e*x))*(a + b*cos(d + e*x) + c*cos(d + e*x)**S(2))**n, x), x)


def replacement3225(c, d, e, f, m, x):
    return Dist(d*m/f, Int((c + d*x)**(m + S(-1))*cos(e + f*x), x), x) - Simp((c + d*x)**m*cos(e + f*x)/f, x)


def replacement3226(c, d, e, f, m, x):
    return -Dist(d*m/f, Int((c + d*x)**(m + S(-1))*sin(e + f*x), x), x) + Simp((c + d*x)**m*sin(e + f*x)/f, x)


def replacement3227(c, d, e, f, m, x):
    return -Dist(f/(d*(m + S(1))), Int((c + d*x)**(m + S(1))*cos(e + f*x), x), x) + Simp((c + d*x)**(m + S(1))*sin(e + f*x)/(d*(m + S(1))), x)


def replacement3228(c, d, e, f, m, x):
    return Dist(f/(d*(m + S(1))), Int((c + d*x)**(m + S(1))*sin(e + f*x), x), x) + Simp((c + d*x)**(m + S(1))*cos(e + f*x)/(d*(m + S(1))), x)


def replacement3229(c, d, e, f, x):
    return Simp(SinIntegral(e + f*x)/d, x)


def replacement3230(c, d, e, f, x):
    return Simp(CosIntegral(e + f*x)/d, x)


def replacement3231(c, d, e, f, x):
    return Dist(sin((-c*f + d*e)/d), Int(cos(c*f/d + f*x)/(c + d*x), x), x) + Dist(cos((-c*f + d*e)/d), Int(sin(c*f/d + f*x)/(c + d*x), x), x)


def replacement3232(c, d, e, f, x):
    return -Dist(sin((-c*f + d*e)/d), Int(sin(c*f/d + f*x)/(c + d*x), x), x) + Dist(cos((-c*f + d*e)/d), Int(cos(c*f/d + f*x)/(c + d*x), x), x)


def replacement3233(c, d, e, f, x):
    return Dist(S(2)/d, Subst(Int(sin(f*x**S(2)/d), x), x, sqrt(c + d*x)), x)


def replacement3234(c, d, e, f, x):
    return Dist(S(2)/d, Subst(Int(cos(f*x**S(2)/d), x), x, sqrt(c + d*x)), x)


def replacement3235(c, d, e, f, x):
    return Dist(sin((-c*f + d*e)/d), Int(cos(c*f/d + f*x)/sqrt(c + d*x), x), x) + Dist(cos((-c*f + d*e)/d), Int(sin(c*f/d + f*x)/sqrt(c + d*x), x), x)


def replacement3236(c, d, e, f, x):
    return -Dist(sin((-c*f + d*e)/d), Int(sin(c*f/d + f*x)/sqrt(c + d*x), x), x) + Dist(cos((-c*f + d*e)/d), Int(cos(c*f/d + f*x)/sqrt(c + d*x), x), x)


def replacement3237(c, d, e, f, m, x):
    return Dist(I/S(2), Int((c + d*x)**m*exp(-I*(e + f*x)), x), x) - Dist(I/S(2), Int((c + d*x)**m*exp(I*(e + f*x)), x), x)


def replacement3238(c, d, e, f, m, x):
    return Dist(S(1)/2, Int((c + d*x)**m*exp(-I*(e + f*x)), x), x) + Dist(S(1)/2, Int((c + d*x)**m*exp(I*(e + f*x)), x), x)


def replacement3239(b, c, d, e, f, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*sin(e + f*x))**(n + S(-2))*(c + d*x), x), x) + Simp(d*(b*sin(e + f*x))**n/(f**S(2)*n**S(2)), x) - Simp(b*(b*sin(e + f*x))**(n + S(-1))*(c + d*x)*cos(e + f*x)/(f*n), x)


def replacement3240(b, c, d, e, f, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*cos(e + f*x))**(n + S(-2))*(c + d*x), x), x) + Simp(d*(b*cos(e + f*x))**n/(f**S(2)*n**S(2)), x) + Simp(b*(b*cos(e + f*x))**(n + S(-1))*(c + d*x)*sin(e + f*x)/(f*n), x)


def replacement3241(b, c, d, e, f, m, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*sin(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) - Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b*sin(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) - Simp(b*(b*sin(e + f*x))**(n + S(-1))*(c + d*x)**m*cos(e + f*x)/(f*n), x) + Simp(d*m*(b*sin(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)


def replacement3242(b, c, d, e, f, m, n, x):
    return Dist(b**S(2)*(n + S(-1))/n, Int((b*cos(e + f*x))**(n + S(-2))*(c + d*x)**m, x), x) - Dist(d**S(2)*m*(m + S(-1))/(f**S(2)*n**S(2)), Int((b*cos(e + f*x))**n*(c + d*x)**(m + S(-2)), x), x) + Simp(b*(b*cos(e + f*x))**(n + S(-1))*(c + d*x)**m*sin(e + f*x)/(f*n), x) + Simp(d*m*(b*cos(e + f*x))**n*(c + d*x)**(m + S(-1))/(f**S(2)*n**S(2)), x)


def replacement3243(c, d, e, f, m, n, x):
    return Int(ExpandTrigReduce((c + d*x)**m, sin(e + f*x)**n, x), x)


def replacement3244(c, d, e, f, m, n, x):
    return Int(ExpandTrigReduce((c + d*x)**m, cos(e + f*x)**n, x), x)


def replacement3245(c, d, e, f, m, n, x):
    return -Dist(f*n/(d*(m + S(1))), Int(ExpandTrigReduce((c + d*x)**(m + S(1)), sin(e + f*x)**(n + S(-1))*cos(e + f*x), x), x), x) + Simp((c + d*x)**(m + S(1))*sin(e + f*x)**n/(d*(m + S(1))), x)


def replacement3246(c, d, e, f, m, n, x):
    return Dist(f*n/(d*(m + S(1))), Int(ExpandTrigReduce((c + d*x)**(m + S(1)), sin(e + f*x)*cos(e + f*x)**(n + S(-1)), x), x), x) + Simp((c + d*x)**(m + S(1))*cos(e + f*x)**n/(d*(m + S(1))), x)


def replacement3247(b, c, d, e, f, m, n, x):
    return -Dist(f**S(2)*n**S(2)/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*sin(e + f*x))**n*(c + d*x)**(m + S(2)), x), x) + Dist(b**S(2)*f**S(2)*n*(n + S(-1))/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*sin(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x), x) + Simp((b*sin(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x) - Simp(b*f*n*(b*sin(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*cos(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))), x)


def replacement3248(b, c, d, e, f, m, n, x):
    return -Dist(f**S(2)*n**S(2)/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*cos(e + f*x))**n*(c + d*x)**(m + S(2)), x), x) + Dist(b**S(2)*f**S(2)*n*(n + S(-1))/(d**S(2)*(m + S(1))*(m + S(2))), Int((b*cos(e + f*x))**(n + S(-2))*(c + d*x)**(m + S(2)), x), x) + Simp((b*cos(e + f*x))**n*(c + d*x)**(m + S(1))/(d*(m + S(1))), x) + Simp(b*f*n*(b*cos(e + f*x))**(n + S(-1))*(c + d*x)**(m + S(2))*sin(e + f*x)/(d**S(2)*(m + S(1))*(m + S(2))), x)


def replacement3249(b, c, d, e, f, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*sin(e + f*x))**(n + S(2))*(c + d*x), x), x) - Simp(d*(b*sin(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x) + Simp((b*sin(e + f*x))**(n + S(1))*(c + d*x)*cos(e + f*x)/(b*f*(n + S(1))), x)


def replacement3250(b, c, d, e, f, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*cos(e + f*x))**(n + S(2))*(c + d*x), x), x) - Simp(d*(b*cos(e + f*x))**(n + S(2))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x) - Simp((b*cos(e + f*x))**(n + S(1))*(c + d*x)*sin(e + f*x)/(b*f*(n + S(1))), x)


def replacement3251(b, c, d, e, f, m, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*sin(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), Int((b*sin(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x), x) + Simp((b*sin(e + f*x))**(n + S(1))*(c + d*x)**m*cos(e + f*x)/(b*f*(n + S(1))), x) - Simp(d*m*(b*sin(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x)


def replacement3252(b, c, d, e, f, m, n, x):
    return Dist((n + S(2))/(b**S(2)*(n + S(1))), Int((b*cos(e + f*x))**(n + S(2))*(c + d*x)**m, x), x) + Dist(d**S(2)*m*(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), Int((b*cos(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-2)), x), x) - Simp((b*cos(e + f*x))**(n + S(1))*(c + d*x)**m*sin(e + f*x)/(b*f*(n + S(1))), x) - Simp(d*m*(b*cos(e + f*x))**(n + S(2))*(c + d*x)**(m + S(-1))/(b**S(2)*f**S(2)*(n + S(1))*(n + S(2))), x)


def replacement3253(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a + b*sin(e + f*x))**n, x), x)


def replacement3254(a, b, c, d, e, f, m, n, x):
    return Int(ExpandIntegrand((c + d*x)**m, (a + b*cos(e + f*x))**n, x), x)


def replacement3255(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**n, Int((c + d*x)**m*cos(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3256(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**IntPart(n)*(a + b*sin(e + f*x))**FracPart(n)*cos(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*cos(-Pi*a/(S(4)*b) + e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3257(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**n, Int((c + d*x)**m*cos(e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3258(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**n, Int((c + d*x)**m*sin(e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3259(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**IntPart(n)*(a + b*cos(e + f*x))**FracPart(n)*cos(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*cos(e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3260(a, b, c, d, e, f, m, n, x):
    return Dist((S(2)*a)**IntPart(n)*(a + b*cos(e + f*x))**FracPart(n)*sin(e/S(2) + f*x/S(2))**(-S(2)*FracPart(n)), Int((c + d*x)**m*sin(e/S(2) + f*x/S(2))**(S(2)*n), x), x)


def replacement3261(a, b, c, d, e, f, m, x):
    return Dist(S(2), Int((c + d*x)**m*exp(I*(e + f*x))/(S(2)*a*exp(I*(e + f*x)) - I*b*exp(S(2)*I*(e + f*x)) + I*b), x), x)


def replacement3262(a, b, c, d, e, f, m, x):
    return Dist(S(2), Int((c + d*x)**m*exp(I*(e + f*x))/(S(2)*a*exp(I*(e + f*x)) + b*exp(S(2)*I*(e + f*x)) + b), x), x)


def replacement3263(a, b, c, d, e, f, m, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((c + d*x)**m/(a + b*sin(e + f*x)), x), x) - Dist(b*d*m/(f*(a**S(2) - b**S(2))), Int((c + d*x)**(m + S(-1))*cos(e + f*x)/(a + b*sin(e + f*x)), x), x) + Simp(b*(c + d*x)**m*cos(e + f*x)/(f*(a + b*sin(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement3264(a, b, c, d, e, f, m, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((c + d*x)**m/(a + b*cos(e + f*x)), x), x) + Dist(b*d*m/(f*(a**S(2) - b**S(2))), Int((c + d*x)**(m + S(-1))*sin(e + f*x)/(a + b*cos(e + f*x)), x), x) - Simp(b*(c + d*x)**m*sin(e + f*x)/(f*(a + b*cos(e + f*x))*(a**S(2) - b**S(2))), x)


def replacement3265(a, b, c, d, e, f, m, n, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((a + b*sin(e + f*x))**(n + S(1))*(c + d*x)**m, x), x) - Dist(b*(n + S(2))/((a**S(2) - b**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(n + S(1))*(c + d*x)**m*sin(e + f*x), x), x) + Dist(b*d*m/(f*(a**S(2) - b**S(2))*(n + S(1))), Int((a + b*sin(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*cos(e + f*x), x), x) - Simp(b*(a + b*sin(e + f*x))**(n + S(1))*(c + d*x)**m*cos(e + f*x)/(f*(a**S(2) - b**S(2))*(n + S(1))), x)


def replacement3266(a, b, c, d, e, f, m, n, x):
    return Dist(a/(a**S(2) - b**S(2)), Int((a + b*cos(e + f*x))**(n + S(1))*(c + d*x)**m, x), x) - Dist(b*(n + S(2))/((a**S(2) - b**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(n + S(1))*(c + d*x)**m*cos(e + f*x), x), x) - Dist(b*d*m/(f*(a**S(2) - b**S(2))*(n + S(1))), Int((a + b*cos(e + f*x))**(n + S(1))*(c + d*x)**(m + S(-1))*sin(e + f*x), x), x) + Simp(b*(a + b*cos(e + f*x))**(n + S(1))*(c + d*x)**m*sin(e + f*x)/(f*(a**S(2) - b**S(2))*(n + S(1))), x)


def replacement3267(a, b, m, n, u, v, x):
    return Int((a + b*sin(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)


def replacement3268(a, b, m, n, u, v, x):
    return Int((a + b*cos(ExpandToSum(v, x)))**n*ExpandToSum(u, x)**m, x)


def replacement3269(a, b, c, d, e, f, m, n, x):
    return Int((a + b*sin(e + f*x))**n*(c + d*x)**m, x)


def replacement3270(a, b, c, d, e, f, m, n, x):
    return Int((a + b*cos(e + f*x))**n*(c + d*x)**m, x)


def replacement3271(a, b, c, d, n, p, x):
    return Int(ExpandIntegrand(sin(c + d*x), (a + b*x**n)**p, x), x)


def replacement3272(a, b, c, d, n, p, x):
    return Int(ExpandIntegrand(cos(c + d*x), (a + b*x**n)**p, x), x)


def replacement3273(a, b, c, d, n, p, x):
    return -Dist(d/(b*n*(p + S(1))), Int(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*cos(c + d*x), x), x) - Dist((S(1) - n)/(b*n*(p + S(1))), Int(x**(-n)*(a + b*x**n)**(p + S(1))*sin(c + d*x), x), x) + Simp(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*sin(c + d*x)/(b*n*(p + S(1))), x)


def replacement3274(a, b, c, d, n, p, x):
    return Dist(d/(b*n*(p + S(1))), Int(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*sin(c + d*x), x), x) - Dist((S(1) - n)/(b*n*(p + S(1))), Int(x**(-n)*(a + b*x**n)**(p + S(1))*cos(c + d*x), x), x) + Simp(x**(S(1) - n)*(a + b*x**n)**(p + S(1))*cos(c + d*x)/(b*n*(p + S(1))), x)


def replacement3275(a, b, c, d, n, p, x):
    return Int(ExpandIntegrand(sin(c + d*x), (a + b*x**n)**p, x), x)


def replacement3276(a, b, c, d, n, p, x):
    return Int(ExpandIntegrand(cos(c + d*x), (a + b*x**n)**p, x), x)


def replacement3277(a, b, c, d, n, p, x):
    return Int(x**(n*p)*(a*x**(-n) + b)**p*sin(c + d*x), x)


def replacement3278(a, b, c, d, n, p, x):
    return Int(x**(n*p)*(a*x**(-n) + b)**p*cos(c + d*x), x)


def replacement3279(a, b, c, d, n, p, x):
    return Int((a + b*x**n)**p*sin(c + d*x), x)


def replacement3280(a, b, c, d, n, p, x):
    return Int((a + b*x**n)**p*cos(c + d*x), x)


def replacement3281(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(sin(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x)


def replacement3282(a, b, c, d, e, m, n, p, x):
    return Int(ExpandIntegrand(cos(c + d*x), (e*x)**m*(a + b*x**n)**p, x), x)


def replacement3283(a, b, c, d, e, m, n, p, x):
    return -Dist(d*e**m/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*cos(c + d*x), x), x) + Simp(e**m*(a + b*x**n)**(p + S(1))*sin(c + d*x)/(b*n*(p + S(1))), x)


def replacement3284(a, b, c, d, e, m, n, p, x):
    return Dist(d*e**m/(b*n*(p + S(1))), Int((a + b*x**n)**(p + S(1))*sin(c + d*x), x), x) + Simp(e**m*(a + b*x**n)**(p + S(1))*cos(c + d*x)/(b*n*(p + S(1))), x)


def replacement3285(a, b, c, d, m, n, p, x):
    return -Dist(d/(b*n*(p + S(1))), Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*cos(c + d*x), x), x) - Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*(a + b*x**n)**(p + S(1))*sin(c + d*x), x), x) + Simp(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sin(c + d*x)/(b*n*(p + S(1))), x)


def replacement3286(a, b, c, d, m, n, p, x):
    return Dist(d/(b*n*(p + S(1))), Int(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*sin(c + d*x), x), x) - Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*(a + b*x**n)**(p + S(1))*cos(c + d*x), x), x) + Simp(x**(m - n + S(1))*(a + b*x**n)**(p + S(1))*cos(c + d*x)/(b*n*(p + S(1))), x)


def replacement3287(a, b, c, d, m, n, p, x):
    return Int(ExpandIntegrand(sin(c + d*x), x**m*(a + b*x**n)**p, x), x)


def replacement3288(a, b, c, d, m, n, p, x):
    return Int(ExpandIntegrand(cos(c + d*x), x**m*(a + b*x**n)**p, x), x)


def replacement3289(a, b, c, d, m, n, p, x):
    return Int(x**(m + n*p)*(a*x**(-n) + b)**p*sin(c + d*x), x)


def replacement3290(a, b, c, d, m, n, p, x):
    return Int(x**(m + n*p)*(a*x**(-n) + b)**p*cos(c + d*x), x)


def replacement3291(a, b, c, d, e, m, n, p, x):
    return Int((e*x)**m*(a + b*x**n)**p*sin(c + d*x), x)


def replacement3292(a, b, c, d, e, m, n, p, x):
    return Int((e*x)**m*(a + b*x**n)**p*cos(c + d*x), x)


def replacement3293(d, x):
    return Simp(sqrt(S(2))*sqrt(Pi)*FresnelS(sqrt(S(2))*x*sqrt(S(1)/Pi)*Rt(d, S(2)))/(S(2)*Rt(d, S(2))), x)


def replacement3294(d, x):
    return Simp(sqrt(S(2))*sqrt(Pi)*FresnelC(sqrt(S(2))*x*sqrt(S(1)/Pi)*Rt(d, S(2)))/(S(2)*Rt(d, S(2))), x)


def replacement3295(c, d, x):
    return Dist(sin(c), Int(cos(d*x**S(2)), x), x) + Dist(cos(c), Int(sin(d*x**S(2)), x), x)


def replacement3296(c, d, x):
    return -Dist(sin(c), Int(sin(d*x**S(2)), x), x) + Dist(cos(c), Int(cos(d*x**S(2)), x), x)


def replacement3297(c, d, n, x):
    return Dist(I/S(2), Int(exp(-I*c - I*d*x**n), x), x) - Dist(I/S(2), Int(exp(I*c + I*d*x**n), x), x)


def replacement3298(c, d, n, x):
    return Dist(S(1)/2, Int(exp(-I*c - I*d*x**n), x), x) + Dist(S(1)/2, Int(exp(I*c + I*d*x**n), x), x)


def replacement3299(a, b, c, d, n, p, x):
    return Int(ExpandTrigReduce((a + b*sin(c + d*x**n))**p, x), x)


def replacement3300(a, b, c, d, n, p, x):
    return Int(ExpandTrigReduce((a + b*cos(c + d*x**n))**p, x), x)


def replacement3301(a, b, c, d, n, p, x):
    return -Subst(Int((a + b*sin(c + d*x**(-n)))**p/x**S(2), x), x, S(1)/x)


def replacement3302(a, b, c, d, n, p, x):
    return -Subst(Int((a + b*cos(c + d*x**(-n)))**p/x**S(2), x), x, S(1)/x)


def With3303(a, b, c, d, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*sin(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)


def With3304(a, b, c, d, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k + S(-1))*(a + b*cos(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)


def replacement3305(c, d, n, x):
    return Dist(I/S(2), Int(exp(-I*c - I*d*x**n), x), x) - Dist(I/S(2), Int(exp(I*c + I*d*x**n), x), x)


def replacement3306(c, d, n, x):
    return Dist(S(1)/2, Int(exp(-I*c - I*d*x**n), x), x) + Dist(S(1)/2, Int(exp(I*c + I*d*x**n), x), x)


def replacement3307(a, b, c, d, n, p, x):
    return Int(ExpandTrigReduce((a + b*sin(c + d*x**n))**p, x), x)


def replacement3308(a, b, c, d, n, p, x):
    return Int(ExpandTrigReduce((a + b*cos(c + d*x**n))**p, x), x)


def replacement3309(a, b, c, d, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*sin(c + d*x**n))**p, x), x, u), x)


def replacement3310(a, b, c, d, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*cos(c + d*x**n))**p, x), x, u), x)


def replacement3311(a, b, c, d, n, p, u, x):
    return Int((a + b*sin(c + d*u**n))**p, x)


def replacement3312(a, b, c, d, n, p, u, x):
    return Int((a + b*cos(c + d*u**n))**p, x)


def replacement3313(a, b, p, u, x):
    return Int((a + b*sin(ExpandToSum(u, x)))**p, x)


def replacement3314(a, b, p, u, x):
    return Int((a + b*cos(ExpandToSum(u, x)))**p, x)


def replacement3315(d, n, x):
    return Simp(SinIntegral(d*x**n)/n, x)


def replacement3316(d, n, x):
    return Simp(CosIntegral(d*x**n)/n, x)


def replacement3317(c, d, n, x):
    return Dist(sin(c), Int(cos(d*x**n)/x, x), x) + Dist(cos(c), Int(sin(d*x**n)/x, x), x)


def replacement3318(c, d, n, x):
    return -Dist(sin(c), Int(sin(d*x**n)/x, x), x) + Dist(cos(c), Int(cos(d*x**n)/x, x), x)


def With3319(a, b, c, d, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    mn = (m + S(1))/n
    if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
        return True
    return False


def replacement3319(a, b, c, d, m, n, p, x):

    mn = (m + S(1))/n
    return Dist(S(1)/n, Subst(Int(x**(mn + S(-1))*(a + b*sin(c + d*x))**p, x), x, x**n), x)


def With3320(a, b, c, d, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    mn = (m + S(1))/n
    if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
        return True
    return False


def replacement3320(a, b, c, d, m, n, p, x):

    mn = (m + S(1))/n
    return Dist(S(1)/n, Subst(Int(x**(mn + S(-1))*(a + b*cos(c + d*x))**p, x), x, x**n), x)


def With3321(a, b, c, d, e, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    mn = (m + S(1))/n
    if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
        return True
    return False


def replacement3321(a, b, c, d, e, m, n, p, x):

    mn = (m + S(1))/n
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sin(c + d*x**n))**p, x), x)


def With3322(a, b, c, d, e, m, n, p, x):
    if isinstance(x, (int, Integer, float, Float)):
        return False
    mn = (m + S(1))/n
    if And(IntegerQ(mn), Or(Equal(p, S(1)), Greater(mn, S(0)))):
        return True
    return False


def replacement3322(a, b, c, d, e, m, n, p, x):

    mn = (m + S(1))/n
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cos(c + d*x**n))**p, x), x)


def replacement3323(a, b, m, n, x):
    return Dist(S(2)/n, Subst(Int(sin(a + b*x**S(2)), x), x, x**(n/S(2))), x)


def replacement3324(a, b, m, n, x):
    return Dist(S(2)/n, Subst(Int(cos(a + b*x**S(2)), x), x, x**(n/S(2))), x)


def replacement3325(c, d, e, m, n, x):
    return Dist(e**n*(m - n + S(1))/(d*n), Int((e*x)**(m - n)*cos(c + d*x**n), x), x) - Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*cos(c + d*x**n)/(d*n), x)


def replacement3326(c, d, e, m, n, x):
    return -Dist(e**n*(m - n + S(1))/(d*n), Int((e*x)**(m - n)*sin(c + d*x**n), x), x) + Simp(e**(n + S(-1))*(e*x)**(m - n + S(1))*sin(c + d*x**n)/(d*n), x)


def replacement3327(c, d, e, m, n, x):
    return -Dist(d*e**(-n)*n/(m + S(1)), Int((e*x)**(m + n)*cos(c + d*x**n), x), x) + Simp((e*x)**(m + S(1))*sin(c + d*x**n)/(e*(m + S(1))), x)


def replacement3328(c, d, e, m, n, x):
    return Dist(d*e**(-n)*n/(m + S(1)), Int((e*x)**(m + n)*sin(c + d*x**n), x), x) + Simp((e*x)**(m + S(1))*cos(c + d*x**n)/(e*(m + S(1))), x)


def replacement3329(c, d, e, m, n, x):
    return Dist(I/S(2), Int((e*x)**m*exp(-I*c - I*d*x**n), x), x) - Dist(I/S(2), Int((e*x)**m*exp(I*c + I*d*x**n), x), x)


def replacement3330(c, d, e, m, n, x):
    return Dist(S(1)/2, Int((e*x)**m*exp(-I*c - I*d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(I*c + I*d*x**n), x), x)


def replacement3331(a, b, m, n, p, x):
    return Dist(b*n*p/(n + S(-1)), Int(sin(a + b*x**n)**(p + S(-1))*cos(a + b*x**n), x), x) - Simp(x**(S(1) - n)*sin(a + b*x**n)**p/(n + S(-1)), x)


def replacement3332(a, b, m, n, p, x):
    return -Dist(b*n*p/(n + S(-1)), Int(sin(a + b*x**n)*cos(a + b*x**n)**(p + S(-1)), x), x) - Simp(x**(S(1) - n)*cos(a + b*x**n)**p/(n + S(-1)), x)


def replacement3333(a, b, m, n, p, x):
    return Dist((p + S(-1))/p, Int(x**m*sin(a + b*x**n)**(p + S(-2)), x), x) + Simp(sin(a + b*x**n)**p/(b**S(2)*n*p**S(2)), x) - Simp(x**n*sin(a + b*x**n)**(p + S(-1))*cos(a + b*x**n)/(b*n*p), x)


def replacement3334(a, b, m, n, p, x):
    return Dist((p + S(-1))/p, Int(x**m*cos(a + b*x**n)**(p + S(-2)), x), x) + Simp(cos(a + b*x**n)**p/(b**S(2)*n*p**S(2)), x) + Simp(x**n*sin(a + b*x**n)*cos(a + b*x**n)**(p + S(-1))/(b*n*p), x)


def replacement3335(a, b, m, n, p, x):
    return Dist((p + S(-1))/p, Int(x**m*sin(a + b*x**n)**(p + S(-2)), x), x) - Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*p**S(2)), Int(x**(m - S(2)*n)*sin(a + b*x**n)**p, x), x) + Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*sin(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)), x) - Simp(x**(m - n + S(1))*sin(a + b*x**n)**(p + S(-1))*cos(a + b*x**n)/(b*n*p), x)


def replacement3336(a, b, m, n, p, x):
    return Dist((p + S(-1))/p, Int(x**m*cos(a + b*x**n)**(p + S(-2)), x), x) - Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*p**S(2)), Int(x**(m - S(2)*n)*cos(a + b*x**n)**p, x), x) + Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*cos(a + b*x**n)**p/(b**S(2)*n**S(2)*p**S(2)), x) + Simp(x**(m - n + S(1))*sin(a + b*x**n)*cos(a + b*x**n)**(p + S(-1))/(b*n*p), x)


def replacement3337(a, b, m, n, p, x):
    return -Dist(b**S(2)*n**S(2)*p**S(2)/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*sin(a + b*x**n)**p, x), x) + Dist(b**S(2)*n**S(2)*p*(p + S(-1))/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*sin(a + b*x**n)**(p + S(-2)), x), x) + Simp(x**(m + S(1))*sin(a + b*x**n)**p/(m + S(1)), x) - Simp(b*n*p*x**(m + n + S(1))*sin(a + b*x**n)**(p + S(-1))*cos(a + b*x**n)/((m + S(1))*(m + n + S(1))), x)


def replacement3338(a, b, m, n, p, x):
    return -Dist(b**S(2)*n**S(2)*p**S(2)/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*cos(a + b*x**n)**p, x), x) + Dist(b**S(2)*n**S(2)*p*(p + S(-1))/((m + S(1))*(m + n + S(1))), Int(x**(m + S(2)*n)*cos(a + b*x**n)**(p + S(-2)), x), x) + Simp(x**(m + S(1))*cos(a + b*x**n)**p/(m + S(1)), x) + Simp(b*n*p*x**(m + n + S(1))*sin(a + b*x**n)*cos(a + b*x**n)**(p + S(-1))/((m + S(1))*(m + n + S(1))), x)


def With3339(a, b, c, d, e, m, n, p, x):
    k = Denominator(m)
    return Dist(k/e, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sin(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(S(1)/k)), x)


def With3340(a, b, c, d, e, m, n, p, x):
    k = Denominator(m)
    return Dist(k/e, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*cos(c + d*e**(-n)*x**(k*n)))**p, x), x, (e*x)**(S(1)/k)), x)


def replacement3341(a, b, c, d, e, m, n, p, x):
    return Int(ExpandTrigReduce((e*x)**m, (a + b*sin(c + d*x**n))**p, x), x)


def replacement3342(a, b, c, d, e, m, n, p, x):
    return Int(ExpandTrigReduce((e*x)**m, (a + b*cos(c + d*x**n))**p, x), x)


def replacement3343(a, b, m, n, p, x):
    return Dist((p + S(2))/(p + S(1)), Int(x**m*sin(a + b*x**n)**(p + S(2)), x), x) - Simp(sin(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))), x) + Simp(x**n*sin(a + b*x**n)**(p + S(1))*cos(a + b*x**n)/(b*n*(p + S(1))), x)


def replacement3344(a, b, m, n, p, x):
    return Dist((p + S(2))/(p + S(1)), Int(x**m*cos(a + b*x**n)**(p + S(2)), x), x) - Simp(cos(a + b*x**n)**(p + S(2))/(b**S(2)*n*(p + S(1))*(p + S(2))), x) - Simp(x**n*sin(a + b*x**n)*cos(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement3345(a, b, m, n, p, x):
    return Dist((p + S(2))/(p + S(1)), Int(x**m*sin(a + b*x**n)**(p + S(2)), x), x) + Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**(m - S(2)*n)*sin(a + b*x**n)**(p + S(2)), x), x) + Simp(x**(m - n + S(1))*sin(a + b*x**n)**(p + S(1))*cos(a + b*x**n)/(b*n*(p + S(1))), x) - Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*sin(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)


def replacement3346(a, b, m, n, p, x):
    return Dist((p + S(2))/(p + S(1)), Int(x**m*cos(a + b*x**n)**(p + S(2)), x), x) + Dist((m - S(2)*n + S(1))*(m - n + S(1))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), Int(x**(m - S(2)*n)*cos(a + b*x**n)**(p + S(2)), x), x) - Simp(x**(m - n + S(1))*sin(a + b*x**n)*cos(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x) - Simp(x**(m - S(2)*n + S(1))*(m - n + S(1))*cos(a + b*x**n)**(p + S(2))/(b**S(2)*n**S(2)*(p + S(1))*(p + S(2))), x)


def replacement3347(a, b, c, d, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*sin(c + d*x**(-n)))**p, x), x, S(1)/x)


def replacement3348(a, b, c, d, m, n, p, x):
    return -Subst(Int(x**(-m + S(-2))*(a + b*cos(c + d*x**(-n)))**p, x), x, S(1)/x)


def With3349(a, b, c, d, e, m, n, p, x):
    k = Denominator(m)
    return -Dist(k/e, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*sin(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k)), x)


def With3350(a, b, c, d, e, m, n, p, x):
    k = Denominator(m)
    return -Dist(k/e, Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*cos(c + d*e**(-n)*x**(-k*n)))**p, x), x, (e*x)**(-S(1)/k)), x)


def replacement3351(a, b, c, d, e, m, n, p, x):
    return -Dist((e*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*sin(c + d*x**(-n)))**p, x), x, S(1)/x), x)


def replacement3352(a, b, c, d, e, m, n, p, x):
    return -Dist((e*x)**m*(S(1)/x)**m, Subst(Int(x**(-m + S(-2))*(a + b*cos(c + d*x**(-n)))**p, x), x, S(1)/x), x)


def With3353(a, b, c, d, m, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*sin(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)


def With3354(a, b, c, d, m, n, p, x):
    k = Denominator(n)
    return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*cos(c + d*x**(k*n)))**p, x), x, x**(S(1)/k)), x)


def replacement3355(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sin(c + d*x**n))**p, x), x)


def replacement3356(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cos(c + d*x**n))**p, x), x)


def replacement3357(a, b, c, d, m, n, p, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*sin(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1))), x)


def replacement3358(a, b, c, d, m, n, p, x):
    return Dist(S(1)/(m + S(1)), Subst(Int((a + b*cos(c + d*x**(n/(m + S(1)))))**p, x), x, x**(m + S(1))), x)


def replacement3359(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*sin(c + d*x**n))**p, x), x)


def replacement3360(a, b, c, d, e, m, n, p, x):
    return Dist(e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m), Int(x**m*(a + b*cos(c + d*x**n))**p, x), x)


def replacement3361(c, d, e, m, n, x):
    return Dist(I/S(2), Int((e*x)**m*exp(-I*c - I*d*x**n), x), x) - Dist(I/S(2), Int((e*x)**m*exp(I*c + I*d*x**n), x), x)


def replacement3362(c, d, e, m, n, x):
    return Dist(S(1)/2, Int((e*x)**m*exp(-I*c - I*d*x**n), x), x) + Dist(S(1)/2, Int((e*x)**m*exp(I*c + I*d*x**n), x), x)


def replacement3363(a, b, c, d, e, m, n, p, x):
    return Int(ExpandTrigReduce((e*x)**m, (a + b*sin(c + d*x**n))**p, x), x)


def replacement3364(a, b, c, d, e, m, n, p, x):
    return Int(ExpandTrigReduce((e*x)**m, (a + b*cos(c + d*x**n))**p, x), x)


def replacement3365(a, b, c, d, m, n, p, u, x):
    return Dist(Coefficient(u, x, S(1))**(-m + S(-1)), Subst(Int((a + b*sin(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u), x)


def replacement3366(a, b, c, d, m, n, p, u, x):
    return Dist(Coefficient(u, x, S(1))**(-m + S(-1)), Subst(Int((a + b*cos(c + d*x**n))**p*(x - Coefficient(u, x, S(0)))**m, x), x, u), x)


def replacement3367(a, b, c, d, e, m, n, p, u, x):
    return Int((e*x)**m*(a + b*sin(c + d*u**n))**p, x)


def replacement3368(a, b, c, d, e, m, n, p, u, x):
    return Int((e*x)**m*(a + b*cos(c + d*u**n))**p, x)


def replacement3369(a, b, e, m, p, u, x):
    return Int((e*x)**m*(a + b*sin(ExpandToSum(u, x)))**p, x)


def replacement3370(a, b, e, m, p, u, x):
    return Int((e*x)**m*(a + b*cos(ExpandToSum(u, x)))**p, x)


def replacement3371(a, b, m, n, p, x):
    return Simp(sin(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement3372(a, b, m, n, p, x):
    return -Simp(cos(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement3373(a, b, m, n, p, x):
    return -Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*sin(a + b*x**n)**(p + S(1)), x), x) + Simp(x**(m - n + S(1))*sin(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement3374(a, b, m, n, p, x):
    return Dist((m - n + S(1))/(b*n*(p + S(1))), Int(x**(m - n)*cos(a + b*x**n)**(p + S(1)), x), x) - Simp(x**(m - n + S(1))*cos(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))), x)


def replacement3375(a, b, c, x):
    return Int(sin((b + S(2)*c*x)**S(2)/(S(4)*c)), x)


def replacement3376(a, b, c, x):
    return Int(cos((b + S(2)*c*x)**S(2)/(S(4)*c)), x)


def replacement3377(a, b, c, x):
    return -Dist(sin((-S(4)*a*c + b**S(2))/(S(4)*c)), Int(cos((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x) + Dist(cos((-S(4)*a*c + b**S(2))/(S(4)*c)), Int(sin((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x)


def replacement3378(a, b, c, x):
    return Dist(sin((-S(4)*a*c + b**S(2))/(S(4)*c)), Int(sin((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x) + Dist(cos((-S(4)*a*c + b**S(2))/(S(4)*c)), Int(cos((b + S(2)*c*x)**S(2)/(S(4)*c)), x), x)


def replacement3379(a, b, c, n, x):
    return Int(ExpandTrigReduce(sin(a + b*x + c*x**S(2))**n, x), x)


def replacement3380(a, b, c, n, x):
    return Int(ExpandTrigReduce(cos(a + b*x + c*x**S(2))**n, x), x)


def replacement3381(n, v, x):
    return Int(sin(ExpandToSum(v, x))**n, x)


def replacement3382(n, v, x):
    return Int(cos(ExpandToSum(v, x))**n, x)


def replacement3383(a, b, c, d, e, x):
    return -Simp(e*cos(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3384(a, b, c, d, e, x):
    return Simp(e*sin(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3385(a, b, c, d, e, x):
    return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(sin(a + b*x + c*x**S(2)), x), x) - Simp(e*cos(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3386(a, b, c, d, e, x):
    return Dist((-b*e + S(2)*c*d)/(S(2)*c), Int(cos(a + b*x + c*x**S(2)), x), x) + Simp(e*sin(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3387(a, b, c, d, e, m, x):
    return Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*cos(a + b*x + c*x**S(2)), x), x) - Simp(e*(d + e*x)**(m + S(-1))*cos(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3388(a, b, c, d, e, m, x):
    return -Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*sin(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*sin(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3389(a, b, c, d, e, m, x):
    return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int((d + e*x)**(m + S(-1))*sin(a + b*x + c*x**S(2)), x), x) + Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*cos(a + b*x + c*x**S(2)), x), x) - Simp(e*(d + e*x)**(m + S(-1))*cos(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3390(a, b, c, d, e, m, x):
    return -Dist((b*e - S(2)*c*d)/(S(2)*c), Int((d + e*x)**(m + S(-1))*cos(a + b*x + c*x**S(2)), x), x) - Dist(e**S(2)*(m + S(-1))/(S(2)*c), Int((d + e*x)**(m + S(-2))*sin(a + b*x + c*x**S(2)), x), x) + Simp(e*(d + e*x)**(m + S(-1))*sin(a + b*x + c*x**S(2))/(S(2)*c), x)


def replacement3391(a, b, c, d, e, m, x):
    return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*cos(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*sin(a + b*x + c*x**S(2))/(e*(m + S(1))), x)


def replacement3392(a, b, c, d, e, m, x):
    return Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*sin(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*cos(a + b*x + c*x**S(2))/(e*(m + S(1))), x)


def replacement3393(a, b, c, d, e, m, x):
    return -Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*cos(a + b*x + c*x**S(2)), x), x) - Dist((b*e - S(2)*c*d)/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(1))*cos(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*sin(a + b*x + c*x**S(2))/(e*(m + S(1))), x)


def replacement3394(a, b, c, d, e, m, x):
    return Dist(S(2)*c/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(2))*sin(a + b*x + c*x**S(2)), x), x) + Dist((b*e - S(2)*c*d)/(e**S(2)*(m + S(1))), Int((d + e*x)**(m + S(1))*sin(a + b*x + c*x**S(2)), x), x) + Simp((d + e*x)**(m + S(1))*cos(a + b*x + c*x**S(2))/(e*(m + S(1))), x)


def replacement3395(a, b, c, d, e, m, x):
    return Int((d + e*x)**m*sin(a + b*x + c*x**S(2)), x)


def replacement3396(a, b, c, d, e, m, x):
    return Int((d + e*x)**m*cos(a + b*x + c*x**S(2)), x)


def replacement3397(a, b, c, d, e, m, n, x):
    return Int(ExpandTrigReduce((d + e*x)**m, sin(a + b*x + c*x**S(2))**n, x), x)


def replacement3398(a, b, c, d, e, m, n, x):
    return Int(ExpandTrigReduce((d + e*x)**m, cos(a + b*x + c*x**S(2))**n, x), x)


def replacement3399(m, n, u, v, x):
    return Int(ExpandToSum(u, x)**m*sin(ExpandToSum(v, x))**n, x)


def replacement3400(m, n, u, v, x):
    return Int(ExpandToSum(u, x)**m*cos(ExpandToSum(v, x))**n, x)
