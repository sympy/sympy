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
        PolynomialRemainder, Factor, PolyLog, CosIntegral, SinIntegral, LogIntegral, SinhIntegral,
        CoshIntegral, Rule, Erf, PolyGamma, ExpIntegralEi, ExpIntegralE, LogGamma , UtilityOperator,
        DerivativeDivides
    )
    from sympy import Integral, S, sqrt, And, Or, Integer, Float, Mod, I, Abs 
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii , Pqq, Q, R, r, C = symbols('i ii Pqq Q R r C_1')
    _UseGamma = False

def miscellaneous_integration(rubi):
    from sympy.integrals.rubi.constraints import cons147, cons1297, cons2, cons3, cons7, cons4, cons5, cons386, cons27, cons50, cons1298, cons1299, cons1300, cons1301, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1302, cons66, cons21, cons84, cons1037, cons1036, cons38, cons1303, cons10, cons1304, cons1305, cons1306, cons209, cons1256, cons1240, cons1307, cons46, cons1308, cons1309, cons1310, cons1311, cons52, cons1312, cons800, cons1313, cons17, cons1314, cons586, cons1315, cons1316, cons1317, cons1318, cons1319, cons1320, cons1321, cons1322, cons1323, cons667, cons196, cons1324, cons840, cons1325, cons18, cons1326, cons148, cons45, cons1327, cons1328, cons1243, cons261, cons1329, cons367, cons1330, cons67, cons1331, cons744, cons1332, cons165, cons1333, cons1334, cons1335, cons1257, cons1336, cons347

    pattern2337 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons147, cons1297)
    def replacement2337(u, c, b, x, p, n, a):
        rubi.append(2337)
        return Dist(c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p)), Int(u*(a + b*x)**(n*p), x), x)
    rule2337 = ReplacementRule(pattern2337, replacement2337)
    pattern2338 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons147, cons386)
    def replacement2338(u, c, b, d, x, p, a, q):
        rubi.append(2338)
        return Dist((c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q), Int(u*(a + b*x)**(p*q), x), x)
    rule2338 = ReplacementRule(pattern2338, replacement2338)
    pattern2339 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons147, cons386)
    def replacement2339(u, c, d, b, x, p, n, a, q):
        rubi.append(2339)
        return Dist((c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q), Int(u*(a + b*x)**(n*p*q), x), x)
    rule2339 = ReplacementRule(pattern2339, replacement2339)
    pattern2340 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1298, cons1299, cons1300, cons1301)
    def replacement2340(F, c, b, B, d, A, e, x, f, C, n, g, a):
        rubi.append(2340)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule2340 = ReplacementRule(pattern2340, replacement2340)
    pattern2341 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1091, cons1298, cons1302)
    def replacement2341(F, c, b, e, A, x, C, n, g, a):
        rubi.append(2341)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule2341 = ReplacementRule(pattern2341, replacement2341)
    pattern2342 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1298, cons1299, cons1300, cons1301)
    def replacement2342(F, c, b, d, e, A, B, x, f, C, n, g, a):
        rubi.append(2342)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule2342 = ReplacementRule(pattern2342, replacement2342)
    pattern2343 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1091, cons1298, cons1302)
    def replacement2343(F, c, b, e, A, x, C, n, g, a):
        rubi.append(2343)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule2343 = ReplacementRule(pattern2343, replacement2343)
    def With2344(u, y, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2344 = Pattern(Integral(u_/y_, x_), CustomConstraint(With2344))
    def replacement2344(u, y, x):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2344)
        return Simp(q*log(RemoveContent(y, x)), x)
    rule2344 = ReplacementRule(pattern2344, replacement2344)
    def With2345(u, y, w, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(w*y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2345 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(With2345))
    def replacement2345(u, y, w, x):
        
        q = DerivativeDivides(w*y, u, x)
        rubi.append(2345)
        return Simp(q*log(RemoveContent(w*y, x)), x)
    rule2345 = ReplacementRule(pattern2345, replacement2345)
    def With2346(u, y, m, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2346 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons21, cons66, CustomConstraint(With2346))
    def replacement2346(u, y, m, x):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2346)
        return Simp(q*y**(m + S(1))/(m + S(1)), x)
    rule2346 = ReplacementRule(pattern2346, replacement2346)
    def With2347(u, z, x, y, n, m):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2347 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons21, cons4, cons66, CustomConstraint(With2347))
    def replacement2347(u, z, x, y, n, m):
        
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        rubi.append(2347)
        return Simp(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), x)
    rule2347 = ReplacementRule(pattern2347, replacement2347)
    def With2348(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = SimplifyIntegrand(u, x)
        if SimplerIntegrandQ(v, u, x):
            return True
        return False
    pattern2348 = Pattern(Integral(u_, x_), CustomConstraint(With2348))
    def replacement2348(u, x):
        
        v = SimplifyIntegrand(u, x)
        rubi.append(2348)
        return Int(v, x)
    rule2348 = ReplacementRule(pattern2348, replacement2348)
    pattern2349 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1037)
    def replacement2349(u, c, b, e, d, x, f, n, m, a):
        rubi.append(2349)
        return Dist((a*e**S(2) - c*f**S(2))**m, Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule2349 = ReplacementRule(pattern2349, replacement2349)
    pattern2350 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1036)
    def replacement2350(u, c, b, e, d, x, f, n, m, a):
        rubi.append(2350)
        return Dist((b*e**S(2) - d*f**S(2))**m, Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule2350 = ReplacementRule(pattern2350, replacement2350)
    pattern2351 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), cons2, cons21, cons4, cons38, cons1303, cons10)
    def replacement2351(u, x, p, w, v, n, m, a):
        rubi.append(2351)
        return Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x)
    rule2351 = ReplacementRule(pattern2351, replacement2351)
    def With2352(u, c, b, d, x, y, v, n, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2352 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1304, CustomConstraint(With2352))
    def replacement2352(u, c, b, d, x, y, v, n, m, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2352)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), x)
    rule2352 = ReplacementRule(pattern2352, replacement2352)
    def With2353(u, c, d, e, b, x, f, p, w, y, v, n, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2353 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons1304, cons1305, CustomConstraint(With2353))
    def replacement2353(u, c, d, e, b, x, f, p, w, y, v, n, m, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2353)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), x)
    rule2353 = ReplacementRule(pattern2353, replacement2353)
    def With2354(u, c, d, e, b, z, g, x, f, p, w, y, v, n, h, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern2354 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons1304, cons1305, cons1306, CustomConstraint(With2354))
    def replacement2354(u, c, d, e, b, z, g, x, f, p, w, y, v, n, h, m, a, q):
        
        r = DerivativeDivides(y, u, x)
        rubi.append(2354)
        return Dist(r, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), x)
    rule2354 = ReplacementRule(pattern2354, replacement2354)
    def With2355(u, b, x, y, n, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2355 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1256, CustomConstraint(With2355))
    def replacement2355(u, b, x, y, n, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2355)
        return Dist(a, Int(u, x), x) + Dist(b*q, Subst(Int(x**n, x), x, y), x)
    rule2355 = ReplacementRule(pattern2355, replacement2355)
    def With2356(u, b, x, y, p, n, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2356 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1240, CustomConstraint(With2356))
    def replacement2356(u, b, x, y, p, n, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2356)
        return Dist(q, Subst(Int((a + b*x**n)**p, x), x, y), x)
    rule2356 = ReplacementRule(pattern2356, replacement2356)
    def With2357(u, b, x, y, p, v, n, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern2357 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons1307, CustomConstraint(With2357))
    def replacement2357(u, b, x, y, p, v, n, m, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(2357)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n)**p, x), x, y), x)
    rule2357 = ReplacementRule(pattern2357, replacement2357)
    def With2358(u, c, b, x, y, p, v, n, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2358 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons1304, CustomConstraint(With2358))
    def replacement2358(u, c, b, x, y, p, v, n, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2358)
        return Dist(q, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule2358 = ReplacementRule(pattern2358, replacement2358)
    def With2359(u, c, b, B, A, x, y, p, w, v, n, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2359 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons46, cons1304, cons1305, CustomConstraint(With2359))
    def replacement2359(u, c, b, B, A, x, y, p, w, v, n, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2359)
        return Dist(q, Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule2359 = ReplacementRule(pattern2359, replacement2359)
    def With2360(u, c, B, A, x, y, p, w, n, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2360 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons4, cons5, cons46, cons1305, CustomConstraint(With2360))
    def replacement2360(u, c, B, A, x, y, p, w, n, n2, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2360)
        return Dist(q, Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule2360 = ReplacementRule(pattern2360, replacement2360)
    def With2361(u, c, b, x, y, p, w, v, n, n2, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern2361 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons1305, CustomConstraint(With2361))
    def replacement2361(u, c, b, x, y, p, w, v, n, n2, m, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(2361)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule2361 = ReplacementRule(pattern2361, replacement2361)
    def With2362(u, c, b, B, z, A, x, y, p, w, v, n, n2, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern2362 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1304, cons1305, CustomConstraint(With2362))
    def replacement2362(u, c, b, B, z, A, x, y, p, w, v, n, n2, m, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(2362)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule2362 = ReplacementRule(pattern2362, replacement2362)
    def With2363(u, c, B, z, A, x, y, p, w, n, n2, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern2363 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1305, CustomConstraint(With2363))
    def replacement2363(u, c, B, z, A, x, y, p, w, n, n2, m, a):
        
        q = Symbol('q')
        r = Symbol('r')
        rubi.append(2363)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule2363 = ReplacementRule(pattern2363, replacement2363)
    def With2364(u, c, d, b, x, y, p, v, n, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2364 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons1304, CustomConstraint(With2364))
    def replacement2364(u, c, d, b, x, y, p, v, n, m, a):
        
        q = DerivativeDivides(y, u, x)
        rubi.append(2364)
        return Dist(q, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), x)
    rule2364 = ReplacementRule(pattern2364, replacement2364)
    def With2365(u, c, d, e, b, x, f, p, w, y, v, n, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern2365 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons1304, cons1305, CustomConstraint(With2365))
    def replacement2365(u, c, d, e, b, x, f, p, w, y, v, n, m, a, q):
        
        r = DerivativeDivides(y, u, x)
        rubi.append(2365)
        return Dist(r, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), x)
    rule2365 = ReplacementRule(pattern2365, replacement2365)
    def With2366(F, v, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2366 = Pattern(Integral(F_**v_*u_, x_), cons1091, cons1091, CustomConstraint(With2366))
    def replacement2366(F, v, u, x):
        
        q = DerivativeDivides(v, u, x)
        rubi.append(2366)
        return Simp(F**v*q/log(F), x)
    rule2366 = ReplacementRule(pattern2366, replacement2366)
    def With2367(F, u, x, w, v, m):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern2367 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), cons1091, cons21, cons1308, CustomConstraint(With2367))
    def replacement2367(F, u, x, w, v, m):
        
        q = DerivativeDivides(v, u, x)
        rubi.append(2367)
        return Dist(q, Subst(Int(F**x*x**m, x), x, v), x)
    rule2367 = ReplacementRule(pattern2367, replacement2367)
    def With2368(u, b, x, p, w, v, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(v*D(w, x) + w*D(v, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2368 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons38, CustomConstraint(With2368))
    def replacement2368(u, b, x, p, w, v, m, a):
        
        c = u/(v*D(w, x) + w*D(v, x))
        rubi.append(2368)
        return Dist(c, Subst(Int((a + b*x**p)**m, x), x, v*w), x)
    rule2368 = ReplacementRule(pattern2368, replacement2368)
    def With2369(u, b, r, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2369 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1309, cons1310, cons1311, CustomConstraint(With2369))
    def replacement2369(u, b, r, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(2369)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w), x)
    rule2369 = ReplacementRule(pattern2369, replacement2369)
    def With2370(u, s, b, r, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2370 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1312, cons1310, cons1311, CustomConstraint(With2370))
    def replacement2370(u, s, b, r, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(2370)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1))), x)
    rule2370 = ReplacementRule(pattern2370, replacement2370)
    def With2371(u, b, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2371 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons1313, cons38, cons17, CustomConstraint(With2371))
    def replacement2371(u, b, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(2371)
        return Dist(c*p, Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))), x)
    rule2371 = ReplacementRule(pattern2371, replacement2371)
    def With2372(u, b, r, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2372 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1314, cons586, cons17, CustomConstraint(With2372))
    def replacement2372(u, b, r, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(2372)
        return -Dist(c*q, Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w), x)
    rule2372 = ReplacementRule(pattern2372, replacement2372)
    def With2373(u, s, b, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2373 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons800, cons1315, cons1316, cons1317, cons17, CustomConstraint(With2373))
    def replacement2373(u, s, b, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(2373)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1))), x)
    rule2373 = ReplacementRule(pattern2373, replacement2373)
    def With2374(u, s, b, r, x, p, w, v, m, a, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern2374 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1318, cons1316, cons1317, cons17, CustomConstraint(With2374))
    def replacement2374(u, s, b, r, x, p, w, v, m, a, q):
        
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(2374)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1))), x)
    rule2374 = ReplacementRule(pattern2374, replacement2374)
    pattern2375 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons21, cons66, cons1319)
    def replacement2375(u, m, x):
        rubi.append(2375)
        return Dist(S(1)/(m + S(1)), Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1))), x)
    rule2375 = ReplacementRule(pattern2375, replacement2375)
    def With2376(u, x):
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
    pattern2376 = Pattern(Integral(u_, x_), CustomConstraint(With2376))
    def replacement2376(u, x):
        
        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(2376)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule2376 = ReplacementRule(pattern2376, replacement2376)
    def With2377(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern2377 = Pattern(Integral(u_, x_), CustomConstraint(With2377))
    def replacement2377(u, x):
        
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        rubi.append(2377)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule2377 = ReplacementRule(pattern2377, replacement2377)
    pattern2378 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons50, cons147, cons10, cons1320, cons1321)
    def replacement2378(u, z, x, p, w, v, n, m, a, q):
        rubi.append(2378)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p), Int(u*v**(m*p)*w**(n*p)*z**(p*q), x), x)
    rule2378 = ReplacementRule(pattern2378, replacement2378)
    pattern2379 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons147, cons10, cons1320)
    def replacement2379(u, x, p, w, v, n, m, a):
        rubi.append(2379)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p), Int(u*v**(m*p)*w**(n*p), x), x)
    rule2379 = ReplacementRule(pattern2379, replacement2379)
    pattern2380 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons5, cons147, cons10, cons1322, cons1323)
    def replacement2380(u, x, p, v, m, a):
        rubi.append(2380)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p), Int(u*v**(m*p), x), x)
    rule2380 = ReplacementRule(pattern2380, replacement2380)
    pattern2381 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons667, cons196, cons1324)
    def replacement2381(u, b, x, p, n, a):
        rubi.append(2381)
        return Dist(FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b)), Int(u*x**(n*p)*(a*x**(-n) + b)**p, x), x)
    rule2381 = ReplacementRule(pattern2381, replacement2381)
    pattern2382 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons147, cons196, cons840, cons1325)
    def replacement2382(u, b, x, p, v, n, a):
        rubi.append(2382)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b)**p, x), x)
    rule2382 = ReplacementRule(pattern2382, replacement2382)
    pattern2383 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons21, cons5, cons147, cons196, cons840)
    def replacement2383(u, b, x, p, v, n, m, a):
        rubi.append(2383)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x), x)
    rule2383 = ReplacementRule(pattern2383, replacement2383)
    def With2384(u, s, b, r, x, m, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        if Not(EqQ(v, S(1))):
            return True
        return False
    pattern2384 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons21, cons52, cons800, cons18, cons1326, CustomConstraint(With2384))
    def replacement2384(u, s, b, r, x, m, a):
        
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        rubi.append(2384)
        return Dist(v, Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), x)
    rule2384 = ReplacementRule(pattern2384, replacement2384)
    def With2385(u, b, x, n, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        if SumQ(v):
            return True
        return False
    pattern2385 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons148, CustomConstraint(With2385))
    def replacement2385(u, b, x, n, a):
        
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        rubi.append(2385)
        return Int(v, x)
    rule2385 = ReplacementRule(pattern2385, replacement2385)
    pattern2386 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38, cons1327)
    def replacement2386(u, c, b, x, p, n, n2, a):
        rubi.append(2386)
        return Dist(S(4)**(-p)*c**(-p), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule2386 = ReplacementRule(pattern2386, replacement2386)
    pattern2387 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons1327)
    def replacement2387(u, c, b, x, p, n, n2, a):
        rubi.append(2387)
        return Dist((b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p, Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule2387 = ReplacementRule(pattern2387, replacement2387)
    def With2388(u, c, b, x, n, n2, a):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        if SumQ(v):
            return True
        return False
    pattern2388 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons46, cons148, CustomConstraint(With2388))
    def replacement2388(u, c, b, x, n, n2, a):
        
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        rubi.append(2388)
        return Int(v, x)
    rule2388 = ReplacementRule(pattern2388, replacement2388)
    pattern2389 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons1328)
    def replacement2389(u, c, b, x, n, m, a):
        rubi.append(2389)
        return Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x)
    rule2389 = ReplacementRule(pattern2389, replacement2389)
    def With2390(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfLinear(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern2390 = Pattern(Integral(u_, x_), CustomConstraint(With2390))
    def replacement2390(u, x):
        
        lst = FunctionOfLinear(u, x)
        rubi.append(2390)
        return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)
    rule2390 = ReplacementRule(pattern2390, replacement2390)
    def With2391(u, x):
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
    pattern2391 = Pattern(Integral(u_/x_, x_), cons1243, cons1324, CustomConstraint(With2391))
    def replacement2391(u, x):
        
        lst = PowerVariableExpn(u, S(0), x)
        rubi.append(2391)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule2391 = ReplacementRule(pattern2391, replacement2391)
    def With2392(u, m, x):
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
    pattern2392 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons17, cons261, cons1243, cons1329, CustomConstraint(With2392))
    def replacement2392(u, m, x):
        
        lst = PowerVariableExpn(u, m + S(1), x)
        rubi.append(2392)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule2392 = ReplacementRule(pattern2392, replacement2392)
    def With2393(u, m, x):
        k = Denominator(m)
        rubi.append(2393)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(S(1)/k)), x)
    pattern2393 = Pattern(Integral(u_*x_**m_, x_), cons367)
    rule2393 = ReplacementRule(pattern2393, With2393)
    def With2394(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        try:
            res = Not(FalseQ(lst))
        except TypeError:
            return False
        if res:
            return True
        return False
    pattern2394 = Pattern(Integral(u_, x_), cons1330, CustomConstraint(With2394))
    def replacement2394(u, x):
        
        lst = FunctionOfSquareRootOfQuadratic(u, x)
        rubi.append(2394)
        return Dist(S(2), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), x)
    rule2394 = ReplacementRule(pattern2394, replacement2394)
    pattern2395 = Pattern(Integral(S(1)/(a_ + v_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement2395(v, b, a, x):
        rubi.append(2395)
        return Dist(S(1)/(S(2)*a), Int(Together(S(1)/(-v/Rt(-a/b, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*a), Int(Together(S(1)/(v/Rt(-a/b, S(2)) + S(1))), x), x)
    rule2395 = ReplacementRule(pattern2395, replacement2395)
    pattern2396 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1331, cons744)
    def replacement2396(b, x, v, n, a):
        rubi.append(2396)
        return Dist(S(2)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x)
    rule2396 = ReplacementRule(pattern2396, replacement2396)
    pattern2397 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1332, cons165)
    def replacement2397(b, x, v, n, a):
        rubi.append(2397)
        return Dist(S(1)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x)
    rule2397 = ReplacementRule(pattern2397, replacement2397)
    pattern2398 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons148, cons1333)
    def replacement2398(u, b, x, v, n, a):
        rubi.append(2398)
        return Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x)
    rule2398 = ReplacementRule(pattern2398, replacement2398)
    def With2399(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = NormalizeIntegrand(u, x)
        if UnsameQ(v, u):
            return True
        return False
    pattern2399 = Pattern(Integral(u_, x_), CustomConstraint(With2399))
    def replacement2399(u, x):
        
        v = NormalizeIntegrand(u, x)
        rubi.append(2399)
        return Int(v, x)
    rule2399 = ReplacementRule(pattern2399, replacement2399)
    def With2400(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = ExpandIntegrand(u, x)
        if SumQ(v):
            return True
        return False
    pattern2400 = Pattern(Integral(u_, x_), CustomConstraint(With2400))
    def replacement2400(u, x):
        
        v = ExpandIntegrand(u, x)
        rubi.append(2400)
        return Int(v, x)
    rule2400 = ReplacementRule(pattern2400, replacement2400)
    pattern2401 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons1334, cons1335, cons1257, cons1336)
    def replacement2401(u, c, b, d, x, p, n, m, a, q):
        rubi.append(2401)
        return Dist(x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q, Int(u*x**(m*p), x), x)
    rule2401 = ReplacementRule(pattern2401, replacement2401)
    pattern2402 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons347)
    def replacement2402(u, c, b, x, p, n, n2, a):
        rubi.append(2402)
        return Dist((S(4)*c)**(-p + S(1)/2)*sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule2402 = ReplacementRule(pattern2402, replacement2402)
    def With2403(u, x):
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
    pattern2403 = Pattern(Integral(u_, x_), CustomConstraint(With2403))
    def replacement2403(u, x):
        
        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(2403)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule2403 = ReplacementRule(pattern2403, replacement2403)
    pattern2404 = Pattern(Integral(u_, x_))
    # def replacement2404(u, x):
    #     rubi.append(2404)
    #     return Int(u, x)
    # rule2404 = ReplacementRule(pattern2404, replacement2404)
    return [rule2337, rule2338, rule2339, rule2340, rule2341, rule2342, rule2343, rule2344, rule2345, rule2346, rule2347, rule2348, rule2349, rule2350, rule2351, rule2352, rule2353, rule2354, rule2355, rule2356, rule2357, rule2358, rule2359, rule2360, rule2361, rule2362, rule2363, rule2364, rule2365, rule2366, rule2367, rule2368, rule2369, rule2370, rule2371, rule2372, rule2373, rule2374, rule2375, rule2376, rule2377, rule2378, rule2379, rule2380, rule2381, rule2382, rule2383, rule2384, rule2385, rule2386, rule2387, rule2388, rule2389, rule2390, rule2391, rule2392, rule2393, rule2394, rule2395, rule2396, rule2397, rule2398, rule2399, rule2400, rule2401, rule2402, rule2403, ]
