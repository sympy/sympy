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
        Zeta, ProductLog, DerivativeDivides, HypergeometricPFQ, IntHide, OneQ, exp, log
    )
    from sympy import (Integral, S, sqrt, And, Or, Integer, Float, Mod, I, Abs, simplify, Mul, Add, Pow)
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (sin, cos, tan, cot, csc, sec, sqrt, erf)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]
    i, ii , Pqq, Q, R, r, C = symbols('i ii Pqq Q R r C_1')
    _UseGamma = False

def miscellaneous_integration(rubi):
    from sympy.integrals.rubi.constraints import cons147, cons1774, cons2, cons3, cons7, cons4, cons5, cons386, cons27, cons50, cons1775, cons1776, cons1777, cons1778, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1779, cons66, cons21, cons84, cons1037, cons1036, cons38, cons1780, cons10, cons1781, cons1782, cons1783, cons209, cons1736, cons1236, cons1784, cons46, cons1785, cons1786, cons1787, cons1788, cons52, cons1789, cons800, cons1790, cons17, cons1791, cons586, cons1792, cons1793, cons1794, cons1795, cons1796, cons1797, cons1798, cons1799, cons1800, cons667, cons196, cons1801, cons840, cons1802, cons18, cons1803, cons148, cons45, cons1804, cons1805, cons1239, cons261, cons1806, cons367, cons1807, cons67, cons1471, cons744, cons1474, cons165, cons1808, cons1809, cons1668, cons1247, cons1810, cons347

    pattern5204 = Pattern(Integral(u_*((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons147, cons1774)
    def replacement5204(a, n, p, b, u, x, c):
        rubi.append(5204)
        return Dist(c**IntPart(p)*(c*(a + b*x)**n)**FracPart(p)*(a + b*x)**(-n*FracPart(p)), Int(u*(a + b*x)**(n*p), x), x)
    rule5204 = ReplacementRule(pattern5204, replacement5204)
    pattern5205 = Pattern(Integral(((d_*(x_*WC('b', S(1)) + WC('a', S(0))))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons5, cons50, cons147, cons386)
    def replacement5205(a, q, p, b, u, x, d, c):
        rubi.append(5205)
        return Dist((c*(d*(a + b*x))**p)**q*(a + b*x)**(-p*q), Int(u*(a + b*x)**(p*q), x), x)
    rule5205 = ReplacementRule(pattern5205, replacement5205)
    pattern5206 = Pattern(Integral((((x_*WC('b', S(1)) + WC('a', S(0)))**n_*WC('d', S(1)))**p_*WC('c', S(1)))**q_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons5, cons50, cons147, cons386)
    def replacement5206(a, n, q, p, b, u, x, d, c):
        rubi.append(5206)
        return Dist((c*(d*(a + b*x)**n)**p)**q*(a + b*x)**(-n*p*q), Int(u*(a + b*x)**(n*p*q), x), x)
    rule5206 = ReplacementRule(pattern5206, replacement5206)
    pattern5207 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1775, cons1776, cons1777, cons1778)
    def replacement5207(e, a, n, B, g, b, F, A, f, x, d, C, c):
        rubi.append(5207)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule5207 = ReplacementRule(pattern5207, replacement5207)
    pattern5208 = Pattern(Integral((F_*sqrt(x_*WC('e', S(1)) + S(1))*WC('b', S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1091, cons1775, cons1779)
    def replacement5208(e, a, n, g, b, F, A, x, C, c):
        rubi.append(5208)
        return Dist(g/C, Subst(Int((a + b*F(c*x))**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule5208 = ReplacementRule(pattern5208, replacement5208)
    pattern5209 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + WC('f', S(0))))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons34, cons35, cons36, cons1091, cons1775, cons1776, cons1777, cons1778)
    def replacement5209(e, a, n, B, g, b, F, A, f, x, d, C, c):
        rubi.append(5209)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(d + e*x)/sqrt(f + g*x)), x)
    rule5209 = ReplacementRule(pattern5209, replacement5209)
    pattern5210 = Pattern(Integral((F_**(sqrt(x_*WC('e', S(1)) + S(1))*WC('c', S(1))/sqrt(x_*WC('g', S(1)) + S(1)))*WC('b', S(1)) + WC('a', S(0)))**WC('n', S(1))/(x_**S(2)*WC('C', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons7, cons48, cons208, cons34, cons36, cons1091, cons1775, cons1779)
    def replacement5210(e, a, n, g, b, F, A, x, C, c):
        rubi.append(5210)
        return Dist(g/C, Subst(Int((F**(c*x)*b + a)**n/x, x), x, sqrt(e*x + S(1))/sqrt(g*x + S(1))), x)
    rule5210 = ReplacementRule(pattern5210, replacement5210)
    def With5211(u, y, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5211 = Pattern(Integral(u_/y_, x_), CustomConstraint(With5211))
    def replacement5211(u, y, x):

        q = DerivativeDivides(y, u, x)
        rubi.append(5211)
        return Simp(q*log(RemoveContent(y, x)), x)
    rule5211 = ReplacementRule(pattern5211, replacement5211)
    def With5212(w, u, y, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(w*y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5212 = Pattern(Integral(u_/(w_*y_), x_), CustomConstraint(With5212))
    def replacement5212(w, u, y, x):

        q = DerivativeDivides(w*y, u, x)
        rubi.append(5212)
        return Simp(q*log(RemoveContent(w*y, x)), x)
    rule5212 = ReplacementRule(pattern5212, replacement5212)
    def With5213(x, u, y, m):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5213 = Pattern(Integral(u_*y_**WC('m', S(1)), x_), cons21, cons66, CustomConstraint(With5213))
    def replacement5213(x, u, y, m):

        q = DerivativeDivides(y, u, x)
        rubi.append(5213)
        return Simp(q*y**(m + S(1))/(m + S(1)), x)
    rule5213 = ReplacementRule(pattern5213, replacement5213)
    def With5214(n, y, m, z, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5214 = Pattern(Integral(u_*y_**WC('m', S(1))*z_**WC('n', S(1)), x_), cons21, cons4, cons66, CustomConstraint(With5214))
    def replacement5214(n, y, m, z, u, x):

        q = DerivativeDivides(y*z, u*z**(-m + n), x)
        rubi.append(5214)
        return Simp(q*y**(m + S(1))*z**(m + S(1))/(m + S(1)), x)
    rule5214 = ReplacementRule(pattern5214, replacement5214)
    def With5215(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = SimplifyIntegrand(u, x)
        if SimplerIntegrandQ(v, u, x):
            return True
        return False
    pattern5215 = Pattern(Integral(u_, x_), CustomConstraint(With5215))
    def replacement5215(u, x):

        v = SimplifyIntegrand(u, x)
        rubi.append(5215)
        return Int(v, x)
    rule5215 = ReplacementRule(pattern5215, replacement5215)
    pattern5216 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1037)
    def replacement5216(e, a, n, m, f, b, u, x, d, c):
        rubi.append(5216)
        return Dist((a*e**S(2) - c*f**S(2))**m, Int(ExpandIntegrand(u*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule5216 = ReplacementRule(pattern5216, replacement5216)
    pattern5217 = Pattern(Integral((sqrt(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))*WC('e', S(1)) + sqrt(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))*WC('f', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons84, cons1036)
    def replacement5217(e, a, n, m, f, b, u, x, d, c):
        rubi.append(5217)
        return Dist((b*e**S(2) - d*f**S(2))**m, Int(ExpandIntegrand(u*x**(m*n)*(e*sqrt(a + b*x**n) - f*sqrt(c + d*x**n))**(-m), x), x), x)
    rule5217 = ReplacementRule(pattern5217, replacement5217)
    pattern5218 = Pattern(Integral(u_**WC('m', S(1))*w_*(u_**n_*WC('a', S(1)) + v_)**WC('p', S(1)), x_), cons2, cons21, cons4, cons38, cons1780, cons10)
    def replacement5218(w, a, n, m, v, p, u, x):
        rubi.append(5218)
        return Int(u**(m + n*p)*w*(a + u**(-n)*v)**p, x)
    rule5218 = ReplacementRule(pattern5218, replacement5218)
    def With5219(n, a, y, m, v, b, u, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5219 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons1781, CustomConstraint(With5219))
    def replacement5219(n, a, y, m, v, b, u, x, d, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5219)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, y), x)
    rule5219 = ReplacementRule(pattern5219, replacement5219)
    def With5220(w, e, n, a, u, y, m, v, p, b, f, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5220 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons1781, cons1782, CustomConstraint(With5220))
    def replacement5220(w, e, n, a, u, y, m, v, p, b, f, x, d, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5220)
        return Dist(q, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, y), x)
    rule5220 = ReplacementRule(pattern5220, replacement5220)
    def With5221(w, e, n, a, h, q, u, y, m, z, v, g, p, b, f, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern5221 = Pattern(Integral(u_*(v_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(w_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(y_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(z_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons1781, cons1782, cons1783, CustomConstraint(With5221))
    def replacement5221(w, e, n, a, h, q, u, y, m, z, v, g, p, b, f, x, d, c):

        r = DerivativeDivides(y, u, x)
        rubi.append(5221)
        return Dist(r, Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, y), x)
    rule5221 = ReplacementRule(pattern5221, replacement5221)
    def With5222(a, n, y, b, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5222 = Pattern(Integral((a_ + y_**n_*WC('b', S(1)))*WC('u', S(1)), x_), cons2, cons3, cons4, cons1736, CustomConstraint(With5222))
    def replacement5222(a, n, y, b, u, x):

        q = DerivativeDivides(y, u, x)
        rubi.append(5222)
        return Dist(a, Int(u, x), x) + Dist(b*q, Subst(Int(x**n, x), x, y), x)
    rule5222 = ReplacementRule(pattern5222, replacement5222)
    def With5223(a, n, y, p, b, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5223 = Pattern(Integral((y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons4, cons5, cons1236, CustomConstraint(With5223))
    def replacement5223(a, n, y, p, b, u, x):

        q = DerivativeDivides(y, u, x)
        rubi.append(5223)
        return Dist(q, Subst(Int((a + b*x**n)**p, x), x, y), x)
    rule5223 = ReplacementRule(pattern5223, replacement5223)
    def With5224(a, n, y, m, v, p, b, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern5224 = Pattern(Integral(v_**WC('m', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons21, cons4, cons5, cons1784, CustomConstraint(With5224))
    def replacement5224(a, n, y, m, v, p, b, u, x):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(5224)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n)**p, x), x, y), x)
    rule5224 = ReplacementRule(pattern5224, replacement5224)
    def With5225(a, n, y, n2, v, p, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5225 = Pattern(Integral((v_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons7, cons4, cons5, cons46, cons1781, CustomConstraint(With5225))
    def replacement5225(a, n, y, n2, v, p, b, u, x, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5225)
        return Dist(q, Subst(Int((a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule5225 = ReplacementRule(pattern5225, replacement5225)
    def With5226(w, a, n, B, y, n2, v, p, b, A, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5226 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons4, cons5, cons46, cons1781, cons1782, CustomConstraint(With5226))
    def replacement5226(w, a, n, B, y, n2, v, p, b, A, u, x, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5226)
        return Dist(q, Subst(Int((A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule5226 = ReplacementRule(pattern5226, replacement5226)
    def With5227(w, a, n, B, y, n2, p, A, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5227 = Pattern(Integral((A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons4, cons5, cons46, cons1782, CustomConstraint(With5227))
    def replacement5227(w, a, n, B, y, n2, p, A, u, x, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5227)
        return Dist(q, Subst(Int((A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule5227 = ReplacementRule(pattern5227, replacement5227)
    def With5228(w, a, n, y, m, n2, v, p, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        if And(Not(FalseQ(Set(r, Divides(y**m, v**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern5228 = Pattern(Integral(v_**WC('m', S(1))*(w_**WC('n2', S(1))*WC('c', S(1)) + y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons21, cons4, cons5, cons46, cons1782, CustomConstraint(With5228))
    def replacement5228(w, a, n, y, m, n2, v, p, b, u, x, c):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, v**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(5228)
        return Dist(q*r, Subst(Int(x**m*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule5228 = ReplacementRule(pattern5228, replacement5228)
    def With5229(w, a, n, B, m, n2, v, y, z, p, b, A, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern5229 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(v_**n_*WC('b', S(1)) + w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1781, cons1782, CustomConstraint(With5229))
    def replacement5229(w, a, n, B, m, n2, v, y, z, p, b, A, u, x, c):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(5229)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + b*x**n + c*x**(S(2)*n))**p, x), x, y), x)
    rule5229 = ReplacementRule(pattern5229, replacement5229)
    def With5230(w, a, n, B, m, n2, y, z, p, A, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        if And(Not(FalseQ(Set(r, Divides(y**m, z**m, x)))), Not(FalseQ(Set(q, DerivativeDivides(y, u, x))))):
            return True
        return False
    pattern5230 = Pattern(Integral(z_**WC('m', S(1))*(A_ + y_**n_*WC('B', S(1)))*(w_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*WC('u', S(1)), x_), cons2, cons7, cons34, cons35, cons21, cons4, cons5, cons46, cons1782, CustomConstraint(With5230))
    def replacement5230(w, a, n, B, m, n2, y, z, p, A, u, x, c):

        q = Symbol('q')
        r = Symbol('r')
        r = Divides(y**m, z**m, x)
        q = DerivativeDivides(y, u, x)
        rubi.append(5230)
        return Dist(q*r, Subst(Int(x**m*(A + B*x**n)*(a + c*x**(S(2)*n))**p, x), x, y), x)
    rule5230 = ReplacementRule(pattern5230, replacement5230)
    def With5231(a, n, y, m, v, p, b, u, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(y, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5231 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons1781, CustomConstraint(With5231))
    def replacement5231(a, n, y, m, v, p, b, u, x, d, c):

        q = DerivativeDivides(y, u, x)
        rubi.append(5231)
        return Dist(q, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p, x), x, y), x)
    rule5231 = ReplacementRule(pattern5231, replacement5231)
    def With5232(w, e, a, n, q, y, m, v, f, p, b, u, x, d, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        r = DerivativeDivides(y, u, x)
        if Not(FalseQ(r)):
            return True
        return False
    pattern5232 = Pattern(Integral((v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('p', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('q', S(1))*(y_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons50, cons1781, cons1782, CustomConstraint(With5232))
    def replacement5232(w, e, a, n, q, y, m, v, f, p, b, u, x, d, c):

        r = DerivativeDivides(y, u, x)
        rubi.append(5232)
        return Dist(r, Subst(Int((a + b*x**n)**m*(c + d*x**n)**p*(e + f*x**n)**q, x), x, y), x)
    rule5232 = ReplacementRule(pattern5232, replacement5232)
    def With5233(v, x, u, F):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5233 = Pattern(Integral(F_**v_*u_, x_), cons1091, cons1091, CustomConstraint(With5233))
    def replacement5233(v, x, u, F):

        q = DerivativeDivides(v, u, x)
        rubi.append(5233)
        return Simp(F**v*q/log(F), x)
    rule5233 = ReplacementRule(pattern5233, replacement5233)
    def With5234(u, m, v, F, w, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        q = DerivativeDivides(v, u, x)
        if Not(FalseQ(q)):
            return True
        return False
    pattern5234 = Pattern(Integral(F_**v_*u_*w_**WC('m', S(1)), x_), cons1091, cons21, cons1785, CustomConstraint(With5234))
    def replacement5234(u, m, v, F, w, x):

        q = DerivativeDivides(v, u, x)
        rubi.append(5234)
        return Dist(q, Subst(Int(F**x*x**m, x), x, v), x)
    rule5234 = ReplacementRule(pattern5234, replacement5234)
    def With5235(w, a, m, v, p, b, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(v*D(w, x) + w*D(v, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5235 = Pattern(Integral(u_*(a_ + v_**WC('p', S(1))*w_**WC('p', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons38, CustomConstraint(With5235))
    def replacement5235(w, a, m, v, p, b, u, x):

        c = u/(v*D(w, x) + w*D(v, x))
        rubi.append(5235)
        return Dist(c, Subst(Int((a + b*x**p)**m, x), x, v*w), x)
    rule5235 = ReplacementRule(pattern5235, replacement5235)
    def With5236(w, a, m, v, p, b, r, u, x, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5236 = Pattern(Integral(u_*v_**WC('r', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1786, cons1787, cons1788, CustomConstraint(With5236))
    def replacement5236(w, a, m, v, p, b, r, u, x, q):

        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(5236)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w), x)
    rule5236 = ReplacementRule(pattern5236, replacement5236)
    def With5237(w, a, m, v, p, b, r, u, x, q, s):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) + q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5237 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(a_ + v_**WC('p', S(1))*w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1789, cons1787, cons1788, CustomConstraint(With5237))
    def replacement5237(w, a, m, v, p, b, r, u, x, q, s):

        c = u/(p*w*D(v, x) + q*v*D(w, x))
        rubi.append(5237)
        return Dist(c*p/(r + S(1)), Subst(Int((a + b*x**(p/(r + S(1))))**m, x), x, v**(r + S(1))*w**(s + S(1))), x)
    rule5237 = ReplacementRule(pattern5237, replacement5237)
    def With5238(w, a, m, v, p, b, u, x, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5238 = Pattern(Integral(u_*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons1790, cons38, cons17, CustomConstraint(With5238))
    def replacement5238(w, a, m, v, p, b, u, x, q):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(5238)
        return Dist(c*p, Subst(Int((a*x**p + b)**m, x), x, v*w**(m*q + S(1))), x)
    rule5238 = ReplacementRule(pattern5238, replacement5238)
    def With5239(w, a, m, v, p, b, r, u, x, q):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5239 = Pattern(Integral(u_*v_**WC('r', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons1791, cons586, cons17, CustomConstraint(With5239))
    def replacement5239(w, a, m, v, p, b, r, u, x, q):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(5239)
        return -Dist(c*q, Subst(Int((a + b*x**q)**m, x), x, v**(m*p + r + S(1))*w), x)
    rule5239 = ReplacementRule(pattern5239, replacement5239)
    def With5240(w, a, m, v, p, b, u, x, q, s):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5240 = Pattern(Integral(u_*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons800, cons1792, cons1793, cons1794, cons17, CustomConstraint(With5240))
    def replacement5240(w, a, m, v, p, b, u, x, q, s):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(5240)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + S(1))*w**(s + S(1))), x)
    rule5240 = ReplacementRule(pattern5240, replacement5240)
    def With5241(w, a, m, v, p, b, r, u, x, q, s):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        c = u/(p*w*D(v, x) - q*v*D(w, x))
        if FreeQ(c, x):
            return True
        return False
    pattern5241 = Pattern(Integral(u_*v_**WC('r', S(1))*w_**WC('s', S(1))*(v_**WC('p', S(1))*WC('a', S(1)) + w_**WC('q', S(1))*WC('b', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons21, cons5, cons50, cons52, cons800, cons1795, cons1793, cons1794, cons17, CustomConstraint(With5241))
    def replacement5241(w, a, m, v, p, b, r, u, x, q, s):

        c = u/(p*w*D(v, x) - q*v*D(w, x))
        rubi.append(5241)
        return -Dist(c*q/(s + S(1)), Subst(Int((a + b*x**(q/(s + S(1))))**m, x), x, v**(m*p + r + S(1))*w**(s + S(1))), x)
    rule5241 = ReplacementRule(pattern5241, replacement5241)
    pattern5242 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons21, cons66, cons1796)
    def replacement5242(x, u, m):
        rubi.append(5242)
        return Dist(S(1)/(m + S(1)), Subst(Int(SubstFor(x**(m + S(1)), u, x), x), x, x**(m + S(1))), x)
    rule5242 = ReplacementRule(pattern5242, replacement5242)
    def With5243(u, x):
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
    pattern5243 = Pattern(Integral(u_, x_), CustomConstraint(With5243))
    def replacement5243(u, x):

        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(5243)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule5243 = ReplacementRule(pattern5243, replacement5243)
    def With5244(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern5244 = Pattern(Integral(u_, x_), CustomConstraint(With5244))
    def replacement5244(u, x):

        lst = SubstForFractionalPowerOfQuotientOfLinears(u, x)
        rubi.append(5244)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule5244 = ReplacementRule(pattern5244, replacement5244)
    pattern5245 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*z_**WC('q', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons50, cons147, cons10, cons1797, cons1798)
    def replacement5245(w, a, n, z, m, v, p, u, x, q):
        rubi.append(5245)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*z**(-q*FracPart(p))*(a*v**m*w**n*z**q)**FracPart(p), Int(u*v**(m*p)*w**(n*p)*z**(p*q), x), x)
    rule5245 = ReplacementRule(pattern5245, replacement5245)
    pattern5246 = Pattern(Integral((v_**WC('m', S(1))*w_**WC('n', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons4, cons5, cons147, cons10, cons1797)
    def replacement5246(w, a, n, m, v, p, u, x):
        rubi.append(5246)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*w**(-n*FracPart(p))*(a*v**m*w**n)**FracPart(p), Int(u*v**(m*p)*w**(n*p), x), x)
    rule5246 = ReplacementRule(pattern5246, replacement5246)
    pattern5247 = Pattern(Integral((v_**WC('m', S(1))*WC('a', S(1)))**p_*WC('u', S(1)), x_), cons2, cons21, cons5, cons147, cons10, cons1799, cons1800)
    def replacement5247(a, m, v, p, u, x):
        rubi.append(5247)
        return Dist(a**IntPart(p)*v**(-m*FracPart(p))*(a*v**m)**FracPart(p), Int(u*v**(m*p), x), x)
    rule5247 = ReplacementRule(pattern5247, replacement5247)
    pattern5248 = Pattern(Integral((x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons667, cons196, cons1801)
    def replacement5248(a, n, p, b, u, x):
        rubi.append(5248)
        return Dist(FullSimplify(x**(-n/S(2))*sqrt(a + b*x**n)/sqrt(a*x**(-n) + b)), Int(u*x**(n*p)*(a*x**(-n) + b)**p, x), x)
    rule5248 = ReplacementRule(pattern5248, replacement5248)
    pattern5249 = Pattern(Integral((v_**n_*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons5, cons147, cons196, cons840, cons1802)
    def replacement5249(a, n, v, p, b, u, x):
        rubi.append(5249)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n)**FracPart(p)*(a*v**(-n) + b)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b)**p, x), x)
    rule5249 = ReplacementRule(pattern5249, replacement5249)
    pattern5250 = Pattern(Integral((v_**n_*x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**p_*WC('u', S(1)), x_), cons2, cons3, cons21, cons5, cons147, cons196, cons840)
    def replacement5250(a, n, m, v, p, b, u, x):
        rubi.append(5250)
        return Dist(v**(-n*FracPart(p))*(a + b*v**n*x**m)**FracPart(p)*(a*v**(-n) + b*x**m)**(-FracPart(p)), Int(u*v**(n*p)*(a*v**(-n) + b*x**m)**p, x), x)
    rule5250 = ReplacementRule(pattern5250, replacement5250)
    def With5251(a, m, b, r, u, x, s):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        if Not(EqQ(v, S(1))):
            return True
        return False
    pattern5251 = Pattern(Integral((x_**WC('r', S(1))*WC('a', S(1)) + x_**WC('s', S(1))*WC('b', S(1)))**m_*WC('u', S(1)), x_), cons2, cons3, cons21, cons52, cons800, cons18, cons1803, CustomConstraint(With5251))
    def replacement5251(a, m, b, r, u, x, s):

        v = x**(-r*FracPart(m))*(a + b*x**(-r + s))**(-FracPart(m))*(a*x**r + b*x**s)**FracPart(m)
        rubi.append(5251)
        return Dist(v, Int(u*x**(m*r)*(a + b*x**(-r + s))**m, x), x)
    rule5251 = ReplacementRule(pattern5251, replacement5251)
    def With5252(a, n, b, u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n), x)
        if SumQ(v):
            return True
        return False
    pattern5252 = Pattern(Integral(u_/(a_ + x_**n_*WC('b', S(1))), x_), cons2, cons3, cons148, CustomConstraint(With5252))
    def replacement5252(a, n, b, u, x):

        v = RationalFunctionExpand(u/(a + b*x**n), x)
        rubi.append(5252)
        return Int(v, x)
    rule5252 = ReplacementRule(pattern5252, replacement5252)
    pattern5253 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons4, cons46, cons45, cons38, cons1804)
    def replacement5253(a, n, n2, p, b, u, x, c):
        rubi.append(5253)
        return Dist(S(4)**(-p)*c**(-p), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule5253 = ReplacementRule(pattern5253, replacement5253)
    pattern5254 = Pattern(Integral(u_*(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons147, cons1804)
    def replacement5254(a, n, n2, p, b, u, x, c):
        rubi.append(5254)
        return Dist((b + S(2)*c*x**n)**(-S(2)*p)*(a + b*x**n + c*x**(S(2)*n))**p, Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule5254 = ReplacementRule(pattern5254, replacement5254)
    def With5255(a, n, n2, b, u, x, c):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        if SumQ(v):
            return True
        return False
    pattern5255 = Pattern(Integral(u_/(x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons46, cons148, CustomConstraint(With5255))
    def replacement5255(a, n, n2, b, u, x, c):

        v = RationalFunctionExpand(u/(a + b*x**n + c*x**(S(2)*n)), x)
        rubi.append(5255)
        return Int(v, x)
    rule5255 = ReplacementRule(pattern5255, replacement5255)
    pattern5256 = Pattern(Integral(WC('u', S(1))/(x_**WC('m', S(1))*WC('a', S(1)) + sqrt(x_**n_*WC('c', S(1)))*WC('b', S(1))), x_), cons2, cons3, cons7, cons21, cons4, cons1805)
    def replacement5256(a, n, m, b, u, x, c):
        rubi.append(5256)
        return Int(u*(a*x**m - b*sqrt(c*x**n))/(a**S(2)*x**(S(2)*m) - b**S(2)*c*x**n), x)
    rule5256 = ReplacementRule(pattern5256, replacement5256)
    def With5257(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        lst = FunctionOfLinear(u, x)
        if Not(FalseQ(lst)):
            return True
        return False
    pattern5257 = Pattern(Integral(u_, x_), CustomConstraint(With5257))
    def replacement5257(u, x):

        lst = FunctionOfLinear(u, x)
        rubi.append(5257)
        return Dist(S(1)/Part(lst, S(3)), Subst(Int(Part(lst, S(1)), x), x, x*Part(lst, S(3)) + Part(lst, S(2))), x)
    rule5257 = ReplacementRule(pattern5257, replacement5257)
    def With5258(u, x):
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
    pattern5258 = Pattern(Integral(u_/x_, x_), cons1239, cons1801, CustomConstraint(With5258))
    def replacement5258(u, x):

        lst = PowerVariableExpn(u, S(0), x)
        rubi.append(5258)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule5258 = ReplacementRule(pattern5258, replacement5258)
    def With5259(x, u, m):
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
    pattern5259 = Pattern(Integral(u_*x_**WC('m', S(1)), x_), cons17, cons261, cons1239, cons1806, CustomConstraint(With5259))
    def replacement5259(x, u, m):

        lst = PowerVariableExpn(u, m + S(1), x)
        rubi.append(5259)
        return Dist(S(1)/Part(lst, S(2)), Subst(Int(NormalizeIntegrand(Part(lst, S(1))/x, x), x), x, (x*Part(lst, S(3)))**Part(lst, S(2))), x)
    rule5259 = ReplacementRule(pattern5259, replacement5259)
    def With5260(u, x, m):
        k = Denominator(m)
        rubi.append(5260)
        return Dist(k, Subst(Int(x**(k*(m + S(1)) + S(-1))*ReplaceAll(u, Rule(x, x**k)), x), x, x**(S(1)/k)), x)
    pattern5260 = Pattern(Integral(u_*x_**m_, x_), cons367)
    rule5260 = ReplacementRule(pattern5260, With5260)
    def With5261(u, x):
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
    pattern5261 = Pattern(Integral(u_, x_), cons1807, CustomConstraint(With5261))
    def replacement5261(u, x):

        lst = FunctionOfSquareRootOfQuadratic(u, x)
        rubi.append(5261)
        return Dist(S(2), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(2))), x)
    rule5261 = ReplacementRule(pattern5261, replacement5261)
    pattern5262 = Pattern(Integral(S(1)/(a_ + v_**S(2)*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement5262(v, a, x, b):
        rubi.append(5262)
        return Dist(S(1)/(S(2)*a), Int(Together(S(1)/(-v/Rt(-a/b, S(2)) + S(1))), x), x) + Dist(S(1)/(S(2)*a), Int(Together(S(1)/(v/Rt(-a/b, S(2)) + S(1))), x), x)
    rule5262 = ReplacementRule(pattern5262, replacement5262)
    pattern5263 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1471, cons744)
    def replacement5263(a, n, v, b, x):
        rubi.append(5263)
        return Dist(S(2)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(4)*k/n)*v**S(2)/Rt(-a/b, n/S(2)))), x), List(k, S(1), n/S(2))), x)
    rule5263 = ReplacementRule(pattern5263, replacement5263)
    pattern5264 = Pattern(Integral(S(1)/(a_ + v_**n_*WC('b', S(1))), x_), cons2, cons3, cons1474, cons165)
    def replacement5264(a, n, v, b, x):
        rubi.append(5264)
        return Dist(S(1)/(a*n), Sum_doit(Int(Together(S(1)/(S(1) - (S(-1))**(-S(2)*k/n)*v/Rt(-a/b, n))), x), List(k, S(1), n)), x)
    rule5264 = ReplacementRule(pattern5264, replacement5264)
    pattern5265 = Pattern(Integral(v_/(a_ + u_**WC('n', S(1))*WC('b', S(1))), x_), cons2, cons3, cons148, cons1808)
    def replacement5265(n, a, v, b, u, x):
        rubi.append(5265)
        return Int(ReplaceAll(ExpandIntegrand(PolynomialInSubst(v, u, x)/(a + b*x**n), x), Rule(x, u)), x)
    rule5265 = ReplacementRule(pattern5265, replacement5265)
    def With5266(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = NormalizeIntegrand(u, x)
        if UnsameQ(v, u):
            return True
        return False
    pattern5266 = Pattern(Integral(u_, x_), CustomConstraint(With5266))
    def replacement5266(u, x):

        v = NormalizeIntegrand(u, x)
        rubi.append(5266)
        return Int(v, x)
    rule5266 = ReplacementRule(pattern5266, replacement5266)
    def With5267(u, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        v = ExpandIntegrand(u, x)
        if SumQ(v):
            return True
        return False
    pattern5267 = Pattern(Integral(u_, x_), CustomConstraint(With5267))
    def replacement5267(u, x):

        v = ExpandIntegrand(u, x)
        rubi.append(5267)
        return Int(v, x)
    rule5267 = ReplacementRule(pattern5267, replacement5267)
    pattern5268 = Pattern(Integral((x_**WC('m', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**WC('n', S(1))*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*WC('u', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons5, cons50, cons1809, cons1668, cons1247, cons1810)
    def replacement5268(a, n, q, m, p, b, u, x, d, c):
        rubi.append(5268)
        return Dist(x**(-m*p)*(a + b*x**m)**p*(c + d*x**n)**q, Int(u*x**(m*p), x), x)
    rule5268 = ReplacementRule(pattern5268, replacement5268)
    pattern5269 = Pattern(Integral(u_*(a_ + x_**WC('n', S(1))*WC('b', S(1)) + x_**WC('n2', S(1))*WC('c', S(1)))**p_, x_), cons2, cons3, cons7, cons4, cons5, cons46, cons45, cons347)
    def replacement5269(n, a, n2, p, b, u, x, c):
        rubi.append(5269)
        return Dist((S(4)*c)**(-p + S(1)/2)*sqrt(a + b*x**n + c*x**(S(2)*n))/(b + S(2)*c*x**n), Int(u*(b + S(2)*c*x**n)**(S(2)*p), x), x)
    rule5269 = ReplacementRule(pattern5269, replacement5269)
    def With5270(u, x):
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
    pattern5270 = Pattern(Integral(u_, x_), CustomConstraint(With5270))
    def replacement5270(u, x):

        lst = SubstForFractionalPowerOfLinear(u, x)
        rubi.append(5270)
        return Dist(Part(lst, S(2))*Part(lst, S(4)), Subst(Int(Part(lst, S(1)), x), x, Part(lst, S(3))**(S(1)/Part(lst, S(2)))), x)
    rule5270 = ReplacementRule(pattern5270, replacement5270)
    pattern5271 = Pattern(Integral(u_, x_))
    # def replacement5271(u, x):
    #     rubi.append(5271)
    #     return Int(u, x)
    # rule5271 = ReplacementRule(pattern5271, replacement5271)
    return [rule5204, rule5205, rule5206, rule5207, rule5208, rule5209, rule5210, rule5211, rule5212, rule5213, rule5214, rule5215, rule5216, rule5217, rule5218, rule5219, rule5220, rule5221, rule5222, rule5223, rule5224, rule5225, rule5226, rule5227, rule5228, rule5229, rule5230, rule5231, rule5232, rule5233, rule5234, rule5235, rule5236, rule5237, rule5238, rule5239, rule5240, rule5241, rule5242, rule5243, rule5244, rule5245, rule5246, rule5247, rule5248, rule5249, rule5250, rule5251, rule5252, rule5253, rule5254, rule5255, rule5256, rule5257, rule5258, rule5259, rule5260, rule5261, rule5262, rule5263, rule5264, rule5265, rule5266, rule5267, rule5268, rule5269, rule5270, ]
