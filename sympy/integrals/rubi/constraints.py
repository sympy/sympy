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


    def cons_f1(a):
        return ZeroQ(a)

    cons1 = CustomConstraint(lambda a: cons_f1(a))

    def cons_f2(a, x):
        return FreeQ(a, x)

    cons2 = CustomConstraint(lambda a, x: cons_f2(a, x))

    def cons_f3(b, x):
        return FreeQ(b, x)

    cons3 = CustomConstraint(lambda b, x: cons_f3(b, x))

    def cons_f4(n, x):
        return FreeQ(n, x)

    cons4 = CustomConstraint(lambda n, x: cons_f4(n, x))

    def cons_f5(p, x):
        return FreeQ(p, x)

    cons5 = CustomConstraint(lambda p, x: cons_f5(p, x))

    def cons_f6(n, j):
        return ZeroQ(j - S(2)*n)

    cons6 = CustomConstraint(lambda n, j: cons_f6(n, j))

    def cons_f7(c, x):
        return FreeQ(c, x)

    cons7 = CustomConstraint(lambda c, x: cons_f7(c, x))

    def cons_f8(b):
        return ZeroQ(b)

    cons8 = CustomConstraint(lambda b: cons_f8(b))

    def cons_f9(c):
        return ZeroQ(c)

    cons9 = CustomConstraint(lambda c: cons_f9(c))

    def cons_f10(v, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return NFreeQ(v, x)

    cons10 = CustomConstraint(lambda v, x: cons_f10(v, x))

    def cons_f11(x, Pm):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return PolyQ(Pm, x)

    cons11 = CustomConstraint(lambda x, Pm: cons_f11(x, Pm))

    def cons_f12(p):
        return Not(RationalQ(p))

    cons12 = CustomConstraint(lambda p: cons_f12(p))

    def cons_f13(p):
        return RationalQ(p)

    cons13 = CustomConstraint(lambda p: cons_f13(p))

    def cons_f14(x, c, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c), x)

    cons14 = CustomConstraint(lambda x, c, a, b: cons_f14(x, c, a, b))

    def cons_f15(a):
        return EqQ(a**S(2), S(1))

    cons15 = CustomConstraint(lambda a: cons_f15(a))

    def cons_f16(u):
        return SumQ(u)

    cons16 = CustomConstraint(lambda u: cons_f16(u))

    def cons_f17(m):
        return IntegerQ(m)

    cons17 = CustomConstraint(lambda m: cons_f17(m))

    def cons_f18(m):
        return Not(IntegerQ(m))

    cons18 = CustomConstraint(lambda m: cons_f18(m))

    def cons_f19(n):
        return PositiveIntegerQ(n + S(1)/2)

    cons19 = CustomConstraint(lambda n: cons_f19(n))

    def cons_f20(n, m):
        return IntegerQ(m + n)

    cons20 = CustomConstraint(lambda n, m: cons_f20(n, m))

    def cons_f21(m, x):
        return FreeQ(m, x)

    cons21 = CustomConstraint(lambda m, x: cons_f21(m, x))

    def cons_f22(n):
        return NegativeIntegerQ(n + S(-1)/2)

    cons22 = CustomConstraint(lambda n: cons_f22(n))

    def cons_f23(n):
        return Not(IntegerQ(n))

    cons23 = CustomConstraint(lambda n: cons_f23(n))

    def cons_f24(n, m):
        return Not(IntegerQ(m + n))

    cons24 = CustomConstraint(lambda n, m: cons_f24(n, m))

    def cons_f25(d, c, a, b):
        return ZeroQ(-a*d + b*c)

    cons25 = CustomConstraint(lambda d, c, a, b: cons_f25(d, c, a, b))

    def cons_f26(c, a, x, b, n, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Or(Not(IntegerQ(n)), SimplerQ(c + d*x, a + b*x))

    cons26 = CustomConstraint(lambda c, a, x, b, n, d: cons_f26(c, a, x, b, n, d))

    def cons_f27(d, x):
        return FreeQ(d, x)

    cons27 = CustomConstraint(lambda d, x: cons_f27(d, x))

    def cons_f28(d, b):
        return PositiveQ(b/d)

    cons28 = CustomConstraint(lambda d, b: cons_f28(d, b))

    def cons_f29(n, m):
        return Not(Or(IntegerQ(m), IntegerQ(n)))

    cons29 = CustomConstraint(lambda n, m: cons_f29(n, m))

    def cons_f30(n, d, m, b):
        return Not(Or(IntegerQ(m), IntegerQ(n), PositiveQ(b/d)))

    cons30 = CustomConstraint(lambda n, d, m, b: cons_f30(n, d, m, b))

    def cons_f31(m):
        return RationalQ(m)

    cons31 = CustomConstraint(lambda m: cons_f31(m))

    def cons_f32(m):
        return LessEqual(m, S(-1))

    cons32 = CustomConstraint(lambda m: cons_f32(m))

    def cons_f33(C, B, a, A, b):
        return ZeroQ(A*b**S(2) - B*a*b + C*a**S(2))

    cons33 = CustomConstraint(lambda C, B, a, A, b: cons_f33(C, B, a, A, b))

    def cons_f34(A, x):
        return FreeQ(A, x)

    cons34 = CustomConstraint(lambda A, x: cons_f34(A, x))

    def cons_f35(B, x):
        return FreeQ(B, x)

    cons35 = CustomConstraint(lambda B, x: cons_f35(B, x))

    def cons_f36(C, x):
        return FreeQ(C, x)

    cons36 = CustomConstraint(lambda C, x: cons_f36(C, x))

    def cons_f37(n, q):
        return ZeroQ(n + q)

    cons37 = CustomConstraint(lambda n, q: cons_f37(n, q))

    def cons_f38(p):
        return IntegerQ(p)

    cons38 = CustomConstraint(lambda p: cons_f38(p))

    def cons_f39(d, c, a, b):
        return ZeroQ(a*c - b*d)

    cons39 = CustomConstraint(lambda d, c, a, b: cons_f39(d, c, a, b))

    def cons_f40(n, m):
        return Not(And(IntegerQ(m), NegQ(n)))

    cons40 = CustomConstraint(lambda n, m: cons_f40(n, m))

    def cons_f41(p, m):
        return ZeroQ(m + p)

    cons41 = CustomConstraint(lambda p, m: cons_f41(p, m))

    def cons_f42(d, c, a, b):
        return ZeroQ(a**S(2)*d + b**S(2)*c)

    cons42 = CustomConstraint(lambda d, c, a, b: cons_f42(d, c, a, b))

    def cons_f43(a):
        return PositiveQ(a)

    cons43 = CustomConstraint(lambda a: cons_f43(a))

    def cons_f44(d):
        return NegativeQ(d)

    cons44 = CustomConstraint(lambda d: cons_f44(d))

    def cons_f45(c, a, b):
        return ZeroQ(-S(4)*a*c + b**S(2))

    cons45 = CustomConstraint(lambda c, a, b: cons_f45(c, a, b))

    def cons_f46(n, n2):
        return ZeroQ(-S(2)*n + n2)

    cons46 = CustomConstraint(lambda n, n2: cons_f46(n, n2))

    def cons_f47(e, c, b, d):
        return ZeroQ(-b*e + S(2)*c*d)

    cons47 = CustomConstraint(lambda e, c, b, d: cons_f47(e, c, b, d))

    def cons_f48(e, x):
        return FreeQ(e, x)

    cons48 = CustomConstraint(lambda e, x: cons_f48(e, x))

    def cons_f49(q, p):
        return PosQ(-p + q)

    cons49 = CustomConstraint(lambda q, p: cons_f49(q, p))

    def cons_f50(q, x):
        return FreeQ(q, x)

    cons50 = CustomConstraint(lambda q, x: cons_f50(q, x))

    def cons_f51(r, p):
        return PosQ(-p + r)

    cons51 = CustomConstraint(lambda r, p: cons_f51(r, p))

    def cons_f52(r, x):
        return FreeQ(r, x)

    cons52 = CustomConstraint(lambda r, x: cons_f52(r, x))

    def cons_f53(n, m):
        return ZeroQ(m - n + S(1))

    cons53 = CustomConstraint(lambda n, m: cons_f53(n, m))

    def cons_f54(p):
        return NonzeroQ(p + S(1))

    cons54 = CustomConstraint(lambda p: cons_f54(p))

    def cons_f55(a1, b2, b1, a2):
        return ZeroQ(a1*b2 + a2*b1)

    cons55 = CustomConstraint(lambda a1, b2, b1, a2: cons_f55(a1, b2, b1, a2))

    def cons_f56(n, m):
        return ZeroQ(m - S(2)*n + S(1))

    cons56 = CustomConstraint(lambda n, m: cons_f56(n, m))

    def cons_f57(a1, x):
        return FreeQ(a1, x)

    cons57 = CustomConstraint(lambda a1, x: cons_f57(a1, x))

    def cons_f58(b1, x):
        return FreeQ(b1, x)

    cons58 = CustomConstraint(lambda b1, x: cons_f58(b1, x))

    def cons_f59(a2, x):
        return FreeQ(a2, x)

    cons59 = CustomConstraint(lambda a2, x: cons_f59(a2, x))

    def cons_f60(b2, x):
        return FreeQ(b2, x)

    cons60 = CustomConstraint(lambda b2, x: cons_f60(b2, x))

    def cons_f61(x, Qm):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return PolyQ(Qm, x)

    cons61 = CustomConstraint(lambda x, Qm: cons_f61(x, Qm))
    def With33(a, x, p, n, Pm, Qm, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        m = Expon(Pm, x)
        if And(Equal(Expon(Qm, x), m + S(-1)), ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x))):
            return True
        return False
    cons_with_33 = CustomConstraint(lambda a, x, p, n, Pm, Qm, b: With33(a, x, p, n, Pm, Qm, b))
    def With34(c, a, x, p, n, Pm, n2, Qm, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        m = Expon(Pm, x)
        if And(Equal(Expon(Qm, x), m + S(-1)), ZeroQ(-Qm*m*Coeff(Pm, x, m) + Coeff(Qm, x, m + S(-1))*D(Pm, x))):
            return True
        return False
    cons_with_34 = CustomConstraint(lambda c, a, x, p, n, Pm, n2, Qm, b: With34(c, a, x, p, n, Pm, n2, Qm, b))

    def cons_f62(m):
        return PositiveIntegerQ(m)

    cons62 = CustomConstraint(lambda m: cons_f62(m))

    def cons_f63(p):
        return NegativeIntegerQ(p)

    cons63 = CustomConstraint(lambda p: cons_f63(p))

    def cons_f64(x, Pq):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return PolyQ(Pq, x)

    cons64 = CustomConstraint(lambda x, Pq: cons_f64(x, Pq))

    def cons_f65(Qr, x):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return PolyQ(Qr, x)

    cons65 = CustomConstraint(lambda Qr, x: cons_f65(Qr, x))
    def With35(u, m, Qr, x, p, Pq):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        gcd = PolyGCD(Pq, Qr, x)
        if NonzeroQ(gcd + S(-1)):
            return True
        return False
    cons_with_35 = CustomConstraint(lambda u, m, Qr, x, p, Pq: With35(u, m, Qr, x, p, Pq))
    def With36(u, Qr, x, p, Pq):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        gcd = PolyGCD(Pq, Qr, x)
        if NonzeroQ(gcd + S(-1)):
            return True
        return False
    cons_with_36 = CustomConstraint(lambda u, Qr, x, p, Pq: With36(u, Qr, x, p, Pq))

    def cons_f66(m):
        return NonzeroQ(m + S(1))

    cons66 = CustomConstraint(lambda m: cons_f66(m))

    def cons_f67(x, a, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b), x)

    cons67 = CustomConstraint(lambda x, a, b: cons_f67(x, a, b))

    def cons_f68(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return LinearQ(u, x)

    cons68 = CustomConstraint(lambda x, u: cons_f68(x, u))

    def cons_f69(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return NonzeroQ(u - x)

    cons69 = CustomConstraint(lambda x, u: cons_f69(x, u))

    def cons_f70(d, c, a, b):
        return ZeroQ(a*d + b*c)

    cons70 = CustomConstraint(lambda d, c, a, b: cons_f70(d, c, a, b))

    def cons_f71(d, c, a, b):
        return NonzeroQ(-a*d + b*c)

    cons71 = CustomConstraint(lambda d, c, a, b: cons_f71(d, c, a, b))

    def cons_f72(n, m):
        return ZeroQ(m + n + S(2))

    cons72 = CustomConstraint(lambda n, m: cons_f72(n, m))

    def cons_f73(m):
        return PositiveIntegerQ(m + S(1)/2)

    cons73 = CustomConstraint(lambda m: cons_f73(m))

    def cons_f74(m):
        return NegativeIntegerQ(m + S(3)/2)

    cons74 = CustomConstraint(lambda m: cons_f74(m))

    def cons_f75(a, m, c):
        return Or(IntegerQ(m), And(PositiveQ(a), PositiveQ(c)))

    cons75 = CustomConstraint(lambda a, m, c: cons_f75(a, m, c))

    def cons_f76(c, a):
        return ZeroQ(a + c)

    cons76 = CustomConstraint(lambda c, a: cons_f76(c, a))

    def cons_f77(m):
        return Not(IntegerQ(S(2)*m))

    cons77 = CustomConstraint(lambda m: cons_f77(m))

    def cons_f78(d, c, a, b):
        return PosQ(b*d/(a*c))

    cons78 = CustomConstraint(lambda d, c, a, b: cons_f78(d, c, a, b))

    def cons_f79(m):
        return IntegerQ(m + S(1)/2)

    cons79 = CustomConstraint(lambda m: cons_f79(m))

    def cons_f80(n):
        return IntegerQ(n + S(1)/2)

    cons80 = CustomConstraint(lambda n: cons_f80(n))

    def cons_f81(n, m):
        return Less(S(0), m, n)

    cons81 = CustomConstraint(lambda n, m: cons_f81(n, m))

    def cons_f82(n, m):
        return Less(m, n, S(0))

    cons82 = CustomConstraint(lambda n, m: cons_f82(n, m))

    def cons_f83(n, c, m):
        return Or(Not(IntegerQ(n)), And(ZeroQ(c), LessEqual(S(7)*m + S(4)*n, S(0))), Less(S(9)*m + S(5)*n + S(5), S(0)), Greater(m + n + S(2), S(0)))

    cons83 = CustomConstraint(lambda n, c, m: cons_f83(n, c, m))

    def cons_f84(m):
        return NegativeIntegerQ(m)

    cons84 = CustomConstraint(lambda m: cons_f84(m))

    def cons_f85(n):
        return IntegerQ(n)

    cons85 = CustomConstraint(lambda n: cons_f85(n))

    def cons_f86(n, m):
        return Not(And(PositiveIntegerQ(n), Less(m + n + S(2), S(0))))

    cons86 = CustomConstraint(lambda n, m: cons_f86(n, m))

    def cons_f87(n):
        return RationalQ(n)

    cons87 = CustomConstraint(lambda n: cons_f87(n))

    def cons_f88(n):
        return Greater(n, S(0))

    cons88 = CustomConstraint(lambda n: cons_f88(n))

    def cons_f89(n):
        return Less(n, S(-1))

    cons89 = CustomConstraint(lambda n: cons_f89(n))

    def cons_f90(d, c, a, b):
        return PosQ((-a*d + b*c)/b)

    cons90 = CustomConstraint(lambda d, c, a, b: cons_f90(d, c, a, b))

    def cons_f91(d, c, a, b):
        return NegQ((-a*d + b*c)/b)

    cons91 = CustomConstraint(lambda d, c, a, b: cons_f91(d, c, a, b))

    def cons_f92(n):
        return Less(S(-1), n, S(0))

    cons92 = CustomConstraint(lambda n: cons_f92(n))

    def cons_f93(n, m):
        return RationalQ(m, n)

    cons93 = CustomConstraint(lambda n, m: cons_f93(n, m))

    def cons_f94(m):
        return Less(m, S(-1))

    cons94 = CustomConstraint(lambda m: cons_f94(m))

    def cons_f95(n, m):
        return Not(And(IntegerQ(n), Not(IntegerQ(m))))

    cons95 = CustomConstraint(lambda n, m: cons_f95(n, m))

    def cons_f96(n, m):
        return Not(And(IntegerQ(m + n), LessEqual(m + n + S(2), S(0)), Or(FractionQ(m), GreaterEqual(m + S(2)*n + S(1), S(0)))))

    cons96 = CustomConstraint(lambda n, m: cons_f96(n, m))

    def cons_f97(c, m, a, x, b, n, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return IntLinearcQ(a, b, c, d, m, n, x)

    cons97 = CustomConstraint(lambda c, m, a, x, b, n, d: cons_f97(c, m, a, x, b, n, d))

    def cons_f98(n, c, a, m):
        return Not(And(Less(n, S(-1)), Or(ZeroQ(a), And(NonzeroQ(c), Less(m, n), IntegerQ(n)))))

    cons98 = CustomConstraint(lambda n, c, a, m: cons_f98(n, c, a, m))

    def cons_f99(n, m):
        return Unequal(m + n + S(1), S(0))

    cons99 = CustomConstraint(lambda n, m: cons_f99(n, m))

    def cons_f100(n, m):
        return Not(And(PositiveIntegerQ(m), Or(Not(IntegerQ(n)), Less(S(0), m, n))))

    cons100 = CustomConstraint(lambda n, m: cons_f100(n, m))

    def cons_f101(n, m):
        return Not(And(IntegerQ(m + n), Less(m + n + S(2), S(0))))

    cons101 = CustomConstraint(lambda n, m: cons_f101(n, m))

    def cons_f102(d, b):
        return ZeroQ(b + d)

    cons102 = CustomConstraint(lambda d, b: cons_f102(d, b))

    def cons_f103(c, a):
        return PositiveQ(a + c)

    cons103 = CustomConstraint(lambda c, a: cons_f103(c, a))

    def cons_f104(d, c, a, b):
        return PositiveQ(-a*d + b*c)

    cons104 = CustomConstraint(lambda d, c, a, b: cons_f104(d, c, a, b))

    def cons_f105(b):
        return PositiveQ(b)

    cons105 = CustomConstraint(lambda b: cons_f105(b))

    def cons_f106(d, b):
        return ZeroQ(b - d)

    cons106 = CustomConstraint(lambda d, b: cons_f106(d, b))

    def cons_f107(m):
        return Less(S(-1), m, S(0))

    cons107 = CustomConstraint(lambda m: cons_f107(m))

    def cons_f108(m):
        return LessEqual(S(3), Denominator(m), S(4))

    cons108 = CustomConstraint(lambda m: cons_f108(m))

    def cons_f109(b, d):
        return PosQ(d/b)

    cons109 = CustomConstraint(lambda b, d: cons_f109(b, d))

    def cons_f110(b, d):
        return NegQ(d/b)

    cons110 = CustomConstraint(lambda b, d: cons_f110(b, d))

    def cons_f111(n, m):
        return Equal(m + n + S(1), S(0))

    cons111 = CustomConstraint(lambda n, m: cons_f111(n, m))

    def cons_f112(n, m):
        return LessEqual(Denominator(n), Denominator(m))

    cons112 = CustomConstraint(lambda n, m: cons_f112(n, m))

    def cons_f113(n, m):
        return NegativeIntegerQ(m + n + S(2))

    cons113 = CustomConstraint(lambda n, m: cons_f113(n, m))

    def cons_f114(n, m):
        return Or(SumSimplerQ(m, S(1)), Not(SumSimplerQ(n, S(1))))

    cons114 = CustomConstraint(lambda n, m: cons_f114(n, m))

    def cons_f115(n, c, b, d):
        return Or(IntegerQ(n), And(PositiveQ(c), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2)), PositiveQ(-d/(b*c))))))

    cons115 = CustomConstraint(lambda n, c, b, d: cons_f115(n, c, b, d))

    def cons_f116(b, c, m, d):
        return Or(IntegerQ(m), PositiveQ(-d/(b*c)))

    cons116 = CustomConstraint(lambda b, c, m, d: cons_f116(b, c, m, d))

    def cons_f117(c):
        return Not(PositiveQ(c))

    cons117 = CustomConstraint(lambda c: cons_f117(c))

    def cons_f118(c, b, d):
        return Not(PositiveQ(-d/(b*c)))

    cons118 = CustomConstraint(lambda c, b, d: cons_f118(c, b, d))

    def cons_f119(n, c, m, d):
        return Or(And(RationalQ(m), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2))))), Not(RationalQ(n)))

    cons119 = CustomConstraint(lambda n, c, m, d: cons_f119(n, c, m, d))

    def cons_f120(d, c, a, b):
        return PositiveQ(b/(-a*d + b*c))

    cons120 = CustomConstraint(lambda d, c, a, b: cons_f120(d, c, a, b))

    def cons_f121(c, m, a, d, n, b):
        return Or(RationalQ(m), Not(And(RationalQ(n), PositiveQ(-d/(-a*d + b*c)))))

    cons121 = CustomConstraint(lambda c, m, a, d, n, b: cons_f121(c, m, a, d, n, b))

    def cons_f122(n, m):
        return Or(RationalQ(m), Not(SimplerQ(n + S(1), m + S(1))))

    cons122 = CustomConstraint(lambda n, m: cons_f122(n, m))

    def cons_f123(x, u):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return NonzeroQ(Coefficient(u, x, S(0)))

    cons123 = CustomConstraint(lambda x, u: cons_f123(x, u))

    def cons_f124(n, m):
        return ZeroQ(m - n)

    cons124 = CustomConstraint(lambda n, m: cons_f124(n, m))

    def cons_f125(f, x):
        return FreeQ(f, x)

    cons125 = CustomConstraint(lambda f, x: cons_f125(f, x))

    def cons_f126(n, p):
        return NonzeroQ(n + p + S(2))

    cons126 = CustomConstraint(lambda n, p: cons_f126(n, p))

    def cons_f127(e, c, a, f, p, b, n, d):
        return ZeroQ(a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))

    cons127 = CustomConstraint(lambda e, c, a, f, p, b, n, d: cons_f127(e, c, a, f, p, b, n, d))

    def cons_f128(p):
        return PositiveIntegerQ(p)

    cons128 = CustomConstraint(lambda p: cons_f128(p))

    def cons_f129(f, e, a, b):
        return ZeroQ(a*f + b*e)

    cons129 = CustomConstraint(lambda f, e, a, b: cons_f129(f, e, a, b))

    def cons_f130(n, p):
        return Not(And(NegativeIntegerQ(n + p + S(2)), Greater(n + S(2)*p, S(0))))

    cons130 = CustomConstraint(lambda n, p: cons_f130(n, p))

    def cons_f131(n, p):
        return Or(NonzeroQ(n + S(1)), Equal(p, S(1)))

    cons131 = CustomConstraint(lambda n, p: cons_f131(n, p))

    def cons_f132(f, e, a, b):
        return NonzeroQ(a*f + b*e)

    cons132 = CustomConstraint(lambda f, e, a, b: cons_f132(f, e, a, b))

    def cons_f133(e, a, f, p, b, n, d):
        return Or(Not(IntegerQ(n)), Less(S(5)*n + S(9)*p, S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, d, e, f)))

    cons133 = CustomConstraint(lambda e, a, f, p, b, n, d: cons_f133(e, a, f, p, b, n, d))

    def cons_f134(e, c, a, f, d, p, n, b):
        return Or(NegativeIntegerQ(n, p), ZeroQ(p + S(-1)), And(PositiveIntegerQ(p), Or(Not(IntegerQ(n)), LessEqual(S(5)*n + S(9)*p + S(10), S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, c, d, e, f)))))

    cons134 = CustomConstraint(lambda e, c, a, f, d, p, n, b: cons_f134(e, c, a, f, d, p, n, b))

    def cons_f135(n, p):
        return ZeroQ(n + p + S(2))

    cons135 = CustomConstraint(lambda n, p: cons_f135(n, p))

    def cons_f136(n, p):
        return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1)))))

    cons136 = CustomConstraint(lambda n, p: cons_f136(n, p))

    def cons_f137(p):
        return Less(p, S(-1))

    cons137 = CustomConstraint(lambda p: cons_f137(p))

    def cons_f138(n, e, c, p):
        return Or(Not(And(RationalQ(n), Less(n, S(-1)))), IntegerQ(p), Not(Or(IntegerQ(n), Not(Or(ZeroQ(e), Not(Or(ZeroQ(c), Less(p, n))))))))

    cons138 = CustomConstraint(lambda n, e, c, p: cons_f138(n, e, c, p))

    def cons_f139(p):
        return SumSimplerQ(p, S(1))

    cons139 = CustomConstraint(lambda p: cons_f139(p))

    def cons_f140(n, p):
        return NonzeroQ(n + p + S(3))

    cons140 = CustomConstraint(lambda n, p: cons_f140(n, p))

    def cons_f141(e, c, a, f, p, b, n, d):
        return ZeroQ(-b*(c*f*(p + S(1)) + d*e*(n + S(1)))*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) + d*f*(a**S(2)*d*f*(n + p + S(3)) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e))*(n + p + S(2)))

    cons141 = CustomConstraint(lambda e, c, a, f, p, b, n, d: cons_f141(e, c, a, f, p, b, n, d))

    def cons_f142(n, m):
        return ZeroQ(m - n + S(-1))

    cons142 = CustomConstraint(lambda n, m: cons_f142(n, m))

    def cons_f143(m):
        return Not(PositiveIntegerQ(m))

    cons143 = CustomConstraint(lambda m: cons_f143(m))

    def cons_f144(n, p, m):
        return NonzeroQ(m + n + p + S(2))

    cons144 = CustomConstraint(lambda n, p, m: cons_f144(n, p, m))

    def cons_f145(p):
        return Less(S(0), p, S(1))

    cons145 = CustomConstraint(lambda p: cons_f145(p))

    def cons_f146(p):
        return Greater(p, S(1))

    cons146 = CustomConstraint(lambda p: cons_f146(p))

    def cons_f147(p):
        return Not(IntegerQ(p))

    cons147 = CustomConstraint(lambda p: cons_f147(p))

    def cons_f148(n):
        return PositiveIntegerQ(n)

    cons148 = CustomConstraint(lambda n: cons_f148(n))

    def cons_f149(p):
        return FractionQ(p)

    cons149 = CustomConstraint(lambda p: cons_f149(p))

    def cons_f150(n, m):
        return IntegersQ(m, n)

    cons150 = CustomConstraint(lambda n, m: cons_f150(n, m))

    def cons_f151(n, p, m):
        return Or(IntegerQ(p), And(Greater(m, S(0)), GreaterEqual(n, S(-1))))

    cons151 = CustomConstraint(lambda n, p, m: cons_f151(n, p, m))

    def cons_f152(n, p):
        return Or(And(RationalQ(n), Less(n, S(-1))), And(ZeroQ(n + p + S(3)), NonzeroQ(n + S(1)), Or(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1))))))

    cons152 = CustomConstraint(lambda n, p: cons_f152(n, p))

    def cons_f153(e, c, a, f, x, b, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c, d, e, f), x)

    cons153 = CustomConstraint(lambda e, c, a, f, x, b, d: cons_f153(e, c, a, f, x, b, d))

    def cons_f154(e, c, a, f, b, d):
        return ZeroQ(S(2)*b*d*e - f*(a*d + b*c))

    cons154 = CustomConstraint(lambda e, c, a, f, b, d: cons_f154(e, c, a, f, b, d))

    def cons_f155(n, m):
        return ZeroQ(m + n + S(1))

    cons155 = CustomConstraint(lambda n, m: cons_f155(n, m))

    def cons_f156(c, a, x, b, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return SimplerQ(a + b*x, c + d*x)

    cons156 = CustomConstraint(lambda c, a, x, b, d: cons_f156(c, a, x, b, d))

    def cons_f157(n, p, m):
        return ZeroQ(m + n + p + S(2))

    cons157 = CustomConstraint(lambda n, p, m: cons_f157(n, p, m))

    def cons_f158(p, m):
        return Not(And(SumSimplerQ(p, S(1)), Not(SumSimplerQ(m, S(1)))))

    cons158 = CustomConstraint(lambda p, m: cons_f158(p, m))

    def cons_f159(n, p, m):
        return ZeroQ(m + n + p + S(3))

    cons159 = CustomConstraint(lambda n, p, m: cons_f159(n, p, m))

    def cons_f160(e, c, m, a, f, d, p, n, b):
        return ZeroQ(a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))

    cons160 = CustomConstraint(lambda e, c, m, a, f, d, p, n, b: cons_f160(e, c, m, a, f, d, p, n, b))

    def cons_f161(m):
        return Or(And(RationalQ(m), Less(m, S(-1))), SumSimplerQ(m, S(1)))

    cons161 = CustomConstraint(lambda m: cons_f161(m))

    def cons_f162(n, p, m):
        return RationalQ(m, n, p)

    cons162 = CustomConstraint(lambda n, p, m: cons_f162(n, p, m))

    def cons_f163(p):
        return Greater(p, S(0))

    cons163 = CustomConstraint(lambda p: cons_f163(p))

    def cons_f164(n, p, m):
        return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), IntegersQ(m, n + p), IntegersQ(p, m + n))

    cons164 = CustomConstraint(lambda n, p, m: cons_f164(n, p, m))

    def cons_f165(n):
        return Greater(n, S(1))

    cons165 = CustomConstraint(lambda n: cons_f165(n))

    def cons_f166(m):
        return Greater(m, S(1))

    cons166 = CustomConstraint(lambda m: cons_f166(m))

    def cons_f167(n, p, m):
        return NonzeroQ(m + n + p + S(1))

    cons167 = CustomConstraint(lambda n, p, m: cons_f167(n, p, m))

    def cons_f168(m):
        return Greater(m, S(0))

    cons168 = CustomConstraint(lambda m: cons_f168(m))

    def cons_f169(n, p, m):
        return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), Or(IntegersQ(m, n + p), IntegersQ(p, m + n)))

    cons169 = CustomConstraint(lambda n, p, m: cons_f169(n, p, m))

    def cons_f170(n, p, m):
        return IntegersQ(S(2)*m, S(2)*n, S(2)*p)

    cons170 = CustomConstraint(lambda n, p, m: cons_f170(n, p, m))

    def cons_f171(n, p):
        return Or(IntegerQ(n), IntegersQ(S(2)*n, S(2)*p))

    cons171 = CustomConstraint(lambda n, p: cons_f171(n, p))

    def cons_f172(n, m):
        return PositiveIntegerQ(m + n + S(1))

    cons172 = CustomConstraint(lambda n, m: cons_f172(n, m))

    def cons_f173(n, m):
        return Or(And(RationalQ(m), Greater(m, S(0))), And(Not(RationalQ(m)), Or(SumSimplerQ(m, S(-1)), Not(SumSimplerQ(n, S(-1))))))

    cons173 = CustomConstraint(lambda n, m: cons_f173(n, m))

    def cons_f174(d, e, c, f):
        return PositiveQ(-f/(-c*f + d*e))

    cons174 = CustomConstraint(lambda d, e, c, f: cons_f174(d, e, c, f))

    def cons_f175(d, e, c, f):
        return Not(PositiveQ(-f/(-c*f + d*e)))

    cons175 = CustomConstraint(lambda d, e, c, f: cons_f175(d, e, c, f))

    def cons_f176(f, e, c, d):
        return NonzeroQ(-c*f + d*e)

    cons176 = CustomConstraint(lambda f, e, c, d: cons_f176(f, e, c, d))

    def cons_f177(c):
        return PositiveQ(c)

    cons177 = CustomConstraint(lambda c: cons_f177(c))

    def cons_f178(e):
        return PositiveQ(e)

    cons178 = CustomConstraint(lambda e: cons_f178(e))

    def cons_f179(d, b):
        return Not(NegativeQ(-b/d))

    cons179 = CustomConstraint(lambda d, b: cons_f179(d, b))

    def cons_f180(d, b):
        return NegativeQ(-b/d)

    cons180 = CustomConstraint(lambda d, b: cons_f180(d, b))

    def cons_f181(e, c):
        return Not(And(PositiveQ(c), PositiveQ(e)))

    cons181 = CustomConstraint(lambda e, c: cons_f181(e, c))

    def cons_f182(f, e, a, b):
        return PositiveQ(b/(-a*f + b*e))

    cons182 = CustomConstraint(lambda f, e, a, b: cons_f182(f, e, a, b))

    def cons_f183(d, c, a, b):
        return Not(NegativeQ(-(-a*d + b*c)/d))

    cons183 = CustomConstraint(lambda d, c, a, b: cons_f183(d, c, a, b))

    def cons_f184(e, c, a, f, x, d, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(And(SimplerQ(c + d*x, a + b*x), PositiveQ(-d/(-a*d + b*c)), PositiveQ(d/(-c*f + d*e)), Not(NegativeQ((-a*d + b*c)/b))))

    cons184 = CustomConstraint(lambda e, c, a, f, x, d, b: cons_f184(e, c, a, f, x, d, b))

    def cons_f185(e, c, a, f, b, d):
        return Not(And(PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e))))

    cons185 = CustomConstraint(lambda e, c, a, f, b, d: cons_f185(e, c, a, f, b, d))

    def cons_f186(f, d, b):
        return Or(PositiveQ(-b/d), NegativeQ(-b/f))

    cons186 = CustomConstraint(lambda f, d, b: cons_f186(f, d, b))

    def cons_f187(f, d, b):
        return Or(PosQ(-b/d), NegQ(-b/f))

    cons187 = CustomConstraint(lambda f, d, b: cons_f187(f, d, b))

    def cons_f188(e, a, f, x, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return SimplerQ(a + b*x, e + f*x)

    cons188 = CustomConstraint(lambda e, a, f, x, b: cons_f188(e, a, f, x, b))

    def cons_f189(e, c, a, f, b, d):
        return Or(PositiveQ(-(-a*d + b*c)/d), NegativeQ(-(-a*f + b*e)/f))

    cons189 = CustomConstraint(lambda e, c, a, f, b, d: cons_f189(e, c, a, f, b, d))

    def cons_f190(e, c, a, f, b, d):
        return Or(PosQ(-(-a*d + b*c)/d), NegQ(-(-a*f + b*e)/f))

    cons190 = CustomConstraint(lambda e, c, a, f, b, d: cons_f190(e, c, a, f, b, d))

    def cons_f191(e, c, a, f, b, d):
        return ZeroQ(-a*d*f - b*c*f + S(2)*b*d*e)

    cons191 = CustomConstraint(lambda e, c, a, f, b, d: cons_f191(e, c, a, f, b, d))

    def cons_f192(n, m):
        return PositiveIntegerQ(m - n)

    cons192 = CustomConstraint(lambda n, m: cons_f192(n, m))

    def cons_f193(n, m):
        return Or(PositiveIntegerQ(m), NegativeIntegerQ(m, n))

    cons193 = CustomConstraint(lambda n, m: cons_f193(n, m))

    def cons_f194(n, p, m):
        return NegativeIntegerQ(m + n + p + S(2))

    cons194 = CustomConstraint(lambda n, p, m: cons_f194(n, p, m))

    def cons_f195(n, p, m):
        return Or(SumSimplerQ(m, S(1)), And(Not(And(NonzeroQ(n + S(1)), SumSimplerQ(n, S(1)))), Not(And(NonzeroQ(p + S(1)), SumSimplerQ(p, S(1))))))

    cons195 = CustomConstraint(lambda n, p, m: cons_f195(n, p, m))

    def cons_f196(n):
        return NegativeIntegerQ(n)

    cons196 = CustomConstraint(lambda n: cons_f196(n))

    def cons_f197(e, p):
        return Or(IntegerQ(p), PositiveQ(e))

    cons197 = CustomConstraint(lambda e, p: cons_f197(e, p))

    def cons_f198(c, b, d):
        return PositiveQ(-d/(b*c))

    cons198 = CustomConstraint(lambda c, b, d: cons_f198(c, b, d))

    def cons_f199(e, c, f, p, d):
        return Or(IntegerQ(p), PositiveQ(d/(-c*f + d*e)))

    cons199 = CustomConstraint(lambda e, c, f, p, d: cons_f199(e, c, f, p, d))

    def cons_f200(c, a, x, d, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(And(PositiveQ(d/(a*d - b*c)), SimplerQ(c + d*x, a + b*x)))

    cons200 = CustomConstraint(lambda c, a, x, d, b: cons_f200(c, a, x, d, b))

    def cons_f201(d, c, a, b):
        return Not(PositiveQ(b/(-a*d + b*c)))

    cons201 = CustomConstraint(lambda d, c, a, b: cons_f201(d, c, a, b))

    def cons_f202(c, a, x, d, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(SimplerQ(c + d*x, a + b*x))

    cons202 = CustomConstraint(lambda c, a, x, d, b: cons_f202(c, a, x, d, b))

    def cons_f203(e, c, a, f, d, x, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(And(PositiveQ(d/(a*d - b*c)), PositiveQ(d/(-c*f + d*e)), SimplerQ(c + d*x, a + b*x)))

    cons203 = CustomConstraint(lambda e, c, a, f, d, x, b: cons_f203(e, c, a, f, d, x, b))

    def cons_f204(e, c, a, f, d, x, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(And(PositiveQ(f/(a*f - b*e)), PositiveQ(f/(c*f - d*e)), SimplerQ(e + f*x, a + b*x)))

    cons204 = CustomConstraint(lambda e, c, a, f, d, x, b: cons_f204(e, c, a, f, d, x, b))

    def cons_f205(f, e, a, b):
        return Not(PositiveQ(b/(-a*f + b*e)))

    cons205 = CustomConstraint(lambda f, e, a, b: cons_f205(f, e, a, b))

    def cons_f206(e, a, f, x, b):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return Not(SimplerQ(e + f*x, a + b*x))

    cons206 = CustomConstraint(lambda e, a, f, x, b: cons_f206(e, a, f, x, b))

    def cons_f207(n, m):
        return Or(PositiveIntegerQ(m), IntegersQ(m, n))

    cons207 = CustomConstraint(lambda n, m: cons_f207(n, m))

    def cons_f208(g, x):
        return FreeQ(g, x)

    cons208 = CustomConstraint(lambda g, x: cons_f208(g, x))

    def cons_f209(h, x):
        return FreeQ(h, x)

    cons209 = CustomConstraint(lambda h, x: cons_f209(h, x))

    def cons_f210(n, m):
        return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(m, S(1)))))

    cons210 = CustomConstraint(lambda n, m: cons_f210(n, m))

    def cons_f211(n, m):
        return Or(And(RationalQ(m), Less(m, S(-2))), And(ZeroQ(m + n + S(3)), Not(And(RationalQ(n), Less(n, S(-2))))))

    cons211 = CustomConstraint(lambda n, m: cons_f211(n, m))

    def cons_f212(m):
        return Or(And(RationalQ(m), Inequality(S(-2), LessEqual, m, Less, S(-1))), SumSimplerQ(m, S(1)))

    cons212 = CustomConstraint(lambda m: cons_f212(m))

    def cons_f213(n, m):
        return NonzeroQ(m + n + S(3))

    cons213 = CustomConstraint(lambda n, m: cons_f213(n, m))

    def cons_f214(n, m):
        return NonzeroQ(m + n + S(2))

    cons214 = CustomConstraint(lambda n, m: cons_f214(n, m))

    def cons_f215(n, p, m):
        return Or(IntegersQ(m, n, p), PositiveIntegerQ(n, p))

    cons215 = CustomConstraint(lambda n, p, m: cons_f215(n, p, m))

    def cons_f216(e, c, a, h, f, x, b, g, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c, d, e, f, g, h), x)

    cons216 = CustomConstraint(lambda e, c, a, h, f, x, b, g, d: cons_f216(e, c, a, h, f, x, b, g, d))

    def cons_f217(e, c, a, h, f, x, b, p, g, n, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c, d, e, f, g, h, n, p), x)

    cons217 = CustomConstraint(lambda e, c, a, h, f, x, b, p, g, n, d: cons_f217(e, c, a, h, f, x, b, p, g, n, d))

    def cons_f218(e, c, f, x, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return SimplerQ(c + d*x, e + f*x)

    cons218 = CustomConstraint(lambda e, c, f, x, d: cons_f218(e, c, f, x, d))

    def cons_f219(n, p, m):
        return Or(SumSimplerQ(m, S(1)), And(Not(SumSimplerQ(n, S(1))), Not(SumSimplerQ(p, S(1)))))

    cons219 = CustomConstraint(lambda n, p, m: cons_f219(n, p, m))

    def cons_f220(q, p):
        return IntegersQ(p, q)

    cons220 = CustomConstraint(lambda q, p: cons_f220(q, p))

    def cons_f221(q):
        return PositiveIntegerQ(q)

    cons221 = CustomConstraint(lambda q: cons_f221(q))

    def cons_f222(e, q, c, a, h, f, m, x, b, p, g, n, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c, d, e, f, g, h, m, n, p, q), x)

    cons222 = CustomConstraint(lambda e, q, c, a, h, f, m, x, b, p, g, n, d: cons_f222(e, q, c, a, h, f, m, x, b, p, g, n, d))

    def cons_f223(e, q, c, a, h, f, m, x, b, p, g, i, n, r, d):
        if isinstance(x, (int, Integer, float, Float)):
            return False
        return FreeQ(List(a, b, c, d, e, f, g, h, i, m, n, p, q, r), x)

    cons223 = CustomConstraint(lambda e, q, c, a, h, f, m, x, b, p, g, i, n, r, d: cons_f223(e, q, c, a, h, f, m, x, b, p, g, i, n, r, d))

    def cons_f224(i, x):
        return FreeQ(i, x)

    cons224 = CustomConstraint(lambda i, x: cons_f224(i, x))
