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
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist
    )
    from sympy import Integral, S, sqrt, And, Or

def cons_f1(a):
    return ZeroQ(a)

cons1 = CustomConstraint(cons_f1)

def cons_f2(a, x):
    return FreeQ(a, x)

cons2 = CustomConstraint(cons_f2)

def cons_f3(b, x):
    return FreeQ(b, x)

cons3 = CustomConstraint(cons_f3)

def cons_f4(n, x):
    return FreeQ(n, x)

cons4 = CustomConstraint(cons_f4)

def cons_f5(p, x):
    return FreeQ(p, x)

cons5 = CustomConstraint(cons_f5)

def cons_f6(n, j):
    return ZeroQ(j - S(2)*n)

cons6 = CustomConstraint(cons_f6)

def cons_f7(c, x):
    return FreeQ(c, x)

cons7 = CustomConstraint(cons_f7)

def cons_f8(b):
    return ZeroQ(b)

cons8 = CustomConstraint(cons_f8)

def cons_f9(c):
    return ZeroQ(c)

cons9 = CustomConstraint(cons_f9)

def cons_f10(v, x):
    return NFreeQ(v, x)

cons10 = CustomConstraint(cons_f10)

def cons_f11(x, Pm):
    return PolyQ(Pm, x)

cons11 = CustomConstraint(cons_f11)

def cons_f12(p):
    return Not(RationalQ(p))

cons12 = CustomConstraint(cons_f12)

def cons_f13(p):
    return RationalQ(p)

cons13 = CustomConstraint(cons_f13)

def cons_f14(x, c, b, a):
    return FreeQ(List(a, b, c), x)

cons14 = CustomConstraint(cons_f14)

def cons_f15(a):
    return EqQ(a**S(2), S(1))

cons15 = CustomConstraint(cons_f15)

def cons_f16(u):
    return SumQ(u)

cons16 = CustomConstraint(cons_f16)

def cons_f17(m):
    return IntegerQ(m)

cons17 = CustomConstraint(cons_f17)

def cons_f18(m):
    return Not(IntegerQ(m))

cons18 = CustomConstraint(cons_f18)

def cons_f19(n):
    return PositiveIntegerQ(n + S(1)/2)

cons19 = CustomConstraint(cons_f19)

def cons_f20(m, n):
    return IntegerQ(m + n)

cons20 = CustomConstraint(cons_f20)

def cons_f21(m, x):
    return FreeQ(m, x)

cons21 = CustomConstraint(cons_f21)

def cons_f22(n):
    return NegativeIntegerQ(n + S(-1)/2)

cons22 = CustomConstraint(cons_f22)

def cons_f23(n):
    return Not(IntegerQ(n))

cons23 = CustomConstraint(cons_f23)

def cons_f24(m, n):
    return Not(IntegerQ(m + n))

cons24 = CustomConstraint(cons_f24)

def cons_f25(d, c, b, a):
    return ZeroQ(-a*d + b*c)

cons25 = CustomConstraint(cons_f25)

def cons_f26(n, x, b, d, c, a):
    return Or(Not(IntegerQ(n)), SimplerQ(c + d*x, a + b*x))

cons26 = CustomConstraint(cons_f26)

def cons_f27(d, x):
    return FreeQ(d, x)

cons27 = CustomConstraint(cons_f27)

def cons_f28(d, b):
    return PositiveQ(b/d)

cons28 = CustomConstraint(cons_f28)

def cons_f29(m, n):
    return Not(Or(IntegerQ(m), IntegerQ(n)))

cons29 = CustomConstraint(cons_f29)

def cons_f30(m, n, d, b):
    return Not(Or(IntegerQ(m), IntegerQ(n), PositiveQ(b/d)))

cons30 = CustomConstraint(cons_f30)

def cons_f31(m):
    return RationalQ(m)

cons31 = CustomConstraint(cons_f31)

def cons_f32(m):
    return LessEqual(m, S(-1))

cons32 = CustomConstraint(cons_f32)

def cons_f33(B, A, b, C, a):
    return ZeroQ(A*b**S(2) - B*a*b + C*a**S(2))

cons33 = CustomConstraint(cons_f33)

def cons_f34(A, x):
    return FreeQ(A, x)

cons34 = CustomConstraint(cons_f34)

def cons_f35(B, x):
    return FreeQ(B, x)

cons35 = CustomConstraint(cons_f35)

def cons_f36(C, x):
    return FreeQ(C, x)

cons36 = CustomConstraint(cons_f36)

def cons_f37(n, q):
    return ZeroQ(n + q)

cons37 = CustomConstraint(cons_f37)

def cons_f38(p):
    return IntegerQ(p)

cons38 = CustomConstraint(cons_f38)

def cons_f39(d, c, b, a):
    return ZeroQ(a*c - b*d)

cons39 = CustomConstraint(cons_f39)

def cons_f40(m, n):
    return Not(And(IntegerQ(m), NegQ(n)))

cons40 = CustomConstraint(cons_f40)

def cons_f41(m, p):
    return ZeroQ(m + p)

cons41 = CustomConstraint(cons_f41)

def cons_f42(d, c, b, a):
    return ZeroQ(a**S(2)*d + b**S(2)*c)

cons42 = CustomConstraint(cons_f42)

def cons_f43(a):
    return PositiveQ(a)

cons43 = CustomConstraint(cons_f43)

def cons_f44(d):
    return NegativeQ(d)

cons44 = CustomConstraint(cons_f44)

def cons_f45(c, b, a):
    return ZeroQ(-S(4)*a*c + b**S(2))

cons45 = CustomConstraint(cons_f45)

def cons_f46(n2, n):
    return ZeroQ(-S(2)*n + n2)

cons46 = CustomConstraint(cons_f46)

def cons_f47(d, e, c, b):
    return ZeroQ(-b*e + S(2)*c*d)

cons47 = CustomConstraint(cons_f47)

def cons_f48(e, x):
    return FreeQ(e, x)

cons48 = CustomConstraint(cons_f48)

def cons_f49(q, p):
    return PosQ(-p + q)

cons49 = CustomConstraint(cons_f49)

def cons_f50(q, x):
    return FreeQ(q, x)

cons50 = CustomConstraint(cons_f50)

def cons_f51(p, r):
    return PosQ(-p + r)

cons51 = CustomConstraint(cons_f51)

def cons_f52(r, x):
    return FreeQ(r, x)

cons52 = CustomConstraint(cons_f52)

def cons_f53(m, n):
    return ZeroQ(m - n + S(1))

cons53 = CustomConstraint(cons_f53)

def cons_f54(p):
    return NonzeroQ(p + S(1))

cons54 = CustomConstraint(cons_f54)

def cons_f55(a1, a2, b1, b2):
    return ZeroQ(a1*b2 + a2*b1)

cons55 = CustomConstraint(cons_f55)

def cons_f56(m, n):
    return ZeroQ(m - S(2)*n + S(1))

cons56 = CustomConstraint(cons_f56)

def cons_f57(a1, x):
    return FreeQ(a1, x)

cons57 = CustomConstraint(cons_f57)

def cons_f58(b1, x):
    return FreeQ(b1, x)

cons58 = CustomConstraint(cons_f58)

def cons_f59(a2, x):
    return FreeQ(a2, x)

cons59 = CustomConstraint(cons_f59)

def cons_f60(b2, x):
    return FreeQ(b2, x)

cons60 = CustomConstraint(cons_f60)

def cons_f61(Qm, x):
    return PolyQ(Qm, x)

cons61 = CustomConstraint(cons_f61)

def cons_f62(m):
    return PositiveIntegerQ(m)

cons62 = CustomConstraint(cons_f62)

def cons_f63(p):
    return NegativeIntegerQ(p)

cons63 = CustomConstraint(cons_f63)

def cons_f64(x, Pq):
    return PolyQ(Pq, x)

cons64 = CustomConstraint(cons_f64)

def cons_f65(x, Qr):
    return PolyQ(Qr, x)

cons65 = CustomConstraint(cons_f65)

def cons_f66(m):
    return NonzeroQ(m + S(1))

cons66 = CustomConstraint(cons_f66)

def cons_f67(x, b, a):
    return FreeQ(List(a, b), x)

cons67 = CustomConstraint(cons_f67)

def cons_f68(x, u):
    return LinearQ(u, x)

cons68 = CustomConstraint(cons_f68)

def cons_f69(x, u):
    return NonzeroQ(u - x)

cons69 = CustomConstraint(cons_f69)

def cons_f70(d, c, b, a):
    return ZeroQ(a*d + b*c)

cons70 = CustomConstraint(cons_f70)

def cons_f71(d, c, b, a):
    return NonzeroQ(-a*d + b*c)

cons71 = CustomConstraint(cons_f71)

def cons_f72(m, n):
    return ZeroQ(m + n + S(2))

cons72 = CustomConstraint(cons_f72)

def cons_f73(m):
    return PositiveIntegerQ(m + S(1)/2)

cons73 = CustomConstraint(cons_f73)

def cons_f74(m):
    return NegativeIntegerQ(m + S(3)/2)

cons74 = CustomConstraint(cons_f74)

def cons_f75(m, c, a):
    return Or(IntegerQ(m), And(PositiveQ(a), PositiveQ(c)))

cons75 = CustomConstraint(cons_f75)

def cons_f76(c, a):
    return ZeroQ(a + c)

cons76 = CustomConstraint(cons_f76)

def cons_f77(m):
    return Not(IntegerQ(S(2)*m))

cons77 = CustomConstraint(cons_f77)

def cons_f78(d, c, b, a):
    return PosQ(b*d/(a*c))

cons78 = CustomConstraint(cons_f78)

def cons_f79(m):
    return IntegerQ(m + S(1)/2)

cons79 = CustomConstraint(cons_f79)

def cons_f80(n):
    return IntegerQ(n + S(1)/2)

cons80 = CustomConstraint(cons_f80)

def cons_f81(m, n):
    return Less(S(0), m, n)

cons81 = CustomConstraint(cons_f81)

def cons_f82(m, n):
    return Less(m, n, S(0))

cons82 = CustomConstraint(cons_f82)

def cons_f83(m, n, c):
    return Or(Not(IntegerQ(n)), And(ZeroQ(c), LessEqual(S(7)*m + S(4)*n, S(0))), Less(S(9)*m + S(5)*n + S(5), S(0)), Greater(m + n + S(2), S(0)))

cons83 = CustomConstraint(cons_f83)

def cons_f84(m):
    return NegativeIntegerQ(m)

cons84 = CustomConstraint(cons_f84)

def cons_f85(n):
    return IntegerQ(n)

cons85 = CustomConstraint(cons_f85)

def cons_f86(m, n):
    return Not(And(PositiveIntegerQ(n), Less(m + n + S(2), S(0))))

cons86 = CustomConstraint(cons_f86)

def cons_f87(n):
    return RationalQ(n)

cons87 = CustomConstraint(cons_f87)

def cons_f88(n):
    return Greater(n, S(0))

cons88 = CustomConstraint(cons_f88)

def cons_f89(n):
    return Less(n, S(-1))

cons89 = CustomConstraint(cons_f89)

def cons_f90(d, c, b, a):
    return PosQ((-a*d + b*c)/b)

cons90 = CustomConstraint(cons_f90)

def cons_f91(d, c, b, a):
    return NegQ((-a*d + b*c)/b)

cons91 = CustomConstraint(cons_f91)

def cons_f92(n):
    return Less(S(-1), n, S(0))

cons92 = CustomConstraint(cons_f92)

def cons_f93(m, n):
    return RationalQ(m, n)

cons93 = CustomConstraint(cons_f93)

def cons_f94(m):
    return Less(m, S(-1))

cons94 = CustomConstraint(cons_f94)

def cons_f95(m, n):
    return Not(And(IntegerQ(n), Not(IntegerQ(m))))

cons95 = CustomConstraint(cons_f95)

def cons_f96(m, n):
    return Not(And(IntegerQ(m + n), LessEqual(m + n + S(2), S(0)), Or(FractionQ(m), GreaterEqual(m + S(2)*n + S(1), S(0)))))

cons96 = CustomConstraint(cons_f96)

def cons_f97(m, n, x, b, d, c, a):
    return IntLinearcQ(a, b, c, d, m, n, x)

cons97 = CustomConstraint(cons_f97)

def cons_f98(m, n, c, a):
    return Not(And(Less(n, S(-1)), Or(ZeroQ(a), And(NonzeroQ(c), Less(m, n), IntegerQ(n)))))

cons98 = CustomConstraint(cons_f98)

def cons_f99(m, n):
    return Unequal(m + n + S(1), S(0))

cons99 = CustomConstraint(cons_f99)

def cons_f100(m, n):
    return Not(And(PositiveIntegerQ(m), Or(Not(IntegerQ(n)), Less(S(0), m, n))))

cons100 = CustomConstraint(cons_f100)

def cons_f101(m, n):
    return Not(And(IntegerQ(m + n), Less(m + n + S(2), S(0))))

cons101 = CustomConstraint(cons_f101)

def cons_f102(d, b):
    return ZeroQ(b + d)

cons102 = CustomConstraint(cons_f102)

def cons_f103(c, a):
    return PositiveQ(a + c)

cons103 = CustomConstraint(cons_f103)

def cons_f104(d, c, b, a):
    return PositiveQ(-a*d + b*c)

cons104 = CustomConstraint(cons_f104)

def cons_f105(b):
    return PositiveQ(b)

cons105 = CustomConstraint(cons_f105)

def cons_f106(d, b):
    return ZeroQ(b - d)

cons106 = CustomConstraint(cons_f106)

def cons_f107(m):
    return Less(S(-1), m, S(0))

cons107 = CustomConstraint(cons_f107)

def cons_f108(m):
    return LessEqual(S(3), Denominator(m), S(4))

cons108 = CustomConstraint(cons_f108)

def cons_f109(d, b):
    return PosQ(d/b)

cons109 = CustomConstraint(cons_f109)

def cons_f110(d, b):
    return NegQ(d/b)

cons110 = CustomConstraint(cons_f110)

def cons_f111(m, n):
    return Equal(m + n + S(1), S(0))

cons111 = CustomConstraint(cons_f111)

def cons_f112(m, n):
    return LessEqual(Denominator(n), Denominator(m))

cons112 = CustomConstraint(cons_f112)

def cons_f113(m, n):
    return NegativeIntegerQ(m + n + S(2))

cons113 = CustomConstraint(cons_f113)

def cons_f114(m, n):
    return Or(SumSimplerQ(m, S(1)), Not(SumSimplerQ(n, S(1))))

cons114 = CustomConstraint(cons_f114)

def cons_f115(d, n, c, b):
    return Or(IntegerQ(n), And(PositiveQ(c), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2)), PositiveQ(-d/(b*c))))))

cons115 = CustomConstraint(cons_f115)

def cons_f116(m, d, c, b):
    return Or(IntegerQ(m), PositiveQ(-d/(b*c)))

cons116 = CustomConstraint(cons_f116)

def cons_f117(c):
    return Not(PositiveQ(c))

cons117 = CustomConstraint(cons_f117)

def cons_f118(d, c, b):
    return Not(PositiveQ(-d/(b*c)))

cons118 = CustomConstraint(cons_f118)

def cons_f119(m, n, d, c):
    return Or(And(RationalQ(m), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2))))), Not(RationalQ(n)))

cons119 = CustomConstraint(cons_f119)

def cons_f120(d, c, b, a):
    return PositiveQ(b/(-a*d + b*c))

cons120 = CustomConstraint(cons_f120)

def cons_f121(m, n, b, d, c, a):
    return Or(RationalQ(m), Not(And(RationalQ(n), PositiveQ(-d/(-a*d + b*c)))))

cons121 = CustomConstraint(cons_f121)

def cons_f122(m, n):
    return Or(RationalQ(m), Not(SimplerQ(n + S(1), m + S(1))))

cons122 = CustomConstraint(cons_f122)

def cons_f123(x, u):
    return NonzeroQ(Coefficient(u, x, S(0)))

cons123 = CustomConstraint(cons_f123)

def cons_f124(m, n):
    return ZeroQ(m - n)

cons124 = CustomConstraint(cons_f124)

def cons_f125(f, x):
    return FreeQ(f, x)

cons125 = CustomConstraint(cons_f125)

def cons_f126(n, p):
    return NonzeroQ(n + p + S(2))

cons126 = CustomConstraint(cons_f126)

def cons_f127(n, p, e, b, d, f, c, a):
    return ZeroQ(a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))

cons127 = CustomConstraint(cons_f127)

def cons_f128(p):
    return PositiveIntegerQ(p)

cons128 = CustomConstraint(cons_f128)

def cons_f129(f, e, b, a):
    return ZeroQ(a*f + b*e)

cons129 = CustomConstraint(cons_f129)

def cons_f130(n, p):
    return Not(And(NegativeIntegerQ(n + p + S(2)), Greater(n + S(2)*p, S(0))))

cons130 = CustomConstraint(cons_f130)

def cons_f131(n, p):
    return Or(NonzeroQ(n + S(1)), Equal(p, S(1)))

cons131 = CustomConstraint(cons_f131)

def cons_f132(f, e, b, a):
    return NonzeroQ(a*f + b*e)

cons132 = CustomConstraint(cons_f132)

def cons_f133(n, p, e, b, d, f, a):
    return Or(Not(IntegerQ(n)), Less(S(5)*n + S(9)*p, S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, d, e, f)))

cons133 = CustomConstraint(cons_f133)

def cons_f134(n, p, e, b, d, f, c, a):
    return Or(NegativeIntegerQ(n, p), ZeroQ(p + S(-1)), And(PositiveIntegerQ(p), Or(Not(IntegerQ(n)), LessEqual(S(5)*n + S(9)*p + S(10), S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, c, d, e, f)))))

cons134 = CustomConstraint(cons_f134)

def cons_f135(n, p):
    return ZeroQ(n + p + S(2))

cons135 = CustomConstraint(cons_f135)

def cons_f136(n, p):
    return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1)))))

cons136 = CustomConstraint(cons_f136)

def cons_f137(p):
    return Less(p, S(-1))

cons137 = CustomConstraint(cons_f137)

def cons_f138(c, n, e, p):
    return Or(Not(And(RationalQ(n), Less(n, S(-1)))), IntegerQ(p), Not(Or(IntegerQ(n), Not(Or(ZeroQ(e), Not(Or(ZeroQ(c), Less(p, n))))))))

cons138 = CustomConstraint(cons_f138)

def cons_f139(p):
    return SumSimplerQ(p, S(1))

cons139 = CustomConstraint(cons_f139)

def cons_f140(n, p):
    return NonzeroQ(n + p + S(3))

cons140 = CustomConstraint(cons_f140)

def cons_f141(n, p, e, b, d, f, c, a):
    return ZeroQ(-b*(c*f*(p + S(1)) + d*e*(n + S(1)))*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) + d*f*(a**S(2)*d*f*(n + p + S(3)) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e))*(n + p + S(2)))

cons141 = CustomConstraint(cons_f141)

def cons_f142(m, n):
    return ZeroQ(m - n + S(-1))

cons142 = CustomConstraint(cons_f142)

def cons_f143(m):
    return Not(PositiveIntegerQ(m))

cons143 = CustomConstraint(cons_f143)

def cons_f144(m, n, p):
    return NonzeroQ(m + n + p + S(2))

cons144 = CustomConstraint(cons_f144)

def cons_f145(p):
    return Less(S(0), p, S(1))

cons145 = CustomConstraint(cons_f145)

def cons_f146(p):
    return Greater(p, S(1))

cons146 = CustomConstraint(cons_f146)

def cons_f147(p):
    return Not(IntegerQ(p))

cons147 = CustomConstraint(cons_f147)

def cons_f148(n):
    return PositiveIntegerQ(n)

cons148 = CustomConstraint(cons_f148)

def cons_f149(p):
    return FractionQ(p)

cons149 = CustomConstraint(cons_f149)

def cons_f150(m, n):
    return IntegersQ(m, n)

cons150 = CustomConstraint(cons_f150)

def cons_f151(m, p, n):
    return Or(IntegerQ(p), And(Greater(m, S(0)), GreaterEqual(n, S(-1))))

cons151 = CustomConstraint(cons_f151)

def cons_f152(n, p):
    return Or(And(RationalQ(n), Less(n, S(-1))), And(ZeroQ(n + p + S(3)), NonzeroQ(n + S(1)), Or(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1))))))

cons152 = CustomConstraint(cons_f152)

def cons_f153(e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f), x)

cons153 = CustomConstraint(cons_f153)

def cons_f154(e, b, d, f, c, a):
    return ZeroQ(S(2)*b*d*e - f*(a*d + b*c))

cons154 = CustomConstraint(cons_f154)

def cons_f155(m, n):
    return ZeroQ(m + n + S(1))

cons155 = CustomConstraint(cons_f155)

def cons_f156(x, b, d, c, a):
    return SimplerQ(a + b*x, c + d*x)

cons156 = CustomConstraint(cons_f156)

def cons_f157(m, n, p):
    return ZeroQ(m + n + p + S(2))

cons157 = CustomConstraint(cons_f157)

def cons_f158(m, p):
    return Not(And(SumSimplerQ(p, S(1)), Not(SumSimplerQ(m, S(1)))))

cons158 = CustomConstraint(cons_f158)

def cons_f159(m, n, p):
    return ZeroQ(m + n + p + S(3))

cons159 = CustomConstraint(cons_f159)

def cons_f160(m, n, p, e, b, d, f, c, a):
    return ZeroQ(a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))

cons160 = CustomConstraint(cons_f160)

def cons_f161(m):
    return Or(And(RationalQ(m), Less(m, S(-1))), SumSimplerQ(m, S(1)))

cons161 = CustomConstraint(cons_f161)

def cons_f162(m, n, p):
    return RationalQ(m, n, p)

cons162 = CustomConstraint(cons_f162)

def cons_f163(p):
    return Greater(p, S(0))

cons163 = CustomConstraint(cons_f163)

def cons_f164(m, n, p):
    return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), IntegersQ(m, n + p), IntegersQ(p, m + n))

cons164 = CustomConstraint(cons_f164)

def cons_f165(n):
    return Greater(n, S(1))

cons165 = CustomConstraint(cons_f165)

def cons_f166(m):
    return Greater(m, S(1))

cons166 = CustomConstraint(cons_f166)

def cons_f167(m, n, p):
    return NonzeroQ(m + n + p + S(1))

cons167 = CustomConstraint(cons_f167)

def cons_f168(m):
    return Greater(m, S(0))

cons168 = CustomConstraint(cons_f168)

def cons_f169(m, n, p):
    return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), Or(IntegersQ(m, n + p), IntegersQ(p, m + n)))

cons169 = CustomConstraint(cons_f169)

def cons_f170(m, n, p):
    return IntegersQ(S(2)*m, S(2)*n, S(2)*p)

cons170 = CustomConstraint(cons_f170)

def cons_f171(n, p):
    return Or(IntegerQ(n), IntegersQ(S(2)*n, S(2)*p))

cons171 = CustomConstraint(cons_f171)

def cons_f172(m, n):
    return PositiveIntegerQ(m + n + S(1))

cons172 = CustomConstraint(cons_f172)

def cons_f173(m, n):
    return Or(And(RationalQ(m), Greater(m, S(0))), And(Not(RationalQ(m)), Or(SumSimplerQ(m, S(-1)), Not(SumSimplerQ(n, S(-1))))))

cons173 = CustomConstraint(cons_f173)

def cons_f174(f, d, c, e):
    return PositiveQ(-f/(-c*f + d*e))

cons174 = CustomConstraint(cons_f174)

def cons_f175(f, d, c, e):
    return Not(PositiveQ(-f/(-c*f + d*e)))

cons175 = CustomConstraint(cons_f175)

def cons_f176(d, f, e, c):
    return NonzeroQ(-c*f + d*e)

cons176 = CustomConstraint(cons_f176)

def cons_f177(c):
    return PositiveQ(c)

cons177 = CustomConstraint(cons_f177)

def cons_f178(e):
    return PositiveQ(e)

cons178 = CustomConstraint(cons_f178)

def cons_f179(d, b):
    return Not(NegativeQ(-b/d))

cons179 = CustomConstraint(cons_f179)

def cons_f180(d, b):
    return NegativeQ(-b/d)

cons180 = CustomConstraint(cons_f180)

def cons_f181(e, c):
    return Not(And(PositiveQ(c), PositiveQ(e)))

cons181 = CustomConstraint(cons_f181)

def cons_f182(f, e, b, a):
    return PositiveQ(b/(-a*f + b*e))

cons182 = CustomConstraint(cons_f182)

def cons_f183(d, c, b, a):
    return Not(NegativeQ(-(-a*d + b*c)/d))

cons183 = CustomConstraint(cons_f183)

def cons_f184(x, e, b, d, f, c, a):
    return Not(And(SimplerQ(c + d*x, a + b*x), PositiveQ(-d/(-a*d + b*c)), PositiveQ(d/(-c*f + d*e)), Not(NegativeQ((-a*d + b*c)/b))))

cons184 = CustomConstraint(cons_f184)

def cons_f185(e, b, d, f, c, a):
    return Not(And(PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e))))

cons185 = CustomConstraint(cons_f185)

def cons_f186(d, f, b):
    return Or(PositiveQ(-b/d), NegativeQ(-b/f))

cons186 = CustomConstraint(cons_f186)

def cons_f187(d, f, b):
    return Or(PosQ(-b/d), NegQ(-b/f))

cons187 = CustomConstraint(cons_f187)

def cons_f188(e, x, b, f, a):
    return SimplerQ(a + b*x, e + f*x)

cons188 = CustomConstraint(cons_f188)

def cons_f189(e, b, d, f, c, a):
    return Or(PositiveQ(-(-a*d + b*c)/d), NegativeQ(-(-a*f + b*e)/f))

cons189 = CustomConstraint(cons_f189)

def cons_f190(e, b, d, f, c, a):
    return Or(PosQ(-(-a*d + b*c)/d), NegQ(-(-a*f + b*e)/f))

cons190 = CustomConstraint(cons_f190)

def cons_f191(e, b, d, f, c, a):
    return ZeroQ(-a*d*f - b*c*f + S(2)*b*d*e)

cons191 = CustomConstraint(cons_f191)

def cons_f192(m, n):
    return PositiveIntegerQ(m - n)

cons192 = CustomConstraint(cons_f192)

def cons_f193(m, n):
    return Or(PositiveIntegerQ(m), NegativeIntegerQ(m, n))

cons193 = CustomConstraint(cons_f193)

def cons_f194(m, n, p):
    return NegativeIntegerQ(m + n + p + S(2))

cons194 = CustomConstraint(cons_f194)

def cons_f195(m, n, p):
    return Or(SumSimplerQ(m, S(1)), And(Not(And(NonzeroQ(n + S(1)), SumSimplerQ(n, S(1)))), Not(And(NonzeroQ(p + S(1)), SumSimplerQ(p, S(1))))))

cons195 = CustomConstraint(cons_f195)

def cons_f196(n):
    return NegativeIntegerQ(n)

cons196 = CustomConstraint(cons_f196)

def cons_f197(p, e):
    return Or(IntegerQ(p), PositiveQ(e))

cons197 = CustomConstraint(cons_f197)

def cons_f198(d, c, b):
    return PositiveQ(-d/(b*c))

cons198 = CustomConstraint(cons_f198)

def cons_f199(p, e, d, f, c):
    return Or(IntegerQ(p), PositiveQ(d/(-c*f + d*e)))

cons199 = CustomConstraint(cons_f199)

def cons_f200(x, b, d, c, a):
    return Not(And(PositiveQ(d/(a*d - b*c)), SimplerQ(c + d*x, a + b*x)))

cons200 = CustomConstraint(cons_f200)

def cons_f201(d, c, b, a):
    return Not(PositiveQ(b/(-a*d + b*c)))

cons201 = CustomConstraint(cons_f201)

def cons_f202(x, b, d, c, a):
    return Not(SimplerQ(c + d*x, a + b*x))

cons202 = CustomConstraint(cons_f202)

def cons_f203(e, x, b, d, f, c, a):
    return Not(And(PositiveQ(d/(a*d - b*c)), PositiveQ(d/(-c*f + d*e)), SimplerQ(c + d*x, a + b*x)))

cons203 = CustomConstraint(cons_f203)

def cons_f204(d, e, x, b, f, c, a):
    return Not(And(PositiveQ(f/(a*f - b*e)), PositiveQ(f/(c*f - d*e)), SimplerQ(e + f*x, a + b*x)))

cons204 = CustomConstraint(cons_f204)

def cons_f205(f, e, b, a):
    return Not(PositiveQ(b/(-a*f + b*e)))

cons205 = CustomConstraint(cons_f205)

def cons_f206(x, e, b, f, a):
    return Not(SimplerQ(e + f*x, a + b*x))

cons206 = CustomConstraint(cons_f206)

def cons_f207(m, n):
    return Or(PositiveIntegerQ(m), IntegersQ(m, n))

cons207 = CustomConstraint(cons_f207)

def cons_f208(g, x):
    return FreeQ(g, x)

cons208 = CustomConstraint(cons_f208)

def cons_f209(h, x):
    return FreeQ(h, x)

cons209 = CustomConstraint(cons_f209)

def cons_f210(m, n):
    return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(m, S(1)))))

cons210 = CustomConstraint(cons_f210)

def cons_f211(m, n):
    return Or(And(RationalQ(m), Less(m, S(-2))), And(ZeroQ(m + n + S(3)), Not(And(RationalQ(n), Less(n, S(-2))))))

cons211 = CustomConstraint(cons_f211)

def cons_f212(m):
    return Or(And(RationalQ(m), Inequality(S(-2), LessEqual, m, Less, S(-1))), SumSimplerQ(m, S(1)))

cons212 = CustomConstraint(cons_f212)

def cons_f213(m, n):
    return NonzeroQ(m + n + S(3))

cons213 = CustomConstraint(cons_f213)

def cons_f214(m, n):
    return NonzeroQ(m + n + S(2))

cons214 = CustomConstraint(cons_f214)

def cons_f215(m, n, p):
    return Or(IntegersQ(m, n, p), PositiveIntegerQ(n, p))

cons215 = CustomConstraint(cons_f215)

def cons_f216(h, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, h), x)

cons216 = CustomConstraint(cons_f216)

def cons_f217(n, h, p, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, h, n, p), x)

cons217 = CustomConstraint(cons_f217)

def cons_f218(x, e, d, f, c):
    return SimplerQ(c + d*x, e + f*x)

cons218 = CustomConstraint(cons_f218)

def cons_f219(m, n, p):
    return Or(SumSimplerQ(m, S(1)), And(Not(SumSimplerQ(n, S(1))), Not(SumSimplerQ(p, S(1)))))

cons219 = CustomConstraint(cons_f219)

def cons_f220(p, q):
    return IntegersQ(p, q)

cons220 = CustomConstraint(cons_f220)

def cons_f221(q):
    return PositiveIntegerQ(q)

cons221 = CustomConstraint(cons_f221)

def cons_f222(m, n, h, p, q, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, h, m, n, p, q), x)

cons222 = CustomConstraint(cons_f222)

def cons_f223(m, n, h, p, q, i, g, e, x, b, d, f, r, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, h, i, m, n, p, q, r), x)

cons223 = CustomConstraint(cons_f223)

def cons_f224(i, x):
    return FreeQ(i, x)

cons224 = CustomConstraint(cons_f224)

def cons_f225(p):
    return NonzeroQ(S(2)*p + S(1))

cons225 = CustomConstraint(cons_f225)

def cons_f226(c, b, a):
    return NonzeroQ(-S(4)*a*c + b**S(2))

cons226 = CustomConstraint(cons_f226)

def cons_f227(c, b, a):
    return PerfectSquareQ(-S(4)*a*c + b**S(2))

cons227 = CustomConstraint(cons_f227)

def cons_f228(c, b, a):
    return Not(PerfectSquareQ(-S(4)*a*c + b**S(2)))

cons228 = CustomConstraint(cons_f228)

def cons_f229(p):
    return IntegerQ(S(4)*p)

cons229 = CustomConstraint(cons_f229)

def cons_f230(p):
    return Unequal(p, S(-3)/2)

cons230 = CustomConstraint(cons_f230)

def cons_f231(c, b, a):
    return PosQ(-S(4)*a*c + b**S(2))

cons231 = CustomConstraint(cons_f231)

def cons_f232(c, b, a):
    return PositiveQ(S(4)*a - b**S(2)/c)

cons232 = CustomConstraint(cons_f232)

def cons_f233(x, c, b):
    return FreeQ(List(b, c), x)

cons233 = CustomConstraint(cons_f233)

def cons_f234(p):
    return LessEqual(S(3), Denominator(p), S(4))

cons234 = CustomConstraint(cons_f234)

def cons_f235(p):
    return Not(IntegerQ(S(4)*p))

cons235 = CustomConstraint(cons_f235)

def cons_f236(m):
    return IntegerQ(m/S(2) + S(1)/2)

cons236 = CustomConstraint(cons_f236)

def cons_f237(m, p):
    return ZeroQ(m + S(2)*p + S(1))

cons237 = CustomConstraint(cons_f237)

def cons_f238(m, p):
    return NonzeroQ(m + S(2)*p + S(1))

cons238 = CustomConstraint(cons_f238)

def cons_f239(d, e, c, b):
    return NonzeroQ(-b*e + S(2)*c*d)

cons239 = CustomConstraint(cons_f239)

def cons_f240(m, p):
    return ZeroQ(m + S(2)*p + S(2))

cons240 = CustomConstraint(cons_f240)

def cons_f241(m):
    return NonzeroQ(m + S(2))

cons241 = CustomConstraint(cons_f241)

def cons_f242(m, p):
    return ZeroQ(m + S(2)*p + S(3))

cons242 = CustomConstraint(cons_f242)

def cons_f243(p):
    return NonzeroQ(p + S(3)/2)

cons243 = CustomConstraint(cons_f243)

def cons_f244(m, p):
    return RationalQ(m, p)

cons244 = CustomConstraint(cons_f244)

def cons_f245(m):
    return Inequality(S(-2), LessEqual, m, Less, S(-1))

cons245 = CustomConstraint(cons_f245)

def cons_f246(p):
    return IntegerQ(S(2)*p)

cons246 = CustomConstraint(cons_f246)

def cons_f247(m):
    return Less(m, S(-2))

cons247 = CustomConstraint(cons_f247)

def cons_f248(m, p):
    return Not(And(NegativeIntegerQ(m + S(2)*p + S(3)), Greater(m + S(3)*p + S(3), S(0))))

cons248 = CustomConstraint(cons_f248)

def cons_f249(m, p):
    return NonzeroQ(m + S(2)*p)

cons249 = CustomConstraint(cons_f249)

def cons_f250(m):
    return Not(And(RationalQ(m), Less(m, S(-2))))

cons250 = CustomConstraint(cons_f250)

def cons_f251(m, p):
    return Not(And(IntegerQ(m), Less(S(0), m, S(2)*p)))

cons251 = CustomConstraint(cons_f251)

def cons_f252(m):
    return Inequality(S(0), Less, m, LessEqual, S(1))

cons252 = CustomConstraint(cons_f252)

def cons_f253(m, p):
    return NonzeroQ(m + p + S(1))

cons253 = CustomConstraint(cons_f253)

def cons_f254(m, p):
    return Or(Not(RationalQ(p)), Inequality(S(-1), LessEqual, p, Less, S(0)), And(IntegerQ(m), Less(S(0), m, S(2)*p)), And(Equal(m, S(1)/2), Less(p, S(0))))

cons254 = CustomConstraint(cons_f254)

def cons_f255(m, p):
    return Or(IntegerQ(m), IntegerQ(S(2)*p))

cons255 = CustomConstraint(cons_f255)

def cons_f256(e, b, d, c, a):
    return ZeroQ(a*e**S(2) - b*d*e + c*d**S(2))

cons256 = CustomConstraint(cons_f256)

def cons_f257(d, e, c, a):
    return ZeroQ(a*e**S(2) + c*d**S(2))

cons257 = CustomConstraint(cons_f257)

def cons_f258(d, p, m, a):
    return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), IntegerQ(m + p)))

cons258 = CustomConstraint(cons_f258)

def cons_f259(m, p):
    return Or(Less(S(0), -m, p), Less(p, -m, S(0)))

cons259 = CustomConstraint(cons_f259)

def cons_f260(m):
    return Unequal(m, S(2))

cons260 = CustomConstraint(cons_f260)

def cons_f261(m):
    return Unequal(m, S(-1))

cons261 = CustomConstraint(cons_f261)

def cons_f262(m, p):
    return PositiveIntegerQ(m + p)

cons262 = CustomConstraint(cons_f262)

def cons_f263(m, p):
    return NegativeIntegerQ(m + S(2)*p + S(2))

cons263 = CustomConstraint(cons_f263)

def cons_f264(m, p):
    return Or(Less(m, S(-2)), ZeroQ(m + S(2)*p + S(1)))

cons264 = CustomConstraint(cons_f264)

def cons_f265(m, p):
    return Or(Inequality(S(-2), LessEqual, m, Less, S(0)), Equal(m + p + S(1), S(0)))

cons265 = CustomConstraint(cons_f265)

def cons_f266(m):
    return GreaterEqual(m, S(1))

cons266 = CustomConstraint(cons_f266)

def cons_f267(m):
    return Less(m, S(0))

cons267 = CustomConstraint(cons_f267)

def cons_f268(d):
    return PositiveQ(d)

cons268 = CustomConstraint(cons_f268)

def cons_f269(m, p):
    return Not(And(ZeroQ(m + S(-3)), Unequal(p, S(1))))

cons269 = CustomConstraint(cons_f269)

def cons_f270(m, p):
    return NonzeroQ(m + S(2)*p + S(3))

cons270 = CustomConstraint(cons_f270)

def cons_f271(m, p):
    return Not(And(EvenQ(m), Less(m + S(2)*p + S(3), S(0))))

cons271 = CustomConstraint(cons_f271)

def cons_f272(m):
    return Not(And(RationalQ(m), Less(m, S(-1))))

cons272 = CustomConstraint(cons_f272)

def cons_f273(m, p):
    return Not(And(PositiveIntegerQ(m/S(2) + S(-1)/2), Or(Not(IntegerQ(p)), Less(m, S(2)*p))))

cons273 = CustomConstraint(cons_f273)

def cons_f274(m):
    return Not(And(RationalQ(m), Greater(m, S(1))))

cons274 = CustomConstraint(cons_f274)

def cons_f275(c, b, a):
    return NegativeQ(c/(-S(4)*a*c + b**S(2)))

cons275 = CustomConstraint(cons_f275)

def cons_f276(m):
    return EqQ(m**S(2), S(1)/4)

cons276 = CustomConstraint(cons_f276)

def cons_f277(m, p):
    return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), OddQ(m))

cons277 = CustomConstraint(cons_f277)

def cons_f278(m, p):
    return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), IntegerQ(m/S(2) + p + S(3)/2))

cons278 = CustomConstraint(cons_f278)

def cons_f279(e, b, d, c, a):
    return NonzeroQ(a*e**S(2) - b*d*e + c*d**S(2))

cons279 = CustomConstraint(cons_f279)

def cons_f280(d, e, c, a):
    return NonzeroQ(a*e**S(2) + c*d**S(2))

cons280 = CustomConstraint(cons_f280)

def cons_f281(m, p):
    return Not(And(ZeroQ(m + S(-1)), Greater(p, S(1))))

cons281 = CustomConstraint(cons_f281)

def cons_f282(c, b, a):
    return NiceSqrtQ(-S(4)*a*c + b**S(2))

cons282 = CustomConstraint(cons_f282)

def cons_f283(c, a):
    return NiceSqrtQ(-a*c)

cons283 = CustomConstraint(cons_f283)

def cons_f284(c, b, a):
    return Not(NiceSqrtQ(-S(4)*a*c + b**S(2)))

cons284 = CustomConstraint(cons_f284)

def cons_f285(c, a):
    return Not(NiceSqrtQ(-a*c))

cons285 = CustomConstraint(cons_f285)

def cons_f286(d, m):
    return Or(NonzeroQ(d), Greater(m, S(2)))

cons286 = CustomConstraint(cons_f286)

def cons_f287(p):
    return Not(And(RationalQ(p), LessEqual(p, S(-1))))

cons287 = CustomConstraint(cons_f287)

def cons_f288(d, e, b, a):
    return ZeroQ(a*e + b*d)

cons288 = CustomConstraint(cons_f288)

def cons_f289(d, e, c, b):
    return ZeroQ(b*e + c*d)

cons289 = CustomConstraint(cons_f289)

def cons_f290(m, p):
    return PositiveIntegerQ(m - p + S(1))

cons290 = CustomConstraint(cons_f290)

def cons_f291(d, e, c, b):
    return NonzeroQ(-b*e + c*d)

cons291 = CustomConstraint(cons_f291)

def cons_f292(m):
    return Equal(m**S(2), S(1)/4)

cons292 = CustomConstraint(cons_f292)

def cons_f293(c):
    return NegativeQ(c)

cons293 = CustomConstraint(cons_f293)

def cons_f294(b):
    return RationalQ(b)

cons294 = CustomConstraint(cons_f294)

def cons_f295(m):
    return ZeroQ(m**S(2) + S(-1)/4)

cons295 = CustomConstraint(cons_f295)

def cons_f296(m, p):
    return Equal(m + S(2)*p + S(2), S(0))

cons296 = CustomConstraint(cons_f296)

def cons_f297(x, e, d, c, a):
    return FreeQ(List(a, c, d, e), x)

cons297 = CustomConstraint(cons_f297)

def cons_f298(m, p):
    return Or(IntegerQ(p), And(RationalQ(m), Less(m, S(-1))))

cons298 = CustomConstraint(cons_f298)

def cons_f299(m, p):
    return Not(NegativeIntegerQ(m + S(2)*p + S(1)))

cons299 = CustomConstraint(cons_f299)

def cons_f300(m, p, e, x, b, d, c, a):
    return IntQuadraticQ(a, b, c, d, e, m, p, x)

cons300 = CustomConstraint(cons_f300)

def cons_f301(m, p, e, x, d, c, a):
    return IntQuadraticQ(a, S(0), c, d, e, m, p, x)

cons301 = CustomConstraint(cons_f301)

def cons_f302(m):
    return Or(Not(RationalQ(m)), Less(m, S(1)))

cons302 = CustomConstraint(cons_f302)

def cons_f303(m, p):
    return Not(NegativeIntegerQ(m + S(2)*p))

cons303 = CustomConstraint(cons_f303)

def cons_f304(m, p):
    return Or(Less(m, S(1)), And(NegativeIntegerQ(m + S(2)*p + S(3)), Unequal(m, S(2))))

cons304 = CustomConstraint(cons_f304)

def cons_f305(m):
    return If(RationalQ(m), Greater(m, S(1)), SumSimplerQ(m, S(-2)))

cons305 = CustomConstraint(cons_f305)

def cons_f306(m, p, e, x, b, d, c, a):
    return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, b, c, d, e, m, p, x)), And(SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

cons306 = CustomConstraint(cons_f306)

def cons_f307(m, p, e, x, d, c, a):
    return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, S(0), c, d, e, m, p, x)), And(SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

cons307 = CustomConstraint(cons_f307)

def cons_f308(e, b, d, c, a):
    return ZeroQ(-S(3)*a*c*e**S(2) + b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

cons308 = CustomConstraint(cons_f308)

def cons_f309(d, e, c, b):
    return PosQ(c*e**S(2)*(-b*e + S(2)*c*d))

cons309 = CustomConstraint(cons_f309)

def cons_f310(d, e, c, a):
    return ZeroQ(-S(3)*a*e**S(2) + c*d**S(2))

cons310 = CustomConstraint(cons_f310)

def cons_f311(d, e, c, b):
    return NegQ(c*e**S(2)*(-b*e + S(2)*c*d))

cons311 = CustomConstraint(cons_f311)

def cons_f312(e, b, d, c, a):
    return ZeroQ(S(9)*a*c*e**S(2) - S(2)*b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

cons312 = CustomConstraint(cons_f312)

def cons_f313(c, b, a):
    return Not(PositiveQ(S(4)*a - b**S(2)/c))

cons313 = CustomConstraint(cons_f313)

def cons_f314(p):
    return Not(IntegerQ(S(2)*p))

cons314 = CustomConstraint(cons_f314)

def cons_f315(f, g, e, d):
    return NonzeroQ(-d*g + e*f)

cons315 = CustomConstraint(cons_f315)

def cons_f316(f, g, c, b):
    return ZeroQ(-b*g + S(2)*c*f)

cons316 = CustomConstraint(cons_f316)

def cons_f317(m):
    return Not(And(RationalQ(m), Greater(m, S(0))))

cons317 = CustomConstraint(cons_f317)

def cons_f318(m, p):
    return Or(Not(RationalQ(p)), And(Greater(p, S(0)), Or(Not(IntegerQ(m)), GreaterEqual(m, -S(2)*p + S(-2)), Less(m, -S(4)*p + S(-4)))))

cons318 = CustomConstraint(cons_f318)

def cons_f319(m, p):
    return NonzeroQ(m + S(2)*p + S(2))

cons319 = CustomConstraint(cons_f319)

def cons_f320(m, p):
    return Or(Not(RationalQ(p)), Less(m, S(2)*p + S(2)))

cons320 = CustomConstraint(cons_f320)

def cons_f321(f, g, c, b):
    return NonzeroQ(-b*g + S(2)*c*f)

cons321 = CustomConstraint(cons_f321)

def cons_f322(p):
    return Less(p, S(0))

cons322 = CustomConstraint(cons_f322)

def cons_f323(m, p, e, b, d, c):
    return Or(And(ZeroQ(m + S(2)*p + S(2)), NonzeroQ(m + S(1))), And(ZeroQ(-b*e + S(2)*c*d), NonzeroQ(m + S(-1))))

cons323 = CustomConstraint(cons_f323)

def cons_f324(m, d, g, x, e, f):
    return Not(And(ZeroQ(m + S(-1)), SimplerQ(f + g*x, d + e*x)))

cons324 = CustomConstraint(cons_f324)

def cons_f325(d, p, m, a):
    return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), ZeroQ(m + p)))

cons325 = CustomConstraint(cons_f325)

def cons_f326(m, p, g, e, b, d, f, c):
    return ZeroQ(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))

cons326 = CustomConstraint(cons_f326)

def cons_f327(m, p, g, e, d, f):
    return ZeroQ(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))

cons327 = CustomConstraint(cons_f327)

def cons_f328(m):
    return SumSimplerQ(m, S(-1))

cons328 = CustomConstraint(cons_f328)

def cons_f329(m, p):
    return Or(And(RationalQ(m), Less(m, S(-1)), Not(PositiveIntegerQ(m + p + S(1)))), And(RationalQ(m, p), Less(m, S(0)), Less(p, S(-1))), ZeroQ(m + S(2)*p + S(2)))

cons329 = CustomConstraint(cons_f329)

def cons_f330(f, g, c, a):
    return ZeroQ(a*g**S(2) + c*f**S(2))

cons330 = CustomConstraint(cons_f330)

def cons_f331(p):
    return Less(p, S(-2))

cons331 = CustomConstraint(cons_f331)

def cons_f332(m, p):
    return Or(Less(S(0), -m, p + S(1)), Less(p, -m, S(0)))

cons332 = CustomConstraint(cons_f332)

def cons_f333(n, p):
    return NegativeIntegerQ(n + S(2)*p)

cons333 = CustomConstraint(cons_f333)

def cons_f334(d, g, e, b, f, c):
    return ZeroQ(-b*e*g + c*d*g + c*e*f)

cons334 = CustomConstraint(cons_f334)

def cons_f335(m, n):
    return NonzeroQ(m - n + S(-1))

cons335 = CustomConstraint(cons_f335)

def cons_f336(f, g, e, d):
    return ZeroQ(d*g + e*f)

cons336 = CustomConstraint(cons_f336)

def cons_f337(m, n):
    return ZeroQ(m - n + S(-2))

cons337 = CustomConstraint(cons_f337)

def cons_f338(n, p):
    return RationalQ(n, p)

cons338 = CustomConstraint(cons_f338)

def cons_f339(n, p):
    return Not(And(IntegerQ(n + p), LessEqual(n + p + S(2), S(0))))

cons339 = CustomConstraint(cons_f339)

def cons_f340(n):
    return Not(PositiveIntegerQ(n))

cons340 = CustomConstraint(cons_f340)

def cons_f341(n, p):
    return Not(And(IntegerQ(n + p), Less(n + p + S(2), S(0))))

cons341 = CustomConstraint(cons_f341)

def cons_f342(p, n):
    return Or(IntegerQ(S(2)*p), IntegerQ(n))

cons342 = CustomConstraint(cons_f342)

def cons_f343(m, p):
    return ZeroQ(m + p + S(-1))

cons343 = CustomConstraint(cons_f343)

def cons_f344(n, p, d, g, e, b, f, c):
    return ZeroQ(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))

cons344 = CustomConstraint(cons_f344)

def cons_f345(n, p, d, g, e, f):
    return ZeroQ(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))

cons345 = CustomConstraint(cons_f345)

def cons_f346(n):
    return Not(And(RationalQ(n), Less(n, S(-1))))

cons346 = CustomConstraint(cons_f346)

def cons_f347(p):
    return IntegerQ(p + S(-1)/2)

cons347 = CustomConstraint(cons_f347)

def cons_f348(m, p):
    return Not(And(Less(m, S(0)), Less(p, S(0))))

cons348 = CustomConstraint(cons_f348)

def cons_f349(p):
    return Unequal(p, S(1)/2)

cons349 = CustomConstraint(cons_f349)

def cons_f350(d, g, e, b, f, c, a):
    return ZeroQ(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)

cons350 = CustomConstraint(cons_f350)

def cons_f351(g, e, d, f, c, a):
    return ZeroQ(a*e*g + c*d*f)

cons351 = CustomConstraint(cons_f351)

def cons_f352(m, e, b, d, c):
    return Not(And(Equal(m, S(1)), Or(ZeroQ(d), ZeroQ(-b*e + S(2)*c*d))))

cons352 = CustomConstraint(cons_f352)

def cons_f353(m, d):
    return Not(And(Equal(m, S(1)), ZeroQ(d)))

cons353 = CustomConstraint(cons_f353)

def cons_f354(p, g, e, b, d, f, c, a):
    return ZeroQ(-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))

cons354 = CustomConstraint(cons_f354)

def cons_f355(p, g, e, d, f, c, a):
    return ZeroQ(a*e*g - c*d*f*(S(2)*p + S(3)))

cons355 = CustomConstraint(cons_f355)

def cons_f356(m):
    return Not(RationalQ(m))

cons356 = CustomConstraint(cons_f356)

def cons_f357(p):
    return Not(PositiveIntegerQ(p))

cons357 = CustomConstraint(cons_f357)

def cons_f358(m, p):
    return ZeroQ(m - p)

cons358 = CustomConstraint(cons_f358)

def cons_f359(m, p):
    return Less(m + S(2)*p, S(0))

cons359 = CustomConstraint(cons_f359)

def cons_f360(m, p):
    return Not(NegativeIntegerQ(m + S(2)*p + S(3)))

cons360 = CustomConstraint(cons_f360)

def cons_f361(m, p):
    return Or(And(RationalQ(m), Less(m, S(-1))), Equal(p, S(1)), And(IntegerQ(p), Not(RationalQ(m))))

cons361 = CustomConstraint(cons_f361)

def cons_f362(m, p):
    return Or(IntegerQ(m), IntegerQ(p), IntegersQ(S(2)*m, S(2)*p))

cons362 = CustomConstraint(cons_f362)

def cons_f363(m, p):
    return Or(IntegerQ(p), Not(RationalQ(m)), Inequality(S(-1), LessEqual, m, Less, S(0)))

cons363 = CustomConstraint(cons_f363)

def cons_f364(m, p, g, e, b, d, f, c, a):
    return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, b, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

cons364 = CustomConstraint(cons_f364)

def cons_f365(m, p, g, e, d, f, c, a):
    return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

cons365 = CustomConstraint(cons_f365)

def cons_f366(m, d, g, x, e, f):
    return Not(And(Equal(m, S(1)), SimplerQ(d + e*x, f + g*x)))

cons366 = CustomConstraint(cons_f366)

def cons_f367(m):
    return FractionQ(m)

cons367 = CustomConstraint(cons_f367)

def cons_f368(m, d, g, x, e, f):
    return Not(And(Equal(m, S(1)), SimplerQ(f + g*x, d + e*x)))

cons368 = CustomConstraint(cons_f368)

def cons_f369(m, p):
    return NegativeIntegerQ(m + S(2)*p + S(3))

cons369 = CustomConstraint(cons_f369)

def cons_f370(e, b, d, c, a):
    return ZeroQ(S(4)*c*(a - d) - (b - e)**S(2))

cons370 = CustomConstraint(cons_f370)

def cons_f371(d, g, e, b, f, a):
    return ZeroQ(e*f*(b - e) - S(2)*g*(-a*e + b*d))

cons371 = CustomConstraint(cons_f371)

def cons_f372(d, e, b, a):
    return NonzeroQ(-a*e + b*d)

cons372 = CustomConstraint(cons_f372)

def cons_f373(g, x, f, c, a):
    return FreeQ(List(a, c, f, g), x)

cons373 = CustomConstraint(cons_f373)

def cons_f374(g, e, x, f, c, a):
    return FreeQ(List(a, c, e, f, g), x)

cons374 = CustomConstraint(cons_f374)

def cons_f375(m, n, p):
    return IntegersQ(m, n, p)

cons375 = CustomConstraint(cons_f375)

def cons_f376(n, p):
    return IntegersQ(n, p)

cons376 = CustomConstraint(cons_f376)

def cons_f377(m, f, d):
    return Or(IntegerQ(m), And(PositiveQ(d), PositiveQ(f)))

cons377 = CustomConstraint(cons_f377)

def cons_f378(m, p, n):
    return Or(IntegerQ(p), IntegersQ(m, n))

cons378 = CustomConstraint(cons_f378)

def cons_f379(m, p, g, e, x, f, c, a):
    return FreeQ(List(a, c, e, f, g, m, p), x)

cons379 = CustomConstraint(cons_f379)

def cons_f380(m, n, p, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p), x)

cons380 = CustomConstraint(cons_f380)

def cons_f381(m, n, p, g, e, x, d, f, c, a):
    return FreeQ(List(a, c, d, e, f, g, m, n, p), x)

cons381 = CustomConstraint(cons_f381)

def cons_f382(d, f, c, a):
    return ZeroQ(-a*f + c*d)

cons382 = CustomConstraint(cons_f382)

def cons_f383(d, e, b, a):
    return ZeroQ(-a*e + b*d)

cons383 = CustomConstraint(cons_f383)

def cons_f384(f, p, c):
    return Or(IntegerQ(p), PositiveQ(c/f))

cons384 = CustomConstraint(cons_f384)

def cons_f385(q, x, e, b, d, f, c, a):
    return Or(Not(IntegerQ(q)), LessEqual(LeafCount(d + e*x + f*x**S(2)), LeafCount(a + b*x + c*x**S(2))))

cons385 = CustomConstraint(cons_f385)

def cons_f386(q):
    return Not(IntegerQ(q))

cons386 = CustomConstraint(cons_f386)

def cons_f387(f, c):
    return Not(PositiveQ(c/f))

cons387 = CustomConstraint(cons_f387)

def cons_f388(q, d, e, b, f, c, a):
    return ZeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

cons388 = CustomConstraint(cons_f388)

def cons_f389(q):
    return NonzeroQ(q + S(1))

cons389 = CustomConstraint(cons_f389)

def cons_f390(q):
    return NonzeroQ(S(2)*q + S(3))

cons390 = CustomConstraint(cons_f390)

def cons_f391(q, d, e, f, c, a):
    return ZeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

cons391 = CustomConstraint(cons_f391)

def cons_f392(q, d, f, c, a):
    return ZeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

cons392 = CustomConstraint(cons_f392)

def cons_f393(q):
    return PositiveIntegerQ(q + S(2))

cons393 = CustomConstraint(cons_f393)

def cons_f394(d, f, e):
    return NonzeroQ(-S(4)*d*f + e**S(2))

cons394 = CustomConstraint(cons_f394)

def cons_f395(q):
    return RationalQ(q)

cons395 = CustomConstraint(cons_f395)

def cons_f396(q):
    return Less(q, S(-1))

cons396 = CustomConstraint(cons_f396)

def cons_f397(q, d, e, b, f, c, a):
    return NonzeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

cons397 = CustomConstraint(cons_f397)

def cons_f398(q, d, e, f, c, a):
    return NonzeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

cons398 = CustomConstraint(cons_f398)

def cons_f399(q, d, f, c, a):
    return NonzeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

cons399 = CustomConstraint(cons_f399)

def cons_f400(q):
    return Not(PositiveIntegerQ(q))

cons400 = CustomConstraint(cons_f400)

def cons_f401(q):
    return Not(And(RationalQ(q), LessEqual(q, S(-1))))

cons401 = CustomConstraint(cons_f401)

def cons_f402(p, q):
    return RationalQ(p, q)

cons402 = CustomConstraint(cons_f402)

def cons_f403(q):
    return Greater(q, S(0))

cons403 = CustomConstraint(cons_f403)

def cons_f404(e, b, d, f, c, a):
    return NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))

cons404 = CustomConstraint(cons_f404)

def cons_f405(p, q):
    return Not(And(Not(IntegerQ(p)), IntegerQ(q), Less(q, S(-1))))

cons405 = CustomConstraint(cons_f405)

def cons_f406(b, d, f, c, a):
    return NonzeroQ(b**S(2)*d*f + (-a*f + c*d)**S(2))

cons406 = CustomConstraint(cons_f406)

def cons_f407(e, d, f, c, a):
    return NonzeroQ(a*c*e**S(2) + (-a*f + c*d)**S(2))

cons407 = CustomConstraint(cons_f407)

def cons_f408(p, q):
    return NonzeroQ(p + q)

cons408 = CustomConstraint(cons_f408)

def cons_f409(p, q):
    return NonzeroQ(S(2)*p + S(2)*q + S(1))

cons409 = CustomConstraint(cons_f409)

def cons_f410(f, e, c, b):
    return ZeroQ(-b*f + c*e)

cons410 = CustomConstraint(cons_f410)

def cons_f411(f, e, c, b):
    return NonzeroQ(-b*f + c*e)

cons411 = CustomConstraint(cons_f411)

def cons_f412(c, a):
    return PosQ(-a*c)

cons412 = CustomConstraint(cons_f412)

def cons_f413(c, b, a):
    return NegQ(-S(4)*a*c + b**S(2))

cons413 = CustomConstraint(cons_f413)

def cons_f414(c, a):
    return NegQ(-a*c)

cons414 = CustomConstraint(cons_f414)

def cons_f415(p, q, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, p, q), x)

cons415 = CustomConstraint(cons_f415)

def cons_f416(p, q, e, x, d, f, c, a):
    return FreeQ(List(a, c, d, e, f, p, q), x)

cons416 = CustomConstraint(cons_f416)

def cons_f417(h, g, b, c, a):
    return ZeroQ(a*h**S(2) - b*g*h + c*g**S(2))

cons417 = CustomConstraint(cons_f417)

def cons_f418(h, g, e, d, f, c, a):
    return ZeroQ(a**S(2)*f*h**S(2) - a*c*e*g*h + c**S(2)*d*g**S(2))

cons418 = CustomConstraint(cons_f418)

def cons_f419(h, g, c, a):
    return ZeroQ(a*h**S(2) + c*g**S(2))

cons419 = CustomConstraint(cons_f419)

def cons_f420(h, g, d, f, c, a):
    return ZeroQ(a**S(2)*f*h**S(2) + c**S(2)*d*g**S(2))

cons420 = CustomConstraint(cons_f420)

def cons_f421(e, b, f, c, a):
    return ZeroQ(a*f**S(2) - b*e*f + c*e**S(2))

cons421 = CustomConstraint(cons_f421)

def cons_f422(f, e, c, a):
    return ZeroQ(a*f**S(2) + c*e**S(2))

cons422 = CustomConstraint(cons_f422)

def cons_f423(m, h, p, g, e, b, f, c):
    return ZeroQ(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

cons423 = CustomConstraint(cons_f423)

def cons_f424(m, h, p, d, g, b, f, c, a):
    return ZeroQ(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

cons424 = CustomConstraint(cons_f424)

def cons_f425(m, h, p, g, e, f, c):
    return ZeroQ(c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

cons425 = CustomConstraint(cons_f425)

def cons_f426(m, h, p, d, f, c, a):
    return ZeroQ(h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

cons426 = CustomConstraint(cons_f426)

def cons_f427(m, h, p, g, b, f, c):
    return ZeroQ(b*f*h*(m + p + S(2)) + S(2)*c*f*g*(p + S(1)))

cons427 = CustomConstraint(cons_f427)

def cons_f428(m, p):
    return Or(IntegersQ(m, p), PositiveIntegerQ(p))

cons428 = CustomConstraint(cons_f428)

def cons_f429(h, g, b, c, a):
    return NonzeroQ(a*h**S(2) - b*g*h + c*g**S(2))

cons429 = CustomConstraint(cons_f429)

def cons_f430(h, g, c, a):
    return NonzeroQ(a*h**S(2) + c*g**S(2))

cons430 = CustomConstraint(cons_f430)

def cons_f431(h, g, b, c, a):
    return NonzeroQ(c*g**S(2) - h*(-a*h + b*g))

cons431 = CustomConstraint(cons_f431)

def cons_f432(p, q):
    return Or(Greater(p, S(0)), Greater(q, S(0)))

cons432 = CustomConstraint(cons_f432)

def cons_f433(p, q):
    return NonzeroQ(p + q + S(1))

cons433 = CustomConstraint(cons_f433)

def cons_f434(c, a):
    return PositiveQ(a*c)

cons434 = CustomConstraint(cons_f434)

def cons_f435(c, a):
    return Not(PositiveQ(a*c))

cons435 = CustomConstraint(cons_f435)

def cons_f436(h, g, e, f):
    return ZeroQ(e*h - S(2)*f*g)

cons436 = CustomConstraint(cons_f436)

def cons_f437(h, g, e, f):
    return NonzeroQ(e*h - S(2)*f*g)

cons437 = CustomConstraint(cons_f437)

def cons_f438(h, g, d, e):
    return ZeroQ(S(2)*d*h - e*g)

cons438 = CustomConstraint(cons_f438)

def cons_f439(h, g, d, e):
    return NonzeroQ(S(2)*d*h - e*g)

cons439 = CustomConstraint(cons_f439)

def cons_f440(h, g, e, b, d, f, c, a):
    return ZeroQ(g**S(2)*(-b*f + c*e) - S(2)*g*h*(-a*f + c*d) + h**S(2)*(-a*e + b*d))

cons440 = CustomConstraint(cons_f440)

def cons_f441(h, g, e, d, f, c, a):
    return ZeroQ(a*e*h**S(2) - c*e*g**S(2) + S(2)*g*h*(-a*f + c*d))

cons441 = CustomConstraint(cons_f441)

def cons_f442(h, g, b, d, f, c, a):
    return ZeroQ(b*d*h**S(2) - b*f*g**S(2) - S(2)*g*h*(-a*f + c*d))

cons442 = CustomConstraint(cons_f442)

def cons_f443(b, d, f, c, a):
    return ZeroQ(c**S(2)*d - f*(-S(3)*a*c + b**S(2)))

cons443 = CustomConstraint(cons_f443)

def cons_f444(h, g, b, c, a):
    return ZeroQ(S(9)*a*c*h**S(2) - S(2)*b**S(2)*h**S(2) - b*c*g*h + c**S(2)*g**S(2))

cons444 = CustomConstraint(cons_f444)

def cons_f445(h, g, c, b):
    return PositiveQ(-S(9)*c*h**S(2)/(-b*h + S(2)*c*g)**S(2))

cons445 = CustomConstraint(cons_f445)

def cons_f446(d, f, c, a):
    return ZeroQ(S(3)*a*f + c*d)

cons446 = CustomConstraint(cons_f446)

def cons_f447(h, g, c, a):
    return ZeroQ(S(9)*a*h**S(2) + c*g**S(2))

cons447 = CustomConstraint(cons_f447)

def cons_f448(a):
    return Not(PositiveQ(a))

cons448 = CustomConstraint(cons_f448)

def cons_f449(h, p, q, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, h, p, q), x)

cons449 = CustomConstraint(cons_f449)

def cons_f450(h, p, q, g, e, x, d, f, c, a):
    return FreeQ(List(a, c, d, e, f, g, h, p, q), x)

cons450 = CustomConstraint(cons_f450)

def cons_f451(z, x):
    return LinearQ(z, x)

cons451 = CustomConstraint(cons_f451)

def cons_f452(v, x, u):
    return QuadraticQ(List(u, v), x)

cons452 = CustomConstraint(cons_f452)

def cons_f453(v, z, x, u):
    return Not(And(LinearMatchQ(z, x), QuadraticMatchQ(List(u, v), x)))

cons453 = CustomConstraint(cons_f453)

def cons_f454(p, q):
    return NonzeroQ(S(2)*p + S(2)*q + S(3))

cons454 = CustomConstraint(cons_f454)

def cons_f455(C, p, q, B, A, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, A, B, C, p, q), x)

cons455 = CustomConstraint(cons_f455)

def cons_f456(C, p, q, A, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, A, C, p, q), x)

cons456 = CustomConstraint(cons_f456)

def cons_f457(C, p, q, B, A, e, x, d, f, c, a):
    return FreeQ(List(a, c, d, e, f, A, B, C, p, q), x)

cons457 = CustomConstraint(cons_f457)

def cons_f458(C, p, q, A, e, x, d, f, c, a):
    return FreeQ(List(a, c, d, e, f, A, C, p, q), x)

cons458 = CustomConstraint(cons_f458)

def cons_f459(n, x, p, b):
    return FreeQ(List(b, n, p), x)

cons459 = CustomConstraint(cons_f459)

def cons_f460(n, p):
    return ZeroQ(p + S(1) + S(1)/n)

cons460 = CustomConstraint(cons_f460)

def cons_f461(n, p):
    return NegativeIntegerQ(p + S(1) + S(1)/n)

cons461 = CustomConstraint(cons_f461)

def cons_f462(n):
    return NonzeroQ(S(3)*n + S(1))

cons462 = CustomConstraint(cons_f462)

def cons_f463(n):
    return Less(n, S(0))

cons463 = CustomConstraint(cons_f463)

def cons_f464(n, p):
    return PositiveIntegerQ(n, p)

cons464 = CustomConstraint(cons_f464)

def cons_f465(p, n):
    return Or(IntegerQ(S(2)*p), And(Equal(n, S(2)), IntegerQ(S(4)*p)), And(Equal(n, S(2)), IntegerQ(S(3)*p)), Less(Denominator(p + S(1)/n), Denominator(p)))

cons465 = CustomConstraint(cons_f465)

def cons_f466(b, a):
    return PosQ(b/a)

cons466 = CustomConstraint(cons_f466)

def cons_f467(n):
    return PositiveIntegerQ(n/S(2) + S(-3)/2)

cons467 = CustomConstraint(cons_f467)

def cons_f468(b, a):
    return PosQ(a/b)

cons468 = CustomConstraint(cons_f468)

def cons_f469(b, a):
    return NegQ(a/b)

cons469 = CustomConstraint(cons_f469)

def cons_f470(b, a):
    return Or(PositiveQ(a), PositiveQ(b))

cons470 = CustomConstraint(cons_f470)

def cons_f471(b, a):
    return Or(NegativeQ(a), NegativeQ(b))

cons471 = CustomConstraint(cons_f471)

def cons_f472(b, a):
    return Or(PositiveQ(a), NegativeQ(b))

cons472 = CustomConstraint(cons_f472)

def cons_f473(b, a):
    return Or(NegativeQ(a), PositiveQ(b))

cons473 = CustomConstraint(cons_f473)

def cons_f474(n):
    return PositiveIntegerQ(n/S(4) + S(-1)/2)

cons474 = CustomConstraint(cons_f474)

def cons_f475(b, a):
    return Or(PositiveQ(a/b), And(PosQ(a/b), AtomQ(SplitProduct(SumBaseQ, a)), AtomQ(SplitProduct(SumBaseQ, b))))

cons475 = CustomConstraint(cons_f475)

def cons_f476(b, a):
    return Not(PositiveQ(a/b))

cons476 = CustomConstraint(cons_f476)

def cons_f477(n):
    return PositiveIntegerQ(n/S(4) + S(-1))

cons477 = CustomConstraint(cons_f477)

def cons_f478(b, a):
    return PositiveQ(a/b)

cons478 = CustomConstraint(cons_f478)

def cons_f479(b):
    return PosQ(b)

cons479 = CustomConstraint(cons_f479)

def cons_f480(b):
    return NegQ(b)

cons480 = CustomConstraint(cons_f480)

def cons_f481(a):
    return PosQ(a)

cons481 = CustomConstraint(cons_f481)

def cons_f482(a):
    return NegQ(a)

cons482 = CustomConstraint(cons_f482)

def cons_f483(b, a):
    return NegQ(b/a)

cons483 = CustomConstraint(cons_f483)

def cons_f484(a):
    return NegativeQ(a)

cons484 = CustomConstraint(cons_f484)

def cons_f485(p):
    return Less(S(-1), p, S(0))

cons485 = CustomConstraint(cons_f485)

def cons_f486(p):
    return Unequal(p, S(-1)/2)

cons486 = CustomConstraint(cons_f486)

def cons_f487(p, n):
    return IntegerQ(p + S(1)/n)

cons487 = CustomConstraint(cons_f487)

def cons_f488(p, n):
    return Less(Denominator(p + S(1)/n), Denominator(p))

cons488 = CustomConstraint(cons_f488)

def cons_f489(n):
    return FractionQ(n)

cons489 = CustomConstraint(cons_f489)

def cons_f490(n):
    return Not(IntegerQ(S(1)/n))

cons490 = CustomConstraint(cons_f490)

def cons_f491(n, p):
    return Not(NegativeIntegerQ(p + S(1)/n))

cons491 = CustomConstraint(cons_f491)

def cons_f492(p, a):
    return Or(IntegerQ(p), PositiveQ(a))

cons492 = CustomConstraint(cons_f492)

def cons_f493(p, a):
    return Not(Or(IntegerQ(p), PositiveQ(a)))

cons493 = CustomConstraint(cons_f493)

def cons_f494(p, a2, a1):
    return Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2)))

cons494 = CustomConstraint(cons_f494)

def cons_f495(n):
    return PositiveIntegerQ(S(2)*n)

cons495 = CustomConstraint(cons_f495)

def cons_f496(p, n):
    return Or(IntegerQ(S(2)*p), Less(Denominator(p + S(1)/n), Denominator(p)))

cons496 = CustomConstraint(cons_f496)

def cons_f497(n):
    return NegativeIntegerQ(S(2)*n)

cons497 = CustomConstraint(cons_f497)

def cons_f498(n):
    return FractionQ(S(2)*n)

cons498 = CustomConstraint(cons_f498)

def cons_f499(m, c):
    return Or(IntegerQ(m), PositiveQ(c))

cons499 = CustomConstraint(cons_f499)

def cons_f500(m, n):
    return IntegerQ((m + S(1))/n)

cons500 = CustomConstraint(cons_f500)

def cons_f501(m, n):
    return Not(IntegerQ((m + S(1))/n))

cons501 = CustomConstraint(cons_f501)

def cons_f502(n):
    return NegQ(n)

cons502 = CustomConstraint(cons_f502)

def cons_f503(m, n, p):
    return ZeroQ(p + S(1) + (m + S(1))/n)

cons503 = CustomConstraint(cons_f503)

def cons_f504(m, n, p):
    return ZeroQ(p + S(1) + (m + S(1))/(S(2)*n))

cons504 = CustomConstraint(cons_f504)

def cons_f505(m, n):
    return IntegerQ((m + S(1))/(S(2)*n))

cons505 = CustomConstraint(cons_f505)

def cons_f506(m, n, p):
    return NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n)

cons506 = CustomConstraint(cons_f506)

def cons_f507(m, n, p):
    return NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n))

cons507 = CustomConstraint(cons_f507)

def cons_f508(m, n, p):
    return Not(NegativeIntegerQ((m + n*p + n + S(1))/n))

cons508 = CustomConstraint(cons_f508)

def cons_f509(m, n, p, x, b, c, a):
    return IntBinomialQ(a, b, c, n, m, p, x)

cons509 = CustomConstraint(cons_f509)

def cons_f510(m, n, p):
    return NonzeroQ(m + S(2)*n*p + S(1))

cons510 = CustomConstraint(cons_f510)

def cons_f511(m, n, a1, p, b2, a2, b1, x, c):
    return IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)

cons511 = CustomConstraint(cons_f511)

def cons_f512(m, n, p):
    return NonzeroQ(m + n*p + S(1))

cons512 = CustomConstraint(cons_f512)

def cons_f513(m):
    return PositiveIntegerQ(m/S(4) + S(-1)/2)

cons513 = CustomConstraint(cons_f513)

def cons_f514(m):
    return NegativeIntegerQ(m/S(4) + S(-1)/2)

cons514 = CustomConstraint(cons_f514)

def cons_f515(m):
    return IntegerQ(S(2)*m)

cons515 = CustomConstraint(cons_f515)

def cons_f516(m):
    return Greater(m, S(3)/2)

cons516 = CustomConstraint(cons_f516)

def cons_f517(m, n):
    return Greater(m + S(1), n)

cons517 = CustomConstraint(cons_f517)

def cons_f518(m, n, p):
    return Not(NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n))

cons518 = CustomConstraint(cons_f518)

def cons_f519(m, n):
    return Greater(m + S(1), S(2)*n)

cons519 = CustomConstraint(cons_f519)

def cons_f520(m, n, p):
    return Not(NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)))

cons520 = CustomConstraint(cons_f520)

def cons_f521(n):
    return PositiveIntegerQ(n/S(2) + S(-1)/2)

cons521 = CustomConstraint(cons_f521)

def cons_f522(m, n):
    return Less(m, n + S(-1))

cons522 = CustomConstraint(cons_f522)

def cons_f523(m, n):
    return PositiveIntegerQ(m, n/S(2) + S(-1)/2)

cons523 = CustomConstraint(cons_f523)

def cons_f524(m, n):
    return PositiveIntegerQ(m, n/S(4) + S(-1)/2)

cons524 = CustomConstraint(cons_f524)

def cons_f525(m, n):
    return PositiveIntegerQ(m, n/S(4))

cons525 = CustomConstraint(cons_f525)

def cons_f526(m, n):
    return Less(m, n/S(2))

cons526 = CustomConstraint(cons_f526)

def cons_f527(m, n):
    return Inequality(n/S(2), LessEqual, m, Less, n)

cons527 = CustomConstraint(cons_f527)

def cons_f528(m, n):
    return PositiveIntegerQ(m, n)

cons528 = CustomConstraint(cons_f528)

def cons_f529(m, n):
    return Greater(m, S(2)*n + S(-1))

cons529 = CustomConstraint(cons_f529)

def cons_f530(m, n):
    return Greater(m, n + S(-1))

cons530 = CustomConstraint(cons_f530)

def cons_f531(m, n):
    return SumSimplerQ(m, -n)

cons531 = CustomConstraint(cons_f531)

def cons_f532(m, n, p):
    return NegativeIntegerQ((m + n*p + S(1))/n)

cons532 = CustomConstraint(cons_f532)

def cons_f533(m, n):
    return SumSimplerQ(m, -S(2)*n)

cons533 = CustomConstraint(cons_f533)

def cons_f534(m, n, p):
    return NegativeIntegerQ((m + S(2)*n*p + S(1))/(S(2)*n))

cons534 = CustomConstraint(cons_f534)

def cons_f535(m, n):
    return SumSimplerQ(m, n)

cons535 = CustomConstraint(cons_f535)

def cons_f536(m, n):
    return SumSimplerQ(m, S(2)*n)

cons536 = CustomConstraint(cons_f536)

def cons_f537(m, p, n):
    return IntegersQ(m, p + (m + S(1))/n)

cons537 = CustomConstraint(cons_f537)

def cons_f538(m, p, n):
    return IntegersQ(m, p + (m + S(1))/(S(2)*n))

cons538 = CustomConstraint(cons_f538)

def cons_f539(m, p, n):
    return Less(Denominator(p + (m + S(1))/n), Denominator(p))

cons539 = CustomConstraint(cons_f539)

def cons_f540(m, p, n):
    return Less(Denominator(p + (m + S(1))/(S(2)*n)), Denominator(p))

cons540 = CustomConstraint(cons_f540)

def cons_f541(m, n):
    return IntegerQ(n/(m + S(1)))

cons541 = CustomConstraint(cons_f541)

def cons_f542(m, n):
    return IntegerQ(S(2)*n/(m + S(1)))

cons542 = CustomConstraint(cons_f542)

def cons_f543(n):
    return Not(IntegerQ(S(2)*n))

cons543 = CustomConstraint(cons_f543)

def cons_f544(m, n, p):
    return ZeroQ(p + (m + S(1))/n)

cons544 = CustomConstraint(cons_f544)

def cons_f545(m, n, p):
    return ZeroQ(p + (m + S(1))/(S(2)*n))

cons545 = CustomConstraint(cons_f545)

def cons_f546(m, p, n):
    return IntegerQ(p + (m + S(1))/n)

cons546 = CustomConstraint(cons_f546)

def cons_f547(m, p, n):
    return IntegerQ(p + (m + S(1))/(S(2)*n))

cons547 = CustomConstraint(cons_f547)

def cons_f548(m, n):
    return FractionQ((m + S(1))/n)

cons548 = CustomConstraint(cons_f548)

def cons_f549(m, n):
    return Or(SumSimplerQ(m, n), SumSimplerQ(m, -n))

cons549 = CustomConstraint(cons_f549)

def cons_f550(p, a):
    return Or(NegativeIntegerQ(p), PositiveQ(a))

cons550 = CustomConstraint(cons_f550)

def cons_f551(p, a):
    return Not(Or(NegativeIntegerQ(p), PositiveQ(a)))

cons551 = CustomConstraint(cons_f551)

def cons_f552(v, x):
    return LinearQ(v, x)

cons552 = CustomConstraint(cons_f552)

def cons_f553(v, x):
    return NonzeroQ(v - x)

cons553 = CustomConstraint(cons_f553)

def cons_f554(v, x, u):
    return LinearPairQ(u, v, x)

cons554 = CustomConstraint(cons_f554)

def cons_f555(p, q):
    return PositiveIntegerQ(p, q)

cons555 = CustomConstraint(cons_f555)

def cons_f556(n, p):
    return ZeroQ(n*p + S(1))

cons556 = CustomConstraint(cons_f556)

def cons_f557(n, q, p):
    return ZeroQ(n*(p + q + S(1)) + S(1))

cons557 = CustomConstraint(cons_f557)

def cons_f558(n, q, p):
    return ZeroQ(n*(p + q + S(2)) + S(1))

cons558 = CustomConstraint(cons_f558)

def cons_f559(p, q, b, d, c, a):
    return ZeroQ(a*d*(p + S(1)) + b*c*(q + S(1)))

cons559 = CustomConstraint(cons_f559)

def cons_f560(p, q):
    return Or(And(RationalQ(p), Less(p, S(-1))), Not(And(RationalQ(q), Less(q, S(-1)))))

cons560 = CustomConstraint(cons_f560)

def cons_f561(n, p, b, d, c, a):
    return ZeroQ(a*d - b*c*(n*(p + S(1)) + S(1)))

cons561 = CustomConstraint(cons_f561)

def cons_f562(p, n):
    return Or(And(RationalQ(p), Less(p, S(-1))), NegativeIntegerQ(p + S(1)/n))

cons562 = CustomConstraint(cons_f562)

def cons_f563(n, p):
    return NonzeroQ(n*(p + S(1)) + S(1))

cons563 = CustomConstraint(cons_f563)

def cons_f564(q):
    return NegativeIntegerQ(q)

cons564 = CustomConstraint(cons_f564)

def cons_f565(p, q):
    return GreaterEqual(p, -q)

cons565 = CustomConstraint(cons_f565)

def cons_f566(d, c, b, a):
    return ZeroQ(S(3)*a*d + b*c)

cons566 = CustomConstraint(cons_f566)

def cons_f567(p):
    return Or(Equal(p, S(1)/2), Equal(Denominator(p), S(4)))

cons567 = CustomConstraint(cons_f567)

def cons_f568(p):
    return Equal(Denominator(p), S(4))

cons568 = CustomConstraint(cons_f568)

def cons_f569(p):
    return Or(Equal(p, S(-5)/4), Equal(p, S(-7)/4))

cons569 = CustomConstraint(cons_f569)

def cons_f570(b, a):
    return PosQ(a*b)

cons570 = CustomConstraint(cons_f570)

def cons_f571(b, a):
    return NegQ(a*b)

cons571 = CustomConstraint(cons_f571)

def cons_f572(p):
    return Or(Equal(p, S(3)/4), Equal(p, S(5)/4))

cons572 = CustomConstraint(cons_f572)

def cons_f573(d, c):
    return PosQ(d/c)

cons573 = CustomConstraint(cons_f573)

def cons_f574(q):
    return Less(S(0), q, S(1))

cons574 = CustomConstraint(cons_f574)

def cons_f575(n, p, q, x, b, d, c, a):
    return IntBinomialQ(a, b, c, d, n, p, q, x)

cons575 = CustomConstraint(cons_f575)

def cons_f576(q):
    return Greater(q, S(1))

cons576 = CustomConstraint(cons_f576)

def cons_f577(p, q):
    return Greater(p + q, S(0))

cons577 = CustomConstraint(cons_f577)

def cons_f578(n, q, p):
    return NonzeroQ(n*(p + q) + S(1))

cons578 = CustomConstraint(cons_f578)

def cons_f579(p):
    return Not(And(IntegerQ(p), Greater(p, S(1))))

cons579 = CustomConstraint(cons_f579)

def cons_f580(d, c, b, a):
    return Not(SimplerSqrtQ(b/a, d/c))

cons580 = CustomConstraint(cons_f580)

def cons_f581(d, c):
    return NegQ(d/c)

cons581 = CustomConstraint(cons_f581)

def cons_f582(d, c, b, a):
    return Not(And(NegQ(b/a), SimplerSqrtQ(-b/a, -d/c)))

cons582 = CustomConstraint(cons_f582)

def cons_f583(d, c, b, a):
    return PositiveQ(a - b*c/d)

cons583 = CustomConstraint(cons_f583)

def cons_f584(n):
    return NonzeroQ(n + S(1))

cons584 = CustomConstraint(cons_f584)

def cons_f585(n, mn):
    return EqQ(mn, -n)

cons585 = CustomConstraint(cons_f585)

def cons_f586(q):
    return IntegerQ(q)

cons586 = CustomConstraint(cons_f586)

def cons_f587(n, p):
    return Or(PosQ(n), Not(IntegerQ(p)))

cons587 = CustomConstraint(cons_f587)

def cons_f588(v, x, u):
    return PseudoBinomialPairQ(u, v, x)

cons588 = CustomConstraint(cons_f588)

def cons_f589(m, p):
    return IntegersQ(p, m/p)

cons589 = CustomConstraint(cons_f589)

def cons_f590(m, p, v, x, u):
    return PseudoBinomialPairQ(u*x**(m/p), v, x)

cons590 = CustomConstraint(cons_f590)

def cons_f591(m, e):
    return Or(IntegerQ(m), PositiveQ(e))

cons591 = CustomConstraint(cons_f591)

def cons_f592(m, n, p, b, d, c, a):
    return ZeroQ(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))

cons592 = CustomConstraint(cons_f592)

def cons_f593(n, non2):
    return ZeroQ(-n/S(2) + non2)

cons593 = CustomConstraint(cons_f593)

def cons_f594(m, n, a1, p, b2, a2, b1, d, c):
    return ZeroQ(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))

cons594 = CustomConstraint(cons_f594)

def cons_f595(m, n, p):
    return ZeroQ(m + n*(p + S(1)) + S(1))

cons595 = CustomConstraint(cons_f595)

def cons_f596(n, e):
    return Or(IntegerQ(n), PositiveQ(e))

cons596 = CustomConstraint(cons_f596)

def cons_f597(m, n):
    return Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1))))

cons597 = CustomConstraint(cons_f597)

def cons_f598(p):
    return Not(And(IntegerQ(p), Less(p, S(-1))))

cons598 = CustomConstraint(cons_f598)

def cons_f599(m):
    return PositiveIntegerQ(m/S(2))

cons599 = CustomConstraint(cons_f599)

def cons_f600(m, p):
    return Or(IntegerQ(p), Equal(m + S(2)*p + S(1), S(0)))

cons600 = CustomConstraint(cons_f600)

def cons_f601(m):
    return NegativeIntegerQ(m/S(2))

cons601 = CustomConstraint(cons_f601)

def cons_f602(m, p, n):
    return Or(IntegerQ(p), Not(RationalQ(m)), And(PositiveIntegerQ(n), NegativeIntegerQ(p + S(1)/2), LessEqual(S(-1), m, -n*(p + S(1)))))

cons602 = CustomConstraint(cons_f602)

def cons_f603(m, n, p):
    return NonzeroQ(m + n*(p + S(1)) + S(1))

cons603 = CustomConstraint(cons_f603)

def cons_f604(m):
    return Or(IntegerQ(m), PositiveIntegerQ(S(2)*m + S(2)), Not(RationalQ(m)))

cons604 = CustomConstraint(cons_f604)

def cons_f605(m, n, p):
    return NonzeroQ(m + n*(p + S(2)) + S(1))

cons605 = CustomConstraint(cons_f605)

def cons_f606(m, p, q):
    return RationalQ(m, p, q)

cons606 = CustomConstraint(cons_f606)

def cons_f607(m, n):
    return Greater(m - n + S(1), S(0))

cons607 = CustomConstraint(cons_f607)

def cons_f608(m, n, p, q, e, x, b, d, c, a):
    return IntBinomialQ(a, b, c, d, e, m, n, p, q, x)

cons608 = CustomConstraint(cons_f608)

def cons_f609(m, n):
    return Greater(m - n + S(1), n)

cons609 = CustomConstraint(cons_f609)

def cons_f610(m, n):
    return Inequality(n, GreaterEqual, m - n + S(1), Greater, S(0))

cons610 = CustomConstraint(cons_f610)

def cons_f611(m, q):
    return RationalQ(m, q)

cons611 = CustomConstraint(cons_f611)

def cons_f612(m, n):
    return LessEqual(n, m, S(2)*n + S(-1))

cons612 = CustomConstraint(cons_f612)

def cons_f613(m, n):
    return IntegersQ(m/S(2), n/S(2))

cons613 = CustomConstraint(cons_f613)

def cons_f614(m, n):
    return Less(S(0), m - n + S(1), n)

cons614 = CustomConstraint(cons_f614)

def cons_f615(n):
    return LessEqual(n, S(4))

cons615 = CustomConstraint(cons_f615)

def cons_f616(d, c, b, a):
    return ZeroQ(-a*d + S(4)*b*c)

cons616 = CustomConstraint(cons_f616)

def cons_f617(m):
    return PositiveIntegerQ(m/S(3) + S(-1)/3)

cons617 = CustomConstraint(cons_f617)

def cons_f618(m):
    return NegativeIntegerQ(m/S(3) + S(-1)/3)

cons618 = CustomConstraint(cons_f618)

def cons_f619(m):
    return IntegerQ(m/S(3) + S(-1)/3)

cons619 = CustomConstraint(cons_f619)

def cons_f620(n):
    return Or(EqQ(n, S(2)), EqQ(n, S(4)))

cons620 = CustomConstraint(cons_f620)

def cons_f621(n, b, d, c, a):
    return Not(And(EqQ(n, S(2)), SimplerSqrtQ(-b/a, -d/c)))

cons621 = CustomConstraint(cons_f621)

def cons_f622(m, p, q, n):
    return IntegersQ(p + (m + S(1))/n, q)

cons622 = CustomConstraint(cons_f622)

def cons_f623(m, n):
    return Or(ZeroQ(m - n), ZeroQ(m - S(2)*n + S(1)))

cons623 = CustomConstraint(cons_f623)

def cons_f624(m, p, q):
    return IntegersQ(m, p, q)

cons624 = CustomConstraint(cons_f624)

def cons_f625(p):
    return GreaterEqual(p, S(-2))

cons625 = CustomConstraint(cons_f625)

def cons_f626(m, q):
    return Or(GreaterEqual(q, S(-2)), And(Equal(q, S(-3)), IntegerQ(m/S(2) + S(-1)/2)))

cons626 = CustomConstraint(cons_f626)

def cons_f627(m, n):
    return NonzeroQ(m - n + S(1))

cons627 = CustomConstraint(cons_f627)

def cons_f628(p, r, q):
    return PositiveIntegerQ(p, q, r)

cons628 = CustomConstraint(cons_f628)

def cons_f629(n, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, n), x)

cons629 = CustomConstraint(cons_f629)

def cons_f630(n, b, d, c, a):
    return Not(And(ZeroQ(n + S(-2)), Or(And(PosQ(b/a), PosQ(d/c)), And(NegQ(b/a), Or(PosQ(d/c), And(PositiveQ(a), Or(Not(PositiveQ(c)), SimplerSqrtQ(-b/a, -d/c))))))))

cons630 = CustomConstraint(cons_f630)

def cons_f631(n, q, p):
    return NonzeroQ(n*(p + q + S(1)) + S(1))

cons631 = CustomConstraint(cons_f631)

def cons_f632(n, p, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, p, n), x)

cons632 = CustomConstraint(cons_f632)

def cons_f633(n, p, q, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, n, p, q), x)

cons633 = CustomConstraint(cons_f633)

def cons_f634(d, c):
    return PositiveQ(d/c)

cons634 = CustomConstraint(cons_f634)

def cons_f635(f, e):
    return PositiveQ(f/e)

cons635 = CustomConstraint(cons_f635)

def cons_f636(d, f, e, c):
    return Not(SimplerSqrtQ(d/c, f/e))

cons636 = CustomConstraint(cons_f636)

def cons_f637(f, e, c, d):
    return Not(SimplerSqrtQ(-f/e, -d/c))

cons637 = CustomConstraint(cons_f637)

def cons_f638(f, e):
    return PosQ(f/e)

cons638 = CustomConstraint(cons_f638)

def cons_f639(f, e, c, d):
    return Not(And(NegQ(f/e), SimplerSqrtQ(-f/e, -d/c)))

cons639 = CustomConstraint(cons_f639)

def cons_f640(q, r):
    return RationalQ(q, r)

cons640 = CustomConstraint(cons_f640)

def cons_f641(r):
    return Greater(r, S(1))

cons641 = CustomConstraint(cons_f641)

def cons_f642(q):
    return LessEqual(q, S(-1))

cons642 = CustomConstraint(cons_f642)

def cons_f643(d, f, e, c):
    return PosQ((-c*f + d*e)/c)

cons643 = CustomConstraint(cons_f643)

def cons_f644(d, f, e, c):
    return NegQ((-c*f + d*e)/c)

cons644 = CustomConstraint(cons_f644)

def cons_f645(n, p, q, e, x, b, d, f, r, c, a):
    return FreeQ(List(a, b, c, d, e, f, n, p, q, r), x)

cons645 = CustomConstraint(cons_f645)

def cons_f646(v, u):
    return ZeroQ(u - v)

cons646 = CustomConstraint(cons_f646)

def cons_f647(w, u):
    return ZeroQ(u - w)

cons647 = CustomConstraint(cons_f647)

def cons_f648(r):
    return IntegerQ(r)

cons648 = CustomConstraint(cons_f648)

def cons_f649(n2, n):
    return ZeroQ(-n/S(2) + n2)

cons649 = CustomConstraint(cons_f649)

def cons_f650(e1, f1, e2, f2):
    return ZeroQ(e1*f2 + e2*f1)

cons650 = CustomConstraint(cons_f650)

def cons_f651(e2, r, e1):
    return Or(IntegerQ(r), And(PositiveQ(e1), PositiveQ(e2)))

cons651 = CustomConstraint(cons_f651)

def cons_f652(e1, x):
    return FreeQ(e1, x)

cons652 = CustomConstraint(cons_f652)

def cons_f653(f1, x):
    return FreeQ(f1, x)

cons653 = CustomConstraint(cons_f653)

def cons_f654(e2, x):
    return FreeQ(e2, x)

cons654 = CustomConstraint(cons_f654)

def cons_f655(f2, x):
    return FreeQ(f2, x)

cons655 = CustomConstraint(cons_f655)

def cons_f656(m, g):
    return Or(IntegerQ(m), PositiveQ(g))

cons656 = CustomConstraint(cons_f656)

def cons_f657(p, r, q):
    return PositiveIntegerQ(p + S(2), q, r)

cons657 = CustomConstraint(cons_f657)

def cons_f658(p, r, q):
    return IntegersQ(p, q, r)

cons658 = CustomConstraint(cons_f658)

def cons_f659(q, e, b, d, f, c, a):
    return Not(And(Equal(q, S(1)), SimplerQ(-a*d + b*c, -a*f + b*e)))

cons659 = CustomConstraint(cons_f659)

def cons_f660(n, q, d, x, e, f, c):
    return Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n)))

cons660 = CustomConstraint(cons_f660)

def cons_f661(r):
    return PositiveIntegerQ(r)

cons661 = CustomConstraint(cons_f661)

def cons_f662(m, n, p, q, g, e, x, b, d, f, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q), x)

cons662 = CustomConstraint(cons_f662)

def cons_f663(m, n, p, q, g, e, x, b, d, f, r, c, a):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q, r), x)

cons663 = CustomConstraint(cons_f663)
