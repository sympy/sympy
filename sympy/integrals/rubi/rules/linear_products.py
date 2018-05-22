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
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

def linear_products(rubi):
    from sympy.integrals.rubi.constraints import cons66, cons21, cons67, cons2, cons3, cons68, cons69, cons70, cons7, cons27, cons71, cons72, cons4, cons73, cons74, cons75, cons43, cons76, cons77, cons78, cons79, cons80, cons81, cons82, cons62, cons83, cons84, cons85, cons86, cons87, cons88, cons89, cons90, cons91, cons92, cons23, cons93, cons94, cons95, cons96, cons97, cons98, cons99, cons100, cons101, cons102, cons103, cons104, cons105, cons106, cons31, cons107, cons108, cons109, cons110, cons111, cons112, cons113, cons114, cons18, cons115, cons116, cons117, cons118, cons119, cons120, cons121, cons122, cons123, cons124, cons17, cons48, cons125, cons5, cons126, cons127, cons128, cons129, cons130, cons131, cons132, cons133, cons134, cons135, cons54, cons136, cons13, cons137, cons138, cons12, cons139, cons140, cons141, cons142, cons143, cons144, cons38, cons145, cons146, cons147, cons148, cons149, cons150, cons151, cons152, cons153, cons154, cons155, cons156, cons157, cons158, cons159, cons160, cons161, cons162, cons163, cons164, cons165, cons166, cons167, cons168, cons169, cons170, cons171, cons172, cons173, cons174, cons175, cons176, cons177, cons178, cons179, cons180, cons181, cons182, cons183, cons184, cons185, cons186, cons187, cons188, cons189, cons190, cons191, cons192, cons193, cons194, cons195, cons196, cons197, cons198, cons199, cons200, cons201, cons202, cons203, cons204, cons205, cons206, cons207, cons208, cons209, cons210, cons211, cons212, cons213, cons214, cons215, cons216, cons217, cons218, cons219, cons50, cons220, cons221, cons222, cons223, cons224, cons52

    pattern37 = Pattern(Integral(S(1)/x_, x_))
    def replacement37(x):
        rubi.append(37)
        return log(x)
    rule37 = ReplacementRule(pattern37, replacement37)
    pattern38 = Pattern(Integral(x_**WC('m', S(1)), x_), cons21, cons66)
    def replacement38(m, x):
        rubi.append(38)
        return x**(m + S(1))/(m + S(1))
    rule38 = ReplacementRule(pattern38, replacement38)
    pattern39 = Pattern(Integral(S(1)/(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons67)
    def replacement39(b, x, a):
        rubi.append(39)
        return log(RemoveContent(a + b*x, x))/b
    rule39 = ReplacementRule(pattern39, replacement39)
    pattern40 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons21, cons66)
    def replacement40(b, m, x, a):
        rubi.append(40)
        return (a + b*x)**(m + S(1))/(b*(m + S(1)))
    rule40 = ReplacementRule(pattern40, replacement40)
    pattern41 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons21, cons68, cons69)
    def replacement41(m, u, x, b, a):
        rubi.append(41)
        return Subst(Int((a + b*x)**m, x), x, u)/Coefficient(u, x, S(1))
    rule41 = ReplacementRule(pattern41, replacement41)
    pattern42 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons70)
    def replacement42(d, x, b, c, a):
        rubi.append(42)
        return Int(S(1)/(a*c + b*d*x**S(2)), x)
    rule42 = ReplacementRule(pattern42, replacement42)
    pattern43 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement43(d, x, b, c, a):
        rubi.append(43)
        return b*Int(S(1)/(a + b*x), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x), x)/(-a*d + b*c)
    rule43 = ReplacementRule(pattern43, replacement43)
    pattern44 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons72, cons66)
    def replacement44(n, d, m, x, b, c, a):
        rubi.append(44)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c))
    rule44 = ReplacementRule(pattern44, replacement44)
    pattern45 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons7, cons27, cons70, cons73)
    def replacement45(d, m, x, b, c, a):
        rubi.append(45)
        return S(2)*a*c*m*Int((a + b*x)**(m + S(-1))*(c + d*x)**(m + S(-1)), x)/(S(2)*m + S(1)) + x*(a + b*x)**m*(c + d*x)**m/(S(2)*m + S(1))
    rule45 = ReplacementRule(pattern45, replacement45)
    pattern46 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(3)/2)*(c_ + x_*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons7, cons27, cons70)
    def replacement46(d, x, b, c, a):
        rubi.append(46)
        return x/(a*c*sqrt(a + b*x)*sqrt(c + d*x))
    rule46 = ReplacementRule(pattern46, replacement46)
    pattern47 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons7, cons27, cons70, cons74)
    def replacement47(d, m, x, b, c, a):
        rubi.append(47)
        return -x*(a + b*x)**(m + S(1))*(c + d*x)**(m + S(1))/(S(2)*a*c*(m + S(1))) + (S(2)*m + S(3))*Int((a + b*x)**(m + S(1))*(c + d*x)**(m + S(1)), x)/(S(2)*a*c*(m + S(1)))
    rule47 = ReplacementRule(pattern47, replacement47)
    pattern48 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons70, cons75)
    def replacement48(d, m, x, b, c, a):
        rubi.append(48)
        return Int((a*c + b*d*x**S(2))**m, x)
    rule48 = ReplacementRule(pattern48, replacement48)
    pattern49 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons70, cons43, cons76)
    def replacement49(d, x, b, c, a):
        rubi.append(49)
        return acosh(b*x/a)/b
    rule49 = ReplacementRule(pattern49, replacement49)
    pattern50 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons70)
    def replacement50(d, x, b, c, a):
        rubi.append(50)
        return S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x))
    rule50 = ReplacementRule(pattern50, replacement50)
    pattern51 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons7, cons27, cons21, cons70, cons77)
    def replacement51(d, m, x, b, c, a):
        rubi.append(51)
        return (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((a*c + b*d*x**S(2))**m, x)
    rule51 = ReplacementRule(pattern51, replacement51)
    pattern52 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(5)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons2, cons3, cons7, cons27, cons70, cons78)
    def replacement52(d, x, b, c, a):
        rubi.append(52)
        return (-a*d + b*c)*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(2)*b) - S(2)/(b*(a + b*x)**(S(1)/4)*(c + d*x)**(S(1)/4))
    rule52 = ReplacementRule(pattern52, replacement52)
    pattern53 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(9)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons2, cons3, cons7, cons27, cons70, cons78)
    def replacement53(d, x, b, c, a):
        rubi.append(53)
        return -d*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(5)*b) - S(4)/(S(5)*b*(a + b*x)**(S(5)/4)*(c + d*x)**(S(1)/4))
    rule53 = ReplacementRule(pattern53, replacement53)
    pattern54 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons7, cons27, cons70, cons79, cons80, cons81)
    def replacement54(n, d, m, x, b, c, a):
        rubi.append(54)
        return S(2)*c*n*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(m + n + S(1)) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1)))
    rule54 = ReplacementRule(pattern54, replacement54)
    pattern55 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons7, cons27, cons70, cons79, cons80, cons82)
    def replacement55(n, d, m, x, b, c, a):
        rubi.append(55)
        return (m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(S(2)*a*(m + S(1))) - (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/(S(2)*a*d*(m + S(1)))
    rule55 = ReplacementRule(pattern55, replacement55)
    pattern56 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons62, cons83)
    def replacement56(n, d, m, x, b, c, a):
        rubi.append(56)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)
    rule56 = ReplacementRule(pattern56, replacement56)
    pattern57 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons84, cons85, cons86)
    def replacement57(n, d, m, x, b, c, a):
        rubi.append(57)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)
    rule57 = ReplacementRule(pattern57, replacement57)
    pattern58 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons71, cons87, cons88)
    def replacement58(n, d, x, b, c, a):
        rubi.append(58)
        return (-a*d + b*c)*Int((c + d*x)**(n + S(-1))/(a + b*x), x)/b + (c + d*x)**n/(b*n)
    rule58 = ReplacementRule(pattern58, replacement58)
    pattern59 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons71, cons87, cons89)
    def replacement59(n, d, x, b, c, a):
        rubi.append(59)
        return b*Int((c + d*x)**(n + S(1))/(a + b*x), x)/(-a*d + b*c) - (c + d*x)**(n + S(1))/((n + S(1))*(-a*d + b*c))
    rule59 = ReplacementRule(pattern59, replacement59)
    def With60(d, x, b, c, a):
        q = Rt((-a*d + b*c)/b, S(3))
        rubi.append(60)
        return S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    pattern60 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons90)
    rule60 = ReplacementRule(pattern60, With60)
    def With61(d, x, b, c, a):
        q = Rt(-(-a*d + b*c)/b, S(3))
        rubi.append(61)
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    pattern61 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons91)
    rule61 = ReplacementRule(pattern61, With61)
    def With62(d, x, b, c, a):
        q = Rt((-a*d + b*c)/b, S(3))
        rubi.append(62)
        return -S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    pattern62 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons7, cons27, cons90)
    rule62 = ReplacementRule(pattern62, With62)
    def With63(d, x, b, c, a):
        q = Rt(-(-a*d + b*c)/b, S(3))
        rubi.append(63)
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    pattern63 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons7, cons27, cons91)
    rule63 = ReplacementRule(pattern63, With63)
    def With64(n, d, x, b, c, a):
        p = Denominator(n)
        rubi.append(64)
        return p*Subst(Int(x**(p*(n + S(1)) + S(-1))/(a*d - b*c + b*x**p), x), x, (c + d*x)**(S(1)/p))
    pattern64 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons71, cons87, cons92)
    rule64 = ReplacementRule(pattern64, With64)
    pattern65 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_/x_, x_), cons7, cons27, cons4, cons23)
    def replacement65(c, n, d, x):
        rubi.append(65)
        return -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), S(1) + d*x/c)/(c*(n + S(1)))
    rule65 = ReplacementRule(pattern65, replacement65)
    pattern66 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons7, cons27, cons4, cons71, cons23)
    def replacement66(n, d, x, b, c, a):
        rubi.append(66)
        return -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), TogetherSimplify(b*(c + d*x)/(-a*d + b*c)))/((n + S(1))*(-a*d + b*c))
    rule66 = ReplacementRule(pattern66, replacement66)
    pattern67 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons71, cons93, cons94, cons88, cons95, cons96, cons97)
    def replacement67(n, d, m, x, b, c, a):
        rubi.append(67)
        return -d*n*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1)), x)/(b*(m + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + S(1)))
    rule67 = ReplacementRule(pattern67, replacement67)
    pattern68 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons71, cons93, cons94, cons98, cons97)
    def replacement68(n, d, m, x, b, c, a):
        rubi.append(68)
        return -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c))
    rule68 = ReplacementRule(pattern68, replacement68)
    pattern69 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons71, cons93, cons88, cons99, cons100, cons101, cons97)
    def replacement69(n, d, m, x, b, c, a):
        rubi.append(69)
        return n*(-a*d + b*c)*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(b*(m + n + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1)))
    rule69 = ReplacementRule(pattern69, replacement69)
    pattern70 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons102, cons103)
    def replacement70(d, x, b, c, a):
        rubi.append(70)
        return Int(S(1)/sqrt(a*c - b**S(2)*x**S(2) - b*x*(a - c)), x)
    rule70 = ReplacementRule(pattern70, replacement70)
    pattern71 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons104, cons105)
    def replacement71(d, x, b, c, a):
        rubi.append(71)
        return S(2)*Subst(Int(S(1)/sqrt(-a*d + b*c + d*x**S(2)), x), x, sqrt(a + b*x))/sqrt(b)
    rule71 = ReplacementRule(pattern71, replacement71)
    pattern72 = Pattern(Integral(S(1)/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons71, cons106)
    def replacement72(d, x, b, c, a):
        rubi.append(72)
        return S(2)*Subst(Int(S(1)/sqrt(-a + c + x**S(2)), x), x, sqrt(a + b*x))/b
    rule72 = ReplacementRule(pattern72, replacement72)
    pattern73 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons71)
    def replacement73(d, x, b, c, a):
        rubi.append(73)
        return S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x))
    rule73 = ReplacementRule(pattern73, replacement73)
    pattern74 = Pattern(Integral((c_ + x_*WC('d', S(1)))**m_*(x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons7, cons27, cons71, cons31, cons107, cons108)
    def replacement74(d, m, x, b, c, a):
        rubi.append(74)
        return (a + b*x)**m*(c + d*x)**m*(a*c + b*d*x**S(2) + x*(a*d + b*c))**(-m)*Int((a*c + b*d*x**S(2) + x*(a*d + b*c))**m, x)
    rule74 = ReplacementRule(pattern74, replacement74)
    def With75(d, x, b, c, a):
        q = Rt(d/b, S(3))
        rubi.append(75)
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d - q*log(c + d*x)/(S(2)*d) - S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(-1))/(S(2)*d)
    pattern75 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons7, cons27, cons71, cons109)
    rule75 = ReplacementRule(pattern75, With75)
    def With76(d, x, b, c, a):
        q = Rt(-d/b, S(3))
        rubi.append(76)
        return sqrt(S(3))*q*ArcTan(-S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d + q*log(c + d*x)/(S(2)*d) + S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(1))/(S(2)*d)
    pattern76 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons7, cons27, cons71, cons110)
    rule76 = ReplacementRule(pattern76, With76)
    def With77(n, d, m, x, b, c, a):
        p = Denominator(m)
        rubi.append(77)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))/(b - d*x**p), x), x, (a + b*x)**(S(1)/p)*(c + d*x)**(-S(1)/p))
    pattern77 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons71, cons93, cons107, cons111)
    rule77 = ReplacementRule(pattern77, With77)
    def With78(n, d, m, x, b, c, a):
        p = Denominator(m)
        rubi.append(78)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))*(-a*d/b + c + d*x**p/b)**n, x), x, (a + b*x)**(S(1)/p))/b
    pattern78 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons71, cons93, cons107, cons92, cons112, cons97)
    rule78 = ReplacementRule(pattern78, With78)
    pattern79 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons113, cons66, cons114)
    def replacement79(n, d, m, x, b, c, a):
        rubi.append(79)
        return -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c))
    rule79 = ReplacementRule(pattern79, replacement79)
    pattern80 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons7, cons27, cons21, cons4, cons18, cons115)
    def replacement80(n, d, m, x, b, c):
        rubi.append(80)
        return c**n*(b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*x/c)/(b*(m + S(1)))
    rule80 = ReplacementRule(pattern80, replacement80)
    pattern81 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons7, cons27, cons21, cons4, cons23, cons116)
    def replacement81(n, d, m, x, b, c):
        rubi.append(81)
        return (-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*Hypergeometric2F1(-m, n + S(1), n + S(2), S(1) + d*x/c)/(d*(n + S(1)))
    rule81 = ReplacementRule(pattern81, replacement81)
    pattern82 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons7, cons27, cons21, cons4, cons18, cons23, cons117, cons118, cons119)
    def replacement82(n, d, m, x, b, c):
        rubi.append(82)
        return c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n, x)
    rule82 = ReplacementRule(pattern82, replacement82)
    pattern83 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons7, cons27, cons21, cons4, cons18, cons23, cons117, cons118)
    def replacement83(n, d, m, x, b, c):
        rubi.append(83)
        return (b*x)**FracPart(m)*(-b*c/d)**IntPart(m)*(-d*x/c)**(-FracPart(m))*Int((-d*x/c)**m*(c + d*x)**n, x)
    rule83 = ReplacementRule(pattern83, replacement83)
    pattern84 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons7, cons27, cons21, cons71, cons18, cons85)
    def replacement84(n, d, m, x, b, c, a):
        rubi.append(84)
        return b**(-n + S(-1))*(a + b*x)**(m + S(1))*(-a*d + b*c)**n*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(m + S(1))
    rule84 = ReplacementRule(pattern84, replacement84)
    pattern85 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons18, cons23, cons120, cons121)
    def replacement85(n, d, m, x, b, c, a):
        print(rubi)
        rubi.append(85)
        return (b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(b*(m + S(1)))
    rule85 = ReplacementRule(pattern85, replacement85)
    pattern86 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons7, cons27, cons21, cons4, cons71, cons18, cons23, cons122)
    def replacement86(n, d, m, x, b, c, a):
        rubi.append(86)
        return (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x)
    rule86 = ReplacementRule(pattern86, replacement86)
    pattern87 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons21, cons4, cons68, cons123)
    def replacement87(n, d, m, u, x, b, c, a):
        rubi.append(87)
        return Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, u)/Coefficient(u, x, S(1))
    rule87 = ReplacementRule(pattern87, replacement87)
    pattern88 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons70, cons124, cons17)
    def replacement88(n, d, e, p, m, x, f, b, c, a):
        rubi.append(88)
        return Int((e + f*x)**p*(a*c + b*d*x**S(2))**m, x)
    rule88 = ReplacementRule(pattern88, replacement88)
    pattern89 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons126, cons127)
    def replacement89(n, d, e, p, x, f, b, c, a):
        rubi.append(89)
        return b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2)))
    rule89 = ReplacementRule(pattern89, replacement89)
    pattern90 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons27, cons48, cons125, cons4, cons128, cons129, cons130)
    def replacement90(n, d, p, e, x, f, b, a):
        rubi.append(90)
        return Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x)
    rule90 = ReplacementRule(pattern90, replacement90)
    pattern91 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons27, cons48, cons125, cons4, cons128, cons131, cons132, cons133)
    def replacement91(n, d, p, e, x, f, b, a):
        rubi.append(91)
        return Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x)
    rule91 = ReplacementRule(pattern91, replacement91)
    pattern92 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons71, cons134)
    def replacement92(n, d, e, p, x, f, b, c, a):
        rubi.append(92)
        return Int(ExpandIntegrand((a + b*x)*(c + d*x)**n*(e + f*x)**p, x), x)
    rule92 = ReplacementRule(pattern92, replacement92)
    pattern93 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons135, cons54, cons136)
    def replacement93(n, d, e, p, x, f, b, c, a):
        rubi.append(93)
        return b*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/f - (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e))
    rule93 = ReplacementRule(pattern93, replacement93)
    pattern94 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons126, cons13, cons137, cons138)
    def replacement94(n, d, e, p, x, f, b, c, a):
        rubi.append(94)
        return -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e))
    rule94 = ReplacementRule(pattern94, replacement94)
    pattern95 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons126, cons12, cons139)
    def replacement95(n, d, e, p, x, f, b, c, a):
        rubi.append(95)
        return -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e))
    rule95 = ReplacementRule(pattern95, replacement95)
    pattern96 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons126)
    def replacement96(n, d, e, p, x, f, b, c, a):
        rubi.append(96)
        return b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))) + (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**p, x)/(d*f*(n + p + S(2)))
    rule96 = ReplacementRule(pattern96, replacement96)
    pattern97 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons126, cons140, cons141)
    def replacement97(n, d, e, p, x, f, b, c, a):
        rubi.append(97)
        return b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(S(2)*a*d*f*(n + p + S(3)) + b*d*f*x*(n + p + S(2)) - b*(c*f*(p + S(2)) + d*e*(n + S(2))))/(d**S(2)*f**S(2)*(n + p + S(2))*(n + p + S(3)))
    rule97 = ReplacementRule(pattern97, replacement97)
    pattern98 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons21, cons4, cons5, cons70, cons142, cons12, cons143, cons144)
    def replacement98(n, d, p, m, x, f, b, c, a):
        rubi.append(98)
        return a*Int((f*x)**p*(a + b*x)**n*(c + d*x)**n, x) + b*Int((f*x)**(p + S(1))*(a + b*x)**n*(c + d*x)**n, x)/f
    rule98 = ReplacementRule(pattern98, replacement98)
    pattern99 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons38)
    def replacement99(d, e, p, x, f, b, c, a):
        rubi.append(99)
        return Int(ExpandIntegrand((e + f*x)**p/((a + b*x)*(c + d*x)), x), x)
    rule99 = ReplacementRule(pattern99, replacement99)
    pattern100 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons13, cons145)
    def replacement100(d, e, p, x, f, b, c, a):
        rubi.append(100)
        return (-a*f + b*e)*Int((e + f*x)**(p + S(-1))/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))/(c + d*x), x)/(-a*d + b*c)
    rule100 = ReplacementRule(pattern100, replacement100)
    pattern101 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons13, cons146)
    def replacement101(d, e, p, x, f, b, c, a):
        rubi.append(101)
        return f*(e + f*x)**(p + S(-1))/(b*d*(p + S(-1))) + Int((e + f*x)**(p + S(-2))*(-a*c*f**S(2) + b*d*e**S(2) + f*x*(-a*d*f - b*c*f + S(2)*b*d*e))/((a + b*x)*(c + d*x)), x)/(b*d)
    rule101 = ReplacementRule(pattern101, replacement101)
    pattern102 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons13, cons137)
    def replacement102(d, e, p, x, f, b, c, a):
        rubi.append(102)
        return f*(e + f*x)**(p + S(1))/((p + S(1))*(-a*f + b*e)*(-c*f + d*e)) + Int((e + f*x)**(p + S(1))*(-a*d*f - b*c*f + b*d*e - b*d*f*x)/((a + b*x)*(c + d*x)), x)/((-a*f + b*e)*(-c*f + d*e))
    rule102 = ReplacementRule(pattern102, replacement102)
    pattern103 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons147)
    def replacement103(d, e, p, x, f, b, c, a):
        rubi.append(103)
        return b*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - d*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c)
    rule103 = ReplacementRule(pattern103, replacement103)
    pattern104 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons148, cons149, cons137)
    def replacement104(n, d, e, p, x, f, b, c, a):
        rubi.append(104)
        return Int(ExpandIntegrand((e + f*x)**FractionalPart(p), (c + d*x)**n*(e + f*x)**IntegerPart(p)/(a + b*x), x), x)
    rule104 = ReplacementRule(pattern104, replacement104)
    pattern105 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons150, cons151)
    def replacement105(n, d, e, p, m, x, f, b, c, a):
        rubi.append(105)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x)
    rule105 = ReplacementRule(pattern105, replacement105)
    pattern106 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons152)
    def replacement106(n, d, e, p, x, f, b, c, a):
        rubi.append(106)
        return (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*d + b*c)**S(2)/(d**S(2)*(n + S(1))*(-c*f + d*e)) - Int((c + d*x)**(n + S(1))*(e + f*x)**p*Simp(a**S(2)*d**S(2)*f*(n + p + S(2)) - S(2)*a*b*d*(c*f*(p + S(1)) + d*e*(n + S(1))) + b**S(2)*c*(c*f*(p + S(1)) + d*e*(n + S(1))) - b**S(2)*d*x*(n + S(1))*(-c*f + d*e), x), x)/(d**S(2)*(n + S(1))*(-c*f + d*e))
    rule106 = ReplacementRule(pattern106, replacement106)
    pattern107 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons140)
    def replacement107(n, d, e, p, x, f, b, c, a):
        rubi.append(107)
        return b*(a + b*x)*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(3))) + Int((c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(n + p + S(3)) + b*x*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e), x), x)/(d*f*(n + p + S(3)))
    rule107 = ReplacementRule(pattern107, replacement107)
    def With108(d, e, x, f, b, c, a):
        q = Rt((-c*f + d*e)/(-a*f + b*e), S(3))
        rubi.append(108)
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(-c*f + d*e) + q*log(e + f*x)/(-S(2)*c*f + S(2)*d*e) - S(3)*q*log(q*(a + b*x)**(S(1)/3) - (c + d*x)**(S(1)/3))/(-S(2)*c*f + S(2)*d*e)
    pattern108 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons153)
    rule108 = ReplacementRule(pattern108, With108)
    pattern109 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons154)
    def replacement109(d, e, x, f, b, c, a):
        rubi.append(109)
        return b*f*Subst(Int(S(1)/(b*f**S(2)*x**S(2) + d*(-a*f + b*e)**S(2)), x), x, sqrt(a + b*x)*sqrt(c + d*x))
    rule109 = ReplacementRule(pattern109, replacement109)
    def With110(n, d, e, m, x, f, b, c, a):
        q = Denominator(m)
        rubi.append(110)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))/(-a*f + b*e - x**q*(-c*f + d*e)), x), x, (a + b*x)**(S(1)/q)*(c + d*x)**(-S(1)/q))
    pattern110 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons155, cons93, cons107, cons156)
    rule110 = ReplacementRule(pattern110, With110)
    pattern111 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons157, cons87, cons88, cons158)
    def replacement111(n, d, e, p, m, x, f, b, c, a):
        rubi.append(111)
        return -n*(-c*f + d*e)*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p, x)/((m + S(1))*(-a*f + b*e)) + (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e))
    rule111 = ReplacementRule(pattern111, replacement111)
    pattern112 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons159, cons160, cons66)
    def replacement112(n, d, e, p, m, x, f, b, c, a):
        rubi.append(112)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule112 = ReplacementRule(pattern112, replacement112)
    pattern113 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons159, cons161)
    def replacement113(n, d, e, p, m, x, f, b, c, a):
        rubi.append(113)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + (a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule113 = ReplacementRule(pattern113, replacement113)
    pattern114 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons162, cons94, cons88, cons163, cons164)
    def replacement114(n, d, e, p, m, x, f, b, c, a):
        rubi.append(114)
        return (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p/(b*(m + S(1))) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(-1))*Simp(c*f*p + d*e*n + d*f*x*(n + p), x), x)/(b*(m + S(1)))
    rule114 = ReplacementRule(pattern114, replacement114)
    pattern115 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons162, cons94, cons165, cons164)
    def replacement115(n, d, e, p, m, x, f, b, c, a):
        rubi.append(115)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(1))*(-a*d + b*c)/(b*(m + S(1))*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-2))*(e + f*x)**p*Simp(a*d*(c*f*(p + S(1)) + d*e*(n + S(-1))) + b*c*(-c*f*(m + p + S(2)) + d*e*(m - n + S(2))) + d*x*(a*d*f*(n + p) + b*(-c*f*(m + n + p + S(1)) + d*e*(m + S(1)))), x), x)/(b*(m + S(1))*(-a*f + b*e))
    rule115 = ReplacementRule(pattern115, replacement115)
    pattern116 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons5, cons162, cons94, cons88, cons164)
    def replacement116(n, d, e, p, m, x, f, b, c, a):
        rubi.append(116)
        return (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(c*f*(m + p + S(2)) + d*e*n + d*f*x*(m + n + p + S(2)), x), x)/((m + S(1))*(-a*f + b*e))
    rule116 = ReplacementRule(pattern116, replacement116)
    pattern117 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons166, cons167, cons17)
    def replacement117(n, d, e, p, m, x, f, b, c, a):
        rubi.append(117)
        return b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1)))
    rule117 = ReplacementRule(pattern117, replacement117)
    pattern118 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons162, cons168, cons88, cons167, cons169)
    def replacement118(n, d, e, p, m, x, f, b, c, a):
        rubi.append(118)
        return (a + b*x)**m*(c + d*x)**n*(e + f*x)**(p + S(1))/(f*(m + n + p + S(1))) - Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(a*n*(-c*f + d*e) + c*m*(-a*f + b*e) + x*(b*n*(-c*f + d*e) + d*m*(-a*f + b*e)), x), x)/(f*(m + n + p + S(1)))
    rule118 = ReplacementRule(pattern118, replacement118)
    pattern119 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons166, cons167, cons170)
    def replacement119(n, d, e, p, m, x, f, b, c, a):
        rubi.append(119)
        return b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1)))
    rule119 = ReplacementRule(pattern119, replacement119)
    pattern120 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons94, cons17, cons171)
    def replacement120(n, d, e, p, m, x, f, b, c, a):
        rubi.append(120)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule120 = ReplacementRule(pattern120, replacement120)
    pattern121 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons31, cons94, cons170)
    def replacement121(n, d, e, p, m, x, f, b, c, a):
        rubi.append(121)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule121 = ReplacementRule(pattern121, replacement121)
    pattern122 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons172, cons173)
    def replacement122(n, d, e, m, x, f, b, c, a):
        rubi.append(122)
        return b*Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x)/f - (-a*f + b*e)*Int((a + b*x)**(m + S(-1))*(c + d*x)**n/(e + f*x), x)/f
    rule122 = ReplacementRule(pattern122, replacement122)
    pattern123 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons174)
    def replacement123(d, e, x, f, b, c, a):
        rubi.append(123)
        return -S(4)*Subst(Int(x**S(2)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4))
    rule123 = ReplacementRule(pattern123, replacement123)
    pattern124 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons175)
    def replacement124(d, e, x, f, b, c, a):
        rubi.append(124)
        return sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(1)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x)
    rule124 = ReplacementRule(pattern124, replacement124)
    pattern125 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons174)
    def replacement125(d, e, x, f, b, c, a):
        rubi.append(125)
        return -S(4)*Subst(Int(S(1)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4))
    rule125 = ReplacementRule(pattern125, replacement125)
    pattern126 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons175)
    def replacement126(d, e, x, f, b, c, a):
        rubi.append(126)
        return sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(3)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x)
    rule126 = ReplacementRule(pattern126, replacement126)
    pattern127 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons176, cons177, cons178, cons179)
    def replacement127(d, e, x, f, b, c):
        rubi.append(127)
        return S(2)*sqrt(e)*EllipticE(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/b
    rule127 = ReplacementRule(pattern127, replacement127)
    pattern128 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons176, cons177, cons178, cons180)
    def replacement128(d, e, x, f, b, c):
        rubi.append(128)
        return sqrt(-b*x)*Int(sqrt(e + f*x)/(sqrt(-b*x)*sqrt(c + d*x)), x)/sqrt(b*x)
    rule128 = ReplacementRule(pattern128, replacement128)
    pattern129 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons176, cons181)
    def replacement129(d, e, x, f, b, c):
        rubi.append(129)
        return sqrt(S(1) + d*x/c)*sqrt(e + f*x)*Int(sqrt(S(1) + f*x/e)/(sqrt(b*x)*sqrt(S(1) + d*x/c)), x)/(sqrt(S(1) + f*x/e)*sqrt(c + d*x))
    rule129 = ReplacementRule(pattern129, replacement129)
    pattern130 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons120, cons182, cons183, cons184)
    def replacement130(d, e, x, f, b, c, a):
        rubi.append(130)
        return S(2)*EllipticE(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*f + b*e)/d, S(2))/b
    rule130 = ReplacementRule(pattern130, replacement130)
    pattern131 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons185, cons183)
    def replacement131(d, e, x, f, b, c, a):
        rubi.append(131)
        return sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*Int(sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))), x)/(sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x))
    rule131 = ReplacementRule(pattern131, replacement131)
    pattern132 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons177, cons178, cons186)
    def replacement132(d, e, x, f, b, c):
        rubi.append(132)
        return S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e))
    rule132 = ReplacementRule(pattern132, replacement132)
    pattern133 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons177, cons178, cons187)
    def replacement133(d, e, x, f, b, c):
        rubi.append(133)
        return S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e))
    rule133 = ReplacementRule(pattern133, replacement133)
    pattern134 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons7, cons27, cons48, cons125, cons181)
    def replacement134(d, e, x, f, b, c):
        rubi.append(134)
        return sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)*Int(S(1)/(sqrt(b*x)*sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)), x)/(sqrt(c + d*x)*sqrt(e + f*x))
    rule134 = ReplacementRule(pattern134, replacement134)
    pattern135 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons120, cons182, cons156, cons188, cons189)
    def replacement135(d, e, x, f, b, c, a):
        rubi.append(135)
        return S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b
    rule135 = ReplacementRule(pattern135, replacement135)
    pattern136 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons120, cons182, cons156, cons188, cons190)
    def replacement136(d, e, x, f, b, c, a):
        rubi.append(136)
        return S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b
    rule136 = ReplacementRule(pattern136, replacement136)
    pattern137 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons185, cons156, cons188)
    def replacement137(d, e, x, f, b, c, a):
        rubi.append(137)
        return sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*Int(S(1)/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))*sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))), x)/(sqrt(c + d*x)*sqrt(e + f*x))
    rule137 = ReplacementRule(pattern137, replacement137)
    def With138(d, e, x, f, b, c, a):
        q = Rt(b*(-a*f + b*e)/(-a*d + b*c)**S(2), S(3))
        rubi.append(138)
        return -sqrt(S(3))*ArcTan(S(2)*sqrt(S(3))*q*(c + d*x)**(S(2)/3)/(S(3)*(e + f*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*q*(-a*d + b*c)) - log(a + b*x)/(S(2)*q*(-a*d + b*c)) + S(3)*log(q*(c + d*x)**(S(2)/3) - (e + f*x)**(S(1)/3))/(S(4)*q*(-a*d + b*c))
    pattern138 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons191)
    rule138 = ReplacementRule(pattern138, With138)
    pattern139 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_/((x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons191, cons17, cons94)
    def replacement139(d, e, m, x, f, b, c, a):
        rubi.append(139)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(S(2)/3)*(e + f*x)**(S(2)/3)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + f*Int((a + b*x)**(m + S(1))*(a*d*(S(3)*m + S(1)) - S(3)*b*c*(S(3)*m + S(5)) - S(2)*b*d*x*(S(3)*m + S(7)))/((c + d*x)**(S(1)/3)*(e + f*x)**(S(1)/3)), x)/(S(6)*(m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule139 = ReplacementRule(pattern139, replacement139)
    pattern140 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons21, cons4, cons5, cons70, cons124, cons43, cons177)
    def replacement140(n, d, p, m, x, f, b, c, a):
        rubi.append(140)
        return Int((f*x)**p*(a*c + b*d*x**S(2))**m, x)
    rule140 = ReplacementRule(pattern140, replacement140)
    pattern141 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons21, cons4, cons5, cons70, cons124)
    def replacement141(n, d, p, m, x, f, b, c, a):
        rubi.append(141)
        return (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((f*x)**p*(a*c + b*d*x**S(2))**m, x)
    rule141 = ReplacementRule(pattern141, replacement141)
    pattern142 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons125, cons21, cons4, cons5, cons70, cons192, cons144)
    def replacement142(n, d, p, m, x, f, b, c, a):
        rubi.append(142)
        return Int(ExpandIntegrand((f*x)**p*(a + b*x)**n*(c + d*x)**n, (a + b*x)**(m - n), x), x)
    rule142 = ReplacementRule(pattern142, replacement142)
    pattern143 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons4, cons5, cons193)
    def replacement143(n, d, e, p, m, x, f, b, c, a):
        rubi.append(143)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x)
    rule143 = ReplacementRule(pattern143, replacement143)
    pattern144 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons194, cons66, cons195)
    def replacement144(n, d, e, p, m, x, f, b, c, a):
        rubi.append(144)
        return b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule144 = ReplacementRule(pattern144, replacement144)
    pattern145 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons5, cons157, cons196)
    def replacement145(n, d, e, p, m, x, f, b, c, a):
        rubi.append(145)
        return (a + b*x)**(m + S(1))*(e + f*x)**(-m + S(-1))*(-a*d + b*c)**n*(-a*f + b*e)**(-n + S(-1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/(m + S(1))
    rule145 = ReplacementRule(pattern145, replacement145)
    pattern146 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons157, cons23)
    def replacement146(n, d, e, p, m, x, f, b, c, a):
        rubi.append(146)
        return ((c + d*x)*(-a*f + b*e)/((e + f*x)*(-a*d + b*c)))**(-n)*(a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/((m + S(1))*(-a*f + b*e))
    rule146 = ReplacementRule(pattern146, replacement146)
    pattern147 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons177, cons197)
    def replacement147(n, d, e, p, m, x, f, b, c):
        rubi.append(147)
        return c**n*e**p*(b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*x/c, -f*x/e)/(b*(m + S(1)))
    rule147 = ReplacementRule(pattern147, replacement147)
    pattern148 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons198, cons199)
    def replacement148(n, d, e, p, m, x, f, b, c):
        rubi.append(148)
        return (d/(-c*f + d*e))**(-p)*(-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*AppellF1(n + S(1), -m, -p, n + S(2), S(1) + d*x/c, -f*(c + d*x)/(-c*f + d*e))/(d*(n + S(1)))
    rule148 = ReplacementRule(pattern148, replacement148)
    pattern149 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons117)
    def replacement149(n, d, e, p, m, x, f, b, c):
        rubi.append(149)
        return c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n*(e + f*x)**p, x)
    rule149 = ReplacementRule(pattern149, replacement149)
    pattern150 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons18, cons23, cons38, cons120, cons200)
    def replacement150(n, d, e, p, m, x, f, b, c, a):
        rubi.append(150)
        return b**(-p + S(-1))*(b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*(-a*f + b*e)**p*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(m + S(1))
    rule150 = ReplacementRule(pattern150, replacement150)
    pattern151 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons18, cons23, cons38, cons201, cons202)
    def replacement151(n, d, e, p, m, x, f, b, c, a):
        rubi.append(151)
        return (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x)
    rule151 = ReplacementRule(pattern151, replacement151)
    pattern152 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons147, cons120, cons182, cons203, cons204)
    def replacement152(n, d, e, p, m, x, f, b, c, a):
        rubi.append(152)
        return (b/(-a*d + b*c))**(-n)*(b/(-a*f + b*e))**(-p)*(a + b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(b*(m + S(1)))
    rule152 = ReplacementRule(pattern152, replacement152)
    pattern153 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons147, cons120, cons205)
    def replacement153(n, d, e, p, m, x, f, b, c, a):
        rubi.append(153)
        return (b/(-a*f + b*e))**(-IntPart(p))*(b*(e + f*x)/(-a*f + b*e))**(-FracPart(p))*(e + f*x)**FracPart(p)*Int((a + b*x)**m*(c + d*x)**n*(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))**p, x)
    rule153 = ReplacementRule(pattern153, replacement153)
    pattern154 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons18, cons23, cons147, cons201, cons202, cons206)
    def replacement154(n, d, e, p, m, x, f, b, c, a):
        rubi.append(154)
        return (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x)
    rule154 = ReplacementRule(pattern154, replacement154)
    pattern155 = Pattern(Integral((e_ + u_*WC('f', S(1)))**WC('p', S(1))*(u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons21, cons4, cons5, cons68, cons69)
    def replacement155(n, d, p, e, m, u, x, f, b, c, a):
        rubi.append(155)
        return Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, u)/Coefficient(u, x, S(1))
    rule155 = ReplacementRule(pattern155, replacement155)
    pattern156 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons207)
    def replacement156(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(156)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)*(g + h*x), x), x)
    rule156 = ReplacementRule(pattern156, replacement156)
    pattern157 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons72, cons66, cons210)
    def replacement157(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(157)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(2)*d*f*h*m - a*b*(-c*f*h*(m + S(1)) + d*(e*h + f*g)) + b**S(2)*d*e*g + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)) + (a*d*f*h*m + b*(-c*f*h*(m + S(2)) + d*(e*h + f*g)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d)
    rule157 = ReplacementRule(pattern157, replacement157)
    pattern158 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons93, cons94, cons89)
    def replacement158(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(158)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*c*d*f*h*(n + S(1)) + a*b*(c**S(2)*f*h*(m + S(1)) - c*d*(e*h + f*g)*(m + n + S(2)) + d**S(2)*e*g*(m + S(1))) + b**S(2)*c*d*e*g*(n + S(1)) + x*(a**S(2)*d**S(2)*f*h*(n + S(1)) - a*b*d**S(2)*(n + S(1))*(e*h + f*g) + b**S(2)*(c**S(2)*f*h*(m + S(1)) - c*d*(m + S(1))*(e*h + f*g) + d**S(2)*e*g*(m + n + S(2)))))/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)) - (a**S(2)*d**S(2)*f*h*(n**S(2) + S(3)*n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m**S(2) + S(3)*m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m**S(2) + m*(S(2)*n + S(5)) + n**S(2) + S(5)*n + S(6))))*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1)), x)/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2))
    rule158 = ReplacementRule(pattern158, replacement158)
    pattern159 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons211)
    def replacement159(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(159)
        return (-d*(m + n + S(3))*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)) + f*h/b**S(2))*Int((a + b*x)**(m + S(2))*(c + d*x)**n, x) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(3)*d*f*h*(n + S(2)) - a**S(2)*b*(c*f*h*m - d*(e*h + f*g)*(m + n + S(3))) - a*b**S(2)*(c*(e*h + f*g) + d*e*g*(S(2)*m + n + S(4))) + b**S(3)*c*e*g*(m + S(2)) + b*x*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2)))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2))
    rule159 = ReplacementRule(pattern159, replacement159)
    pattern160 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons212, cons66, cons213)
    def replacement160(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(160)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*d*f*h*(n + S(2)) + a*b*(c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*d*e*g*(m + n + S(3)) + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))) - (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3)))
    rule160 = ReplacementRule(pattern160, replacement160)
    pattern161 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons214, cons213)
    def replacement161(h, n, d, e, m, g, x, f, b, c, a):
        rubi.append(161)
        return -(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a*d*f*h*(n + S(2)) + b*c*f*h*(m + S(2)) - b*d*f*h*x*(m + n + S(2)) - b*d*(e*h + f*g)*(m + n + S(3)))/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))) + (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**m*(c + d*x)**n, x)/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3)))
    rule161 = ReplacementRule(pattern161, replacement161)
    pattern162 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons215)
    def replacement162(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(162)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x), x), x)
    rule162 = ReplacementRule(pattern162, replacement162)
    pattern163 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons93, cons94, cons88, cons17)
    def replacement163(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(163)
        return (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e))
    rule163 = ReplacementRule(pattern163, replacement163)
    pattern164 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons5, cons93, cons94, cons88, cons170)
    def replacement164(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(164)
        return (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e))
    rule164 = ReplacementRule(pattern164, replacement164)
    pattern165 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons31, cons94, cons17)
    def replacement165(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(165)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule165 = ReplacementRule(pattern165, replacement165)
    pattern166 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons31, cons94, cons170)
    def replacement166(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(166)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule166 = ReplacementRule(pattern166, replacement166)
    pattern167 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons31, cons168, cons144, cons17)
    def replacement167(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(167)
        return h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2)))
    rule167 = ReplacementRule(pattern167, replacement167)
    pattern168 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons31, cons168, cons144, cons170)
    def replacement168(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(168)
        return h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2)))
    rule168 = ReplacementRule(pattern168, replacement168)
    pattern169 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons194, cons66, cons195)
    def replacement169(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(169)
        return (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e))
    rule169 = ReplacementRule(pattern169, replacement169)
    pattern170 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement170(h, d, e, p, g, x, f, b, c, a):
        rubi.append(170)
        return (-a*h + b*g)*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - (-c*h + d*g)*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c)
    rule170 = ReplacementRule(pattern170, replacement170)
    pattern171 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons4, cons5, cons217)
    def replacement171(h, n, d, e, p, g, x, f, b, c, a):
        rubi.append(171)
        return h*Int((c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((c + d*x)**n*(e + f*x)**p/(a + b*x), x)/b
    rule171 = ReplacementRule(pattern171, replacement171)
    pattern172 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons188, cons218)
    def replacement172(h, d, e, a, x, f, b, c, g):
        rubi.append(172)
        return h*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)), x)/f + (-e*h + f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)), x)/f
    rule172 = ReplacementRule(pattern172, replacement172)
    pattern173 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons219)
    def replacement173(h, n, d, e, p, m, g, x, f, b, c, a):
        rubi.append(173)
        return h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x)/b
    rule173 = ReplacementRule(pattern173, replacement173)
    pattern174 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons50, cons13, cons145)
    def replacement174(h, d, q, e, p, g, x, f, b, c, a):
        rubi.append(174)
        return (-a*f + b*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(c + d*x), x)/(-a*d + b*c)
    rule174 = ReplacementRule(pattern174, replacement174)
    pattern175 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement175(h, d, e, g, x, f, b, c, a):
        rubi.append(175)
        return -S(2)*sqrt(d*(e + f*x)/(-c*f + d*e))*sqrt(d*(g + h*x)/(-c*h + d*g))*EllipticPi(-b*(-c*f + d*e)/(f*(-a*d + b*c)), asin(sqrt(-f/(-c*f + d*e))*sqrt(c + d*x)), h*(-c*f + d*e)/(f*(-c*h + d*g)))/(sqrt(-f/(-c*f + d*e))*sqrt(e + f*x)*sqrt(g + h*x)*(-a*d + b*c))
    rule175 = ReplacementRule(pattern175, replacement175)
    pattern176 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons80)
    def replacement176(h, n, d, e, g, x, f, b, c, a):
        rubi.append(176)
        return Int(ExpandIntegrand(S(1)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), (c + d*x)**(n + S(1)/2)/(a + b*x), x), x)
    rule176 = ReplacementRule(pattern176, replacement176)
    pattern177 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement177(h, d, e, g, x, f, b, c, a):
        rubi.append(177)
        return (-a*f + b*e)*(-a*h + b*g)*Int(S(1)/((a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2) + Int((-a*f*h + b*e*h + b*f*g + b*f*h*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2)
    rule177 = ReplacementRule(pattern177, replacement177)
    pattern178 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement178(h, d, e, g, x, f, b, c, a):
        rubi.append(178)
        return -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g))
    rule178 = ReplacementRule(pattern178, replacement178)
    pattern179 = Pattern(Integral(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement179(h, d, e, g, x, f, b, c, a):
        rubi.append(179)
        return -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*(-c*h + d*g)*Subst(Int(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))/sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1)), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)**S(2))
    rule179 = ReplacementRule(pattern179, replacement179)
    pattern180 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement180(h, d, e, g, x, f, b, c, a):
        rubi.append(180)
        return S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/((-b*x**S(2) + h)*sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x))
    rule180 = ReplacementRule(pattern180, replacement180)
    pattern181 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement181(h, d, e, g, x, f, b, c, a):
        rubi.append(181)
        return b*Int(sqrt(c + d*x)/((a + b*x)**(S(3)/2)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c) - d*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c)
    rule181 = ReplacementRule(pattern181, replacement181)
    pattern182 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement182(h, d, e, g, x, f, b, c, a):
        rubi.append(182)
        return sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)/(h*sqrt(e + f*x)) - (-c*f + d*e)*(-e*h + f*g)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*(e + f*x)**(S(3)/2)*sqrt(g + h*x)), x)/(S(2)*f*h) + (-c*f + d*e)*(-S(2)*a*f*h + b*e*h + b*f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h) + (a*d*f*h - b*(-c*f*h + d*e*h + d*f*g))*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h)
    rule182 = ReplacementRule(pattern182, replacement182)
    pattern183 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons216)
    def replacement183(h, d, e, g, x, f, b, c, a):
        rubi.append(183)
        return b*Int(sqrt(a + b*x)*sqrt(c + d*x)/(sqrt(e + f*x)*sqrt(g + h*x)), x)/d - (-a*d + b*c)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/d
    rule183 = ReplacementRule(pattern183, replacement183)
    pattern184 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons220)
    def replacement184(h, n, d, q, e, p, m, g, x, f, b, c, a):
        rubi.append(184)
        return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x)
    rule184 = ReplacementRule(pattern184, replacement184)
    pattern185 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons221, cons219)
    def replacement185(h, n, d, q, e, p, m, g, x, f, b, c, a):
        rubi.append(185)
        return h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b
    rule185 = ReplacementRule(pattern185, replacement185)
    pattern186 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons222)
    def replacement186(h, n, d, q, e, p, m, g, x, f, b, c, a):
        rubi.append(186)
        return Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x)
    rule186 = ReplacementRule(pattern186, replacement186)
    pattern187 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(u_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(u_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons21, cons4, cons5, cons50, cons68, cons69)
    def replacement187(h, n, d, q, e, p, m, u, g, x, f, b, c, a):
        rubi.append(187)
        return Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, u)/Coefficient(u, x, S(1))
    rule187 = ReplacementRule(pattern187, replacement187)
    pattern188 = Pattern(Integral(((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_*WC('i', S(1)))**r_, x_), cons2, cons3, cons7, cons27, cons48, cons125, cons208, cons209, cons224, cons21, cons4, cons5, cons50, cons52, cons223)
    def replacement188(h, n, d, q, e, p, m, r, g, x, f, b, c, i, a):
        rubi.append(188)
        return (i*(a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q)**r*(a + b*x)**(-m*r)*(c + d*x)**(-n*r)*(e + f*x)**(-p*r)*(g + h*x)**(-q*r)*Int((a + b*x)**(m*r)*(c + d*x)**(n*r)*(e + f*x)**(p*r)*(g + h*x)**(q*r), x)
    rule188 = ReplacementRule(pattern188, replacement188)
    return [rule37, rule38, rule39, rule40, rule41, rule42, rule43, rule44, rule45, rule46, rule47, rule48, rule49, rule50, rule51, rule52, rule53, rule54, rule55, rule56, rule57, rule58, rule59, rule60, rule61, rule62, rule63, rule64, rule65, rule66, rule67, rule68, rule69, rule70, rule71, rule72, rule73, rule74, rule75, rule76, rule77, rule78, rule79, rule80, rule81, rule82, rule83, rule84, rule85, rule86, rule87, rule88, rule89, rule90, rule91, rule92, rule93, rule94, rule95, rule96, rule97, rule98, rule99, rule100, rule101, rule102, rule103, rule104, rule105, rule106, rule107, rule108, rule109, rule110, rule111, rule112, rule113, rule114, rule115, rule116, rule117, rule118, rule119, rule120, rule121, rule122, rule123, rule124, rule125, rule126, rule127, rule128, rule129, rule130, rule131, rule132, rule133, rule134, rule135, rule136, rule137, rule138, rule139, rule140, rule141, rule142, rule143, rule144, rule145, rule146, rule147, rule148, rule149, rule150, rule151, rule152, rule153, rule154, rule155, rule156, rule157, rule158, rule159, rule160, rule161, rule162, rule163, rule164, rule165, rule166, rule167, rule168, rule169, rule170, rule171, rule172, rule173, rule174, rule175, rule176, rule177, rule178, rule179, rule180, rule181, rule182, rule183, rule184, rule185, rule186, rule187, rule188, ]
