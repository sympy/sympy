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


def linear_products():
    from sympy.integrals.rubi.constraints import cons2, cons68, cons19, cons69, cons3, cons70, cons71, cons72, cons8, cons29, cons73, cons74, cons4, cons75, cons76, cons77, cons45, cons78, cons79, cons80, cons81, cons82, cons83, cons84, cons64, cons85, cons86, cons87, cons88, cons89, cons90, cons91, cons92, cons93, cons94, cons25, cons95, cons96, cons97, cons98, cons99, cons100, cons101, cons102, cons103, cons104, cons105, cons106, cons107, cons108, cons33, cons109, cons110, cons111, cons112, cons113, cons114, cons115, cons116, cons21, cons117, cons118, cons119, cons120, cons121, cons122, cons123, cons124, cons125, cons126, cons20, cons50, cons127, cons5, cons128, cons129, cons130, cons131, cons132, cons133, cons134, cons135, cons136, cons137, cons56, cons138, cons13, cons139, cons140, cons12, cons141, cons142, cons143, cons144, cons145, cons146, cons40, cons147, cons148, cons149, cons150, cons151, cons152, cons153, cons154, cons155, cons156, cons157, cons158, cons159, cons160, cons161, cons162, cons163, cons164, cons165, cons166, cons167, cons168, cons169, cons170, cons171, cons172, cons173, cons174, cons175, cons176, cons177, cons178, cons179, cons180, cons181, cons182, cons183, cons184, cons185, cons186, cons187, cons188, cons189, cons190, cons191, cons192, cons193, cons194, cons195, cons196, cons197, cons198, cons199, cons200, cons201, cons202, cons203, cons204, cons205, cons206, cons207, cons208, cons209, cons210, cons211, cons212, cons213, cons214, cons215, cons216, cons217, cons218, cons219, cons220, cons221, cons52, cons222, cons223, cons224, cons225, cons226, cons54


    pattern39 = Pattern(Integral(a_, x_), cons2, cons2)
    rule39 = ReplacementRule(pattern39, replacement39)

    pattern40 = Pattern(Integral(S(1)/x_, x_))
    rule40 = ReplacementRule(pattern40, replacement40)

    pattern41 = Pattern(Integral(x_**WC('m', S(1)), x_), cons19, cons68)
    rule41 = ReplacementRule(pattern41, replacement41)

    pattern42 = Pattern(Integral(S(1)/(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons69)
    rule42 = ReplacementRule(pattern42, replacement42)

    pattern43 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons19, cons68)
    rule43 = ReplacementRule(pattern43, replacement43)

    pattern44 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons19, cons70, cons71)
    rule44 = ReplacementRule(pattern44, replacement44)

    pattern45 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons72)
    rule45 = ReplacementRule(pattern45, replacement45)

    pattern46 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule46 = ReplacementRule(pattern46, replacement46)

    pattern47 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons74, cons68)
    rule47 = ReplacementRule(pattern47, replacement47)

    pattern48 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons8, cons29, cons72, cons75)
    rule48 = ReplacementRule(pattern48, replacement48)

    pattern49 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(3)/2)*(c_ + x_*WC('d', S(1)))**(S(3)/2)), x_), cons2, cons3, cons8, cons29, cons72)
    rule49 = ReplacementRule(pattern49, replacement49)

    pattern50 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons8, cons29, cons72, cons76)
    rule50 = ReplacementRule(pattern50, replacement50)

    pattern51 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('m', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons72, cons77)
    rule51 = ReplacementRule(pattern51, replacement51)

    pattern52 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons72, cons45, cons78)
    rule52 = ReplacementRule(pattern52, replacement52)

    pattern53 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons72)
    rule53 = ReplacementRule(pattern53, replacement53)

    pattern54 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons2, cons3, cons8, cons29, cons19, cons72, cons79)
    rule54 = ReplacementRule(pattern54, replacement54)

    pattern55 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(5)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons2, cons3, cons8, cons29, cons72, cons80)
    rule55 = ReplacementRule(pattern55, replacement55)

    pattern56 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(9)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons2, cons3, cons8, cons29, cons72, cons80)
    rule56 = ReplacementRule(pattern56, replacement56)

    pattern57 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons8, cons29, cons72, cons81, cons82, cons83)
    rule57 = ReplacementRule(pattern57, replacement57)

    pattern58 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons8, cons29, cons72, cons81, cons82, cons84)
    rule58 = ReplacementRule(pattern58, replacement58)

    pattern59 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons64, cons85)
    rule59 = ReplacementRule(pattern59, replacement59)

    pattern60 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons86, cons87, cons88)
    rule60 = ReplacementRule(pattern60, replacement60)

    pattern61 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons73, cons89, cons90)
    rule61 = ReplacementRule(pattern61, replacement61)

    pattern62 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons73, cons89, cons91)
    rule62 = ReplacementRule(pattern62, replacement62)

    pattern63 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons92)
    rule63 = ReplacementRule(pattern63, With63)

    pattern64 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons93)
    rule64 = ReplacementRule(pattern64, With64)

    pattern65 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons8, cons29, cons92)
    rule65 = ReplacementRule(pattern65, With65)

    pattern66 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons8, cons29, cons93)
    rule66 = ReplacementRule(pattern66, With66)

    pattern67 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons73, cons89, cons94)
    rule67 = ReplacementRule(pattern67, With67)

    pattern68 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_/x_, x_), cons8, cons29, cons4, cons25)
    rule68 = ReplacementRule(pattern68, replacement68)

    pattern69 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + x_*WC('b', S(1))), x_), cons2, cons3, cons8, cons29, cons4, cons73, cons25)
    rule69 = ReplacementRule(pattern69, replacement69)

    pattern70 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons73, cons95, cons96, cons90, cons97, cons98, cons99)
    rule70 = ReplacementRule(pattern70, replacement70)

    pattern71 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons73, cons95, cons96, cons100, cons99)
    rule71 = ReplacementRule(pattern71, replacement71)

    pattern72 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons73, cons95, cons90, cons101, cons102, cons103, cons99)
    rule72 = ReplacementRule(pattern72, replacement72)

    pattern73 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons104, cons105)
    rule73 = ReplacementRule(pattern73, replacement73)

    pattern74 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons106, cons107)
    rule74 = ReplacementRule(pattern74, replacement74)

    pattern75 = Pattern(Integral(S(1)/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons73, cons108)
    rule75 = ReplacementRule(pattern75, replacement75)

    pattern76 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons73)
    rule76 = ReplacementRule(pattern76, replacement76)

    pattern77 = Pattern(Integral((c_ + x_*WC('d', S(1)))**m_*(x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons2, cons3, cons8, cons29, cons73, cons33, cons109, cons110)
    rule77 = ReplacementRule(pattern77, replacement77)

    pattern78 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons8, cons29, cons73, cons111)
    rule78 = ReplacementRule(pattern78, With78)

    pattern79 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons2, cons3, cons8, cons29, cons73, cons112)
    rule79 = ReplacementRule(pattern79, With79)

    pattern80 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons73, cons95, cons109, cons113)
    rule80 = ReplacementRule(pattern80, With80)

    pattern81 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons73, cons95, cons109, cons94, cons114, cons99)
    rule81 = ReplacementRule(pattern81, With81)

    pattern82 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons115, cons68, cons116)
    rule82 = ReplacementRule(pattern82, replacement82)

    pattern83 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons8, cons29, cons19, cons4, cons21, cons117)
    rule83 = ReplacementRule(pattern83, replacement83)

    pattern84 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons8, cons29, cons19, cons4, cons25, cons118)
    rule84 = ReplacementRule(pattern84, replacement84)

    pattern85 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons8, cons29, cons19, cons4, cons21, cons25, cons119, cons120, cons121)
    rule85 = ReplacementRule(pattern85, replacement85)

    pattern86 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons3, cons8, cons29, cons19, cons4, cons21, cons25, cons119, cons120)
    rule86 = ReplacementRule(pattern86, replacement86)

    pattern87 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons73, cons21, cons87)
    rule87 = ReplacementRule(pattern87, replacement87)

    pattern88 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons21, cons25, cons122, cons123)
    rule88 = ReplacementRule(pattern88, replacement88)

    pattern89 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons2, cons3, cons8, cons29, cons19, cons4, cons73, cons21, cons25, cons124)
    rule89 = ReplacementRule(pattern89, replacement89)

    pattern90 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons19, cons4, cons70, cons125)
    rule90 = ReplacementRule(pattern90, replacement90)

    pattern91 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons72, cons126, cons20)
    rule91 = ReplacementRule(pattern91, replacement91)

    pattern92 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons128, cons129)
    rule92 = ReplacementRule(pattern92, replacement92)

    pattern93 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons130, cons131, cons132)
    rule93 = ReplacementRule(pattern93, replacement93)

    pattern94 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons29, cons50, cons127, cons4, cons130, cons133, cons134, cons135)
    rule94 = ReplacementRule(pattern94, replacement94)

    pattern95 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons73, cons136)
    rule95 = ReplacementRule(pattern95, replacement95)

    pattern96 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons137, cons56, cons138)
    rule96 = ReplacementRule(pattern96, replacement96)

    pattern97 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons128, cons13, cons139, cons140)
    rule97 = ReplacementRule(pattern97, replacement97)

    pattern98 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons128, cons12, cons141)
    rule98 = ReplacementRule(pattern98, replacement98)

    pattern99 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons128)
    rule99 = ReplacementRule(pattern99, replacement99)

    pattern100 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons128, cons142, cons143)
    rule100 = ReplacementRule(pattern100, replacement100)

    pattern101 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons19, cons4, cons5, cons72, cons144, cons12, cons145, cons146)
    rule101 = ReplacementRule(pattern101, replacement101)

    pattern102 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons40)
    rule102 = ReplacementRule(pattern102, replacement102)

    pattern103 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons13, cons147)
    rule103 = ReplacementRule(pattern103, replacement103)

    pattern104 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons13, cons148)
    rule104 = ReplacementRule(pattern104, replacement104)

    pattern105 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons13, cons139)
    rule105 = ReplacementRule(pattern105, replacement105)

    pattern106 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons149)
    rule106 = ReplacementRule(pattern106, replacement106)

    pattern107 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons150, cons151, cons139)
    rule107 = ReplacementRule(pattern107, replacement107)

    pattern108 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons152, cons153)
    rule108 = ReplacementRule(pattern108, replacement108)

    pattern109 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons154)
    rule109 = ReplacementRule(pattern109, replacement109)

    pattern110 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons142)
    rule110 = ReplacementRule(pattern110, replacement110)

    pattern111 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons155)
    rule111 = ReplacementRule(pattern111, With111)

    pattern112 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons156)
    rule112 = ReplacementRule(pattern112, replacement112)

    pattern113 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons157, cons95, cons109, cons158)
    rule113 = ReplacementRule(pattern113, With113)

    pattern114 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons159, cons89, cons90, cons160)
    rule114 = ReplacementRule(pattern114, replacement114)

    pattern115 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons161, cons162, cons68)
    rule115 = ReplacementRule(pattern115, replacement115)

    pattern116 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons161, cons163)
    rule116 = ReplacementRule(pattern116, replacement116)

    pattern117 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons164, cons96, cons90, cons165, cons166)
    rule117 = ReplacementRule(pattern117, replacement117)

    pattern118 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons164, cons96, cons167, cons166)
    rule118 = ReplacementRule(pattern118, replacement118)

    pattern119 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons5, cons164, cons96, cons90, cons166)
    rule119 = ReplacementRule(pattern119, replacement119)

    pattern120 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons33, cons168, cons169, cons20)
    rule120 = ReplacementRule(pattern120, replacement120)

    pattern121 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons164, cons170, cons90, cons169, cons171)
    rule121 = ReplacementRule(pattern121, replacement121)

    pattern122 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons33, cons168, cons169, cons172)
    rule122 = ReplacementRule(pattern122, replacement122)

    pattern123 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons33, cons96, cons20, cons173)
    rule123 = ReplacementRule(pattern123, replacement123)

    pattern124 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons33, cons96, cons172)
    rule124 = ReplacementRule(pattern124, replacement124)

    pattern125 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons174, cons175)
    rule125 = ReplacementRule(pattern125, replacement125)

    pattern126 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons176)
    rule126 = ReplacementRule(pattern126, replacement126)

    pattern127 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons177)
    rule127 = ReplacementRule(pattern127, replacement127)

    pattern128 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons176)
    rule128 = ReplacementRule(pattern128, replacement128)

    pattern129 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons177)
    rule129 = ReplacementRule(pattern129, replacement129)

    pattern130 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons178, cons179, cons180, cons181)
    rule130 = ReplacementRule(pattern130, replacement130)

    pattern131 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons178, cons179, cons180, cons182)
    rule131 = ReplacementRule(pattern131, replacement131)

    pattern132 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons178, cons183)
    rule132 = ReplacementRule(pattern132, replacement132)

    pattern133 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons122, cons184, cons185, cons186)
    rule133 = ReplacementRule(pattern133, replacement133)

    pattern134 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons187, cons185)
    rule134 = ReplacementRule(pattern134, replacement134)

    pattern135 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons179, cons180, cons188)
    rule135 = ReplacementRule(pattern135, replacement135)

    pattern136 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons179, cons180, cons189)
    rule136 = ReplacementRule(pattern136, replacement136)

    pattern137 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons3, cons8, cons29, cons50, cons127, cons183)
    rule137 = ReplacementRule(pattern137, replacement137)

    pattern138 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons122, cons184, cons158, cons190, cons191)
    rule138 = ReplacementRule(pattern138, replacement138)

    pattern139 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons122, cons184, cons158, cons190, cons192)
    rule139 = ReplacementRule(pattern139, replacement139)

    pattern140 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons187, cons158, cons190)
    rule140 = ReplacementRule(pattern140, replacement140)

    pattern141 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons193)
    rule141 = ReplacementRule(pattern141, With141)

    pattern142 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_/((x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons193, cons20, cons96)
    rule142 = ReplacementRule(pattern142, replacement142)

    pattern143 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons19, cons4, cons5, cons72, cons126, cons45, cons179)
    rule143 = ReplacementRule(pattern143, replacement143)

    pattern144 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons19, cons4, cons5, cons72, cons126)
    rule144 = ReplacementRule(pattern144, replacement144)

    pattern145 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons127, cons19, cons4, cons5, cons72, cons194, cons146)
    rule145 = ReplacementRule(pattern145, replacement145)

    pattern146 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons4, cons5, cons195)
    rule146 = ReplacementRule(pattern146, replacement146)

    pattern147 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons196, cons68, cons197)
    rule147 = ReplacementRule(pattern147, replacement147)

    pattern148 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons5, cons159, cons198)
    rule148 = ReplacementRule(pattern148, replacement148)

    pattern149 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons159, cons25)
    rule149 = ReplacementRule(pattern149, replacement149)

    pattern150 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons179, cons199)
    rule150 = ReplacementRule(pattern150, replacement150)

    pattern151 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons200, cons201)
    rule151 = ReplacementRule(pattern151, replacement151)

    pattern152 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons119)
    rule152 = ReplacementRule(pattern152, replacement152)

    pattern153 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons21, cons25, cons40, cons122, cons202)
    rule153 = ReplacementRule(pattern153, replacement153)

    pattern154 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons21, cons25, cons40, cons203, cons204)
    rule154 = ReplacementRule(pattern154, replacement154)

    pattern155 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons149, cons122, cons184, cons205, cons206)
    rule155 = ReplacementRule(pattern155, replacement155)

    pattern156 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons149, cons122, cons207)
    rule156 = ReplacementRule(pattern156, replacement156)

    pattern157 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons21, cons25, cons149, cons203, cons204, cons208)
    rule157 = ReplacementRule(pattern157, replacement157)

    pattern158 = Pattern(Integral((e_ + u_*WC('f', S(1)))**WC('p', S(1))*(u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons19, cons4, cons5, cons70, cons71)
    rule158 = ReplacementRule(pattern158, replacement158)

    pattern159 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons209)
    rule159 = ReplacementRule(pattern159, replacement159)

    pattern160 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons74, cons68, cons212)
    rule160 = ReplacementRule(pattern160, replacement160)

    pattern161 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons95, cons96, cons91)
    rule161 = ReplacementRule(pattern161, replacement161)

    pattern162 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons213)
    rule162 = ReplacementRule(pattern162, replacement162)

    pattern163 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons214, cons68, cons215)
    rule163 = ReplacementRule(pattern163, replacement163)

    pattern164 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons216, cons215)
    rule164 = ReplacementRule(pattern164, replacement164)

    pattern165 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons217)
    rule165 = ReplacementRule(pattern165, replacement165)

    pattern166 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons95, cons96, cons90, cons20)
    rule166 = ReplacementRule(pattern166, replacement166)

    pattern167 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons5, cons95, cons96, cons90, cons172)
    rule167 = ReplacementRule(pattern167, replacement167)

    pattern168 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons33, cons96, cons20)
    rule168 = ReplacementRule(pattern168, replacement168)

    pattern169 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons33, cons96, cons172)
    rule169 = ReplacementRule(pattern169, replacement169)

    pattern170 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons33, cons170, cons146, cons20)
    rule170 = ReplacementRule(pattern170, replacement170)

    pattern171 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons33, cons170, cons146, cons172)
    rule171 = ReplacementRule(pattern171, replacement171)

    pattern172 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons196, cons68, cons197)
    rule172 = ReplacementRule(pattern172, replacement172)

    pattern173 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule173 = ReplacementRule(pattern173, replacement173)

    pattern174 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons4, cons5, cons219)
    rule174 = ReplacementRule(pattern174, replacement174)

    pattern175 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons190, cons220)
    rule175 = ReplacementRule(pattern175, replacement175)

    pattern176 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons221)
    rule176 = ReplacementRule(pattern176, replacement176)

    pattern177 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons52, cons13, cons147)
    rule177 = ReplacementRule(pattern177, replacement177)

    pattern178 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule178 = ReplacementRule(pattern178, replacement178)

    pattern179 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons82)
    rule179 = ReplacementRule(pattern179, replacement179)

    pattern180 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule180 = ReplacementRule(pattern180, replacement180)

    pattern181 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule181 = ReplacementRule(pattern181, replacement181)

    pattern182 = Pattern(Integral(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule182 = ReplacementRule(pattern182, replacement182)

    pattern183 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule183 = ReplacementRule(pattern183, replacement183)

    pattern184 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule184 = ReplacementRule(pattern184, replacement184)

    pattern185 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule185 = ReplacementRule(pattern185, replacement185)

    pattern186 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons218)
    rule186 = ReplacementRule(pattern186, replacement186)

    pattern187 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons222)
    rule187 = ReplacementRule(pattern187, replacement187)

    pattern188 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons223, cons221)
    rule188 = ReplacementRule(pattern188, replacement188)

    pattern189 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons52, cons224)
    rule189 = ReplacementRule(pattern189, replacement189)

    pattern190 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(u_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(u_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons19, cons4, cons5, cons52, cons70, cons71)
    rule190 = ReplacementRule(pattern190, replacement190)

    pattern191 = Pattern(Integral(((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_*WC('i', S(1)))**r_, x_), cons2, cons3, cons8, cons29, cons50, cons127, cons210, cons211, cons226, cons19, cons4, cons5, cons52, cons54, cons225)
    rule191 = ReplacementRule(pattern191, replacement191)
    return [rule39, rule40, rule41, rule42, rule43, rule44, rule45, rule46, rule47, rule48, rule49, rule50, rule51, rule52, rule53, rule54, rule55, rule56, rule57, rule58, rule59, rule60, rule61, rule62, rule63, rule64, rule65, rule66, rule67, rule68, rule69, rule70, rule71, rule72, rule73, rule74, rule75, rule76, rule77, rule78, rule79, rule80, rule81, rule82, rule83, rule84, rule85, rule86, rule87, rule88, rule89, rule90, rule91, rule92, rule93, rule94, rule95, rule96, rule97, rule98, rule99, rule100, rule101, rule102, rule103, rule104, rule105, rule106, rule107, rule108, rule109, rule110, rule111, rule112, rule113, rule114, rule115, rule116, rule117, rule118, rule119, rule120, rule121, rule122, rule123, rule124, rule125, rule126, rule127, rule128, rule129, rule130, rule131, rule132, rule133, rule134, rule135, rule136, rule137, rule138, rule139, rule140, rule141, rule142, rule143, rule144, rule145, rule146, rule147, rule148, rule149, rule150, rule151, rule152, rule153, rule154, rule155, rule156, rule157, rule158, rule159, rule160, rule161, rule162, rule163, rule164, rule165, rule166, rule167, rule168, rule169, rule170, rule171, rule172, rule173, rule174, rule175, rule176, rule177, rule178, rule179, rule180, rule181, rule182, rule183, rule184, rule185, rule186, rule187, rule188, rule189, rule190, rule191, ]





def replacement39(x):
    return Simp(a_*x, x)


def replacement40(x):
    return Simp(log(x), x)


def replacement41(m, x):
    return Simp(x**(m + S(1))/(m + S(1)), x)


def replacement42(a, b, x):
    return Simp(log(RemoveContent(a + b*x, x))/b, x)


def replacement43(a, b, m, x):
    return Simp((a + b*x)**(m + S(1))/(b*(m + S(1))), x)


def replacement44(a, b, m, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x)**m, x), x, u), x)


def replacement45(a, b, c, d, x):
    return Int(S(1)/(a*c + b*d*x**S(2)), x)


def replacement46(a, b, c, d, x):
    return Dist(b/(-a*d + b*c), Int(S(1)/(a + b*x), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(c + d*x), x), x)


def replacement47(a, b, c, d, m, n, x):
    return Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)), x)


def replacement48(a, b, c, d, m, x):
    return Dist(S(2)*a*c*m/(S(2)*m + S(1)), Int((a + b*x)**(m + S(-1))*(c + d*x)**(m + S(-1)), x), x) + Simp(x*(a + b*x)**m*(c + d*x)**m/(S(2)*m + S(1)), x)


def replacement49(a, b, c, d, x):
    return Simp(x/(a*c*sqrt(a + b*x)*sqrt(c + d*x)), x)


def replacement50(a, b, c, d, m, x):
    return Dist((S(2)*m + S(3))/(S(2)*a*c*(m + S(1))), Int((a + b*x)**(m + S(1))*(c + d*x)**(m + S(1)), x), x) - Simp(x*(a + b*x)**(m + S(1))*(c + d*x)**(m + S(1))/(S(2)*a*c*(m + S(1))), x)


def replacement51(a, b, c, d, m, x):
    return Int((a*c + b*d*x**S(2))**m, x)


def replacement52(a, b, c, d, x):
    return Simp(acosh(b*x/a)/b, x)


def replacement53(a, b, c, d, x):
    return Dist(S(2), Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)), x)


def replacement54(a, b, c, d, m, x):
    return Dist((a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m)), Int((a*c + b*d*x**S(2))**m, x), x)


def replacement55(a, b, c, d, x):
    return Dist((-a*d + b*c)/(S(2)*b), Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x), x) + Simp(-S(2)/(b*(a + b*x)**(S(1)/4)*(c + d*x)**(S(1)/4)), x)


def replacement56(a, b, c, d, x):
    return -Dist(d/(S(5)*b), Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x), x) + Simp(-S(4)/(S(5)*b*(a + b*x)**(S(5)/4)*(c + d*x)**(S(1)/4)), x)


def replacement57(a, b, c, d, m, n, x):
    return Dist(S(2)*c*n/(m + n + S(1)), Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))), x)


def replacement58(a, b, c, d, m, n, x):
    return Dist((m + n + S(2))/(S(2)*a*(m + S(1))), Int((a + b*x)**(m + S(1))*(c + d*x)**n, x), x) - Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/(S(2)*a*d*(m + S(1))), x)


def replacement59(a, b, c, d, m, n, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)


def replacement60(a, b, c, d, m, n, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)


def replacement61(a, b, c, d, n, x):
    return Dist((-a*d + b*c)/b, Int((c + d*x)**(n + S(-1))/(a + b*x), x), x) + Simp((c + d*x)**n/(b*n), x)


def replacement62(a, b, c, d, n, x):
    return Dist(b/(-a*d + b*c), Int((c + d*x)**(n + S(1))/(a + b*x), x), x) - Simp((c + d*x)**(n + S(1))/((n + S(1))*(-a*d + b*c)), x)


def With63(a, b, c, d, x):
    q = Rt((-a*d + b*c)/b, S(3))
    return Dist(S(3)/(S(2)*b), Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3)), x) - Dist(S(3)/(S(2)*b*q), Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3)), x) - Simp(log(RemoveContent(a + b*x, x))/(S(2)*b*q), x)


def With64(a, b, c, d, x):
    q = Rt(-(-a*d + b*c)/b, S(3))
    return Dist(S(3)/(S(2)*b), Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3)), x) - Dist(S(3)/(S(2)*b*q), Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3)), x) + Simp(log(RemoveContent(a + b*x, x))/(S(2)*b*q), x)


def With65(a, b, c, d, x):
    q = Rt((-a*d + b*c)/b, S(3))
    return -Dist(S(3)/(S(2)*b*q**S(2)), Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3)), x) - Dist(S(3)/(S(2)*b*q), Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3)), x) - Simp(log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2)), x)


def With66(a, b, c, d, x):
    q = Rt(-(-a*d + b*c)/b, S(3))
    return Dist(S(3)/(S(2)*b*q**S(2)), Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3)), x) + Dist(S(3)/(S(2)*b*q), Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3)), x) - Simp(log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2)), x)


def With67(a, b, c, d, n, x):
    p = Denominator(n)
    return Dist(p, Subst(Int(x**(p*(n + S(1)) + S(-1))/(a*d - b*c + b*x**p), x), x, (c + d*x)**(S(1)/p)), x)


def replacement68(c, d, n, x):
    return -Simp((c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), S(1) + d*x/c)/(c*(n + S(1))), x)


def replacement69(a, b, c, d, n, x):
    return -Simp((c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), TogetherSimplify(b*(c + d*x)/(-a*d + b*c)))/((n + S(1))*(-a*d + b*c)), x)


def replacement70(a, b, c, d, m, n, x):
    return -Dist(d*n/(b*(m + S(1))), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + S(1))), x)


def replacement71(a, b, c, d, m, n, x):
    return -Dist(d*(m + n + S(2))/((m + S(1))*(-a*d + b*c)), Int((a + b*x)**(m + S(1))*(c + d*x)**n, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)), x)


def replacement72(a, b, c, d, m, n, x):
    return Dist(n*(-a*d + b*c)/(b*(m + n + S(1))), Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))), x)


def replacement73(a, b, c, d, x):
    return Int(S(1)/sqrt(a*c - b**S(2)*x**S(2) - b*x*(a - c)), x)


def replacement74(a, b, c, d, x):
    return Dist(S(2)/sqrt(b), Subst(Int(S(1)/sqrt(-a*d + b*c + d*x**S(2)), x), x, sqrt(a + b*x)), x)


def replacement75(a, b, c, d, x):
    return Dist(S(2)/b, Subst(Int(S(1)/sqrt(-a + c + x**S(2)), x), x, sqrt(a + b*x)), x)


def replacement76(a, b, c, d, x):
    return Dist(S(2), Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)), x)


def replacement77(a, b, c, d, m, x):
    return Dist((a + b*x)**m*(c + d*x)**m*(a*c + b*d*x**S(2) + x*(a*d + b*c))**(-m), Int((a*c + b*d*x**S(2) + x*(a*d + b*c))**m, x), x)


def With78(a, b, c, d, x):
    q = Rt(d/b, S(3))
    return -Simp(q*log(c + d*x)/(S(2)*d), x) - Simp(S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(-1))/(S(2)*d), x) - Simp(sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d, x)


def With79(a, b, c, d, x):
    q = Rt(-d/b, S(3))
    return Simp(q*log(c + d*x)/(S(2)*d), x) + Simp(S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(1))/(S(2)*d), x) + Simp(sqrt(S(3))*q*ArcTan(-S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d, x)


def With80(a, b, c, d, m, n, x):
    p = Denominator(m)
    return Dist(p, Subst(Int(x**(p*(m + S(1)) + S(-1))/(b - d*x**p), x), x, (a + b*x)**(S(1)/p)*(c + d*x)**(-S(1)/p)), x)


def With81(a, b, c, d, m, n, x):
    p = Denominator(m)
    return Dist(p/b, Subst(Int(x**(p*(m + S(1)) + S(-1))*(-a*d/b + c + d*x**p/b)**n, x), x, (a + b*x)**(S(1)/p)), x)


def replacement82(a, b, c, d, m, n, x):
    return -Dist(d*(m + n + S(2))/((m + S(1))*(-a*d + b*c)), Int((a + b*x)**(m + S(1))*(c + d*x)**n, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)), x)


def replacement83(b, c, d, m, n, x):
    return Simp(c**n*(b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*x/c)/(b*(m + S(1))), x)


def replacement84(b, c, d, m, n, x):
    return Simp((-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*Hypergeometric2F1(-m, n + S(1), n + S(2), S(1) + d*x/c)/(d*(n + S(1))), x)


def replacement85(b, c, d, m, n, x):
    return Dist(c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n), Int((b*x)**m*(S(1) + d*x/c)**n, x), x)


def replacement86(b, c, d, m, n, x):
    return Dist((b*x)**FracPart(m)*(-b*c/d)**IntPart(m)*(-d*x/c)**(-FracPart(m)), Int((-d*x/c)**m*(c + d*x)**n, x), x)


def replacement87(a, b, c, d, m, n, x):
    return Simp(b**(-n + S(-1))*(a + b*x)**(m + S(1))*(-a*d + b*c)**n*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(m + S(1)), x)


def replacement88(a, b, c, d, m, n, x):
    return Simp((b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(b*(m + S(1))), x)


def replacement89(a, b, c, d, m, n, x):
    return Dist((b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n), Int((a + b*x)**m*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x), x)


def replacement90(a, b, c, d, m, n, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, u), x)


def replacement91(a, b, c, d, e, f, m, n, p, x):
    return Int((e + f*x)**p*(a*c + b*d*x**S(2))**m, x)


def replacement92(a, b, c, d, e, f, n, p, x):
    return Simp(b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))), x)


def replacement93(a, b, d, e, f, n, p, x):
    return Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x)


def replacement94(a, b, d, e, f, n, p, x):
    return Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x)


def replacement95(a, b, c, d, e, f, n, p, x):
    return Int(ExpandIntegrand((a + b*x)*(c + d*x)**n*(e + f*x)**p, x), x)


def replacement96(a, b, c, d, e, f, n, p, x):
    return Dist(b/f, Int((c + d*x)**n*(e + f*x)**(p + S(1)), x), x) - Simp((c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)), x)


def replacement97(a, b, c, d, e, f, n, p, x):
    return -Dist((a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))/(f*(p + S(1))*(c*f - d*e)), Int((c + d*x)**n*(e + f*x)**(p + S(1)), x), x) - Simp((c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)), x)


def replacement98(a, b, c, d, e, f, n, p, x):
    return -Dist((a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))/(f*(p + S(1))*(c*f - d*e)), Int((c + d*x)**n*(e + f*x)**(p + S(1)), x), x) - Simp((c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)), x)


def replacement99(a, b, c, d, e, f, n, p, x):
    return Dist((a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))/(d*f*(n + p + S(2))), Int((c + d*x)**n*(e + f*x)**p, x), x) + Simp(b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))), x)


def replacement100(a, b, c, d, e, f, n, p, x):
    return Simp(b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(S(2)*a*d*f*(n + p + S(3)) + b*d*f*x*(n + p + S(2)) - b*(c*f*(p + S(2)) + d*e*(n + S(2))))/(d**S(2)*f**S(2)*(n + p + S(2))*(n + p + S(3))), x)


def replacement101(a, b, c, d, f, m, n, p, x):
    return Dist(a, Int((f*x)**p*(a + b*x)**n*(c + d*x)**n, x), x) + Dist(b/f, Int((f*x)**(p + S(1))*(a + b*x)**n*(c + d*x)**n, x), x)


def replacement102(a, b, c, d, e, f, p, x):
    return Int(ExpandIntegrand((e + f*x)**p/((a + b*x)*(c + d*x)), x), x)


def replacement103(a, b, c, d, e, f, p, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int((e + f*x)**(p + S(-1))/(a + b*x), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int((e + f*x)**(p + S(-1))/(c + d*x), x), x)


def replacement104(a, b, c, d, e, f, p, x):
    return Dist(S(1)/(b*d), Int((e + f*x)**(p + S(-2))*(-a*c*f**S(2) + b*d*e**S(2) + f*x*(-a*d*f - b*c*f + S(2)*b*d*e))/((a + b*x)*(c + d*x)), x), x) + Simp(f*(e + f*x)**(p + S(-1))/(b*d*(p + S(-1))), x)


def replacement105(a, b, c, d, e, f, p, x):
    return Dist(S(1)/((-a*f + b*e)*(-c*f + d*e)), Int((e + f*x)**(p + S(1))*(-a*d*f - b*c*f + b*d*e - b*d*f*x)/((a + b*x)*(c + d*x)), x), x) + Simp(f*(e + f*x)**(p + S(1))/((p + S(1))*(-a*f + b*e)*(-c*f + d*e)), x)


def replacement106(a, b, c, d, e, f, p, x):
    return Dist(b/(-a*d + b*c), Int((e + f*x)**p/(a + b*x), x), x) - Dist(d/(-a*d + b*c), Int((e + f*x)**p/(c + d*x), x), x)


def replacement107(a, b, c, d, e, f, n, p, x):
    return Int(ExpandIntegrand((e + f*x)**FractionalPart(p), (c + d*x)**n*(e + f*x)**IntegerPart(p)/(a + b*x), x), x)


def replacement108(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x)


def replacement109(a, b, c, d, e, f, n, p, x):
    return -Dist(S(1)/(d**S(2)*(n + S(1))*(-c*f + d*e)), Int((c + d*x)**(n + S(1))*(e + f*x)**p*Simp(a**S(2)*d**S(2)*f*(n + p + S(2)) - S(2)*a*b*d*(c*f*(p + S(1)) + d*e*(n + S(1))) + b**S(2)*c*(c*f*(p + S(1)) + d*e*(n + S(1))) - b**S(2)*d*x*(n + S(1))*(-c*f + d*e), x), x), x) + Simp((c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*d + b*c)**S(2)/(d**S(2)*(n + S(1))*(-c*f + d*e)), x)


def replacement110(a, b, c, d, e, f, n, p, x):
    return Dist(S(1)/(d*f*(n + p + S(3))), Int((c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(n + p + S(3)) + b*x*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e), x), x), x) + Simp(b*(a + b*x)*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(3))), x)


def With111(a, b, c, d, e, f, x):
    q = Rt((-c*f + d*e)/(-a*f + b*e), S(3))
    return Simp(q*log(e + f*x)/(-S(2)*c*f + S(2)*d*e), x) - Simp(S(3)*q*log(q*(a + b*x)**(S(1)/3) - (c + d*x)**(S(1)/3))/(-S(2)*c*f + S(2)*d*e), x) - Simp(sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(-c*f + d*e), x)


def replacement112(a, b, c, d, e, f, x):
    return Dist(b*f, Subst(Int(S(1)/(b*f**S(2)*x**S(2) + d*(-a*f + b*e)**S(2)), x), x, sqrt(a + b*x)*sqrt(c + d*x)), x)


def With113(a, b, c, d, e, f, m, n, x):
    q = Denominator(m)
    return Dist(q, Subst(Int(x**(q*(m + S(1)) + S(-1))/(-a*f + b*e - x**q*(-c*f + d*e)), x), x, (a + b*x)**(S(1)/q)*(c + d*x)**(-S(1)/q)), x)


def replacement114(a, b, c, d, e, f, m, n, p, x):
    return -Dist(n*(-c*f + d*e)/((m + S(1))*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)), x)


def replacement115(a, b, c, d, e, f, m, n, p, x):
    return Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement116(a, b, c, d, e, f, m, n, p, x):
    return Dist((a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x), x) + Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement117(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(1)/(b*(m + S(1))), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(-1))*Simp(c*f*p + d*e*n + d*f*x*(n + p), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p/(b*(m + S(1))), x)


def replacement118(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/(b*(m + S(1))*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-2))*(e + f*x)**p*Simp(a*d*(c*f*(p + S(1)) + d*e*(n + S(-1))) + b*c*(-c*f*(m + p + S(2)) + d*e*(m - n + S(2))) + d*x*(a*d*f*(n + p) + b*(-c*f*(m + n + p + S(1)) + d*e*(m + S(1)))), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(1))*(-a*d + b*c)/(b*(m + S(1))*(-a*f + b*e)), x)


def replacement119(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(1)/((m + S(1))*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(c*f*(m + p + S(2)) + d*e*n + d*f*x*(m + n + p + S(2)), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)), x)


def replacement120(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/(d*f*(m + n + p + S(1))), Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x), x) + Simp(b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))), x)


def replacement121(a, b, c, d, e, f, m, n, p, x):
    return -Dist(S(1)/(f*(m + n + p + S(1))), Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(a*n*(-c*f + d*e) + c*m*(-a*f + b*e) + x*(b*n*(-c*f + d*e) + d*m*(-a*f + b*e)), x), x), x) + Simp((a + b*x)**m*(c + d*x)**n*(e + f*x)**(p + S(1))/(f*(m + n + p + S(1))), x)


def replacement122(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/(d*f*(m + n + p + S(1))), Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x), x) + Simp(b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))), x)


def replacement123(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x), x) + Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement124(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x), x) + Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement125(a, b, c, d, e, f, m, n, x):
    return Dist(b/f, Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x), x) - Dist((-a*f + b*e)/f, Int((a + b*x)**(m + S(-1))*(c + d*x)**n/(e + f*x), x), x)


def replacement126(a, b, c, d, e, f, x):
    return Dist(S(-4), Subst(Int(x**S(2)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)), x)


def replacement127(a, b, c, d, e, f, x):
    return Dist(sqrt(-f*(c + d*x)/(-c*f + d*e))/sqrt(c + d*x), Int(S(1)/((a + b*x)*(e + f*x)**(S(1)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x), x)


def replacement128(a, b, c, d, e, f, x):
    return Dist(S(-4), Subst(Int(S(1)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)), x)


def replacement129(a, b, c, d, e, f, x):
    return Dist(sqrt(-f*(c + d*x)/(-c*f + d*e))/sqrt(c + d*x), Int(S(1)/((a + b*x)*(e + f*x)**(S(3)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x), x)


def replacement130(b, c, d, e, f, x):
    return Simp(S(2)*sqrt(e)*EllipticE(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/b, x)


def replacement131(b, c, d, e, f, x):
    return Dist(sqrt(-b*x)/sqrt(b*x), Int(sqrt(e + f*x)/(sqrt(-b*x)*sqrt(c + d*x)), x), x)


def replacement132(b, c, d, e, f, x):
    return Dist(sqrt(S(1) + d*x/c)*sqrt(e + f*x)/(sqrt(S(1) + f*x/e)*sqrt(c + d*x)), Int(sqrt(S(1) + f*x/e)/(sqrt(b*x)*sqrt(S(1) + d*x/c)), x), x)


def replacement133(a, b, c, d, e, f, x):
    return Simp(S(2)*EllipticE(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*f + b*e)/d, S(2))/b, x)


def replacement134(a, b, c, d, e, f, x):
    return Dist(sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)/(sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)), Int(sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))), x), x)


def replacement135(b, c, d, e, f, x):
    return Simp(S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)), x)


def replacement136(b, c, d, e, f, x):
    return Simp(S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)), x)


def replacement137(b, c, d, e, f, x):
    return Dist(sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)/(sqrt(c + d*x)*sqrt(e + f*x)), Int(S(1)/(sqrt(b*x)*sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)), x), x)


def replacement138(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b, x)


def replacement139(a, b, c, d, e, f, x):
    return Simp(S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b, x)


def replacement140(a, b, c, d, e, f, x):
    return Dist(sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))/(sqrt(c + d*x)*sqrt(e + f*x)), Int(S(1)/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))*sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))), x), x)


def With141(a, b, c, d, e, f, x):
    q = Rt(b*(-a*f + b*e)/(-a*d + b*c)**S(2), S(3))
    return -Simp(log(a + b*x)/(S(2)*q*(-a*d + b*c)), x) + Simp(S(3)*log(q*(c + d*x)**(S(2)/3) - (e + f*x)**(S(1)/3))/(S(4)*q*(-a*d + b*c)), x) - Simp(sqrt(S(3))*ArcTan(S(2)*sqrt(S(3))*q*(c + d*x)**(S(2)/3)/(S(3)*(e + f*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*q*(-a*d + b*c)), x)


def replacement142(a, b, c, d, e, f, m, x):
    return Dist(f/(S(6)*(m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(a*d*(S(3)*m + S(1)) - S(3)*b*c*(S(3)*m + S(5)) - S(2)*b*d*x*(S(3)*m + S(7)))/((c + d*x)**(S(1)/3)*(e + f*x)**(S(1)/3)), x), x) + Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(S(2)/3)*(e + f*x)**(S(2)/3)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement143(a, b, c, d, f, m, n, p, x):
    return Int((f*x)**p*(a*c + b*d*x**S(2))**m, x)


def replacement144(a, b, c, d, f, m, n, p, x):
    return Dist((a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m)), Int((f*x)**p*(a*c + b*d*x**S(2))**m, x), x)


def replacement145(a, b, c, d, f, m, n, p, x):
    return Int(ExpandIntegrand((f*x)**p*(a + b*x)**n*(c + d*x)**n, (a + b*x)**(m - n), x), x)


def replacement146(a, b, c, d, e, f, m, n, p, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x)


def replacement147(a, b, c, d, e, f, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x), x) + Simp(b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement148(a, b, c, d, e, f, m, n, p, x):
    return Simp((a + b*x)**(m + S(1))*(e + f*x)**(-m + S(-1))*(-a*d + b*c)**n*(-a*f + b*e)**(-n + S(-1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/(m + S(1)), x)


def replacement149(a, b, c, d, e, f, m, n, p, x):
    return Simp(((c + d*x)*(-a*f + b*e)/((e + f*x)*(-a*d + b*c)))**(-n)*(a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/((m + S(1))*(-a*f + b*e)), x)


def replacement150(b, c, d, e, f, m, n, p, x):
    return Simp(c**n*e**p*(b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*x/c, -f*x/e)/(b*(m + S(1))), x)


def replacement151(b, c, d, e, f, m, n, p, x):
    return Simp((d/(-c*f + d*e))**(-p)*(-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*AppellF1(n + S(1), -m, -p, n + S(2), S(1) + d*x/c, -f*(c + d*x)/(-c*f + d*e))/(d*(n + S(1))), x)


def replacement152(b, c, d, e, f, m, n, p, x):
    return Dist(c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n), Int((b*x)**m*(S(1) + d*x/c)**n*(e + f*x)**p, x), x)


def replacement153(a, b, c, d, e, f, m, n, p, x):
    return Simp(b**(-p + S(-1))*(b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*(-a*f + b*e)**p*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(m + S(1)), x)


def replacement154(a, b, c, d, e, f, m, n, p, x):
    return Dist((b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n), Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x), x)


def replacement155(a, b, c, d, e, f, m, n, p, x):
    return Simp((b/(-a*d + b*c))**(-n)*(b/(-a*f + b*e))**(-p)*(a + b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(b*(m + S(1))), x)


def replacement156(a, b, c, d, e, f, m, n, p, x):
    return Dist((b/(-a*f + b*e))**(-IntPart(p))*(b*(e + f*x)/(-a*f + b*e))**(-FracPart(p))*(e + f*x)**FracPart(p), Int((a + b*x)**m*(c + d*x)**n*(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))**p, x), x)


def replacement157(a, b, c, d, e, f, m, n, p, x):
    return Dist((b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n), Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x), x)


def replacement158(a, b, c, d, e, f, m, n, p, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, u), x)


def replacement159(a, b, c, d, e, f, g, h, m, n, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)*(g + h*x), x), x)


def replacement160(a, b, c, d, e, f, g, h, m, n, x):
    return Dist((a*d*f*h*m + b*(-c*f*h*(m + S(2)) + d*(e*h + f*g)))/(b**S(2)*d), Int((a + b*x)**(m + S(1))*(c + d*x)**n, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(2)*d*f*h*m - a*b*(-c*f*h*(m + S(1)) + d*(e*h + f*g)) + b**S(2)*d*e*g + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)), x)


def replacement161(a, b, c, d, e, f, g, h, m, n, x):
    return -Dist((a**S(2)*d**S(2)*f*h*(n**S(2) + S(3)*n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m**S(2) + S(3)*m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m**S(2) + m*(S(2)*n + S(5)) + n**S(2) + S(5)*n + S(6))))/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1)), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*c*d*f*h*(n + S(1)) + a*b*(c**S(2)*f*h*(m + S(1)) - c*d*(e*h + f*g)*(m + n + S(2)) + d**S(2)*e*g*(m + S(1))) + b**S(2)*c*d*e*g*(n + S(1)) + x*(a**S(2)*d**S(2)*f*h*(n + S(1)) - a*b*d**S(2)*(n + S(1))*(e*h + f*g) + b**S(2)*(c**S(2)*f*h*(m + S(1)) - c*d*(m + S(1))*(e*h + f*g) + d**S(2)*e*g*(m + n + S(2)))))/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)), x)


def replacement162(a, b, c, d, e, f, g, h, m, n, x):
    return Dist(-d*(m + n + S(3))*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)) + f*h/b**S(2), Int((a + b*x)**(m + S(2))*(c + d*x)**n, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(3)*d*f*h*(n + S(2)) - a**S(2)*b*(c*f*h*m - d*(e*h + f*g)*(m + n + S(3))) - a*b**S(2)*(c*(e*h + f*g) + d*e*g*(S(2)*m + n + S(4))) + b**S(3)*c*e*g*(m + S(2)) + b*x*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2)))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)), x)


def replacement163(a, b, c, d, e, f, g, h, m, n, x):
    return -Dist((a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))), Int((a + b*x)**(m + S(1))*(c + d*x)**n, x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*d*f*h*(n + S(2)) + a*b*(c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*d*e*g*(m + n + S(3)) + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))), x)


def replacement164(a, b, c, d, e, f, g, h, m, n, x):
    return Dist((a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))), Int((a + b*x)**m*(c + d*x)**n, x), x) - Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a*d*f*h*(n + S(2)) + b*c*f*h*(m + S(2)) - b*d*f*h*x*(m + n + S(2)) - b*d*(e*h + f*g)*(m + n + S(3)))/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))), x)


def replacement165(a, b, c, d, e, f, g, h, m, n, p, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x), x), x)


def replacement166(a, b, c, d, e, f, g, h, m, n, p, x):
    return -Dist(S(1)/(b*(m + S(1))*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)), x)


def replacement167(a, b, c, d, e, f, g, h, m, n, p, x):
    return -Dist(S(1)/(b*(m + S(1))*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)), x)


def replacement168(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement169(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement170(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(S(1)/(d*f*(m + n + p + S(2))), Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x), x) + Simp(h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))), x)


def replacement171(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(S(1)/(d*f*(m + n + p + S(2))), Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x), x) + Simp(h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))), x)


def replacement172(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(S(1)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x), x) + Simp((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)), x)


def replacement173(a, b, c, d, e, f, g, h, p, x):
    return Dist((-a*h + b*g)/(-a*d + b*c), Int((e + f*x)**p/(a + b*x), x), x) - Dist((-c*h + d*g)/(-a*d + b*c), Int((e + f*x)**p/(c + d*x), x), x)


def replacement174(a, b, c, d, e, f, g, h, n, p, x):
    return Dist(h/b, Int((c + d*x)**n*(e + f*x)**p, x), x) + Dist((-a*h + b*g)/b, Int((c + d*x)**n*(e + f*x)**p/(a + b*x), x), x)


def replacement175(a, b, c, d, e, f, g, h, x):
    return Dist(h/f, Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)), x), x) + Dist((-e*h + f*g)/f, Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)), x), x)


def replacement176(a, b, c, d, e, f, g, h, m, n, p, x):
    return Dist(h/b, Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x), x) + Dist((-a*h + b*g)/b, Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x)


def replacement177(a, b, c, d, e, f, g, h, p, q, x):
    return Dist((-a*f + b*e)/(-a*d + b*c), Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(a + b*x), x), x) - Dist((-c*f + d*e)/(-a*d + b*c), Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(c + d*x), x), x)


def replacement178(a, b, c, d, e, f, g, h, x):
    return Simp(-S(2)*sqrt(d*(e + f*x)/(-c*f + d*e))*sqrt(d*(g + h*x)/(-c*h + d*g))*EllipticPi(-b*(-c*f + d*e)/(f*(-a*d + b*c)), asin(sqrt(-f/(-c*f + d*e))*sqrt(c + d*x)), h*(-c*f + d*e)/(f*(-c*h + d*g)))/(sqrt(-f/(-c*f + d*e))*sqrt(e + f*x)*sqrt(g + h*x)*(-a*d + b*c)), x)


def replacement179(a, b, c, d, e, f, g, h, n, x):
    return Int(ExpandIntegrand(S(1)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), (c + d*x)**(n + S(1)/2)/(a + b*x), x), x)


def replacement180(a, b, c, d, e, f, g, h, x):
    return Dist(b**(S(-2)), Int((-a*f*h + b*e*h + b*f*g + b*f*h*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x), x) + Dist((-a*f + b*e)*(-a*h + b*g)/b**S(2), Int(S(1)/((a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x), x)


def replacement181(a, b, c, d, e, f, g, h, x):
    return Dist(-S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)), Subst(Int(S(1)/(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x)), x)


def replacement182(a, b, c, d, e, f, g, h, x):
    return Dist(-S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*(-c*h + d*g)/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)**S(2)), Subst(Int(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))/sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1)), x), x, sqrt(g + h*x)/sqrt(a + b*x)), x)


def replacement183(a, b, c, d, e, f, g, h, x):
    return Dist(S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)), Subst(Int(S(1)/((-b*x**S(2) + h)*sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x)), x)


def replacement184(a, b, c, d, e, f, g, h, x):
    return Dist(b/(-a*d + b*c), Int(sqrt(c + d*x)/((a + b*x)**(S(3)/2)*sqrt(e + f*x)*sqrt(g + h*x)), x), x) - Dist(d/(-a*d + b*c), Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x), x)


def replacement185(a, b, c, d, e, f, g, h, x):
    return Dist((a*d*f*h - b*(-c*f*h + d*e*h + d*f*g))/(S(2)*f**S(2)*h), Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)), x), x) + Dist((-c*f + d*e)*(-S(2)*a*f*h + b*e*h + b*f*g)/(S(2)*f**S(2)*h), Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x), x) - Dist((-c*f + d*e)*(-e*h + f*g)/(S(2)*f*h), Int(sqrt(a + b*x)/(sqrt(c + d*x)*(e + f*x)**(S(3)/2)*sqrt(g + h*x)), x), x) + Simp(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)/(h*sqrt(e + f*x)), x)


def replacement186(a, b, c, d, e, f, g, h, x):
    return Dist(b/d, Int(sqrt(a + b*x)*sqrt(c + d*x)/(sqrt(e + f*x)*sqrt(g + h*x)), x), x) - Dist((-a*d + b*c)/d, Int(sqrt(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x), x)


def replacement187(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x)


def replacement188(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Dist(h/b, Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x), x) + Dist((-a*h + b*g)/b, Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x), x)


def replacement189(a, b, c, d, e, f, g, h, m, n, p, q, x):
    return Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x)


def replacement190(a, b, c, d, e, f, g, h, m, n, p, q, u, x):
    return Dist(S(1)/Coefficient(u, x, S(1)), Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, u), x)


def replacement191(a, b, c, d, e, f, g, h, i, m, n, p, q, r, x):
    return Dist((i*(a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q)**r*(a + b*x)**(-m*r)*(c + d*x)**(-n*r)*(e + f*x)**(-p*r)*(g + h*x)**(-q*r), Int((a + b*x)**(m*r)*(c + d*x)**(n*r)*(e + f*x)**(p*r)*(g + h*x)**(q*r), x), x)
